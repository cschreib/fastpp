#include "fast++.hpp"

extern "C" {
    #include "tinyexpr.h"
}

double expr_step(double a) {
    return (a >= 0 ? 1.0 : 0.0);
}

double expr_min(double a, double b) {
    return std::min(a, b);
}

double expr_max(double a, double b) {
    return std::max(a, b);
}

bool gridder_t::tinyexpr_wrapper::compile(const std::string& sexpr, const vec1s& params) {
    uint_t nfuncparam = params.size();
    uint_t nfunc = 3;

    vars_glue = new te_variable[nfuncparam+nfunc];
    for (uint_t p : range(nfuncparam+nfunc)) {
        vars_glue[p].context = nullptr;
    }

    vars.resize(nfuncparam);

    // Time
    vars_glue[0].name = "t";
    vars_glue[0].address = &vars[0];
    vars_glue[0].type = TE_VARIABLE;

    // Custom parameters
    for (uint_t p : range(nfuncparam)) {
        vars_glue[p].name = params[p].c_str();
        vars_glue[p].address = &vars[p];
        vars_glue[p].type = TE_VARIABLE;
    }

    // Custom functions
    vars_glue[nfuncparam+0].name = "step";
    vars_glue[nfuncparam+0].address = (void*)(&expr_step);
    vars_glue[nfuncparam+0].type = TE_FUNCTION1 | TE_FLAG_PURE;
    vars_glue[nfuncparam+1].name = "min";
    vars_glue[nfuncparam+1].address = (void*)(&expr_min);
    vars_glue[nfuncparam+1].type = TE_FUNCTION2 | TE_FLAG_PURE;
    vars_glue[nfuncparam+2].name = "max";
    vars_glue[nfuncparam+2].address = (void*)(&expr_max);
    vars_glue[nfuncparam+2].type = TE_FUNCTION2 | TE_FLAG_PURE;

    // Compile expression
    int err = 0;
    expr = te_compile(sexpr.c_str(), vars_glue, nfuncparam+nfunc, &err);
    if (err > 0) {
        std::string head = "could not parse SFH expression: ";
        error(head, sexpr);
        error(std::string(head.size()+err-1, ' ')+'^');
        return false;
    }

    return true;
}

double gridder_t::tinyexpr_wrapper::eval() {
    return te_eval(expr);
}

gridder_t::tinyexpr_wrapper::~tinyexpr_wrapper() {
    if (expr != nullptr) {
        te_free(expr);
        delete[] vars_glue;
    }
}

std::string gridder_t::get_library_file_ssp(uint_t im) const {
    return opts.library_dir+"ssp"+"."+opts.resolution+"/"+
        opts.library+"_"+opts.resolution+"_"+opts.name_imf+
        "_z"+replace(to_string(output.grid[grid_id::metal][im]), "0.", "");
}

void gridder_t::evaluate_sfh_custom(const vec1u& idm, const vec1d& t, vec1d& sfh) const {
    sfh_expr.vars[1] = output.grid[grid_id::age][idm[grid_id::age]];
    for (uint_t i : range(opts.custom_params)) {
        sfh_expr.vars[i+2] = output.grid[grid_id::custom+i][idm[grid_id::custom+i]];
    }

    sfh.resize(t.size());
    double nage = e10(output.grid[grid_id::age][idm[grid_id::age]]);
    for (uint_t i : range(t)) {
        sfh_expr.vars[0] = (opts.custom_sfh_lookback ? t.safe[i] : nage - t.safe[i]);
        sfh.safe[i] = sfh_expr.eval();
    }
}

bool gridder_t::build_and_send_custom(fitter_t& fitter) {
    model_id_pair m;
    m.model.flux.resize(input.lambda.size());
    m.model.props.resize(nprop);
    m.idm.resize(nparam);

    ssp_bc03 ssp;

    const vec1f& output_metal = output.grid[grid_id::metal];
    const vec1f& output_age = output.grid[grid_id::age];
    const vec1f& output_z = output.grid[grid_id::z];

    auto pg = progress_start(nmodel);
    for (uint_t im : range(output_metal)) {
        m.idm[_] = 0;
        m.idm[grid_id::metal] = im;

        // Load SSP
        std::string filename = get_library_file_ssp(im);
        if (!ssp.read(filename)) {
            return false;
        }

        // Apply velocity dispersion
        convolve_rest(ssp.lambda, ssp.sed);

        // Pre-compute dust law & IGM absorption (they don't change with SFH)
        vec1d dust_law = build_dust_law(ssp.lambda);
        vec2d igm_abs;
        if (!opts.no_igm) {
            igm_abs = build_igm_absorption(output_z, ssp.lambda);
        }

        // Function to build a model
        auto do_model = [&](model_id_pair& tm) {
            float& model_mass  = tm.model.props[prop_id::mass];
            float& model_mform = tm.model.props[prop_id::mform];
            float& model_sfr   = tm.model.props[prop_id::sfr];
            float& model_ssfr  = tm.model.props[prop_id::ssfr];

            uint_t ia = tm.idm[grid_id::age];

            const double dt = opts.custom_sfh_step;
            // Compute lookback time (t=0 is when the galaxy is observed, t>0 is in the past) in [yr]
            vec1d ltime = dt*indgen<double>(uint_t(ceil(e10(output_age[ia])/dt)+1.0));

            // Build analytic SFH
            vec1d sfh; {
                auto lock = (opts.parallel == parallel_choice::generators ?
                    std::unique_lock<std::mutex>(sfh_mutex) : std::unique_lock<std::mutex>());

                evaluate_sfh_custom(tm.idm, ltime, sfh);
            }

            // Compute SFH quantities
            if (!input.sfh_quant.empty()) {
                compute_sfh_quantities_impl(ltime, sfh, tm.model);
            }

            if (opts.sfr_avg > 0) {
                // Average SFR over the past X yr
                double t1 = min(opts.sfr_avg, ltime.back());
                model_sfr = integrate(ltime, sfh, 0.0, t1)/opts.sfr_avg;
            } else {
                // Use instantaneous SFR
                model_sfr = interpolate(sfh, ltime, 0.0);
            }

            // Integrate SFH for total mass
            double tmodel_mass  = 0.0;
            double tformed_mass = 0.0;
            ssp.integrate(ltime, sfh, [&](uint_t it, double formed) {
                tmodel_mass  += formed*ssp.mass.safe[it];
                tformed_mass += formed;
            });

            model_mass = tmodel_mass;
            model_mform = tformed_mass*ssp.mass.safe[0];

            model_ssfr = model_sfr/model_mass;

            if (opts.differential_a_v) {
                // Need to separate yound and old stars for differential attenuation

                // Integrate SFH for template
                vec1d tpl_flux_young(ssp.lambda.size());
                vec1d tpl_flux_old(ssp.lambda.size());
                double age_bc = e10(opts.log_bc_age_max); // [yr]
                ssp.integrate(ltime, sfh, [&](uint_t it, double formed) {
                    if (ssp.age.safe[it] <= age_bc) {
                        tpl_flux_young += formed*ssp.sed.safe(it,_);
                    } else {
                        tpl_flux_old += formed*ssp.sed.safe(it,_);
                    }
                });

                // The rest is not specific to the SFH, use generic code
                attenuate_and_send(fitter, pg, ssp.lambda, tpl_flux_young, tpl_flux_old,
                    dust_law, igm_abs, output_age[ia], tm.idm, tm.model);
            } else {
                // Treat all stars the same way

                // Integrate SFH for template
                vec1d tpl_flux(ssp.lambda.size());
                ssp.integrate(ltime, sfh, [&](uint_t it, double formed) {
                    tpl_flux += formed*ssp.sed.safe(it,_);
                });

                // The rest is not specific to the SFH, use generic code
                attenuate_and_send(fitter, pg, ssp.lambda, tpl_flux, dust_law, igm_abs,
                    output_age[ia], tm.idm, tm.model);
            }
        };

        thread::worker_pool<model_id_pair> pool;
        if (opts.parallel == parallel_choice::generators) {
            pool.start(opts.n_thread, do_model);
        }

        // Iterate over all models
        for (uint_t ic = 0; ic < ncustom; ++ic) {
            for (uint_t ia : range(output_age)) {
                m.idm[grid_id::age] = ia;

                if (opts.parallel == parallel_choice::generators) {
                    // Parallel
                    while (opts.max_queued_fits > 0 &&
                        pool.remaining() > opts.max_queued_fits) {
                        thread::sleep_for(1e-6);
                    }

                    pool.process(m);
                } else {
                    // Single threaded
                    do_model(m);
                }
            }

            // Go to next model
            increment_index_list(m.idm, grid_dims);
        }

        if (opts.parallel == parallel_choice::generators) {
            while (pool.remaining() != 0) {
                thread::sleep_for(1e-6);
            }

            pool.join();
        }
    }

    return true;
}

bool gridder_t::build_template_custom(uint_t iflat, vec1f& lam, vec1f& flux) const {
    vec1u idm = grid_ids(iflat);
    uint_t ia = idm[grid_id::age];
    uint_t im = idm[grid_id::metal];
    const vec1f& output_age = output.grid[grid_id::age];

    // Build analytic SFH
    double dt = opts.custom_sfh_step;
    vec1d ltime = dt*indgen<double>(uint_t(ceil(e10(output_age[ia])/dt)+1.0));
    vec1d sfh; {
        auto lock = (opts.n_thread > 1 ?
            std::unique_lock<std::mutex>(sfh_mutex) : std::unique_lock<std::mutex>());

        evaluate_sfh_custom(idm, ltime, sfh);
    }

    auto lock = (opts.n_thread > 1 ?
        std::unique_lock<std::mutex>(sed_mutex) : std::unique_lock<std::mutex>());

    if (!cached_ssp_bc03) {
        cached_ssp_bc03 = std::unique_ptr<ssp_bc03>(new ssp_bc03());
    }

    ssp_bc03* ssp = cached_ssp_bc03.get();

    // Load SSP if not cached
    std::string filename = get_library_file_ssp(im);
    if (filename != cached_library) {
        if (!ssp->read(filename)) {
            return false;
        }

        cached_library = filename;

        // Apply velocity dispersion
        convolve_rest(ssp->lambda, ssp->sed);
    }

    // Integrate SFH on local time grid
    vec1d tpl_flux(ssp->lambda.size());
    ssp->integrate(ltime, sfh, [&](uint_t it, double formed) {
        tpl_flux += formed*ssp->sed.safe(it,_);
    });

    lam = ssp->lambda;
    flux = tpl_flux;

    return true;
}

bool gridder_t::build_template_custom(uint_t iflat, vec1f& lam, vec1f& flux_young,
    vec1f& flux_old) const {

    vec1u idm = grid_ids(iflat);
    uint_t ia = idm[grid_id::age];
    uint_t im = idm[grid_id::metal];
    const vec1f& output_age = output.grid[grid_id::age];

    // Build analytic SFH
    double dt = opts.custom_sfh_step;
    vec1d ltime = dt*indgen<double>(uint_t(ceil(e10(output_age[ia])/dt)+1.0));
    vec1d sfh; {
        auto lock = (opts.n_thread > 1 ?
            std::unique_lock<std::mutex>(sfh_mutex) : std::unique_lock<std::mutex>());

        evaluate_sfh_custom(idm, ltime, sfh);
    }

    auto lock = (opts.n_thread > 1 ?
        std::unique_lock<std::mutex>(sed_mutex) : std::unique_lock<std::mutex>());

    if (!cached_ssp_bc03) {
        cached_ssp_bc03 = std::unique_ptr<ssp_bc03>(new ssp_bc03());
    }

    ssp_bc03* ssp = cached_ssp_bc03.get();

    // Load SSP if not cached
    std::string filename = get_library_file_ssp(im);
    if (filename != cached_library) {
        if (!ssp->read(filename)) {
            return false;
        }

        cached_library = filename;

        // Apply velocity dispersion
        convolve_rest(ssp->lambda, ssp->sed);
    }

    // Integrate SFH for template
    vec1d tpl_flux_young(ssp->lambda.size());
    vec1d tpl_flux_old(ssp->lambda.size());
    double age_bc = e10(opts.log_bc_age_max);
    ssp->integrate(ltime, sfh, [&](uint_t it, double formed) {
        if (ssp->age.safe[it] <= age_bc) {
            tpl_flux_young += formed*ssp->sed.safe(it,_);
        } else {
            tpl_flux_old += formed*ssp->sed.safe(it,_);
        }
    });

    lam = ssp->lambda;
    flux_young = tpl_flux_young;
    flux_old = tpl_flux_old;

    return true;
}

bool gridder_t::get_sfh_custom(uint_t iflat, const vec1d& t, vec1d& sfh,
    const std::string& type) const {

    vec1u idm = grid_ids(iflat);
    double nage = e10(output.grid[grid_id::age][idm[grid_id::age]]);
    double t_obs = e10(auniv[idm[grid_id::z]]);

    vec1d tlb = reverse(t_obs - t);

    // Evaluate SFH
    uint_t i1 = upper_bound(t, nage);
    if (i1 == npos) {
        i1 = t.size()-1;
    }

    sfh = replicate(0.0, t.size()); {
        vec1d tsfh;
        evaluate_sfh_custom(idm, t[_-i1], tsfh);
        sfh[_-i1] = tsfh;
    }

    // Load SSP (only extras) to get mass
    std::string filename = get_library_file_ssp(idm[grid_id::metal]);

    ssp_bc03 ssp;
    if (!ssp.read(filename, true)) {
        return false;
    }

    if (type == "sfr") {
        // Compute total mass at epoch of observation and normalize
        double mass = 0.0;
        ssp.integrate(tlb, sfh, [&](uint_t it, double formed) {
            mass += formed*ssp.mass.safe[it];
        });

        sfh /= mass;
    } else if (type == "mass") {
        // Integrate mass, including mass loss
        vec1d mass(t.size());

        for (uint_t i : range(t)) {
            ssp.integrate(tlb[i-_] - tlb.safe[i], sfh[i-_], [&](uint_t it, double formed) {
                mass.safe[i] += formed*ssp.mass.safe[it];
            });
        }

        // Normalize to unit mass at observation
        mass /= mass[0];

        // Return mass instead of SFR
        std::swap(mass, sfh);
    } else {
        error("unknown SFH type '", type, "'");
        return false;
    }

    sfh = reverse(sfh);

    return true;
}
