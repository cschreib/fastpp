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
    uint_t nparam = params.size();
    uint_t nfunc = 3;

    vars_glue = new te_variable[nparam+nfunc];
    for (uint_t p : range(nparam+nfunc)) {
        vars_glue[p].context = nullptr;
    }

    vars.resize(nparam);

    // Time
    vars_glue[0].name = "t";
    vars_glue[0].address = &vars[0];
    vars_glue[0].type = TE_VARIABLE;

    // Custom parameters
    for (uint_t p : range(nparam)) {
        vars_glue[p].name = params[p].c_str();
        vars_glue[p].address = &vars[p];
        vars_glue[p].type = TE_VARIABLE;
    }

    // Custom functions
    vars_glue[nparam+0].name = "step";
    vars_glue[nparam+0].address = (void*)(&expr_step);
    vars_glue[nparam+0].type = TE_FUNCTION1 | TE_FLAG_PURE;
    vars_glue[nparam+1].name = "min";
    vars_glue[nparam+1].address = (void*)(&expr_min);
    vars_glue[nparam+1].type = TE_FUNCTION2 | TE_FLAG_PURE;
    vars_glue[nparam+2].name = "max";
    vars_glue[nparam+2].address = (void*)(&expr_max);
    vars_glue[nparam+2].type = TE_FUNCTION2 | TE_FLAG_PURE;

    // Compile expression
    int err = 0;
    expr = te_compile(sexpr.c_str(), vars_glue, nparam+nfunc, &err);
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

bool ssp_bc03::read_ascii(std::string filename) {
    std::string state = "";

    try {
        std::ifstream in(filename);
        in.exceptions(in.failbit);

        state = "read number of time steps";
        uint_t ntime = 0;
        in >> ntime;

        state = "read time steps";
        age.resize(ntime);
        for (uint_t i : range(age)) {
            in >> age[i];
        }

        state = "read IMF and other parameters";
        double ml, mu;
        uint_t iseg;
        in >> ml >> mu >> iseg;
        for (uint_t i = 0; i < iseg; ++i) {
            double xx, lm, um, baux, cn, cc;
            in >> xx >> lm >> um >> baux >> cn >> cc;
        }

        state = "read additional parameters";
        double totm, totn, avs, jo, tau, id, tau1, tau2, tau3, tau4;
        in >> totm >> totn >> avs >> jo >> tau >> id >> tau1 >> tau2 >> tau3 >> tau4;

        char id2;
        in >> id2;
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        std::string id3, iop, stelib;
        std::getline(in, id3);
        std::getline(in, iop);
        std::getline(in, stelib);

        state = "read number of wavelength elements";
        uint_t nlam;
        in >> nlam;

        state = "read wavelength elements";
        lambda.resize(nlam);
        for (uint_t i : range(lambda)) {
            in >> lambda[i];
        }

        state = "read SEDs";
        sed.resize(ntime, nlam);
        for (uint_t it : range(ntime)) {
            uint_t nstep = 0;
            in >> nstep;

            if (nstep != nlam) {
                error("corrupted file, wavelength step mismatch: ", nstep, " vs. ", nlam);
                error("reading time step ", it, " of ", ntime);
                return 1;
            }

            for (uint_t il : range(lambda)) {
                in >> sed.safe(it,il);
            }

            // Read extra information
            uint_t nfunc = 0;
            in >> nfunc;

            for (uint_t i = 0; i < nfunc; ++i) {
                double f;
                in >> f;
            }
        }

        state = "read extras";
        uint_t nextra = 12;
        for (uint_t ie : range(nextra)) {
            uint_t nstep = 0;
            in >> nstep;

            vec1d extra(nstep);
            for (uint_t i : range(nstep)) {
                in >> extra[i];
            }

            if (ie == 1) {
                mass = extra/totm;
            }
        }
    } catch (...) {
        print("");
        error("could not read data in library file '", filename, "'");
        error("could not ", state);
        error("the file is probably corrupted, try re-downloading it");
        print("");
        return false;
    }

    // Write FITS file for faster reading next time
    fits::write_table(filename+".fits", ftable(age, mass, lambda, sed));

    return true;
}

bool ssp_bc03::read_fits(std::string filename, bool noflux) {
    fits::input_table itbl(filename);
    if (!itbl.read_column("age", age)) {
        print("");
        error("could not read column 'AGE' (time) from FITS file '", filename, "'");
        print("");
        return false;
    }
    if (!itbl.read_column("mass", mass)) {
        print("");
        error("could not read column 'MASS' from FITS file '", filename, "'");
        print("");
        return false;
    }
    if (!noflux) {
        if (!itbl.read_column("lambda", lambda)) {
            print("");
            error("could not read column 'LAMBDA' from FITS file '", filename, "'");
            print("");
            return false;
        }
        if (!itbl.read_column("sed", sed)) {
            print("");
            error("could not read column 'SED' from FITS file '", filename, "'");
            print("");
            return false;
        }
    }
    return true;
}

bool ssp_bc03::read(std::string filename, bool noflux) {
    if (file::exists(filename+".ised_ASCII.fits")) {
        return read_fits(filename+".ised_ASCII.fits", noflux);
    } else if (file::exists(filename+".ised_ASCII")) {
        return read_ascii(filename+".ised_ASCII");
    } else {
        print("");
        error("could not find library: '", filename, "'");
        error("expected extensions *.fits or *.ised_ASCII");
        print("");
        return false;
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
        sfh_expr.vars[0] = (opts.custom_sfh_lookback ? nage - t.safe[i] : t.safe[i]);
        sfh.safe[i] = sfh_expr.eval();
    }
}

template<typename F>
void integrate_ssp(const vec1d& age, const vec1d& sfr, const vec1d& ssp_age, F&& func) {
    double t2 = 0.0;
    uint_t ihint = npos;
    for (uint_t it : range(ssp_age)) {
        double t1 = t2;
        if (it < ssp_age.size()-1) {
            t2 = 0.5*(ssp_age.safe[it] + ssp_age.safe[it+1]);
        } else {
            t2 = ssp_age.back();
        }

        if (t2 <= age.front()) {
            continue;
        }

        t1 = max(t1, age.front());
        t2 = min(t2, age.back());

        func(it, integrate_hinted(age, sfr, ihint, t1, t2));

        if (t2 >= age.back()) {
            break;
        }
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
    const vec1f& output_av = output.grid[grid_id::av];

    // Compute "cosmic" time (t=0 is when the galaxy is born)
    const double dt = opts.custom_sfh_step;
    const vec1d ctime = reverse(dt*dindgen(uint_t(ceil(e10(max(output_age))/dt)+1.0)));
    // NB: age array is sorted from largest to smallest

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
        if (is_finite(opts.apply_vdisp)) {
            ssp.sed = convolve_vdisp(ssp.lambda, ssp.sed, opts.apply_vdisp);
        }

        // Pre-compute dust law & IGM absorption (they don't change with SFH)
        vec2d dust_law = build_dust_law(output_av, ssp.lambda);
        vec2d igm_abs = build_igm_absorption(output_z, ssp.lambda);

        // Function to build a model
        auto do_model = [&](model_id_pair& tm) {
            float& model_mass  = tm.model.props[prop_id::mass];
            float& model_mform = tm.model.props[prop_id::mform];
            float& model_sfr   = tm.model.props[prop_id::sfr];
            float& model_ssfr  = tm.model.props[prop_id::ssfr];

            // Compute lookback time (t=0 is when the galaxy is observed, t>0 is in the past)
            uint_t ia = tm.idm[grid_id::age];
            vec1d ltime = e10(output_age[ia]) - ctime;

            // Build analytic SFH
            vec1d sfh; {
                auto lock = (opts.parallel == parallel_choice::generators ?
                    std::unique_lock<std::mutex>(sfh_mutex) : std::unique_lock<std::mutex>());

                evaluate_sfh_custom(tm.idm, ctime, sfh);
            }

            // Integrate SFH on local time grid
            vec1d tpl_flux(ssp.lambda.size());
            double tmodel_mass  = 0.0;
            double tformed_mass = 0.0;
            integrate_ssp(ltime, sfh, ssp.age, [&](uint_t it, double formed) {
                tmodel_mass  += formed*ssp.mass.safe[it];
                tformed_mass += formed*ssp.mass.safe[0];
                tpl_flux     += formed*ssp.sed.safe(it,_);
            });

            model_mass = tmodel_mass;
            model_mform = tformed_mass;

            if (opts.sfr_avg > 0) {
                // Average SFR over the past X yr
                double t1 = min(opts.sfr_avg, ltime.back());
                model_sfr = integrate(ltime, sfh, 0.0, t1)/opts.sfr_avg;
            } else {
                // Use instantaneous SFR
                model_sfr = interpolate(sfh, ltime, 0.0);
            }

            model_ssfr = model_sfr/model_mass;

            // The rest is not specific to the SFH, use generic code
            build_and_send_impl(fitter, pg, ssp.lambda, tpl_flux, dust_law, igm_abs,
                output_age[ia], tm.idm, tm.model);
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

    ssp_bc03* ssp = cached_ssp_bc03.get();

    {
        auto lock = (opts.n_thread > 1 ?
            std::unique_lock<std::mutex>(sed_mutex) : std::unique_lock<std::mutex>());

        if (!cached_ssp_bc03) {
            cached_ssp_bc03 = std::unique_ptr<ssp_bc03>(new ssp_bc03());
            ssp = cached_ssp_bc03.get();
        }

        // Load SSP
        std::string filename = get_library_file_ssp(im);
        if (filename != cached_library) {
            if (!ssp->read(filename)) {
                return false;
            }

            cached_library = filename;

            // Apply velocity dispersion
            if (is_finite(opts.apply_vdisp)) {
                ssp->sed = convolve_vdisp(ssp->lambda, ssp->sed, opts.apply_vdisp);
            }
        }
    }

    // Build analytic SFH
    double dt = opts.custom_sfh_step;
    vec1d ctime = reverse(dt*dindgen(uint_t(ceil(e10(max(output_age))/dt)+1.0)));
    // NB: age array is sorted from largest to smallest
    vec1d sfh; {
        auto lock = (opts.n_thread > 1 ?
            std::unique_lock<std::mutex>(sfh_mutex) : std::unique_lock<std::mutex>());

        evaluate_sfh_custom(idm, ctime, sfh);
    }

    // Integrate SFH on local time grid
    vec1d tpl_flux(ssp->lambda.size());
    integrate_ssp(e10(output_age[ia]) - ctime, sfh, ssp->age,
        [&](uint_t it, double formed) {
            tpl_flux += formed*ssp->sed.safe(it,_);
        }
    );

    lam = ssp->lambda;
    flux = tpl_flux;

    return true;
}

bool gridder_t::get_sfh_custom(uint_t iflat, const vec1d& t, vec1d& sfh,
    const std::string& type) const {

    vec1u idm = grid_ids(iflat);
    double nage = e10(output.grid[grid_id::age][idm[grid_id::age]]);
    double age_obs = e10(auniv[grid_id::z]);
    double age_born = age_obs - nage;

    // Evaluate SFH
    uint_t i0 = upper_bound(t, age_born);
    uint_t i1 = upper_bound(t, age_obs);

    if (i1 == npos) {
        i1 = t.size()-1;
    }
    if (i0 == npos) {
        i0 = 0;
    }

    sfh = replicate(0.0, t.size()); {
        vec1d tsfh;
        evaluate_sfh_custom(idm, t[i0-_] - age_born, tsfh);
        sfh[i0-_] = tsfh;
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
        integrate_ssp(nage + age_born - reverse(t), reverse(sfh), ssp.age,
            [&](uint_t it, double formed) {
                mass += formed*ssp.mass.safe[it];
            }
        );

        sfh /= mass;
    } else if (type == "mass") {
        // Integrate mass, including mass loss
        vec1d mass(t.size());
        vec1d lsfh = reverse(sfh);

        for (uint_t i : range(t)) {
            integrate_ssp(t.safe[i] - reverse(t), lsfh, ssp.age,
                [&](uint_t it, double formed) {
                    mass.safe[i] += formed*ssp.mass.safe[it];
                }
            );
        }

        // Normalize to unit mass at observation
        mass /= interpolate(mass, t, nage + age_born);

        // Return mass instead of SFR
        std::swap(mass, sfh);
    } else {
        error("unknown SFH type '", type, "'");
        return false;
    }

    return true;
}
