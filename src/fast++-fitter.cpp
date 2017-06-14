#include "fast++.hpp"
#include <phypp/utility/thread.hpp>

fitter_t::fitter_t(const options_t& opt, const input_state_t& inp, const gridder_t& gri,
    output_state_t& out) : opts(opt), input(inp), gridder(gri), output(out) {

    vec1f& output_z = output.grid[grid_id::z];

    // Pre-compute template error
    if (!opts.temp_err_file.empty()) {
        if (opts.verbose) note("initializing template error function...");
        double l0 = median(input.tplerr_lam);

        tpl_err.resize(output_z.size(), input.lambda.size());
        for (uint_t iz : range(output_z))
        for (uint_t il : range(input.lambda)) {
            tpl_err.safe(iz,il) = sqr(astro::sed2flux(
                input.filters[il].wl, input.filters[il].tr,
                input.tplerr_lam*(1.0 + output_z[iz]), input.tplerr_err
            ));

            if (!is_finite(tpl_err.safe(iz,il))) {
                // The filter goes out of the template error function, extrapolate
                if (input.lambda.safe[il] < l0) {
                    tpl_err.safe(iz,il) = sqr(input.tplerr_err.front());
                } else {
                    tpl_err.safe(iz,il) = sqr(input.tplerr_err.back());
                }
            }
        }
    }

    // Pre-identify zgrid for galaxies with zspec or zphot
    if (opts.verbose) note("identify zspecs and redshift bounds in the grid...");
    idz = replicate(npos, input.id.size());
    idzp = replicate(npos, input.id.size());
    for (uint_t is : range(input.id)) {
        uint_t izp = npos;
        if (!input.zphot.empty()) {
            // Use zphot if we have one
            if (is_finite(input.zphot.safe(is,0))) {
                izp = min_id(abs(output_z - input.zphot.safe(is,0)));
            }
        }

        idzp.safe[is] = izp;
        if (opts.force_zphot) {
            idz.safe[is] = izp;
        }

        // Override with zspec if we have one
        if (is_finite(input.zspec.safe[is])) {
            idzp.safe[is] = idz.safe[is] = min_id(abs(output_z - input.zspec.safe[is]));
        }
    }

    // Pre-identify zgrid boundaries for galaxies with zphot confidence interval
    idzl = replicate(0, input.id.size());
    idzu = replicate(output_z.size()-1, input.id.size());

    if (is_finite(opts.zphot_conf) && !input.zphot.empty()) {
        // Use chosen confidence interval from EAzY
        uint_t ic = min_id(abs(input.zphot_conf - opts.zphot_conf));
        uint_t ilow = 1 + 2*ic + 0;
        uint_t iup  = 1 + 2*ic + 1;

        for (uint_t is : range(input.id)) {
            // Get range from zlow and zup
            if (is_finite(input.zphot.safe(is,ilow)) && is_finite(input.zphot.safe(is,iup))) {
                idzl.safe[is] = min_id(abs(output_z - input.zphot.safe(is,ilow)));
                idzu.safe[is] = min_id(abs(output_z - input.zphot.safe(is,iup)));
            }
        }
    }

    // Pre-generate random fluctuations
    if (opts.n_sim > 0) {
        auto seed = make_seed(42);
        sim_rnd = randomn(seed, opts.n_sim, input.lambda.size()+1);
    }

    // Initialize chi2 grid if asked
    save_chi2 = opts.save_chi_grid;
    if (save_chi2) {
        // Create chi2 grid on disk
        if (opts.verbose) {
            double expsize = input.id.size()*double(gridder.nmodel)*(1+gridder.nprop)*sizeof(float);

            std::string unit = "B";
            vec1s units = vec1s{"k", "M", "G", "T", "P"}+"B";
            for (uint_t i : range(units)) {
                if (expsize > 1024) {
                    expsize /= 1024;
                    unit = units[i];
                } else {
                    break;
                }
            }
            note("initializing chi2 grid on disk... (expected size ", expsize, " ", unit, ")");
        }

        ochi2.out_filename = opts.output_dir+"chi2.grid";
        ochi2.out_file.open(ochi2.out_filename, std::ios::binary | std::ios::out);

        // Write header
        // Format:
        // uint32: size of header in bytes (to skip it)
        // uint32: number of galaxies
        // uint32: number of grid axis
        // for each grid axis:
        //     uint32: number of values
        //     float[*]: grid values
        file::write_as<std::uint32_t>(ochi2.out_file, 0);
        file::write_as<std::uint32_t>(ochi2.out_file, input.id.size());
        file::write_as<std::uint32_t>(ochi2.out_file, gridder.grid_dims.size());
        for (uint_t i : range(gridder.grid_dims)) {
            file::write_as<std::uint32_t>(ochi2.out_file, gridder.grid_dims[i]);
            file::write(ochi2.out_file, output.grid[i]);
        }

        ochi2.hpos = ochi2.out_file.tellp();

        ochi2.out_file.seekp(0);
        file::write_as<std::uint32_t>(ochi2.out_file, ochi2.hpos);
        ochi2.out_file.seekp(ochi2.hpos);

        // Populate the file with empty data now
        // For each point of the grid, we store chi2 and properties
        uint_t nmodel1 = gridder.nmodel;
        uint_t idm = max_id(gridder.grid_dims);
        nmodel1 /= gridder.grid_dims[idm];
        vec1f chunk = replicate(fnan, gridder.grid_dims[idm]*input.id.size()*(1+gridder.nprop));
        for (uint_t im = 0; im < nmodel1; ++im) {
            if (!file::write(ochi2.out_file, chunk)) {
                warning("the chi2 grid could not be initialized");
                ochi2.out_file.close();
                file::remove(ochi2.out_filename);

                warning("in case you ran out of disk space, the file was deleted "
                    "and the grid will not be saved");
                save_chi2 = false;
                break;
            }
        }
    }

    // Start threads if using parallel execution
    if (opts.parallel == parallel_choice::none) {
        // nothing to do
    } else if (opts.parallel == parallel_choice::models) {
        workers_multi_model = std::unique_ptr<workers_multi_model_t>(
            new workers_multi_model_t(*this)
        );
    } else if (opts.parallel == parallel_choice::sources) {
        workers_multi_source = std::unique_ptr<workers_multi_source_t>(
            new workers_multi_source_t(*this)
        );
    }
}

fitter_t::workers_multi_source_t::workers_multi_source_t(fitter_t& f) : fitter(f) {
    workers.start(fitter.opts.n_thread, [this](const model_source_pair& p) {
        fitter.fit_galaxies(p.model, p.i0, p.i1);
    });
}

void fitter_t::workers_multi_source_t::process(const model_t& model) {
    while (fitter.opts.max_queued_fits > 0 && workers.remaining() > fitter.opts.max_queued_fits) {
        thread::sleep_for(1e-6);
    }

    uint_t i0 = 0;
    uint_t dn = fitter.input.id.size()/fitter.opts.n_thread;
    for (uint_t iw : range(fitter.opts.n_thread)) {
        uint_t i1 = (iw == fitter.opts.n_thread-1 ? fitter.input.id.size() : i0 + dn);
        workers.process(iw, model_source_pair(model, i0, i1));
        i0 = i1;
    }
}

fitter_t::workers_multi_model_t::workers_multi_model_t(fitter_t& f) : fitter(f) {
    workers.start(fitter.opts.n_thread, [this](const model_t& model) {
        fitter.fit_galaxies(model, 0, fitter.input.id.size());
    });
}

void fitter_t::workers_multi_model_t::process(const model_t& model) {
    while (fitter.opts.max_queued_fits > 0 && workers.remaining() > fitter.opts.max_queued_fits) {
        thread::sleep_for(1e-6);
    }

    workers.process(model);
}

void fitter_t::write_chi2(uint_t igrid, const vec1f& chi2, const vec2f& props, uint_t i0) {
    // TODO: consider putting this in a worker thread if this is slowing down too much
    auto p0 = ochi2.hpos + igrid*input.id.size()*(1+gridder.nprop)*sizeof(float);

    ochi2.out_file.seekp(p0 + i0*sizeof(float));
    file::write(ochi2.out_file, chi2);

    for (uint_t p : range(gridder.nprop)) {
        ochi2.out_file.seekp(p0 + (p*input.id.size() + i0)*sizeof(float));
        file::write(ochi2.out_file, props(_,p));
    }
}

struct fitter_workspace {
    // Loop local, per source
    char*   pool     = nullptr;
    double* weight   = nullptr;
    double* wflux    = nullptr;
    double* wmodel   = nullptr;
    double* rflux    = nullptr;
    vec1f  mc_chi2;
    vec2f  mc_props;

    // For all sources
    vec1f chi2;
    vec2f props;

    fitter_workspace(uint_t ngal, uint_t nflux, uint_t nprop, uint_t nsim) {
        uint_t buffer_size = sizeof(double)*nflux*(3 + (nsim > 0 ? 1 : 0));
        pool = new char[buffer_size];

        std::ptrdiff_t off = 0;
        weight = reinterpret_cast<double*>(pool + off); off += nflux*sizeof(double);
        wflux  = reinterpret_cast<double*>(pool + off); off += nflux*sizeof(double);
        wmodel = reinterpret_cast<double*>(pool + off); off += nflux*sizeof(double);

        if (nsim > 0) {
            rflux = reinterpret_cast<double*>(pool + off); off += nflux*sizeof(double);
        }

        phypp_check(uint_t(off) == buffer_size, "mismatch in buffer and arrays, please report!");

        mc_chi2.resize(nsim);
        mc_props.resize(nprop, nsim);

        chi2.resize(ngal);
        props.resize(ngal, nprop);
    }

    fitter_workspace(const fitter_workspace&) = delete;
    fitter_workspace(fitter_workspace&&) = delete;
    fitter_workspace& operator=(const fitter_workspace&) = delete;
    fitter_workspace& operator=(fitter_workspace&&) = delete;

    ~fitter_workspace() {
        delete pool;
    }
};

void fitter_t::fit_galaxies(const model_t& model, uint_t i0, uint_t i1) {
    fitter_workspace wsp(i1-i0, input.lambda.size()+1, model.props.size(), opts.n_sim);

    uint_t iz = gridder.grid_ids(model.igrid)[grid_id::z];

    for (uint_t i : range(i1-i0)) {
        uint_t is = i + i0;

        // Apply constraints on redshift
        if ((idz.safe[is] != npos && iz != idz.safe[is]) ||
            iz < idzl.safe[is] || iz > idzu.safe[is]) {
            wsp.chi2.safe[i]  = finf;
            wsp.props.safe(i,_) = fnan;
            continue;
        }

        // Compute weights and scaling factor
        double wfm = 0, wmm = 0;

        uint_t ndata = input.lambda.size();
        if (opts.temp_err_file.empty()) {
            for (uint_t il : range(input.lambda)) {
                wsp.weight[il] = 1.0/input.eflux.safe(is,il);

                wsp.wflux[il] = input.flux.safe(is,il)*wsp.weight[il];
                wsp.wmodel[il] = model.flux.safe[il]*wsp.weight[il];

                wfm += wsp.wmodel[il]*wsp.wflux[il];
                wmm += sqr(wsp.wmodel[il]);
            }
        } else {
            for (uint_t il : range(input.lambda)) {
                wsp.weight[il] = 1.0/sqrt((sqr(input.eflux.safe(is,il)) +
                    tpl_err.safe(iz,il)*sqr(input.flux.safe(is,il))));

                wsp.wflux[il] = input.flux.safe(is,il)*wsp.weight[il];
                wsp.wmodel[il] = model.flux.safe[il]*wsp.weight[il];

                wfm += wsp.wmodel[il]*wsp.wflux[il];
                wmm += sqr(wsp.wmodel[il]);
            }
        }

        if (!input.lir.empty() && is_finite(input.lir.safe[is])) {
            // Add LIR as a data point in the fit
            wsp.weight[ndata] = 1.0/input.lir_err.safe[is];
            wsp.wflux[ndata] = input.lir.safe[is]*wsp.weight[ndata];
            wsp.wmodel[ndata] = model.props.safe[prop_id::ldust]*wsp.weight[ndata];

            wfm += wsp.wmodel[ndata]*wsp.wflux[ndata];
            wmm += sqr(wsp.wmodel[ndata]);

            ++ndata;
        }

        double scale = wfm/wmm;

        // Compute chi2
        double tchi2 = 0;
        for (uint_t il : range(ndata)) {
            tchi2 += sqr(wsp.wflux[il] - scale*wsp.wmodel[il]);
        }

        wsp.chi2.safe[i] = tchi2;
        for (uint_t ip : range(model.props)) {
            wsp.props.safe(i,ip) = (output.param_scale.safe[gridder.nparam+ip] ?
                scale*model.props.safe[ip] : model.props.safe[ip]
            );
        }

        if (!opts.best_from_sim && opts.parallel == parallel_choice::none) {
            // Compare to best
            // WARNING: read/modify shared resource
            if (output.best_chi2.safe[is]  > wsp.chi2.safe[i] &&
                (!opts.best_at_zphot || idzp.safe[is] == npos || iz == idzp.safe[is])) {
                output.best_chi2.safe[is]  = wsp.chi2.safe[i];
                output.best_model.safe[is] = model.igrid;
                for (uint_t ip : range(model.props)) {
                    output.best_params.safe(is,gridder.nparam+ip,0) = wsp.props.safe(i,ip);
                }
            }
        }

        // Do MC simulation
        if (opts.n_sim > 0) {
            for (uint_t im : range(opts.n_sim)) {
                // Generate "random" fluxes and compute scaling factor
                // NB: since we create the randomness just once at the beginning of the fit
                // all models (and each galaxy) will use the same random numbers
                wfm = 0;
                for (uint_t il : range(ndata)) {
                    // In weighted units, the random perturbations have a sigma of unity
                    wsp.rflux[il] = wsp.wflux[il] + sim_rnd.safe(im,il);
                    wfm += wsp.wmodel[il]*wsp.rflux[il];
                }

                scale = wfm/wmm;

                // Compute chi2
                tchi2 = 0;
                for (uint_t il : range(ndata)) {
                    tchi2 += sqr(wsp.rflux[il] - scale*wsp.wmodel[il]);
                }

                if (opts.parallel == parallel_choice::none) {
                    // Compare to best
                    // WARNING: read/modify shared resource
                    if (output.mc_best_chi2.safe(is,im)  > tchi2) {
                        output.mc_best_chi2.safe(is,im)  = tchi2;
                        output.mc_best_model.safe(is,im) = model.igrid;
                        for (uint_t ip : range(model.props)) {
                            output.mc_best_props.safe(is,ip,im) = (
                                output.param_scale.safe[gridder.nparam+ip] ?
                                scale*model.props.safe[ip] : model.props.safe[ip]
                            );
                        }
                    }
                } else {
                    // If multithreaded, accumulate simulated values and commit
                    // them to the shared array in one batch
                    wsp.mc_chi2.safe[im] = tchi2;
                    for (uint_t ip : range(model.props)) {
                        wsp.mc_props.safe(ip,im) = (output.param_scale.safe[gridder.nparam+ip] ?
                            scale*model.props.safe[ip] : model.props.safe[ip]
                        );
                    }
                }
            }

            if (opts.parallel != parallel_choice::none) {
                // Compare to best
                // WARNING: read/modify shared resource
                std::lock_guard<std::mutex> lock(output.fit_result_mutex);

                for (uint_t im : range(opts.n_sim)) {
                    if (output.mc_best_chi2.safe(is,im)  > wsp.mc_chi2.safe[im]) {
                        output.mc_best_chi2.safe(is,im)  = wsp.mc_chi2.safe[im];
                        output.mc_best_model.safe(is,im) = model.igrid;
                        for (uint_t ip : range(model.props)) {
                            output.mc_best_props.safe(is,ip,im) = wsp.mc_props.safe(ip,im);
                        }
                    }
                }
            }
        }
    }

    if (save_chi2) {
        if (opts.parallel == parallel_choice::none) {
            write_chi2(model.igrid, wsp.chi2, wsp.props, i0);
        } else {
            // For thread safety
            std::lock_guard<std::mutex> lock(ochi2.write_mutex);
            write_chi2(model.igrid, wsp.chi2, wsp.props, i0);
        }
    }

    if (!opts.best_from_sim && opts.parallel != parallel_choice::none) {
        // Compare to best
        // WARNING: read/modify shared resource
        std::lock_guard<std::mutex> lock(output.fit_result_mutex);

        for (uint_t i : range(i1-i0)) {
            uint_t is = i + i0;
            if (output.best_chi2.safe[is]    > wsp.chi2[i] &&
                (!opts.best_at_zphot || idzp.safe[is] == npos || iz == idzp.safe[is])) {
                output.best_chi2.safe[is]    = wsp.chi2[i];
                output.best_model.safe[is]   = model.igrid;
                for (uint_t ip : range(model.props)) {
                    output.best_params.safe(is,gridder.nparam+ip,0) = wsp.props.safe(i,ip);
                }
            }
        }
    }
}

void fitter_t::fit(const model_t& model) {
    if (opts.parallel == parallel_choice::none) {
        fit_galaxies(model, 0, input.id.size());
    } else if (opts.parallel == parallel_choice::sources) {
        workers_multi_source->process(model);
    } else if (opts.parallel == parallel_choice::models) {
        workers_multi_model->process(model);
    }
}

void fitter_t::find_best_fits() {
    if (opts.parallel == parallel_choice::models) {
        if (opts.verbose) note("waiting for all models to finish...");
        workers_multi_model->workers.join();
    } else if (opts.parallel == parallel_choice::sources) {
        if (opts.verbose) note("waiting for all models to finish...");
        workers_multi_source->workers.join();
    }

    if (opts.verbose) note("finding best fits...");
    for (uint_t is : range(input.id)) {
        if (!opts.best_from_sim) {
            vec1u ids = gridder.grid_ids(output.best_model[is]);
            for (uint_t ip : range(gridder.nparam)) {
                output.best_params(is,ip,0) = output.grid[ip][ids[ip]];
            }
        }

        if (opts.n_sim > 0) {
            vec1u bmodel = output.mc_best_model(is,_);

            // Gather simulated parameters (+properties) for this source
            vec2f bparams(gridder.nparam+gridder.nprop, opts.n_sim);
            for (uint_t im : range(opts.n_sim)) {
                vec1u ids = gridder.grid_ids(bmodel[im]);
                for (uint_t ip : range(gridder.nparam)) {
                    bparams.safe(ip,im) = output.grid[ip][ids[ip]];
                }
            }

            for (uint_t ip : range(gridder.nprop)) {
                bparams.safe(gridder.nparam+ip,_) = output.mc_best_props.safe(is,ip,_);
            }

            if (opts.save_sim) {
                std::string odir = opts.output_dir+"best_fits/";
                if (!file::mkdir(odir)) {
                    warning("could not save simulations");
                    warning("the output directory '", odir, "' could not be created");
                } else {
                    fits::write_table(odir+opts.catalog+"_"+input.id[is]+".sims.fits",
                        "result", bparams, "params", output.param_names,
                        "chi2", output.mc_best_chi2(is,_)
                    );
                }
            }

            if (opts.best_from_sim) {
                for (uint_t ip : range(bparams.dims[0])) {
                    output.best_params(is,ip,0) = median(bparams.safe(ip,_));
                }
            }

            for (uint_t ip : range(bparams.dims[0])) {
                vec1d bp = bparams.safe(ip,_);
                for (uint_t ic : range(input.conf_interval)) {
                    output.best_params(is,ip,1+ic) = inplace_percentile(bp, input.conf_interval[ic]);
                }
            }
        }
    }
}
