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

        bool error_shown = false;

        for (uint_t is : range(input.id)) {
            if (is_finite(input.zphot.safe(is,ilow)) && is_finite(input.zphot.safe(is,iup))) {
                // Get range from zlow and zup
                idzl.safe[is] = min_id(abs(output_z - input.zphot.safe(is,ilow)));
                idzu.safe[is] = min_id(abs(output_z - input.zphot.safe(is,iup)));

                // Check that zphot falls inside confidence interval
                if (opts.best_at_zphot && idzp.safe[is] != npos && idz.safe[is] == npos &&
                    (idzp.safe[is] < idzl.safe[is] || idzp.safe[is] > idzu.safe[is])) {
                    warning("for galaxy ", input.id.safe[is], " the photo-z (",
                        input.zphot.safe(is,0), ") falls outside of the chosen confidence "
                        "interval (", input.zphot.safe(is,ilow), " to ", input.zphot.safe(is,iup),
                        ")");
                    error_shown = true;
                }
            }
        }

        if (error_shown) {
            warning("for these galaxies, the confidence intervals may not be consistent with "
                "the best fit values");
        }
    }

    for (uint_t is : range(input.id)) {
        if (idz.safe[is] != npos) {
            // Range is limited to z_spec
            idzl.safe[is] = idz.safe[is];
            idzu.safe[is] = idz.safe[is];
        }
    }

    has_spec = replicate(false, input.id.size());
    for (uint_t is : range(input.id)) {
        for (uint_t il : range(input.spec_start, input.spec_end)) {
            if (is_finite(input.eflux.safe(is,il))) {
                has_spec.safe[is] = true;
                break;
            }
        }
    }

    // Pre-generate random fluctuations
    if (opts.n_sim > 0) {
        auto seed = make_seed(42);
        sim_rnd = randomn(seed, opts.n_sim, input.lambda.size()+1);
    }

    // Initialize chi2 grid if asked
    save_chi2 = opts.save_chi_grid || opts.save_bestchi > 0;

    if (opts.save_chi_grid) {
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
        ochi2.out_file.open(ochi2.out_filename, std::ios::binary | std::ios::out | std::ios::trunc);

        // Write header
        // Format:
        // char:   file type ('C' = chi2 grid)
        // uint32: size of header in bytes (to skip it)
        // uint32: number of galaxies
        // uint32: number of properties
        // for each property:
        //     char[*]: name
        // uint32: number of grid axis
        // for each grid axis:
        //     char[*]: name
        //     uint32: number of values
        //     float[*]: grid values
        const unsigned char ftype_chi2grid = 'C';
        file::write_as<std::uint32_t>(ochi2.out_file, 0);
        file::write(ochi2.out_file, ftype_chi2grid);
        file::write_as<std::uint32_t>(ochi2.out_file, input.id.size());

        file::write_as<std::uint32_t>(ochi2.out_file, gridder.nprop);
        for (uint_t i : range(gridder.nprop)) {
            file::write(ochi2.out_file, output.param_names[gridder.nparam+i]);
        }

        file::write_as<std::uint32_t>(ochi2.out_file, gridder.grid_dims.size());
        for (uint_t i : range(gridder.grid_dims)) {
            file::write(ochi2.out_file, output.param_names[i]);
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

    if (opts.save_bestchi > 0) {
        if (save_chi2) {
            std::string odir = opts.output_dir+"best_chi2/";
            file::mkdir(odir);

            chi2_filename.resize(input.id.size());
            best_chi2 = replicate(finf, input.id.size());

            std::ofstream out;
            for (uint_t is : range(input.id)) {
                chi2_filename[is] = odir+opts.catalog+"_"+input.id[is]+".chi2.grid";
                out.open(chi2_filename[is], std::ios::binary | std::ios::out | std::ios::trunc);

                if (is == 0) {
                    // Format:
                    // char:   file type ('B' = best chi2)
                    // uint32: size of header in bytes (to skip it)
                    // uint32: number of properties
                    // for each property:
                    //     char[*]: name
                    // uint32: number of grid axis
                    // for each grid axis:
                    //     char[*]: name
                    //     uint32: number of values
                    //     float[*]: grid values
                    const unsigned char ftype_bestchi2 = 'B';
                    file::write_as<std::uint32_t>(out, 0);
                    file::write(out, ftype_bestchi2);
                    file::write_as<std::uint32_t>(out, gridder.nprop);
                    for (uint_t i : range(gridder.nprop)) {
                        file::write(out, output.param_names[gridder.nparam+i]);
                    }

                    file::write_as<std::uint32_t>(out, gridder.grid_dims.size());
                    for (uint_t i : range(gridder.grid_dims)) {
                        file::write(out, output.param_names[i]);
                        file::write_as<std::uint32_t>(out, gridder.grid_dims[i]);
                        file::write(out, output.grid[i]);
                    }

                    obchi2.hpos = out.tellp();

                    out.seekp(0);
                    file::write_as<std::uint32_t>(out, obchi2.hpos);
                    out.close();

                    // Read back header for faster writes
                    std::ifstream in(chi2_filename[is], std::ios::binary | std::ios::in);
                    obchi2.header.resize(obchi2.hpos);
                    file::read(in, obchi2.header);
                } else {
                    file::write(out, obchi2.header);
                    out.close();
                }
            }
        } else {
            warning("best chi2 file will also not be saved because of lack of disk space");
        }
    }

    // Start threads if using parallel execution
    if (opts.parallel == parallel_choice::models) {
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

    if (opts.save_bestchi > 0) {
        std::ofstream out;
        std::ifstream in;

        for (uint_t cis : range(chi2)) {
            uint_t is = cis + i0;

            if (chi2.safe[cis] > best_chi2.safe[is] + opts.save_bestchi) continue;

            if (chi2.safe[cis] < best_chi2.safe[is]) {
                best_chi2.safe[is] = chi2.safe[cis];

                struct datum {
                    uint32_t id;
                    float chi2;
                    vec1f p;
                };

                // Read the saved data and write simultaneously
                file::move(chi2_filename.safe[is], chi2_filename.safe[is]+".old");
                in.open(chi2_filename.safe[is]+".old", std::ios::binary | std::ios::in);
                in.seekg(obchi2.hpos);

                out.open(chi2_filename.safe[is], std::ios::binary | std::ios::out);
                file::write(out, obchi2.header);

                while (in) {
                    uint32_t id;
                    float chi2;
                    vec1f p(gridder.nprop);

                    if (file::read(in, id) && file::read(in, chi2) && file::read(in, p)) {
                        if (chi2 > best_chi2.safe[is] + opts.save_bestchi) continue;
                        file::write(out, id);
                        file::write(out, chi2);
                        file::write(out, p);
                    }
                }

                in.close();
                file::remove(chi2_filename.safe[is]+".old");

                file::write_as<std::uint32_t>(out, igrid);
                file::write(out, chi2.safe[cis]);
                file::write(out, vec1f(props.safe(cis,_)));
                out.close();
            } else {
                out.open(chi2_filename.safe[is], std::ios::binary | std::ios::out | std::ios::app);
                file::write_as<std::uint32_t>(out, igrid);
                file::write(out, chi2.safe[cis]);
                file::write(out, vec1f(props.safe(cis,_)));
                out.close();
            }
        }
    }

    if (opts.save_chi_grid) {
        auto p0 = ochi2.hpos + igrid*input.id.size()*(1+gridder.nprop)*sizeof(float);

        ochi2.out_file.seekp(p0 + i0*sizeof(float));
        file::write(ochi2.out_file, chi2);

        for (uint_t p : range(gridder.nprop)) {
            ochi2.out_file.seekp(p0 + ((1+p)*input.id.size() + i0)*sizeof(float));
            file::write(ochi2.out_file, vec1f(props.safe(_,p)));
        }
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

        chi2 = replicate(finf, ngal);
        props = replicate(fnan, ngal, nprop);
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

    if (opts.debug) {
        // DEBUG: check that model has all finite values
        if (count(!is_finite(model.flux)) > 0) {
            fits::write_table("debug.fits",
                "flux", model.flux, "grid_id", gridder.grid_ids(model.igrid));
            phypp_check(false, "model has invalid values, saved in debug.fits for inspection");
        }
        // DEBUG: check that model is not only zero values
        if (count(model.flux > 0.0) == 0) {
            fits::write_table("debug.fits",
                "flux", model.flux, "grid_id", gridder.grid_ids(model.igrid));
            phypp_check(false, "model has all zero values, saved in debug.fits for inspection");
        }
    }

    for (uint_t i : range(i1-i0)) {
        uint_t is = i + i0;

        // Apply constraints on redshift
        bool dofit = (idzl.safe[is] <= iz && iz <= idzu.safe[is]);
        bool keepfit = true;
        if (opts.best_at_zphot && idzp.safe[is] != npos) {
            if (idzp.safe[is] == iz) {
                dofit = true;
                keepfit = true;
            } else {
                if (opts.n_sim == 0) {
                    dofit = false;
                }
                keepfit = false;
            }
        }

        if (!dofit) {
            // Skip this model
            continue;
        }

        // Compute weights and scaling factor
        double wfm = 0, wmm = 0;
        uint_t ndata = input.lambda.size();
        uint_t iflx = 0;

        uint_t nscale = ndata;
        bool spec_auto_scale = opts.auto_scale && has_spec.safe[is];
        if (spec_auto_scale) {
            nscale -= input.spec_end - input.spec_start;
        }

        // Add LIR as a data point in the fit
        if (!input.lir.empty() && is_finite(input.lir.safe[is])) {
            wsp.weight[0] = 1.0/input.lir_err.safe[is];
            wsp.wflux[0] = input.lir.safe[is]*wsp.weight[0];
            wsp.wmodel[0] = model.props.safe[prop_id::ldust]*wsp.weight[0];

            wfm += wsp.wmodel[0]*wsp.wflux[0];
            wmm += sqr(wsp.wmodel[0]);

            ++ndata;
            ++nscale;
            ++iflx;
        }

        for (uint_t il : range(iflx, nscale)) {
            wsp.weight[il] = (opts.temp_err_file.empty() ?
                1.0/input.eflux.safe(is,il-iflx) :
                1.0/sqrt((sqr(input.eflux.safe(is,il-iflx)) +
                    tpl_err.safe(iz,il-iflx)*sqr(input.flux.safe(is,il-iflx))))
            );

            wsp.wflux[il] = input.flux.safe(is,il-iflx)*wsp.weight[il];
            wsp.wmodel[il] = model.flux.safe[il-iflx]*wsp.weight[il];

            wfm += wsp.wmodel[il]*wsp.wflux[il];
            wmm += sqr(wsp.wmodel[il]);
        }

        double scale = wfm/wmm;
        double spec_scale = scale;

        // Auto scaling for spectral data
        // When enabled, spectral data doesn't participate in the global normalization,
        // just to the chi2. It is separately rescaled to match the model, which is equivalent
        // to assuming we don't know the absolute flux calibration but we are confident on
        // the shape of the spectrum.
        double swfm = 0, swmm = 0;
        if (spec_auto_scale) {
            for (uint_t il : range(nscale, ndata)) {
                wsp.weight[il] = 1.0/input.eflux.safe(is,il-iflx);
                wsp.wflux[il] = input.flux.safe(is,il-iflx)*wsp.weight[il];
                wsp.wmodel[il] = model.flux.safe[il-iflx]*wsp.weight[il];

                swfm += wsp.wmodel[il]*wsp.wflux[il];
                swmm += sqr(wsp.wmodel[il]);
            }

            spec_scale = swfm/swmm;

            if (!is_finite(spec_scale)) {
                // No spectral data
                spec_scale = scale;
            }
        }

        // Compute chi2
        double tchi2 = 0;
        for (uint_t il : range(nscale)) {
            tchi2 += sqr(wsp.wflux[il] - scale*wsp.wmodel[il]);
        }
        for (uint_t il : range(nscale, ndata)) {
            tchi2 += sqr(wsp.wflux[il] - spec_scale*wsp.wmodel[il]);
        }

        if (keepfit) {
            // Save chi2 and properties
            wsp.chi2.safe[i] = tchi2;
            for (uint_t ip : range(model.props)) {
                double tscale = (ip == prop_id::spec_scale ? scale/spec_scale : scale);
                wsp.props.safe(i,ip) = (output.param_scale.safe[gridder.nparam+ip] ?
                    tscale*model.props.safe[ip] : model.props.safe[ip]
                );
            }
        }

        if (opts.parallel == parallel_choice::none) {
            // Compare to best
            // WARNING: read/modify shared resource
            ++output.num_models[is];
            if (output.best_chi2.safe[is]  > wsp.chi2.safe[i]) {
                output.best_chi2.safe[is]  = wsp.chi2.safe[i];
                output.best_model.safe[is] = model.igrid;
                if (!opts.best_from_sim) {
                    for (uint_t ip : range(model.props)) {
                        output.best_params.safe(is,gridder.nparam+ip,0) = wsp.props.safe(i,ip);
                    }
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
                swfm = 0;
                for (uint_t il : range(nscale)) {
                    // In weighted units, the random perturbations have a sigma of unity
                    wsp.rflux[il] = wsp.wflux[il] + sim_rnd.safe(im,il);
                    wfm += wsp.wmodel[il]*wsp.rflux[il];
                }
                for (uint_t il : range(nscale, ndata)) {
                    // When auto scale is ON, spectral data doesn't participate in scale factor
                    wsp.rflux[il] = wsp.wflux[il] + sim_rnd.safe(im,il);
                    swfm += wsp.wmodel[il]*wsp.rflux[il];
                }

                scale = wfm/wmm;
                if (!spec_auto_scale) {
                    spec_scale = scale;
                } else {
                    spec_scale = swfm/swmm;
                }

                // Compute chi2
                tchi2 = 0;
                for (uint_t il : range(nscale)) {
                    tchi2 += sqr(wsp.rflux[il] - scale*wsp.wmodel[il]);
                }
                for (uint_t il : range(nscale, ndata)) {
                    tchi2 += sqr(wsp.rflux[il] - spec_scale*wsp.wmodel[il]);
                }

                if (opts.parallel == parallel_choice::none) {
                    // Compare to best
                    // WARNING: read/modify shared resource
                    if (output.mc_best_chi2.safe(is,im)  > tchi2) {
                        output.mc_best_chi2.safe(is,im)  = tchi2;
                        output.mc_best_model.safe(is,im) = model.igrid;
                        for (uint_t ip : range(model.props)) {
                            double tscale = (ip == prop_id::spec_scale ? scale/spec_scale : scale);
                            output.mc_best_props.safe(is,ip,im) = (
                                output.param_scale.safe[gridder.nparam+ip] ?
                                tscale*model.props.safe[ip] : model.props.safe[ip]
                            );
                        }
                    }
                } else {
                    // If multithreaded, accumulate simulated values and commit
                    // them to the shared array in one batch
                    wsp.mc_chi2.safe[im] = tchi2;
                    for (uint_t ip : range(model.props)) {
                        double tscale = (ip == prop_id::spec_scale ? scale/spec_scale : scale);
                        wsp.mc_props.safe(ip,im) = (output.param_scale.safe[gridder.nparam+ip] ?
                            tscale*model.props.safe[ip] : model.props.safe[ip]
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

    if (opts.parallel != parallel_choice::none) {
        // Compare to best
        // WARNING: read/modify shared resource
        std::lock_guard<std::mutex> lock(output.fit_result_mutex);

        for (uint_t i : range(i1-i0)) {
            uint_t is = i + i0;
            ++output.num_models[is];
            if (output.best_chi2.safe[is]  > wsp.chi2[i]) {
                output.best_chi2.safe[is]  = wsp.chi2[i];
                output.best_model.safe[is] = model.igrid;
                if (!opts.best_from_sim) {
                    for (uint_t ip : range(model.props)) {
                        output.best_params.safe(is,gridder.nparam+ip,0) = wsp.props.safe(i,ip);
                    }
                }
            }
        }
    }
}

void fitter_t::fit(const model_t& model) {
    if (opts.parallel == parallel_choice::sources) {
        workers_multi_source->process(model);
    } else if (opts.parallel == parallel_choice::models) {
        workers_multi_model->process(model);
    } else {
        fit_galaxies(model, 0, input.id.size());
    }
}

vec2d make_grid_bins(const vec1d& grid) {
    vec2d bins(2, grid.size());

    if (grid.size() == 1) {
        // Single bin, just make it large enough to account for potential numerical jitter
        bins(0,1) = grid.safe[0]*0.9 - 0.1;
        bins(1,1) = grid.safe[0]*1.1 + 0.1;
    } else {
        // Build bins surrounding each grid value
        for (uint_t k : range(grid)) {
            if (k == 0) {
                bins(0,k) = grid.safe[k] - 0.5*(grid.safe[k+1] - grid.safe[k]);
            } else {
                bins(0,k) = bins(1,k-1);
            }

            if (k == grid.size()-1) {
                bins(1,k) = grid.safe[k] + 0.5*(grid.safe[k] - grid.safe[k-1]);
            } else {
                bins(1,k) = 0.5*(grid.safe[k] + grid.safe[k+1]);
            }
        }
    }

    return bins;
}

void fitter_t::find_best_fits() {
    if (opts.parallel == parallel_choice::models) {
        if (opts.verbose) note("waiting for all models to finish...");
        workers_multi_model->workers.join();
    } else if (opts.parallel == parallel_choice::sources) {
        if (opts.verbose) note("waiting for all models to finish...");
        workers_multi_source->workers.join();
    }

    bool silence_invalid_chi2 = false;

    if (opts.verbose) note("finding best fits...");
    for (uint_t is : range(input.id)) {
        if (!silence_invalid_chi2 && !is_finite(output.best_chi2[is])) {
            warning("galaxy ", input.id[is], " has no best fit solution");
            warning("there is probably a problem with the models, please re-run FAST++ with DEBUG=1");
            warning("(further occurences of this warning for other galaxies will be suppressed)");
            silence_invalid_chi2 = true;
        }

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

            // For grid parameters, use cumulative distribution
            for (uint_t ip : range(gridder.nparam)) {
                auto& grid = output.grid[ip];

                if (grid.size() == 1) {
                    if (opts.best_from_sim) {
                        output.best_params(is,ip,0) = grid[0];
                    }

                    for (uint_t ic : range(input.conf_interval)) {
                        output.best_params(is,ip,1+ic) = grid[0];
                    }
                } else {
                    vec2d bins = make_grid_bins(grid);
                    vec1d cnt = cumul(histogram(bparams.safe(ip,_), bins));
                    cnt /= cnt.back();

                    if (opts.best_from_sim) {
                        if (cnt[0] < 0.5) {
                            output.best_params(is,ip,0) = interpolate(grid, cnt, 0.5);
                        } else {
                            output.best_params(is,ip,0) = grid[0];
                        }
                    }

                    for (uint_t ic : range(input.conf_interval)) {
                        double c = input.conf_interval[ic];
                        if (cnt[0] < c) {
                            output.best_params(is,ip,1+ic) = interpolate(grid, cnt, c);
                        } else {
                            output.best_params(is,ip,1+ic) = grid[0];
                        }
                    }
                }
            }

            // For properties, use percentiles
            for (uint_t ip : range(gridder.nparam, bparams.dims[0])) {
                vec1d bp = bparams.safe(ip,_);

                if (opts.best_from_sim) {
                    output.best_params(is,ip,0) = inplace_median(bp);
                }

                for (uint_t ic : range(input.conf_interval)) {
                    output.best_params(is,ip,1+ic) = inplace_percentile(bp, input.conf_interval[ic]);
                }
            }
        }
    }
}
