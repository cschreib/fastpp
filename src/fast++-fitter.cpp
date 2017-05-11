#include "fast++.hpp"
#include <phypp/utility/thread.hpp>

fitter_t::fitter_t(const options_t& opt, const input_state_t& inp, const gridder_t& gri,
    output_state_t& out) : opts(opt), input(inp), gridder(gri), output(out) {

    // Pre-compute template error
    if (!opts.temp_err_file.empty()) {
        if (opts.verbose) note("initializing template error function...");
        double l0 = median(input.tplerr_lam);

        tpl_err.resize(output.z.size(), input.lambda.size());
        for (uint_t iz : range(output.z))
        for (uint_t il : range(input.lambda)) {
            tpl_err.safe(iz,il) = sqr(astro::sed2flux(
                input.filters[il].wl, input.filters[il].tr,
                input.tplerr_lam*(1.0 + output.z[iz]), input.tplerr_err
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
    for (uint_t is : range(input.id)) {
        if (opts.force_zphot && !input.zphot.empty()) {
            // Use zphot if we have one
            if (is_finite(input.zphot.safe(is,0))) {
                idz.safe[is] = min_id(abs(output.z - input.zphot.safe(is,0)));
            }
        }
        // Override with zspec if we have one
        if (is_finite(input.zspec.safe[is])) {
            idz.safe[is] = min_id(abs(output.z - input.zspec.safe[is]));
        }
    }

    // Pre-identify zgrid boundaries for galaxies with zphot confidence interval
    idzl = replicate(npos, input.id.size());
    idzu = replicate(npos, input.id.size());

    if (is_finite(opts.zphot_conf) && !input.zphot.empty()) {
        // Use chosen confidence interval from EAzY
        uint_t ic = min_id(abs(input.zphot_conf - opts.zphot_conf));
        uint_t ilow = 1 + 2*ic + 0;
        uint_t iup  = 1 + 2*ic + 1;

        for (uint_t is : range(input.id)) {
            // Get range from zlow and zup
            if (is_finite(input.zphot.safe(is,ilow)) && is_finite(input.zphot.safe(is,iup))) {
                idzl.safe[is] = min_id(abs(output.z - input.zphot.safe(is,ilow)));
                idzu.safe[is] = min_id(abs(output.z - input.zphot.safe(is,iup)));
            }
        }
    }

    // Pre-generate random fluctuations
    if (opts.n_sim > 0) {
        auto seed = make_seed(42);
        sim_rnd = randomn(seed, opts.n_sim, input.lambda.size());
    }

    // Initialize chi2 grid if asked
    save_chi2 = opts.save_chi_grid;
    if (save_chi2) {
        // Create chi2 grid on disk
        if (opts.verbose) note("initializing chi2 grid on disk...");
        ochi2.out_filename = opts.output_dir+"chi2.grid";
        ochi2.out_file.open(ochi2.out_filename);

        // Write header
        // Format:
        // uint32: size of header in bytes (to skip it)
        // uint32: number of grid axis (5)
        // uint32: number of metallicities
        // float[*]: metallicities
        // uint32: number of tau
        // float[*]: tau
        // uint32: number of ages
        // float[*]: ages
        // uint32: number of av
        // float[*]: av
        // uint32: number of z
        // float[*]: z
        ochi2.out_file << gridder.dims.size();
        ochi2.out_file << std::uint32_t(0);
        ochi2.out_file << std::uint32_t(gridder.dims.size());
        for (uint_t i : range(gridder.dims)) {
            ochi2.out_file << std::uint32_t(gridder.dims[i]);

            vec1f* gp = nullptr;
            switch (i) {
            case 0: gp = &output.metal; break;
            case 1: gp = &output.tau;   break;
            case 2: gp = &output.age;   break;
            case 3: gp = &output.av;    break;
            case 4: gp = &output.z;     break;
            }

            for (float v : *gp) {
                ochi2.out_file << v;
            }
        }

        ochi2.hpos = ochi2.out_file.tellp();
        ochi2.out_file.seekp(0);
        ochi2.out_file << std::uint32_t(ochi2.hpos);
        ochi2.out_file.seekp(ochi2.hpos);

        // Populate the file with empty data now
        // For each point of the grid, we store 3 values: chi2, mass and sfr
        uint_t nmodel1 = gridder.dims[0]*gridder.dims[1]*gridder.dims[2];
        vec1f empty = replicate(fnan, gridder.dims[3]*gridder.dims[4]*3);
        for (uint_t im = 0; im < nmodel1; ++im) {
            if (!ochi2.out_file.write(reinterpret_cast<const char*>(empty.data.data()),
                empty.size()*sizeof(float))) {
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
}

void fitter_t::chi2_output_manager_t::write_chi2(uint_t igrid, float chi2, float mass, float sfr) {
    out_file.seekp(hpos + igrid*3*sizeof(float));
    out_file << chi2 << mass << sfr;
}

void fitter_t::fit(model_t model) {
    vec1f chi2 = replicate(finf, input.id.size());
    vec1f mass = replicate(fnan, input.id.size());
    vec1f sfr = replicate(fnan, input.id.size());

    vec1d weight(input.lambda.size());
    vec1d wflux(input.lambda.size());
    vec1d wmodel(input.lambda.size());
    vec1d rflux;
    vec1f mc_chi2, mc_mass, mc_sfr;
    if (opts.n_sim > 0) {
        rflux.resize(input.lambda.size());
        mc_chi2 = mc_mass = mc_sfr = replicate(0.0f, opts.n_sim);
    }

    for (uint_t is : range(input.id)) {
        // Apply constraints on redshift
        if ((idz.safe[is] != npos && model.iz != idz.safe[is]) ||
            (idzl.safe[is] != npos && model.iz < idzl.safe[is]) ||
            (idzu.safe[is] != npos && model.iz > idzu.safe[is])) {
            continue;
        }

        // Compute weights and scaling factor
        double wfm = 0, wmm = 0;
        if (opts.temp_err_file.empty()) {
            for (uint_t il : range(input.lambda)) {
                weight.safe[il] = 1.0/input.eflux.safe(is,il);

                wflux.safe[il] = input.flux.safe(is,il)*weight.safe[il];
                wmodel.safe[il] = model.flux.safe[il]*weight.safe[il];

                wfm += wmodel.safe[il]*wflux.safe[il];
                wmm += sqr(wmodel.safe[il]);
            }
        } else {
            for (uint_t il : range(input.lambda)) {
                weight.safe[il] = 1.0/sqrt((sqr(input.eflux.safe(is,il)) +
                    tpl_err.safe(model.iz,il)*sqr(input.flux.safe(is,il))));

                wflux.safe[il] = input.flux.safe(is,il)*weight.safe[il];
                wmodel.safe[il] = model.flux.safe[il]*weight.safe[il];

                wfm += wmodel.safe[il]*wflux.safe[il];
                wmm += sqr(wmodel.safe[il]);
            }
        }

        double scale = wfm/wmm;

        // Compute chi2
        double tchi2 = 0;
        for (uint_t il : range(input.lambda)) {
            tchi2 += sqr(wflux.safe[il] - scale*wmodel.safe[il]);
        }

        if (save_chi2) {
            ochi2.write_chi2(model.igrid, tchi2, scale*model.mass, scale*model.sfr);
        }

        chi2.safe[is] = tchi2;
        mass.safe[is] = scale*model.mass;
        sfr.safe[is] = scale*model.sfr;

        // Do MC simulation
        if (opts.n_sim > 0) {
            for (uint_t im : range(opts.n_sim)) {
                // Generate "random" fluxes and compute scaling factor
                // NB: since we create the randomness just once at the beginning of the fit
                // all models (and each galaxy) will use the same random numbers
                double wfm = 0;
                for (uint_t il : range(input.lambda)) {
                    rflux.safe[il] = wflux.safe[il] + sim_rnd.safe(im,il);

                    wfm += wmodel.safe[il]*rflux.safe[il];
                }

                scale = wfm/wmm;

                // Compute chi2
                tchi2 = 0;
                for (uint_t il : range(input.lambda)) {
                    tchi2 += sqr(rflux.safe[il] - scale*wmodel.safe[il]);
                }

                mc_chi2.safe[im] = tchi2;
                mc_mass.safe[im] = scale*model.mass;
                mc_sfr.safe[im] = scale*model.sfr;
            }

            {
                // Compare to best
                // WARNING: read/modify shared resource
                std::lock_guard<std::mutex> lock(output.fit_result_mutex);

                for (uint_t im : range(opts.n_sim)) {
                    if (output.mc_best_chi2.safe(is,im)  > mc_chi2.safe[im]) {
                        output.mc_best_chi2.safe(is,im)  = mc_chi2.safe[im];
                        output.mc_best_mass.safe(is,im)  = mc_mass.safe[im];
                        output.mc_best_sfr.safe(is,im)   = mc_sfr.safe[im];
                        output.mc_best_model.safe(is,im) = model.igrid;
                    }
                }
            }
        }
    }

    {
        // Compare to best
        // WARNING: read/modify shared resource
        std::lock_guard<std::mutex> lock(output.fit_result_mutex);

        for (uint_t is : range(input.id)) {
            if (output.best_chi2.safe[is]   > chi2.safe[is]) {
                output.best_chi2.safe[is]   = chi2.safe[is];
                output.best_mass.safe(is,0) = mass.safe[is];
                output.best_sfr.safe(is,0)  = sfr.safe[is];
                output.best_model.safe[is]  = model.igrid;
            }
        }
    }
}

void fitter_t::find_best_fits() {
    for (uint_t is : range(input.id)) {
        vec1u ids = gridder.grid_ids(output.best_model[is]);
        output.best_ssfr(is,0)  = output.best_sfr(is,0)/output.best_mass(is,0);
        output.best_metal(is,0) = output.metal[ids[0]];
        output.best_tau(is,0)   = output.tau[ids[1]];
        output.best_age(is,0)   = output.age[ids[2]];
        output.best_av(is,0)    = output.av[ids[3]];
        output.best_z(is,0)     = output.z[ids[4]];

        if (opts.n_sim > 0) {
            vec1f bmass, bsfr, bssfr, bmetal, btau, bage, bav, bz;

            bmass = output.mc_best_mass(is,_);
            bsfr = output.mc_best_sfr(is,_);
            bssfr = bsfr/bmass;

            vec1u bmodel = output.mc_best_model(is,_);
            bmetal = btau = bage = bav = bz = replicate(fnan, opts.n_sim);
            for (uint_t im : range(opts.n_sim)) {
                vec1u ids = gridder.grid_ids(bmodel[im]);
                bmetal[im] = output.metal[ids[0]];
                btau[im]   = output.tau[ids[1]];
                bage[im]   = output.age[ids[2]];
                bav[im]    = output.av[ids[3]];
                bz[im]     = output.z[ids[4]];
            }

            if (opts.save_sim) {
                std::string odir = opts.output_dir+"best_fits/";
                if (!file::mkdir(odir)) {
                    warning("could not save simulations");
                    warning("the output directory '", odir, "' could not be created");
                } else {
                    fits::write_table(odir+opts.catalog+"_"+input.id[is]+".sims.fits",
                        "mass", bmass, "sfr", bsfr, "ssfr", bssfr, "metal", bmetal,
                        "tau", btau, "age", bage, "av", bav, "z", bz, "chi2", output.mc_best_chi2
                    );
                }
            }

            for (uint_t ic : range(input.conf_interval)) {
                output.best_mass(is,1+ic)  = inplace_percentile(bmass,  input.conf_interval[ic]);
                output.best_sfr(is,1+ic)   = inplace_percentile(bsfr,   input.conf_interval[ic]);
                output.best_ssfr(is,1+ic)  = inplace_percentile(bssfr,  input.conf_interval[ic]);
                output.best_metal(is,1+ic) = inplace_percentile(bmetal, input.conf_interval[ic]);
                output.best_tau(is,1+ic)   = inplace_percentile(btau,   input.conf_interval[ic]);
                output.best_age(is,1+ic)   = inplace_percentile(bage,   input.conf_interval[ic]);
                output.best_av(is,1+ic)    = inplace_percentile(bav,    input.conf_interval[ic]);
                output.best_z(is,1+ic)     = inplace_percentile(bz,     input.conf_interval[ic]);
            }
        }
    }
}
