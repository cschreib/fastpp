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
        if (opts.verbose) {
            double expsize = input.id.size()*3*sizeof(float);
            for (uint_t i : range(5)) {
                expsize *= gridder.dims[i];
            }

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
        file::write_as<std::uint32_t>(ochi2.out_file, 0);
        file::write_as<std::uint32_t>(ochi2.out_file, input.id.size());
        file::write_as<std::uint32_t>(ochi2.out_file, gridder.dims.size());
        for (uint_t i : range(gridder.dims)) {
            file::write_as<std::uint32_t>(ochi2.out_file, gridder.dims[i]);

            switch (i) {
            case 0: file::write(ochi2.out_file, output.metal); break;
            case 1: file::write(ochi2.out_file, output.tau);   break;
            case 2: file::write(ochi2.out_file, output.age);   break;
            case 3: file::write(ochi2.out_file, output.av);    break;
            case 4: file::write(ochi2.out_file, output.z);     break;
            default: break;
            }
        }

        ochi2.hpos = ochi2.out_file.tellp();

        ochi2.out_file.seekp(0);
        file::write_as<std::uint32_t>(ochi2.out_file, ochi2.hpos);
        ochi2.out_file.seekp(ochi2.hpos);

        // Populate the file with empty data now
        // For each point of the grid, we store 3 values: chi2, mass and sfr
        uint_t nmodel1 = gridder.dims[0]*gridder.dims[1]*gridder.dims[2]*gridder.dims[3];
        vec1f junk = replicate(fnan, gridder.dims[4]*input.id.size()*3);
        for (uint_t im = 0; im < nmodel1; ++im) {
            if (!file::write(ochi2.out_file, junk)) {
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

void fitter_t::write_chi2(uint_t igrid, const vec1f& chi2, const vec1f& mass, const vec1f& sfr) {
    // TODO: configer putting this in a worker thread if this is slowing down too much
    ochi2.out_file.seekp(ochi2.hpos + igrid*input.id.size()*3*sizeof(float));

    file::write(ochi2.out_file, chi2);
    file::write(ochi2.out_file, mass);
    file::write(ochi2.out_file, sfr);
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
        if ((idz.safe[is]  != npos && model.iz != idz.safe[is]) ||
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

        chi2.safe[is] = tchi2;
        mass.safe[is] = scale*model.mass;
        sfr.safe[is]  = scale*model.sfr;

        // Do MC simulation
        if (opts.n_sim > 0) {
            for (uint_t im : range(opts.n_sim)) {
                // Generate "random" fluxes and compute scaling factor
                // NB: since we create the randomness just once at the beginning of the fit
                // all models (and each galaxy) will use the same random numbers
                double wfm = 0;
                for (uint_t il : range(input.lambda)) {
                    // In weighted units, the random perturbations have a sigma of unity
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
                mc_sfr.safe[im]  = scale*model.sfr;
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

    if (save_chi2) {
        write_chi2(model.igrid, chi2, mass, sfr);
    }

    if (!opts.best_from_sim) {
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
        if (!opts.best_from_sim) {
            vec1u ids = gridder.grid_ids(output.best_model[is]);
            output.best_ssfr(is,0)  = output.best_sfr(is,0)/output.best_mass(is,0);
            output.best_metal(is,0) = output.metal[ids[0]];
            output.best_tau(is,0)   = output.tau[ids[1]];
            output.best_age(is,0)   = output.age[ids[2]];
            output.best_av(is,0)    = output.av[ids[3]];
            output.best_z(is,0)     = output.z[ids[4]];
        }

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

            if (opts.best_from_sim) {
                output.best_mass(is,0)  = inplace_median(bmass);
                output.best_sfr(is,0)   = inplace_median(bsfr);
                output.best_ssfr(is,0)  = inplace_median(bssfr);
                output.best_metal(is,0) = inplace_median(bmetal);
                output.best_tau(is,0)   = inplace_median(btau);
                output.best_age(is,0)   = inplace_median(bage);
                output.best_av(is,0)    = inplace_median(bav);
                output.best_z(is,0)     = inplace_median(bz);
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
