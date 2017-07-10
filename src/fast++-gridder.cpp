#include "fast++.hpp"

void gridder_t::cache_manager_t::write_model(const model_t& model) {
    if (!cache_file.is_open()) return;

    // TODO: Consider doing this in a worker thread if it slows down execution
    file::write_as<std::uint32_t>(cache_file, model.igrid);
    file::write(cache_file, model.props);
    file::write(cache_file, model.flux);

    if (!cache_file) {
        print("");
        warning("could not write to cache file anymore");
        warning("in case you ran out of disk space, the cache file has been removed");
        print("");
        cache_file.close();
        file::remove(cache_filename);
    }
}

bool gridder_t::cache_manager_t::read_model(model_t& model) {
    if (!cache_file.is_open()) return false;

    file::read_as<std::uint32_t>(cache_file, model.igrid);
    file::read(cache_file, model.props);
    file::read(cache_file, model.flux);

    if (!cache_file) {
        cache_file.close();
        return false;
    }

    return true;
}

gridder_t::gridder_t(const options_t& opt, const input_state_t& inp, output_state_t& out) :
    opts(opt), input(inp), output(out) {

    if (opts.verbose) note("define grid...");

    switch (opts.sfh) {
    case sfh_type::gridded:
        nparam = 5; // [z,av,age,tau,metal]
        nprop = 5;  // [mass,sfr,ssfr,ldust,a2t]
        break;
    case sfh_type::single:
        nparam = 4; // [z,av,age,metal]
        nprop = 4;  // [mass,sfr,ssfr,ldust]
        break;
    case sfh_type::custom:
        nparam = 4+opts.custom_params.size(); // [z,av,age,metal,...]
        nprop = 4;  // [mass,sfr,ssfr,ldust]
        break;
    }

    output.grid.resize(nparam);
    output.param_names.resize(nparam+nprop);
    output.param_descriptions.resize(nparam+nprop);
    output.param_scale.resize(nparam+nprop);
    output.param_log.resize(nparam+nprop);
    output.param_precision.resize(nparam+nprop);

    auto set_param = [&](uint_t id, std::string name, std::string desc,
        bool scale, bool log, float precision) {

        output.param_names[id]        = name;
        output.param_descriptions[id] = desc;
        output.param_scale[id]        = scale;
        output.param_log[id]          = log;
        output.param_precision[id]    = precision;
    };

    auto set_prop = [&](uint_t id, std::string name, std::string desc,
        bool scale, bool log, float precision) {
        set_param(nparam+id, name, desc, scale, log, precision);
    };

    // Custom grid & properties
    switch (opts.sfh) {
    case sfh_type::gridded:
        output.grid[grid_id::custom] = rgen_step(opts.log_tau_min, opts.log_tau_max, opts.log_tau_step);

        set_param(grid_id::custom, "ltau", "log[tau/yr]",  false, false, 1e-2);
        set_prop(prop_id::custom,  "la2t", "log[age/tau]", false, false, 1e-2);
        break;
    case sfh_type::single:
        break;
    case sfh_type::custom:
        for (uint_t ig : range(opts.custom_params)) {
            output.grid[grid_id::custom+ig] = rgen_step(
                opts.custom_params_min[ig], opts.custom_params_max[ig], opts.custom_params_step[ig]
            );

            set_param(grid_id::custom+ig, opts.custom_params[ig], "", false, false, 1e-2);
        }
        break;
    }

    // Base grid parameters
    output.grid[grid_id::av]    = rgen_step(opts.a_v_min,     opts.a_v_max,     opts.a_v_step);
    output.grid[grid_id::age]   = rgen_step(opts.log_age_min, opts.log_age_max, opts.log_age_step);
    output.grid[grid_id::metal] = opts.metal;

    set_param(grid_id::z,     "z",      "",                   false, false, 1e-4);
    set_param(grid_id::av,    "Av",     "",                   false, false, 1e-2);
    set_param(grid_id::age,   "lage",   "log[age/yr]",        false, false, 1e-2);
    set_param(grid_id::metal, "metal",  "",                   false, false, 1e-4);
    set_prop(prop_id::mass,   "lmass",  "log[mass/Msol]",     true,  true,  1e-2);
    set_prop(prop_id::sfr,    "lsfr",   "log[sfr/(Msol/yr)]", true,  true,  1e-2);
    set_prop(prop_id::ssfr,   "lssfr",  "log[ssfr*yr]",       false, true,  1e-2);
    set_prop(prop_id::ldust,  "lldust", "log[ldust/Lsol]",    true,  true,  1e-2);

    // Redshift grid
    vec1f& output_z = output.grid[grid_id::z];

    // First build a new array from the provided parameters
    if (opts.z_step == 0.0 || abs(opts.z_max - opts.z_min) < 0.5*opts.z_step) {
        // A single redshift
        output_z = {opts.z_min};
    } else {
        if (opts.z_step_type == 0) {
            // "dz = cte" grid
            output_z = rgen_step(opts.z_min, opts.z_max, opts.z_step);
        } else {
            // "dz ~ (1+z)" grid
            output_z = e10(rgen_step(
                log10(1.0 + opts.z_min), log10(1.0 + opts.z_max), opts.z_step
            )) - 1.0;
        }
    }

    // If we have fewer galaxies to fit than the size of the above grid,
    // we use the individual zspec/zphot as the grid.
    if (opts.spectrum.empty() && opts.n_sim == 0 && !input.zphot.empty() &&
        input.zphot.size() < output_z.size()) {

        // First compile valid zphot & zspecs
        vec1f cz = input.zspec;
        vec1u idzp = where(!is_finite(cz));
        cz[idzp] = input.zphot[idzp];

        // ... and only keep valid and unique values
        cz = cz[where(is_finite(cz))];
        cz = unique_values(cz, sort(cz));

        output_z = max(cz, 0.00001);
    }

    // Pre-compute distances (galaxev templates are in Lsun/A)
    {
        const double dist_Mpc_to_cgs = 3.0856e24; // [cm/Mpc]
        const double lum_sol_to_cgs  = 3.839e33;  // [erg/s/Lsol]
        const double factor = 1e19*lum_sol_to_cgs/sqr(dist_Mpc_to_cgs);
        lum2fl = factor/(4.0*dpi*(1.0+output_z)*sqr(astro::lumdist(output_z, opts.cosmo)));
    }

    // Pre-compute age of Universe
    auniv = log10(lookback_time(1000.0, opts.cosmo) - lookback_time(output_z, opts.cosmo)) + 9;

    // Pre-compute grid properties and number of parameters
    grid_dims = replicate(0u, output.grid.size());
    grid_dims_pitch = replicate(1u, output.grid.size());
    nmodel = 1;
    nparam = 0;
    ncustom = 1;
    nfreeparam = 0; // TODO: change, this should be 1 (normalization *is* a degree of freedom)
    for (uint_t i : range(grid_dims)) {
        grid_dims[i] = output.grid[i].size();

        for (uint_t j : range(i+1, grid_dims.size())) {
            grid_dims_pitch[i] *= output.grid[j].size();
        }

        nmodel *= grid_dims[i];
        if (i >= grid_id::custom) {
            ncustom *= grid_dims[i];
        }

        ++nparam;
        if (grid_dims[i] > 1) ++nfreeparam;
    }

    if (opts.verbose) {
        std::string grid_common = "nmetal="+strn(output.grid[grid_id::metal].size())+
            ",nage="+strn(output.grid[grid_id::age].size())+
            ",nav="+strn(output.grid[grid_id::av].size())+",nz="+strn(output_z.size());

        switch (opts.sfh) {
        case sfh_type::gridded:
            note("fitting a grid of ", nmodel,
                " templates (ntau=", output.grid[grid_id::custom].size(), ",", grid_common, ")");
            break;
        case sfh_type::single:
            note("fitting a grid of ", nmodel, " templates (", grid_common, ")");
            break;
        case sfh_type::custom:
            std::string grid_custom;
            for (uint_t ig : range(opts.custom_params)) {
                grid_custom += "n"+opts.custom_params[ig]+"="+
                    strn(output.grid[grid_id::custom+ig].size())+",";
            }
            note("fitting a grid of ", nmodel, " templates (", grid_custom, grid_common, ")");
            break;
        }
    }

    output.best_chi2 = replicate(finf, input.id.size());
    output.best_model = replicate(npos, input.id.size());
    output.best_params = replicate(fnan, input.id.size(), nparam+nprop, 1+input.conf_interval.size());

    if (opts.n_sim > 0) {
        out.mc_best_chi2 = replicate(finf, input.id.size(), opts.n_sim);
        out.mc_best_model = replicate(npos, input.id.size(), opts.n_sim);
        out.mc_best_props = replicate(fnan, input.id.size(), nprop, opts.n_sim);
    }

    if (!opts.no_cache) {
        // Base library properties
        cache.cache_filename = opts.output_dir+opts.library+"_"+opts.resolution+"_"+
            opts.name_imf+"_"+opts.name_sfh+"_"+opts.dust_law+"_";

        std::string grid_hash = hash(output.grid[_-(grid_id::custom-1)], input.lambda,
            opts.dust_noll_eb, opts.dust_noll_delta, opts.sfr_avg,
            opts.cosmo.H0, opts.cosmo.wm, opts.cosmo.wL);

        // Additional grid parameter
        switch (opts.sfh) {
        case sfh_type::gridded:
            grid_hash = hash(grid_hash, output.grid[grid_id::custom]);
            break;
        case sfh_type::single:
            break;
        case sfh_type::custom:
            grid_hash = hash(grid_hash, opts.custom_sfh,
                output.grid[grid_id::custom+uindgen(nparam-grid_id::custom)]);
            break;
        }

        cache.cache_filename += grid_hash+".grid";

        if (opts.verbose) {
            note("cache file is '", cache.cache_filename, "'");
        }

        // Open grid file
        if (file::exists(cache.cache_filename)) {
            if (opts.verbose) note("checking cache integrity...");
            cache.cache_file.open(cache.cache_filename, std::ios::binary | std::ios::in);

            cache.cache_file.seekg(0, std::ios_base::end);
            uint_t size = cache.cache_file.tellg();
            uint_t size_expected = nmodel*(
                sizeof(std::uint32_t)+sizeof(float)*(nprop+input.lambda.size())
            );

            if (size != size_expected) {
                warning("cache file is corrupted or invalid, will overwrite it");
                warning("found ", size, " bytes, expected ", size_expected, " (one model is ",
                    size_expected/nmodel, " bytes)");
                cache.cache_file.close();
                read_from_cache = false;
            } else {
                if (opts.verbose) note("cache file exists and seems valid, will use it");
                cache.cache_file.seekg(0, std::ios_base::beg);
            }
        } else {
            read_from_cache = false;
        }

        if (!read_from_cache) {
            cache.cache_file.open(cache.cache_filename, std::ios::binary | std::ios::out);

            if (!cache.cache_file.is_open()) {
                warning("cache file could not be created");
                warning("the program will not use the cache");
            }
        }
    } else {
        read_from_cache = false;
    }
}

bool gridder_t::check_options() const {
    // Check that the requested column names make sense
    bool bad = false;
    vec1s base_names = {"id", "chi2"};
    for (auto c : opts.output_columns) {
        c = tolower(c);
        if (!is_any_of(c, tolower(output.param_names)) && !is_any_of(c, tolower(output.param_names))
            && !is_any_of(c, base_names)) {
            error("unknown column '", c, "'");
            bad = true;
        }
    }

    if (bad) {
        error("some of the requested columns do not exist, cannot proceed");
        return false;
    }

    // If we use a custom SFH, compile it and check that the expression is valid
    if (opts.sfh == sfh_type::custom) {
        vec1s p = {"t", "lage"};
        append(p, opts.custom_params);

        if (!sfh_expr.compile(opts.custom_sfh, p)) {
            return false;
        }
    }

    // If we use an exclude pattern, compile it and check that the expression is valid
    if (!opts.grid_exclude.empty()) {
        vec1s p(nparam);
        for (uint_t ip : range(nparam)) {
            p[ip] = output.param_names[ip];
        }

        if (!exclude_expr.compile(opts.grid_exclude, p)) {
            return false;
        }
    }

    return true;
}

namespace dust {
    auto calzetti2000 = vectorize_lambda([](double l) {
        // http://adsabs.harvard.edu/abs/2000ApJ...533..682C

        const double iRv = 1.0/4.05;

        l *= 1e-4; // Angstrom to um
        if (l <= 0.63) {
            l = (2.659*iRv)*(-2.156 + 1.509/l - 0.198*pow(l, -2) + 0.011*pow(l, -3)) + 1.0;
        } else {
            l = (2.659*iRv)*(-1.857 + 1.040/l) + 1.0;
        }

        return l;
    });

    auto milky_way = vectorize_lambda([](double l) {
        // Reference?

        const double iRv = 1.0/3.1;

        l = 1e-4/l; // Angstrom to 1/um
        if (l <= 1.1) {
            l = (0.574 - 0.527*iRv)*pow(l, 1.61);
        } else if (l <= 3.3) {
            l -= 1.82;
            l = 1.0 + (0.17699 + 1.41338*iRv)*l - (0.50447 - 2.28305*iRv)*pow(l,2) -
                (0.02427 - 1.07233*iRv)*pow(l,3) + (0.72085 - 5.38434*iRv)*pow(l,4) +
                (0.01979 - 0.62251*iRv)*pow(l,5) - (0.77530 - 5.3026*iRv)*pow(l,6) +
                (0.32999 - 2.09002*iRv)*pow(l,7);
        } else if (l <= 8.0) {
            l = 1.752 - 3.09*iRv - (0.316 - 1.825)*l - 0.104/(pow(l - 4.67, 2) + 0.341) +
                1.206*iRv/(pow(l - 4.62, 2) + 0.263) + (l <= 5.9 ? 0 :
                    -0.04473*pow(l - 5.9, 2) - 0.009779*pow(l - 5.9, 3) +
                    (0.2130*pow(l - 5.9, 2) + 0.1207*pow(l - 5.9, 3))*iRv
                );
        } else {
            l = -1.073 + 13.67*iRv - (0.628 - 4.257*iRv)*(l - 8.0) +
                (0.137 - 0.42*iRv)*pow(l - 8.0, 2) + 0.374*iRv*pow(l - 8.0, 3);
        }

        return l;
    });

    auto noll2009 = vectorize_lambda([](double l, double eb, double delta) {
        // http://adsabs.harvard.edu/abs/2009A%26A...507.1793N

        const double iRv = 1.0/4.05;
        const double width2 = pow(350.0, 2);
        const double clam2 = pow(2175.0, 2);

        double l2 = pow(l, 2);
        return (calzetti2000(l) + iRv*eb*l2*width2/((l2 - clam2) + l2*width2))*pow(l/5500.0, delta);
    });
}

namespace igm {
    vec1d madau1995(double z, vec1d lam) {
        // http://adsabs.harvard.edu/abs/1995ApJ...441...18M
        // TODO: check this implementation someday, I suspect this is wrong or
        // very approximate (taken directly from FAST)

        double da; {
            double l0 = 1050.0*(1.0 + z);
            double l1 = 1170.0*(1.0 + z);
            uint_t nstep = 100;
            vec1d tl = rgen(l0, l1, nstep);
            vec1d ptau = exp(-3.6e-3*pow(tl/1216.0, 3.46));
            da = total(ptau)*(l1-l0)/nstep/(120.0*(1.0 + z));
        }

        double db; {
            double l0 = 920.0*(1.0 + z);
            double l1 = 1015.0*(1.0 + z);
            uint_t nstep = 100;
            vec1d tl = rgen(l0, l1, nstep);
            vec1d ptau = exp(-1.7e-3*pow(tl/1026.0, 3.46) - 1.2e-3*pow(tl/972.5, 3.46) -
                9.3e-4*pow(tl/950.0, 3.46));
            db = total(ptau)*(l1-l0)/nstep/(95.0*(1.0 + z));
        }

        for (auto& l : lam) {
            if (l < 912) {
                l = 0.0;
            } else if (l < 1026) {
                l = db;
            } else if (l < 1216) {
                l = da;
            } else {
                l = 1.0;
            }
        }

        return lam;
    }
}

vec1d gridder_t::build_dust_law(const vec1f& lambda) const {
    vec1d dust_law;
    if (opts.dust_law == "calzetti") {
        dust_law = dust::calzetti2000(lambda);
    } else if (opts.dust_law == "mw") {
        dust_law = dust::milky_way(lambda);
    } else if (opts.dust_law == "noll") {
        dust_law = dust::noll2009(lambda, opts.dust_noll_eb, opts.dust_noll_delta);
    } else if (opts.dust_law == "kc") {
        dust_law = dust::noll2009(lambda, 1.0, -0.1);
    }

    return dust_law;
}

vec2d gridder_t::build_igm_absorption(const vec1f& z, const vec1f& lambda) const {
    vec2d igm_abs(z.size(), lambda.size());
    for (uint_t iz : range(z)) {
        igm_abs(iz,_) = igm::madau1995(z.safe[iz], lambda);
    }
    return igm_abs;
}

void gridder_t::build_and_send_impl(fitter_t& fitter, progress_t& pg,
    const vec1d& lam, const vec1d& tpl_flux, const vec1d& dust_law, const vec2d& igm_abs,
    float lage, vec1u& idm, model_t& model) {

    vec1f& output_av = output.grid[grid_id::av];
    vec1f& output_z = output.grid[grid_id::z];
    float& model_ldust = model.props[prop_id::ldust];

    // Pre-compute bolometric luminosity
    double lbol = integrate(lam, tpl_flux);

    for (uint_t id : range(output_av)) {
        idm[grid_id::av] = id;

        // Apply dust reddening
        vec1f tpl_att_flux = tpl_flux;
        if (output_av[id] > 0) {
            for (uint_t il : range(tpl_att_flux)) {
                tpl_att_flux.safe[il] *= e10(-0.4*output_av[id]*dust_law.safe[il]);
            }

            // Compute absorbed energy
            double lobs = integrate(lam, tpl_att_flux);
            model_ldust = lbol - lobs;
        } else {
            model_ldust = 0;
        }

        for (uint_t iz : range(output_z)) {
            idm[grid_id::z] = iz;
            model.igrid = model_id(idm);

            vec1f tpl_att_z_lam = lam;
            vec1f tpl_att_z_flux = tpl_att_flux;

            for (uint_t il : range(tpl_att_z_flux)) {
                // Apply IGM absorption & redshift
                tpl_att_z_flux.safe[il] *= lum2fl.safe[iz]*igm_abs.safe(iz,il);
                tpl_att_z_lam.safe[il] *= (1.0 + output_z.safe[iz]);
            }

            // Integrate
            for (uint_t il : range(input.lambda)) {
                model.flux.safe[il] = astro::sed2flux(
                    input.filters.safe[il].wl, input.filters.safe[il].tr,
                    tpl_att_z_lam, tpl_att_z_flux
                );

                if (!is_finite(model.flux.safe[il])) {
                    // Filter goes out of model coverage, assume zero
                    model.flux.safe[il] = 0;
                }
            }

            // See if we want to use this model or not
            bool nofit = false;
            if (!opts.no_max_age && lage > auniv[iz]) {
                // Age greater than the age of the universe
                nofit = true;
            }

            if (!nofit && !opts.grid_exclude.empty()) {
                // Custom exclude function
                for (uint_t i : range(nparam)) {
                    exclude_expr.vars.safe[i] = output.grid.safe[i].safe[idm.safe[i]];
                }

                nofit = abs(exclude_expr.eval()) > 0;
            }

            // Send to fitter
            if (!nofit) {
                fitter.fit(model);
            }

            // Cache
            cache.write_model(model);

            if (opts.verbose) progress(pg, 131);
        }
    }
}

bool gridder_t::build_and_send(fitter_t& fitter) {
    if (read_from_cache) {
        model_t model;
        model.flux.resize(input.lambda.size());
        model.props.resize(nprop);

        auto pg = progress_start(nmodel);
        for (uint_t i = 0; i < nmodel; ++i) {
            if (cache.read_model(model)) {
                bool nofit = false;
                if (!opts.no_max_age) {
                    vec1u idm = grid_ids(model.igrid);
                    if (output.grid[grid_id::age][idm[grid_id::age]] > auniv[idm[grid_id::z]]) {
                        nofit = true;
                    }
                }

                if (!nofit && !opts.grid_exclude.empty()) {
                    // Custom exclude function
                    vec1u idm = grid_ids(model.igrid);
                    for (uint_t i : range(nparam)) {
                        exclude_expr.vars.safe[i] = output.grid.safe[i].safe[idm.safe[i]];
                    }

                    nofit = abs(exclude_expr.eval()) > 0;
                }

                if (!nofit) {
                    // Send to fitter
                    fitter.fit(model);
                }

                if (opts.verbose) progress(pg, 131);
            } else {
                print("");
                error("could not read data from cache file");
                error("the cache is probably corrupted, please remove it and try again");
                print("");
                return false;
            }
        }
    } else {
        switch (opts.sfh) {
        case sfh_type::gridded:
            return build_and_send_ised(fitter);
        case sfh_type::custom:
            return build_and_send_custom(fitter);
        default:
            error("this SFH is not implemented yet");
            return false;
        }

        // Make sure we flush all the cache out
        cache.cache_file.close();
    }

    return true;
}

bool gridder_t::build_template_impl(uint_t iflat, bool nodust,
    vec1f& lam, vec1f& flux, vec1f& iflux) const {

    vec1u idm = grid_ids(iflat);
    uint_t iz = idm[grid_id::z];
    uint_t id = idm[grid_id::av];

    float z = output.grid[grid_id::z][iz];
    float av = output.grid[grid_id::av][id];

    switch (opts.sfh) {
    case sfh_type::gridded:
        if (!build_template_ised(iflat, lam, flux)) {
            return false;
        }
        break;
    case sfh_type::custom:
        if (!build_template_custom(iflat, lam, flux)) {
            return false;
        }
        break;
    default:
        error("this SFH is not implemented yet");
        return false;
    }

    // Apply dust reddening
    if (av > 0 && !nodust) {
        vec1d dust_law = build_dust_law(lam);

        for (uint_t il : range(flux)) {
            flux.safe[il] *= e10(-0.4*av*dust_law.safe[il]);
        }
    }

    // Compute IGM absorption
    vec2d igm_abs = build_igm_absorption({z}, lam);

    for (uint_t il : range(flux)) {
        // Apply IGM absorption & redshift
        flux.safe[il] *= lum2fl[iz]*igm_abs.safe(0,il);
        lam.safe[il] *= (1.0 + z);
    }

    // Integrate
    iflux.resize(input.lambda.size());
    for (uint_t il : range(input.lambda)) {
        iflux[il] = astro::sed2flux(
            input.filters[il].wl, input.filters[il].tr,
            lam, flux
        );
    }

    return true;
}

bool gridder_t::build_template(uint_t igrid, vec1f& lam, vec1f& flux, vec1f& iflux) const {
    return build_template_impl(igrid, false, lam, flux, iflux);
}

bool gridder_t::build_template_nodust(uint_t igrid, vec1f& lam, vec1f& flux, vec1f& iflux) const {
    return build_template_impl(igrid, true, lam, flux, iflux);
}

bool gridder_t::get_sfh(uint_t igrid, const vec1d& t, vec1d& sfh) const {
    switch (opts.sfh) {
    case sfh_type::gridded:
        return get_sfh_ised(igrid, t, sfh, opts.sfh_output);
    case sfh_type::custom:
        return get_sfh_custom(igrid, t, sfh, opts.sfh_output);
    default:
        error("this SFH is not implemented yet");
        return false;
    }
}


bool gridder_t::write_seds() const {
    if (opts.make_seds.empty()) return true;

    if (opts.verbose) {
        note("writing list of SEDs to "+opts.output_dir+"seds directory");
    }

    std::string odir = opts.output_dir+"seds/";
    file::mkdir(odir);

    std::ifstream in(opts.make_seds);
    std::string line;
    uint_t l = 0;
    uint_t written = 0;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split_any_of(line, " \t");
        if (spl.size() != 3) {
            error("ill formed line at ", opts.make_seds, ":", l, ": expected three values, got ",
                spl.size());
            return false;
        }

        std::string id = spl[0];
        uint_t igrid; double mass;
        if (!from_string(spl[1], igrid)) {
            error("could not convert '", spl[1], "' to integer at ", opts.make_seds, ":", l);
            return false;
        }
        if (!from_string(spl[2], mass)) {
            error("could not convert '", spl[2], "' to double at ", opts.make_seds, ":", l);
            return false;
        }

        vec1f lam, sed, flx;
        if (!build_template(igrid, lam, sed, flx)) {
            return false;
        }

        sed *= mass;
        flx *= mass;

        // Save model
        // TODO: use a more efficient serialization code
        std::ofstream fout(odir+id+".fit");
        fout << "# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        for (uint_t il : range(lam)) {
            fout << align_right(strn(lam.safe[il]), 13)
                 << align_right(strn(sed.safe[il]), 13) << "\n";
        }
        fout.close();

        // Save fluxes
        fout.open(odir+id+".input_res.fit");
        fout << "# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        for (uint_t il : range(input.lambda)) {
            fout << align_right(strn(float(input.lambda[il])), 13)
                 << align_right(strn(flx[il]), 13) << "\n";
        }
        fout.close();

        ++written;
    }

    if (opts.verbose) {
        note("wrote ", written, " SEDs");
    }

    return true;
}

uint_t gridder_t::model_id(const vec1u& ids) const {
    uint_t fid = 0;

    phypp_check(ids.size() == grid_dims.size(), "mismatching IDs with grid");

    for (uint_t i : range(ids)) {
        fid += grid_dims_pitch[i]*ids[i];
    }

    return fid;
}

vec1u gridder_t::grid_ids(uint_t iflat) const {
    vec1u idm(grid_dims_pitch.size());
    for (uint_t i : range(grid_dims)) {
        uint_t j = grid_dims.size()-1-i;
        idm.safe[j] = iflat % grid_dims.safe[j];
        iflat /= grid_dims.safe[j];
    }

    return idm;
}
