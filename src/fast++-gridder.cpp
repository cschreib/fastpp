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

    const uint_t nparam_core = grid_id::custom; // [z,av,av_bc,age,metal]

    nprop = prop_id::custom;
    uint_t nparam_sfh = 0;

    switch (opts.sfh) {
    case sfh_type::gridded:
        nparam_sfh = 1; // [tau]
        nprop += 1;  // [a2t]
        break;
    case sfh_type::single:
        nparam_sfh = 0; // []
        break;
    case sfh_type::custom:
        nparam_sfh = opts.custom_params.size(); // [...]
        break;
    }

    nparam = nparam_core + nparam_sfh;

    output.ifirst_rlum  = nprop;
    nprop += opts.rest_mag.size();
    output.ifirst_abs   = nprop;
    nprop += input.abs_lines.size();
    output.ifirst_ratio = nprop;
    nprop += input.cont_ratios.size();
    output.ifirst_sfhq = nprop;
    nprop += input.sfh_quant.size();

    output.grid.resize(nparam);
    output.param_names.resize(nparam+nprop);
    output.param_descriptions.resize(nparam+nprop);
    output.param_scale.resize(nparam+nprop);
    output.param_log.resize(nparam+nprop);
    output.param_precision.resize(nparam+nprop);

    auto set_param = [&](uint_t id, std::string name, std::string desc,
        bool scale, uint_t lstyle, float precision) {

        output.param_names[id]        = name;
        output.param_descriptions[id] = desc;
        output.param_scale[id]        = scale;
        output.param_log[id]          = lstyle;
        output.param_precision[id]    = precision;
    };

    auto set_prop = [&](uint_t id, std::string name, std::string desc,
        bool scale, uint_t lstyle, float precision) {
        set_param(nparam+id, name, desc, scale, lstyle, precision);
    };

    // Custom grid & properties
    switch (opts.sfh) {
    case sfh_type::gridded:
        output.grid[grid_id::custom] = rgen_step(opts.log_tau_min, opts.log_tau_max, opts.log_tau_step);

        set_param(grid_id::custom, "ltau", "log[tau/yr]",  false, log_style::none, 1e-2);
        set_prop(prop_id::custom,  "la2t", "log[age/tau]", false, log_style::none, 1e-2);
        break;
    case sfh_type::single:
        break;
    case sfh_type::custom:
        for (uint_t ig : range(opts.custom_params)) {
            output.grid[grid_id::custom+ig] = rgen_step(
                opts.custom_params_min[ig], opts.custom_params_max[ig], opts.custom_params_step[ig]
            );

            set_param(grid_id::custom+ig, opts.custom_params[ig], "", false, log_style::none, 1e-2);
        }
        break;
    }

    // Base grid parameters
    output.grid[grid_id::av]    = rgen_step(opts.a_v_min,     opts.a_v_max,     opts.a_v_step);
    output.grid[grid_id::av_bc] = rgen_step(opts.a_v_bc_min,  opts.a_v_bc_max,  opts.a_v_bc_step);
    output.grid[grid_id::age]   = rgen_step(opts.log_age_min, opts.log_age_max, opts.log_age_step);
    output.grid[grid_id::metal] = opts.metal;

    set_param(grid_id::z,       "z",         "",                   false, log_style::none,    1e-4);
    set_param(grid_id::av,      "Av",        "",                   false, log_style::none,    1e-2);
    set_param(grid_id::av_bc,   "Av_bc",     "",                   false, log_style::none,    1e-2);
    set_param(grid_id::age,     "lage",      "log[age/yr]",        false, log_style::none,    1e-2);
    set_param(grid_id::metal,   "metal",     "",                   false, log_style::none,    1e-4);
    set_prop(prop_id::scale,    "lscale",    "log[scale]",         true,  log_style::decimal, 1e-2);
    set_prop(prop_id::spec_scale, "sscale",  "spec_rescale",       true,  log_style::none,    1e-3);
    set_prop(prop_id::mass,     "lmass",     "log[mass/Msol]",     true,  log_style::decimal, 1e-2);
    set_prop(prop_id::sfr,      "lsfr",      "log[sfr/(Msol/yr)]", true,  log_style::decimal, 1e-2);
    set_prop(prop_id::ssfr,     "lssfr",     "log[ssfr*yr]",       false, log_style::decimal, 1e-2);
    set_prop(prop_id::ldust,    "lldust",    "log[lum/Lsol]",      true,  log_style::decimal, 1e-2);
    set_prop(prop_id::ldust_bc, "lldust_bc", "log[lum/Lsol]",      true,  log_style::decimal, 1e-2);
    set_prop(prop_id::lion,     "llion",     "log[lum/Lsol]",      true,  log_style::decimal, 1e-2);
    set_prop(prop_id::mform,    "lmform",    "log[mass/Msol]",     true,  log_style::decimal, 1e-2);

    // Rest luminosities
    for (uint_t i : range(opts.rest_mag)) {
        set_prop(output.ifirst_rlum+i, "M"+to_string(opts.rest_mag[i]),
            "[ABmag]", true, log_style::abmag, 1e-2);
    }

    // Continuum indices
    for (uint_t i : range(input.abs_lines)) {
        auto& l = input.abs_lines[i];
        set_prop(output.ifirst_abs+i, l.name, "[A]", false, log_style::none, 1e-2);
    }
    for (uint_t i : range(input.cont_ratios)) {
        auto& r = input.cont_ratios[i];
        set_prop(output.ifirst_ratio+i, r.name, "", false, log_style::none, 1e-2);
    }

    // SFH parameters
    for (uint_t i : range(input.sfh_quant)) {
        auto& q = input.sfh_quant[i];
        set_prop(output.ifirst_sfhq+i, q.name, q.unit, q.scale, log_style::decimal, 1e-2);
    }

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

    if (opts.fast_spec_convolve) {
        double dz = (max(output_z) - min(output_z))/(1.0 + mean(output_z));
        if (dz > 0.05) {
            warning("delta_z/(1+z) is large (", dz, "); consider turning off FAST_SPEC_CONVOLVE");
        }
    }

    // Pre-compute distances (galaxev templates are in Lsun/A)
    {
        const double dist_Mpc_to_cgs = 3.0856e24; // [cm/Mpc]
        const double lum_sol_to_cgs  = 3.839e33;  // [erg/s/Lsol]
        const double factor = 1e19*lum_sol_to_cgs/sqr(dist_Mpc_to_cgs);
        lum2fl = factor/(4.0*dpi*(1.0+output_z)*sqr(astro::lumdist(output_z, opts.cosmo)));

        const double dist_Mpc_to_si = 3.0856e22; // [m/Mpc]
        const double lum_sol_to_si  = 3.839e26;  // [W/Lsol]
        const double flux_to_uJy    = 1.0e32;    // [uJy/(W/m2/Hz)]
        const double speed_of_light = 2.9979e18; // [A/s]
        const double rffactor = flux_to_uJy*lum_sol_to_si/(speed_of_light*sqr(dist_Mpc_to_si));
        rflum2fl = rffactor/(4.0*dpi*sqr(1e-5));
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
        std::string grid_common = "nmetal="+to_string(output.grid[grid_id::metal].size())+
            ",nage="+to_string(output.grid[grid_id::age].size())+
            ",nav="+to_string(output.grid[grid_id::av].size())+
            ",navbc="+to_string(output.grid[grid_id::av_bc].size())+
            ",nz="+to_string(output_z.size());

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
                    to_string(output.grid[grid_id::custom+ig].size())+",";
            }
            note("fitting a grid of ", nmodel, " templates (", grid_custom, grid_common, ")");
            break;
        }
    }

    output.best_chi2 = replicate(finf, input.id.size());
    output.best_model = replicate(npos, input.id.size());
    output.best_params = replicate(fnan, input.id.size(), nparam+nprop, 1+input.conf_interval.size());

    output.num_models = replicate(0, input.id.size());

    if (opts.n_sim > 0) {
        out.mc_best_chi2 = replicate(finf, input.id.size(), opts.n_sim);
        out.mc_best_model = replicate(npos, input.id.size(), opts.n_sim);
        out.mc_best_props = replicate(fnan, input.id.size(), nprop, opts.n_sim);
    }

    if (!opts.no_cache) {
        // Base library properties
        cache.cache_filename = opts.output_dir+opts.library+"_"+opts.resolution+"_"+
            opts.name_imf+"_"+opts.name_sfh+"_"+opts.dust_law+"_";

        std::string grid_hash = hash(output.grid[_-(grid_id::custom-1)], output.param_names,
            input.lambda, opts.dust_noll_eb, opts.dust_noll_delta, opts.sfr_avg, opts.lambda_ion,
            opts.cosmo.H0, opts.cosmo.wm, opts.cosmo.wL, opts.apply_vdisp, opts.no_igm,
            opts.log_bc_age_max);

        // Additional grid parameter
        switch (opts.sfh) {
        case sfh_type::gridded:
            grid_hash = hash(grid_hash, output.grid[grid_id::custom]);
            break;
        case sfh_type::single:
            break;
        case sfh_type::custom:
            grid_hash = hash(grid_hash, opts.custom_sfh,
                output.grid[grid_id::custom+indgen(nparam-grid_id::custom)]);
            break;
        }

        // Continuum indices
        for (uint_t i : range(input.abs_lines)) {
            auto& l = input.abs_lines[i];
            grid_hash = hash(grid_hash, l.name, l.line_low, l.line_up, l.cont_low, l.cont_up);
        }
        for (uint_t i : range(input.cont_ratios)) {
            auto& r = input.cont_ratios[i];
            grid_hash = hash(grid_hash, r.name, r.cont1_low, r.cont1_up, r.cont2_low, r.cont2_up);
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
    vec1s base_names = {"id", "chi2", "nmodel"};
    for (auto c : opts.output_columns) {
        c = to_lower(c);
        if (!is_any_of(c, to_lower(output.param_names)) && !is_any_of(c, to_lower(output.param_names))
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

        if (l < 0) l = 0;

        return l;
    });

    auto milky_way = vectorize_lambda([](double l) {
        // Reference?

        const double iRv = 1.0/3.1;

        l = 1e4/l; // Angstrom to 1/um
        if (l <= 1.1) {
            l = (0.574 - 0.527*iRv)*pow(l, 1.61);
        } else if (l <= 3.3) {
            l -= 1.82;
            l = 1.0 + (0.17699 + 1.41338*iRv)*l - (0.50447 - 2.28305*iRv)*pow(l,2) -
                (0.02427 - 1.07233*iRv)*pow(l,3) + (0.72085 - 5.38434*iRv)*pow(l,4) +
                (0.01979 - 0.62251*iRv)*pow(l,5) - (0.77530 - 5.3026*iRv)*pow(l,6) +
                (0.32999 - 2.09002*iRv)*pow(l,7);
        } else if (l <= 8.0) {
            double fafb = 0.0;
            if (l > 5.9) {
                fafb = -0.04473*pow(l - 5.9, 2) - 0.009779*pow(l - 5.9, 3)
                + (0.2130*pow(l - 5.9, 2) + 0.1207*pow(l - 5.9, 3))*iRv;
            }

            l = 1.752 - 3.09*iRv + (-0.316 + 1.825*iRv)*l - 0.104/(pow(l - 4.67, 2) + 0.341) +
                1.206*iRv/(pow(l - 4.62, 2) + 0.263) + fafb;
        } else {
            l -= 8.0;
            l = -1.073 + 13.67*iRv - (0.628 - 4.257*iRv)*l +
                (0.137 - 0.42*iRv)*pow(l, 2) + 0.374*iRv*pow(l, 3);
        }

        if (l < 0) l = 0;

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
            da = mean(ptau);
        }

        double db; {
            double l0 = 920.0*(1.0 + z);
            double l1 = 1015.0*(1.0 + z);
            uint_t nstep = 100;
            vec1d tl = rgen(l0, l1, nstep);
            vec1d ptau = exp(-1.7e-3*pow(tl/1026.0, 3.46) - 1.2e-3*pow(tl/972.5, 3.46) -
                9.3e-4*pow(tl/950.0, 3.46));
            db = mean(ptau);
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

void gridder_t::attenuate_and_send(fitter_t& fitter, progress_t& pg,
    const vec1d& lam, const vec1d& tpl_flux, const vec1d& dust_law, const vec2d& igm_abs,
    float lage, vec1u& idm, model_t& model) {

    vec1f& output_av = output.grid[grid_id::av];
    float& model_ldust = model.props[prop_id::ldust];
    float& model_ldust_bc = model.props[prop_id::ldust_bc];
    float& model_lion = model.props[prop_id::lion];

    // Pre-compute bolometric luminosity
    double lbol = integrate(lam, tpl_flux);
    model_ldust = 0.0;
    model_ldust_bc = 0.0;
    model_lion = integrate(lam, tpl_flux, lam.front(), opts.lambda_ion);

    for (uint_t id : range(output_av)) {
        idm[grid_id::av] = id;

        // Apply dust reddening
        vec1f tpl_att_flux = tpl_flux;
        if (output_av[id] > 0) {
            for (uint_t il : range(tpl_att_flux)) {
                tpl_att_flux.safe[il] *= e10(-0.4*output_av.safe[id]*dust_law.safe[il]);
            }

            // Compute absorbed energy
            double lobs = integrate(lam, tpl_att_flux);
            model_ldust = lbol - lobs;
        }

        // Generic treatment
        redshift_and_send(fitter, pg, lam, tpl_att_flux, igm_abs, lage, idm, model);
    }
}

void gridder_t::attenuate_and_send(fitter_t& fitter, progress_t& pg,
    const vec1d& lam, const vec1d& tpl_flux_young, const vec1d& tpl_flux_old,
    const vec1d& dust_law, const vec2d& igm_abs,
    float lage, vec1u& idm, model_t& model) {

    vec1f& output_av = output.grid[grid_id::av];
    vec1f& output_av_bc = output.grid[grid_id::av_bc];
    float& model_ldust = model.props[prop_id::ldust];
    float& model_ldust_bc = model.props[prop_id::ldust_bc];
    float& model_lion = model.props[prop_id::lion];

    // Pre-compute bolometric luminosity
    double lbol_young = integrate(lam, tpl_flux_young);
    double lbol_old = integrate(lam, tpl_flux_old);

    model_ldust = 0.0;
    model_ldust_bc = 0.0;
    model_lion = integrate(lam, tpl_flux_young, lam.front(), opts.lambda_ion)
               + integrate(lam, tpl_flux_old, lam.front(), opts.lambda_ion);

    for (uint_t id : range(output_av)) {
        idm[grid_id::av] = id;

        // Apply dust reddening to old stars
        vec1f tpl_att_flux_old = tpl_flux_old;
        double ldust_old = 0.0;
        if (output_av[id] > 0) {
            for (uint_t il : range(tpl_att_flux_old)) {
                tpl_att_flux_old.safe[il] *= e10(-0.4*output_av.safe[id]*dust_law.safe[il]);
            }

            // Compute absorbed energy
            double lobs = integrate(lam, tpl_att_flux_old);
            ldust_old = lbol_old - lobs;
        }

        for (uint_t id_bc : range(output_av_bc)) {
            idm[grid_id::av_bc] = id_bc;

            // Apply dust reddening to young stars
            vec1f tpl_att_flux_young = tpl_flux_young;
            if (output_av_bc[id_bc] > 0 || output_av[id] > 0) {
                double bc_av = output_av.safe[id] + output_av_bc[id_bc];
                for (uint_t il : range(tpl_att_flux_young)) {
                    tpl_att_flux_young.safe[il] *= e10(-0.4*bc_av*dust_law.safe[il]);
                }

                // Compute absorbed energy
                double lobs = integrate(lam, tpl_att_flux_young);
                double ldust_young = lbol_young - lobs;
                model_ldust = ldust_young + ldust_old;
                model_ldust_bc = ldust_young;
            }

            // Generic treatment
            redshift_and_send(
                fitter, pg, lam, tpl_att_flux_young + tpl_att_flux_old, igm_abs, lage, idm, model
            );
        }
    }
}

void gridder_t::redshift_and_send(fitter_t& fitter, progress_t& pg,
    const vec1d& lam, const vec1d& tpl_att_flux, const vec2d& igm_abs,
    float lage, vec1u& idm, model_t& model) {

    vec1f& output_z = output.grid[grid_id::z];
    float& model_scale = model.props[prop_id::scale];
    float& model_sscale = model.props[prop_id::spec_scale];

    model_scale = 1.0; // by definition
    model_sscale = 1.0; // by definition

    // Compute rest-frame luminosities
    for (uint_t i : range(opts.rest_mag)) {
        model.props[output.ifirst_rlum+i] = rflum2fl*sqr(input.rf_lambda[i])*astro::sed2flux(
            input.rf_filters.safe[i].wl, input.rf_filters.safe[i].tr,
            lam, tpl_att_flux
        );

        if (!is_finite(model.props[output.ifirst_rlum+i])) {
            // Filter goes out of model coverage, assume zero
            model.props[output.ifirst_rlum+i] = 0;
        }
    }

    // Compute absorption lines EW
    for (uint_t i : range(input.abs_lines)) {
        auto& l = input.abs_lines[i];

        double fl = integrate(lam, tpl_att_flux, l.line_low, l.line_up);

        double fc = 0.0;
        if (l.cont_low.size() == 1) {
            // Single window, assume continuum is constant
            fc = integrate(lam, tpl_att_flux, l.cont_low[0], l.cont_up[0])/
                (l.cont_up[0] - l.cont_low[0]);

            // Subtract continuum
            fl -= fc*(l.line_up - l.line_low);
        } else {
            // Two windows, assume continuum is linear
            double l1 = 0.5*(l.cont_low[0] + l.cont_up[0]);
            double l2 = 0.5*(l.cont_low[1] + l.cont_up[1]);

            double f1 = integrate(lam, tpl_att_flux, l.cont_low[0], l.cont_up[0])/
                (l.cont_up[0] - l.cont_low[0]);
            double f2 = integrate(lam, tpl_att_flux, l.cont_low[1], l.cont_up[1])/
                (l.cont_up[1] - l.cont_low[1]);

            fc = interpolate(f1, f2, l1, l2, 0.5*(l.line_low + l.line_up));

            // Subtract continuum
            double a = (f1*l2 - f2*l1)/(l2 - l1);
            double b = (f2 - f1)/(l2 - l1)/2.0;
            fl -= a*(l.line_up - l.line_low) + b*(sqr(l.line_up) - sqr(l.line_low));
        }

        model.props[output.ifirst_abs+i] = -fl/fc;
    }

    // Compute continuum indices
    for (uint_t i : range(input.cont_ratios)) {
        auto& r = input.cont_ratios[i];

        double fc1 = integrate(lam, tpl_att_flux, r.cont1_low, r.cont1_up)/
            (r.cont1_up - r.cont1_low);
        double fc2 = integrate(lam, tpl_att_flux, r.cont2_low, r.cont2_up)/
            (r.cont2_up - r.cont2_low);

        model.props[output.ifirst_ratio+i] = fc2/fc1;
    }

    // Redshift, integrate, and send to fitter
    for (uint_t iz : range(output_z)) {
        idm[grid_id::z] = iz;
        model.igrid = model_id(idm);

        vec1f tpl_att_z_lam = lam;
        vec1f tpl_att_z_flux = tpl_att_flux;

        for (uint_t il : range(tpl_att_z_flux)) {
            // Apply IGM absorption & redshift
            if (opts.no_igm) {
                tpl_att_z_flux.safe[il] *= lum2fl.safe[iz];
            } else {
                tpl_att_z_flux.safe[il] *= lum2fl.safe[iz]*igm_abs.safe(iz,il);
            }
            tpl_att_z_lam.safe[il] *= (1.0 + output_z.safe[iz]);
        }

        // Apply LSF
        convolve_obs(tpl_att_z_lam, tpl_att_z_flux);

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
        // NB we still want to compute the model fluxes for the cache, this is just to
        // decide whether the model is used in the fit
        bool nofit = false;
        if (!opts.no_max_age && lage > auniv[iz]) {
            // Age greater than the age of the universe
            nofit = true;
        }

        if (!nofit && !opts.grid_exclude.empty()) {
            // Custom exclude function
            auto lock = (opts.parallel == parallel_choice::generators ?
                std::unique_lock<std::mutex>(exclude_mutex) : std::unique_lock<std::mutex>());

            for (uint_t i : range(nparam)) {
                exclude_expr.vars.safe[i] = output.grid.safe[i].safe[idm.safe[i]];
            }

            nofit = abs(exclude_expr.eval()) > 0;
        }

        // Send to fitter
        if (!nofit) {
            fitter.fit(model);
        }

        // Cache and print progress
        {
            auto lock = (opts.parallel == parallel_choice::generators ?
                std::unique_lock<std::mutex>(progress_mutex) : std::unique_lock<std::mutex>());

            cache.write_model(model);
            if (opts.verbose) progress_tick(pg, 0.5);
        }
    }
}

void gridder_t::compute_sfh_quantities_impl(const vec1d& ltime, const vec1d& sfh, model_t& model) {
    // Pre-compute some useful things
    vec1d csfh = cumul(ltime, sfh);         // Cumulative SFH
    csfh /= csfh.back();
    uint_t imax = max_id(sfh);              // Index of peak SFR
    double tot_sfr = integrate(ltime, sfh); // Total sum of SFR

    for (uint_t i : range(input.sfh_quant)) {
        auto& q = input.sfh_quant[i];
        auto& o = model.props[output.ifirst_sfhq+i];

        switch (q.type) {
        case sfh_quantity_type::sfr : {
            // Average SFR over the past X yr in [Msun/yr]
            double t1 = min(q.param, ltime.back());
            o = integrate(ltime, sfh, 0.0, t1)/q.param;
            break;
        }
        case sfh_quantity_type::tform : {
            // Time of formation in [yr]
            o = interpolate(ltime, csfh, q.param);
            break;
        }
        case sfh_quantity_type::tsf : {
            // Duration of star formation in [yr]

            // Start from peak SFR and expand in the direction of max SFR until
            // summed SFR equals 68% of the total summed SFR.
            uint_t im0 = imax;
            uint_t im1 = imax;
            double sum = 0.0;
            double prev_length = 0.0;
            double prev_sum = 0.0;
            bool set = false;

            for (uint_t i = 1; i < sfh.size(); ++i) {
                bool left = true;
                if (im0 == 0) {
                    left = false;
                } else if (im1 == sfh.size()-1) {
                    left = true;
                } else {
                    if (sfh.safe[im1+1] > sfh.safe[im0-1]) {
                        left = false;
                    } else {
                        left = true;
                    }
                }

                if (left) {
                    sum += 0.5*(sfh.safe[im0] + sfh.safe[im0-1])*
                                (ltime.safe[im0] - ltime.safe[im0-1]);
                    --im0;
                } else {
                    sum += 0.5*(sfh.safe[im1+1] + sfh.safe[im1])*
                                (ltime.safe[im1+1] - ltime.safe[im1]);
                    ++im1;
                }

                double length = ltime.safe[im1] - ltime.safe[im0];
                if (sum > 0.68*tot_sfr) {
                    o = interpolate(prev_length, length, prev_sum, sum, 0.68*tot_sfr);
                    set = true;
                    break;
                }

                prev_length = length;
                prev_sum = sum;
            }

            if (!set) {
                o = 0.68*ltime.back();
            }

            break;
        }
        case sfh_quantity_type::past_sfr : {
            // Mean SFR during main star formation episode in [Msun/yr]
            double tsf = model.props[output.ifirst_sfhq+q.depends[0]];
            o = 0.68*tot_sfr/tsf;
            break;
        }
        case sfh_quantity_type::tquench : {
            // Time spent since quenching in [yr]
            // Quenching defined as time when SFR is X times lower than the mean past SFR
            double msfr = model.props[output.ifirst_sfhq+q.depends[0]];
            for (uint_t i : range(sfh)) {
                if (sfh.safe[i] > msfr/q.param) {
                    if (i == 0) {
                        o = 0;
                    } else {
                        o = interpolate(ltime.safe[i-1], ltime.safe[i],
                            sfh.safe[i-1], sfh.safe[i], msfr/q.param);
                    }
                    break;
                }
            }
            break;
        }
        case sfh_quantity_type::brate : {
            // Birth rate parameter: SFR_now/SFR_past
            double msfr = model.props[output.ifirst_sfhq+q.depends[0]];
            double sfr = model.props[output.ifirst_sfhq+q.depends[1]];
            o = sfr/msfr;
            break;
        }
        }
    }
}

bool gridder_t::build_and_send(fitter_t& fitter) {
    if (opts.verbose) {
        note("start fitting...");
    }

    if (read_from_cache) {
        model_t model;
        model.flux.resize(input.lambda.size());
        model.props.resize(nprop);

        auto pg = progress_start(nmodel);
        for (uint_t m = 0; m < nmodel; ++m) {
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
        bool ret = false;

        switch (opts.sfh) {
        case sfh_type::gridded:
            ret = build_and_send_ised(fitter);
            break;
        case sfh_type::custom:
            ret = build_and_send_custom(fitter);
            break;
        default:
            error("this SFH is not implemented yet");
            break;
        }

        // Make sure we flush all the cache out
        cache.cache_file.close();

        return ret;
    }

    return true;
}

bool gridder_t::build_template_impl(uint_t iflat, bool nodust,
    vec1f& lam, vec1f& flux, vec1f& iflux) const {

    vec1u idm = grid_ids(iflat);
    uint_t iz = idm[grid_id::z];
    uint_t id = idm[grid_id::av];
    uint_t id_bc = idm[grid_id::av_bc];

    float z = output.grid[grid_id::z][iz];
    float av = output.grid[grid_id::av][id];
    float av_bc = output.grid[grid_id::av_bc][id_bc];

    if (av_bc == 0.0 || nodust) {
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
    } else {
        vec1f flux_young, flux_old;

        switch (opts.sfh) {
        case sfh_type::custom:
            if (!build_template_custom(iflat, lam, flux_young, flux_old)) {
                return false;
            }
            break;
        default:
            error("this SFH is not implemented yet");
            return false;
        }

        vec1d dust_law = build_dust_law(lam);

        // Apply dust reddening
        if (av > 0) {
            for (uint_t il : range(flux_old)) {
                flux_old.safe[il] *= e10(-0.4*av*dust_law.safe[il]);
            }
        }
        if (av+av_bc > 0) {
            for (uint_t il : range(flux_young)) {
                flux_young.safe[il] *= e10(-0.4*(av+av_bc)*dust_law.safe[il]);
            }
        }

        // Combine young and old
        flux = flux_young + flux_old;
    }

    // Compute IGM absorption
    vec2d igm_abs;
    if (!opts.no_igm) {
        igm_abs = build_igm_absorption({z}, lam);
    }

    for (uint_t il : range(flux)) {
        // Apply IGM absorption & redshift
        if (opts.no_igm) {
            flux.safe[il] *= lum2fl[iz];
        } else {
            flux.safe[il] *= lum2fl[iz]*igm_abs.safe(0,il);
        }
        lam.safe[il] *= (1.0 + z);
    }

    // Apply LSF
    convolve_obs(lam, flux);

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

vec2d gridder_t::convolve_function(const vec1d& lam, const vec2d& osed, uint_t n_thread,
    std::function<double(double)> sigma_fun) const {

    vec2d sed = replicate(0.0, osed.dims);

    const uint_t nlam = lam.size();
    const double max_sigma = 5.0;

    auto do_convolve = [&](uint_t lbegin, uint_t lend) {
        for (uint_t l : range(lbegin, lend)) {
            // Get sigma in wavelength units
            double sigma = sigma_fun(lam.safe[l]);

            // Find bounds
            uint_t i0 = l, i1 = l;
            double l0 = lam.safe[l] - max_sigma*sigma;
            while (i0 > 1 && lam.safe[i0] > l0) --i0;
            double l1 = lam.safe[l] + max_sigma*sigma;
            while (i1 < nlam-2 && lam.safe[i1] < l1) ++i1;

            // Integrate
            for (uint_t tl = i0; tl <= i1; ++tl) {
                l0 = 0.5*(lam.safe[tl] + lam.safe[tl-1]);
                l1 = 0.5*(lam.safe[tl] + lam.safe[tl+1]);

                double k = (l1 - l0)*integrate_gauss(l0, l1, lam.safe[l], sigma);
                for (uint_t s : range(osed.dims[0])) {
                    sed.safe(s,l) += osed.safe(s,tl)*k;
                }
            }
        }
    };

    if (n_thread <= 1) {
        // Single-threaded
        do_convolve(0, lam.size());
    } else {
        // Multi-threaded
        auto tp = thread::pool(n_thread);
        uint_t dl = lam.size()/n_thread + 1;
        uint_t l0 = 0;
        uint_t l1 = dl;
        for (uint_t it : range(n_thread)) {
            tp[it].start(do_convolve, l0, l1);
            l0 += dl;
            l1 += dl;
            if (l1 > lam.size()) {
                l1 = lam.size();
            }
        }

        for (uint_t it : range(n_thread)) {
            tp[it].join();
        }
    }

    return sed;
}

void gridder_t::convolve_obs(const vec1f& lam, vec1f& osed) const {
    // We can't use multi-threading if 'generators' is chosen as parallel mode, since
    // we will already be generating multiple models in different threads.
    uint_t nthread = (opts.parallel == parallel_choice::generators ? 1 : opts.n_thread);

    if (!input.tpllsf_sigma.empty() && !opts.fast_spec_convolve) {
        osed = flatten(convolve_function(lam, reform(osed, 1, osed.size()), nthread, [&](double lobs) {
            if (lobs < input.tpllsf_lam.front()) {
                return input.tpllsf_sigma.front();
            } else if (lobs > input.tpllsf_lam.back()) {
                return input.tpllsf_sigma.back();
            } else {
                return interpolate(input.tpllsf_sigma, input.tpllsf_lam, lobs);
            }
        }));
    }
}

void gridder_t::convolve_rest(const vec1d& lam, vec2d& osed) const {
    if (!input.tpllsf_sigma.empty() && opts.fast_spec_convolve) {
        // Approximation: convolve LSF + velocity dispersion in one go for speed
        double vdisp = is_finite(opts.apply_vdisp) ? opts.apply_vdisp : 0.0;
        double z = mean(output.grid[grid_id::z]);
        osed = convolve_function(lam, osed, opts.n_thread, [&](double lrest) {
            double lobs = lrest*(1.0 + z);
            double lsf = 0.0;
            if (lobs < input.tpllsf_lam.front()) {
                lsf = input.tpllsf_sigma.front();
            } else if (lobs > input.tpllsf_lam.back()) {
                lsf = input.tpllsf_sigma.back();
            } else {
                lsf = interpolate(input.tpllsf_sigma, input.tpllsf_lam, lobs);
            }

            return std::sqrt(sqr(lsf/(1.0 + z)) + sqr(lrest*vdisp/2.99792e5));
        });
    } else {
        if (is_finite(opts.apply_vdisp)) {
            osed = convolve_function(lam, osed, opts.n_thread, [&](double lrest) {
                return lrest*opts.apply_vdisp/2.99792e5;
            });
        }
    }
}

void gridder_t::convolve_rest(const vec1f& lam, vec2f& osed) const {
    vec2d dsed = osed;
    convolve_rest(lam, dsed);
    osed = dsed;
}

struct sed_id_pair {
    std::string id;
    uint_t igrid;
    double scale;
};

bool gridder_t::write_seds() const {
    if (opts.make_seds.empty()) return true;

    if (opts.verbose) {
        note("writing list of SEDs to "+opts.output_dir+"seds directory...");
    }

    std::string odir = opts.output_dir+"seds/";
    file::mkdir(odir);

    uint_t nseds = 0; {
        std::ifstream in(opts.make_seds);
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;
            ++nseds;
        }
    }

    auto write_sed = [&](sed_id_pair p) {
        vec1f lam, sed, flx;
        if (!build_template(p.igrid, lam, sed, flx)) {
            return;
        }

        sed *= p.scale;
        flx *= p.scale;

        // Save model
        std::ofstream fout(odir+p.id+".fit");
        fout << "# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        for (uint_t il : range(lam)) {
            fout << std::setw(13) << lam.safe[il] << std::setw(13) << sed.safe[il] << "\n";
        }

        // Save fluxes
        fout.close();
        fout.open(odir+p.id+".input_res.fit");
        fout << "# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        for (uint_t il : range(input.lambda)) {
            fout << std::setw(13) << float(input.lambda.safe[il]) << std::setw(13) << flx.safe[il] << "\n";
        }
    };

    thread::worker_pool<sed_id_pair> pool;
    if (opts.n_thread > 1) {
        pool.start(opts.n_thread, write_sed);
    }

    std::ifstream in(opts.make_seds);
    std::string line;
    uint_t l = 0;
    uint_t written = 0;
    auto pg = progress_start(nseds);
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
        uint_t igrid; double scale;
        if (!from_string(spl[1], igrid)) {
            error("could not convert '", spl[1], "' to integer at ", opts.make_seds, ":", l);
            return false;
        }
        if (!from_string(spl[2], scale)) {
            error("could not convert '", spl[2], "' to double at ", opts.make_seds, ":", l);
            return false;
        }

        if (opts.n_thread > 1) {
            // Parallel
            while (opts.max_queued_fits > 0 &&
                pool.remaining() > opts.max_queued_fits) {
                thread::sleep_for(1e-6);
            }

            pool.process({id, igrid, scale});
        } else {
            // Single threaded
            write_sed({id, igrid, scale});
        }

        ++written;
        if (opts.verbose) progress(pg);
    }

    if (opts.n_thread > 1) {
        while (pool.remaining() != 0) {
            thread::sleep_for(1e-6);
        }

        pool.join();
    }

    if (opts.verbose) {
        note("wrote ", written, " SEDs");
    }

    return true;
}

uint_t gridder_t::model_id(const vec1u& ids) const {
    uint_t fid = 0;

    vif_check(ids.size() == grid_dims.size(), "mismatching IDs with grid");

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
