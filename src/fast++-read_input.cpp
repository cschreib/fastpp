#include "fast++.hpp"

extern const char* fastpp_version;

std::string remove_first_last(std::string val, std::string charlist) {
    if (val.empty()) return val;
    uint_t p0 = 0, n = val.size();
    if (charlist.find_first_of(val[0]) != charlist.npos) {
        ++p0; --n;
    }
    if (n > 0 && charlist.find_first_of(val[val.size()-1]) != charlist.npos) {
        --n;
    }

    return val.substr(p0, n);
}

template <typename T>
bool parse_value_impl(const std::string& val, T& out) {
    return from_string(val, out);
}

bool parse_value_impl(const std::string& val, std::string& out) {
    out = remove_first_last(val, "\'\"");
    return true;
}

bool parse_value_impl(std::string val, parallel_choice& out) {
    val = trim(to_lower(remove_first_last(val, "\'\"")));
    if (val == "none" || val.empty()) {
        out = parallel_choice::none;
    } else if (val == "sources") {
        out = parallel_choice::sources;
    } else if (val == "models") {
        out = parallel_choice::models;
    } else if (val == "generators") {
        out = parallel_choice::generators;
    } else {
        error("unknown parallelization choice '", val, "'");
        error("must be one of 'none', sources', 'models' or 'generators'");
        return false;
    }
    return true;
}

template <typename T>
bool parse_value_impl(std::string val, vec<1,T>& out) {
    if (val.empty()) return true;
    val = remove_first_last(val, "[]");
    vec1s spl = trim(split(val, ","));

    if (spl.size() == 1 && spl[0].empty()) {
        out.clear();
        return true;
    }

    out.resize(spl.size());
    for (uint_t i : range(spl)) {
        if (!parse_value_impl(spl[i], out[i])) {
            return false;
        }
    }

    return true;
}

template <typename T>
bool parse_value(const std::string& key, const std::string& val, T& out) {
    if (!parse_value_impl(val, out)) {
        error("could not parse value of parameter ", key);
        note("could not convert '", val, "' into ", pretty_type_t(T));
        return false;
    }

    return true;
}

bool read_params(options_t& opts, input_state_t& state, const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open()) {
        error("could not open parameter file '", filename, "'");
        return false;
    }

    opts.cosmo = astro::get_cosmo("std");

    vec1s unparsed_key, unparsed_val;

    auto do_parse = [&](const std::string& key, const std::string& val) {
        #define PARSE_OPTION(name) if (key == #name) { return parse_value(key, val, opts.name); }
        #define PARSE_OPTION_RENAME(opt, name) if (key == name) { return parse_value(key, val, opts.opt); }

        PARSE_OPTION(catalog)
        PARSE_OPTION(ab_zeropoint)
        PARSE_OPTION(filters_res)
        PARSE_OPTION_RENAME(filters_format, "filter_format")
        PARSE_OPTION(temp_err_file)
        PARSE_OPTION(name_zphot)
        PARSE_OPTION(spectrum)
        PARSE_OPTION(auto_scale)
        PARSE_OPTION(output_dir)
        PARSE_OPTION(output_file)
        PARSE_OPTION(n_sim)
        PARSE_OPTION(c_interval)
        PARSE_OPTION(best_fit)
        PARSE_OPTION(library_dir)
        PARSE_OPTION(library)
        PARSE_OPTION(resolution)
        PARSE_OPTION_RENAME(name_imf, "imf")
        PARSE_OPTION_RENAME(name_sfh, "sfh")
        PARSE_OPTION(dust_law)
        PARSE_OPTION_RENAME(dust_noll_eb, "e_b")
        PARSE_OPTION_RENAME(dust_noll_delta, "delta")
        PARSE_OPTION(my_sfh)
        PARSE_OPTION(log_tau_min)
        PARSE_OPTION(log_tau_max)
        PARSE_OPTION(log_tau_step)
        PARSE_OPTION(log_age_min)
        PARSE_OPTION(log_age_max)
        PARSE_OPTION(log_age_step)
        PARSE_OPTION(no_max_age)
        PARSE_OPTION(z_min)
        PARSE_OPTION(z_max)
        PARSE_OPTION(z_step)
        PARSE_OPTION(z_step_type)
        PARSE_OPTION(a_v_min)
        PARSE_OPTION(a_v_max)
        PARSE_OPTION(a_v_step)
        PARSE_OPTION(metal)
        PARSE_OPTION_RENAME(cosmo.H0, "h0")
        PARSE_OPTION_RENAME(cosmo.wm, "omega_m")
        PARSE_OPTION_RENAME(cosmo.wL, "omega_l")
        PARSE_OPTION(save_chi_grid)

        // Not in original FAST
        PARSE_OPTION(force_zphot)
        PARSE_OPTION(best_at_zphot)
        PARSE_OPTION(zphot_conf)
        PARSE_OPTION(save_sim)
        PARSE_OPTION(best_from_sim)
        PARSE_OPTION(no_cache)
        PARSE_OPTION(parallel)
        PARSE_OPTION(n_thread)
        PARSE_OPTION(max_queued_fits)
        PARSE_OPTION(verbose)
        PARSE_OPTION(debug)
        PARSE_OPTION(sfr_avg)
        PARSE_OPTION(intrinsic_best_fit)
        PARSE_OPTION(best_sfhs)
        PARSE_OPTION(sfh_output_step)
        PARSE_OPTION(sfh_output)
        PARSE_OPTION(use_lir)
        PARSE_OPTION(output_columns)
        PARSE_OPTION(custom_sfh)
        PARSE_OPTION(custom_sfh_step)
        PARSE_OPTION(custom_sfh_lookback)
        PARSE_OPTION(custom_params)
        PARSE_OPTION(grid_exclude)
        PARSE_OPTION(no_igm)
        PARSE_OPTION(make_seds)
        PARSE_OPTION(lambda_ion)
        PARSE_OPTION(save_bestchi)
        PARSE_OPTION(apply_vdisp)
        PARSE_OPTION(rest_mag)
        PARSE_OPTION(continuum_indices)
        PARSE_OPTION(sfh_quantities)

        #undef  PARSE_OPTION
        #undef  PARSE_OPTION_RENAME

        unparsed_key.push_back(key);
        unparsed_val.push_back(val);

        return true;
    };

    // Read the parameter file line by line
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        // Remove inline comments
        auto cp = line.find_first_of('#');
        if (cp != line.npos) {
            line = line.substr(0, cp);
        }

        // Split on '=' and parse key and values
        auto eqp = line.find_first_of('=');
        if (eqp == line.npos) {
            error("ill formed line in configuration file");
            note(line);
            return false;
        }

        std::string key = to_lower(trim(line.substr(0, eqp)));
        std::string val = trim(line.substr(eqp+1));

        if (!do_parse(key, val)) {
            return false;
        }
    }

    // Initialize custom parameter grid
    opts.custom_params_min = opts.custom_params_max = opts.custom_params_step =
        replicate(fnan, opts.custom_params.size());

    // Check unparsed key/value pairs for additional options we couldn't parse before
    for (uint_t p : range(unparsed_key)) {
        std::string key = unparsed_key[p];
        std::string val = unparsed_val[p];
        bool read = false;

        // Read custom grid parameters
        for (uint_t c : range(opts.custom_params)) {
            if (key == to_lower(opts.custom_params[c])+"_min") {
                read = true;
                if (!parse_value(key, val, opts.custom_params_min[c])) {
                    return false;
                }
            } else if (key == to_lower(opts.custom_params[c])+"_max") {
                read = true;
                if (!parse_value(key, val, opts.custom_params_max[c])) {
                    return false;
                }
            } else if (key == to_lower(opts.custom_params[c])+"_step") {
                read = true;
                if (!parse_value(key, val, opts.custom_params_step[c])) {
                    return false;
                }
            }
        }

        if (!read) {
            warning("unknown parameter '", to_upper(key), "'");
        }
    }

    if (opts.verbose) {
        note("this is FAST++ version '", fastpp_version, "'");
    }

    // Create output directory, if it doesn't exist
    if (!opts.output_dir.empty()) {
        opts.output_dir = file::directorize(opts.output_dir);
        if (!file::mkdir(opts.output_dir)) {
            error("could not create output directory '", opts.output_dir, "'");
            error("make sure you have the proper rights");
            return false;
        }
    }

    if (opts.output_file.empty()) opts.output_file = opts.catalog;

    // Now check for the consistency of the input and make corrections when necessary

    vec1s possible_dust_laws = {"calzetti", "mz", "noll", "kc"};
    if (!is_any_of(opts.dust_law, possible_dust_laws)) {
        error("unknown dust law '", opts.dust_law, "'");
        error("possible values: ", possible_dust_laws);
        return false;
    }

    if (opts.sfr_avg < 0 || !is_finite(opts.sfr_avg)) {
        opts.sfr_avg = 0;
    }

    auto check_sfh_step = [&](float& step, std::string which) {
        if (step < 0) {
            error("step of ", which, " SFH cannot be negative");
            return false;
        } else if (step > 1e4) {
            error("step of ", which, " SFH must be given in Myr");
            error("(you gave ", step, " which is larger than the age of the Universe)");
            return false;
        }

        step *= 1e6; // convert to yr

        return true;
    };

    check_sfh_step(opts.sfh_output_step, "output");
    check_sfh_step(opts.custom_sfh_step, "custom");

    if (opts.verbose) {
        if (opts.sfr_avg == 0) {
            note("using instantaneous SFRs");
        } else {
            note("averaging SFRs over ", opts.sfr_avg, " Myr");
        }

        if (opts.no_igm) {
            note("not applying IGM absorption");
        } else {
            note("applying Madau+1995 IGM absorption");
        }
    }

    opts.sfh_output = to_lower(opts.sfh_output);

    if (!(opts.sfh_output == "sfr" || opts.sfh_output == "mass")) {
        error("'SFH_OUTPUT' must be either 'sfr' or 'mass'");
        return false;
    }

    opts.sfr_avg *= 1e6;

    if (opts.intrinsic_best_fit && !opts.best_fit) {
        warning("'INTRINSIC_BEST_FIT' has no effect if 'BEST_SIM' is not set to 1");
    }

    if (opts.best_sfhs && !opts.my_sfh.empty()) {
        warning("cannot output best fit SFH when using 'MY_SFH'");
        opts.best_sfhs = false;
    }

    if (opts.best_from_sim && opts.n_sim == 0) {
        error("cannot use the option 'BEST_FROM_SIM' if simulations are not enabled (set N_SIM > 0)");
        return false;
    }

    if (!opts.my_sfh.empty()) {
        opts.sfh = sfh_type::single;
    } else if (!opts.custom_sfh.empty()) {
        opts.sfh = sfh_type::custom;
    } else {
        opts.sfh = sfh_type::gridded;
    }

    if (opts.spectrum.empty()) {
        opts.auto_scale = false;
    }

    if (opts.parallel != parallel_choice::none && opts.n_thread <= 1) {
        warning("parallelization is enabled but the number of thread is set to zero");
        warning("parallelization will be disabled unless you set N_THREAD > 1");
        opts.parallel = parallel_choice::none;
    }

    if (opts.best_at_zphot && opts.force_zphot) {
        note("BEST_AT_ZPHOT=1 is automatically true if FORCE_ZPHOT=1, "
            "so you do not have to specify both");
    }

    for (double c : opts.c_interval) {
        if (abs(c - 68.0) > 0.01 && abs(c - 95.0) > 0.01 && abs(c - 99.0) > 0.01) {
            error("confidence interval must be one of 68, 95 or 99% (got ", c, ")");
            return false;
        }
    }

    if (opts.n_sim != 0) {
        state.conf_interval = 0.5*(1.0 - opts.c_interval/100.0);
        inplace_sort(opts.c_interval);
        vec1f cint = state.conf_interval;
        state.conf_interval.clear();
        for (float c : cint) {
            state.conf_interval.push_back(c);
            state.conf_interval.push_back(1.0 - c);
        }
    } else {
        opts.c_interval.clear();
    }

    if (opts.force_zphot || abs(opts.zphot_conf) < 0.1) {
        opts.zphot_conf = fnan;
    } else {
        opts.zphot_conf /= 100.0;
    }

    if (opts.apply_vdisp <= 0.0) {
        opts.apply_vdisp = fnan;
    }

    if (is_finite(opts.apply_vdisp) && opts.verbose) {
        note("convolving templates with sigma=", opts.apply_vdisp, " km/s");
    }

    if (opts.metal.empty()) {
        if (opts.library == "bc03" || opts.library == "ma05") {
            opts.metal = {0.02};
        } else if (opts.library == "co11") {
            opts.metal = {0.019};
        } else {
            // Unknown stellar library, simply guess what could be a good default value...
            opts.metal = {0.02};
        }
    }

    // SFH parameters
    if (!opts.sfh_quantities.empty()) {
        opts.sfh_quantities = to_lower(opts.sfh_quantities);

        // Scan the list and add dependencies

        // tquench and brate need past_sfr
        if (count(opts.sfh_quantities == "tquench") > 0 ||
            count(begins_with(opts.sfh_quantities, "brate")) > 0) {
            opts.sfh_quantities.push_back("past_sfr");
        }

        // brateX needs sfrX
        for (auto c : opts.sfh_quantities[where(begins_with(opts.sfh_quantities, "brate"))]) {
            opts.sfh_quantities.push_back("sfr"+erase_begin(c, "brate"));
        }

        // past_sfr needs tsf
        if (count(opts.sfh_quantities == "past_sfr") > 0) {
            opts.sfh_quantities.push_back("tsf");
        }

        // Remove duplicates
        opts.sfh_quantities = unique_values(opts.sfh_quantities);

        // Sort by dependency
        vec1u dependency = vectorize_lambda([](const std::string& q) {
            if (begins_with(q, "brate"))   return 5;
            if (begins_with(q, "sfr"))     return 4;
            if (begins_with(q, "tquench")) return 3;
            if (q == "past_sfr")           return 2;
            if (q == "tsf")                return 1;
            return 0;
        })(opts.sfh_quantities);

        opts.sfh_quantities = opts.sfh_quantities[sort(dependency)];

        // Now add the parameters
        for (auto c : opts.sfh_quantities) {
            sfh_quantity_t q;
            q.name = c;

            if (c == "tsf") {
                q.type = sfh_quantity_type::tsf;
                q.unit = "log[yr]";
            } else if (begins_with(c, "past_sfr")) {
                q.type = sfh_quantity_type::past_sfr;
                q.unit = "log[Msol/yr]";
                q.scale = true;
                q.depends.push_back(where_first(opts.sfh_quantities == "tsf"));
            } else if (begins_with(c, "tquench")) {
                c = erase_begin(c, "tquench");
                if (!from_string(c, q.param)) {
                    error("could not read 'tquench' parameter from '", q.name, "'");
                    return false;
                }
                if (q.param <= 0) {
                    error("'tquench' parameter must be > 0 (got ", q.param,
                        " from '", q.name, "')");
                    return false;
                }

                q.type = sfh_quantity_type::tquench;
                q.unit = "log[yr]";
                q.depends.push_back(where_first(opts.sfh_quantities == "past_sfr"));
            } else if (begins_with(c, "tform")) {
                c = erase_begin(c, "tform");
                if (!from_string(c, q.param)) {
                    error("could not read 'tform' parameter from '", q.name, "'");
                    return false;
                }
                if (q.param <= 0 || q.param > 100) {
                    error("'tform' parameter must be > 0 and <= 100 (got ", q.param,
                        " from '", q.name, "')");
                    return false;
                }

                q.param = 1.0 - q.param*1e-2; // [%] to fraction

                q.type = sfh_quantity_type::tform;
                q.unit = "log[yr]";
            } else if (begins_with(c, "sfr")) {
                c = erase_begin(c, "sfr");
                if (!from_string(c, q.param)) {
                    error("could not read 'sfr' parameter from '", q.name, "'");
                    return false;
                }
                if (q.param <= 0) {
                    error("'sfr' parameter must be > 0 (got ", q.param, " from '", q.name, "')");
                    return false;
                }

                q.param *= 1e6; // [Myr] to [yr]

                q.type = sfh_quantity_type::sfr;
                q.unit = "log[Msol/yr]";
                q.scale = true;
            } else if (begins_with(c, "brate")) {
                c = erase_begin(c, "brate");
                if (!from_string(c, q.param)) {
                    error("could not read 'brate' parameter from '", q.name, "'");
                    return false;
                }
                if (q.param <= 0) {
                    error("'brate' parameter must be > 0 (got ", q.param, " from '", q.name, "')");
                    return false;
                }

                q.param *= 1e6; // [Myr] to [yr]

                q.type = sfh_quantity_type::brate;
                q.depends.push_back(where_first(opts.sfh_quantities == "past_sfr"));
                q.depends.push_back(where_first(opts.sfh_quantities == "sfr"+c));
            } else {
                error("unknown SFH quantity '", q.name, "'");
                return false;
            }

            vif_check(count(q.depends >= state.sfh_quant.size()) == 0,
                "dependent SFH quantity will be processed after, this is a bug "
                "(quantities: ", opts.sfh_quantities,")");

            state.sfh_quant.push_back(q);
        }

        if (opts.verbose) note("computing SFH quantities: ", opts.sfh_quantities);
    }

    if (opts.output_columns.empty()) {
        // Default columns to display
        switch (opts.sfh) {
        case sfh_type::gridded:
            opts.output_columns = {
                "id", "z", "ltau", "metal", "lage", "Av", "lmass", "lsfr", "lssfr", "la2t"
            };
            break;
        case sfh_type::single:
            opts.output_columns = {
                "id", "z", "metal", "lage", "Av", "lmass", "lsfr", "lssfr"
            };
            break;
        case sfh_type::custom:
            opts.output_columns = {
                "id", "z", "metal", "lage", "Av", "lmass", "lsfr", "lssfr"
            };
            append(opts.output_columns, opts.custom_params);
            break;
        }

        if (!opts.rest_mag.empty()) {
            append(opts.output_columns, "M"+to_string_vector(opts.rest_mag));
        }

        if (!opts.sfh_quantities.empty()) {
            append(opts.output_columns, to_lower(opts.sfh_quantities));
        }

        opts.output_columns.push_back("chi2");
    } else {
        // Check that column names are OK
        if (count(opts.output_columns == "") != 0) {
            error("OUTPUT_COLUMNS contains empty column names, please fix and re-run");
            error("if you want to use the default columns, set this option to an empty array: []");
            return false;
        }
    }

    if (!opts.custom_params.empty()) {
        // Add tau to the custom parameter grid if the user used it
        uint_t idp = where_first(to_lower(opts.custom_params) == "log_tau");
        if (idp != npos) {
            opts.custom_params_min[idp] = opts.log_tau_min;
            opts.custom_params_max[idp] = opts.log_tau_max;
            opts.custom_params_step[idp] = opts.log_tau_step;
        }
    }

    if (!opts.custom_sfh.empty()) {
        // Check that all parameters in the grid have been properly defined
        bool bad = false;
        for (uint_t c : range(opts.custom_params)) {
            if (!is_finite(opts.custom_params_min[c])) {
                error("missing ", to_upper(opts.custom_params[c])+"_MIN value");
                bad = true;
            }
            if (!is_finite(opts.custom_params_max[c])) {
                error("missing ", to_upper(opts.custom_params[c])+"_MAX value");
                bad = true;
            }
            if (!is_finite(opts.custom_params_step[c])) {
                error("missing ", to_upper(opts.custom_params[c])+"_STEP value");
                bad = true;
            }
        }

        if (bad) return false;
    }

    if (!opts.make_seds.empty() && !file::exists(opts.make_seds)) {
        error("file '", opts.make_seds, "' is empty (MAKE_SEDS=...)");
        error("if you do not want to generate these SEDs, remove the MAKE_SEDS=... option");
        return false;
    }

    // Use the default 'share' directory if nothing is provided
    if (opts.filters_res.empty()) {
        opts.filters_res = std::string(FASTPP_SHARE_DIR)+"/FILTER.RES.latest";
    }

    if (opts.temp_err_file.empty()) {
        opts.temp_err_file = std::string(FASTPP_SHARE_DIR)+"/TEMPLATE_ERROR.fast.v0.2";
    }

    if (opts.library_dir.empty()) {
        opts.library_dir = std::string(FASTPP_SHARE_DIR)+"/libraries/";
    } else {
        opts.library_dir = file::directorize(opts.library_dir);
    }

    return true;
}

bool read_header(const std::string& filename, vec1s& header) {
    std::ifstream in(filename);
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);

        // Find the first non-empty comment line
        if (line.empty() || line[0] != '#') continue;

        line = trim(line.substr(1));
        if (line.empty()) continue;

        // Split column names by spaces
        header = split_any_of(line, " \t\n\r");
        return true;
    }

    error("missing header in '", filename, "'");
    note("the header line must start with # and list the column names");

    return false;
}

bool adjust_filter(const options_t& opts, fast_filter_t& f) {
    // Convert filter to the right definition
    if (opts.filters_format == 1) {
        f.tr *= f.wl;
    } else {
        f.tr *= sqr(f.wl);
    }

    // Compute filter integral
    double ttot = integrate(f.wl, f.tr);
    if (!is_finite(ttot) || ttot == 0) {
        error("filter ", f.id, " has zero or invalid througput");
        return false;
    }

    // Normalize filter to unit integral
    f.tr /= ttot;

    // Check the units are OK
    double wl0 = integrate(f.wl, f.wl*f.tr);
    if (wl0 < 100.0) {
        error("the cental wavelength of filter ", f.id, " is unexpectedly low (", wl0, ")");
        error("please check the wavelenth axis is in Angstroms");
        return false;
    }

    return true;
}

bool read_filters(const options_t& opts, input_state_t& state) {
    if (opts.filters_res.empty()) {
        // No filter, maybe just trying to fit a spectrum?
        return true;
    }

    if (opts.verbose) {
        note("reading filters from '", opts.filters_res, "'");
    }

    std::ifstream in(opts.filters_res);
    if (!in.is_open()) {
        error("could not open filter file '", opts.filters_res, "'");
        return false;
    }

    state.rf_filters.resize(opts.rest_mag.size());

    // Read the required filters from the database
    std::string line; uint_t l = 0;
    vec1u idcat, idfil;
    uint_t ntotfilt = 0;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty()) continue;

        vec1s spl = split_any_of(line, " \t\n\r");
        uint_t npts;
        if (!from_string(spl[0], npts)) {
            warning("could not understand l.", l, " in ", opts.filters_res);
            continue;
        }

        // Start of a new filter
        fast_filter_t filt;

        ++ntotfilt;
        filt.id = ntotfilt;

        // Determine if this filter is used in the catalog
        bool doread = false;
        bool cat_filt = false;
        bool rf_filt = false;

        vec1u idused = where(state.no_filt == filt.id);
        if (!idused.empty()) {
            doread = true;
            cat_filt = true;
        }

        uint_t idrf = where_first(opts.rest_mag == filt.id);
        if (idrf != npos) {
            doread = true;
            rf_filt = true;
        }

        uint_t lcnt = 0;
        while (lcnt < npts && std::getline(in, line)) {
            ++l;

            if (line.empty()) continue;

            if (doread) {
                // Reading the filter response line by line
                spl = split_any_of(line, " \t\n\r");
                float wl, tr;
                if (spl.size() != 3 || !from_string(spl[1], wl) || !from_string(spl[2], tr)) {
                    error("could not parse values from line ", l);
                    error("expected '[id] [wavelength] [throughput]', got '", line, "'");
                    note("reading '", opts.filters_res, "'");
                    return false;
                }

                filt.wl.push_back(wl);
                filt.tr.push_back(tr);
            }

            ++lcnt;
        }

        if (cat_filt) {
            // Keep the ID aside for later sorting
            append(idcat, idused);
            append(idfil, replicate(state.filters.size(), idused.size()));

            // Add filter to database
            state.filters.push_back(std::move(filt));
        }

        if (rf_filt) {
            // Add filter to database
            state.rf_filters[idrf] = std::move(filt);
        }
    }

    if (opts.verbose) {
        note("found ", ntotfilt, " filters in the database, will use ",
            state.filters.size(), " of them");
        if (!state.rf_filters.empty()) {
            note("... plus ", state.rf_filters.size(), " for rest-frame magnitudes");
        }
    }

    // Sort the filter list as in the catalog
    state.filters = state.filters[idfil[sort(idcat)]];

    // Make sure we are not missing any
    vec1u notfound;
    if (state.filters.size() != state.no_filt.size()) {
        for (uint_t b : state.no_filt) {
            bool found = false;
            for (auto& f : state.filters) {
                if (f.id == b) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                notfound.push_back(b);
            }
        }
    }

    for (uint_t f : range(opts.rest_mag)) {
        if (state.rf_filters[f].tr.empty()) {
            notfound.push_back(opts.rest_mag[f]);
        }
    }

    if (!notfound.empty()) {
        error("filters ", notfound, " are missing from the filter library");
        return false;
    }

    // Normalize filters and compute the central wavelength
    for (auto& f : state.filters) {
        if (!adjust_filter(opts, f)) {
            return false;
        }

        state.lambda.push_back(integrate(f.wl, f.wl*f.tr));
    }

    for (auto& f : state.rf_filters) {
        if (!adjust_filter(opts, f)) {
            return false;
        }

        state.rf_lambda.push_back(integrate(f.wl, f.wl*f.tr));
    }

    return true;
}

bool read_fluxes(const options_t& opts, input_state_t& state) {
    std::string catalog_file = opts.catalog+".cat";

    if (opts.verbose) {
        note("reading fluxes from '", catalog_file, "'");
    }

    if (!file::exists(catalog_file)) {
        error("could not open photometric catalog '", catalog_file, "'");
        return false;
    }

    // Read the header of the photometric catalog
    vec1s header;
    if (!read_header(catalog_file, header)) {
        return false;
    }

    // Read and apply the column translation file, if it exists
    vec1s header_trans = header;
    std::string translate_file = opts.catalog+".translate";
    if (file::exists(translate_file)) {
        if (opts.verbose) {
            note("using column translation file '", translate_file, "'");
        }

        vec1s tr_from, tr_to;
        ascii::read_table(translate_file, tr_from, tr_to);

        vec1u idh, idt;
        match(header, tr_from, idh, idt);
        header_trans[idh] = tr_to[idt];
    }

    header_trans = to_upper(header_trans);

    // Parse the filter IDs from the (translated) header
    vec1u col_flux, col_eflux;
    for (uint_t ih : range(header)) {
        vec2s ext = regex_extract(header_trans[ih], "^F([0-9]+)$");
        if (ext.empty()) continue;

        uint_t ihe = where_first(header_trans == "E"+ext[0]);
        if (ihe == npos) {
            warning("flux column ", header[ih], " has no uncertainty (E"+ext[0]+") "
                "and will be ignored");
            continue;
        }

        uint_t tid;
        from_string(ext[0], tid);
        state.no_filt.push_back(tid);
        col_flux.push_back(ih);
        col_eflux.push_back(ihe);
    }

    // Check that we have all the columns we need
    uint_t col_id = where_first(header_trans == "ID");
    uint_t col_zspec = where_first(header_trans == "Z_SPEC");
    uint_t col_tot = where_first(is_any_of(header_trans, "TOT"+to_string_vector(state.no_filt)));

    if (col_id == npos) {
        error("missing ID column in photometric catalog");
        return false;
    }

    // Read filters from the database
    if (!read_filters(opts, state)) {
        return false;
    }

    if (opts.verbose) {
        note("reading fluxes...");
    }

    // Read all lines to determine the number of galaxies
    uint_t ngal = 0; {
        std::ifstream in(catalog_file);
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;

            ++ngal;
        }
    }

    // Resize arrays now to avoid reallocation later
    state.id.resize(ngal);
    state.zspec = replicate(fnan, ngal);
    state.flux = replicate(fnan, ngal, state.no_filt.size());
    state.eflux = replicate(fnan, ngal, state.no_filt.size());

    // Now read the catalog itself, only keeping the columns we are interested in
    uint_t l = 0;
    uint_t gid = 0;
    std::ifstream in(catalog_file);
    std::string line;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split_any_of(line, " \t\n\r");

        if (spl.size() != header_trans.size()) {
            error("line ", l, " has ", spl.size(), " columns while header has ", header_trans.size());
            return false;
        }

        // Read the ID
        state.id.safe[gid] = spl[col_id];

        // Read the zspec if any
        if (col_zspec != npos) {
            float tz;
            if (!from_string(spl[col_zspec], tz)) {
                error("could not read z_spec from line ", l);
                note("must be a floating point number, got: '", spl[col_zspec], "'");
                return false;
            }

            // Negative means "no zspec"
            if (tz < 0) {
                tz = fnan;
            }

            state.zspec.safe[gid] = tz;
        }

        // Read the fluxes and uncertainties
        vec1f flx, err;
        vec1b good = from_string(spl[col_flux], flx);
        if (count(!good) != 0) {
            for (uint_t b : where(!good)) {
                error("could not read flux (", header[col_flux[b]], ") from line ", l);
                note("must be a floating point number, got: '", spl[col_flux[b]], "'");
            }

            return false;
        }

        good = from_string(spl[col_eflux], err);
        if (count(!good) != 0) {
            for (uint_t b : where(!good)) {
                error("could not read uncertainty (", header[col_eflux[b]], ") from line ", l);
                note("must be a floating point number, got: '", spl[col_eflux[b]], "'");
            }

            return false;
        }

        // Apply scaling to total fluxes if requested
        if (col_tot != npos) {
            float totcor;
            if (!from_string(spl[col_tot], totcor)) {
                error("could not read total flux scaling (", header[col_tot], ") from line ", l);
                note("must be a floating point number, got: '", spl[col_tot], "'");
                return false;
            }

            flx *= totcor;
            err *= totcor;
        }

        // Flag bad values
        vec1u idb = where(err < 0 || !is_finite(flx) || !is_finite(err));
        err.safe[idb] = finf; flx.safe[idb] = 0;

        if (idb.size() == flx.size()) {
            warning("object ", state.id.safe[gid], " (l.", l, ") has no valid photometry");
        }

        // Save flux and uncertainties in the input state
        state.flux.safe(gid,_) = flx;
        state.eflux.safe(gid,_) = err;

        ++gid;
    }

    // Convert photometry from [catalog unit] to [uJy]
    float abzp = e10(0.4*(23.9 - opts.ab_zeropoint));
    // Convert photometry from fnu [uJy] to flambda [1e-19 x erg/s/cm2/A]
    vec1u idbb = indgen(state.no_filt.size());
    for (uint_t i : range(state.id)) {
        state.flux(i,idbb) = 1e19*abzp*astro::uJy2cgs(state.lambda*1e-4, state.flux(i,idbb));
        state.eflux(i,idbb) = 1e19*abzp*astro::uJy2cgs(state.lambda*1e-4, state.eflux(i,idbb));
    }

    if (opts.verbose) {
        note("fitting ", state.flux.dims[0], " source", (state.flux.dims[0] > 1 ? "s" : ""),
            " with ", state.flux.dims[1], " fluxes", (state.flux.dims[0] > 1 ? " each" : ""));
    }

    return true;
}

bool read_spectra(const options_t& opts, input_state_t& state) {
    if (opts.spectrum.empty()) {
        return true;
    }

    if (opts.verbose) {
        note("reading spectra from '", opts.spectrum, "'");
    }

    if (!file::exists(opts.spectrum+".spec")) {
        error("could not open spectral catalog '", opts.spectrum, ".spec'");
        return false;
    }

    // Read the header
    vec1s spec_header;
    if (!read_header(opts.spectrum+".spec", spec_header)) {
        return false;
    }

    spec_header = to_upper(spec_header);

    // Check we have the right columns there
    uint_t col_wl0 = where_first(spec_header == "WL_LOW");
    uint_t col_wl1 = where_first(spec_header == "WL_UP");
    uint_t col_wl  = where_first(spec_header == "WL");
    uint_t col_tr  = where_first(spec_header == "TR");
    uint_t col_bin = where_first(spec_header == "BIN");

    // First check if the user is using the FAST-IDL format
    if ((col_wl0 == npos || col_wl1 == npos) && col_wl != npos) {
        col_wl0 = npos;
        col_wl1 = npos;
    } else {
        col_wl = npos;

        // Then check for missing wavelength columns
        if (col_wl0 == npos) {
            error("missing lower wavelength column (WL_LOW) in spectral catalog");
            return false;
        }
        if (col_wl1 == npos) {
            error("missing upper wavelength column (WL_UP) in spectral catalog");
            return false;
        }
    }

    vec1u sid;
    vec1u col_flux, col_eflux;
    for (uint_t b : range(spec_header)) {
        if (spec_header[b][0] != 'F') continue;

        std::string id = erase_begin(spec_header[b], "F");
        uint_t cid = where_first(to_upper(state.id) == id);
        if (cid == npos) {
            warning("spectrum for source ", id, " has no corresponding photometry "
                "and will be ignored");
            continue;
        }

        uint_t ce = where_first(spec_header == "E"+id);
        if (ce == npos) {
            warning("spectral flux column ", spec_header[b], " has no uncertainty (E"+id+") "
                "and will be ignored");
            continue;
        }

        sid.push_back(cid);
        col_flux.push_back(b);
        col_eflux.push_back(ce);
    }

    // Read the catalog
    vec1u bin;
    vec2f sflx, serr;
    vec1f slam0, slam1;

    uint_t l = 0;
    std::ifstream in(opts.spectrum+".spec");
    std::string line;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split_any_of(line, " \t\n\r");

        if (spl.size() != spec_header.size()) {
            error("line ", l, " has ", spl.size(), " columns while header has ", spec_header.size());
            return false;
        }

        if (col_tr != npos) {
            bool tr;
            if (!from_string(spl[col_tr], tr)) {
                error("could not read spectral transmission from line ", l);
                note("must be a boolean (0 or 1), got: '", spl[col_tr], "'");
                return false;
            }

            // Transmission is binary: 0 = ignore, 1 = use
            if (!tr) continue;
        }

        if (col_bin != npos) {
            uint_t tbin;
            if (!from_string(spl[col_bin], tbin)) {
                error("could not read spectral binning from line ", l);
                note("must be a positive integer (>= 0), got: '", spl[col_bin], "'");
                return false;
            }

            bin.push_back(tbin);
        }

        if (col_wl0 != npos) {
            double wl0;
            if (!from_string(spl[col_wl0], wl0)) {
                error("could not read spectral lower wavelength from line ", l);
                note("must be a floating point number, got: '", spl[col_wl0], "'");
                return false;
            }

            slam0.push_back(wl0);
        }

        if (col_wl1 != npos) {
            double wl1;
            if (!from_string(spl[col_wl1], wl1)) {
                error("could not read spectral upper wavelength from line ", l);
                note("must be a floating point number, got: '", spl[col_wl1], "'");
                return false;
            }

            slam1.push_back(wl1);
        }

        if (col_wl != npos) {
            double wl;
            if (!from_string(spl[col_wl], wl)) {
                error("could not read spectral wavelength from line ", l);
                note("must be a floating point number, got: '", spl[col_wl], "'");
                return false;
            }

            slam0.push_back(wl);
        }

        // Read the fluxes and uncertainties
        vec1f flx, err;
        vec1b good = from_string(spl[col_flux], flx);
        if (count(!good) != 0) {
            for (uint_t b : where(!good)) {
                error("could not read spectral flux (", spec_header[col_flux[b]],
                    ") from line ", l);
                note("must be a floating point number, got: '", spl[col_flux[b]], "'");
            }

            return false;
        }

        good = from_string(spl[col_eflux], err);
        if (count(!good) != 0) {
            for (uint_t b : where(!good)) {
                error("could not read spectral uncertainty (", spec_header[col_eflux[b]],
                    ") from line ", l);
                note("must be a floating point number, got: '", spl[col_eflux[b]], "'");
            }

            return false;
        }

        // Flag bad values
        vec1u idb = where(err < 0 || !is_finite(flx) || !is_finite(err));
        err[idb] = finf; flx[idb] = 0;

        // Add fluxes to the list
        append<1>(sflx, reform(flx, flx.size(), 1));
        append<1>(serr, reform(err, err.size(), 1));
    }

    if (opts.verbose) {
        note("found ", sflx.dims[0], " spectr", (sflx.dims[0] > 1 ? "a" : "um"),
            " with ", sflx.dims[1], " spectral elements", (sflx.dims[0] > 1 ? " each" : ""));
    }

    // Compute upper/lower wavelength if only WL was provided, assuming contiguous coverage and
    // mostly uniform binning as advised in the FAST documentation
    if (col_wl != npos) {
        double lam0 = slam0.front();
        double lam1 = slam0.back();
        double dlam_whole = (slam0.back()-slam0.front())/(slam0.size()-1);

        for (uint_t il : range(1, slam0.size())) {
            if (abs((slam0[il]-slam0[il-1])/dlam_whole - 1.0) > 0.05) {
                error("the wavelength grid of the input spectrum is not uniform");
                error("if this was intended, please use the new WL_LOW/WL_UP synthax to "
                    "specify each spectral element unambiguously");
                return false;
            }
        }

        slam0 = rgen(lam0, lam1, slam0.size()) - dlam_whole*0.5;
        slam1 = slam0 + dlam_whole;
    }

    // Flag bad values
    vec1u idb = where(serr < 0 || !is_finite(sflx) || !is_finite(serr));
    serr[idb] = finf; sflx[idb] = 0;

    // Apply binning if required
    if (!bin.empty()) {
        // First sort by bin ID
        vec1u sortid = sort(bin);
        sflx = sflx(_,sortid);
        serr = serr(_,sortid);
        bin = bin[sortid];
        slam0 = slam0[sortid];
        slam1 = slam1[sortid];

        // Now combine photometry by inverse variance weighting
        vec2f oflx = sflx, oerr = serr;
        sflx.clear(); serr.clear();
        vec1f oslam0 = slam0, oslam1 = slam1;
        slam0.clear(); slam1.clear();

        vec1f tw(oflx.dims[0]), tf(oflx.dims[0]);
        uint_t lastb = npos;
        uint_t nb = 0;
        double min_lam = finf, max_lam = -finf;
        for (uint_t b : range(oslam0)) {
            if (bin[b] != lastb) {
                if (nb > 0) {
                    // Save the previous bin
                    append<1>(sflx, reform(tf/tw, tf.size(), 1));
                    append<1>(serr, reform(1/sqrt(tw), tw.size(), 1));
                    slam0.push_back(min_lam);
                    slam1.push_back(max_lam);
                }

                // Cleanup for the next one
                tw[_] = 0;
                tf[_] = 0;
                nb = 0;
                lastb = bin[b];
                min_lam = finf;
                max_lam = -finf;
            }

            // Accumulate the flux with weighting
            vec1f w = 1/sqr(oerr(_,b));
            tf += oflx(_,b)*w;
            tw += w;
            ++nb;

            // The width of the effective spectral bin is taken
            // from the extrema of the binned elements, even if
            // they are not contiguous... be warned
            min_lam = min(min_lam, oslam0[b]);
            max_lam = max(max_lam, oslam1[b]);
        }

        if (nb > 0) {
            // Save the last bin, if any
            append<1>(sflx, reform(tf/tw, tf.size(), 1));
            append<1>(serr, reform(1/sqrt(tw), tw.size(), 1));
            slam0.push_back(min_lam);
            slam1.push_back(max_lam);
        }

        vec1u idb = where(!is_finite(serr) || !is_finite(sflx));
        serr[idb] = finf; sflx[idb] = 0;

        if (opts.verbose) {
            if (oflx.dims[1] == sflx.dims[1]) {
                note("binned and original spectra resolution are the same (", sflx.dims[1],
                    " spectral elements)");
            } else {
                note("rebinned spectra from ", oflx.dims[1], " to ", sflx.dims[1],
                    " spectral elements");
            }
        }
    } else {
        if (opts.verbose) {
            note("fitting spectra at original resolution (", sflx.dims[1], " spectral elements)");
        }
    }

    // Create synthetic filters
    for (uint_t b : range(sflx.dims[1])) {
        fast_filter_t f;
        f.spectral = true;
        if (state.no_filt.size() > 0){
          f.id = max(state.no_filt)+1 + b;
        }
        else{
          f.id = 1 + b;
        }
        f.wl = {slam0[b], slam1[b]};
        f.tr = {1.0f, 1.0f};

        double ttot = integrate(f.wl, f.tr);
        if (!is_finite(ttot) || ttot == 0) {
            error("synthetic filter ", b, " for spectral data (wavelength ", slam0[b], " to ",
                slam1[b], " A) has zero or invalid througput");
            return false;
        }

        f.tr /= ttot;

        state.filters.push_back(f);
        state.lambda.push_back(0.5f*(slam0[b] + slam1[b]));
    }

    // Merge the spectra into the flux catalog
    vec2f pflx = std::move(state.flux);
    vec2f perr = std::move(state.eflux);
    uint_t nplam = pflx.dims[1];
    uint_t nslam = sflx.dims[1];
    uint_t ngal = pflx.dims[0];

    state.flux = replicate(0.0f, ngal, nplam+nslam);
    state.eflux = replicate(finf, ngal, nplam+nslam);

    if (nplam > 0){
      for (uint_t i : range(ngal)) {
          state.flux(i,_-(nplam-1)) = pflx(i,_);
          state.eflux(i,_-(nplam-1)) = perr(i,_);
      }
    }

    for (uint_t i : range(sid)) {
        state.flux(sid[i],nplam-_) = sflx(i,_);
        state.eflux(sid[i],nplam-_) = serr(i,_);
    }

    state.spec_start = nplam;
    state.spec_end = nplam+nslam;

    return true;
}

bool read_photoz(const options_t& opts, input_state_t& state) {
    std::string catalog_file = opts.catalog+".zout";
    if (!file::exists(catalog_file)) {
        return true;
    }

    if (opts.verbose) {
        note("reading photometric redshifts from '", catalog_file, "'");
    }

    // Read the header
    vec1s header;
    if (!read_header(catalog_file, header)) {
        return false;
    }

    header = to_upper(header);

    // Check we have all the columns we need
    uint_t col_zphot = where_first(header == to_upper(opts.name_zphot));
    uint_t col_zspec = where_first(header == "Z_SPEC");
    uint_t col_id    = where_first(header == "ID");

    vec1u col_up, col_low;
    if (is_finite(opts.zphot_conf)) {
        for (float c : {0.68, 0.95, 0.99}) {
            uint_t cl = where_first(header == "L"+to_string(round(100.0*c)));
            uint_t cu = where_first(header == "U"+to_string(round(100.0*c)));
            if (cl == npos) {
                error("could not find redshift confidence interval column L"+to_string(round(100.0*c)));
                continue;
            } else if (cu == npos) {
                error("could not find redshift confidence interval column U"+to_string(round(100.0*c)));
                continue;
            }

            col_up.push_back(cu);
            col_low.push_back(cl);
            state.zphot_conf.push_back(c);
        }

        {
            vec1u ids = sort(state.zphot_conf);
            state.zphot_conf = state.zphot_conf[ids];
            col_up = col_up[ids];
            col_low = col_low[ids];
        }

        if (!is_any_of(round(100*opts.zphot_conf), round(100*state.zphot_conf))) {
            error("missing zphot confidence intervals (", round(100*opts.zphot_conf), "th "
                "percentiles)");
            return false;
        }
    }

    if (col_zphot == npos) {
        error("missing zphot (", to_upper(opts.name_zphot), ") column in photometric redshift file");
        return false;
    }
    if (col_id == npos) {
        error("missing ID column in photometric redshift file");
        return false;
    }

    // Initialize the zphot columns
    state.zphot = replicate(fnan, state.id.size(), 1 + 2*state.zphot_conf.size());

    // Read the catalog
    uint_t l = 0;
    uint_t i = 0;
    std::ifstream in(catalog_file);
    std::string line;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split_any_of(line, " \t\n\r");

        std::string id = spl[col_id];

        float zp;
        if (!from_string(spl[col_zphot], zp)) {
            error("could not read photometric redshift from line ", l);
            note("must be a floating point number, got: '", spl[col_zphot], "'");
            return false;
        }

        if (id != state.id[i]) {
            error("photometric redshift and photometry catalogs do not match");
            return false;
        }

        state.zphot(i,0) = (zp > 0 ? zp : fnan);

        if (col_zspec != npos && !is_finite(state.zspec[i])) {
            if (!from_string(spl[col_zspec], zp)) {
                error("could not read spectroscopic redshift from line ", l);
                note("must be a floating point number, got: '", spl[col_zspec], "'");
                return false;
            }

            state.zspec[i] = (zp > 0 ? zp : fnan);
        }

        for (uint_t c : range(col_low)) {
            if (!from_string(spl[col_low[c]], zp)) {
                error("could not read photometric redshift lower ", round(100.0*state.zphot_conf[c]),
                    "th confidence boundary from line ", l);
                note("must be a floating point number, got: '", spl[col_low[c]], "'");
                return false;
            }

            state.zphot(i,1+2*c+0) = (zp > 0 ? zp : fnan);
        }

        for (uint_t c : range(col_up)) {
            if (!from_string(spl[col_up[c]], zp)) {
                error("could not read photometric redshift upper ", round(100.0*state.zphot_conf[c]),
                    "th confidence boundary from line ", l);
                note("must be a floating point number, got: '", spl[col_up[c]], "'");
                return false;
            }

            state.zphot(i,1+2*c+1) = (zp > 0 ? zp : fnan);
        }

        ++i;
    }

    if (i != state.id.size()) {
        error("photometric redshift and photometry catalogs do not match (", i, " vs. ",
            state.id.size(), ")");
        return false;
    }

    return true;
}


bool read_lir(const options_t& opts, input_state_t& state) {
    if (!opts.use_lir) return true;

    std::string catalog_file = opts.catalog+".lir";
    if (!file::exists(catalog_file)) {
        return true;
    }

    if (opts.verbose) {
        note("reading infrared luminosities from '", catalog_file, "'");
    }

    // Read the header
    vec1s header;
    if (!read_header(catalog_file, header)) {
        return false;
    }

    header = to_upper(header);

    // Check we have all the columns we need
    uint_t col_lir = where_first(header == "LIR");
    uint_t col_err = where_first(header == "ELIR");
    uint_t col_id  = where_first(header == "ID");

    if (col_lir == npos) {
        error("missing LIR column in infrared luminosity file");
        return false;
    }
    if (col_err == npos) {
        error("missing ELIR column in infrared luminosity file");
        return false;
    }
    if (col_id == npos) {
        error("missing ID column in infrared luminosity file");
        return false;
    }

    // Initialize the lir columns
    state.lir = replicate(fnan, state.id.size());
    state.lir_err = replicate(fnan, state.id.size());

    // Read the catalog
    uint_t l = 0;
    uint_t i = 0;
    std::ifstream in(catalog_file);
    std::string line;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split_any_of(line, " \t\n\r");

        std::string id = spl[col_id];

        float lir;
        if (!from_string(spl[col_lir], lir)) {
            error("could not read infrared luminosity from line ", l);
            note("must be a floating point number, got: '", spl[col_lir], "'");
            return false;
        }

        float err;
        if (!from_string(spl[col_err], err)) {
            error("could not read infrared luminosity uncertainty from line ", l);
            note("must be a floating point number, got: '", spl[col_err], "'");
            return false;
        }

        if (i >= state.id.size() || id != state.id[i]) {
            error("infrared luminosity and photometry catalogs do not match");
            return false;
        }

        if (err > 0) {
            state.lir[i] = lir;
            state.lir_err[i] = err;
        }

        ++i;
    }

    if (i != state.id.size()) {
        error("infrared luminosity and photometry catalogs do not match (", i, " vs. ",
            state.id.size(), ")");
        return false;
    }

    return true;
}

bool read_template_error(const options_t& opts, input_state_t& state) {
    if (opts.temp_err_file.empty()) {
        return true;
    }

    if (opts.verbose) {
        note("apply template library error function to photometry");
    }

    if (!file::exists(opts.temp_err_file)) {
        error("could not open template error function file '", opts.temp_err_file, "'");
        return false;
    }

    ascii::read_table(opts.temp_err_file, state.tplerr_lam, state.tplerr_err);

    if (state.tplerr_lam.empty()) {
        error("template error function is empty, something must be wrong in the file");
        return false;
    }

    return true;
}

bool check_input(options_t& opts, input_state_t& state) {
    if (!state.zspec.empty()) {
        // Check that zspecs are covered by the redshift grid
        vec1u idb = where(state.zspec < opts.z_min);
        if (!idb.empty()) {
            warning("found ", idb.size(), " galaxies with z_spec lower than "
                "the minimum value of the grid (z < ", opts.z_min, ")");
            warning("will put them at z=", opts.z_min);

            state.zspec[idb] = opts.z_min;
        }

        idb = where(state.zspec > opts.z_max);
        if (!idb.empty()) {
            warning("found ", idb.size(), " galaxies with z_spec larger than "
                "the maximum value of the grid (z > ", opts.z_max, ")");
            warning("will put them at z=", opts.z_max);

            state.zspec[idb] = opts.z_max;
        }
    }

    if (!state.zphot.empty()) {
        // Check that zphots are covered by the redshift grid
        for (uint_t c : range(state.zphot.dims[1])) {
            vec1u idb = where(state.zphot(_,c) < opts.z_min);
            if (!idb.empty()) {
                if (c == 0) {
                    warning("found ", idb.size(), " galaxies with z_phot lower than "
                        "the minimum value of the grid (z < ", opts.z_min, ")");
                    warning("will put them at z=", opts.z_min);
                }

                state.zphot(idb,c) = opts.z_min;
            }

            idb = where(state.zphot(_,c) > opts.z_max);
            if (!idb.empty()) {
                if (c == 0) {
                    warning("found ", idb.size(), " galaxies with z_phot larger than "
                        "the maximum value of the grid (z > ", opts.z_max, ")");
                    warning("will put them at z=", opts.z_max);
                }

                state.zphot(idb,c) = opts.z_max;
            }
        }
    }

    return true;
}

bool read_continuum_indices(const options_t& opts, input_state_t& state) {
    if (opts.continuum_indices.empty()) {
        return true;
    }

    if (opts.verbose) {
        note("using continuum indices defined in '", opts.continuum_indices, "'");
    }

    std::ifstream in(opts.continuum_indices);
    if (!in.is_open()) {
        error("could not open continuum indices file '", opts.continuum_indices, "'");
        return false;
    }

    std::string line;
    uint_t l = 0;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split(line, "=");
        std::string name = trim(spl[0]);

        spl = trim(split(spl[1], ","));
        std::string type = spl[0];

        if (type == "abs") {
            if ((spl.size()-1) != 4 && (spl.size()-1) != 6) {
                error("too few parameters for 'abs' continuum index (need at 4 or 6, got ",
                    spl.size()-1, ")");
                return false;
            }

            state.abs_lines.push_back(absorption_line_t{});
            auto& tline = state.abs_lines.back();
            tline.name = name;

            if (!from_string(spl[1], tline.line_low)) {
                error("could not read lower wavelength of absorption line from '", spl[1], "' in",
                    opts.continuum_indices, ":", l);
                return false;
            }
            if (!from_string(spl[2], tline.line_up)) {
                error("could not read upper wavelength of absorption line from '", spl[2], "' in",
                    opts.continuum_indices, ":", l);
                return false;
            }

            uint_t ncont = (spl.size()-3)/2;
            tline.cont_low.resize(ncont);
            tline.cont_up.resize(ncont);
            for (uint_t c : range(ncont)) {
                if (!from_string(spl[c*2+3], tline.cont_low[c])) {
                    error("could not read lower wavelength of continuum from '", spl[c*2+3], "' in",
                        opts.continuum_indices, ":", l);
                    return false;
                }
                if (!from_string(spl[c*2+4], tline.cont_up[c])) {
                    error("could not read upper wavelength of continuum from '", spl[c*2+4], "' in",
                        opts.continuum_indices, ":", l);
                    return false;
                }
            }
        } else if (type == "ratio") {
            if ((spl.size()-1) != 4) {
                error("incorrect number of parameters for 'ratio' continuum index (need 4, got ",
                    spl.size()-1, ")");
                return false;
            }

            state.cont_ratios.push_back(continuum_ratio_t{});
            auto& ratio = state.cont_ratios.back();
            ratio.name = name;

            if (!from_string(spl[1], ratio.cont1_low)) {
                error("could not read lower wavelength of continuum ratio from '", spl[1], "' in",
                    opts.continuum_indices, ":", l);
                return false;
            }
            if (!from_string(spl[2], ratio.cont1_up)) {
                error("could not read upper wavelength of continuum ratio from '", spl[2], "' in",
                    opts.continuum_indices, ":", l);
                return false;
            }
            if (!from_string(spl[3], ratio.cont2_low)) {
                error("could not read lower wavelength of continuum ratio from '", spl[3], "' in",
                    opts.continuum_indices, ":", l);
                return false;
            }
            if (!from_string(spl[4], ratio.cont2_up)) {
                error("could not read upper wavelength of continuum ratio from '", spl[4], "' in",
                    opts.continuum_indices, ":", l);
                return false;
            }
        } else {
            error("unknown continuum index type '", type, "' in ", opts.continuum_indices, ":", l);
            return false;
        }
    }

    return true;
}

bool read_input(options_t& opts, input_state_t& state, const std::string& filename) {
    // First read options from the parameter file
    if (!read_params(opts, state, filename)) {
        return false;
    }

    // Read the photometry + filters
    if (!read_fluxes(opts, state)) {
        return false;
    }

    // Read the spectra, if any
    if (!read_spectra(opts, state)) {
        return false;
    }

    // Read the photometric redshift catalog from EAzY, if any
    if (!read_photoz(opts, state)) {
        return false;
    }

    // Read infrared luminosities, if any
    if (!read_lir(opts, state)) {
        return false;
    }

    // Read the template error function, if any
    if (!read_template_error(opts, state)) {
        return false;
    }

    // Read the continuum indices definitions, if any
    if (!read_continuum_indices(opts, state)) {
        return false;
    }

    // Check and adjust input
    if (!check_input(opts, state)) {
        return false;
    }

    return true;
}
