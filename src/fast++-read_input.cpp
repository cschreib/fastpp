#include "fast++.hpp"

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

template <typename T>
bool parse_value_impl(std::string val, vec<1,T>& out) {
    if (val.empty()) return true;
    val = remove_first_last(val, "[]");
    vec1s spl = split(val, ",");
    return count(!from_string(spl, out)) == 0;
}

bool parse_value_impl(const std::string& val, std::string& out) {
    out = remove_first_last(val, "\'\"");
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

bool read_params(options_t& opts, const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open()) {
        error("could not open parameter file '", filename, "'");
        return false;
    }

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

        std::string key = toupper(trim(line.substr(0, eqp)));
        std::string val = trim(line.substr(eqp+1));

        if      (key == "CATALOG")       { if (!parse_value(key, val, opts.catalog))        return false; }
        else if (key == "AB_ZEROPOINT")  { if (!parse_value(key, val, opts.ab_zeropoint))   return false; }
        else if (key == "FILTERS_RES")   { if (!parse_value(key, val, opts.filters_res))    return false; }
        else if (key == "FILTER_FORMAT") { if (!parse_value(key, val, opts.filters_format)) return false; }
        else if (key == "TEMP_ERR_FILE") { if (!parse_value(key, val, opts.temp_err_file))  return false; }
        else if (key == "NAME_ZPHOT")    { if (!parse_value(key, val, opts.name_zphot))     return false; }
        else if (key == "SPECTRUM")      { if (!parse_value(key, val, opts.spectrum))       return false; }
        else if (key == "AUTO_SCALE")    { if (!parse_value(key, val, opts.auto_scale))     return false; }
        else if (key == "OUTPUT_DIR")    { if (!parse_value(key, val, opts.output_dir))     return false; }
        else if (key == "OUTPUT_FILE")   { if (!parse_value(key, val, opts.output_file))    return false; }
        else if (key == "N_SIM")         { if (!parse_value(key, val, opts.n_sim))          return false; }
        else if (key == "C_INTERVAL")    { if (!parse_value(key, val, opts.c_interval))     return false; }
        else if (key == "BEST_FIT")      { if (!parse_value(key, val, opts.best_fit))       return false; }
        else if (key == "LIBRARY_DIR")   { if (!parse_value(key, val, opts.library_dir))    return false; }
        else if (key == "LIBRARY")       { if (!parse_value(key, val, opts.library))        return false; }
        else if (key == "RESOLUTION")    { if (!parse_value(key, val, opts.resolution))     return false; }
        else if (key == "IMF")           { if (!parse_value(key, val, opts.imf))            return false; }
        else if (key == "SFH")           { if (!parse_value(key, val, opts.sfh))            return false; }
        else if (key == "DUST_LAW")      { if (!parse_value(key, val, opts.dust_law))       return false; }
        else if (key == "MY_SFH")        { if (!parse_value(key, val, opts.my_sfh))         return false; }
        else if (key == "LOG_TAU_MIN")   { if (!parse_value(key, val, opts.log_tau_min))    return false; }
        else if (key == "LOG_TAU_MAX")   { if (!parse_value(key, val, opts.log_tau_max))    return false; }
        else if (key == "LOG_TAU_STEP")  { if (!parse_value(key, val, opts.log_tau_step))   return false; }
        else if (key == "LOG_AGE_MIN")   { if (!parse_value(key, val, opts.log_age_min))    return false; }
        else if (key == "LOG_AGE_MAX")   { if (!parse_value(key, val, opts.log_age_max))    return false; }
        else if (key == "LOG_AGE_STEP")  { if (!parse_value(key, val, opts.log_age_step))   return false; }
        else if (key == "NO_MAX_AGE")    { if (!parse_value(key, val, opts.no_max_age))     return false; }
        else if (key == "Z_MIN")         { if (!parse_value(key, val, opts.z_min))          return false; }
        else if (key == "Z_MAX")         { if (!parse_value(key, val, opts.z_max))          return false; }
        else if (key == "Z_STEP")        { if (!parse_value(key, val, opts.z_step))         return false; }
        else if (key == "Z_STEP_TYPE")   { if (!parse_value(key, val, opts.z_step_type))    return false; }
        else if (key == "A_V_MIN")       { if (!parse_value(key, val, opts.a_v_min))        return false; }
        else if (key == "A_V_MAX")       { if (!parse_value(key, val, opts.a_v_max))        return false; }
        else if (key == "A_V_STEP")      { if (!parse_value(key, val, opts.a_v_step))       return false; }
        else if (key == "METAL")         { if (!parse_value(key, val, opts.metal))          return false; }
        else if (key == "H0")            { if (!parse_value(key, val, opts.h0))             return false; }
        else if (key == "OMEGA_M")       { if (!parse_value(key, val, opts.omega_m))        return false; }
        else if (key == "OMEGA_L")       { if (!parse_value(key, val, opts.omega_l))        return false; }
        else if (key == "SAVE_CHI_GRID") { if (!parse_value(key, val, opts.save_chi_grid))  return false; }
        // Not in original FAST
        else if (key == "VERBOSE")       { if (!parse_value(key, val, opts.verbose))        return false; }
        else {
            warning("unknown parameter '", key, "'");
        }
    }

    // Now check for the consistency of the output and make corrections when necessary
    if (opts.spectrum.empty()) {
        opts.auto_scale = false;
    }

    if (opts.metal.empty()) {
        if (opts.library == "bc03" || opts.library == "ma05") {
            opts.metal = {0.02};
        } else if (opts.library == "co11") {
            opts.metal = {0.019};
        } else {
            // Uknown stellar library, simply guess what could be a good default value...
            opts.metal = {0.02};
        }
    }

    // Use the default 'share' directory if nothing is provided
    if (opts.filters_res.empty()) {
        opts.filters_res = std::string(FASTPP_SHARE_DIR)+"/FILTER.RES.latest";
    }

    if (opts.temp_err_file.empty()) {
        opts.temp_err_file = std::string(FASTPP_SHARE_DIR)+"/TEMPLATE_ERROR.fast.v0.2";
    }

    if (opts.library_dir.empty()) {
        opts.library_dir = std::string(FASTPP_SHARE_DIR)+"/Libraries/";
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

bool read_filters(const options_t& opts, input_state_t& state) {
    if (opts.filters_res.empty()) {
        // No filter, maybe just trying to fit a spectrum?
        return true;
    }

    if (opts.verbose) {
        note("reading filters from ", opts.filters_res);
    }

    std::ifstream in(opts.filters_res);
    if (!in.is_open()) {
        error("could not open filter file '", opts.filters_res, "'");
        return false;
    }

    // Read the required filters from the database
    std::string line; uint_t l = 0;
    fast_filter_t filt;
    vec1u idcat;
    bool doread = false;
    uint_t ntotfilt = 0;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty()) continue;

        vec1s spl = split_any_of(line, " \t\n\r");
        if (spl.size() > 3) {
            // Start of a new filter

            // Save the previous one, if any
            if (!filt.wl.empty()) {
                state.filters.push_back(filt);

                // Cleanup for next filter
                filt.wl.clear();
                filt.tr.clear();
                filt.id = npos;
            }

            ++ntotfilt;
            filt.id = ntotfilt;

            // Find if this filter is used in the catalog
            uint_t idfil = where_first(state.no_filt == filt.id);
            if (idfil != npos) {
                // It is there, keep the ID aside for later sorting
                idcat.push_back(idfil);
                doread = true;
            } else {
                // If not, discard its values
                doread = false;
            }

        } else if (doread && spl.size() == 3) {
            // Reading the filter response
            double wl, tr;
            if (!from_string(spl[1], wl) || !from_string(spl[2], tr)) {
                error("could not parse values from line ", l);
                note("reading '", opts.filters_res, "'");
                return false;
            }

            filt.wl.push_back(wl);
            filt.tr.push_back(tr);
        }
    }

    // Save the last filter, if any
    if (!filt.wl.empty()) {
        state.filters.push_back(filt);
    }

    if (opts.verbose) {
        note("found ", ntotfilt, " filters in the database, will use ",
            state.filters.size(), " of them");
    }

    // Sort the filter list as in the catalog
    vec1u sid = sort(idcat);
    state.filters = state.filters[sid];

    // Make sure we are not missing any
    if (state.filters.size() != state.no_filt.size()) {
        vec1u nfound;
        for (uint_t b : state.no_filt) {
            bool found = false;
            for (auto& f : state.filters) {
                if (f.id == b) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                nfound.push_back(b);
            }
        }

        error("filters ", nfound, " are missing from the filter library");
        return false;
    }

    // Normalize filter and compute the central wavelength
    for (auto f : state.filters) {
        f.tr /= integrate(f.wl, f.tr);

        state.lambda.push_back(integrate(f.wl, f.wl*f.tr));
    }

    return true;
}

bool read_fluxes(const options_t& opts, input_state_t& state) {
    std::string catalog_file = opts.catalog+".cat";

    if (opts.verbose) {
        note("reading fluxes from ", catalog_file);
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
            note("using column translation file ", catalog_file);
        }

        vec1s tr_from, tr_to;
        file::read_table(translate_file, file::find_skip(translate_file), tr_from, tr_to);

        vec1u idh, idt;
        match(header, tr_from, idh, idt);
        header_trans[idh] = tr_to[idt];
    }

    header_trans = toupper(header_trans);

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
    uint_t col_tot = where_first(is_any_of(header_trans, "TOT"+strna(state.no_filt)));

    if (col_id == npos) {
        error("missing ID column from photometric catalog");
        return false;
    }

    // Read filters from the database
    if (!read_filters(opts, state)) {
        return false;
    }

    // Now read the catalog itself, only keeping the columns we are interested in
    uint_t l = 0;
    std::ifstream in(catalog_file);
    std::string line;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split_any_of(line, " \t\n\r");

        // Read the ID
        uint_t tid;
        if (!from_string(spl[col_id], tid)) {
            error("could not read ID from line ", l);
            note("must be a positive integer (>= 0), got: '", spl[col_id], "'");
            return false;
        }

        state.id.push_back(tid);

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

            state.zspec.push_back(tz);
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

        // Flag bad values
        vec1u idb = where(err < 0 || !is_finite(flx) || !is_finite(err));
        err[idb] = fnan; flx[idb] = 0;

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

        // Save flux and uncertainties in the input state
        append<0>(state.flux, reform(flx, 1, flx.size()));
        append<0>(state.eflux, reform(err, 1, err.size()));
    }

    if (col_zspec == npos) {
        // If no zspec column, give no zspec to all galaxies
        state.zspec = replicate(fnan, state.id.size());
    }

    // Convert photometry from [catalog unit] to [uJy]
    float abzp = e10(0.4*(23.9 - opts.ab_zeropoint));
    // Convert photometry from fnu [uJy] to flambda [1e-19 x erg/s/cm2/A]
    vec1u idbb = uindgen(state.no_filt.size());
    for (uint_t i : range(state.id)) {
        state.flux(i,idbb) = 1e19*abzp*uJy2cgs(state.lambda*1e-4, state.flux(i,idbb));
        state.eflux(i,idbb) = 1e19*abzp*uJy2cgs(state.lambda*1e-4, state.eflux(i,idbb));
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
        note("reading spectra from ", opts.spectrum);
    }

    if (!file::exists(opts.spectrum)) {
        error("could not open spectral catalog '", opts.spectrum, "'");
        return false;
    }

    // Read the header
    vec1s spec_header;
    if (!read_header(opts.spectrum, spec_header)) {
        return false;
    }

    spec_header = toupper(spec_header);

    // Check we have the right columns there
    uint_t col_wl0 = where_first(spec_header == "WL_LOW");
    uint_t col_wl1 = where_first(spec_header == "WL_UP");
    uint_t col_tr = where_first(spec_header == "TR");
    uint_t col_bin = where_first(spec_header == "BIN");

    // First check if the user is not using the FAST-IDL format
    if ((col_wl0 == npos || col_wl1 == npos) && count(spec_header == "WL") != 0) {
        error("FAST++ requires WL_LOW and WL_UP columns for spectra instead of WL");
        note("WL_LOW and WL_UP define the width of the spectral element");
        return false;
    }

    // Then check for missing wavelength columns
    if (col_wl0 == npos) {
        error("missing lower wavelength column (WL_LOW) in spectral catalog");
        return false;
    }
    if (col_wl1 == npos) {
        error("missing upper wavelength column (WL_UP) in spectral catalog");
        return false;
    }

    vec1u sid;
    vec1u col_flux, col_eflux;
    for (uint_t b : range(spec_header)) {
        vec2s ext = regex_extract(spec_header[b], "F([0-9]+)");
        if (ext.empty()) continue;
        uint_t id;
        from_string(ext[0], id);

        uint_t cid = where_first(state.id == id);
        if (cid == npos) {
            warning("spectrum for source ", id, " has no corresponding photometry "
                "and will be ignored");
            continue;
        }

        uint_t ce = where_first(spec_header == "E"+ext[0]);
        if (ce == npos) {
            warning("spectral flux column ", spec_header[b], " has no uncertainty (E"+ext[0]+") "
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
    std::ifstream in(opts.spectrum);
    std::string line;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split_any_of(line, " \t\n\r");

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

        double wl0;
        if (!from_string(spl[col_wl0], wl0)) {
            error("could not read spectral lower wavelength from line ", l);
            note("must be a floating point number, got: '", spl[col_wl0], "'");
            return false;
        }

        slam0.push_back(wl0);

        double wl1;
        if (!from_string(spl[col_wl1], wl1)) {
            error("could not read spectral upper wavelength from line ", l);
            note("must be a floating point number, got: '", spl[col_wl1], "'");
            return false;
        }

        slam1.push_back(wl1);

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
        err[idb] = fnan; flx[idb] = 0;

        // Add fluxes to the list
        append<1>(sflx, reform(flx, flx.size(), 1));
        append<1>(serr, reform(err, err.size(), 1));
    }

    if (opts.verbose) {
        note("found ", sflx.dims[0], " spectr", (sflx.dims[0] > 1 ? "a" : "um"),
            " with ", sflx.dims[1], " spectral elements", (sflx.dims[0] > 1 ? " each" : ""));
    }

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
            w[where(!is_finite(w))] = 0;
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
        serr[idb] = fnan; sflx[idb] = 0;

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
    for (uint_t b : range(sflx)) {
        fast_filter_t f;
        f.spectral = true;
        f.id = max(state.no_filt)+1 + b;
        double dl = slam1[b] - slam0[b];
        f.wl = {slam0[b]-dl, slam0[b]-1e-4*dl, slam0[b], slam1[b], slam1[b]+1e-4*dl, slam1[b]+dl};
        f.tr = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
        f.tr /= integrate(f.wl, f.tr);

        state.filters.push_back(f);
        state.lambda.push_back(0.5*(slam0[b] + slam1[b]));
    }

    // Merge the spectra into the flux catalog
    vec2f pflx = std::move(state.flux);
    vec2f perr = std::move(state.eflux);
    uint_t nplam = pflx.dims[1];
    uint_t nslam = sflx.dims[1];
    uint_t ngal = pflx.dims[0];

    state.flux = replicate(fnan, ngal, nplam+nslam);
    state.eflux = replicate(fnan, ngal, nplam+nslam);

    for (uint_t i : range(ngal)) {
        state.flux(i,_-(nplam-1)) = pflx(i,_);
        state.eflux(i,_-(nplam-1)) = perr(i,_);
    }

    for (uint_t i : range(sid)) {
        state.flux(sid[i],nplam-_) = sflx(i,_);
        state.eflux(sid[i],nplam-_) = serr(i,_);
    }

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

    header = toupper(header);

    // Check we have all the columns we need
    uint_t col_zphot = where_first(header == toupper(opts.name_zphot));
    uint_t col_zspec = where_first(header == "Z_SPEC");
    uint_t col_id = where_first(header == "ID");

    vec1u col_up, col_low;
    for (uint_t i : range(opts.c_interval)) {
        uint_t cl = where_first(header == "L"+strn(opts.c_interval[i]));
        uint_t cu = where_first(header == "U"+strn(opts.c_interval[i]));
        if (cl == npos) {
            error("could not find redshift confidence interval column L"+strn(opts.c_interval[i]));
            return false;
        } else if (cu == npos) {
            error("could not find redshift confidence interval column U"+strn(opts.c_interval[i]));
            return false;
        }

        col_up.push_back(cu);
        col_low.push_back(cl);
    }

    if (col_zphot == npos) {
        error("missing zphot (", toupper(opts.name_zphot), ") column in photometric redshift file");
        return false;
    }
    if (col_id == npos) {
        error("missing ID column in photometric redshift file");
        return false;
    }

    // Initialize the zphot columns
    state.zphot = replicate(fnan, state.id.size());
    state.zphot_low = replicate(fnan, state.id.size(), opts.c_interval.size());
    state.zphot_up = replicate(fnan, state.id.size(), opts.c_interval.size());

    // Read the catalog
    uint_t l = 0;
    uint_t i = 0;
    std::ifstream in(opts.spectrum);
    std::string line;
    while (std::getline(in, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = split_any_of(line, " \t\n\r");

        uint_t id;
        if (!from_string(spl[col_id], id)) {
            error("could not read source ID from line ", l);
            note("must be a positive integer (>= 0), got: '", spl[col_id], "'");
            return false;
        }

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

        state.zphot[i] = (zp > 0 ? zp : fnan);

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
                error("could not read photometric redshift ", opts.c_interval[c],
                    "th confidence interval from line ", l);
                note("must be a floating point number, got: '", spl[col_low[c]], "'");
                return false;
            }

            state.zphot_low(i,c) = (zp > 0 ? zp : fnan);
        }

        for (uint_t c : range(col_up)) {
            if (!from_string(spl[col_up[c]], zp)) {
                error("could not read photometric redshift ", opts.c_interval[c],
                    "th confidence interval from line ", l);
                note("must be a floating point number, got: '", spl[col_up[c]], "'");
                return false;
            }

            state.zphot_up(i,c) = (zp > 0 ? zp : fnan);
        }

        ++i;
    }

    if (i != state.id.size()) {
        error("photometric redshift and photometry catalogs do not match");
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

    file::read_table(opts.temp_err_file, file::find_skip(opts.temp_err_file),
        state.tplerr_lam, state.tplerr_err
    );

    if (state.tplerr_lam.empty()) {
        error("template error function is empty, something must be wrong in the file");
        return false;
    }

    return true;
}

bool read_input(options_t& opts, input_state_t& state, const std::string& filename) {
    // First read options from the parameter file
    if (!read_params(opts, filename)) {
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

    // Read the template error function, if any
    if (!read_template_error(opts, state)) {
        return false;
    }

    return true;
}