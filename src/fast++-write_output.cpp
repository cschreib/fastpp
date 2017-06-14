#include "fast++.hpp"

std::string pretty_library(const std::string& lib) {
    if (lib == "bc03") return "Bruzual & Charlot (2003)";
    if (lib == "cb07") return "Charlot & Bruzual (2007)";
    if (lib == "ma05") return "Maraston (2005)         ";
    if (lib == "ma11") return "Maraston (2011)         ";
    if (lib == "co11") return "FSPS (Conroy et al.)    ";
    return lib;
}

std::string pretty_sfh(const std::string& sfh) {
    if (sfh == "exp") return "Exponentially declining SFH: SFR ~ exp(-t/tau)";
    if (sfh == "del") return "Delayed exponential SFH: SFR ~ t exp(-t/tau)  ";
    if (sfh == "tru") return "Truncated SFH: single burst of length tau     ";
    return sfh;
}

std::string pretty_imf(const std::string& imf) {
    if (imf == "ch") return "Chabrier";
    if (imf == "sa") return "Salpeter";
    if (imf == "kr") return "Kroupa  ";
    return imf;
}

std::string pretty_dust_law(const std::string& law, float eb, float delta) {
    if (law == "mw")       return "Milky Way dust attenuation law                    ";
    if (law == "calzetti") return "Calzetti (2000) dust attenuation law              ";
    if (law == "kc")       return "Kriek & Conroy (2013) average dust attenuation law";
    if (law == "noll")     return "Noll et al. dust law, E_B="+strn(eb)+", delta="+strn(delta);
    return law;
}

template<typename T>
std::string pretty_grid(T x0, T x1, T step, uint_t ndecimal) {
    double dd = e10(ndecimal);
    uint_t cwidth = 7;
    return align_left(strn(round(dd*x0)/dd), cwidth)+" - "+
           align_left(strn(round(dd*x1)/dd), cwidth)+" in steps of "+
           strn(round(dd*step)/dd);
}

template<typename T>
std::string pretty_strn(T v) {
    return strn(v);
}

template<typename T>
std::string pretty_strn_float(T f) {
    if (!is_nan(f) && !is_finite(f)) {
        if (f < 0) {
            return "-99";
        } else {
            return "99";
        }
    } else {
        return strn(f);
    }
}

std::string pretty_strn(float f) {
    return pretty_strn_float(f);
}

std::string pretty_strn(double f) {
    return pretty_strn_float(f);
}

void write_catalog(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
    const output_state_t& output) {

    if (opts.verbose) note("saving catalog");

    std::ofstream fout(opts.output_dir+opts.output_file+".fout");

    // Print header
    fout << "# FAST++ version: " << fastpp_version << std::endl;
    fout << "# Photometric catalog file: " << opts.catalog << ".cat" << std::endl;
    if (!input.zphot.empty()) {
    fout << "# Photometric redshift file: " << opts.catalog << ".zout" << std::endl;
    }
    if (!opts.spectrum.empty()) {
    fout << "# Spectrum file: " << opts.spectrum << std::endl;
    }
    if (!opts.temp_err_file.empty()) {
    fout << "# Template error function: " << opts.temp_err_file << std::endl;
    }
    fout << "# AB ZP:       " << opts.ab_zeropoint << std::endl;
    fout << "# Library:     " << pretty_library(opts.library) << std::endl;
    if (opts.sfh == sfh_type::single) {
    fout << "# SFH:         " << opts.my_sfh;
    } else {
    fout << "# SFH:         " << pretty_sfh(opts.name_sfh);
    }
    if (opts.sfr_avg > 0) {
        fout << " (<SFR> over " << opts.sfr_avg/1e6 << " Myr)" << std::endl;
    } else {
        fout << " (inst. SFR)" << std::endl;
    }
    fout << "# Stellar IMF: " << pretty_imf(opts.name_imf) << std::endl;
    fout << "# Dust law:    " <<
        pretty_dust_law(opts.dust_law, opts.dust_noll_eb, opts.dust_noll_delta) << std::endl;
    fout << "# metallicity: " << collapse(strna(opts.metal), "  ") << std::endl;
    if (opts.sfh == sfh_type::gridded) {
    fout << "# log(tau/yr): " <<
        pretty_grid(opts.log_tau_min, opts.log_tau_max, opts.log_tau_step, 2) << std::endl;
    }
    fout << "# log(age/yr): " <<
        pretty_grid(opts.log_age_min, opts.log_age_max, opts.log_age_step, 2) << std::endl;
    fout << "# A_V:         " <<
        pretty_grid(opts.a_v_min, opts.a_v_max, opts.a_v_step, 2) << std::endl;
    fout << "# z:           " <<
        pretty_grid(opts.z_min, opts.z_max, opts.z_step, 4) <<
        " (" << (opts.z_step_type == 0 ? "linear" : "logarithmic") << ")" << std::endl;
    fout << "# Filters:     " << collapse(strna(input.no_filt), "  ") << std::endl;

    vec1u idp(opts.output_columns.size());
    for (uint_t ic : range(opts.output_columns)) {
        idp[ic] = where_first(tolower(output.param_names) == tolower(opts.output_columns[ic]));
    }

    std::string abbrev = "# ";
    for (uint_t ic : range(idp)) {
        if (idp[ic] == npos || output.param_descriptions[idp[ic]].empty()) continue;
        if (abbrev != "# ") abbrev += ", ";
        abbrev += output.param_names[idp[ic]]+": "+output.param_descriptions[idp[ic]];
    }

    fout << abbrev << std::endl;
    fout << "# For value=0. log[value] is set to -99" << std::endl;

    vec1s param;
    vec1u cwidth;
    vec1u iparam;
    vec1u iconf;

    for (uint_t ic : range(idp)) {
        std::string cname = opts.output_columns[ic];
        param.push_back(cname);
        iparam.push_back(idp[ic]);
        iconf.push_back(0);

        if (cname == "id") {
            uint_t maxid = 7;
            if (!input.id.empty()) {
                maxid = max(maxid, max(length(input.id)+1));
            }
            cwidth.push_back(maxid);
        } else if (cname == "chi2") {
            cwidth.push_back(15);
        } else {
            uint_t ocwidth = 10;
            // Make sure we use the right format
            cname = output.param_names[idp[ic]];
            param.back() = cname;
            cwidth.push_back(max(ocwidth, cname.size()+1));

            if (!opts.c_interval.empty()) {
                for (uint_t ip : range(input.conf_interval)) {
                    float cc = 100*(1-2*input.conf_interval[ip]);
                    std::string is = (cc < 0.0 ? "u" : "l")+strn(round(abs(cc)));
                    param.push_back(is+"_"+cname);
                    cwidth.push_back(max(ocwidth, param.back().size()+1));
                    iparam.push_back(idp[ic]);
                    iconf.push_back(1+ip);
                }
            }
        }
    }

    fout << "#";
    for (uint_t ic : range(param)) {
        fout << align_right(param[ic], cwidth[ic]);
    }
    fout << std::endl;

    // Print data
    for (uint_t is : range(input.id)) {
        fout << " ";

        for (uint_t ic : range(param)) {
            if (iparam[ic] == npos) {
                std::string cname = param[ic];
                if (cname == "id") {
                    fout << align_right(input.id[is], cwidth[ic]);
                } else if (cname == "chi2") {
                    uint_t nobs = count(is_finite(input.eflux(is,_)));
                    if (!input.lir.empty() && is_finite(input.lir[is])) ++nobs;

                    uint_t ndof = (nobs > gridder.nfreeparam ? nobs - gridder.nfreeparam : 1u);
                    float chi2 = output.best_chi2[is]/ndof;
                    fout << align_right(strn_sci(chi2), cwidth[ic]);
                }
            } else {
                float value = output.best_params.safe(is,iparam[ic],iconf[ic]);
                if (output.param_log.safe[iparam[ic]]) {
                    value = log10(value);
                }

                float precision = output.param_precision.safe[iparam[ic]];
                value = round(value/precision)*precision;
                fout << align_right(pretty_strn(value), cwidth[ic]);
            }
        }

        fout << "\n";
    }
}

void write_best_fits(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
    const output_state_t& output) {

    if (opts.verbose) note("saving best fits");

    std::string odir = opts.output_dir+"best_fits/";
    if (!file::mkdir(odir)) {
        warning("could not save best fit SEDs");
        warning("the output directory '", odir, "' could not be created");
        return;
    }

    vec1f lam, sed, sed_nodust, flx;
    auto pg = progress_start(input.id.size());
    for (uint_t is : range(input.id)) {
        float mass = output.best_params(is,gridder.nparam+prop_id::mass,0);

        // Get model
        if (opts.intrinsic_best_fit) {
            gridder.build_template_nodust(output.best_model[is], lam, sed_nodust, flx);
            sed_nodust *= mass;
        }

        gridder.build_template(output.best_model[is], lam, sed, flx);
        sed *= mass;
        flx *= mass;

        // Save model
        // TODO: use a more efficient serialization code
        std::ofstream fout(odir+opts.catalog+"_"+input.id[is]+".fit");
        if (opts.intrinsic_best_fit) {
            fout << "# wl fl fl_nodust (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        } else {
            fout << "# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        }

        for (uint_t il : range(lam)) {
            fout << align_right(strn(lam.safe[il]), 13)
                 << align_right(strn(sed.safe[il]), 13);

            if (opts.intrinsic_best_fit) {
            fout << align_right(strn(sed_nodust.safe[il]), 13);
            }

            fout << "\n";
        }
        fout.close();

        // Save fluxes
        fout.open(odir+opts.catalog+"_"+input.id[is]+".input_res.fit");
        fout << "# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        for (uint_t il : range(input.lambda)) {
            fout << align_right(strn(float(input.lambda[il])), 13)
                 << align_right(strn(flx[il]), 13)
                 << align_right(strn(input.flux(is,il)), 13)
                 << align_right(strn(input.eflux(is,il)), 13) << "\n";
        }
        fout.close();

        if (opts.verbose) progress(pg, 13);
    }
}

void write_sfhs(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
    const output_state_t& output) {

    if (opts.verbose) note("saving star formation histories");

    std::string odir = opts.output_dir+"best_fits/";
    if (!file::mkdir(odir)) {
        warning("could not save best fit SFHs");
        warning("the output directory '", odir, "' could not be created");
        return;
    }

    std::string quantity = "SFR";
    std::string unit = "Msol/yr";
    if (opts.sfh_output == "mass") {
        quantity = "Mstar";
        unit = "Msol";
    }

    std::string header = "# t "+quantity+"(t) ("+unit+")";
    if (!opts.c_interval.empty()) {
        header += " med_"+quantity;
        for (float c : input.conf_interval) {
            float cc = 100*(1-2*c);
            std::string is = (cc < 0.0 ? "u" : "l")+strn(round(abs(cc)));
            header += " "+is+"_"+quantity;
        }
    }

    header += "\n";

    auto pg = progress_start(input.id.size());
    for (uint_t is : range(input.id)) {
        vec1u idm = gridder.grid_ids(output.best_model[is]);
        vec1f t = rgen_step(1e6, e10(gridder.auniv[idm[grid_id::z]]), opts.sfh_step);

        // Get best fit SFH
        vec2f sfh(t.size(), 1+input.conf_interval.size()+(input.conf_interval.empty() ? 0 : 1)); {
            vec1f tsfh;
            if (!gridder.get_sfh(output.best_model[is], t, tsfh)) {
                return;
            }

            float mass = output.best_params(is,gridder.nparam+prop_id::mass,0);
            sfh(_,0) = tsfh*mass;
        }

        // Get confidence intervals
        if (!opts.c_interval.empty()) {
            vec2f sim_sfh(t.size(), opts.n_sim);
            for (uint_t ir : range(opts.n_sim)) {
                vec1f tsfh;
                if (!gridder.get_sfh(output.mc_best_model(is,ir), t, tsfh)) {
                    return;
                }

                float mass = output.mc_best_props(is,prop_id::mass,ir);
                sim_sfh(_,ir) = tsfh*mass;
            }

            for (uint_t it : range(t)) {
                vec1f tsfh = sim_sfh(it,_);
                sfh(it,1) = inplace_median(tsfh);
                for (uint_t ic : range(input.conf_interval)) {
                    sfh(it,2+ic) = inplace_percentile(tsfh, input.conf_interval[ic]);
                }
            }
        }

        // Save SFH
        // TODO: use a more efficient serialization code
        std::ofstream fout(odir+opts.catalog+"_"+input.id[is]+".sfh");
        fout << header;

        for (uint_t it : range(t)) {
            fout << align_right(strn_sci(t.safe[it]), 13);

            for (uint_t ic : range(sfh.dims[1])) {
                fout << align_right(strn_sci(sfh.safe(it,ic)), 13);
            }

            fout << "\n";
        }

        fout.close();

        if (opts.verbose) progress(pg, 13);
    }
}

void write_output(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
    const output_state_t& output) {

    write_catalog(opts, input, gridder, output);

    if (opts.best_fit) {
        write_best_fits(opts, input, gridder, output);
    }

    if (opts.best_sfhs) {
        write_sfhs(opts, input, gridder, output);
    }
}
