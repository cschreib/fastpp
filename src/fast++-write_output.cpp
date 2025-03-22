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
    if (law == "noll")     return "Noll et al. dust law, E_B="+to_string(eb)+", delta="+to_string(delta);
    return law;
}

template<typename T>
std::string pretty_grid(T x0, T x1, T step, uint_t ndecimal) {
    if (step != 0) {
        uint_t ndecimal_step = max(0u, -floor(log10(step)));
        if (ndecimal_step > ndecimal) {
            ndecimal = ndecimal_step;
        }
    }
    double dd = e10(ndecimal);
    uint_t cwidth = 7;
    return align_left(to_string(round(dd*x0)/dd), cwidth)+" - "+
           align_left(to_string(round(dd*x1)/dd), cwidth)+" in steps of "+
           to_string(round(dd*step)/dd);
}

void write_catalog(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
    const output_state_t& output) {

    std::string output_path = opts.output_dir+opts.output_file+".fout";
    if (opts.verbose) note("saving catalog at ", output_path);

    std::ofstream fout(output_path);

    // Print header
    if (std::string(fastpp_git_hash).empty()) {
        fout << "# FAST++ version: " << fastpp_version << std::endl;
    } else {
        fout << "# FAST++ version: " << fastpp_version << "-" << fastpp_git_hash << std::endl;
    }
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
    if (!opts.temp_err_spec_file.empty()) {
    fout << "# Template spectrum error function: " << opts.temp_err_spec_file << std::endl;
    }
    fout << "# AB ZP:       " << opts.ab_zeropoint << std::endl;
    fout << "# Library:     " << pretty_library(opts.library) << std::endl;
    switch (opts.sfh) {
    case sfh_type::gridded:
        fout << "# SFH:         " << pretty_sfh(opts.name_sfh); break;
    case sfh_type::custom:
        fout << "# SFH:         " << opts.custom_sfh; break;
    case sfh_type::single:
        fout << "# SFH:         " << opts.my_sfh; break;
    }
    if (opts.sfr_avg > 0) {
        fout << " (<SFR> over " << opts.sfr_avg/1e6 << " Myr)" << std::endl;
    } else {
        fout << " (inst. SFR)" << std::endl;
    }
    fout << "# Stellar IMF: " << pretty_imf(opts.name_imf) << std::endl;
    fout << "# Dust law:    " <<
        pretty_dust_law(opts.dust_law, opts.dust_noll_eb, opts.dust_noll_delta);
    if (opts.differential_a_v) {
        fout << " (with differential attenuation)" << std::endl;
    }
    fout << std::endl;
    fout << "# metallicity: " << collapse(to_string_vector(opts.metal), "  ") << std::endl;
    if (opts.sfh == sfh_type::gridded) {
        fout << "# log(tau/yr): " <<
            pretty_grid(opts.log_tau_min, opts.log_tau_max, opts.log_tau_step, 2) << std::endl;
    } else if (opts.sfh == sfh_type::custom) {
        for (uint_t ip : range(opts.custom_params.size())) {
            fout << "# " << align_left(opts.custom_params[ip]+": ", 13) <<
                pretty_grid(opts.custom_params_min[ip], opts.custom_params_max[ip],
                    opts.custom_params_step[ip], -round(log10(output.param_precision[ip])))
                     << std::endl;
        }
    }
    fout << "# log(age/yr): " <<
        pretty_grid(opts.log_age_min, opts.log_age_max, opts.log_age_step, 2) << std::endl;
    fout << "# A_V:         " <<
        pretty_grid(opts.a_v_min, opts.a_v_max, opts.a_v_step, 2) << std::endl;
    if (opts.differential_a_v) {
        fout << "# A_V_BC:      " <<
            pretty_grid(opts.a_v_bc_min, opts.a_v_bc_max, opts.a_v_bc_step, 2) << std::endl;
    }
    fout << "# z:           " <<
        pretty_grid(opts.z_min, opts.z_max, opts.z_step, 4) <<
        " (" << (opts.z_step_type == 0 ? "linear" : "logarithmic") << ")" << std::endl;
    fout << "# Filters:     " << collapse(to_string_vector(input.no_filt), "  ") << std::endl;

    vec1u idp(opts.output_columns.size());
    for (uint_t ic : range(opts.output_columns)) {
        idp[ic] = where_first(to_lower(output.param_names) == to_lower(opts.output_columns[ic]));
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
        } else if (cname == "nmodel") {
            cwidth.push_back(15);
        } else {
            uint_t ocwidth = opts.output_precision > 0 && opts.output_precision < 1e-6 ?
                4 - round(log10(opts.output_precision)) : 10;

            // Make sure we use the right format
            cname = output.param_names[idp[ic]];
            param.back() = cname;
            cwidth.push_back(max(ocwidth, cname.size()+1));

            if (!opts.c_interval.empty()) {
                for (uint_t ip : range(input.conf_interval)) {
                    float cc = 100*(1-2*input.conf_interval[ip]);
                    std::string is = (cc < 0.0 ? "u" : "l")+to_string(round(abs(cc)));
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

    if (opts.output_precision > 0) {
        fout << std::setprecision(-log10(opts.output_precision));
    }

    // Print data
    for (uint_t is : range(input.id)) {
        fout << " ";

        for (uint_t ic : range(param)) {
            if (iparam[ic] == npos) {
                std::string cname = param[ic];
                if (cname == "id") {
                    fout << std::setw(cwidth[ic]) << input.id[is];
                } else if (cname == "chi2") {
                    uint_t nobs = count(is_finite(input.eflux(is,_)));
                    if (!input.lir.empty() && is_finite(input.lir[is])) ++nobs;

                    uint_t ndof = (nobs > gridder.nfreeparam ? nobs - gridder.nfreeparam : 1u);
                    float chi2 = output.best_chi2[is]/ndof;
                    fout << std::setw(cwidth[ic]) << std::scientific << chi2;
                } else if (cname == "nmodel") {
                    fout << std::setw(cwidth[ic]) << output.num_models[is];
                }
            } else {
                float value = output.best_params.safe(is,iparam[ic],iconf[ic]);
                if (output.param_log.safe[iparam[ic]] == log_style::decimal) {
                    value = log10(value);
                } else if (output.param_log.safe[iparam[ic]] == log_style::abmag) {
                    value = -2.5*log10(value) + 23.9;
                }

                float precision = opts.output_precision > 0 ?
                    opts.output_precision : output.param_precision.safe[iparam[ic]];

                value = round(value/precision)*precision;
                if (!is_nan(value) && !is_finite(value)) {
                    if (value < 0) {
                        value = -99.0f;
                    } else {
                        value = 99.0f;
                    }
                }

                fout.unsetf(std::ios_base::floatfield);
                fout << std::setw(cwidth[ic]) << value;
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

    vec1f lam, sed, sed_nodust, sed_lsf, flx;
    auto pg = progress_start(input.id.size());
    for (uint_t is : range(input.id)) {
        std::string output_file_tpl = odir+file::get_basename(opts.catalog)+"_"+input.id[is]+".fit";
        std::string output_file_flx = odir+file::get_basename(opts.catalog)+"_"+input.id[is]+".input_res.fit";

        if (!input.good[is]) {
            file::remove(output_file_tpl);
            file::remove(output_file_flx);
            if (opts.verbose) progress(pg, 13);
            continue;
        }

        float scale = finf;
        uint_t model = npos;

        if (!opts.best_from_sim) {
            // Get the model with smallest chi2
            scale = output.best_params(is,gridder.nparam+prop_id::scale,0);
            model = output.best_model[is];
        } else {
            // Get the "median model" among MC simulations
            // NB: the median is computed on the multidimensional space of the fitting grid,
            // yet there is no unique definition of the median in N-dimensional space. Here
            // we pick the model with the smallest distance to the other models
            double best_dist = dinf;
            for (uint_t imc : range(opts.n_sim)) {
                vec1i idm = vec1i{gridder.grid_ids(output.mc_best_model(is,imc))};
                double dist = 0.0;
                for (uint_t imc2 : range(opts.n_sim)) {
                    if (imc2 == imc) continue;

                    vec1i idm2 = vec1i{gridder.grid_ids(output.mc_best_model(is,imc2))};
                    dist += sqrt(total(sqr(idm - idm2)));
                }

                if (dist < best_dist) {
                    best_dist = dist;
                    scale = output.mc_best_props(is,prop_id::scale,imc);
                    model = output.mc_best_model(is,imc);
                }
            }
        }

        // Get model
        if (opts.intrinsic_best_fit) {
            if (!gridder.build_template_nodust(model, lam, sed_nodust, flx)) {
                return;
            }

            sed_nodust *= scale;
        }

        if (!gridder.build_template(model, lam, sed, flx)) {
            return;
        }

        sed *= scale;
        flx *= scale;

        if (opts.lsf_best_fit) {
            sed_lsf = gridder.apply_lsf(lam, sed);
        }

        // Save model
        std::ofstream fout(output_file_tpl);
        fout << "# wl fl";
        if (opts.intrinsic_best_fit) {
            fout << " fl_nodust";
        }
        if (opts.lsf_best_fit) {
            fout << " fl_lsf";
        }
        fout << " (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";

        for (uint_t il : range(lam)) {
            fout << std::setw(13) << lam.safe[il] << std::setw(13) << sed.safe[il];

            if (opts.intrinsic_best_fit) {
                fout << std::setw(13) << sed_nodust.safe[il];
            }
            if (opts.lsf_best_fit) {
                fout << std::setw(13) << sed_lsf.safe[il];
            }

            fout << "\n";
        }
        fout.close();

        // Save fluxes
        fout.open(output_file_flx);
        fout << "# wl fl_model fl_obs unc_obs (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        for (uint_t il : range(input.lambda)) {
            fout << std::setw(13) << float(input.lambda[il])
                 << std::setw(13) << flx[il]
                 << std::setw(13) << input.flux(is,il)
                 << std::setw(13) << input.eflux(is,il) << "\n";
        }
        fout.close();

        if (opts.verbose) progress(pg, 13);
    }
}

void write_sfhs(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
    const fitter_t& fitter, const output_state_t& output) {

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

    uint_t nconf = input.conf_interval.size();
    bool has_median = false;
    if (!input.conf_interval.empty() && !opts.interval_from_chi2) {
        // Also add the median SFH for MC simulations
        ++nconf;
        has_median = true;
    }

    std::string header = "# t "+quantity+"(t) ("+unit+")";
    if (!opts.c_interval.empty()) {
        if (has_median) {
            header += " med_"+quantity;
        }

        for (float c : input.conf_interval) {
            float cc = 100*(1-2*c);
            std::string is = (cc < 0.0 ? "u" : "l")+to_string(round(abs(cc)));
            header += " "+is+"_"+quantity;
        }
    }

    header += "\n";

    auto pg = progress_start(input.id.size());
    for (uint_t is : range(input.id)) {
        std::string output_file = odir+opts.catalog+"_"+input.id[is]+".sfh";

        if (!input.good[is]) {
            file::remove(output_file);
            if (opts.verbose) progress(pg, 13);
            continue;
        }

        vec1u idm = gridder.grid_ids(output.best_model[is]);
        vec1d t = rgen_step(1e6, e10(gridder.auniv[idm[grid_id::z]]), opts.sfh_output_step);
        vec2f sfh = replicate(fnan, t.size(), 1+nconf);

        // Best fit
        {
            vec1d tsfh;
            if (!gridder.get_sfh(output.best_model[is], t, tsfh)) {
                return;
            }

            float mass = output.best_params(is,gridder.nparam+prop_id::mass,0);
            sfh(_,0) = tsfh*mass;
        }

        // Get confidence intervals
        if (!opts.c_interval.empty()) {
            if (opts.n_sim != 0) {
                // From Monte Carlo simulations
                vec2f sim_sfh(t.size(), opts.n_sim);
                for (uint_t ir : range(opts.n_sim)) {
                    vec1d tsfh;
                    if (!gridder.get_sfh(output.mc_best_model(is,ir), t, tsfh)) {
                        return;
                    }

                    float mass = output.mc_best_props(is,prop_id::mass,ir);
                    sim_sfh(_,ir) = tsfh*mass;
                }

                for (uint_t it : range(t)) {
                    vec1f tsfh = sim_sfh.safe(it,_);
                    sfh.safe(it,1) = inplace_median(tsfh);
                    for (uint_t ic : range(input.conf_interval)) {
                        sfh.safe(it,2+ic) = inplace_percentile(tsfh, input.conf_interval[ic]);
                    }
                }
            } else {
                // From chi2 grid

                fitter.iterate_best_chi2(is, [&](uint_t id, const vec1f& p, float chi2) {
                    if (chi2 - output.best_chi2.safe[is] > input.delta_chi2.back()) return;

                    vec1d tsfh;
                    if (!gridder.get_sfh(id, t, tsfh)) {
                        return;
                    }

                    float mass = p[prop_id::mass];
                    tsfh *= mass;

                    for (uint_t ic : range(input.conf_interval)) {
                        if (chi2 - output.best_chi2.safe[is] < abs(input.delta_chi2[ic])) {
                            for (uint_t it : range(t)) {
                                double v = tsfh.safe[it];
                                float& saved = sfh.safe(it,1+ic);
                                if (input.delta_chi2[ic] < 0.0) {
                                    if (saved > v || !is_finite(saved)) {
                                        saved = v;
                                    }
                                } else {
                                    if (saved < v || !is_finite(saved)) {
                                        saved = v;
                                    }
                                }
                            }
                        }
                    }
                });
            }
        }

        // Save SFH
        std::ofstream fout(output_file);
        fout << header;

        for (uint_t it : range(t)) {
            fout << std::setw(13) << t.safe[it];

            for (uint_t ic : range(sfh.dims[1])) {
                fout << std::setw(13) << sfh.safe(it,ic);
            }

            fout << "\n";
        }

        fout.close();

        if (opts.verbose) progress(pg, 13);
    }
}

void write_output(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
    const fitter_t& fitter, const output_state_t& output) {

    write_catalog(opts, input, gridder, output);

    if (opts.best_fit) {
        write_best_fits(opts, input, gridder, output);
    }

    if (opts.best_sfhs) {
        write_sfhs(opts, input, gridder, fitter, output);
    }

    if (opts.verbose) {
        note("done writing outputs");
    }
}
