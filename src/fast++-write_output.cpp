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
    if (opts.my_sfh.empty()) {
    fout << "# SFH:         " << pretty_sfh(opts.sfh) << std::endl;
    } else {
    fout << "# SFH:         " << "custom SFH: " << opts.my_sfh << std::endl;
    }
    fout << "# Stellar IMF: " << pretty_imf(opts.imf) << std::endl;
    fout << "# Dust law:    " <<
        pretty_dust_law(opts.dust_law, opts.dust_noll_eb, opts.dust_noll_delta) << std::endl;
    fout << "# metallicity: " << collapse(strna(opts.metal), "  ") << std::endl;
    fout << "# log(tau/yr): " <<
        pretty_grid(opts.log_tau_min, opts.log_tau_max, opts.log_tau_step, 2) << std::endl;
    fout << "# log(age/yr): " <<
        pretty_grid(opts.log_age_min, opts.log_age_max, opts.log_age_step, 2) << std::endl;
    fout << "# A_V:         " <<
        pretty_grid(opts.a_v_min, opts.a_v_max, opts.a_v_step, 2) << std::endl;
    fout << "# z:           " <<
        pretty_grid(opts.z_min, opts.z_max, opts.z_step, 4) <<
        " (" << (opts.z_step_type == 0 ? "linear" : "logarithmic") << ")" << std::endl;
    fout << "# Filters:     " << collapse(strna(input.no_filt), "  ") << std::endl;

    std::string additional_abbrev;
    if (opts.output_ldust) {
        additional_abbrev += ", lldust: log[ldust/Lsol]";
    }

    fout << "# ltau: log[tau/yr], lage: log[age/yr], lmass: log[mass/Msol], "
        "lsfr: log[sfr/(Msol/yr)], lssfr: log[ssfr*yr], la2t: log[age/tau]"+additional_abbrev << std::endl;
    fout << "# For sfr=0. lsfr is set to -99" << std::endl;

    vec1s param = {"id"};
    uint_t maxid = 7;
    if (!input.id.empty()) {
        maxid = max(maxid, max(length(input.id)+1));
    }
    vec1u cwidth = {maxid};

    vec1s oparam = {"z", "ltau", "metal", "lage", "Av", "lmass", "lsfr", "lssfr", "la2t"};
    uint_t ocwidth = 10;
    for (uint_t i : range(oparam)) {
        param.push_back(oparam[i]);
        cwidth.push_back(max(ocwidth, param.back().size()+1));
        if (!opts.c_interval.empty()) {
            for (float c : input.conf_interval) {
                float cc = 100*(1-2*c);
                std::string is = (cc < 0.0 ? "u" : "l")+strn(round(abs(cc)));
                param.push_back(is+"_"+oparam[i]);
                cwidth.push_back(max(ocwidth, param.back().size()+1));
            }
        }
    }

    param.push_back("chi2");
    cwidth.push_back(15);

    fout << "#";
    for (uint_t ip : range(param)) {
        fout << align_right(param[ip], cwidth[ip]);
    }
    fout << std::endl;

    // Print data
    for (uint_t is : range(input.id)) {
        fout << " ";
        uint_t c = 0;
        fout << align_right(input.id[is], cwidth[c]); ++c;

        for (uint_t ip : range(output.best_z.dims[1])) {
            fout << align_right(pretty_strn(round(1e4*output.best_z(is,ip))/1e4), cwidth[c]); ++c;
        }
        for (uint_t ip : range(output.best_tau.dims[1])) {
            fout << align_right(pretty_strn(round(1e2*output.best_tau(is,ip))/1e2), cwidth[c]); ++c;
        }
        for (uint_t ip : range(output.best_metal.dims[1])) {
            fout << align_right(pretty_strn(round(1e4*output.best_metal(is,ip))/1e4), cwidth[c]); ++c;
        }
        for (uint_t ip : range(output.best_age.dims[1])) {
            fout << align_right(pretty_strn(round(1e2*output.best_age(is,ip))/1e2), cwidth[c]); ++c;
        }
        for (uint_t ip : range(output.best_av.dims[1])) {
            fout << align_right(pretty_strn(round(1e2*output.best_av(is,ip))/1e2), cwidth[c]); ++c;
        }
        for (uint_t ip : range(output.best_mass.dims[1])) {
            fout << align_right(pretty_strn(round(1e2*log10(output.best_mass(is,ip)))/1e2), cwidth[c]); ++c;
        }
        for (uint_t ip : range(output.best_sfr.dims[1])) {
            fout << align_right(pretty_strn(round(1e2*log10(output.best_sfr(is,ip)))/1e2), cwidth[c]); ++c;
        }
        for (uint_t ip : range(output.best_ssfr.dims[1])) {
            fout << align_right(pretty_strn(round(1e2*log10(output.best_ssfr(is,ip)))/1e2), cwidth[c]); ++c;
        }
        for (uint_t ip : range(output.best_tau.dims[1])) {
            float la2t = output.best_age(is,ip) - output.best_tau(is,ip);
            fout << align_right(pretty_strn(round(1e2*la2t)/1e2), cwidth[c]); ++c;
        }

        uint_t nobs = count(is_finite(input.eflux(is,_)));
        float chi2 = output.best_chi2[is]/max(1u, nobs > gridder.nparam ? nobs - gridder.nparam : 1u);
        fout << align_right(strn_sci(chi2), cwidth[c]);
        fout << "\n";
        ++c;

        if (c != cwidth.size()) {
            error("mismatch in column writing code, please report!");
            return;
        }
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

    vec1f lam, sed, flx;
    auto pg = progress_start(input.id.size());
    for (uint_t is : range(input.id)) {
        // Get model
        gridder.build_template(output.best_model[is], lam, sed, flx);
        sed *= output.best_mass(is,0);
        flx *= output.best_mass(is,0);

        // Save model
        // TODO: use a more efficient serialization code
        std::ofstream fout(odir+opts.catalog+"_"+input.id[is]+".fit");
        fout << "# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
        for (uint_t il : range(lam)) {
            fout << align_right(strn(lam.safe[il]), 13)
                 << align_right(strn(sed.safe[il]), 13)
                 << "\n";
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

void write_output(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
    const output_state_t& output) {

    write_catalog(opts, input, gridder, output);

    if (opts.best_fit) {
        write_best_fits(opts, input, gridder, output);
    }
}
