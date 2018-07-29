#include <phypp.hpp>
#include <iomanip>
#include "fast++-ssp.hpp"

int phypp_main(int argc, char* argv[]) {
    if (argc <= 1) {
        print("usage: fast++-sfh2sed sfh=... tobs=... ssp=... out=... "
            "[Av=... z=... H0=... omega_L=... omega_m=...]");
        return 0;
    }

    // Read command line arguments
    std::string sfh_file, library_file, out_file;
    double tobs = dnan;
    double z = dnan;
    double Av = 0.0; // [mag]
    double sfh_step = 1.0; // [Myr]
    double H0 = 70.0;
    double omega_L = 0.7;
    double omega_m = 0.3;
    bool lookback = false;

    read_args(argc, argv, arg_list(name(sfh_file, "sfh"), name(library_file, "ssp"),
        name(out_file, "out"), tobs, z, Av, sfh_step, H0, omega_m, omega_L, lookback));

    // Check arguments
    if (!is_finite(tobs) && !lookback) {
        error("please provide the time of observation in 'tobs=...' (in Myr), "
            "or specify 'lookback' if time is in lookback unit (tobs=0)");
        return 1;
    }
    if (sfh_file.empty()) {
        error("please provide the path to the file containing the SFH in 'sfh=...'");
        return 1;
    }
    if (library_file.empty()) {
        error("please provide the path to the SSP library in 'ssp=...'");
        return 1;
    }
    if (out_file.empty()) {
        error("please provide the path to the output file in 'out=...'");
        return 1;
    }

    cosmo_t cosmo;
    cosmo.H0 = H0;
    cosmo.wL = omega_L;
    cosmo.wm = omega_m;

    // Read SSP library
    ssp_bc03 ssp;
    ssp.read(file::get_directory(library_file)+file::get_basename(library_file));
    ssp.age /= 1e6; // yr to Myr
    double max_age = max(ssp.age);

    // Read SFH
    vec1d input_t, input_sfr;
    if (ends_with(sfh_file, ".fits")) {
        fits::read_table(sfh_file, "t", input_t, "sfr", input_sfr);
    } else {
        ascii::read_table(sfh_file, input_t, input_sfr);
    }

    // Check input
    double min_t = min(input_t), max_t = max(input_t);
    if (max_t > 1e5) {
        error("time axis must be given in Myr, found too high values (", max_t, ")");
        return 1;
    }
    if (min_t < 0) {
        error("time axis contains negative values");
        return 1;
    }
    if (!lookback && tobs > max_t) {
        error("requested observation time (", tobs, ") is larger than maximum time in SFH (",
            max_t, ")");
        return 1;
    }
    if (lookback && min_t > 0) {
        error("when 'lookback' is set, time axis must reach zero (found minimum: ", min_t, ")");
        return 1;
    }
    if ((lookback && max_t > max_age) || (!lookback && tobs - min_t > max_age)) {
        error("SFH goes beyond the maximum age of the library (", max_age/1e3, " Gyr)");
        return 1;
    }
    if (min(input_sfr) < 0) {
        error("input SFR is negative");
        return 1;
    }
    if (input_t.size() > 1 && input_t[1] < input_t[0]) {
        error("time axis must be sorted by increasing value");
        return 1;
    }

    // Resample SFH on grid
    double tt0 = min_t;
    double tt1 = tobs;
    if (lookback) {
        tt0 = 0.0;
        tt1 = max_t;
    }

    vec1d t = rgen_step(tt0, tt1, sfh_step);
    vec1d sfr;

    input_sfr *= 1e6;

    double orig_step = median(input_t[1-_] - input_t[_-(input_t.size()-2)]);
    if (orig_step < sfh_step) {
        // Average
        sfr.resize(t.size());
        for (uint_t i : range(t)) {
            double t0 = max(min_t, t[i] - 0.5*sfh_step);
            double t1 = min(max_t, t[i] + 0.5*sfh_step);
            if (t0 >= t1) break;
            sfr[i] = integrate(input_t, input_sfr, t0, t1)/(t1 - t0);
        }
    } else {
        // Interpolate
        sfr = interpolate(input_sfr, input_t, t);
    }

    // Put in lookback unit
    if (!lookback) {
        t = reverse(tobs - t);
        sfr = reverse(sfr);
    }

    fits::write_table("tmp.fits", ftable(t, sfr));

    // Sum SSPs
    vec1d tpl_flux(ssp.lambda.size());
    ssp.integrate(t, sfr, [&](uint_t it, double formed) {
        tpl_flux += formed*ssp.sed.safe(it,_);
    });

    if (Av > 0) {
        // Apply dust attenuation
        auto calzetti2000 = [](double l) {
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
        };

        for (uint_t il : range(tpl_flux)) {
            tpl_flux.safe[il] *= e10(-0.4*Av*calzetti2000(ssp.lambda.safe[il]));
        }
    }

    if (is_finite(z)) {
        // Apply redshift
        const double dist_Mpc_to_cgs = 3.0856e24; // [cm/Mpc]
        const double lum_sol_to_cgs  = 3.839e33;  // [erg/s/Lsol]
        const double factor = 1e19*lum_sol_to_cgs/sqr(dist_Mpc_to_cgs);
        double lum2fl = factor/(4.0*dpi*(1.0+z)*sqr(astro::lumdist(z, cosmo)));

        auto madau1995 = [](double tz, vec1d lam) {
            // http://adsabs.harvard.edu/abs/1995ApJ...441...18M
            // TODO: check this implementation someday, I suspect this is wrong or
            // very approximate (taken directly from FAST)

            double da; {
                double l0 = 1050.0*(1.0 + tz);
                double l1 = 1170.0*(1.0 + tz);
                uint_t nstep = 100;
                vec1d tl = rgen(l0, l1, nstep);
                vec1d ptau = exp(-3.6e-3*pow(tl/1216.0, 3.46));
                da = total(ptau)*(l1-l0)/nstep/(120.0*(1.0 + tz));
            }

            double db; {
                double l0 = 920.0*(1.0 + tz);
                double l1 = 1015.0*(1.0 + tz);
                uint_t nstep = 100;
                vec1d tl = rgen(l0, l1, nstep);
                vec1d ptau = exp(-1.7e-3*pow(tl/1026.0, 3.46) - 1.2e-3*pow(tl/972.5, 3.46) -
                    9.3e-4*pow(tl/950.0, 3.46));
                db = total(ptau)*(l1-l0)/nstep/(95.0*(1.0 + tz));
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
        };

        ssp.lambda *= (1.0 + z);
        tpl_flux *= lum2fl*madau1995(z, ssp.lambda);
    }

    // Save SED
    std::ofstream fout(out_file);
    if (is_finite(z)) {
        fout << "# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)\n";
    } else {
        fout << "# wl fl (Lsol Angstrom^-1)\n";
    }

    for (uint_t il : range(ssp.lambda)) {
        fout << std::setw(13) << ssp.lambda.safe[il] << std::setw(13) << tpl_flux.safe[il] << "\n";
    }
    fout.close();

    return 0;
}
