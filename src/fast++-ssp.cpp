#include "fast++-ssp.hpp"
#include <phypp/utility/string.hpp>
#include <phypp/io/fits.hpp>

bool ssp_bc03::read_ascii(std::string filename) {
    std::string state = "";

    try {
        std::ifstream in(filename);
        in.exceptions(in.failbit);

        state = "read number of time steps";
        uint_t ntime = 0;
        in >> ntime;

        state = "read time steps";
        age.resize(ntime);
        for (uint_t i : range(age)) {
            in >> age[i];
        }

        state = "read IMF and other parameters";
        double ml, mu;
        uint_t iseg;
        in >> ml >> mu >> iseg;
        for (uint_t i = 0; i < iseg; ++i) {
            double xx, lm, um, baux, cn, cc;
            in >> xx >> lm >> um >> baux >> cn >> cc;
        }

        state = "read additional parameters";
        double totm, totn, avs, jo, tau, id, tau1, tau2, tau3, tau4;
        in >> totm >> totn >> avs >> jo >> tau >> id >> tau1 >> tau2 >> tau3 >> tau4;

        char id2;
        in >> id2;
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        std::string id3, iop, stelib;
        std::getline(in, id3);
        std::getline(in, iop);
        std::getline(in, stelib);

        state = "read number of wavelength elements";
        uint_t nlam;
        in >> nlam;

        state = "read wavelength elements";
        lambda.resize(nlam);
        for (uint_t i : range(lambda)) {
            in >> lambda[i];
        }

        state = "read SEDs";
        sed.resize(ntime, nlam);
        for (uint_t it : range(ntime)) {
            uint_t nstep = 0;
            in >> nstep;

            if (nstep != nlam) {
                error("corrupted file, wavelength step mismatch: ", nstep, " vs. ", nlam);
                error("reading time step ", it, " of ", ntime);
                return 1;
            }

            for (uint_t il : range(lambda)) {
                in >> sed.safe(it,il);
            }

            // Read extra information
            uint_t nfunc = 0;
            in >> nfunc;

            for (uint_t i = 0; i < nfunc; ++i) {
                double f;
                in >> f;
            }
        }

        state = "read extras";
        uint_t nextra = 12;
        for (uint_t ie : range(nextra)) {
            uint_t nstep = 0;
            in >> nstep;

            vec1d extra(nstep);
            for (uint_t i : range(nstep)) {
                in >> extra[i];
            }

            if (ie == 1) {
                mass = extra/totm;
            }
        }
    } catch (...) {
        print("");
        error("could not read data in library file '", filename, "'");
        error("could not ", state);
        error("the file is probably corrupted, try re-downloading it");
        print("");
        return false;
    }

    // Write FITS file for faster reading next time
    fits::write_table(filename+".fits", ftable(age, mass, lambda, sed));

    return true;
}

bool ssp_bc03::read_fits(std::string filename, bool noflux) {
    fits::input_table itbl(filename);
    if (!itbl.read_column("age", age)) {
        print("");
        error("could not read column 'AGE' (time) from FITS file '", filename, "'");
        print("");
        return false;
    }
    if (!itbl.read_column("mass", mass)) {
        print("");
        error("could not read column 'MASS' from FITS file '", filename, "'");
        print("");
        return false;
    }
    if (!noflux) {
        if (!itbl.read_column("lambda", lambda)) {
            print("");
            error("could not read column 'LAMBDA' from FITS file '", filename, "'");
            print("");
            return false;
        }
        if (!itbl.read_column("sed", sed)) {
            print("");
            error("could not read column 'SED' from FITS file '", filename, "'");
            print("");
            return false;
        }
    }
    return true;
}

bool ssp_bc03::read(std::string filename, bool noflux) {
    if (file::exists(filename+".ised_ASCII.fits")) {
        return read_fits(filename+".ised_ASCII.fits", noflux);
    } else if (file::exists(filename+".ised_ASCII")) {
        return read_ascii(filename+".ised_ASCII");
    } else {
        print("");
        error("could not find library: '", filename, "'");
        error("expected extensions *.fits or *.ised_ASCII");
        print("");
        return false;
    }
}

