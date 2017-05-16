#include "fast++.hpp"
#include <phypp/utility/thread.hpp>

void gridder_t::cache_manager_t::write_model(const model_t& model) {
    if (!cache_file.is_open()) return;

    // TODO: Consider doing this in a worker thread if it slows down execution
    file::write_as<std::uint32_t>(cache_file, model.igrid);
    file::write(cache_file, model.mass);
    file::write(cache_file, model.sfr);
    file::write(cache_file, model.flux);

    if (!cache_file) {
        warning("could not write to cache file anymore");
        warning("in case you ran out of disk space, the cache file has been removed");
        cache_file.close();
        file::remove(cache_filename);
    }
}

bool gridder_t::cache_manager_t::read_model(model_t& model) {
    if (!cache_file.is_open()) return false;

    file::read_as<std::uint32_t>(cache_file, model.igrid);
    file::read(cache_file, model.mass);
    file::read(cache_file, model.sfr);
    file::read(cache_file, model.flux);

    if (!cache_file) {
        cache_file.close();
        return false;
    }

    return true;
}

gridder_t::gridder_t(const options_t& opt, const input_state_t& inp, output_state_t& out) :
    opts(opt), input(inp), output(out) {

    // First build a new array from the provided parameters
    if (opts.verbose) note("define redshift grid...");
    if (opts.z_step == 0.0 || abs(opts.z_max - opts.z_min) < opts.z_step) {
        // A single redshift
        output.z = {opts.z_min};
    } else {
        if (opts.z_step_type == 0) {
            // "dz = cte" grid
            output.z = rgen_step(opts.z_min, opts.z_max, opts.z_step);
        } else {
            // "dz ~ (1+z)" grid
            output.z = e10(rgen_step(
                log10(1.0 + opts.z_min), log10(1.0 + opts.z_max), opts.z_step
            )) - 1.0;
        }
    }

    // If we have fewer galaxies to fit than the size of the above grid,
    // we use the individual zspec/zphot as the grid.
    if (opts.spectrum.empty() && opts.n_sim == 0 && !input.zphot.empty() &&
        input.zphot.size() < output.z.size()) {

        // First compile valid zphot & zspecs
        vec1f cz = input.zspec;
        vec1u idzp = where(!is_finite(cz));
        cz[idzp] = input.zphot[idzp];

        // ... and only keep valid and unique values
        cz = cz[where(is_finite(cz))];
        cz = unique_values(cz, sort(cz));

        output.z = max(cz, 0.00001);
    }

    // Pre-compute distances
    lum2fl = 4018.5161/(4.0*dpi*(1.0+output.z)*
        sqr(astro::lumdist(output.z, opts.cosmo)));

    // Pre-compute age of Universe
    auniv = log10(lookback_time(1000.0, opts.cosmo) - lookback_time(output.z, opts.cosmo)) + 9;

    // Grid the rest of the parameter space
    if (opts.verbose) note("define parameter grid...");
    output.metal = opts.metal;
    output.tau   = rgen_step(opts.log_tau_min, opts.log_tau_max, opts.log_tau_step);
    output.age   = rgen_step(opts.log_age_min, opts.log_age_max, opts.log_age_step);
    output.av    = rgen_step(opts.a_v_min,     opts.a_v_max,     opts.a_v_step);

    dims = {output.metal.size(), output.tau.size(), output.age.size(),
        output.av.size(), output.z.size()};

    uint_t nmodel = dims[0]*dims[1]*dims[2]*dims[3]*dims[4];

    nparam = 0; // TODO: change, this should be 1 (normalization is a degree of freedom)
    for (uint_t i : range(dims)) {
        if (dims[i] > 1) nparam += 1;
    }

    if (opts.verbose) {
        note("fitting a grid of ", nmodel,
            " templates (nmetal=", output.metal.size(), ",ntau=", output.tau.size(),
            ",nage=", output.age.size(), ",nav=", output.av.size(), ",nz=", output.z.size(), ")");
    }

    output.best_mass = output.best_sfr = output.best_z = output.best_metal =
        output.best_av = output.best_age = output.best_tau = output.best_ssfr =
        replicate(fnan, input.id.size(), 1+input.conf_interval.size());

    output.best_chi2 = replicate(finf, input.id.size());
    output.best_model = replicate(npos, input.id.size());

    if (opts.n_sim > 0) {
        out.mc_best_chi2 = replicate(finf, input.id.size(), opts.n_sim);
        out.mc_best_mass = out.mc_best_sfr = replicate(fnan, input.id.size(), opts.n_sim);
        out.mc_best_model = replicate(npos, input.id.size(), opts.n_sim);
    }

    // Base library properties
    cache.cache_filename = opts.output_dir+opts.library+"_"+opts.resolution+"_"+opts.imf+
        "_"+opts.sfh+"_"+opts.dust_law+"_";
    // Grid parameters
    cache.cache_filename += hash(output.z, output.metal, output.av, output.age, output.tau,
        input.lambda, opts.dust_noll_eb, opts.dust_noll_delta, opts.cosmo.H0, opts.cosmo.wm,
        opts.cosmo.wL)+".grid";

    if (opts.verbose) {
        note("cache file is '", cache.cache_filename, "'");
    }

    // Open grid file
    if (file::exists(cache.cache_filename)) {
        if (opts.verbose) note("checking cache integrity...");
        cache.cache_file.open(cache.cache_filename, std::ios::binary | std::ios::in);

        cache.cache_file.seekg(0, std::ios_base::end);
        uint_t size = cache.cache_file.tellg();
        uint_t size_expected = (sizeof(std::uint32_t)+sizeof(float)*(2+input.lambda.size()))*nmodel;

        if (size != size_expected) {
            warning("cache file is corrupted or invalid, will overwrite it");
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
}

struct file_wrapper {
    std::ifstream in;
    bool doswap = false;

    explicit file_wrapper(std::string filename) : in(filename, std::ios::binary) {
        in.exceptions(in.failbit);
    }

    static void swap_endian(char& t) {}

    template<typename T>
    static void swap_endian(T& t) {
        static_assert(std::is_fundamental<T>::value, "cannot swap endian of complex type");

        union {
            T u;
            unsigned char u8[sizeof(T)];
        } source, dest;

        source.u = t;

        for (size_t k : range(sizeof(T))) {
            dest.u8[k] = source.u8[sizeof(T)-k-1];
        }

        t = dest.u;
    }

    template<std::size_t D, typename T>
    static void swap_endian(vec<D,char>& v) {}

    template<std::size_t D, typename T>
    static void swap_endian(vec<D,T>& v) {
        for (auto& val : v) {
            swap_endian(val);
        }
    }

    template<typename T>
    void read(T& val) {
        in.read(reinterpret_cast<char*>(&val), sizeof(T));

        if (doswap) {
            swap_endian(val);
        }
    }

    template<std::size_t D, typename T>
    void read(vec<D,T>& val) {
        in.read(reinterpret_cast<char*>(val.data.data()), sizeof(T)*val.size());

        if (doswap) {
            swap_endian(val);
        }
    }

    template<typename T>
    void read(T* val, uint_t n) {
        in.read(reinterpret_cast<char*>(val), sizeof(T)*n);

        if (doswap) {
            for (uint_t i : range(n)) {
                swap_endian(val[i]);
            }
        }
    }

    template<typename T>
    void seekg(T i) {
        in.seekg(i);
    }

    template<typename T, typename U>
    void seekg(T i, U u) {
        in.seekg(i, u);
    }
};

struct galaxev_ised {
    vec1f age, sfr, mass;
    vec1f lambda;
    vec2f fluxes;

private :
    vec2f extras;

public :
    bool read(std::string filename) {
        std::string state;
        try {
            file_wrapper lib(filename);

            // The file might have been written with a different endianess...
            // As in FAST, we try to read the number of time steps and see if
            // it makes sense, if not we switch endianess.

            {
                // Ignore the first four bytes (FORTRAN convention)
                state = "read first bytes";
                lib.seekg(4);

                state = "determine endianess";
                std::int32_t ntime = 0;
                lib.read(ntime);

                // Check if the number of steps makes sense
                lib.doswap = ntime < 0 || ntime > 10000;
                if (lib.doswap) {
                    // It doesn't, try again swapping the bytes
                    lib.seekg(4);
                    lib.read(ntime);
                }

                // Read time steps
                state = "read time steps";
                age.resize(ntime);
                lib.read(age);
            }

            // Read through info section (variable total size)
            // Get the number of extras data, in passing
            {
                // Skip IMF boundaries
                state = "read IMF";
                lib.seekg(2*4, std::ios::cur);
                // Skip IMF segments
                std::int32_t imf_seg = 0;
                lib.read(imf_seg);
                lib.seekg(6*imf_seg*4, std::ios::cur);
                // Skip some more stuff
                state = "read header";
                lib.seekg(3*4, std::ios::cur);
                float info;
                lib.read(info);
                extras.resize(info == 0 ? 12 : 10, age.size());
                // Skip the end of the info section (fixed size)
                lib.seekg(1*4 + 80 + 4*4 + 80 + 80 + 2 + 2 + 3*4, std::ios::cur);
            }

            // Read wavelength grid
            {
                state = "read wavelength grid";
                std::int32_t nwave = 0;
                if (nwave < 0) {
                    error("invalid library file: number of wavelength values is negative");
                    return false;
                }
                lib.read(nwave);
                lambda.resize(nwave);
                lib.read(lambda);
            }

            // Read fluxes
            {
                state = "read model fluxes";
                fluxes.resize(age.size(), lambda.size());
                for (uint_t i : range(age)) {
                    // Discard extra data
                    lib.seekg(2*4, std::ios::cur);

                    // Read data
                    std::int32_t nlam = 0;
                    lib.read(nlam);
                    if (nlam < 0) {
                        error("invalid library file: number of wavelength values "
                            "is negative for time step id=", i);
                        return false;
                    }
                    if (uint_t(nlam) > lambda.size()) {
                        error("invalid library file: too many wavelength values "
                            "at time step id=", i);
                        return false;
                    }

                    lib.read(&fluxes(i,0), nlam);

                    // Discard extra data
                    std::int32_t nspec = 0;
                    lib.read(nspec);
                    lib.seekg(nspec*4, std::ios::cur);
                }
            }

            // Read extras
            state = "read extra data (mass, sfr)";
            for (uint_t i : range(extras.dims[0])) {
                // Discard extra data
                lib.seekg(2*4, std::ios::cur);

                // Read data
                std::int32_t ntime = 0;
                lib.read(ntime);
                if (ntime < 0) {
                    error("invalid library file: number of time steps is negative "
                        "for extra data id=", i);
                    return false;
                }
                if (uint_t(ntime) > extras.dims[1]) {
                    error("invalid library file: too many time steps for extra data id=", i);
                    return false;
                }

                lib.read(&extras(i,0), ntime);
            }
        } catch (...) {
            error("could not read data in library file '", filename, "'");
            error("could not ", state);
            error("the file is probably corrupted, try re-downloading it");
            return false;
        }

        mass = extras(1,_);
        sfr = extras(2,_);

        return true;
    }
};

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
    vec1d madau1995(vec1d lam, double z) {
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

vec1d gridder_t::build_igm_absorption(const vec1f& lambda, float z) const {
    return igm::madau1995(lambda, z);
}

std::string gridder_t::get_library_file(uint_t im, uint_t it) const {
    std::string stau = strn(output.tau[it]);
    if (stau.find(".") == stau.npos) stau += ".0";

    return opts.library_dir+"ised_"+opts.sfh+"."+opts.resolution+"/"+
        opts.library+"_"+opts.resolution+"_"+opts.imf+
        "_z"+replace(strn(output.metal[im]), "0.", "")+"_ltau"+stau+".ised";
}

bool gridder_t::build_and_send(fitter_t& fitter) {
    model_t model;
    model.flux.resize(input.lambda.size());

    if (read_from_cache) {
        auto pg = progress_start(dims[0]*dims[1]*dims[2]*dims[3]*dims[4]);
        for (uint_t im : range(output.metal))
        for (uint_t it : range(output.tau))
        for (uint_t ia : range(output.age))
        for (uint_t id : range(output.av))
        for (uint_t iz : range(output.z)) {
            if (cache.read_model(model)) {
                if (model.igrid != model_id(im, it, ia, id, iz)) {
                    error("data in the cache file are not what was expected");
                    error("the cache is probably corrupted, please remove it and try again");
                    return false;
                }

                if ((opts.no_max_age || output.age[ia] <= auniv[iz])) {
                    // Send to fitter
                    model.im = im; model.it = it; model.ia = ia; model.id = id; model.iz = iz;
                    fitter.fit(model);
                }

                if (opts.verbose) progress(pg, 131);
            } else {
                error("could not read data from cache file");
                error("the cache is probably corrupted, please remove it and try again");
                return false;
            }
        }
    } else {
        galaxev_ised ised;

        auto pg = progress_start(dims[0]*dims[1]*dims[2]*dims[3]*dims[4]);
        for (uint_t im : range(output.metal))
        for (uint_t it : range(output.tau)) {
            model.im = im; model.it = it;

            // Load SSP in galaxev ised format
            std::string filename = get_library_file(im, it);
            if (!ised.read(filename)) {
                return false;
            }

            // Make sure input is correct
            phypp_check(is_sorted(ised.age), "galaxev age array is not sorted: ", ised.age);

            // Pre-compute dust law
            vec1d dust_law = build_dust_law(ised.lambda);

            // Pre-compute IGM absorption
            vec2d igm_abs(output.z.size(), ised.lambda.size());
            for (uint_t iz : range(output.z)) {
                igm_abs(iz,_) = build_igm_absorption(ised.lambda, output.z[iz]);
            }

            for (uint_t ia : range(output.age)) {
                model.ia = ia;

                // Interpolate the galaxev grid at the requested age
                vec1f tpl_flux;
                double nage = e10(output.age[ia]);
                auto p = bounds(nage, ised.age);
                if (p[0] == npos) {
                    error("requested age is lower than allowed by the template library (",
                        output.age[ia], " vs. ", log10(ised.age.safe[p[1]]), ")");
                    return false;
                } else if (p[1] == npos) {
                    if (nage > ised.age.safe[p[0]]) {
                        error("requested age is larger than allowed by the template library (",
                            output.age[ia], " vs. ", log10(ised.age.safe[p[0]]), ")");
                        return false;
                    }

                    tpl_flux = ised.fluxes.safe(p[0],_);
                    model.sfr = ised.sfr.safe[p[0]];
                    model.mass = ised.mass.safe[p[0]];
                } else {
                    double x = (output.age.safe[ia] - log10(ised.age.safe[p[0]]))/
                        (log10(ised.age.safe[p[1]]) - log10(ised.age.safe[p[0]]));

                    tpl_flux = ised.fluxes.safe(p[0],_)*(1.0 - x) + ised.fluxes.safe(p[1],_)*x;
                    model.sfr = ised.sfr.safe[p[0]]*(1.0 - x) + ised.sfr.safe[p[1]]*x;
                    model.mass = ised.mass.safe[p[0]]*(1.0 - x) + ised.mass.safe[p[1]]*x;
                }

                for (uint_t id : range(output.av)) {
                    model.id = id;

                    // Apply dust reddening
                    vec1f tpl_att_flux = tpl_flux;
                    if (output.av[id] > 0) {
                        for (uint_t il : range(tpl_att_flux)) {
                            tpl_att_flux.safe[il] *= e10(-0.4*output.av[id]*dust_law.safe[il]);
                        }
                    }

                    for (uint_t iz : range(output.z)) {
                        model.iz = iz;
                        model.igrid = model_id(model.im, model.it, model.ia, model.id, model.iz);
                        vec1f tpl_att_z_lam = ised.lambda;
                        vec1f tpl_att_z_flux = tpl_att_flux;

                        for (uint_t il : range(tpl_att_z_flux)) {
                            // Apply IGM absorption & redshift
                            tpl_att_z_flux.safe[il] *= lum2fl.safe[iz]*igm_abs.safe(iz,il);
                            tpl_att_z_lam.safe[il] *= (1.0 + output.z.safe[iz]);
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

                        // Send to fitter
                        if (opts.no_max_age || output.age[ia] <= auniv[iz]) {
                            fitter.fit(model);
                        }

                        // Cache
                        cache.write_model(model);

                        if (opts.verbose) progress(pg, 131);
                    }
                }
            }
        }
    }

    return true;
}

bool gridder_t::build_template(uint_t im, uint_t it, uint_t ia, uint_t id, uint_t iz,
    vec1f& lam, vec1f& flux, vec1f& iflux) const {

    galaxev_ised ised;

    // Load SSP in galaxev ised format
    std::string filename = get_library_file(im, it);

    if (!ised.read(filename)) {
        return false;
    }

    // Make sure input is correct
    phypp_check(is_sorted(ised.age), "galaxev age array is not sorted: ", ised.age);

    // Interpolate the galaxev grid at the requested age
    double nage = e10(output.age[ia]);
    auto p = bounds(nage, ised.age);
    double mass = 0;
    if (p[0] == npos) {
        error("requested age is lower than allowed by the template library (",
            output.age[ia], " vs. ", log10(ised.age[p[1]]), ")");
        return false;
    } else if (p[1] == npos) {
        if (nage > ised.age[p[0]]) {
            error("requested age is larger than allowed by the template library (",
                output.age[ia], " vs. ", log10(ised.age[p[0]]), ")");
            return false;
        }

        flux = ised.fluxes(p[0],_);
        mass = ised.mass[p[0]];
    } else {
        double x = (output.age[ia] - log10(ised.age[p[0]]))/
            (log10(ised.age[p[1]]) - log10(ised.age[p[0]]));

        flux = ised.fluxes(p[0],_)*(1.0 - x) + ised.fluxes(p[1],_)*x;
        mass = ised.mass[p[0]]*(1.0 - x) + ised.mass[p[1]]*x;
    }

    // Apply dust reddening
    if (output.av[id] > 0) {
        vec1d dust_law = build_dust_law(ised.lambda);

        for (uint_t il : range(flux)) {
            flux.safe[il] *= e10(-0.4*output.av[id]*dust_law.safe[il]);
        }
    }

    // Compute IGM absorption
    vec1d igm_abs = build_igm_absorption(ised.lambda, output.z[iz]);

    lam = ised.lambda;
    for (uint_t il : range(flux)) {
        // Apply IGM absorption & redshift & normalize to unit mass
        flux.safe[il] *= (lum2fl[iz]/mass)*igm_abs.safe[il];
        lam.safe[il] *= (1.0 + output.z[iz]);
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
    vec1u ids = grid_ids(igrid);
    return build_template(ids[0], ids[1], ids[2], ids[3], ids[4], lam, flux, iflux);
}

uint_t gridder_t::model_id(uint_t im, uint_t it, uint_t ia, uint_t id, uint_t iz) const {
    return flat_id(dims, im, it, ia, id, iz);
}

vec1u gridder_t::grid_ids(uint_t iflat) const {
    return mult_ids(dims, iflat);
}
