#include "fast++.hpp"

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
    bool read(std::string filename, bool noflux = false) {
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

                // Make sure input is correct
                phypp_check(is_sorted(age), "galaxev age array is not sorted: ", age);
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
            std::int32_t nwave = 0; {
                state = "read wavelength grid";
                if (nwave < 0) {
                    error("invalid library file: number of wavelength values is negative");
                    return false;
                }
                lib.read(nwave);
                if (noflux) {
                    lib.seekg(nwave*4, std::ios::cur);
                } else {
                    lambda.resize(nwave);
                    lib.read(lambda);
                }
            }

            // Read fluxes
            {
                state = "read model fluxes";

                if (!noflux) {
                    fluxes.resize(age.size(), lambda.size());
                }

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
                    if (nlam > nwave) {
                        error("invalid library file: too many wavelength values "
                            "at time step id=", i);
                        return false;
                    }

                    if (noflux) {
                        lib.seekg(nlam*4, std::ios::cur);
                    } else {
                        lib.read(&fluxes(i,0), nlam);
                    }

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

std::string gridder_t::get_library_file_ised(uint_t im, uint_t it) const {
    std::string stau = strn(output.grid[grid_id::custom+0][it]);
    if (stau.find(".") == stau.npos) stau += ".0";

    return opts.library_dir+"ised_"+opts.name_sfh+"."+opts.resolution+"/"+
        opts.library+"_"+opts.resolution+"_"+opts.name_imf+
        "_z"+replace(strn(output.grid[grid_id::metal][im]), "0.", "")+"_ltau"+stau+".ised";
}

bool gridder_t::get_age_bounds(const vec1f& ised_age, float nage,
    std::array<uint_t,2>& p, double& x) const {

    p = bounds(nage, ised_age);

    if (p[0] == npos) {
        error("requested age is lower than allowed by the template library (",
            log10(nage), " vs. ", log10(ised_age.safe[p[1]]), ")");
        return false;
    } else if (p[1] == npos) {
        if (nage > ised_age.safe[p[0]]) {
            error("requested age is larger than allowed by the template library (",
                log10(nage), " vs. ", log10(ised_age.safe[p[0]]), ")");
            return false;
        }

        // We picked exactly the oldest age of the library
        p[1] = p[0];
        p[0] = p[1]-1;
    }

    x = (log10(nage) - log10(ised_age.safe[p[0]]))/
        (log10(ised_age.safe[p[1]]) - log10(ised_age.safe[p[0]]));

    return true;
}

bool gridder_t::build_and_send_ised(fitter_t& fitter) {
    model_t model;
    model.flux.resize(input.lambda.size());
    model.props.resize(nprop);
    vec1u idm(nparam);

    galaxev_ised ised;

    vec1f& output_metal = output.grid[grid_id::metal];
    vec1f& output_tau = output.grid[grid_id::custom+0];
    vec1f& output_age = output.grid[grid_id::age];
    vec1f& output_z = output.grid[grid_id::z];

    float& model_mass = model.props[prop_id::mass];
    float& model_sfr = model.props[prop_id::sfr];
    float& model_ssfr = model.props[prop_id::ssfr];
    float& model_a2t = model.props[prop_id::custom+0];

    auto pg = progress_start(nmodel);
    for (uint_t im : range(output_metal))
    for (uint_t it : range(output_tau)) {
        idm[grid_id::metal] = im;
        idm[grid_id::custom] = it;

        // Load CSP
        std::string filename = get_library_file_ised(im, it);
        if (!ised.read(filename)) {
            return false;
        }

        // Pre-compute dust law & IGM absorption (they don't change with age)
        vec1d dust_law = build_dust_law(ised.lambda);
        vec2d igm_abs = build_igm_absorption(output_z, ised.lambda);

        for (uint_t ia : range(output_age)) {
            idm[grid_id::age] = ia;

            // Interpolate the galaxev grid at the requested age
            vec1f tpl_flux;
            double nage = e10(output_age[ia]);

            std::array<uint_t,2> p;
            double x;
            if (!get_age_bounds(ised.age, nage, p, x)) {
                return false;
            }

            tpl_flux = (1.0 - x)*ised.fluxes.safe(p[0],_) + x*ised.fluxes.safe(p[1],_);
            model_mass = (1.0 - x)*ised.mass.safe[p[0]] + x*ised.mass.safe[p[1]];

            if (opts.sfr_avg > 0) {
                // Average SFR over the past X yr
                double t1 = nage;
                double t0 = max(t1 - opts.sfr_avg, ised.age[0]);
                model_sfr = integrate(ised.age, ised.sfr, t0, t1)/opts.sfr_avg;
            } else {
                // Use instantaneous SFR
                model_sfr = (1.0 - x)*ised.sfr.safe[p[0]] + x*ised.sfr.safe[p[1]];
            }

            model_ssfr = model_sfr/model_mass;
            model_a2t = output_age[ia] - output_tau[it];

            // The rest is not specific to the SFH, use generic code
            build_and_send_impl(fitter, pg, ised.lambda, tpl_flux, dust_law, igm_abs,
                output_age[ia], idm, model);
        }
    }

    return true;
}

bool gridder_t::build_template_ised(uint_t iflat, vec1f& lam, vec1f& flux) const {
    vec1u idm = grid_ids(iflat);
    uint_t ia = idm[grid_id::age];
    uint_t it = idm[grid_id::custom+0];
    uint_t im = idm[grid_id::metal];

    galaxev_ised ised;

    // Load CSP
    std::string filename = get_library_file_ised(im, it);
    if (!ised.read(filename)) {
        return false;
    }

    // Interpolate the galaxev grid at the requested age
    std::array<uint_t,2> p;
    double x;
    double nage = e10(output.grid[grid_id::age][ia]);
    if (!get_age_bounds(ised.age, nage, p, x)) {
        return false;
    }

    lam = ised.lambda;
    flux = ised.fluxes(p[0],_)*(1.0 - x) + ised.fluxes(p[1],_)*x;
    float mass = ised.mass[p[0]]*(1.0 - x) + ised.mass[p[1]]*x;
    flux /= mass;

    return true;
}

bool gridder_t::get_sfh_ised(uint_t iflat, const vec1d& t, vec1d& sfh,
    const std::string& type) const {

    uint_t im, it, ia, iz; {
        vec1u idm = grid_ids(iflat);
        iz = idm[grid_id::z];
        ia = idm[grid_id::age];
        it = idm[grid_id::custom+0];
        im = idm[grid_id::metal];
    }

    galaxev_ised ised;

    // Load CSP (only extras)
    std::string filename = get_library_file_ised(im, it);
    if (!ised.read(filename, true)) {
        return false;
    }

    double nage = e10(output.grid[grid_id::age][ia]);
    double age_obs = e10(auniv[iz]);
    double age_born = age_obs - nage;

    uint_t i0 = upper_bound(age_born, t);
    uint_t i1 = upper_bound(age_obs, t);

    if (i1 == npos) {
        i1 = t.size()-1;
    }
    if (i0 == npos) {
        i0 = 0;
    }

    if (type == "sfr") {
        sfh = replicate(0.0, t.size());
        sfh[i0-_-i1] = interpolate(ised.sfr, ised.age, t[i0-_-i1] - age_born);

        // Interpolate the galaxev grid at the requested age
        std::array<uint_t,2> p;
        double x;
        if (!get_age_bounds(ised.age, nage, p, x)) {
            return false;
        }

        double ised_mass = ised.mass[p[0]]*(1.0 - x) + ised.mass[p[1]]*x;
        sfh /= ised_mass;
    } else if (type == "mass") {
        sfh = replicate(0.0, t.size());
        sfh[i0-_-i1] = interpolate(ised.mass, ised.age, t[i0-_-i1] - age_born);
        sfh /= interpolate(sfh, t, age_obs);
        sfh[i1-_] = 1.0; // NB: ignores mass loss after time of observations!
    } else {
        error("unknown SFH type '", type, "'");
        return false;
    }

    return true;
}
