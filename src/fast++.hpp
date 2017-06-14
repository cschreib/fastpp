#ifndef FASTPP_HPP
#define FASTPP_HPP

#include <phypp/core/typedefs.hpp>
#include <phypp/core/vec.hpp>
#include <phypp/reflex/reflex.hpp>
#include <phypp/core/error.hpp>
#include <phypp/utility/generic.hpp>
#include <phypp/utility/string.hpp>
#include <phypp/utility/time.hpp>
#include <phypp/math/base.hpp>
#include <phypp/astro/astro.hpp>
#include <phypp/io/ascii.hpp>
#include "thread_worker_pool.hpp"

using namespace phypp;

extern const char* fastpp_version;

// Program structures
// ------------------

enum class parallel_choice {
    none, sources, models
};

enum class sfh_type {
    gridded, single
};

// Program options read from parameter file
struct options_t {
    // Input catalog parameters
    std::string catalog;
    float ab_zeropoint = 23.9;
    std::string name_zphot = "z_phot";

    // Input spectrum parameters
    std::string spectrum;
    bool auto_scale = false;

    // Filter database parameters
    std::string filters_res;
    int_t filters_format = 1;

    // Templates parameter
    std::string temp_err_file;
    std::string library_dir;
    std::string library;
    std::string dust_law;
    std::string my_sfh;
    sfh_type    sfh;
    float dust_noll_eb = 1.0;
    float dust_noll_delta = -0.2;

    // Core grid parameters (valid for all models)
    float z_min = 0.001;
    float z_max = 9.0;
    float z_step = 0.02;
    int_t z_step_type = 0;
    float a_v_min = 0;
    float a_v_max = 3;
    float a_v_step = 0.1;
    float log_age_min = 6.3;
    float log_age_max = 10.3;
    float log_age_step = 0.5;
    bool no_max_age = false;

    // Parameters for pre-computed SFHs with galaxev
    std::string resolution;
    std::string name_imf;
    std::string name_sfh;
    float log_tau_min = 6.5;
    float log_tau_max = 10.5;
    float log_tau_step = 0.5;
    vec1f metal;

    // Cosmology parameters
    astro::cosmo_t cosmo;

    // Output parameters
    std::string output_dir;
    std::string output_file;
    bool save_chi_grid = false;
    uint_t n_sim = 0;
    vec1f c_interval = {68.0};
    bool best_fit = false;

    // NB: parameters below were not in original FAST
    // ----------------------------------------------

    // z-phot behavior
    bool force_zphot = false;
    bool best_at_zphot = false;
    float zphot_conf = fnan;

    // LIR prior
    bool use_lir = false;

    // Miscelaneous
    bool verbose = true;

    // Outputs
    float sfr_avg = 0.0;
    bool intrinsic_best_fit = false;
    bool best_sfhs = false;
    std::string sfh_output = "sfr";
    float sfh_step = 10.0;
    vec1s output_columns;

    // Simulations
    bool save_sim = false;
    bool best_from_sim = false;

    // Cache
    bool no_cache = false;

    // Multithreading
    parallel_choice parallel = parallel_choice::none;
    uint_t n_thread = 0;
    uint_t max_queued_fits = 1000;
};

// Filter passband
struct fast_filter_t {
    uint_t id = npos;
    bool spectral = false;
    vec1f wl, tr;
};

// Holds the input state of the program
struct input_state_t {
    // List of filter ID used in the photometric catalog
    vec1u no_filt;                      // [nfilt]
    // List of spectral bins ID in the spectroscopic catalog
    vec1u no_spec;                      // [nspec]
    // Central wavelength of the filters
    vec1d lambda;                       // [nfilt+nspec]

    // Photometry & parameters of all sources in the catalog
    vec1s id;                           // [ngal]
    vec1f zspec;                        // [ngal]
    vec2f zphot;                        // [ngal,1+2*nzconf]
    vec1f zphot_conf;                   // [nzconf]
    vec1f lir, lir_err;                 // [ngal]
    vec2f flux, eflux;                  // [ngal,nfilt+nspec]
    vec1u nobs;                         // [ngal]
    vec1f conf_interval;                // [nconf]

    // Filter database
    // NB: for reference, in FAST: filters = [4,nfilt+nspec], with 4={ID,wl,tr,type}
    // with 'type=1' for photometry and 'type=0' for spectroscopy
    vec<1,fast_filter_t> filters; // [nfilt+nspec]

    // Template error function
    vec1f tplerr_lam, tplerr_err;

    // Baked grid cache name
    std::string name;
};

struct grid_id {
    static const constexpr uint_t z = 0, av = 1, age = 2, metal = 3, custom = 4;
};

struct prop_id {
    static const constexpr uint_t mass = 0, sfr = 1, ssfr = 2, ldust = 3, custom = 4;
};

// Holds the output state of the program
struct output_state_t {
    // Grid parameters
    vec<1,vec1f> grid;               // [ngrid][...]

    // Model properties
    vec1s param_names;                // [ngrid+nprop]
    vec1s param_descriptions;         // [ngrid+nprop]
    vec1b param_scale;                // [ngrid+nprop]
    vec1b param_log;                  // [ngrid+nprop]
    vec1f param_precision;            // [ngrid+nprop]

    // Best fits
    vec3f best_params;               // [ngal,ngrid+nprop,1+nconf]
    vec1f best_chi2;                 // [ngal]
    vec1u best_model;                // [ngal]

    // Monte Carlo simulations
    vec3f mc_best_props;             // [ngal,nprop,nsim]
    vec2f mc_best_chi2;              // [ngal,nsim]
    vec2u mc_best_model;             // [ngal,nsim]

    // For thread safety
    std::mutex fit_result_mutex;
};

// Structure holding the integrated fluxes of a model and associated physical parameters
struct model_t {
    vec1f flux;                      // [nfilt+nspec]
    vec1f props;                     // [nprop]
    uint_t igrid = npos;
};

struct fitter_t;

// Build the grid of models and sends models to the fitter
struct gridder_t {
    const options_t& opts;
    const input_state_t& input;
    output_state_t& output;

    struct cache_manager_t {
        std::fstream cache_file;
        std::string cache_filename;

        void write_model(const model_t& model);
        bool read_model(model_t& model);
    };

    bool read_from_cache = true;
    cache_manager_t cache;

    vec1u grid_dims;                 // [ngrid]
    vec1u grid_dims_pitch;           // [ngrid]

    vec1d lum2fl;                    // [nz]
    vec1d auniv;                     // [nz]
    uint_t nparam = 0, nprop = 0, nfreeparam = 0, nmodel = 0;

    explicit gridder_t(const options_t& opts, const input_state_t& input, output_state_t& output);

    bool check_options() const;

    bool build_and_send(fitter_t& fitter);
    bool build_template(uint_t igrid, vec1f& lam, vec1f& flux, vec1f& iflux) const;
    bool build_template_nodust(uint_t igrid, vec1f& lam, vec1f& flux, vec1f& iflux) const;
    bool get_sfh(uint_t iflat, const vec1f& t, vec1f& sfh) const;

    uint_t model_id(const vec1u& ids) const;
    vec1u grid_ids(uint_t iflat) const;

private :
    bool build_and_send_ised(fitter_t& fitter);
    bool build_template_ised_impl(uint_t im, uint_t it, uint_t iz, double nage, double av,
        vec1f& lam, vec1f& flux, vec1f& iflux) const;
    bool build_template_ised(uint_t iflat, vec1f& lam, vec1f& flux, vec1f& iflux) const;
    bool build_template_ised_nodust(uint_t iflat, vec1f& lam, vec1f& flux, vec1f& iflux) const;
    bool get_sfh_ised(uint_t iflat, const vec1f& t, vec1f& sfh) const;

    std::string get_library_file(uint_t im, uint_t it) const;
    vec1d build_dust_law(const vec1f& lambda) const;
    vec1d build_igm_absorption(const vec1f& lambda, float z) const;
};

// Fit a model to observed fluxes
struct fitter_t {
    const options_t& opts;
    const input_state_t& input;
    const gridder_t& gridder;
    output_state_t& output;

    struct chi2_output_manager_t {
        std::fstream out_file;
        std::string out_filename;
        uint_t hpos = 0;

        // For thread safety
        std::mutex write_mutex;
    };

    bool save_chi2 = false;
    chi2_output_manager_t ochi2;

    struct model_source_pair {
        model_t model;
        uint_t i0 = 0, i1 = 0;

        model_source_pair() = default;
        explicit model_source_pair(const model_t& m, uint_t ti0, uint_t ti1) :
            model(m), i0(ti0), i1(ti1) {}
    };

    struct workers_multi_source_t {
        fitter_t& fitter;
        thread::worker_pool<model_source_pair> workers;

        explicit workers_multi_source_t(fitter_t& f);
        void process(const model_t& model);
    };

    struct workers_multi_model_t {
        fitter_t& fitter;
        thread::worker_pool<model_t> workers;

        explicit workers_multi_model_t(fitter_t& f);
        void process(const model_t& model);
    };

    std::unique_ptr<workers_multi_source_t> workers_multi_source;
    std::unique_ptr<workers_multi_model_t> workers_multi_model;

    vec2d tpl_err;               // [nz,nfilt+nspec]
    vec1u idz, idzp, idzl, idzu; // [ngal]
    vec2d sim_rnd;               // [nsim,nfilt]

    explicit fitter_t(const options_t& opts, const input_state_t& input, const gridder_t& gridder,
        output_state_t& output);

    void fit(const model_t& model);
    void find_best_fits();

private :
    inline void write_chi2(uint_t igrid, const vec1f& chi2, const vec2f& props, uint_t i0);
    inline void fit_galaxies(const model_t& model, uint_t i0, uint_t i1);
};

// Main functions
// --------------

// fast++-read_input.cpp
bool read_input(options_t& opts, input_state_t& state, const std::string& filename);

// fast++-write_output.cpp
void write_output(const options_t& opts, const input_state_t& state, const gridder_t& gridder,
    const output_state_t& output);


// Helper functions
// ----------------

namespace phypp {
namespace file {
    template<typename S, typename T>
    bool write(S& s, const T& t) {
        s.write(reinterpret_cast<const char*>(&t), sizeof(T));
        return !s.fail();
    }

    template<typename R, typename S, typename T>
    bool write_as(S& s, const T& t) {
        R r = t;
        s.write(reinterpret_cast<const char*>(&r), sizeof(R));
        return !s.fail();
    }

    template<typename S, std::size_t D, typename T>
    bool write(S& s, const vec<D,T>& t) {
        s.write(reinterpret_cast<const char*>(t.data.data()), sizeof(T)*t.size());
        return !s.fail();
    }

    template<typename S, typename T>
    bool read(S& s, T& t) {
        s.read(reinterpret_cast<char*>(&t), sizeof(T));
        return !s.fail();
    }

    template<typename R, typename S, typename T>
    bool read_as(S& s, T& t) {
        R r;
        if (s.read(reinterpret_cast<char*>(&r), sizeof(R))) {
            t = r;
        }

        return !s.fail();
    }

    template<typename S, std::size_t D, typename T>
    bool read(S& s, vec<D,T>& t) {
        s.read(reinterpret_cast<char*>(t.data.data()), sizeof(T)*t.size());
        return !s.fail();
    }
}
}

#endif
