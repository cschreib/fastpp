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

using namespace phypp;

// Program structures
// ------------------

// Program options read from parameter file
struct options_t {
    std::string catalog;
    float ab_zeropoint = 23.9;
    std::string filters_res;
    int_t filters_format = 1;
    std::string temp_err_file;
    std::string name_zphot = "z_phot";
    std::string spectrum;
    bool auto_scale = false;
    std::string output_dir;
    std::string output_file;
    uint_t n_sim = 0;
    vec1f c_interval = {68.0};
    bool best_fit = false;
    std::string library_dir;
    std::string library;
    std::string resolution;
    std::string imf;
    std::string sfh;
    std::string dust_law;
    std::string my_sfh;
    float dust_noll_eb = 1.0;
    float dust_noll_delta = -0.2;
    float log_tau_min = 6.5;
    float log_tau_max = 10.5;
    float log_tau_step = 0.5;
    float log_age_min = 6.3;
    float log_age_max = 10.3;
    float log_age_step = 0.5;
    bool no_max_age = false;
    float z_min = 0.001;
    float z_max = 9.0;
    float z_step = 0.02;
    int_t z_step_type = 0;
    float a_v_min = 0;
    float a_v_max = 3;
    float a_v_step = 0.1;
    vec1f metal;
    astro::cosmo_t cosmo;
    bool save_chi_grid = false;

    // Not in original FAST
    bool verbose = false;
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
    vec1u id;                           // [ngal]
    vec1f zphot, zspec;                 // [ngal]
    vec2f zphot_low, zphot_up;          // [ngal,ninterval]
    vec2f flux, eflux;                  // [ngal,nfilt+nspec]

    // Filter database
    // NB: for reference, in FAST: filters = [4,nfilt+nspec], with 4={ID,wl,tr,type}
    // with 'type=1' for photometry and 'type=0' for spectroscopy
    vec<1,fast_filter_t> filters; // [nfilt+nspec]

    // Template error function
    vec1f tplerr_lam, tplerr_err;

    // Baked grid cache name
    std::string name;
};

// Holds the output state of the program
struct output_state_t {
    vec1f z, metal, av, age, tau;
    vec5f chi2;
};

// Structure holding the integrated fluxes of a model and associated physical parameters
struct model_t {
    vec1f flux;
    float mass = 0, sfr = 0;
    uint_t im, it, ia, id, iz;
};

// Fit a model to observed fluxes
struct fitter_t {
    const options_t& opts;
    const input_state_t& input;
    output_state_t& output;

    explicit fitter_t(const options_t& opts, const input_state_t& input, output_state_t& output);

    void fit(model_t model);
};

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

    explicit gridder_t(const options_t& opts, const input_state_t& input, output_state_t& output);

    bool build_and_send(fitter_t& fitter);
};

// Main functions
// --------------

// fast++-read_input.cpp
bool read_input(options_t& opts, input_state_t& state, const std::string& filename);

// fast++-fitter.cpp
bool fit_data(const options_t& opts, const input_state_t& input, output_state_t& output);

#endif
