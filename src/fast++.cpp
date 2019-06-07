#include "fast++.hpp"
#include <vif/core/main.hpp>

#ifndef FASTPP_GIT_HASH
#define FASTPP_GIT_HASH ""
#endif

const char* fastpp_version = "1.2.2";
const char* fastpp_git_hash = FASTPP_GIT_HASH;

const constexpr uint_t grid_id::z;
const constexpr uint_t grid_id::av;
const constexpr uint_t grid_id::age;
const constexpr uint_t grid_id::metal;
const constexpr uint_t grid_id::custom;

const constexpr uint_t prop_id::scale;
const constexpr uint_t prop_id::spec_scale;
const constexpr uint_t prop_id::mass;
const constexpr uint_t prop_id::sfr;
const constexpr uint_t prop_id::ssfr;
const constexpr uint_t prop_id::ldust;
const constexpr uint_t prop_id::lion;
const constexpr uint_t prop_id::mform;
const constexpr uint_t prop_id::custom;

const constexpr uint_t log_style::none;
const constexpr uint_t log_style::decimal;
const constexpr uint_t log_style::abmag;


int vif_main(int argc, char* argv[]) {
    std::string param_file = (argc >= 2 ? argv[1] : "fast.param");

    // Read input data
    options_t opts;
    input_state_t input;
    if (!read_input(opts, input, param_file)) {
        return 1;
    }

    // Initialize the grid
    output_state_t output;
    gridder_t gridder(opts, input, output);
    if (!gridder.check_options()) {
        return 1;
    }

    if (!opts.make_seds.empty()) {
        // Write SEDs if asked, then terminate
        gridder.write_seds();
        return 0;
    }

    if (gridder.read_from_cache && opts.parallel == parallel_choice::generators &&
        opts.n_thread > 0) {
        if (opts.verbose) {
            note("using cache, switched parallel execution from 'generators' to 'models'");
        }

        opts.parallel = parallel_choice::models;
    }

    // Initizalize the fitter
    fitter_t fitter(opts, input, gridder, output);

    // Build/read the grid and fit galaxies
    if (!gridder.build_and_send(fitter)) {
        return 1;
    }

    // Compile results
    fitter.find_best_fits();

    // Write output to disk
    write_output(opts, input, gridder, fitter, output);

    return 0;
}
