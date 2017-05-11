#include "fast++.hpp"
#include <phypp/core/main.hpp>

const char* fastpp_version = "0.1";

int phypp_main(int argc, char* argv[]) {
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

    // Initizalize the fitter
    fitter_t fitter(opts, input, gridder, output);

    // Build/read the grid for this redshift
    if (!gridder.build_and_send(fitter)) {
        return 1;
    }

    fitter.find_best_fits();

    // Write output to disk
    write_output(opts, input, gridder, output);

    return 0;
}
