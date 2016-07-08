#include "fast++.hpp"
#include <phypp/main.hpp>

int phypp_main(int argc, char* argv[]) {
    std::string param_file = (argc >= 2 ? argv[1] : "fast.param");

    // Read input data
    options_t opts;
    input_state_t input;
    if (!read_input(opts, input, param_file)) {
        return 1;
    }

    // Fit them all!
    output_state_t output;
    if (!fit_data(opts, input, output)) {
        return 1;
    }

    return 0;
}
