#include "fast++.hpp"
#include <phypp/utility/thread.hpp>

fitter_t::fitter_t(const options_t& opt, const input_state_t& inp, output_state_t& out) :
    opts(opt), input(inp), output(out) {}

void fitter_t::fit(model_t model) {

}

bool fit_data(const options_t& opts, const input_state_t& input, output_state_t& output) {
    // Initialize the grid
    gridder_t grid(opts, input, output);
    // Initizalize the fitter
    fitter_t fitter(opts, input, output);

    // Build/read the grid for this redshift
    if (!grid.build_and_send(fitter)) {
        return false;
    }

    return true;
}
