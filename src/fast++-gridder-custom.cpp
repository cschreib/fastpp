#include "fast++.hpp"
#include <phypp/utility/thread.hpp>

extern "C" {
    #include "tinyexpr.h"
}

bool gridder_t::build_and_send_custom(fitter_t& fitter) {
    return false;
}

bool gridder_t::build_template_custom(uint_t iflat, vec1f& lam, vec1f& flux) const {
    return false;
}

bool gridder_t::get_sfh_custom(uint_t iflat, const vec1f& t, vec1f& sfh) const {
    return false;
}
