#ifndef FASTPP_SSP_HPP
#define FASTPP_SSP_HPP

#include <phypp/core/typedefs.hpp>
#include <phypp/core/vec.hpp>
#include <phypp/core/range.hpp>
#include <phypp/math/reduce.hpp>

using namespace phypp;

// SED libraries
struct ssp_bc03 {
    vec1d age;
    vec1d mass;
    vec1d lambda;
    vec2d sed;

    bool read_ascii(std::string filename);
    bool read_fits(std::string filename, bool noflux);
    bool read(std::string filename, bool noflux = false);

    template<typename F>
    void integrate(const vec1d& iage, const vec1d& isfr, F&& func) {
        double t2 = 0.0;
        uint_t ihint = npos;
        for (uint_t it : range(age)) {
            double t1 = t2;
            if (it < age.size()-1) {
                t2 = 0.5*(age.safe[it] + age.safe[it+1]);
            } else {
                t2 = age.back();
            }

            if (t2 <= iage.front()) {
                continue;
            }

            t1 = max(t1, iage.front());
            t2 = min(t2, iage.back());

            func(it, integrate_hinted(iage, isfr, ihint, t1, t2));

            if (t2 >= iage.back()) {
                break;
            }
        }
    }
};

#endif
