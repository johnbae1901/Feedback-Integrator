#ifndef GETL0E0_H
#define GETL0E0_H

#include <vector>
#include <tuple>
#include <array>
#include "dynamics_perturbed_kepler.h"

void getL0E0(const std::vector<double>& xi, const PKParams& P)
{
    // initial state
    double x1 = xi[0], x2 = xi[1], x3 = xi[2];
    double v1 = xi[3], v2 = xi[4], v3 = xi[5];

    std::array<double,6> s = {x1,x2,x3, v1,v2,v3};

    // L0, E0 내부 계산
    const std::array<double,3> L0 = {
        x2*v3 - x3*v2,
        x3*v1 - x1*v3,
        x1*v2 - x2*v1
    };
    const double r0  = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    const double v02 = (v1*v1 + v2*v2 + v3*v3);
    const double E0  = 0.5*v02 - P.mu/r0 - P.delta/(r0*r0*r0);

    std::cout << "L0         : " << L0[0] << ", " << L0[1] << ", " << L0[2] << "\n";
    std::cout << "E0         : " << E0 << "\n";
}



#endif