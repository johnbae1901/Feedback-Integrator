#ifndef STORMER_VERLET_B_H
#define STORMER_VERLET_B_H

#include <array>
#include <vector>
#include <tuple>
#include "dynamics_perturbed_kepler.h"

/**
 * @brief Integrates a 6D state (x, y, z, vx, vy, vz) using a Stormer–Verlet (leapfrog) integrator.
 *
 * @param xi   Initial state as an std::array<double,6> (order: [x, y, z, vx, vy, vz]).
 * @param tf   Final time.
 * @param h    Time step.
 * @param mu   Gravitational (or other system) parameter.
 *
 * @return A tuple of 7 vectors: (x1, x2, x3, v1, v2, v3, t), where each vector holds the
 *         time history of that component.
 */
std::tuple<
    std::vector<double>, // x1
    std::vector<double>, // x2
    std::vector<double>, // x3
    std::vector<double>, // v1
    std::vector<double>, // v2
    std::vector<double>, // v3
    std::vector<double>  // t
>
stormer_verlet_B(const std::vector<double>& xi,
                 double tf,
                 double h,
                 const PKParams& P);

std::tuple<double, double, double> 
    stormer_verlet_B_error(const std::vector<double>& xi,
                            double tf,
                            double h,
                            double k1,  // -> kL
                            double k2,  // -> kE
                            const PKParams& P);
#endif // STORMER_VERLET_B_H
