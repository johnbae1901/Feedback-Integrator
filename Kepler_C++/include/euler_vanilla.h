#ifndef EULER_VANILLA_H
#define EULER_VANILLA_H

#include <array>
#include <tuple>
#include <vector>

/** 
 * @brief Dynamics for the Euler integrator: returns [dx, dv] for 6D state 
 */
std::array<double,6> dynamics_Euler(const std::array<double,6>& x, double mu);

/**
 * @brief Euler vanilla integrator
 * @return A tuple of 7 vectors (x1, x2, x3, v1, v2, v3, t)
 */
std::tuple<
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>
> euler_vanilla(const std::vector<double>& xi, double tf, double h, double mu);

std::tuple<double, double, double> 
        euler_vanilla_error(const std::vector<double>& xi, double tf, double h, 
                            const std::array<double,3>& L0,
                            const std::array<double,3>& A0,
                            double k1,
                            double k2,
                            double mu);

#endif // EULER_VANILLA_H
