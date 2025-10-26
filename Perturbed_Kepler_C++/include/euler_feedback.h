#ifndef EULER_FEEDBACK_H
#define EULER_FEEDBACK_H

#include <vector>
#include <tuple>
#include "basic_operations.h"
#include "dynamics_perturbed_kepler.h"

/**
 * @brief EL-feedback RHS for perturbed Kepler.
 *        Uses fixed references (L0, E0) passed in from the integrator.
 *        kL_eff = alpha*kL, kE_eff = alpha*kE are already scaled gains.
 *
 * xdot = v - kE_eff*dE*Ur(r)*x/r - kL_eff*(v × dL)
 * vdot = a(x) - kE_eff*dE*v      - kL_eff*(dL × x)
 */
std::array<double,6> dynamics_Euler_feedback(const std::array<double,6>& x,
                                             const PKParams& P,
                                             double kL_eff,
                                             double kE_eff,
                                             const std::array<double,3>& L0,
                                             double E0);

/**
 * @brief Euler integration of the 6D system with feedback:
 *        x_{n+1} = x_n + h * dynamics_Euler_feedback(x_n, ...).
 *
 * @param xi    Initial 6D state: (x, y, z, vx, vy, vz).
 * @param tf    Final time.
 * @param h     Time step size.
 * @param mu    Gravitational (or system) parameter.
 * @param alpha Feedback parameter.
 * @param k1, k2 Gains for the feedback law.
 * @param L0, A0 Reference vectors for feedback correction.
 *
 * @return A 7-tuple of std::vectors: (x1, x2, x3, v1, v2, v3, t).
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
euler_feedback(const std::vector<double>& xi,
               double tf,
               double h,
               const PKParams& P,
               double alpha,
               double kL,
               double kE);

std::tuple<double, double, double> 
    euler_feedback_error(const std::vector<double>& xi,
                        double tf,
                        double h,
                        const PKParams& P,
                        double alpha,
                        double kL,
                        double kE);
#endif // EULER_FEEDBACK_H
