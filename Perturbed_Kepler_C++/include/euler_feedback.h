#ifndef EULER_FEEDBACK_H
#define EULER_FEEDBACK_H

#include <vector>
#include <tuple>
#include "basic_operations.h"

/**
 * @brief Computes the feedback dynamics for the 6D state x = [r(3), v(3)]:
 *        Fh = [v, -mu * r / ||r||^3] - alpha * [delVx, delVv].
 *
 * @param x      6D state vector: (x, y, z, vx, vy, vz).
 * @param mu     Gravitational parameter (or relevant constant).
 * @param alpha  Scalar feedback parameter.
 * @param k1, k2 Additional feedback gains.
 * @param L0, A0 Reference 3D vectors for feedback calculations.
 *
 * @return A 6-element vector representing [dx, dv].
 */
std::array<double,6> dynamics_Euler_feedback(const std::array<double,6>& x,
                                            double mu,
                                            double alpha,
                                            double k1,
                                            double k2,
                                            const std::array<double,3>& L0,
                                            const std::array<double,3>& A0);

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
               double mu,
               double alpha,
               double k1,
               double k2,
               const std::array<double,3>& L0,
               const std::array<double,3>& A0);

std::tuple<double, double, double> 
    euler_feedback_error(const std::vector<double>& xi,
                        double tf,
                        double h,
                        double mu,
                        double alpha,
                        double k1,
                        double k2,
                        const std::array<double,3>& L0,
                        const std::array<double,3>& A0);
#endif // EULER_FEEDBACK_H
