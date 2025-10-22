#ifndef EULER_FEEDBACK_ADAPTIVE_LIGHT_H
#define EULER_FEEDBACK_ADAPTIVE_LIGHT_H

#include <vector>
#include <tuple>

/**
 * @brief euler_feedback_adaptive:
 *   - Integrates the system from t=0 to t=tf with step h,
 *   - Returns arrays (x1, x2, x3, v1, v2, v3, t).
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
euler_feedback_adaptive_light(const std::vector<double>& xi,
                        double tf,
                        double h,
                        double mu,
                        double c,
                        double k1,
                        double k2,
                        const std::array<double,3>& L0,
                        const std::array<double,3>& A0,
                        const int m,
                        const double lambda,
                        const double Hmin);

/**
 * @brief dynamics_Euler_feedback
 *   - The 6D right-hand side with an adaptive alpha = 1/(h * lipConstant).
 */
// std::array<double,6> dynamics_Euler_feedback_light(const std::array<double,6>& x,
//                                             double mu,
//                                             double c,
//                                             double h,
//                                             double k1,
//                                             double k2,
//                                             const std::array<double,3>& L0,
//                                             const std::array<double,3>& A0);

std::array<double,6> dynamics_Euler_feedback_light(const std::array<double,6>& x,
                                            double mu,
                                            double c,
                                            double h,
                                            double k1,
                                            double k2,
                                            const std::array<double,3>& L0,
                                            const std::array<double,3>& A0,
                                            const double lipConstant,
                                            const double Hmin);

double euler_feedback_adaptive_error_light(const std::vector<double>& xi,
                        double tf,
                        double h,
                        double mu,
                        double c,
                        double k1,
                        double k2,
                        const std::array<double,3>& L0,
                        const std::array<double,3>& A0,
                        const int m,
                        const double lambda,
                        const double Hmin);

#endif // EULER_FEEDBACK_ADAPTIVE_LIGHT_H
