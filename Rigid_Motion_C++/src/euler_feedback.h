#ifndef EULER_FEEDBACK_H
#define EULER_FEEDBACK_H

#include <array>
#include <tuple>

void dynamics_feedback(const double state[3][4],
                       const double I[3][3],
                       const double I_inv[3][3],
                       double k0, double k1, double k2,
                       double E0, const double Pi0[3],
                       double alpha,
                       double out[3][4]);
void euler_feedback(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    double k0, double k1, double k2,
                    double E0, const double Pi0[3],
                    double alpha,
                    double tf, double h,
                    double *&R_out, double *&Omega_out, double *&t_out,
                    int &N);
std::tuple<
    double, 
    double, 
    double, 
    double> 
    euler_feedback(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    double k0, double k1, double k2,
                    double E0, const double Pi0[3],
                    double alpha,
                    double tf, double h);

#endif // EULER_FEEDBACK_H
