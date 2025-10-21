#ifndef EULER_ADAPTIVE_H
#define EULER_ADAPTIVE_H

#include <array>
#include <tuple>

double estimateL(const double R[3][3], const double Omega[3]);

void dynamics_adaptive(const double state[3][4],
                         const double I[3][3],
                         double k0, double k1, double k2,
                         double E0, const double Pi0[3],
                         double h,
                         double out[3][4]);

void dynamics_adaptive_light(const double state[3][4],
                         const double I[3][3],
                         double k0, double k1, double k2,
                         double E0, const double Pi0[3],
                         double h, double lipConstant,
                         double out[3][4]);

void euler_adaptive(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    double k0, double k1, double k2,
                    double E0, const double Pi0[3],
                    double tf, double h, 
                    int m, double lambda,
                    double *&R_out, double *&Omega_out, double *&t_out,
                    int &N);

double euler_adaptive(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    double k0, double k1, double k2,
                    double E0, const double Pi0[3],
                    double tf, double h,
                    int m, double lambda);

#endif // EULER_ADAPTIVE_H
