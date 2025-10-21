#ifndef EULER_VANILLA_H
#define EULER_VANILLA_H

#include <array>
#include <tuple>

void dynamics_euler(const double state[3][4], const double I[3][3], double out[3][4]);
void euler_vanilla(const double R0[3][3],
                   const double Omega0[3],
                   const double I[3][3],
                   double tf, double h,
                   double *&R_out, double *&Omega_out, double *&t_out,
                   int &N);
double euler_vanilla(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    const double E0, 
                    const double Pi0[3], 
                    const double k0,
                    const double k1,
                    const double k2,
                    double tf, double h);

// double euler_vanilla_error(const std::vector<double>& xi, double tf, double h, 
//                             const std::array<double,3>& L0,
//                             const std::array<double,3>& A0,
//                             double k1,
//                             double k2,
//                             double mu);

#endif // EULER_VANILLA_H
