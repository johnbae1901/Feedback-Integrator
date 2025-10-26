#ifndef ESTIMATE_L_H
#define ESTIMATE_L_H

#include <array>
#include "dynamics_perturbed_kepler.h" // PKParams {mu, delta}

/**
 * @brief Spectral (2-)norm of Hessian of V at the state x=[r;v] when (E=E0, L=L0)
 *        V = 0.5*kE*(E-E0)^2 + 0.5*kL*||L-L0||^2
 *        Hess V = kE * (J_E^T J_E) + kL * (J_L^T J_L)
 *
 * @param x   6-vector state [x1,x2,x3, v1,v2,v3]
 * @param P   PKParams {mu, delta} for perturbed Kepler
 * @param kL  weight for angular-momentum term
 * @param kE  weight for energy term
 * @return    || Hess V ||_2  (largest eigenvalue; symmetric PSD)
 */
double estimateL(const std::array<double,6>& x,
                 const PKParams& P,
                 double kL,
                 double kE);

#endif // ESTIMATE_L_H
