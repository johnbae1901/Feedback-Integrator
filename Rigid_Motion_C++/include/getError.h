#ifndef GET_ERROR_H
#define GET_ERROR_H

#include "basic_operations.h"

struct ErrorReport {
    double frob_RT_R_minus_I_sq; // (1) ||R^T R - I||_F^2
    double abs_deltaE;           // (2) |E(Omega,I) - E0|
    double deltaPi_norm_sq;      // (3) ||Pi - Pi0||_2^2
    double weighted_sum;         // (4) (k0/4)*(1) + (k1/2)*|ΔE|^2 + (k2/2)*(3)
};

inline double frobeniusNormRT_R_minus_I_Square(const double R[3][3]) {
    double sum = 0.0;
    // Loop over indices of the 3x3 matrix for (RᵀR)
    // (RᵀR)_{ij} = sum_{k=0}^{2} R[k][i] * R[k][j]
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double val = 0.0;
            for (int k = 0; k < 3; k++) {
                val += R[k][i] * R[k][j];
            }
            // Subtract the identity matrix contribution (only on the diagonal)
            if (i == j) {
                val -= 1.0;
            }
            // Sum up the square of the element.
            sum += val * val;
        }
    }
    return sum;
}

inline ErrorReport makeErrorReport_(const double R[3][3],
                                    const double Omega[3],
                                    const double I[3][3],
                                    const double E0,
                                    const double Pi0[3],
                                    const double k0,
                                    const double k1,
                                    const double k2)
{
    const double frob_sq = frobeniusNormRT_R_minus_I_Square(R);

    const double dE = E(Omega, I) - E0;
    const double abs_dE = std::abs(dE);
    const double dE_sq  = dE * dE;

    double P[3];
    Pi(R, Omega, I, P);
    const double dPi[3] = { P[0] - Pi0[0], P[1] - Pi0[1], P[2] - Pi0[2] };
    const double dPi_sq = norm3Square(dPi);

    ErrorReport out;
    out.frob_RT_R_minus_I_sq = frob_sq;
    out.abs_deltaE           = abs_dE;
    out.deltaPi_norm_sq      = dPi_sq;
    out.weighted_sum         = (k0 * 0.25) * frob_sq + (k1 * 0.5) * dE_sq + (k2 * 0.5) * dPi_sq;
    return out;
}

inline ErrorReport getError(const double R[3][3],
                            const double Omega[3],
                            const double I[3][3],
                            const double E0,
                            const double Pi0[3],
                            const double k0,
                            const double k1,
                            const double k2)
{
    return makeErrorReport_(R, Omega, I, E0, Pi0, k0, k1, k2);
}

inline ErrorReport getError(const double state[3][4],
                            const double I[3][3],
                            const double E0,
                            const double Pi0[3],
                            const double k0,
                            const double k1,
                            const double k2)
{
    double R[3][3];
    double Omega[3];
    // state에서 R, Omega 추출
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) R[i][j] = state[i][j];
        Omega[i] = state[i][3];
    }
    return makeErrorReport_(R, Omega, I, E0, Pi0, k0, k1, k2);
}

#endif