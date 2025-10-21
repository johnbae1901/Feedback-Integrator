#ifndef GET_ERROR_H
#define GET_ERROR_H

#include "basic_operations.h"

// Computes the Frobenius norm of (RᵀR - I)
// where R is a 3x3 matrix.
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

inline double getError(const double R[3][3], 
                    const double Omega[3], 
                    const double I[3][3], 
                    const double E0, 
                    const double Pi0[3], 
                    const double k0,
                    const double k1,
                    const double k2)
{
    const double deltaE = E(Omega, I) - E0;
    double P[3];
    Pi(R, Omega, I, P);
    const double deltaPi[3] = {P[0] - Pi0[0], P[1] - Pi0[1], P[2] - Pi0[2]};
    return k0/4 * frobeniusNormRT_R_minus_I_Square(R) + k1/2 * deltaE * deltaE + k2/2 * norm3Square(deltaPi);
}

inline double getError(const double state[3][4],
                    const double I[3][3], 
                    const double E0, 
                    const double Pi0[3], 
                    const double k0,
                    const double k1,
                    const double k2)
{
    double R[3][3];
    double Omega[3];
    // Extract R and Omega from state.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R[i][j] = state[i][j];
        }
        Omega[i] = state[i][3];
    }
    const double deltaE = E(Omega, I) - E0;
    double P[3];
    Pi(R, Omega, I, P);
    const double deltaPi[3] = {P[0] - Pi0[0], P[1] - Pi0[1], P[2] - Pi0[2]};
    return k0/4 * frobeniusNormRT_R_minus_I_Square(R) + k1/2 * deltaE * deltaE + k2/2 * norm3Square(deltaPi);
}
                    
#endif