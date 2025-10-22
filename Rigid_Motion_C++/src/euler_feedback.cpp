#include <cmath>      // for std::ceil, std::sqrt
#include <vector>
#include <tuple>
#include <array>
#include <iostream>
#include "../include/basic_operations.h"
#include "../include/getError.h"

using namespace std;

// Computes the feedback dynamics.
//  - state: a 3x4 array where the first 3 columns are the rotation matrix R and the 4th column is Omega.
//  - I: 3x3 inertia matrix.
//  - k0, k1, k2: scalar feedback gains.
//  - E0: scalar reference energy.
//  - Pi0: 3-vector reference momentum.
//  - alpha: scalar feedback factor.
// The output is a 3x4 array (first 3 columns give R_dot and the 4th column gives Omega_dot).
void dynamics_feedback(const double state[3][4],
                       const double I[3][3],
                       const double I_inv[3][3],
                       double k0, double k1, double k2,
                       double E0, const double Pi0[3],
                       double alpha,
                       double out[3][4])
{
    // 1. Extract R and Omega from 'state'
    double R[3][3];
    double Omega[3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R[i][j] = state[i][j];
        }
        Omega[i] = state[i][3];
    }
    
    // 2. Compute OmegaHat = getSkew(Omega)
    double OmegaHat[3][3];
    getSkew(Omega, OmegaHat);
    
    // 3. Compute A = R * OmegaHat (a 3x3 matrix)
    double A[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                A[i][j] += R[i][k] * OmegaHat[k][j];
            }
        }
    }
    
    // 4. Compute M = R' * R - I (subtracting the identity matrix)
    double M[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += R[k][i] * R[k][j];  // R' (transpose) means R[k][i]
            }
            M[i][j] = sum - ((i == j) ? 1.0 : 0.0);  // subtract identity
        }
    }
    
    // 5. Compute term1 = k0 * R * M
    double term1[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += R[i][k] * M[k][j];
            }
            term1[i][j] = k0 * sum;
        }
    }
    
    // 6. Compute P = Pi(R, Omega, I)
    double P[3];
    Pi(R, Omega, I, P);
    
    // 7. Compute diffP = P - Pi0 (a 3-element vector)
    double diffP[3];
    for (int i = 0; i < 3; i++) {
        diffP[i] = P[i] - Pi0[i];
    }
    
    // 8. Compute X = Omega' * I (treated as a 1x3 row vector)
    double X[3] = {0};
    for (int j = 0; j < 3; j++) {
        double sum = 0.0;
        for (int i = 0; i < 3; i++) {
            sum += Omega[i] * I[i][j];
        }
        X[j] = sum;
    }
    
    // 9. Compute term2 = k2 * (diffP * X) where the product is the outer product to form a 3x3 matrix.
    double term2[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            term2[i][j] = k2 * diffP[i] * X[j];
        }
    }
    
    // 10. Compute R_dot = A - alpha * (term1 + term2)
    double R_dot[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_dot[i][j] = A[i][j] - alpha * (term1[i][j] + term2[i][j]);
        }
    }
    
    //------ Compute Omega_dot ------
    // 11. Compute I*Omega and store in I_Omega (3-vector)
    double I_Omega[3] = {0};
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += I[i][j] * Omega[j];
        }
        I_Omega[i] = sum;
    }
    
    // 12. Compute crossVec = cross(I*Omega, Omega)
    double crossVec[3];
    cross3(I_Omega, Omega, crossVec);
    
    // 13. Compute I_inv (the inverse of I)
    // double I_inv[3][3];
    // invertSymmetric3x3(I, I_inv);
    
    // 14. Compute vA = I_inv * crossVec (3-vector)
    double vA[3] = {0};
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += I_inv[i][j] * crossVec[j];
        }
        vA[i] = sum;
    }
    
    // 15. Compute termB = k1*(E(Omega, I) - E0)*(I*Omega)
    double E_val = E(Omega, I);
    double termB[3];
    for (int i = 0; i < 3; i++) {
        termB[i] = k1 * (E_val - E0) * I_Omega[i];
    }
    
    // 16. Compute R' (transpose of R)
    // double R_T[3][3];
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         R_T[i][j] = R[j][i];
    //     }
    // }
    
    // 17. Compute Y = R_T * diffP (a 3-vector)
    double Y[3] = {0};
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            // sum += R_T[i][j] * diffP[j];
            sum += R[j][i] * diffP[j];
        }
        Y[i] = sum;
    }
    
    // 18. Compute termC = k2 * (I * Y) (a 3-vector)
    double termC[3] = {0};
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += I[i][j] * Y[j];
        }
        termC[i] = k2 * sum;
    }
    
    // 19. Finally, compute Omega_dot = vA - alpha*(termB + termC)
    double Omega_dot[3];
    for (int i = 0; i < 3; i++) {
        Omega_dot[i] = vA[i] - alpha * (termB[i] + termC[i]);
    }
    
    // 20. Pack the results into the output 'out'
    // The first three columns are R_dot and the fourth column is Omega_dot.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = R_dot[i][j];
        }
        out[i][3] = Omega_dot[i];
    }
}

//------------------------------------------------------------------------------
// Euler integration with feedback control.
// Inputs:
//   - R0: initial 3x3 rotation matrix.
//   - Omega0: initial 3x1 angular velocity vector.
//   - I: 3x3 inertia matrix.
//   - k0, k1, k2: scalar feedback gains.
//   - E0: scalar energy reference.
//   - Pi0: 3x1 momentum reference vector.
//   - alpha: scalar feedback factor.
//   - tf: final time.
//   - h: time step.
// Outputs (allocated inside the function):
//   - R_out: pointer to a dynamically allocated array of size 9 * N (3x3 matrices stored consecutively).
//   - Omega_out: pointer to a dynamically allocated array of size 3 * N (each column is a 3-vector).
//   - t_out: pointer to an array of N time points.
//   - N: number of integration steps.
void euler_feedback(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    double k0, double k1, double k2,
                    double E0, const double Pi0[3],
                    double alpha,
                    double tf, double h,
                    double *&R_out, double *&Omega_out, double *&t_out,
                    int &N)
{
    const long long N_total = static_cast<long long>(std::ceil(tf / h)) + 1;

    const int stride = 1000;
    const long long approx = 1 + ((N_total - 1) + stride - 1) / stride;

    R_out     = new double[9 * approx];   // 3x3 블록 * approx
    Omega_out = new double[3 * approx];   // 3 * approx
    t_out     = new double[approx];

    // I^{-1}
    double I_inv[3][3];
    invertSymmetric3x3(I, I_inv);

    double current_state[3][4];
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) current_state[r][c] = R0[r][c];
        current_state[r][3] = Omega0[r];
    }
    double current_time = 0.0;

    long long k = 0;

    for (int r = 0; r < 3; ++r) {
        R_out[k*9 + r*3 + 0] = current_state[r][0];
        R_out[k*9 + r*3 + 1] = current_state[r][1];
        R_out[k*9 + r*3 + 2] = current_state[r][2];
        Omega_out[k*3 + r]   = current_state[r][3];
    }
    t_out[k] = current_time;
    ++k;

    for (long long i = 0; i < N_total - 1; ++i) {
        double dstate[3][4];
        dynamics_feedback(current_state, I, I_inv, k0, k1, k2, E0, Pi0, alpha, dstate);

        double next_state[3][4];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 4; ++c)
                next_state[r][c] = current_state[r][c] + h * dstate[r][c];

        current_time += h;

        const bool periodic_save = ((i + 1) % stride == 0);
        const bool is_last       = (i == N_total - 2);

        if (periodic_save || (is_last && !periodic_save)) {
            for (int r = 0; r < 3; ++r) {
                R_out[k*9 + r*3 + 0] = next_state[r][0];
                R_out[k*9 + r*3 + 1] = next_state[r][1];
                R_out[k*9 + r*3 + 2] = next_state[r][2];
                Omega_out[k*3 + r]   = next_state[r][3];
            }
            t_out[k] = current_time;
            ++k;
        }

        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 4; ++c)
                current_state[r][c] = next_state[r][c];
    }
    N = static_cast<int>(k);
}

double euler_feedback(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    double k0, double k1, double k2,
                    double E0, const double Pi0[3],
                    double alpha,
                    double tf, double h)
{
    // Calculate the number of steps: N = ceil(tf/h) + 1.
    int N = static_cast<int>(std::ceil(tf/h)) + 1;

    double maxError = 0.0;
    double currentError = 0.0;

    double I_inv[3][3];
    invertSymmetric3x3(I, I_inv);
    
    // Combine the initial rotation R0 and angular velocity Omega0 into current_state (3x4 array).
    double current_state[3][4];
    for (int i = 0; i < 3; i++) {
        // Columns 0-2: rotation matrix R0.
        for (int j = 0; j < 3; j++) {
            current_state[i][j] = R0[i][j];
        }
        // Column 3: angular velocity Omega0.
        current_state[i][3] = Omega0[i];
    }
    
    double current_time = 0.0;
    
    // Euler integration loop.
    for (int i = 0; i < N; i++) {
        // Compute the state derivative using the feedback dynamics.
        double dstate[3][4];
        dynamics_feedback(current_state, I, I_inv, k0, k1, k2, E0, Pi0, alpha, dstate);
        
        // Compute next_state = current_state + h * dstate.
        double next_state[3][4];
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 4; c++) {
                next_state[r][c] = current_state[r][c] + h * dstate[r][c];
            }
        }
        
        // Update time.
        current_time += h;
        
        // Update current_state for the next iteration.
        for (int r = 0; r < 3; r++) 
            for (int c = 0; c < 4; c++) 
                current_state[r][c] = next_state[r][c];
        
        currentError = getError(current_state, I, E0, Pi0, k0, k1, k2);
        if (currentError >= maxError)
            maxError = currentError;
    }
    return maxError;
}