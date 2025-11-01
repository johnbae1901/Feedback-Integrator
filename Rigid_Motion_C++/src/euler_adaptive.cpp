#include <cmath>      // for std::ceil, std::sqrt
#include <vector>
#include <tuple>
#include <array>
#include <iostream>
#include "../include/basic_operations.h"
#include "../include/getError.h"

using namespace std;

double estimateL(const double R[3][3], const double Omega[3])
{
    // Extract R's elements
    double r1_1 = R[0][0], r1_2 = R[0][1], r1_3 = R[0][2];
    double r2_1 = R[1][0], r2_2 = R[1][1], r2_3 = R[1][2];
    double r3_1 = R[2][0], r3_2 = R[2][1], r3_3 = R[2][2];
    // Extract Omega's components:
    double w1 = Omega[0], w2 = Omega[1], w3 = Omega[2];

    // Build the 12×12 matrix H.
    double H[12][12];

    // Row 0
    H[0][0] = 150*r1_1*r1_1 + 50*r1_2*r1_2 + 50*r1_3*r1_3 + 50*r2_1*r2_1 + 50*r3_1*r3_1 + 450*w1*w1 - 50;
    H[0][1] = 100*r1_1*r1_2 + 50*r2_1*r2_2 + 50*r3_1*r3_2 + 300*w1*w2;
    H[0][2] = 100*r1_1*r1_3 + 50*r2_1*r2_3 + 50*r3_1*r3_3 + 150*w1*w3;
    H[0][3] = 100*r1_1*r2_1 + 50*r1_2*r2_2 + 50*r1_3*r2_3;
    H[0][4] = 50*r1_2*r2_1;
    H[0][5] = 50*r1_3*r2_1;
    H[0][6] = 100*r1_1*r3_1 + 50*r1_2*r3_2 + 50*r1_3*r3_3;
    H[0][7] = 50*r1_2*r3_1;
    H[0][8] = 50*r1_3*r3_1;
    H[0][9] = 900*r1_1*w1 + 300*r1_2*w2 + 150*r1_3*w3 - 450;
    H[0][10] = 300*r1_2*w1;
    H[0][11] = 150*r1_3*w1;

    // Row 1
    H[1][0]  = 100*r1_1*r1_2 + 50*r2_1*r2_2 + 50*r3_1*r3_2 + 300*w1*w2;
    H[1][1]  = 50*r1_1*r1_1 + 150*r1_2*r1_2 + 50*r1_3*r1_3 + 50*r2_2*r2_2 + 50*r3_2*r3_2 + 200*w2*w2 - 50;
    H[1][2]  = 100*r1_2*r1_3 + 50*r2_2*r2_3 + 50*r3_2*r3_3 + 100*w2*w3;
    H[1][3]  = 50*r1_1*r2_2;
    H[1][4]  = 50*r1_1*r2_1 + 100*r1_2*r2_2 + 50*r1_3*r2_3;
    H[1][5]  = 50*r1_3*r2_2;
    H[1][6]  = 50*r1_1*r3_2;
    H[1][7]  = 50*r1_1*r3_1 + 100*r1_2*r3_2 + 50*r1_3*r3_3;
    H[1][8]  = 50*r1_3*r3_2;
    H[1][9]  = 300*r1_1*w2;
    H[1][10] = 300*r1_1*w1 + 400*r1_2*w2 + 100*r1_3*w3 - 300;
    H[1][11] = 100*r1_3*w2;
    
    // Row 2
    H[2][0]  = 100*r1_1*r1_3 + 50*r2_1*r2_3 + 50*r3_1*r3_3 + 150*w1*w3;
    H[2][1]  = 100*r1_2*r1_3 + 50*r2_2*r2_3 + 50*r3_2*r3_3 + 100*w2*w3;
    H[2][2]  = 50*r1_1*r1_1 + 50*r1_2*r1_2 + 150*r1_3*r1_3 + 50*r2_3*r2_3 + 50*r3_3*r3_3 + 50*w3*w3 - 50;
    H[2][3]  = 50*r1_1*r2_3;
    H[2][4]  = 50*r1_2*r2_3;
    H[2][5]  = 50*r1_1*r2_1 + 50*r1_2*r2_2 + 100*r1_3*r2_3;
    H[2][6]  = 50*r1_1*r3_3;
    H[2][7]  = 50*r1_2*r3_3;
    H[2][8]  = 50*r1_1*r3_1 + 50*r1_2*r3_2 + 100*r1_3*r3_3;
    H[2][9]  = 150*r1_1*w3;
    H[2][10] = 100*r1_2*w3;
    H[2][11] = 150*r1_1*w1 + 100*r1_2*w2 + 100*r1_3*w3 - 150;
    
    // Row 3
    H[3][0]  = 100*r1_1*r2_1 + 50*r1_2*r2_2 + 50*r1_3*r2_3;
    H[3][1]  = 50*r1_1*r2_2;
    H[3][2]  = 50*r1_1*r2_3;
    H[3][3]  = 50*r1_1*r1_1 + 150*r2_1*r2_1 + 50*r2_2*r2_2 + 50*r2_3*r2_3 + 50*r3_1*r3_1 + 450*w1*w1 - 50;
    H[3][4]  = 50*r1_1*r1_2 + 100*r2_1*r2_2 + 50*r3_1*r3_2 + 300*w1*w2;
    H[3][5]  = 50*r1_1*r1_3 + 100*r2_1*r2_3 + 50*r3_1*r3_3 + 150*w1*w3;
    H[3][6]  = 100*r2_1*r3_1 + 50*r2_2*r3_2 + 50*r2_3*r3_3;
    H[3][7]  = 50*r2_2*r3_1;
    H[3][8]  = 50*r2_3*r3_1;
    H[3][9]  = 900*r2_1*w1 + 300*r2_2*w2 + 150*r2_3*w3 - 300;
    H[3][10] = 300*r2_2*w1;
    H[3][11] = 150*r2_3*w1;
    
    // Row 4
    H[4][0]  = 50*r1_2*r2_1;
    H[4][1]  = 50*r1_1*r2_1 + 100*r1_2*r2_2 + 50*r1_3*r2_3;
    H[4][2]  = 50*r1_2*r2_3;
    H[4][3]  = 50*r1_1*r1_2 + 100*r2_1*r2_2 + 50*r3_1*r3_2 + 300*w1*w2;
    H[4][4]  = 50*r1_2*r1_2 + 50*r2_1*r2_1 + 150*r2_2*r2_2 + 50*r2_3*r2_3 + 50*r3_2*r3_2 + 200*w2*w2 - 50;
    H[4][5]  = 50*r1_2*r1_3 + 100*r2_2*r2_3 + 50*r3_2*r3_3 + 100*w2*w3;
    H[4][6]  = 50*r2_1*r3_2;
    H[4][7]  = 50*r2_1*r3_1 + 100*r2_2*r3_2 + 50*r2_3*r3_3;
    H[4][8]  = 50*r2_3*r3_2;
    H[4][9]  = 300*r2_1*w2;
    H[4][10] = 300*r2_1*w1 + 400*r2_2*w2 + 100*r2_3*w3 - 200;
    H[4][11] = 100*r2_3*w2;
    
    // Row 5
    H[5][0]  = 50*r1_3*r2_1;
    H[5][1]  = 50*r1_3*r2_2;
    H[5][2]  = 50*r1_1*r2_1 + 50*r1_2*r2_2 + 100*r1_3*r2_3;
    H[5][3]  = 50*r1_1*r2_3 + 100*r2_1*r2_3 + 50*r3_1*r3_3 + 150*w1*w3;
    H[5][4]  = 50*r1_2*r1_3 + 100*r2_2*r2_3 + 50*r3_2*r3_3 + 100*w2*w3;
    H[5][5]  = 50*r1_3*r1_3 + 50*r2_1*r2_1 + 50*r2_2*r2_2 + 150*r2_3*r2_3 + 50*r3_3*r3_3 + 50*w3*w3 - 50;
    H[5][6]  = 50*r2_1*r3_3;
    H[5][7]  = 50*r2_2*r3_3;
    H[5][8]  = 50*r2_1*r3_1 + 50*r2_2*r3_2 + 100*r2_3*r3_3;
    H[5][9]  = 150*r2_1*w3;
    H[5][10] = 100*r2_2*w3;
    H[5][11] = 150*r2_1*w1 + 100*r2_2*w2 + 100*r2_3*w3 - 100;
    
    // Row 6
    H[6][0]  = 100*r1_1*r3_1 + 50*r1_2*r3_2 + 50*r1_3*r3_3;
    H[6][1]  = 50*r1_1*r3_2;
    H[6][2]  = 50*r1_1*r3_3;
    H[6][3]  = 100*r2_1*r3_1 + 50*r2_2*r3_2 + 50*r2_3*r3_3;
    H[6][4]  = 50*r2_1*r3_2;
    H[6][5]  = 50*r2_1*r3_3;
    H[6][6]  = 50*r1_1*r1_1 + 50*r2_1*r2_1 + 150*r3_1*r3_1 + 50*r3_2*r3_2 + 50*r3_3*r3_3 + 450*w1*w1 - 50;
    H[6][7]  = 50*r1_1*r1_2 + 50*r2_1*r2_2 + 100*r3_1*r3_2 + 300*w1*w2;
    H[6][8]  = 50*r1_1*r1_3 + 50*r2_1*r2_3 + 100*r3_1*r3_3 + 150*w1*w3;
    H[6][9]  = 900*r3_1*w1 + 300*r3_2*w2 + 150*r3_3*w3 - 150;
    H[6][10] = 300*r3_2*w1;
    H[6][11] = 150*r3_3*w1;
    
    // Row 7
    H[7][0]  = 50*r1_2*r3_1;
    H[7][1]  = 50*r1_1*r3_1 + 100*r1_2*r3_2 + 50*r1_3*r3_3;
    H[7][2]  = 50*r1_2*r3_3;
    H[7][3]  = 50*r2_2*r3_1;
    H[7][4]  = 50*r2_1*r3_1 + 100*r2_2*r3_2 + 50*r2_3*r3_3;
    H[7][5]  = 50*r2_2*r3_3;
    H[7][6]  = 50*r1_1*r1_2 + 50*r2_1*r2_2 + 100*r3_1*r3_2 + 300*w1*w2;
    H[7][7]  = 50*r1_2*r1_2 + 50*r2_2*r2_2 + 50*r3_1*r3_1 + 150*r3_2*r3_2 + 50*r3_3*r3_3 + 200*w2*w2 - 50;
    H[7][8]  = 50*r1_2*r1_3 + 50*r2_2*r2_3 + 100*r3_2*r3_3 + 100*w2*w3;
    H[7][9]  = 300*r3_1*w2;
    H[7][10] = 300*r3_1*w1 + 400*r3_2*w2 + 100*r3_3*w3 - 100;
    H[7][11] = 100*r3_3*w2;
    
    // Row 8
    H[8][0]  = 50*r1_3*r3_1;
    H[8][1]  = 50*r1_3*r3_2;
    H[8][2]  = 50*r1_1*r3_1 + 50*r1_2*r3_2 + 100*r1_3*r3_3;
    H[8][3]  = 50*r2_3*r3_1;
    H[8][4]  = 50*r2_3*r3_2;
    H[8][5]  = 50*r2_1*r3_1 + 50*r2_2*r3_2 + 100*r2_3*r3_3;
    H[8][6]  = 50*r1_1*r1_3 + 50*r2_1*r2_3 + 100*r3_1*r3_3 + 150*w1*w3;
    H[8][7]  = 50*r1_2*r1_3 + 50*r2_2*r2_3 + 100*r3_2*r3_3 + 100*w2*w3;
    H[8][8]  = 50*r1_3*r1_3 + 50*r2_3*r2_3 + 50*r3_1*r3_1 + 50*r3_2*r3_2 + 150*r3_3*r3_3 + 50*w3*w3 - 300;
    H[8][9]  = 150*r3_1*w3;
    H[8][10] = 100*r3_2*w3;
    H[8][11] = 150*r3_1*w1 + 100*r3_2*w2 + 100*r3_3*w3 - 50;
    
    // Row 9
    H[9][0]  = 900*r1_1*w1 + 300*r1_2*w2 + 150*r1_3*w3 - 450;
    H[9][1]  = 300*r1_1*w2;
    H[9][2]  = 150*r1_1*w3;
    H[9][3]  = 900*r2_1*w1 + 300*r2_2*w2 + 150*r2_3*w3 - 300;
    H[9][4]  = 300*r2_1*w2;
    H[9][5]  = 150*r2_1*w3;
    H[9][6]  = 900*r3_1*w1 + 300*r3_2*w2 + 150*r3_3*w3 - 150;
    H[9][7]  = 300*r3_1*w2;
    H[9][8]  = 150*r3_1*w3;
    H[9][9]  = 450*r1_1*r1_1 + 450*r2_1*r2_1 + 450*r3_1*r3_1 + 1350*w1*w1 + 300*w2*w2 + 150*w3*w3 - 900;
    H[9][10] = 300*r1_1*r1_2 + 300*r2_1*r2_2 + 300*r3_1*r3_2 + 600*w1*w2;
    H[9][11] = 150*r1_1*r1_3 + 150*r2_1*r2_3 + 150*r3_1*r3_3 + 300*w1*w3;
    
    // Row 10
    H[10][0]  = 300*r1_2*w1;
    H[10][1]  = 300*r1_1*w1 + 400*r1_2*w2 + 100*r1_3*w3 - 300;
    H[10][2]  = 100*r1_2*w3;
    H[10][3]  = 300*r2_2*w1;
    H[10][4]  = 300*r2_1*w1 + 400*r2_2*w2 + 100*r2_3*w3 - 200;
    H[10][5]  = 100*r2_2*w3;
    H[10][6]  = 300*r3_2*w1;
    H[10][7]  = 300*r3_1*w1 + 400*r3_2*w2 + 100*r3_3*w3 - 100;
    H[10][8]  = 100*r3_2*w3;
    H[10][9]  = 300*r1_1*r1_2 + 300*r2_1*r2_2 + 300*r3_1*r3_2 + 600*w1*w2;
    H[10][10] = 200*r1_2*r1_2 + 200*r2_2*r2_2 + 200*r3_2*r3_2 + 300*w1*w1 + 600*w2*w2 + 100*w3*w3 - 600;
    H[10][11] = 100*r1_2*r1_3 + 100*r2_2*r2_3 + 100*r3_2*r3_3 + 200*w2*w3;
    
    // Row 11
    H[11][0]  = 150*r1_3*w1;
    H[11][1]  = 100*r1_3*w2;
    H[11][2]  = 150*r1_1*w1 + 100*r1_2*w2 + 100*r1_3*w3 - 150;
    H[11][3]  = 150*r2_3*w1;
    H[11][4]  = 100*r2_3*w2;
    H[11][5]  = 150*r2_1*w1 + 100*r2_2*w2 + 100*r2_3*w3 - 100;
    H[11][6]  = 150*r3_3*w1;
    H[11][7]  = 100*r3_3*w2;
    H[11][8]  = 150*r3_1*w1 + 100*r3_2*w2 + 100*r3_3*w3 - 50;
    H[11][9]  = 150*r1_1*r1_3 + 150*r2_1*r2_3 + 150*r3_1*r3_3 + 300*w1*w3;
    H[11][10] = 100*r1_2*r1_3 + 100*r2_2*r2_3 + 100*r3_2*r3_3 + 200*w2*w3;
    H[11][11] = 50*r1_3*r1_3 + 50*r2_3*r2_3 + 50*r3_3*r3_3 + 150*w1*w1 + 100*w2*w2 + 150*w3*w3 - 300;

    double output = 0;
    for (int i = 0; i < 11; i++) {
        for (int j = 0; j < 11; j++) {
            output += H[i][j] * H[i][j];
        }
    }
    
    return sqrt(output);
}

//------------------------------------------------------------------------------
// dynamics_adaptive: Compute adaptive feedback dynamics.
// Inputs:
//   - state: a 3x4 array where columns 0-2 are rotation matrix R and column 3 is Omega.
//   - I: 3x3 inertia matrix.
//   - k0, k1, k2: scalar feedback gains.
//   - E0: scalar reference energy.
//   - Pi0: 3-vector reference momentum.
//   - h: time step (used here in the adaptive gain computation).
// Output:
//   - out: a 3x4 array containing the derivatives [R_dot, Omega_dot].
//------------------------------------------------------------------------------
void dynamics_adaptive(const double state[3][4],
                         const double I[3][3],
                         double k0, double k1, double k2,
                         double E0, const double Pi0[3],
                         double h,
                         double out[3][4],
                         const double Hmin)
{
    // 1. Extract R (3x3) and Omega (3x1) from state.
    double R[3][3];
    double Omega[3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R[i][j] = state[i][j];
        }
        Omega[i] = state[i][3];
    }
    
    // 2. Compute the skew-symmetric matrix OmegaHat from Omega.
    double OmegaHat[3][3];
    getSkew(Omega, OmegaHat);
    
    // 3. Estimate L and compute adaptive gain alpha = 1 / (h * L)
    double L = estimateL(R, Omega);
    // double L = 1986;
    double alpha = 1.0 / (h * max(Hmin, L));
    
    // 4. Compute A = R * OmegaHat  (a 3x3 matrix)
    double A[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                A[i][j] += R[i][k] * OmegaHat[k][j];
            }
        }
    }
    
    // 5. Compute M = R' * R - I
    //    R' is the transpose of R so that M(i,j) = sum_k(R[k][i]*R[k][j]) - delta(i,j)
    double M[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += R[k][i] * R[k][j];
            }
            M[i][j] = sum - ((i == j) ? 1.0 : 0.0);
        }
    }
    
    // 6. Compute term1 = k0 * R * M.
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
    
    // 7. Compute P = Pi(R, Omega, I)
    double P[3];
    Pi(R, Omega, I, P);
    
    // 8. Compute diffP = P - Pi0 (a 3-vector)
    double diffP[3];
    for (int i = 0; i < 3; i++) {
        diffP[i] = P[i] - Pi0[i];
    }
    
    // 9. Compute X = Omega' * I (i.e. a row vector; here stored as a 3-vector):
    double X[3] = {0, 0, 0};
    for (int j = 0; j < 3; j++) {
        double sum = 0.0;
        for (int i = 0; i < 3; i++) {
            sum += Omega[i] * I[i][j];
        }
        X[j] = sum;
    }
    
    // 10. Compute term2 = k2 * outer(diffP, X): term2(i,j) = k2 * diffP[i] * X[j]
    double term2[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            term2[i][j] = k2 * diffP[i] * X[j];
        }
    }
    
    // 11. Compute the derivative of R: R_dot = A - alpha * (term1 + term2)
    double R_dot[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_dot[i][j] = A[i][j] - alpha * (term1[i][j] + term2[i][j]);
        }
    }
    
    //------ Compute Omega_dot ------
    // 12. Compute I*Omega and store in I_Omega.
    double I_Omega[3] = {0};
    for (int i = 0; i < 3; i++) {
        I_Omega[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            I_Omega[i] += I[i][j] * Omega[j];
        }
    }
    
    // 13. Compute crossVec = cross(I_Omega, Omega).
    double crossVec[3];
    cross3(I_Omega, Omega, crossVec);
    
    // 14. Compute I_inv = inverse(I) (using a general 3x3 inversion).
    double I_inv[3][3];
    invert3x3(I, I_inv);
    
    // 15. Compute vA = I_inv * crossVec.
    double vA[3] = {0};
    for (int i = 0; i < 3; i++) {
        vA[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            vA[i] += I_inv[i][j] * crossVec[j];
        }
    }
    
    // 16. Compute E_val = E(Omega, I).
    double E_val = E(Omega, I);
    
    // 17. Compute termB = k1*(E_val - E0)*(I*Omega) (elementwise scaling).
    double termB[3] = {0};
    for (int i = 0; i < 3; i++) {
        termB[i] = k1 * (E_val - E0) * I_Omega[i];
    }
    
    // 18. Compute R_T, the transpose of R.
    double R_T[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_T[i][j] = R[j][i];
        }
    }
    
    // 19. Compute Y = R_T * diffP.
    double Y[3] = {0};
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += R_T[i][j] * diffP[j];
        }
        Y[i] = sum;
    }
    
    // 20. Compute termC = k2 * (I * Y).
    double termC[3] = {0};
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += I[i][j] * Y[j];
        }
        termC[i] = k2 * sum;
    }
    
    // 21. Compute Omega_dot = vA - alpha*(termB + termC).
    double Omega_dot[3] = {0};
    for (int i = 0; i < 3; i++) {
        Omega_dot[i] = vA[i] - alpha * (termB[i] + termC[i]);
    }
    
    // 22. Pack the results into the output.
    //     The first three columns correspond to R_dot and the fourth column to Omega_dot.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = R_dot[i][j];
        }
        out[i][3] = Omega_dot[i];
    }
}

void dynamics_adaptive_light(const double state[3][4],
                         const double I[3][3],
                         double k0, double k1, double k2,
                         double E0, const double Pi0[3],
                         double h, double lipConstant,
                         double out[3][4],
                         const double Hmin)
{
    // 1. Extract R (3x3) and Omega (3x1) from state.
    double R[3][3];
    double Omega[3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R[i][j] = state[i][j];
        }
        Omega[i] = state[i][3];
    }
    
    // 2. Compute the skew-symmetric matrix OmegaHat from Omega.
    double OmegaHat[3][3];
    getSkew(Omega, OmegaHat);
    
    // 3. Estimate L and compute adaptive gain alpha = 1 / (h * L)
    double alpha = 1.0 / (h * max(Hmin, lipConstant));
    
    // 4. Compute A = R * OmegaHat  (a 3x3 matrix)
    double A[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                A[i][j] += R[i][k] * OmegaHat[k][j];
            }
        }
    }
    
    // 5. Compute M = R' * R - I
    //    R' is the transpose of R so that M(i,j) = sum_k(R[k][i]*R[k][j]) - delta(i,j)
    double M[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += R[k][i] * R[k][j];
            }
            M[i][j] = sum - ((i == j) ? 1.0 : 0.0);
        }
    }
    
    // 6. Compute term1 = k0 * R * M.
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
    
    // 7. Compute P = Pi(R, Omega, I)
    double P[3];
    Pi(R, Omega, I, P);
    
    // 8. Compute diffP = P - Pi0 (a 3-vector)
    double diffP[3];
    for (int i = 0; i < 3; i++) {
        diffP[i] = P[i] - Pi0[i];
    }
    
    // 9. Compute X = Omega' * I (i.e. a row vector; here stored as a 3-vector):
    double X[3] = {0, 0, 0};
    for (int j = 0; j < 3; j++) {
        double sum = 0.0;
        for (int i = 0; i < 3; i++) {
            sum += Omega[i] * I[i][j];
        }
        X[j] = sum;
    }
    
    // 10. Compute term2 = k2 * outer(diffP, X): term2(i,j) = k2 * diffP[i] * X[j]
    double term2[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            term2[i][j] = k2 * diffP[i] * X[j];
        }
    }
    
    // 11. Compute the derivative of R: R_dot = A - alpha * (term1 + term2)
    double R_dot[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_dot[i][j] = A[i][j] - alpha * (term1[i][j] + term2[i][j]);
        }
    }
    
    //------ Compute Omega_dot ------
    // 12. Compute I*Omega and store in I_Omega.
    double I_Omega[3] = {0};
    for (int i = 0; i < 3; i++) {
        I_Omega[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            I_Omega[i] += I[i][j] * Omega[j];
        }
    }
    
    // 13. Compute crossVec = cross(I_Omega, Omega).
    double crossVec[3];
    cross3(I_Omega, Omega, crossVec);
    
    // 14. Compute I_inv = inverse(I) (using a general 3x3 inversion).
    double I_inv[3][3];
    invert3x3(I, I_inv);
    
    // 15. Compute vA = I_inv * crossVec.
    double vA[3] = {0};
    for (int i = 0; i < 3; i++) {
        vA[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            vA[i] += I_inv[i][j] * crossVec[j];
        }
    }
    
    // 16. Compute E_val = E(Omega, I).
    double E_val = E(Omega, I);
    
    // 17. Compute termB = k1*(E_val - E0)*(I*Omega) (elementwise scaling).
    double termB[3] = {0};
    for (int i = 0; i < 3; i++) {
        termB[i] = k1 * (E_val - E0) * I_Omega[i];
    }
    
    // 18. Compute R_T, the transpose of R.
    double R_T[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_T[i][j] = R[j][i];
        }
    }
    
    // 19. Compute Y = R_T * diffP.
    double Y[3] = {0};
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += R_T[i][j] * diffP[j];
        }
        Y[i] = sum;
    }
    
    // 20. Compute termC = k2 * (I * Y).
    double termC[3] = {0};
    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += I[i][j] * Y[j];
        }
        termC[i] = k2 * sum;
    }
    
    // 21. Compute Omega_dot = vA - alpha*(termB + termC).
    double Omega_dot[3] = {0};
    for (int i = 0; i < 3; i++) {
        Omega_dot[i] = vA[i] - alpha * (termB[i] + termC[i]);
    }
    
    // 22. Pack the results into the output.
    //     The first three columns correspond to R_dot and the fourth column to Omega_dot.
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
void euler_adaptive(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    double k0, double k1, double k2,
                    double E0, const double Pi0[3],
                    double tf, double h,
                    int m, double lambda,
                    double *&R_out, double *&Omega_out, double *&t_out,
                    int &N,
                    const double Hmin)
{
    const long long steps = static_cast<long long>(std::ceil(tf / h));
    const int stride = 1000;

    const long long approx = 1 + (steps + stride - 1) / stride;

    R_out     = new double[9 * approx];
    Omega_out = new double[3 * approx];
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
    double lipConstant = 0.0;

    long long ksave = 0;

    for (int r = 0; r < 3; ++r) {
        R_out[ksave*9 + r*3 + 0] = current_state[r][0];
        R_out[ksave*9 + r*3 + 1] = current_state[r][1];
        R_out[ksave*9 + r*3 + 2] = current_state[r][2];
        Omega_out[ksave*3 + r]   = current_state[r][3];
    }
    t_out[ksave] = current_time;
    ++ksave;

    for (long long s = 0; s < steps; ++s) {
        if ((s % m) == 0) {
            double Rtmp[3][3];
            double Otmp[3];
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) Rtmp[r][c] = current_state[r][c];
                Otmp[r] = current_state[r][3];
            }
            lipConstant = lambda * estimateL(Rtmp, Otmp);
        }

        double dstate[3][4];
        dynamics_adaptive_light(current_state, I, k0, k1, k2, E0, Pi0,
                                h, lipConstant, dstate, Hmin);

        double next_state[3][4];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 4; ++c)
                next_state[r][c] = current_state[r][c] + h * dstate[r][c];

        current_time += h;

        const bool periodic_save = ((s + 1) % stride == 0);
        const bool is_last       = (s == steps - 1);

        if (periodic_save || (is_last && !periodic_save)) {
            for (int r = 0; r < 3; ++r) {
                R_out[ksave*9 + r*3 + 0] = next_state[r][0];
                R_out[ksave*9 + r*3 + 1] = next_state[r][1];
                R_out[ksave*9 + r*3 + 2] = next_state[r][2];
                Omega_out[ksave*3 + r]   = next_state[r][3];
            }
            t_out[ksave] = current_time;
            ++ksave;
        }

        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 4; ++c)
                current_state[r][c] = next_state[r][c];
    }
    N = static_cast<int>(ksave);
}

std::tuple<
    double, 
    double, 
    double, 
    double> 
    euler_adaptive(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    double k0, double k1, double k2,
                    double E0, const double Pi0[3],
                    double tf, double h,
                    int m, double lambda,
                    const double Hmin)
{
    // Calculate the number of steps: N = ceil(tf/h) + 1.
    int N = static_cast<int>(std::ceil(tf/h)) + 1;

    double maxV = 0.0;
    double maxdE = 0.0;
    double maxdPi_sq = 0.0;
    double maxdDet_sq = 0.0;
    ErrorReport currentError;

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
    
    // L(x)
    double lipConstant;
    // double minL = 1e5;
    // double maxL = 0.0;

    // Euler integration loop.
    for (int i = 0; i < N; i++) {
        if ((i % m) == 0)
        {
            // Extract R (3x3) and Omega (3x1) from state.
            double R[3][3];
            double Omega[3];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    R[i][j] = current_state[i][j];
                }
                Omega[i] = current_state[i][3];
            }
            // Estimate L and compute adaptive gain alpha = 1 / (h * L)
            lipConstant = lambda * estimateL(R, Omega);
            // if (lipConstant < minL) minL = lipConstant; 
            // if (lipConstant > maxL) maxL = lipConstant;
        }

        // Compute the state derivative using the feedback dynamics.
        double dstate[3][4];
        dynamics_adaptive_light(current_state, I, k0, k1, k2, E0, Pi0, h, lipConstant, dstate, Hmin);
        
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
        if (currentError.weighted_sum >= maxV)
            maxV = currentError.weighted_sum;
        if (currentError.abs_deltaE >= maxdE)
            maxdE = currentError.abs_deltaE;
        if (currentError.deltaPi_norm_sq >= maxdPi_sq)
            maxdPi_sq = currentError.deltaPi_norm_sq;
        if (currentError.frob_RT_R_minus_I_sq >= maxdDet_sq)
            maxdDet_sq = currentError.frob_RT_R_minus_I_sq;
    }
    // cout << minL / lambda << endl;
    // cout << maxL / lambda << endl;
    return {maxV, maxdE, sqrt(maxdPi_sq), sqrt(maxdDet_sq)};
}