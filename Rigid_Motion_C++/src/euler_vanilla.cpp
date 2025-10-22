#include <cmath>      // for std::ceil, std::sqrt
#include <vector>
#include <tuple>
#include <array>
#include <iostream>
#include "../include/basic_operations.h"
#include "../include/getError.h"

using namespace std;

// Compute the derivatives of the state (3x4):
//   [ R_dot, Omega_dot ] = [ R * getSkew(Omega), inv(I)*cross(I*Omega, Omega) ]
// Here, 'state' is a 3x4 array with columns 0-2 for R and column 3 for Omega.
void dynamics_euler(const double state[3][4], const double I[3][3], double out[3][4]) {
    double R[3][3];
    double Omega[3];
    // Extract R and Omega from state.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R[i][j] = state[i][j];
        }
        Omega[i] = state[i][3];
    }
    
    // Compute OmegaHat = getSkew(Omega)
    double OmegaHat[3][3];
    getSkew(Omega, OmegaHat);
    
    // Compute R_dot = R * OmegaHat.
    double R_dot[3][3] = {0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            R_dot[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
                R_dot[i][j] += R[i][k] * OmegaHat[k][j];
        }
    
    // Compute I * Omega.
    double I_Omega[3] = {0};
    for (int i = 0; i < 3; i++) {
        I_Omega[i] = 0.0;
        for (int j = 0; j < 3; j++)
            I_Omega[i] += I[i][j] * Omega[j];
    }
    
    // Compute cross(I*Omega, Omega).
    double crossVec[3];
    cross3(I_Omega, Omega, crossVec);
    
    // Solve for Omega_dot: Omega_dot = inv(I) * crossVec.
    double I_inv[3][3];
    invert3x3(I, I_inv);
    double Omega_dot[3] = {0};
    for (int i = 0; i < 3; i++) {
        Omega_dot[i] = 0.0;
        for (int j = 0; j < 3; j++)
            Omega_dot[i] += I_inv[i][j] * crossVec[j];
    }
    
    // Pack the derivatives into out (first three columns: R_dot, fourth column: Omega_dot).
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            out[i][j] = R_dot[i][j];
        out[i][3] = Omega_dot[i];
    }
}

// Perform Euler integration. The inputs are initial R0 (3x3), Omega0 (3),
// the inertia matrix I (3x3), final time tf, and time step h.
// The outputs are pointers to arrays R_out (3x3xN stored in a 1D array of size 9*N),
// Omega_out (3xN stored in a 1D array of size 3*N), and t_out (time array of length N).
// The total number of time steps is returned in N.
void euler_vanilla(const double R0[3][3],
                   const double Omega0[3],
                   const double I[3][3],
                   double tf, double h,
                   double *&R_out, double *&Omega_out, double *&t_out,
                   int &N)
{
    const long long N_total = static_cast<long long>(std::ceil(tf / h)) + 1;

    const int stride = 1000;
    const long long approx = 1 + ((N_total - 1) + stride - 1) / stride;  // = 1 + ceil((N_total-1)/stride)

    R_out     = new double[9 * approx];
    Omega_out = new double[3 * approx];
    t_out     = new double[approx];

    double current_state[3][4];
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) current_state[r][c] = R0[r][c];
        current_state[r][3] = Omega0[r];
    }
    double current_time = 0.0;

    // 저장 인덱스
    long long k = 0;

    // [초기 샘플] t=0 저장
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
        dynamics_euler(current_state, I, dstate);

        double next_state[3][4];
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 4; ++c) {
                next_state[r][c] = current_state[r][c] + h * dstate[r][c];
            }
        }
        current_time += h;

        const bool periodic_save = ((i + 1) % stride == 0);
        const bool is_last       = (i == N_total - 2);

        if (periodic_save || (is_last && !periodic_save)) {
            // next_state를 슬롯 k에 기록
            for (int r = 0; r < 3; ++r) {
                R_out[k*9 + r*3 + 0] = next_state[r][0];
                R_out[k*9 + r*3 + 1] = next_state[r][1];
                R_out[k*9 + r*3 + 2] = next_state[r][2];
                Omega_out[k*3 + r]   = next_state[r][3];
            }
            t_out[k] = current_time;
            ++k;
        }

        // 다음 루프 준비
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 4; ++c)
                current_state[r][c] = next_state[r][c];
    }
    N = static_cast<int>(k);
}

std::tuple<
    double, 
    double, 
    double, 
    double> 
    euler_vanilla(const double R0[3][3],
                    const double Omega0[3],
                    const double I[3][3],
                    const double E0, 
                    const double Pi0[3], 
                    const double k0,
                    const double k1,
                    const double k2,
                    double tf, double h)
{
    // Determine the number of time steps. (Equivalent to: N = ceil(tf/h) + 1)
    int N = static_cast<int>(std::ceil(tf/h)) + 1;

    double maxV = 0.0;
    double maxdE = 0.0;
    double maxdPi_sq = 0.0;
    double maxdDet_sq = 0.0;
    ErrorReport currentError;
    
    // Set up the initial state as a 3x4 array.
    // The first three columns are the initial rotation matrix R0, and
    // the fourth column is the initial angular velocity Omega0.
    double current_state[3][4];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            current_state[i][j] = R0[i][j];
        current_state[i][3] = Omega0[i];
    }
    
    double current_time = 0.0;
    
    // Euler integration loop.
    for (int i = 0; i < N; i++) {
        // Compute the state derivative dstate = dynamicsEuler(current_state, I).
        double dstate[3][4];
        dynamics_euler(current_state, I, dstate);
        
        // Compute next_state = current_state + h * dstate.
        double next_state[3][4];
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 4; c++)
                next_state[r][c] = current_state[r][c] + h * dstate[r][c];
        }
        
        // Update current time.
        current_time += h;

        // Prepare for next iteration.
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
    return {maxV, maxdE, sqrt(maxdPi_sq), sqrt(maxdDet_sq)};
}