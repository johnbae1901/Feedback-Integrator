#pragma once
#include <string>

// Save trajectories into a single CSV with header:
// t, omega_x, omega_y, omega_z, R00,R01,R02,R10,R11,R12,R20,R21,R22
// Arrays are expected to have lengths: t[N], Omega[3*N], R[9*N]
// R is row-major 3x3 per step: [R00 R01 R02 R10 R11 R12 R20 R21 R22] for each row.
bool save_state_csv(const std::string& csv_path,
                    const double* t, const double* Omega, const double* R,
                    int N);
