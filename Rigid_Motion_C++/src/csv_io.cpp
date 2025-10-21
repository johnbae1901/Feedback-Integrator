#include "../include/csv_io.h"
#include <fstream>
#include <filesystem>

bool save_state_csv(const std::string& csv_path,
                    const double* t, const double* Omega, const double* R,
                    int N)
{
    namespace fs = std::filesystem;
    fs::path p(csv_path);
    fs::create_directories(p.parent_path());

    std::ofstream f(csv_path);
    if (!f.is_open()) return false;

    f << "t,omega_x,omega_y,omega_z,R00,R01,R02,R10,R11,R12,R20,R21,R22\n";
    for (int k = 0; k < N; ++k) {
        f << t[k] << ","
          << Omega[3*k+0] << "," << Omega[3*k+1] << "," << Omega[3*k+2] << ",";
        const int rbase = 9*k;
        f << R[rbase+0] << "," << R[rbase+1] << "," << R[rbase+2] << ","
          << R[rbase+3] << "," << R[rbase+4] << "," << R[rbase+5] << ","
          << R[rbase+6] << "," << R[rbase+7] << "," << R[rbase+8] << "\n";
    }
    return true;
}
