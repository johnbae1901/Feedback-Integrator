#include <vector>
#include <array>
#include <cmath>
#include "basic_operations.h"
#include "getError.h"

static double dist3_squared(const std::vector<double>& a,
                            const std::vector<double>& b)
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return dx*dx + dy*dy + dz*dz;
}

static double dist3_squared(const std::array<double,3>& a,
                            const std::array<double,3>& b)
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return dx*dx + dy*dy + dz*dz;
}

// ===== perturbed Kepler 전용 (권장): error = kL*||L-L0||^2 + kE*(E-E0)^2 =====
std::vector<double> getError_PK(
    const std::vector<double>& x1,
    const std::vector<double>& x2,
    const std::vector<double>& x3,
    const std::vector<double>& v1,
    const std::vector<double>& v2,
    const std::vector<double>& v3,
    const std::array<double,3>& L0,
    double E0,
    double kL,
    double kE,
    const PKParams& P)
{
    size_t N = x1.size();
    std::vector<double> error(N, 0.0);
    for(size_t i=0; i<N; ++i){
        double r[3] = { x1[i], x2[i], x3[i] };
        double v[3] = { v1[i], v2[i], v3[i] };
        double L[3]; cross3(r, v, L);
        std::array<double,3> Lvec = {L[0], L[1], L[2]};
        double dL2 = dist3_squared(Lvec, L0);
        double dE  = energy_pk(r, v, P) - E0;
        error[i] = kL * dL2 + kE * (dE*dE);
    }
    return error;
}

LEError getError_PK(
    const double& x1,
    const double& x2,
    const double& x3,
    const double& v1,
    const double& v2,
    const double& v3,
    const std::array<double,3>& L0,
    double E0,
    double kL,
    double kE,
    const PKParams& P)
{
    double r[3] = { x1, x2, x3 };
    double v[3] = { v1, v2, v3 };
    double L[3]; cross3(r, v, L);
    std::array<double,3> Lvec = {L[0], L[1], L[2]};
    double dL2 = dist3_squared(Lvec, L0);
    double dE  = energy_pk(r, v, P) - E0;

    LEError out;
    out.distL_sq = dL2;
    out.distE_sq = dE * dE;
    out.error = kL * dL2 + kE * (dE*dE);
    return out;
}
