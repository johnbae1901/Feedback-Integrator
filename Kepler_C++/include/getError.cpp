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

LAError getError(
    const double& x1,
    const double& x2,
    const double& x3,
    const double& v1,
    const double& v2,
    const double& v3,
    const std::array<double,3>& L0,
    const std::array<double,3>& A0,
    double k1,
    double k2,
    double mu)
{
    double r[3] = { x1, x2, x3 };
    double v[3] = { v1, v2, v3 };

    // L = r × v
    double L[3];
    cross3(r, v, L);

    // A = v × L − μ r / ||r||
    double temp[3];
    cross3(v, L, temp);
    double rNorm = norm3(r);
    double A[3] = {
        temp[0] - mu * r[0] / rNorm,
        temp[1] - mu * r[1] / rNorm,
        temp[2] - mu * r[2] / rNorm
    };

    std::array<double,3> Lvec = {L[0], L[1], L[2]};
    std::array<double,3> Avec = {A[0], A[1], A[2]};

    double distL_sq = dist3_squared(Lvec, L0);
    double distA_sq = dist3_squared(Avec, A0);

    LAError out;
    out.distL_sq = distL_sq;
    out.distA_sq = distA_sq;
    out.error = k1*distL_sq + k2*distA_sq;
    return out;
}