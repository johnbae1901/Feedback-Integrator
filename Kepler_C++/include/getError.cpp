#include <vector>
#include <array>
#include <cmath>
#include "basic_operations.h"

/**
 * @brief Return the squared Euclidean distance between a and b => ||a - b||^2
 */
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

/**
 * @brief getError: replicates the MATLAB function
 *
 * In MATLAB:
 *     N = size(x,1);
 *     for n=1:N
 *        L = cross(x(n,1:3), x(n,4:6));
 *        A = cross(x(n,4:6), L) - mu*x(n,1:3)/norm(x(n,1:3));
 *        error(n) = k1*||L-L0||^2 + k2*||A-A0||^2;
 *     end
 *
 * Here:
 *   - xVec is a vector of length N, where each element is {x, y, z, vx, vy, vz}.
 *   - L0, A0 are each 3D arrays.
 */
std::vector<double> getError(
    const std::vector<double>& x1,
    const std::vector<double>& x2,
    const std::vector<double>& x3,
    const std::vector<double>& v1,
    const std::vector<double>& v2,
    const std::vector<double>& v3,
    const std::vector<double>& L0,
    const std::vector<double>& A0,
    double k1,
    double k2,
    double mu)
{
    // Assume x1.size() == x2.size() == ... == v3.size() for safety
    size_t N = x1.size();
    std::vector<double> error(N, 0.0);

    for(size_t i = 0; i < N; ++i)
    {
        // r = (x1[i], x2[i], x3[i])
        // v = (v1[i], v2[i], v3[i])
        double r[3] = { x1[i], x2[i], x3[i] };
        double v[3] = { v1[i], v2[i], v3[i] };

        // Compute L = cross(r, v)
        double L[3];
        cross3(r, v, L);

        // Compute A = cross(v, L) - mu * r / ||r||
        double temp[3];
        cross3(v, L, temp);
        double rNorm = norm3(r);
        double A[3] = {
            temp[0] - mu * r[0] / rNorm,
            temp[1] - mu * r[1] / rNorm,
            temp[2] - mu * r[2] / rNorm
        };

        // error(i) = k1*||L - L0||^2 + k2*||A - A0||^2
        // Here L0, A0 are assumed to be 3-element vectors
        std::vector<double> Lvec(L, L + 3);
        std::vector<double> Avec(A, A + 3);

        double distL = dist3_squared(Lvec, L0);
        double distA = dist3_squared(Avec, A0);
        error[i] = k1*distL + k2*distA;
    }

    return error;
}


double getError(
    const double& x1,
    const double& x2,
    const double& x3,
    const double& v1,
    const double& v2,
    const double& v3,
    const std::vector<double>& L0,
    const std::vector<double>& A0,
    double k1,
    double k2,
    double mu)
{
    double error;

    double r[3] = { x1, x2, x3 };
    double v[3] = { v1, v2, v3 };

    // Compute L = cross(r, v)
    double L[3];
    cross3(r, v, L);

    // Compute A = cross(v, L) - mu * r / ||r||
    double temp[3];
    cross3(v, L, temp);
    double rNorm = norm3(r);
    double A[3] = {
        temp[0] - mu * r[0] / rNorm,
        temp[1] - mu * r[1] / rNorm,
        temp[2] - mu * r[2] / rNorm
    };

    std::vector<double> Lvec(L, L + 3);
    std::vector<double> Avec(A, A + 3);

    double distL = dist3_squared(Lvec, L0);
    double distA = dist3_squared(Avec, A0);
    error = k1*distL + k2*distA;

    return error;
}


double getError(
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
    double error;

    double r[3] = { x1, x2, x3 };
    double v[3] = { v1, v2, v3 };

    // Compute L = cross(r, v)
    double L[3];
    cross3(r, v, L);

    // Compute A = cross(v, L) - mu * r / ||r||
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

    double distL = dist3_squared(Lvec, L0);
    double distA = dist3_squared(Avec, A0);
    error = k1*distL + k2*distA;

    return error;
}
