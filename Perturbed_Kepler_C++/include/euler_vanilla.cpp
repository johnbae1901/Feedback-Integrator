#include <cmath>      // for std::ceil, std::sqrt
#include <vector>
#include <tuple>
#include <array>
#include <iostream>
#include "getError.h"

using namespace std;

/**
 * @brief Compute the 6x1 dynamics for the simple orbital model:
 *          Fh = [ v1, v2, v3, -mu x1/r^3, -mu x2/r^3, -mu x3/r^3 ],
 *        where r = sqrt(x1^2 + x2^2 + x3^2).
 */
std::array<double,6> dynamics_Euler(const std::array<double,6>& x, double mu)
{
    // x = [x1, x2, x3, v1, v2, v3]
    // Position components
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    // Velocity components
    double v1 = x[3];
    double v2 = x[4];
    double v3 = x[5];

    // Compute norm(r)
    double r = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    double r3 = r * r * r;

    std::array<double,6> F;
    F[0] = v1;
    F[1] = v2;
    F[2] = v3;
    F[3] = -mu * x1 / r3;
    F[4] = -mu * x2 / r3;
    F[5] = -mu * x3 / r3;
    return F;
}

/**
 * @brief Vanilla Euler integrator:
 *        [x1, x2, x3, v1, v2, v3, t] = euler_vanilla(xi, tf, h, mu)
 *
 * @param xi  Initial 6D state: (x1, x2, x3, v1, v2, v3)
 * @param tf  Final time
 * @param h   Time step
 * @param mu  Parameter (e.g. gravitational constant)
 *
 * @return A tuple of 7 std::vectors: (x1, x2, x3, v1, v2, v3, t).
 */
std::tuple<
    std::vector<double>, // x1
    std::vector<double>, // x2
    std::vector<double>, // x3
    std::vector<double>, // v1
    std::vector<double>, // v2
    std::vector<double>, // v3
    std::vector<double>  // t
>
euler_vanilla(const std::vector<double>& xi, double tf, double h, double mu)
{
    // Number of points we will store
    long long N = static_cast<long long>(std::ceil(tf / h)) + 1;

    const int stride = 1000;
    const long long approx = 1 + (N - 1) / stride + 1;
    vector<double> x1; x1.reserve(approx);
    vector<double> x2; x2.reserve(approx);
    vector<double> x3; x3.reserve(approx);
    vector<double> v1; v1.reserve(approx);
    vector<double> v2; v2.reserve(approx);
    vector<double> v3; v3.reserve(approx);
    vector<double> t;  t.reserve(approx);

    // Set initial values
    x1.push_back(xi[0]);  x2.push_back(xi[1]);  x3.push_back(xi[2]);
    v1.push_back(xi[3]);  v2.push_back(xi[4]);  v3.push_back(xi[5]);
    t.push_back(0.0);

    // Current state and time
    std::array<double,6> current_state = {x1[0], x2[0], x3[0], v1[0], v2[0], v3[0]};
    double current_time = 0.0;

    // Euler integration loop
    for(int i = 0; i < N - 1; ++i)
    {
        // Compute the derivative at the current state
        std::array<double,6> F = dynamics_Euler(current_state, mu);
        // Advance one step: next_state = current_state + h * F
        std::array<double,6> next_state;
        for(int j = 0; j < 6; ++j) {
            next_state[j] = current_state[j] + h * F[j];
        }
        double next_time = current_time + h;

        // Update
        current_state = next_state;
        current_time  = next_time;

        const bool periodic_save = ((i + 1) % stride == 0);
        const bool is_last       = (i == N - 2);
        if (periodic_save || (is_last && !periodic_save)) {
            x1.push_back(current_state[0]); x2.push_back(current_state[1]); x3.push_back(current_state[2]);
            v1.push_back(current_state[3]); v2.push_back(current_state[4]); v3.push_back(current_state[5]);
            t .push_back(current_time);
        }
    }

    return { x1, x2, x3, v1, v2, v3, t };
}


std::tuple<double, double, double>
        euler_vanilla_error(const std::vector<double>& xi, double tf, double h, 
                            const std::array<double,3>& L0,
                            const std::array<double,3>& A0,
                            double k1,
                            double k2,
                            double mu)
{
    // Number of points we will store
    long long N = static_cast<long long>(std::ceil(tf / h)) + 1;
    
    double maxV = 0.0;
    double maxdL_sq = 0.0;
    double maxdA_sq = 0.0;
    LAError currentError;
    double x1, x2, x3, v1, v2, v3, t;

    x1 = xi[0];  x2 = xi[1];  x3 = xi[2];
    v1 = xi[3];  v2 = xi[4];  v3 = xi[5];
    t  = 0.0;

    std::array<double,6> current_state = {x1, x2, x3, v1, v2, v3};
    double current_time = 0.0;

    for(int i = 0; i < N - 1; ++i)
    {
        std::array<double,6> F = dynamics_Euler(current_state, mu);
        std::array<double,6> next_state;
        for(int j = 0; j < 6; ++j) {
            next_state[j] = current_state[j] + h * F[j];
        }
        double next_time = current_time + h;

        current_state = next_state;
        current_time  = next_time;

        x1 = current_state[0];
        x2 = current_state[1];
        x3 = current_state[2];
        v1 = current_state[3];
        v2 = current_state[4];
        v3 = current_state[5];
        t  = current_time;

        currentError = getError(x1, x2, x3, v1, v2, v3, L0, A0, k1, k2, mu);
        if ( currentError.error >= maxV) maxV = currentError.error;
        if ( currentError.distL_sq >= maxdL_sq) maxdL_sq = currentError.distL_sq;
        if ( currentError.distA_sq >= maxdA_sq) maxdA_sq = currentError.distA_sq;
    }
    return {maxV, sqrt(maxdL_sq), sqrt(maxdA_sq)};
}