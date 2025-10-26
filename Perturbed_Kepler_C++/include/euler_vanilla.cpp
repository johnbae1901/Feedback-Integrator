#include <cmath>      // for std::ceil, std::sqrt
#include <vector>
#include <tuple>
#include <array>
#include <iostream>
#include "getError.h"
#include "euler_vanilla.h"

using namespace std;

// double getError_PK(
//     const double& x1, const double& x2, const double& x3,
//     const double& v1, const double& v2, const double& v3,
//     const std::array<double,3>& L0,
//     double E0,
//     double kL,
//     double kE,
//     const PKParams& P);

std::array<double,6> dynamics_Euler(const std::array<double,6>& x, const PKParams& P)
{
    const double x1 = x[0], x2 = x[1], x3 = x[2];
    const double v1 = x[3], v2 = x[4], v3 = x[5];
    const double r  = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    const double r3 = r*r*r, r5 = r3*r*r;
    const double c  = -(P.mu/r3 + 3.0*P.delta/r5);

    std::array<double,6> F;
    F[0] = v1; F[1] = v2; F[2] = v3;
    F[3] = c*x1; F[4] = c*x2; F[5] = c*x3;
    return F;
}

std::tuple<
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>
>
euler_vanilla(const std::vector<double>& xi, double tf, double h, const PKParams& P)
{
    // Number of points we will store
    long long N = static_cast<long long>(std::ceil(tf / h)) + 1;

    const int stride = 10;
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

    for(int i = 0; i < N - 1; ++i)
    {
        std::array<double,6> F = dynamics_Euler(current_state, P);
        std::array<double,6> next_state;
        for(int j = 0; j < 6; ++j) {
            next_state[j] = current_state[j] + h * F[j];
        }
        
        double next_time = current_time + h;
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
        euler_vanilla_error(const std::vector<double>& xi,
                           double tf,
                           double h,
                           double kL,
                           double kE,
                           const PKParams& P)
{
    // Number of points we will store
    long long N = static_cast<long long>(std::ceil(tf / h)) + 1;
    
    double maxV = 0.0;
    double maxdE_sq = 0.0;
    double maxdL_sq = 0.0;
    LEError currentError;
    double x1, x2, x3, v1, v2, v3, t;

    x1 = xi[0];  x2 = xi[1];  x3 = xi[2];
    v1 = xi[3];  v2 = xi[4];  v3 = xi[5];
    t  = 0.0;

    std::array<double,6> current_state = {x1, x2, x3, v1, v2, v3};
    double current_time = 0.0;

    const std::array<double,3> L0 = {
        x2*v3 - x3*v2,
        x3*v1 - x1*v3,
        x1*v2 - x2*v1
    };
    const double r0  = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    const double v02 = (v1*v1 + v2*v2 + v3*v3);
    const double E0  = 0.5*v02 - P.mu/r0 - P.delta/(r0*r0*r0);

    for(int i = 0; i < N - 1; ++i)
    {
        std::array<double,6> F = dynamics_Euler(current_state, P);
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

        currentError = getError_PK(x1,x2,x3, v1,v2,v3, L0, E0, kL, kE, P);
        if ( currentError.error >= maxV) maxV = currentError.error;
        if ( currentError.distE_sq >= maxdE_sq) maxdE_sq = currentError.distE_sq;
        if ( currentError.distL_sq >= maxdL_sq) maxdL_sq = currentError.distL_sq;
    }
    return {maxV, sqrt(maxdE_sq), sqrt(maxdL_sq)};
}