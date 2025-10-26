#include <vector>
#include <tuple>
#include <cmath>
#include <array>
#include <iostream>
#include "getError.h"

using namespace std;

// half-kick
static std::array<double,3> kick1(const std::array<double,3>& x,
                                  const std::array<double,3>& v,
                                  double h,
                                  const PKParams& P)
{
    std::array<double,3> v_half = v;
    double r = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    double coeff = -(P.mu/(r*r*r) + 3.0*P.delta/(r*r*r*r*r));
    v_half[0] += 0.5*h*coeff*x[0];
    v_half[1] += 0.5*h*coeff*x[1];
    v_half[2] += 0.5*h*coeff*x[2];
    return v_half;
}

static std::array<double,3> drift(const std::array<double,3>& x,
                                  const std::array<double,3>& v_half,
                                  double h)
{
    return { x[0] + h*v_half[0],
             x[1] + h*v_half[1],
             x[2] + h*v_half[2] };
}

static std::array<double,3> kick2(const std::array<double,3>& x_next,
                                  const std::array<double,3>& v_half,
                                  double h,
                                  const PKParams& P)
{
    std::array<double,3> v_next = v_half;
    double r = std::sqrt(x_next[0]*x_next[0]
                       + x_next[1]*x_next[1]
                       + x_next[2]*x_next[2]);
    double coeff = -(P.mu/(r*r*r) + 3.0*P.delta/(r*r*r*r*r));
    v_next[0] += 0.5*h*coeff*x_next[0];
    v_next[1] += 0.5*h*coeff*x_next[1];
    v_next[2] += 0.5*h*coeff*x_next[2];
    return v_next;
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
stormer_verlet_B(const std::vector<double>& xi,
                 double tf,
                 double h,
                 const PKParams& P)
{
    long long N = static_cast<long long>(std::ceil(tf/h)) + 1;

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

    // current state
    std::array<double,3> x = { xi[0], xi[1], xi[2] };
    std::array<double,3> v = { xi[3], xi[4], xi[5] };
    double current_time = 0.0;

    for(int i=0; i<N; ++i)
    {
        // 1) Kick1: half-step velocity
        std::array<double,3> v_half = kick1(x, v, h, P);

        // 2) Drift: full-step position
        std::array<double,3> x_next = drift(x, v_half, h);

        // 3) Kick2: half-step velocity
        std::array<double,3> v_next = kick2(x_next, v_half, h, P);

        // Build next state
        x = x_next;
        v = v_next;
        double next_time = current_time + h;

        current_time = next_time;

        const bool periodic_save = ((i + 1) % stride == 0);
        const bool is_last       = (i == N - 2);
        if (periodic_save || (is_last && !periodic_save)) {
            x1.push_back(x[0]); x2.push_back(x[1]); x3.push_back(x[2]);
            v1.push_back(v[0]); v2.push_back(v[1]); v3.push_back(v[2]);
            t .push_back(current_time);
        }
    }

    return {x1, x2, x3, v1, v2, v3, t};
}

std::tuple<double, double, double> 
    stormer_verlet_B_error(const std::vector<double>& xi,
                            double tf,
                            double h,
                            double k1,  // -> kL
                            double k2,  // -> kE
                            const PKParams& P)
{
    long long N = static_cast<long long>(std::ceil(tf/h)) + 1;

    double maxV = 0.0;
    double maxdE_sq = 0.0;
    double maxdL_sq = 0.0;
    LEError currentError;
    double x1, x2, x3, v1, v2, v3, t;

    // Set initial values
    x1 = xi[0];  x2 = xi[1];  x3 = xi[2];
    v1 = xi[3];  v2 = xi[4];  v3 = xi[5];
    t  = 0.0;

    // current state
    std::array<double,3> x = { xi[0], xi[1], xi[2] };
    std::array<double,3> v = { xi[3], xi[4], xi[5] };
    double current_time = 0.0;

    // ---- compute L0, E0 internally ----
    const std::array<double,3> L0 = {
        x2*v3 - x3*v2,
        x3*v1 - x1*v3,
        x1*v2 - x2*v1
    };
    const double r0  = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    const double v02 = (v1*v1 + v2*v2 + v3*v3);
    const double E0  = 0.5*v02 - P.mu/r0 - P.delta/(r0*r0*r0);

    for(int i=0; i<N; ++i)
    {
        // 1) Kick1: half-step velocity
        std::array<double,3> v_half = kick1(x, v, h, P);

        // 2) Drift: full-step position
        std::array<double,3> x_next = drift(x, v_half, h);

        // 3) Kick2: half-step velocity
        std::array<double,3> v_next = kick2(x_next, v_half, h, P);

        // Build next state
        x = x_next;
        v = v_next;
        double next_time = current_time + h;

        // Move on
        current_time = next_time;

        // Store
        x1 = x[0];
        x2 = x[1];
        x3 = x[2];
        v1 = v[0];
        v2 = v[1];
        v3 = v[2];
        t  = current_time;

        currentError = getError_PK(x1,x2,x3, v1,v2,v3, L0, E0, k1, k2, P);
        if ( currentError.error >= maxV) maxV = currentError.error;
        if ( currentError.distE_sq >= maxdE_sq) maxdE_sq = currentError.distE_sq;
        if ( currentError.distL_sq >= maxdL_sq) maxdL_sq = currentError.distL_sq;
    }
    return {maxV, sqrt(maxdE_sq), sqrt(maxdL_sq)};
}