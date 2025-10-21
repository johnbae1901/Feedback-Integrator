#include <vector>
#include <tuple>
#include <cmath>
#include <array>
#include <iostream>
#include "getError.h"

using namespace std;

/**
 * @brief Half-step kick1: 
 *        v_temp = v + (h/2) * a(r).
 *        a(r) = -mu * r / ||r||^3
 */
static std::array<double,3> kick1(const std::array<double,3>& x,
                                  const std::array<double,3>& v,
                                  double mu, double h)
{
    std::array<double,3> v_temp;
    // norm(x):
    double r = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    double coeff = -mu / (r*r*r); // acceleration factor
    for(int i=0; i<3; i++){
        v_temp[i] = v[i] + 0.5*h*coeff*x[i];
    }
    return v_temp;
}

/**
 * @brief drift: 
 *        x_next = x + h * v_temp
 */
static std::array<double,3> drift(const std::array<double,3>& x,
                                  const std::array<double,3>& v_temp,
                                  double h)
{
    std::array<double,3> x_next;
    for(int i=0; i<3; i++){
        x_next[i] = x[i] + h*v_temp[i];
    }
    return x_next;
}

/**
 * @brief kick2: 
 *        v_next = v_temp + (h/2) * a(r_next)
 */
static std::array<double,3> kick2(const std::array<double,3>& x_next,
                                  const std::array<double,3>& v_temp,
                                  double mu, double h)
{
    std::array<double,3> v_next;
    double r = std::sqrt(x_next[0]*x_next[0]
                       + x_next[1]*x_next[1]
                       + x_next[2]*x_next[2]);
    double coeff = -mu / (r*r*r);
    for(int i=0; i<3; i++){
        v_next[i] = v_temp[i] + 0.5*h*coeff*x_next[i];
    }
    return v_next;
}

/**
 * @brief stormer_verlet_B integrator for the 6D state:
 *        (x, y, z, vx, vy, vz).
 *
 * @param xi  Initial 6D state: [x, y, z, vx, vy, vz].
 * @param tf  Final time.
 * @param h   Time step.
 * @param mu  Parameter (e.g., gravitational).
 *
 * @return A tuple of 7 vectors:
 *         ( x1, x2, x3, v1, v2, v3, t ).
 */
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
                 double mu)
{
    long long N = static_cast<long long>(std::ceil(tf/h)) + 1;

    // Allocate output arrays
    std::vector<double> x1(N), x2(N), x3(N);
    std::vector<double> v1(N), v2(N), v3(N);
    std::vector<double> t (N);

    // Initialize
    x1[0] = xi[0];  x2[0] = xi[1];  x3[0] = xi[2];
    v1[0] = xi[3];  v2[0] = xi[4];  v3[0] = xi[5];
    t[0]  = 0.0;

    // current state
    std::array<double,3> x = { xi[0], xi[1], xi[2] };
    std::array<double,3> v = { xi[3], xi[4], xi[5] };
    double current_time = 0.0;

    for(int i=0; i<N; ++i)
    {
        // 1) Kick1: half-step velocity
        std::array<double,3> v_temp = kick1(x, v, mu, h);

        // 2) Drift: full-step position
        std::array<double,3> x_next = drift(x, v_temp, h);

        // 3) Kick2: half-step velocity
        std::array<double,3> v_next = kick2(x_next, v_temp, mu, h);

        // Build next state
        x = x_next;
        v = v_next;
        double next_time = current_time + h;

        // Store
        x1[i] = x[0];
        x2[i] = x[1];
        x3[i] = x[2];
        v1[i] = v[0];
        v2[i] = v[1];
        v3[i] = v[2];
        t [i] = current_time;

        // Move on
        current_time = next_time;
    }

    return {x1, x2, x3, v1, v2, v3, t};
}

double stormer_verlet_B_error(const std::vector<double>& xi,
                double tf,
                double h,
                const std::array<double,3>& L0,
                const std::array<double,3>& A0,
                double k1,
                double k2,
                double mu)
{
    long long N = static_cast<long long>(std::ceil(tf/h)) + 1;
    // cout << N << endl;

    double maxError = 0;
    double x1, x2, x3, v1, v2, v3, t;

    // Set initial values
    x1 = xi[0];  x2 = xi[1];  x3 = xi[2];
    v1 = xi[3];  v2 = xi[4];  v3 = xi[5];
    t  = 0.0;

    // current state
    std::array<double,3> x = { xi[0], xi[1], xi[2] };
    std::array<double,3> v = { xi[3], xi[4], xi[5] };
    double current_time = 0.0;

    for(int i=0; i<N; ++i)
    {
        // 1) Kick1: half-step velocity
        std::array<double,3> v_temp = kick1(x, v, mu, h);

        // 2) Drift: full-step position
        std::array<double,3> x_next = drift(x, v_temp, h);

        // 3) Kick2: half-step velocity
        std::array<double,3> v_next = kick2(x_next, v_temp, mu, h);

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

        // double currentError = getError(x1, x2, x3, v1, v2, v3, L0, A0, k1, k2, mu);
        // if ( currentError >= maxError) {
        //     maxError = currentError;
        // }
    }

    return maxError;
}