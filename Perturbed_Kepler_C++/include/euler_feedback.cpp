#include <cmath>      // for ceil, pow, etc.
#include <vector>
#include <iostream>
#include "basic_operations.h"
#include "getError.h"

using namespace std;

std::array<double,6> dynamics_Euler_feedback(const std::array<double,6>& x,
                                             const PKParams& P,
                                             double kL_eff,
                                             double kE_eff,
                                             const std::array<double,3>& L0,
                                             double E0)
{
    const double rVec[3] = { x[0], x[1], x[2] };
    const double vVec[3] = { x[3], x[4], x[5] };
    const double r  = std::max(1e-15, norm3(rVec));
    const double upr_over_r = Ur_scalar(r, P) / r;

    // deviations
    double L[3]; cross3(rVec, vVec, L);
    const double dL[3] = { L[0]-L0[0], L[1]-L0[1], L[2]-L0[2] };
    const double E = 0.5*(vVec[0]*vVec[0]+vVec[1]*vVec[1]+vVec[2]*vVec[2]) + U_scalar(r,P);
    const double dE = E - E0;

    // feedback terms
    double vXdL[3]; cross3(vVec, dL, vXdL);
    double dLxX[3]; cross3(dL, rVec, dLxX);

    // base accel
    double ax[3]; accel(rVec, P, ax);

    std::array<double,6> F;
    // xdot
    F[0] = vVec[0] - kE_eff*dE*upr_over_r*rVec[0] - kL_eff*vXdL[0];
    F[1] = vVec[1] - kE_eff*dE*upr_over_r*rVec[1] - kL_eff*vXdL[1];
    F[2] = vVec[2] - kE_eff*dE*upr_over_r*rVec[2] - kL_eff*vXdL[2];
    // vdot
    F[3] = ax[0]   - kE_eff*dE*vVec[0]           - kL_eff*dLxX[0];
    F[4] = ax[1]   - kE_eff*dE*vVec[1]           - kL_eff*dLxX[1];
    F[5] = ax[2]   - kE_eff*dE*vVec[2]           - kL_eff*dLxX[2];

    return F;
}

/**
 * @brief euler_feedback(xi, tf, h, mu, alpha, k1, k2, L0, A0)
 *        returns {x1, x2, x3, v1, v2, v3, t} in 7 parallel std::vectors.
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
euler_feedback(const std::vector<double>& xi,
               double tf,
               double h,
               const PKParams& P,
               double alpha,
               double kL,
               double kE)
{
    // N = ceil(tf/h) + 1
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

    // current state/time
    std::array<double,6> current_state = {x1[0], x2[0], x3[0], v1[0], v2[0], v3[0]};
    double current_time = 0.0;

    const std::array<double,3> L0 = {
        x2[0]*v3[0] - x3[0]*v2[0],
        x3[0]*v1[0] - x1[0]*v3[0],
        x1[0]*v2[0] - x2[0]*v1[0]
    };
    const double r0  = std::max(1e-15, std::sqrt(x1[0]*x1[0]+x2[0]*x2[0]+x3[0]*x3[0]));
    const double v02 = (v1[0]*v1[0]+v2[0]*v2[0]+v3[0]*v3[0]);
    const double E0  = 0.5*v02 + U_scalar(r0,P);

    const double kL_eff = alpha * kL;
    const double kE_eff = alpha * kE;

    for(int i=0; i<N; ++i)
    {
        // Compute next state using Euler step
        std::array<double,6> F = dynamics_Euler_feedback(current_state, P, kL_eff, kE_eff, L0, E0);

        // next_state = current_state + h*F
        std::array<double,6> next_state;
        for(int j=0; j<6; j++){
            next_state[j] = current_state[j] + h*F[j];
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

    return {x1, x2, x3, v1, v2, v3, t};
}

std::tuple<double, double, double> 
    euler_feedback_error(const std::vector<double>& xi,
                        double tf,
                        double h,
                        const PKParams& P,
                        double alpha,
                        double kL,
                        double kE)
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

    std::array<double,6> current_state = {x1, x2, x3, v1, v2, v3};
    double current_time = 0.0;

    const std::array<double,3> L0 = {
        xi[1]*xi[5] - xi[2]*xi[4],
        xi[2]*xi[3] - xi[0]*xi[5],
        xi[0]*xi[4] - xi[1]*xi[3]
    };
    const double r0  = std::max(1e-15, std::sqrt(xi[0]*xi[0]+xi[1]*xi[1]+xi[2]*xi[2]));
    const double v02 = (xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5]);
    const double E0  = 0.5*v02 + U_scalar(r0,P);

    const double kL_eff = alpha * kL;
    const double kE_eff = alpha * kE;

    for(int i=0; i<N; ++i)
    {
        std::array<double,6> F = dynamics_Euler_feedback(current_state, P, kL_eff, kE_eff, L0, E0);
        std::array<double,6> next_state;
        for(int j=0; j<6; j++){
            next_state[j] = current_state[j] + h*F[j];
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