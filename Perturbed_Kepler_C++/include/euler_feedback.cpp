#include <cmath>      // for ceil, pow, etc.
#include <vector>
#include <iostream>
#include "basic_operations.h"
#include "getError.h"

using namespace std;

/**
 * @brief dynamics_Euler_feedback(x, mu, alpha, k1, k2, L0, A0)
 *
 * Replicates the MATLAB function:
 *
 *     Fh = [dx, dv] - alpha*[delVx; delVv]
 */
std::array<double,6> dynamics_Euler_feedback(const std::array<double,6>& x,
                                            double mu,
                                            double alpha,
                                            double k1,
                                            double k2,
                                            const std::array<double,3>& L0,
                                            const std::array<double,3>& A0)
{
    // x = [ x(0), x(1), x(2), x(3), x(4), x(5) ] = [r, v]
    // We expect x to have size=6
    // L0 and A0 each size=3

    // Extract position & velocity
    double rVec[3] = { x[0], x[1], x[2] };
    double vVec[3] = { x[3], x[4], x[5] };

    // 1) L = cross(r, v)
    double L[3];
    cross3(rVec, vVec, L);

    // 2) A = cross(v, L) - mu * r / norm(r)
    double tempCross[3];
    cross3(vVec, L, tempCross);  // cross(v, L)
    double rNorm = norm3(rVec);
    
    double A[3];
    A[0] = tempCross[0] - mu * rVec[0] / rNorm;
    A[1] = tempCross[1] - mu * rVec[1] / rNorm;
    A[2] = tempCross[2] - mu * rVec[2] / rNorm;

    // deltaL = L - L0
    double deltaL[3] = { L[0] - L0[0], L[1] - L0[1], L[2] - L0[2] };
    // deltaA = A - A0
    double deltaA[3] = { A[0] - A0[0], A[1] - A0[1], A[2] - A0[2] };

    // ---------- delVx ----------
    // delVx = k1 * cross(v, deltaL)
    //       + k2 * [ cross( v, cross(deltaA, v) )
    //                - mu*deltaA / norm(r)
    //                + mu * r*(r'.deltaA)/norm(r)^3 ]
    double cross_v_deltaL[3];
    cross3(vVec, deltaL, cross_v_deltaL);

    // cross(deltaA, v)
    double cross_deltaA_v[3];
    cross3(deltaA, vVec, cross_deltaA_v);

    // cross( v, cross_deltaA_v )
    double cross_v_crossDAv[3];
    cross3(vVec, cross_deltaA_v, cross_v_crossDAv);

    // - mu * deltaA / norm(r)
    double term_mu_deltaA[3] = {
        -mu * deltaA[0] / rNorm,
        -mu * deltaA[1] / rNorm,
        -mu * deltaA[2] / rNorm
    };

    // mu * r*(r'.deltaA)/norm(r)^3
    double dot_r_deltaA = rVec[0]*deltaA[0] + rVec[1]*deltaA[1] + rVec[2]*deltaA[2];
    double factor = mu * dot_r_deltaA / (rNorm*rNorm*rNorm);
    double term_mu_rdot[3] = {
        factor * rVec[0],
        factor * rVec[1],
        factor * rVec[2]
    };

    double delVx[3];
    for(int i=0; i<3; i++){
        delVx[i] = k1*cross_v_deltaL[i]
                   + k2*( cross_v_crossDAv[i]
                          + term_mu_deltaA[i]
                          + term_mu_rdot[i] );
    }

    // ---------- delVv ----------
    // delVv = k1 * cross(deltaL, r)
    //       + k2 * [ cross( cross(r, v), deltaA ) + cross( r, cross(v, deltaA) ) ]
    double cross_deltaL_r[3];
    cross3(deltaL, rVec, cross_deltaL_r);

    // cross(r, v) = L
    // cross(L, deltaA)
    double cross_L_deltaA[3];
    cross3(L, deltaA, cross_L_deltaA);

    // cross(v, deltaA)
    double cross_v_deltaA[3];
    cross3(vVec, deltaA, cross_v_deltaA);

    // cross(r, cross_v_deltaA)
    double cross_r_crossVD[3];
    cross3(rVec, cross_v_deltaA, cross_r_crossVD);

    double delVv[3];
    for(int i=0; i<3; i++){
        delVv[i] = k1*cross_deltaL_r[i]
                   + k2*( cross_L_deltaA[i]
                          + cross_r_crossVD[i] );
    }

    // ---------- Build gravitational acceleration: -mu*r / norm(r)^3 ----------
    double r3 = rNorm*rNorm*rNorm;
    double acc[3] = {
        -mu * rVec[0]/r3,
        -mu * rVec[1]/r3,
        -mu * rVec[2]/r3
    };

    // Fh = [v, acc] - alpha*[delVx, delVv]
    std::array<double,6> F;
    F[0] = vVec[0] - alpha*delVx[0];
    F[1] = vVec[1] - alpha*delVx[1];
    F[2] = vVec[2] - alpha*delVx[2];
    F[3] = acc[0] - alpha*delVv[0];
    F[4] = acc[1] - alpha*delVv[1];
    F[5] = acc[2] - alpha*delVv[2];

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
               double mu,
               double alpha,
               double k1,
               double k2,
               const std::array<double,3>& L0,
               const std::array<double,3>& A0)
{
    // N = ceil(tf/h) + 1
    long long N = static_cast<long long>(std::ceil(tf/h)) + 1;

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

    // current state/time
    std::array<double,6> current_state = {x1[0], x2[0], x3[0], v1[0], v2[0], v3[0]};
    double current_time = 0.0;

    for(int i=0; i<N; ++i)
    {
        // Compute next state using Euler step
        std::array<double,6> F = dynamics_Euler_feedback(current_state, mu, alpha, k1, k2, L0, A0);

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
                        double mu,
                        double alpha,
                        double k1,
                        double k2,
                        const std::array<double,3>& L0,
                        const std::array<double,3>& A0)
    {
    long long N = static_cast<long long>(std::ceil(tf/h)) + 1;

    double maxV = 0.0;
    double maxdL_sq = 0.0;
    double maxdA_sq = 0.0;
    LAError currentError;
    double x1, x2, x3, v1, v2, v3, t;

    // Set initial values
    x1 = xi[0];  x2 = xi[1];  x3 = xi[2];
    v1 = xi[3];  v2 = xi[4];  v3 = xi[5];
    t  = 0.0;

    std::array<double,6> current_state = {x1, x2, x3, v1, v2, v3};
    double current_time = 0.0;

    for(int i=0; i<N; ++i)
    {
        std::array<double,6> F = dynamics_Euler_feedback(current_state, mu, alpha, k1, k2, L0, A0);
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

        currentError = getError(x1, x2, x3, v1, v2, v3, L0, A0, k1, k2, mu);
        if ( currentError.error >= maxV) maxV = currentError.error;
        if ( currentError.distL_sq >= maxdL_sq) maxdL_sq = currentError.distL_sq;
        if ( currentError.distA_sq >= maxdA_sq) maxdA_sq = currentError.distA_sq;
    }
    return {maxV, sqrt(maxdL_sq), sqrt(maxdA_sq)};
}