#include "euler_feedback_adaptive.h"
#include "basic_operations.h"  // for cross3(...) and norm3(...)
#include <cmath>               // for std::ceil, pow
#include <vector>
#include <tuple>
#include <iostream>            // (Optional) for debug prints
#include "getError.h"
#include "estimateL.h"
#include <Eigen/Core>
#include <Eigen/LU>        // only for the (rare) spectral‑norm path
#include <Eigen/Core>
#include <Eigen/Dense>
#include <array>
#include <cassert>

using namespace std;

//------------------------------------
// dynamics_Euler_feedback
//------------------------------------
std::array<double,6> dynamics_Euler_feedback(const std::array<double,6>& x,
                                            double mu,
                                            double c,
                                            double h,
                                            double k1,
                                            double k2,
                                            const std::array<double,3>& L0,
                                            const std::array<double,3>& A0,
                                            const double lambda)
{
    // Expect x.size()=6: x=[x1,x2,x3, v1,v2,v3]
    double r[3] = { x[0], x[1], x[2] };
    double v[3] = { x[3], x[4], x[5] };

    // 1) L = cross(r,v)
    double L[3];
    cross3(r, v, L);

    // 2) A = cross(v,L) - mu*r / ||r||
    double tempCross[3];
    cross3(v, L, tempCross);
    double rNorm = norm3(r);

    double A[3];
    A[0] = tempCross[0] - mu*r[0]/rNorm;
    A[1] = tempCross[1] - mu*r[1]/rNorm;
    A[2] = tempCross[2] - mu*r[2]/rNorm;

    // deltaL = L - L0
    double deltaL[3] = { L[0] - L0[0], L[1] - L0[1], L[2] - L0[2] };
    // deltaA = A - A0
    double deltaA[3] = { A[0] - A0[0], A[1] - A0[1], A[2] - A0[2] };

    // 3) Compute lipschitz constant
    double lipConstant = lambda * estimateL(x, mu, k1, k2);
    // double lipConstant = 1.0 / h;

    // alpha = 1 / (h * lipConstant)
    double alpha = 0.0;
    if(lipConstant > 1.0e-14){
        alpha = 1.0 / (h * lipConstant);
    }

    // delVx
    // = k1 * cross(v, deltaL)
    //   + k2 * [ cross(v, cross(deltaA, v)) - mu*deltaA / ||r|| + mu*r*( (r'.deltaA) / ||r||^3 ) ]
    double cross_v_dL[3]; // cross(v, deltaL)
    cross3(v, deltaL, cross_v_dL);

    double cross_dA_v[3]; // cross(deltaA, v)
    cross3(deltaA, v, cross_dA_v);

    double cross_v_crossDA[3]; // cross(v, cross_dA_v)
    cross3(v, cross_dA_v, cross_v_crossDA);

    // - mu*deltaA/||r||
    double term_mu_deltaA[3] = {
        -mu*deltaA[0]/rNorm,
        -mu*deltaA[1]/rNorm,
        -mu*deltaA[2]/rNorm
    };
    //  mu*r*(r'.deltaA)/||r||^3
    double dot_r_dA = r[0]*deltaA[0] + r[1]*deltaA[1] + r[2]*deltaA[2];
    double factor = mu * dot_r_dA / (rNorm*rNorm*rNorm);
    double term_mu_r[3] = {
        factor*r[0],
        factor*r[1],
        factor*r[2]
    };

    double delVx[3];
    for(int i=0; i<3; i++){
        delVx[i] = k1*cross_v_dL[i]
                   + k2*( cross_v_crossDA[i] + term_mu_deltaA[i] + term_mu_r[i] );
    }

    // delVv
    // = k1 * cross(deltaL, r)
    //   + k2 * [ cross( cross(r,v), deltaA ) + cross(r, cross(v, deltaA)) ]
    double cross_dL_r[3];
    cross3(deltaL, r, cross_dL_r);

    // cross(r,v) => L above
    double cross_L_dA[3];
    cross3(L, deltaA, cross_L_dA);

    double cross_v_dA[3];
    cross3(v, deltaA, cross_v_dA);

    double cross_r_crossv_dA[3];
    cross3(r, cross_v_dA, cross_r_crossv_dA);

    double delVv[3];
    for(int i=0; i<3; i++){
        delVv[i] = k1*cross_dL_r[i]
                   + k2*( cross_L_dA[i] + cross_r_crossv_dA[i] );
    }

    // gravitational acceleration: -mu*r / ||r||^3
    double r3 = rNorm*rNorm*rNorm;
    double acc[3] = {
        -mu*r[0]/r3,
        -mu*r[1]/r3,
        -mu*r[2]/r3
    };

    // Fh = [v, acc] - alpha*[delVx, delVv]
    // We'll store in std::vector<double> of length=6
    std::array<double,6> F;
    // Positions change = v
    F[0] = v[0] - alpha*delVx[0];
    F[1] = v[1] - alpha*delVx[1];
    F[2] = v[2] - alpha*delVx[2];

    // Velocity change = acc
    F[3] = acc[0] - alpha*delVv[0];
    F[4] = acc[1] - alpha*delVv[1];
    F[5] = acc[2] - alpha*delVv[2];
    return F;
}

//------------------------------------
// euler_feedback_adaptive
//------------------------------------
std::tuple<
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>
>
euler_feedback_adaptive(const std::vector<double>& xi,
                        double tf,
                        double h,
                        double mu,
                        double c,
                        double k1,
                        double k2,
                        const std::array<double,3>& L0,
                        const std::array<double,3>& A0,
                        const double lambda)
{
    // N = ceil(tf/h) + 1
    long long N = static_cast<long long>(std::ceil(tf/h)) + 1;
    cout << N << endl;

    std::vector<double> x1(N), x2(N), x3(N);
    std::vector<double> v1(N), v2(N), v3(N);
    std::vector<double> t(N);

    // Initialize
    x1[0] = xi[0];
    x2[0] = xi[1];
    x3[0] = xi[2];
    v1[0] = xi[3];
    v2[0] = xi[4];
    v3[0] = xi[5];
    t [0] = 0.0;

    // Current state/time
    std::array<double,6> current_state = {x1[0], x2[0], x3[0], v1[0], v2[0], v3[0]};
    double current_time = 0.0;

    // Euler loop
    for(int i = 0; i < N; ++i)
    {
        // Evaluate dynamics at current_state
        std::array<double,6> F = dynamics_Euler_feedback(
            current_state, mu, c, h, k1, k2, L0, A0, lambda
        );

        // next_state = current_state + h*F
        std::array<double,6> next_state;
        for(int j=0; j<6; j++){
            next_state[j] = current_state[j] + h*F[j];
        }
        double next_time = current_time + h;

        // Store in x1(i), x2(i), ...
        x1[i] = current_state[0];
        x2[i] = current_state[1];
        x3[i] = current_state[2];
        v1[i] = current_state[3];
        v2[i] = current_state[4];
        v3[i] = current_state[5];
        t [i] = current_time;

        // Update
        current_state = next_state;
        current_time  = next_time;
    }

    return {x1, x2, x3, v1, v2, v3, t};
}


double euler_feedback_adaptive_error(const std::vector<double>& xi,
                        double tf,
                        double h,
                        double mu,
                        double c,
                        double k1,
                        double k2,
                        const std::array<double,3>& L0,
                        const std::array<double,3>& A0,
                        const double lambda)
{
    // N = ceil(tf/h) + 1
    long long N = static_cast<long long>(std::ceil(tf/h)) + 1;

    double maxError = 0;
    double x1, x2, x3, v1, v2, v3, t;

    // Set initial values
    x1 = xi[0];  x2 = xi[1];  x3 = xi[2];
    v1 = xi[3];  v2 = xi[4];  v3 = xi[5];
    t  = 0.0;

    // Current state/time
    std::array<double,6> current_state = {x1, x2, x3, v1, v2, v3};
    double current_time = 0.0;

    // Euler loop
    for(int i = 0; i < N; ++i)
    {
        // Evaluate dynamics at current_state
        std::array<double,6> F = dynamics_Euler_feedback(
            current_state, mu, c, h, k1, k2, L0, A0, lambda
        );

        // next_state = current_state + h*F
        std::array<double,6> next_state;
        for(int j=0; j<6; j++){
            next_state[j] = current_state[j] + h*F[j];
        }
        double next_time = current_time + h;

        // Update
        current_state = next_state;
        current_time  = next_time;

        // Store results for step i+1
        x1 = current_state[0];
        x2 = current_state[1];
        x3 = current_state[2];
        v1 = current_state[3];
        v2 = current_state[4];
        v3 = current_state[5];
        t  = current_time;

        // double currentError = getError(x1, x2, x3, v1, v2, v3, L0, A0, k1, k2, mu);
        // if ( currentError >= maxError) {
        //     maxError = currentError;
        // }
    }

    return maxError;
}
