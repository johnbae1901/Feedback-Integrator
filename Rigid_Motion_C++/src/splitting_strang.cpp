#include "../include/getError.h"
#include "../include/basic_operations.h"
#include <cmath>
#include <algorithm>

static inline void mat3_mul_vec3(const double A[3][3], const double x[3], double y[3]) {
    for (int i = 0; i < 3; ++i)
        y[i] = A[i][0]*x[0] + A[i][1]*x[1] + A[i][2]*x[2];
}

static inline void mat3_mul_mat3(const double A[3][3], const double B[3][3], double C[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j];
}

static inline void mat3_transpose(const double A[3][3], double AT[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            AT[i][j] = A[j][i];
}

static inline void mat3_copy(const double A[3][3], double B[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            B[i][j] = A[i][j];
}

std::tuple<
    double, 
    double, 
    double, 
    double> 
            splitting_strang(
            const double R0[3][3],
            const double Omega0[3],
            const double I[3][3],
            const double E0,
            const double Pi0[3],
            const double k0,
            const double k1,
            const double k2,
            double tf,
            double h
                        ) {
    // Steps follow user's Euler style: N = ceil(tf/h) + 1
    const int N = static_cast<int>(std::ceil(tf / h)) + 1;

    // Assume diagonal inertia in body frame (as in the MATLAB splitting code)
    const double I1 = I[0][0], I2 = I[1][1], I3 = I[2][2];

    // State variables
    double R[3][3]; mat3_copy(R0, R);
    double y[3] = { I1*Omega0[0], I2*Omega0[1], I3*Omega0[2] }; // y = I*Omega

    auto rotX = [](double th, double S[3][3]) {
        const double c = std::cos(th), s = std::sin(th);
        S[0][0]=1.0; S[0][1]=0.0; S[0][2]=0.0;
        S[1][0]=0.0; S[1][1]= c ; S[1][2]= s ;
        S[2][0]=0.0; S[2][1]=-s ; S[2][2]= c ;
    };
    auto rotY = [](double th, double S[3][3]) {
        const double c = std::cos(th), s = std::sin(th);
        S[0][0]= c ; S[0][1]=0.0; S[0][2]=-s ;
        S[1][0]=0.0; S[1][1]=1.0; S[1][2]=0.0;
        S[2][0]= s ; S[2][1]=0.0; S[2][2]= c ;
    };
    auto rotZ = [](double th, double S[3][3]) {
        const double c = std::cos(th), s = std::sin(th);
        S[0][0]= c ; S[0][1]= s ; S[0][2]=0.0;
        S[1][0]=-s ; S[1][1]= c ; S[1][2]=0.0;
        S[2][0]=0.0; S[2][1]=0.0; S[2][2]=1.0;
    };

    double S[3][3], ST[3][3], Rtmp[3][3];
    double maxV = 0.0;
    double maxdE = 0.0;
    double maxdPi_sq = 0.0;
    double maxdDet_sq = 0.0;
    ErrorReport currentError;

    // current_state in 3x4 layout for getError(...)
    auto eval_and_update_max = [&](void){
        // Build current_state = [R | Omega]
        const double Omega[3] = { y[0]/I1, y[1]/I2, y[2]/I3 };
        double current_state[3][4];
        for (int i = 0; i < 3; ++i) {
            current_state[i][0] = R[i][0];
            current_state[i][1] = R[i][1];
            current_state[i][2] = R[i][2];
            current_state[i][3] = Omega[i];
        }
        currentError = getError(current_state, I, E0, Pi0, k0, k1, k2);
        if (currentError.weighted_sum >= maxV)
            maxV = currentError.weighted_sum;
        if (currentError.abs_deltaE >= maxdE)
            maxdE = currentError.abs_deltaE;
        if (currentError.deltaPi_norm_sq >= maxdPi_sq)
            maxdPi_sq = currentError.deltaPi_norm_sq;
        if (currentError.frob_RT_R_minus_I_sq >= maxdDet_sq)
            maxdDet_sq = currentError.frob_RT_R_minus_I_sq;
    };

    for (int k = 0; k < N; ++k) {
        // One Strang step of size h
        // R3 (h/2)
        {
            const double theta = (y[2]/I3)*(h*0.5);
            rotZ(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew); y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }
        // R2 (h/2)
        {
            const double theta = (y[1]/I2)*(h*0.5);
            rotY(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew); y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }
        // R1 (h)
        {
            const double theta = (y[0]/I1)*h;
            rotX(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew); y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }
        // R2 (h/2)
        {
            const double theta = (y[1]/I2)*(h*0.5);
            rotY(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew); y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }
        // R3 (h/2)
        {
            const double theta = (y[2]/I3)*(h*0.5);
            rotZ(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew); y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }

        // Evaluate error after completing the full step
        eval_and_update_max();
    }
    return {maxV, maxdE, sqrt(maxdPi_sq), sqrt(maxdDet_sq)};
}
