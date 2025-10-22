#include "splitting_strang_full.h"
#include <cmath>

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

void splitting_strang_full(
    const double R0[3][3],
    const double Omega0[3],
    const double I[3][3],
    double tf,
    double h,
    double*& R_hist,
    double*& Omega_hist,
    double*& t_hist,
    int& N_out
) {
    const long long N_total = static_cast<long long>(std::ceil(tf / h)) + 1;

    const int stride = 1000;
    const long long approx = 1 + ((N_total - 1) + stride - 1) / stride;

    const double I1 = I[0][0], I2 = I[1][1], I3 = I[2][2];

    R_hist     = new double[9 * approx];
    Omega_hist = new double[3 * approx];
    t_hist     = new double[approx];

    double R[3][3]; mat3_copy(R0, R);
    double y[3] = { I1*Omega0[0], I2*Omega0[1], I3*Omega0[2] }; // y = I * Omega

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

    auto save_R = [&](long long k){
        const long long base = 9*k;
        for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
                R_hist[base + 3*i + j] = R[i][j];
    };

    double S[3][3], ST[3][3], Rtmp[3][3];

    long long k = 0;
    double t = 0.0;

    Omega_hist[3*k+0] = y[0]/I1;
    Omega_hist[3*k+1] = y[1]/I2;
    Omega_hist[3*k+2] = y[2]/I3;
    save_R(k);
    ++k;

    for (long long i = 0; i < N_total - 1; ++i) {
        // Strang step
        {
            const double theta = (y[2]/I3)*(h*0.5);
            rotZ(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew);
            y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }
        {
            const double theta = (y[1]/I2)*(h*0.5);
            rotY(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew);
            y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }
        {
            const double theta = (y[0]/I1)*h;
            rotX(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew);
            y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }
        {
            const double theta = (y[1]/I2)*(h*0.5);
            rotY(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew);
            y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }
        {
            const double theta = (y[2]/I3)*(h*0.5);
            rotZ(theta, S);
            double ynew[3]; mat3_mul_vec3(S, y, ynew);
            y[0]=ynew[0]; y[1]=ynew[1]; y[2]=ynew[2];
            mat3_transpose(S, ST);
            mat3_mul_mat3(R, ST, Rtmp); mat3_copy(Rtmp, R);
        }

        t += h;

        const bool periodic_save = ((i + 1) % stride == 0);
        const bool is_last       = (i == N_total - 2);
        if (periodic_save || (is_last && !periodic_save)) {
            t_hist[k] = t;
            Omega_hist[3*k+0] = y[0]/I1;
            Omega_hist[3*k+1] = y[1]/I2;
            Omega_hist[3*k+2] = y[2]/I3;
            save_R(k);
            ++k;
        }
    }
    N_out = static_cast<int>(k);
}
