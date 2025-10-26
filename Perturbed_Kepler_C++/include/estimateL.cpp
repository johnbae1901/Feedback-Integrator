#include "basic_operations.h"     // crossProdMtx, norm3
#include <array>
#include <cmath>
#include <algorithm>
#include "estimateL.h"

// H = kE*(J_E^T J_E) + kL*(J_L^T J_L)
// J_E = [ gEx^T  gEv^T ] (1x6),  gEx=(mu/r^3+3δ/r^5) x,  gEv=v
// J_L = [ -[v]_x   [x]_x ] (3x6)
// Return ||H||_F = sqrt( sum_ij H_ij^2 )
double estimateL(const std::array<double,6>& x,
                 const PKParams& P,
                 double kL,
                 double kE)
{
    // unpack
    const double r[3] = { x[0], x[1], x[2] };
    const double v[3] = { x[3], x[4], x[5] };
    const double R    = std::max(1e-15, norm3(r)); // avoid 0-div

    // ----- build J_L^T J_L via M = [ -[v]_x  0 ; 0  [x]_x ] -----
    double xMtx[9]; crossProdMtx(r, xMtx); // [x]_x
    double vMtx[9]; crossProdMtx(v, vMtx); // [v]_x

    double M6x6[36] = {0.0};
    // top-left 3x3: -[v]_x
    for (int row=0; row<3; ++row){
        const int base = row*6;
        const int src  = row*3;
        M6x6[base+0] = -vMtx[src+0];
        M6x6[base+1] = -vMtx[src+1];
        M6x6[base+2] = -vMtx[src+2];
    }
    // bottom-right 3x3: [x]_x
    for (int row=0; row<3; ++row){
        const int base = (row+3)*6 + 3;
        const int src  = row*3;
        M6x6[base+0] =  xMtx[src+0];
        M6x6[base+1] =  xMtx[src+1];
        M6x6[base+2] =  xMtx[src+2];
    }

    // B = J_L^T J_L = M^T M
    double B[36];
    for (int i=0;i<6;++i){
        for (int j=0;j<6;++j){
            double s=0.0;
            for (int k=0;k<6;++k) s += M6x6[k*6+i]*M6x6[k*6+j];
            B[i*6+j] = s;
        }
    }

    // ----- build A = J_E^T J_E (rank-1 outer product) -----
    // gEx = (mu/r^3 + 3δ/r^5) * x
    const double coef = P.mu/(R*R*R) + 3.0*P.delta/(R*R*R*R*R);
    const double gEx[3] = { coef*r[0], coef*r[1], coef*r[2] };
    const double gEv[3] = { v[0], v[1], v[2] };
    double JE[6] = { gEx[0], gEx[1], gEx[2], gEv[0], gEv[1], gEv[2] };

    double A[36];
    for (int i=0;i<6;++i)
        for (int j=0;j<6;++j)
            A[i*6+j] = JE[i]*JE[j];

    // ----- H = kE*A + kL*B  → ||H||_F -----
    long double sumsq = 0.0L;
    for (int i=0;i<36;++i){
        const long double hij = (long double)kE * (long double)A[i]
                              + (long double)kL * (long double)B[i];
        sumsq += hij*hij;
    }
    double HF = std::sqrt((double)sumsq);
    return HF;
}
