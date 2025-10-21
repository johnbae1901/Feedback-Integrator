#include "basic_operations.h"  // for cross3(...) and norm3(...)
#include <cmath>               // for std::ceil, pow
#include <vector>
#include <tuple>
#include <iostream>            // (Optional) for debug prints
#include <array>
#include <cassert>
#include "estimateL.h"

using namespace std;


//------------------------------------
// estimateL
//------------------------------------
double estimateL(const std::array<double,6>& x,
                 double mu,
                 double k1,
                 double k2)
{
    // x = [x1, x2, x3, v1, v2, v3]
    double r[3] = { x[0], x[1], x[2] };
    double v[3] = { x[3], x[4], x[5] };

    //-----------------------------------
    // crossProdMtx(r), crossProdMtx(v)
    //-----------------------------------
    double xMtx[9]; // crossProdMtx(r)
    crossProdMtx(r, xMtx);

    double vMtx[9]; // crossProdMtx(v)
    crossProdMtx(v, vMtx);

    // Build the 6x6 matrix M in one step using row-major order:
    // M = [ -vMtx,    0 ]
    //     [   0,    xMtx ]
    double M6x6[36] = {0.0}; // Zero-initialized

    // Top-left 3x3 block: -vMtx
    for (int row = 0; row < 3; ++row)
    {
        int base = row * 6;
        int src = row * 3;
        M6x6[base + 0] = -vMtx[src + 0];
        M6x6[base + 1] = -vMtx[src + 1];
        M6x6[base + 2] = -vMtx[src + 2];
    }
    // Bottom-right 3x3 block: xMtx
    for (int row = 0; row < 3; ++row)
    {
        int base = (row + 3) * 6 + 3;
        int src = row * 3;
        M6x6[base + 0] = xMtx[src + 0];
        M6x6[base + 1] = xMtx[src + 1];
        M6x6[base + 2] = xMtx[src + 2];
    }

    // Instead of doing a separate transpose and multiplication, compute (M^T * M)
    // directly, taking advantage of symmetry and the known matrix structure.
    // Let MTM = M^T * M. We will find DV1 = k1 * MTM.
    // double MTM[36] = {0.0};
    double DV1[36];

    // We are computing for each element MTM[i*6+j] = sum_k M[k*6+i] * M[k*6+j].
    // Unroll the loops for a 6x6 fixed size.
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i; j < 6; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < 6; ++k)
            {
                sum += M6x6[k*6 + i] * M6x6[k*6 + j];
            }
            DV1[i*6 + j] = k1 * sum;
            DV1[j*6 + i] = DV1[i*6 + j]; // use symmetry, assign the mirrored element
        }
    }

    //-----------------------------------
    // Next compute DV2 = k2 * (DA'^T * DA')
    // Where DA is the block [ leftBlock, rightBlock ] (3×6).
    //-----------------------------------
    // leftBlock = v^2 * I - v*v^T - mu/||r|| * (I - (r*r^T)/||r||^2)
    // rightBlock = 2*r*v^T - v*r^T - (v^T*r)*I
    // Then DA is 3×6 => so DA^T is 6×3 => DA'^T * DA' is a 6×6

    // 1) leftBlock (3×3)
    double v2   = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double rNorm = norm3(r);
    double r2   = rNorm * rNorm;
    
    // Identity matrix (3x3) in row-major order.
    const double I3[9] = { 1, 0, 0,
                           0, 1, 0,
                           0, 0, 1 };

    // 1) Compute leftBlock (3×3)
    // Formula: leftBlock = v^2 * I - v*v^T - mu/||r|| * (I - (r*r^T)/||r||^2)
    double leftBlock[9];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            // Compute each element explicitly:
            // v2*δ_ij - (v[i]*v[j]) - (mu/rNorm)*(δ_ij - (r[i]*r[j])/r2)
            leftBlock[i*3 + j] = v2 * I3[i*3 + j]
                                - v[i] * v[j]
                                - (mu / rNorm) * ( I3[i*3 + j] - (r[i]*r[j]) / r2 );
        }
    }

    // 2) Compute rightBlock (3×3)
    // Formula: rightBlock = 2*r*v^T - v*r^T - (v^T*r)*I
    double dotVR = r[0]*v[0] + r[1]*v[1] + r[2]*v[2];
    double rightBlock[9];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            // 2*r[i]*v[j] - v[i]*r[j] - dotVR*δ_ij
            rightBlock[i*3 + j] = 2.0 * r[i] * v[j] - v[i] * r[j] - dotVR * I3[i*3 + j];
        }
    }

    // 3) Build DA (3×6) as [ leftBlock, rightBlock ]
    // We store DA in row-major order: 3 rows and 6 columns.
    double DA[18]; // 3 * 6 = 18 elements.
    for (int row = 0; row < 3; ++row)
    {
        // First half: leftBlock row.
        DA[row*6 + 0] = leftBlock[row*3 + 0];
        DA[row*6 + 1] = leftBlock[row*3 + 1];
        DA[row*6 + 2] = leftBlock[row*3 + 2];
        // Second half: rightBlock row.
        DA[row*6 + 3] = rightBlock[row*3 + 0];
        DA[row*6 + 4] = rightBlock[row*3 + 1];
        DA[row*6 + 5] = rightBlock[row*3 + 2];
    }

    // 4) Compute (DA^T * DA)  (a 6×6 matrix)
    // Instead of explicitly forming the transpose of DA, compute:
    // (DATDA)[i,j] = sum_{k=0}^{2} DA[k, i] * DA[k, j]
    // where DA[k, i] is DA[k*6 + i] in row-major format.
    // double DATDA[36] = { 0.0 };
    // 5) Compute DV2 = k2 * (DA^T * DA) and combine with DV1 into D2V.
    double D2V[36], DV2[36];
    double L = 0.0;

    for (int i = 0; i < 6; ++i)
    {
        for (int j = i; j < 6; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k)
            {
                sum += DA[k*6 + i] * DA[k*6 + j];
            }
            double dv2_val = k2 * sum;
            if (i == j)
            {
                // For diagonal elements, update once.
                DV2[i*6 + i] = dv2_val;
                D2V[i*6 + i] = DV1[i*6 + i] + dv2_val;
                L += D2V[i*6 + i] * D2V[i*6 + i];
            }
            else
            {
                // For off-diagonal elements, assign symmetrically.
                DV2[i*6 + j] = dv2_val;
                DV2[j*6 + i] = dv2_val;

                D2V[i*6 + j] = DV1[i*6 + j] + dv2_val;
                D2V[j*6 + i] = DV1[j*6 + i] + dv2_val;

                L += D2V[i*6 + j] * D2V[i*6 + j] + D2V[j*6 + i] * D2V[j*6 + i];
            }
        }
    }

    // Here we do the *Frobenius* norm for demonstration
    // If you truly need spectral norm, do SVD or use a linear algebra library
    // double L = maxEigenvalue(D2V, 6);
    // for(int i=0; i<36; i++){
    //     L += D2V[i] * D2V[i];
    // }
    return sqrt(L);
}