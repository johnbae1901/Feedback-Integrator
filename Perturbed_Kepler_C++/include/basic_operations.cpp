#include "basic_operations.h"
#include <cmath>
#include <Eigen/Dense>

double norm3(const double v[3])
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

//------------------------------------
// Helper: crossProdMtx
//------------------------------------
/**
 * @brief Build the 3x3 skew-symmetric matrix from a 3D vector v.
 *        M = [  0   -v[2]  v[1]
 *               v[2]   0  -v[0]
 *              -v[1]  v[0]  0  ]
 * We store M in row-major order in a double[9].
 */
void crossProdMtx(const double v[3], double M[9])
{
    // row-major layout
    // M[row*3 + col]
    M[0] = 0.0;    M[1] = -v[2];  M[2] =  v[1];
    M[3] =  v[2];  M[4] =  0.0;   M[5] = -v[0];
    M[6] = -v[1];  M[7] =  v[0];  M[8] =  0.0;
}


// For a fixed-size 6x6 matrix
double matrix2norm(const Eigen::Matrix<double, 6, 6>& M)
{
    // For fixed-size matrices, use full U/V (thin SVD isn't available).
    Eigen::JacobiSVD<Eigen::Matrix<double,6,6>> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // The singular values are sorted in decreasing order, so the first is the largest.
    return svd.singularValues()(0);
}

double maxEigenvalue(const double* M, int n)
{
    // Map the raw pointer to an Eigen matrix.
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mat(M, n, n);
    
    // For symmetric matrices, SelfAdjointEigenSolver is fast and accurate.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(mat);
    if(eigensolver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue computation failed.");
    }
    // The eigenvalues are sorted in ascending order.
    return eigensolver.eigenvalues().maxCoeff();
}
