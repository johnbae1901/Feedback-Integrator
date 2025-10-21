#include "basic_operations.h"
#include <cmath>
#include <Eigen/Dense>
#include <vector>

std::vector<double> linspace(double start, double end, int numPoints) {
    std::vector<double> result;
    result.reserve(numPoints);
    
    // Edge cases: if numPoints <= 0, return empty vector
    if (numPoints <= 0) {
        return result;
    }
    
    // If numPoints == 1, return only the start value
    if (numPoints == 1) {
        result.push_back(start);
        return result;
    }
    
    double step = (end - start) / (static_cast<double>(numPoints) - 1);

    for (int i = 0; i < numPoints; ++i) {
        result.push_back(start + i * step);
    }
    
    // Ensure the last value is exactly 'end' if floating-point inaccuracies occur
    result.back() = end;
    
    return result;
}

double norm3(const double v[3])
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double norm3Square(const double v[3])
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
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

// Create the 3x3 skew-symmetric matrix for a 3-vector Omega.
void getSkew(const double Omega[3], double skew[3][3]) {
    skew[0][0] = 0.0;         skew[0][1] = -Omega[2];   skew[0][2] =  Omega[1];
    skew[1][0] =  Omega[2];     skew[1][1] = 0.0;         skew[1][2] = -Omega[0];
    skew[2][0] = -Omega[1];     skew[2][1] =  Omega[0];   skew[2][2] = 0.0;
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

// Invert a 3x3 matrix A (assumes A is non-singular).
void invert3x3(const double A[3][3], double A_inv[3][3]) {
    double det = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
               - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
               + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
    double invDet = 1.0 / det;
    
    A_inv[0][0] =  (A[1][1]*A[2][2]-A[1][2]*A[2][1]) * invDet;
    A_inv[0][1] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1]) * invDet;
    A_inv[0][2] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1]) * invDet;
    
    A_inv[1][0] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0]) * invDet;
    A_inv[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0]) * invDet;
    A_inv[1][2] = -(A[0][0]*A[1][2]-A[0][2]*A[1][0]) * invDet;
    
    A_inv[2][0] =  (A[1][0]*A[2][1]-A[1][1]*A[2][0]) * invDet;
    A_inv[2][1] = -(A[0][0]*A[2][1]-A[0][1]*A[2][0]) * invDet;
    A_inv[2][2] =  (A[0][0]*A[1][1]-A[0][1]*A[1][0]) * invDet;
}

// Optimized inversion for a symmetric 3x3 matrix.
void invertSymmetric3x3(const double A[3][3], double A_inv[3][3]) {
    // Extract the unique elements of the symmetric matrix.
    double a = A[0][0];
    double b = A[0][1]; // same as A[1][0]
    double c = A[0][2]; // same as A[2][0]
    double d = A[1][1];
    double e = A[1][2]; // same as A[2][1]
    double f = A[2][2];
    
    // Compute the determinant.
    double det = a*(d * f - e * e) - b*(b * f - c * e) + c*(b * e - c * d);
    double invDet = 1.0 / det;
    
    // Compute the cofactors divided by the determinant.
    A_inv[0][0] = (d * f - e * e) * invDet;
    A_inv[0][1] = (c * e - b * f) * invDet;
    A_inv[0][2] = (b * e - c * d) * invDet;
    A_inv[1][0] = A_inv[0][1];   // symmetry
    A_inv[1][1] = (a * f - c * c) * invDet;
    A_inv[1][2] = (b * c - a * e) * invDet;
    A_inv[2][0] = A_inv[0][2];   // symmetry
    A_inv[2][1] = A_inv[1][2];   // symmetry
    A_inv[2][2] = (a * d - b * b) * invDet;
}

// Computes Pi = R * I * Omega.
// R is a 3x3 matrix, I is a 3x3 matrix, Omega is a 3-vector.
// The result (a 3-vector) is stored in the 'out' parameter.
void Pi(const double R[3][3], const double Omega[3], const double I[3][3], double out[3]) {
    double temp[3] = {0};
    // First, compute temp = I * Omega.
    for (int i = 0; i < 3; i++) {
        temp[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            temp[i] += I[i][j] * Omega[j];
        }
    }
    // Then, compute out = R * temp.
    for (int i = 0; i < 3; i++) {
        out[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            out[i] += R[i][j] * temp[j];
        }
    }
}

// Computes E = 1/2 * Omega' * I * Omega.
// Omega is a 3-vector and I is a 3x3 matrix.
double E(const double Omega[3], const double I[3][3]) {
    double temp[3] = {0};
    for (int i = 0; i < 3; i++) {
        temp[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            temp[i] += I[i][j] * Omega[j];
        }
    }
    double dot = 0.0;
    for (int i = 0; i < 3; i++) {
        dot += Omega[i] * temp[i];
    }
    return 0.5 * dot;
}