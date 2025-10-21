#ifndef BASIC_OPERATIONS_H
#define BASIC_OPERATIONS_H

#include <Eigen/Dense>
#include <vector>

std::vector<double> linspace(double start, double end, int numPoints);
double norm3(const double v[3]);
double norm3Square(const double v[3]);
void cross3(const double a[3], const double b[3], double c[3]);
void crossProdMtx(const double v[3], double M[9]);
void getSkew(const double Omega[3], double skew[3][3]);
double matrix2norm(const Eigen::Matrix<double, 6, 6>& M);
double maxEigenvalue(const double* M, int n);
void invert3x3(const double A[3][3], double A_inv[3][3]);
void invertSymmetric3x3(const double A[3][3], double A_inv[3][3]);
void Pi(const double R[3][3], const double Omega[3], const double I[3][3], double out[3]);
double E(const double Omega[3], const double I[3][3]);

#endif