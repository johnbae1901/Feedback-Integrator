#ifndef BASIC_OPERATIONS_H
#define BASIC_OPERATIONS_H

#include <Eigen/Dense>

double norm3(const double v[3]);
void cross3(const double a[3], const double b[3], double c[3]);
void crossProdMtx(const double v[3], double M[9]);
double matrix2norm(const Eigen::Matrix<double, 6, 6>& M);
double maxEigenvalue(const double* M, int n);

#endif