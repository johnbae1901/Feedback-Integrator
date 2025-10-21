#ifndef GET_ERROR_H
#define GET_ERROR_H

#include <array>
#include <tuple>
#include <vector>

std::vector<double> getError(
    const std::vector<double>& x1,
    const std::vector<double>& x2,
    const std::vector<double>& x3,
    const std::vector<double>& v1,
    const std::vector<double>& v2,
    const std::vector<double>& v3,
    const std::vector<double>& L0,
    const std::vector<double>& A0,
    double k1,
    double k2,
    double mu);

double getError(
    const double& x1,
    const double& x2,
    const double& x3,
    const double& v1,
    const double& v2,
    const double& v3,
    const std::vector<double>& L0,
    const std::vector<double>& A0,
    double k1,
    double k2,
    double mu);

double getError(
    const double& x1,
    const double& x2,
    const double& x3,
    const double& v1,
    const double& v2,
    const double& v3,
    const std::array<double,3>& L0,
    const std::array<double,3>& A0,
    double k1,
    double k2,
    double mu);

#endif // GET_ERROR_H
