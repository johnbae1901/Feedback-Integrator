#ifndef GET_ERROR_H
#define GET_ERROR_H

#include <array>
#include <tuple>
#include <vector>

struct LAError { double distL_sq, distA_sq, error; };

LAError getError(
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
