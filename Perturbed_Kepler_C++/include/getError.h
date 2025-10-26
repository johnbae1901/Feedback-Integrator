#ifndef GET_ERROR_H
#define GET_ERROR_H

#include <array>
#include <tuple>
#include <vector>
#include "dynamics_perturbed_kepler.h"

struct LEError { double distL_sq, distE_sq, error; };

std::vector<double> getError_PK(
    const std::vector<double>& x1,
    const std::vector<double>& x2,
    const std::vector<double>& x3,
    const std::vector<double>& v1,
    const std::vector<double>& v2,
    const std::vector<double>& v3,
    const std::array<double,3>& L0,
    double E0,
    double kL,
    double kE,
    const PKParams& P);

LEError getError_PK(
    const double& x1,
    const double& x2,
    const double& x3,
    const double& v1,
    const double& v2,
    const double& v3,
    const std::array<double,3>& L0,
    double E0,
    double kL,
    double kE,
    const PKParams& P);

#endif // GET_ERROR_H
