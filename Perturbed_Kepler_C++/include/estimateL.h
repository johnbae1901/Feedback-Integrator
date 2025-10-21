#ifndef ESTIMATE_L_H
#define ESTIMATE_L_H

#include <vector>
#include <tuple>

double estimateL(const std::array<double,6>& x,
                 double mu,
                 double k1,
                 double k2);

#endif // ESTIMATE_L_H
