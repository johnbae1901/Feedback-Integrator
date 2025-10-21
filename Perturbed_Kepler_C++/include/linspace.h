#ifndef LINSPACE_H
#define LINSPACE_H

#include <vector>

/**
 * @brief Generates a vector of linearly spaced values from start to end.
 * @param start     The starting value.
 * @param end       The ending value.
 * @param numPoints Number of points to generate.
 * @return std::vector<double> containing the generated values.
 */
std::vector<double> linspace(double start, double end, int numPoints);

#endif // LINSPACE_H
