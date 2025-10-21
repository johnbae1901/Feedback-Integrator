#include <iostream>
#include "linspace.h"
#include <vector>
#include <iomanip> // for std::fixed, std::setprecision

/**
 * @brief Generates a vector of linearly spaced values from start to end.
 * @param start     The starting value.
 * @param end       The ending value.
 * @param numPoints Number of points to generate.
 * @return std::vector<double> containing the generated values.
 */
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

// int main() {
//     // Example usage:
//     double start = 0.0;
//     double end  = 10.0;
//     int N        = 5;
    
//     std::vector<double> values = linspace(start, end, N);
    
//     // Print results
//     std::cout << std::fixed << std::setprecision(3);
//     for (const auto& val : values) {
//         std::cout << val << " ";
//     }
//     std::cout << std::endl;
    
//     return 0;
// }
