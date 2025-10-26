#pragma once
#include <string>
#include <vector>

// Save trajectories into a single CSV with header:
// t,x1,x2,x3,v1,v2,v3
// All vectors must have the same length N (or we take the minimum length).
bool save_state_csv(const std::string& csv_path,
                    const std::vector<double>& t,
                    const std::vector<double>& x1,
                    const std::vector<double>& x2,
                    const std::vector<double>& x3,
                    const std::vector<double>& v1,
                    const std::vector<double>& v2,
                    const std::vector<double>& v3);
