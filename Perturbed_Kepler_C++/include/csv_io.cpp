#include "../include/csv_io.h"
#include <fstream>
#include <filesystem>
#include <algorithm>

bool save_state_csv(const std::string& csv_path,
                    const std::vector<double>& t,
                    const std::vector<double>& x1,
                    const std::vector<double>& x2,
                    const std::vector<double>& x3,
                    const std::vector<double>& v1,
                    const std::vector<double>& v2,
                    const std::vector<double>& v3)
{
    namespace fs = std::filesystem;
    fs::path p(csv_path);
    fs::create_directories(p.parent_path());

    const std::size_t N = std::min({ t.size(), x1.size(), x2.size(), x3.size(),
                                     v1.size(), v2.size(), v3.size() });

    std::ofstream f(csv_path);
    if (!f.is_open()) return false;

    f << "t,x1,x2,x3,v1,v2,v3\n";
    for (std::size_t k = 0; k < N; ++k) {
        f << t[k]  << ","
          << x1[k] << "," << x2[k] << "," << x3[k] << ","
          << v1[k] << "," << v2[k] << "," << v3[k] << "\n";
    }
    return true;
}
