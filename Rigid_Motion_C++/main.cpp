#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <filesystem>
#include <fstream>

#include "include/csv_io.h"
#include "include/basic_operations.h"
#include "include/linspace.h"
#include "src/euler_vanilla.h"
#include "src/euler_feedback.h"
#include "src/euler_adaptive.h"
#include "src/splitting_strang.h"
#include "src/splitting_strang_full.h"

using namespace std;

struct SaveOpts {
    bool enable = false;
    double h = -1.0;              // if <0, use the middle h
    std::string methods = "all";  // "vanilla,feedback,adaptive,strang" or "all"
    std::string outdir = "traj";
};
SaveOpts parse_cli(int argc, char** argv){
    SaveOpts o;
    for (int i=1;i<argc;++i){
        std::string s(argv[i]);
        auto eat = [&](const char* key, std::string& dst){
            const std::string k = std::string(key) + "=";
            if (s.rfind(k,0)==0){ dst = s.substr(k.size()); return true; }
            return false;
        };
        std::string v;
        if (eat("--save", v)) { o.enable = (v=="1" || v=="true" || v=="on"); continue; }
        if (eat("--save_h", v)) { o.h = std::stod(v); continue; }
        if (eat("--save_methods", v)) { o.methods = v; continue; }
        if (eat("--outdir", v)) { o.outdir = v; continue; }
    }
    return o;
}

int main(int argc, char** argv)
{
    SaveOpts save = parse_cli(argc, argv);
    const int    numOfIter = 10;
    const vector<double> log_h = linspace(-1.0, -7.0, numOfIter);

    const double tf = 1e3;
    const double L_update_period = 30.0;
    const double lambda = 1.1;

    const double k0 = 50.0, k1 = 100.0, k2 = 50.0;
    const double L0 = 1986.0;

    // Initial R, Omega, inertia
    const double R0[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double Omega0[3] = {1.0, 1.0, 1.0};
    const double I[3][3] = {
        {3.0, 0.0, 0.0},
        {0.0, 2.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    // Precompute invariants from initial state
    double E0 = E(Omega0, I);
    double Pi0[3];
    Pi(R0, Omega0, I, Pi0);

    // Storage for per-h results
    vector<double> err_vanilla(numOfIter), err_feedback(numOfIter), err_feedback_vanilla(numOfIter), err_adaptive(numOfIter), err_split_strang(numOfIter);
    vector<double> t_vanilla(numOfIter), t_feedback(numOfIter), t_feedback_vanilla(numOfIter), t_adaptive(numOfIter), t_split_strang(numOfIter);

    for (int n = 0; n < numOfIter; ++n) {
        const double h = pow(10.0, log_h[n]);
        const double alpha = 1/(h*L0);
        const int m = static_cast<int>(ceil(L_update_period / h));

        cout << "----------------------------------------------------------------------------\n";
        cout << "h : " << h << "  (10^" << log_h[n] << ")\n";

        #pragma omp parallel sections shared(err_vanilla, err_feedback_vanilla, err_feedback, err_adaptive, err_split_strang, t_vanilla, t_feedback_vanilla, t_feedback, t_adaptive, t_split_strang)
        {
            #pragma omp section
            {
                auto start = chrono::system_clock::now();
                double e = euler_vanilla(R0, Omega0, I, E0, Pi0, k0, k1, k2, tf, h);
                auto dur  = chrono::duration<double>(chrono::system_clock::now() - start).count();
                err_vanilla[n] = e;
                t_vanilla[n]   = dur;
                cout << "[h=" << h << "] CPU time (Vanilla Euler): " << dur << " s\n";
            }

            #pragma omp section
            {
                auto start = chrono::system_clock::now();
                // In the original main_rigid, feedback call used a fixed "1.0" parameter.
                // We keep that behavior to avoid changing algorithmic settings.
                double e = euler_feedback(R0, Omega0, I,
                                          k0, k1, k2,
                                          E0, Pi0,
                                          1.0,  // keep as-is (not switching to 1/(h*L0_const) to preserve original behavior)
                                          tf, h);
                auto dur  = chrono::duration<double>(chrono::system_clock::now() - start).count();
                err_feedback_vanilla[n] = e;
                t_feedback_vanilla[n]   = dur;
                cout << "[h=" << h << "] CPU time (Vanilla Feedback Euler): " << dur << " s\n";
            }

            #pragma omp section
            {
                auto start = chrono::system_clock::now();
                // In the original main_rigid, feedback call used a fixed "1.0" parameter.
                // We keep that behavior to avoid changing algorithmic settings.
                double e = euler_feedback(R0, Omega0, I,
                                          k0, k1, k2,
                                          E0, Pi0,
                                          alpha,  // keep as-is (not switching to 1/(h*L0_const) to preserve original behavior)
                                          tf, h);
                auto dur  = chrono::duration<double>(chrono::system_clock::now() - start).count();
                err_feedback[n] = e;
                t_feedback[n]   = dur;
                cout << "[h=" << h << "] CPU time (Feedback Euler): " << dur << " s\n";
            }

            #pragma omp section
            {
                auto start = chrono::system_clock::now();
                double e = euler_adaptive(R0, Omega0, I,
                                          k0, k1, k2,
                                          E0, Pi0,
                                          tf, h,
                                          m, lambda);
                auto dur  = chrono::duration<double>(chrono::system_clock::now() - start).count();
                err_adaptive[n] = e;
                t_adaptive[n]   = dur;
                cout << "[h=" << h << "] CPU time (Adaptive Feedback Euler): " << dur << " s\n";
            }

            #pragma omp section
            {
                auto start = chrono::system_clock::now();
                double e = splitting_strang(R0, Omega0, I, E0, Pi0, k0, k1, k2, tf, h);
                auto dur  = chrono::duration<double>(chrono::system_clock::now() - start).count();
                err_split_strang[n] = e; t_split_strang[n] = dur;
                cout << "[h=" << h << "] CPU time (Strang Splitting): " << dur << " s\n";
            }
        }
    }

    // ---------- Summary ----------
    cout << "============================================================================\n";
    cout << "Errors (Vanilla) : ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << err_vanilla[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "Errors (Vanilla Feedback): ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << err_feedback_vanilla[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "Errors (Feedback): ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << err_feedback[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "Errors (Adaptive): ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << err_adaptive[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "Errors (Strang Splitting): ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << err_split_strang[i] << (i+1==numOfIter?'\n':',');
    }

    cout << "CPU times (s) — Vanilla : ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << t_vanilla[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "CPU times (s) — Vanilla Feedback: ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << t_feedback_vanilla[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "CPU times (s) — Feedback: ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << t_feedback[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "CPU times (s) — Adaptive: ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << t_adaptive[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "CPU times (s) — Strang Splitting: ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << t_split_strang[i] << (i+1==numOfIter?'\n':',');
    }

    // Saving error data
    std::filesystem::create_directories("results/error_data");
    std::ofstream fout("results/error_data/errors_vs_h.csv");
    fout << "h,log10_h,vanilla,vanilla_feedback,feedback,adaptive,strang\n";
    for (int i = 0; i < numOfIter; ++i) {
        double h_i = std::pow(10.0, log_h[i]);
        fout << h_i << ","
            << std::log10(h_i) << ","
            << err_vanilla[i] << ","
            << err_feedback_vanilla[i] << ","
            << err_feedback[i] << ","
            << err_adaptive[i] << ","
            << err_split_strang[i] << "\n";
    }
    fout.close();
    std::cout << "[SAVE] Wrote error_data/errors_vs_h.csv\n";

    // Saving trajectrories
    if (save.enable) {
        // pick h
        double h_save = 0.0;
        int idx_save = 0;
        if (save.h > 0) {
            // choose nearest by absolute difference in log-space
            double best = 1e300;
            for (int i=0;i<numOfIter;++i){
                double hi = pow(10.0, log_h[i]);
                double d = std::abs(std::log10(hi) - std::log10(save.h));
                if (d < best){ best = d; idx_save = i; h_save = hi; }
            }
        } else {
            idx_save = numOfIter/2;
            h_save = pow(10.0, log_h[idx_save]);
        }
        // double h_save = save.h;
        std::cout << "[SAVE] dumping trajectories for h=" << h_save
                << " to '" << save.outdir << "'\n";

        using std::filesystem::path;

        // vanilla
        if (save.methods=="all" || save.methods.find("vanilla")!=std::string::npos){
            double *R=nullptr,*Om=nullptr,*t=nullptr; int N=0;
            euler_vanilla(R0, Omega0, I, tf, h_save, R, Om, t, N);
            save_state_csv( (path(save.outdir)/("vanilla/h="+std::to_string(h_save)+"/traj.csv")).string(),
                            t, Om, R, N);
            delete[] R; delete[] Om; delete[] t;
        }

        // vanilla feedback
        if (save.methods=="all" || save.methods.find("vanilla_feedback")!=std::string::npos){
            double *R=nullptr,*Om=nullptr,*t=nullptr; int N=0;
            const double alpha = 1.0;
            euler_feedback(R0, Omega0, I, k0, k1, k2, E0, Pi0, alpha, tf, h_save, R, Om, t, N);
            save_state_csv( (path(save.outdir)/("vanilla_feedback/h="+std::to_string(h_save)+"/traj.csv")).string(),
                            t, Om, R, N);
            delete[] R; delete[] Om; delete[] t;
        }

        // feedback
        if (save.methods=="all" || save.methods.find("feedback")!=std::string::npos){
            double *R=nullptr,*Om=nullptr,*t=nullptr; int N=0;
            const double alpha = 1/(h_save * L0);
            euler_feedback(R0, Omega0, I, k0, k1, k2, E0, Pi0, alpha, tf, h_save, R, Om, t, N);
            save_state_csv( (path(save.outdir)/("feedback/h="+std::to_string(h_save)+"/traj.csv")).string(),
                            t, Om, R, N);
            delete[] R; delete[] Om; delete[] t;
        }

        // adaptive
        if (save.methods=="all" || save.methods.find("adaptive")!=std::string::npos){
            double *R=nullptr,*Om=nullptr,*t=nullptr; int N=0;
            int m = 30; /* update Lipschitz estimate every m steps; adjust as you like */
            euler_adaptive(R0, Omega0, I, k0, k1, k2, E0, Pi0, tf, h_save, m, 1.1 /*lambda*/,
                        R, Om, t, N);
            save_state_csv( (path(save.outdir)/("adaptive/h="+std::to_string(h_save)+"/traj.csv")).string(),
                            t, Om, R, N);
            delete[] R; delete[] Om; delete[] t;
        }

        // strang
        if (save.methods=="all" || save.methods.find("strang")!=std::string::npos){
            double *R=nullptr,*Om=nullptr,*t=nullptr; int N=0;
            splitting_strang_full(R0, Omega0, I, tf, h_save, R, Om, t, N);
            save_state_csv( (path(save.outdir)/("strang/h="+std::to_string(h_save)+"/traj.csv")).string(),
                            t, Om, R, N);
            delete[] R; delete[] Om; delete[] t;
        }
    }

    return 0;
}
