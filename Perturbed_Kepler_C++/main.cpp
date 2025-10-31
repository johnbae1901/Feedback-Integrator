#include <iostream>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <math.h>
#include <chrono>
#include <filesystem>
#include <fstream>

#include "include/csv_io.h"
#include "include/basic_operations.h"
#include "include/euler_vanilla.h"
#include "include/euler_feedback.h"
#include "include/euler_feedback_adaptive_light.h"
#include "include/stormer_verlet_B.h"
#include "include/getError.h"
#include "include/linspace.h"
#include "include/getL0E0.h"

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
    // Parameters
    const int numOfIter = 10;
    const double L_update_period = 0.1;
    const double lambda = 1.1;
    const vector<double> log_h = linspace(-1, -9, numOfIter);

    const double e = 0.6;
    const vector<double> xi = {1.0 - e, 0.0, 0.0, 0.0, sqrt((1+e)/(1-e)), 0.0}; // initial [x1,x2,x3, v1,v2,v3]
    const double ri_array[3] = { xi[0], xi[1], xi[2] };
    const double vi_array[3] = { xi[3], xi[4], xi[5] };
    const double tf = 200;    // final time
    const double mu = 1.0;
    const double k1 = 2, k2 = 3;
    const double L = 487.69;
    const double Hmin = 1e-10;

    PKParams P; P.mu = 1.0; P.delta = 0.0025;

    getL0E0(xi, P);

    double maxErrorVanilla[numOfIter], maxErrorVanillaFeedback[numOfIter], maxErrorFeedback[numOfIter], maxErrorAdaptive[numOfIter], maxErrorAdaptiveLight[numOfIter], maxErrorSV[numOfIter];
    double maxdE_Vanilla[numOfIter], maxdE_VanillaFeedback[numOfIter], maxdE_Feedback[numOfIter], maxdE_AdaptiveLight[numOfIter], maxdE_SV[numOfIter];
    double maxdL_Vanilla[numOfIter], maxdL_VanillaFeedback[numOfIter], maxdL_Feedback[numOfIter], maxdL_AdaptiveLight[numOfIter], maxdL_SV[numOfIter];
    vector<double> t_vanilla(numOfIter), t_feedback_vanilla(numOfIter), t_feedback(numOfIter), t_adaptive(numOfIter), t_SV(numOfIter);


    for (int n = 0; n < numOfIter; n++) {
        double h = pow(10, log_h[n]);
        double alpha = 1 / (h * L);
        int m = static_cast<int>(ceil(L_update_period/h));
        cout << "----------------------------------------------------------------------------" << endl;
        cout << "h : " << h << endl;
        #pragma omp parallel sections shared( \
                maxErrorVanilla, maxErrorVanillaFeedback, maxErrorFeedback, maxErrorAdaptive, maxErrorAdaptiveLight, maxErrorSV, \
                maxdE_Vanilla, maxdE_VanillaFeedback, maxdE_Feedback, maxdE_AdaptiveLight, maxdE_SV, \
                maxdL_Vanilla, maxdL_VanillaFeedback, maxdL_Feedback, maxdL_AdaptiveLight, maxdL_SV, \
                t_vanilla, t_feedback_vanilla, t_feedback, t_adaptive, t_SV)
        {
            #pragma omp section
            {
            // Run Euler integrator
            chrono::system_clock::time_point start = chrono::system_clock::now();
            auto [error, maxdE, maxdL] = euler_vanilla_error(xi, tf, h, k1, k2, P);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Vanilla Euler method : " << duration.count() << endl;
            t_vanilla[n] = duration.count();
            maxErrorVanilla[n] = error;
            maxdE_Vanilla[n] = maxdE;
            maxdL_Vanilla[n] = maxdL;
            }

            #pragma omp section
            {
            // Run Vanilla Feedback Euler integrator
            chrono::system_clock::time_point start = chrono::system_clock::now();
            auto [error, maxdE, maxdL] = euler_feedback_error(xi, tf, h, P, 1.0, k1, k2);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Vanilla Feedback Euler method : " << duration.count() << endl;
            t_feedback_vanilla[n] = duration.count();
            maxErrorVanillaFeedback[n] = error;
            maxdE_VanillaFeedback[n] = maxdE;
            maxdL_VanillaFeedback[n] = maxdL;
            }

            #pragma omp section
            {
            // Run Feedback Euler integrator
            chrono::system_clock::time_point start = chrono::system_clock::now();
            auto [error, maxdE, maxdL] = euler_feedback_error(xi, tf, h, P, alpha, k1, k2);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Feedback Euler method : " << duration.count() << endl;
            t_feedback[n] = duration.count();
            maxErrorFeedback[n] = error;
            maxdE_Feedback[n] = maxdE;
            maxdL_Feedback[n] = maxdL;
            }

            #pragma omp section
            {
            // Run Adaptive Feedback Integrator (Light)
            chrono::system_clock::time_point start = chrono::system_clock::now();
            auto [error, maxdE, maxdL] = euler_feedback_adaptive_error_light(xi, tf, h, P, k1, k2, m, lambda, Hmin);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Adaptive Feedback Euler method (Light) : " << duration.count() << endl;
            t_adaptive[n] = duration.count();
            maxErrorAdaptiveLight[n] = error;
            maxdE_AdaptiveLight[n] = maxdE;
            maxdL_AdaptiveLight[n] = maxdL;
            }

            #pragma omp section
            {
            // Run Störmer-Verlet-B
            chrono::system_clock::time_point start = chrono::system_clock::now();
            auto [error, maxdE, maxdL] = stormer_verlet_B_error(xi, tf, h, k1, k2, P);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Stormer-Verlet method (B) : " << duration.count() << endl;
            t_SV[n] = duration.count();
            maxErrorSV[n] = error;
            maxdE_SV[n] = maxdE;
            maxdL_SV[n] = maxdL;
            }
        }
    }
    
    // ---------- Summary ----------
    cout << "============================================================================\n";
    cout << "Errors (Vanilla) : ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << maxErrorVanilla[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "Errors (Vanilla Feedback): ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << maxErrorVanillaFeedback[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "Errors (Feedback): ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << maxErrorFeedback[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "Errors (Adaptive): ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << maxErrorAdaptiveLight[i] << (i + 1 == numOfIter ? "\n" : ", ");
    }
    cout << "Errors (Stormer Verlet): ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << maxErrorSV[i] << (i+1==numOfIter?'\n':',');
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
    cout << "CPU times (s) — Stormer Verlet: ";
    for (int i = 0; i < numOfIter; ++i) {
        cout << t_SV[i] << (i+1==numOfIter?'\n':',');
    }

    // Saving error data
    std::filesystem::create_directories("results/error_data");
    std::ofstream fout("results/error_data/errors_vs_h.csv");
    fout << "h,log10_h,vanilla,vanilla_feedback,feedback,adaptive,SV\n";
    for (int i = 0; i < numOfIter; ++i) {
        double h_i = std::pow(10.0, log_h[i]);
        fout << h_i << ","
            << std::log10(h_i) << ","
            << maxErrorVanilla[i] << ","
            << maxErrorVanillaFeedback[i] << ","
            << maxErrorFeedback[i] << ","
            << maxErrorAdaptiveLight[i] << ","
            << maxErrorSV[i] << "\n";
    }
    fout.close();
    std::cout << "[SAVE] Wrote error_data/errors_vs_h.csv\n";

    {
        std::ofstream fE("results/error_data/maxdE_vs_h.csv");
        fE << "h,log10_h,vanilla,vanilla_feedback,feedback,adaptive,SV\n";
        for (int i = 0; i < numOfIter; ++i) {
            const double h_i = std::pow(10.0, log_h[i]);
            fE << h_i << ","
            << std::log10(h_i) << ","
            << maxdE_Vanilla[i] << ","
            << maxdE_VanillaFeedback[i] << ","
            << maxdE_Feedback[i] << ","
            << maxdE_AdaptiveLight[i] << ","
            << maxdE_SV[i] << "\n";
        }
        fE.close();
        std::cout << "[SAVE] Wrote error_data/maxdE_vs_h.csv\n";
    }

    {
        std::ofstream fL("results/error_data/maxdL_vs_h.csv");
        fL << "h,log10_h,vanilla,vanilla_feedback,feedback,adaptive,SV\n";
        for (int i = 0; i < numOfIter; ++i) {
            const double h_i = std::pow(10.0, log_h[i]);
            fL << h_i << ","
            << std::log10(h_i) << ","
            << maxdL_Vanilla[i] << ","
            << maxdL_VanillaFeedback[i] << ","
            << maxdL_Feedback[i] << ","
            << maxdL_AdaptiveLight[i] << ","
            << maxdL_SV[i] << "\n";
        }
        fL.close();
        std::cout << "[SAVE] Wrote error_data/maxdL_vs_h.csv\n";
    }

    // Saving CPU time data
    std::filesystem::create_directories("results/cpu_time");
    std::ofstream f("results/cpu_time/cpu_time_vs_h_all.csv");
    f << "h,log10_h,vanilla,vanilla_feedback,feedback,adaptive,SV\n";
    for (int i = 0; i < numOfIter; ++i) {
        const double h_i = std::pow(10.0, log_h[i]);
        f << h_i << ","
          << std::log10(h_i) << ","
          << t_vanilla[i] << ","
          << t_feedback_vanilla[i] << ","
          << t_feedback[i] << ","
          << t_adaptive[i] << ","
          << t_SV[i] << "\n";
    }
    f.close();
    std::cout << "[SAVE] Wrote cpu_time/cpu_time_vs_h_all.csv\n";

    // Saving trajectrories
    if (save.enable) {
        double h_save = save.h;
        std::cout << "[SAVE] dumping trajectories for h=" << h_save
                << " to '" << save.outdir << "'\n";

        using std::filesystem::path;

        // vanilla
        if (save.methods=="all" || save.methods.find("vanilla")!=std::string::npos){
            auto [x1, x2, x3, v1, v2, v3, tvec] = euler_vanilla(xi, tf, h_save, P);
            save_state_csv(
                (path(save.outdir)/("vanilla/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }

        // vanilla feedback
        if (save.methods=="all" || save.methods.find("vanilla_feedback")!=std::string::npos){
            const double alpha = 1.0;
            auto [x1, x2, x3, v1, v2, v3, tvec] = euler_feedback(xi, tf, h_save, P, alpha, k1, k2);
            save_state_csv(
                (path(save.outdir)/("vanilla_feedback/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }

        // feedback
        if (save.methods=="all" || save.methods.find("feedback")!=std::string::npos){
            const double alpha = 1/(h_save * L);
            auto [x1, x2, x3, v1, v2, v3, tvec] = euler_feedback(xi, tf, h_save, P, alpha, k1, k2);
            save_state_csv(
                (path(save.outdir)/("feedback/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }

        // adaptive
        if (save.methods=="all" || save.methods.find("adaptive")!=std::string::npos){
            int m = static_cast<int>(ceil(L_update_period/h_save));
            auto [x1, x2, x3, v1, v2, v3, tvec] = euler_feedback_adaptive_light(xi, tf, h_save, P, k1, k2, m, lambda, Hmin);
            save_state_csv(
                (path(save.outdir)/("adaptive/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }

        // Stormer Verlet
        if (save.methods=="all" || save.methods.find("SV")!=std::string::npos){
            auto [x1, x2, x3, v1, v2, v3, tvec] = stormer_verlet_B(xi, tf, h_save, P);
            save_state_csv(
                (path(save.outdir)/("SV/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }
    }

    return 0;
}
