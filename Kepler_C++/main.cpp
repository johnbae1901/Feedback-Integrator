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
#include "include/euler_feedback_adaptive.h"
#include "include/euler_feedback_adaptive_light.h"
#include "include/stormer_verlet_B.h"
#include "include/getError.h"
#include "include/linspace.h"

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
    const vector<double> log_h = linspace(-1, -5, numOfIter);

    const vector<double> xi = {1.0, 0.0, 0.0, 0.0, sqrt(1.8), 0.0}; // initial [x1,x2,x3, v1,v2,v3]
    const double ri_array[3] = { xi[0], xi[1], xi[2] };
    const double vi_array[3] = { xi[3], xi[4], xi[5] };
    const double T = 70.2481;
    const double tf = 2000*T;    // final time
    const double mu = 1.0;
    const double k1 = 4, k2 = 2;
    const double L = 487.69;
    const double c = 1.0; // TODO: delete this variable.
    const double Hmin = 1e-10;

    double maxErrorVanilla[numOfIter], maxErrorVanillaFeedback[numOfIter], maxErrorFeedback[numOfIter], maxErrorAdaptive[numOfIter], maxErrorAdaptiveLight[numOfIter], maxErrorSVB[numOfIter];
    vector<double> t_vanilla(numOfIter), t_feedback_vanilla(numOfIter), t_feedback(numOfIter), t_adaptive(numOfIter), t_SV(numOfIter);

    // Compute L0
    double L0_double[3];
    cross3(ri_array, vi_array, L0_double);
    const vector<double> L0(L0_double, L0_double + 3);
    const array<double,3> L0_array = {L0_double[0], L0_double[1], L0_double[2]};

    // Compute A0
    double cross_v_L0[3];
    cross3(vi_array, L0_double, cross_v_L0);
    double riNorm = norm3(ri_array);
    double A0_double[3];
    A0_double[0] = cross_v_L0[0] - mu * ri_array[0] / riNorm;
    A0_double[1] = cross_v_L0[1] - mu * ri_array[1] / riNorm;
    A0_double[2] = cross_v_L0[2] - mu * ri_array[2] / riNorm;
    const vector<double> A0(A0_double, A0_double + 3);
    const array<double,3> A0_array = {A0_double[0], A0_double[1], A0_double[2]};

    for (int n = 0; n < numOfIter; n++) {
        double h = pow(10, log_h[n]);
        double alpha = 1 / (h * L);
        // double alpha = 1.0;
        int m = static_cast<int>(ceil(L_update_period/h));
        cout << "----------------------------------------------------------------------------" << endl;
        cout << "h : " << h << endl;
        #pragma omp parallel sections shared(maxErrorVanilla, maxErrorVanillaFeedback, maxErrorFeedback, maxErrorAdaptive, maxErrorAdaptiveLight, maxErrorSVB, t_vanilla, t_feedback_vanilla, t_feedback, t_adaptive, t_SV)
        {
            #pragma omp section
            {
            // Run Vanilla Euler integrator
            chrono::system_clock::time_point start = chrono::system_clock::now();
            double errorVanilla = euler_vanilla_error(xi, tf, h, L0_array, A0_array, k1, k2, mu);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Vanilla Euler method : " << duration.count() << endl;
            t_vanilla[n] = duration.count();
            maxErrorVanilla[n] = move(errorVanilla);
            }

            #pragma omp section
            {
            // Run Vanilla Feedback Euler integrator
            chrono::system_clock::time_point start = chrono::system_clock::now();
            double errorVanillaFeedback = euler_feedback_error(xi, tf, h, mu, 1.0, k1, k2, L0_array, A0_array);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Vanilla Feedback Euler method : " << duration.count() << endl;
            t_feedback_vanilla[n] = duration.count();
            maxErrorVanillaFeedback[n] = move(errorVanillaFeedback);
            }

            #pragma omp section
            {
            // Run Feedback Euler integrator
            chrono::system_clock::time_point start = chrono::system_clock::now();
            double errorFeedback = euler_feedback_error(xi, tf, h, mu, alpha, k1, k2, L0_array, A0_array);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Feedback Euler method : " << duration.count() << endl;
            t_feedback[n] = duration.count();
            maxErrorFeedback[n] = move(errorFeedback);
            }

            // #pragma omp section
            // {
            // // Run Adaptive Feedback Integrator
            // // auto [X1_adaptive, X2_adaptive, X3_adaptive, V1_adaptive, V2_adaptive, V3_adaptive, T_adaptive] =
            // //     euler_feedback_adaptive(xi, tf, h, mu, c, k1, k2, L0, A0);
            // // vector<double> errorAdaptive = getError(X1_adaptive, X2_adaptive, X3_adaptive, 
            // //                                     V1_adaptive, V2_adaptive, V3_adaptive, 
            // //                                     L0, A0, k1, k2, mu);
            // // maxErrorAdaptive = move(*max_element(errorAdaptive.begin(), errorAdaptive.end()));
            // chrono::system_clock::time_point start = chrono::system_clock::now();
            // double errorAdaptive = euler_feedback_adaptive_error(xi, tf, h, mu, c, k1, k2, L0_array, A0_array);
            // chrono::duration<double> duration = chrono::system_clock::now() - start;
            // cout << "CPU time for Adaptive Feedback Euler method : " << duration.count() << endl;
            // maxErrorAdaptive[n] = move(errorAdaptive);
            // }

            #pragma omp section
            {
            // Run Adaptive Feedback Integrator (Light)
            chrono::system_clock::time_point start = chrono::system_clock::now();
            double errorAdaptiveLight = euler_feedback_adaptive_error_light(xi, tf, h, mu, c, k1, k2, L0_array, A0_array, m, lambda, Hmin);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Adaptive Feedback Euler method (Light) : " << duration.count() << endl;
            t_adaptive[n] = duration.count();
            maxErrorAdaptiveLight[n] = move(errorAdaptiveLight);
            }

            #pragma omp section
            {
            // Run Störmer-Verlet-B
            chrono::system_clock::time_point start = chrono::system_clock::now();
            double errorSVB = stormer_verlet_B_error(xi, tf, h, L0_array, A0_array, k1, k2, mu);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Stormer-Verlet method (B) : " << duration.count() << endl;
            t_SV[n] = duration.count();
            maxErrorSVB[n] = move(errorSVB);
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
        cout << maxErrorSVB[i] << (i+1==numOfIter?'\n':',');
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
            << maxErrorSVB[i] << "\n";
    }
    fout.close();
    std::cout << "[SAVE] Wrote error_data/errors_vs_h.csv\n";

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
        // // pick h
        // double h_save = 0.0;
        // int idx_save = 0;
        // if (save.h > 0) {
        //     // choose nearest by absolute difference in log-space
        //     double best = 1e300;
        //     for (int i=0;i<numOfIter;++i){
        //         double hi = pow(10.0, log_h[i]);
        //         double d = std::abs(std::log10(hi) - std::log10(save.h));
        //         if (d < best){ best = d; idx_save = i; h_save = hi; }
        //     }
        // } else {
        //     idx_save = numOfIter/2;
        //     h_save = pow(10.0, log_h[idx_save]);
        // }
        double h_save = save.h;
        std::cout << "[SAVE] dumping trajectories for h=" << h_save
                << " to '" << save.outdir << "'\n";

        using std::filesystem::path;

        // vanilla
        if (save.methods=="all" || save.methods.find("vanilla")!=std::string::npos){
            auto [x1, x2, x3, v1, v2, v3, tvec] = euler_vanilla(xi, tf, h_save, mu);
            save_state_csv(
                (path(save.outdir)/("vanilla/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }

        // vanilla feedback
        if (save.methods=="all" || save.methods.find("vanilla_feedback")!=std::string::npos){
            const double alpha = 1.0;
            auto [x1, x2, x3, v1, v2, v3, tvec] = euler_feedback(xi, tf, h_save, mu, alpha, k1, k2, L0_array, A0_array);
            save_state_csv(
                (path(save.outdir)/("vanilla_feedback/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }

        // feedback
        if (save.methods=="all" || save.methods.find("feedback")!=std::string::npos){
            const double alpha = 1/(h_save * L);
            auto [x1, x2, x3, v1, v2, v3, tvec] = euler_feedback(xi, tf, h_save, mu, alpha, k1, k2, L0_array, A0_array);
            save_state_csv(
                (path(save.outdir)/("feedback/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }

        // adaptive
        if (save.methods=="all" || save.methods.find("adaptive")!=std::string::npos){
            int m = static_cast<int>(ceil(L_update_period/h_save));
            auto [x1, x2, x3, v1, v2, v3, tvec] = euler_feedback_adaptive_light(xi, tf, h_save, mu, c, k1, k2, L0_array, A0_array, m, lambda, Hmin);
            save_state_csv(
                (path(save.outdir)/("adaptive/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }

        // Stormer Verlet
        if (save.methods=="all" || save.methods.find("SV")!=std::string::npos){
            auto [x1, x2, x3, v1, v2, v3, tvec] = stormer_verlet_B(xi, tf, h_save, mu);
            save_state_csv(
                (path(save.outdir)/("SV/h="+std::to_string(h_save)+"/traj.csv")).string(),
                tvec, x1, x2, x3, v1, v2, v3
            );
        }
    }

    return 0;
}
