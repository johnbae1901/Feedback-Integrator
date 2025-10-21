#include <iostream>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <math.h>
#include <chrono>
#include "include/basic_operations.h"
#include "include/euler_vanilla.h"
#include "include/euler_feedback.h"
#include "include/euler_feedback_adaptive.h"
#include "include/euler_feedback_adaptive_light.h"
#include "include/stormer_verlet_B.h"
#include "include/getError.h"
#include "include/linspace.h"

using namespace std;

int main()
{
    // Parameters
    const int numOfIter = 5;
    const double h_ratio = 1.0;
    // const int m = 20;
    const double L_update_period = 0.1;
    const double lambda = 1.1;
    const vector<double> log_h = linspace(-1, -3, numOfIter);

    // const vector<double> xi = {1.0, 0.0, 0.0, 0.0, sqrt(1.8), 0.0}; // initial [x1,x2,x3, v1,v2,v3]
    const vector<double> xi = {1.0, 0.0, 0.0, 0.0, 0.13, 0.0}; // initial [x1,x2,x3, v1,v2,v3]
    const double ri_array[3] = { xi[0], xi[1], xi[2] };
    const double vi_array[3] = { xi[3], xi[4], xi[5] };
    const double T = 70.2481;
    const double tf = 2000*T;    // final time
    const double mu = 1.0;
    const double k1 = 4, k2 = 2;
    const double L = 487.69;
    const double c = 1.0; // TODO: delete this variable.

    double maxErrorVanilla[numOfIter], maxErrorFeedback[numOfIter], maxErrorAdaptive[numOfIter], maxErrorAdaptiveLight[numOfIter], maxErrorSVB[numOfIter];

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
        int m = static_cast<int>(ceil(L_update_period/h));
        cout << "----------------------------------------------------------------------------" << endl;
        cout << "h : " << h << ", "<< "h_ratio : " << h_ratio << endl;
        #pragma omp parallel sections shared(maxErrorVanilla, maxErrorFeedback, maxErrorAdaptive, maxErrorAdaptiveLight, maxErrorSVB)
        {
            #pragma omp section
            {
            // Run Vanilla Euler integrator
            chrono::system_clock::time_point start = chrono::system_clock::now();
            double errorVanilla = euler_vanilla_error(xi, tf, h_ratio * h, L0_array, A0_array, k1, k2, mu);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Vanilla Euler method : " << duration.count() << endl;
            maxErrorVanilla[n] = move(errorVanilla);
            }

            #pragma omp section
            {
            // Run Feedback Euler integrator
            chrono::system_clock::time_point start = chrono::system_clock::now();
            double errorFeedback = euler_feedback_error(xi, tf, h_ratio * h, mu, alpha, k1, k2, L0_array, A0_array);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Feedback Euler method : " << duration.count() << endl;
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
            double errorAdaptiveLight = euler_feedback_adaptive_error_light(xi, tf, h, mu, c, k1, k2, L0_array, A0_array, m, lambda);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Adaptive Feedback Euler method (Light) : " << duration.count() << endl;
            maxErrorAdaptiveLight[n] = move(errorAdaptiveLight);
            }

            #pragma omp section
            {
            // Run Störmer-Verlet-B
            chrono::system_clock::time_point start = chrono::system_clock::now();
            double errorSVB = stormer_verlet_B_error(xi, tf, h_ratio * h, L0_array, A0_array, k1, k2, mu);
            chrono::duration<double> duration = chrono::system_clock::now() - start;
            cout << "CPU time for Stormer-Verlet method (B) : " << duration.count() << endl;
            maxErrorSVB[n] = move(errorSVB);
            }
        }
    }
    
    // Print Results
    cout << "----------------------------------------------------------------------------" << endl;
    cout << "Max Error Vanilla          : " << maxErrorVanilla[0] << ", " <<  maxErrorVanilla[1] << 
                                            ", " <<  maxErrorVanilla[2] << ", " <<  maxErrorVanilla[3] <<
                                            ", " <<  maxErrorVanilla[4] << "\n";
    cout << "Max Error SVB              : " << maxErrorSVB[0] << ", " << maxErrorSVB[1] << 
                                            ", " << maxErrorSVB[2] << ", " << maxErrorSVB[3] <<
                                            ", " << maxErrorSVB[4] << "\n";
    cout << "Max Error Feedback         : " << maxErrorFeedback[0] << ", " << maxErrorFeedback[1] <<
                                            ", " << maxErrorFeedback[2] << ", " << maxErrorFeedback[3] <<
                                            ", " << maxErrorFeedback[4] << "\n";
    // cout << "Max Error Adaptive :      " << maxErrorAdaptive[0] << ", " << maxErrorAdaptive[1] <<
        //                                 ", " << maxErrorAdaptive[2] << ", " << maxErrorAdaptive[3] <<
        //                                 ", " << maxErrorAdaptive[4] <<"\n";
    cout << "Max Error Adaptive (Light) : " << maxErrorAdaptiveLight[0] << ", " << maxErrorAdaptiveLight[1] <<
                                            ", " << maxErrorAdaptiveLight[2] << ", " << maxErrorAdaptiveLight[3] <<
                                            ", " << maxErrorAdaptiveLight[4] <<"\n";

    return 0;
}
