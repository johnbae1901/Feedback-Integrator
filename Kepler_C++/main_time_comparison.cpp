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
    double h_ratio = 1.0;
    const double lambda = 1.1;
    const int m = 10;
    vector<double> xi = {1.0, 0.0, 0.0, 0.0, sqrt(1.8), 0.0}; // initial [x1,x2,x3, v1,v2,v3]
    double ri_array[3] = { xi[0], xi[1], xi[2] };
    double vi_array[3] = { xi[3], xi[4], xi[5] };
    double T = 70.2481;
    double tf = 1000*T;    // final time
    double h  = 1e-2;    // step size
    double mu = 1.0;
    double k1 = 4, k2 = 2;
    double L = 487.69;
    double alpha = 1 / (h * L);
    double c = 1.0;
    const double Hmin = 1e-10;

    double maxErrorVanilla, maxErrorFeedback, maxErrorAdaptive, maxErrorAdaptiveLight, maxErrorSVB;

    // Compute L0
    double L0_double[3];
    cross3(ri_array, vi_array, L0_double);
    vector<double> L0(L0_double, L0_double + 3);
    array<double,3> L0_array = {L0_double[0], L0_double[1], L0_double[2]};

    // Compute A0
    double cross_v_L0[3];
    cross3(vi_array, L0_double, cross_v_L0);
    double rNorm = norm3(ri_array);
    double A0_double[3];
    A0_double[0] = cross_v_L0[0] - mu * ri_array[0] / rNorm;
    A0_double[1] = cross_v_L0[1] - mu * ri_array[1] / rNorm;
    A0_double[2] = cross_v_L0[2] - mu * ri_array[2] / rNorm;
    vector<double> A0(A0_double, A0_double + 3);
    array<double,3> A0_array = {A0_double[0], A0_double[1], A0_double[2]};

    chrono::system_clock::time_point start = chrono::system_clock::now();
    chrono::duration<double> duration;

    cout << "----------------------------------------------------------------------------" << endl;
    #pragma omp parallel sections
    {
        #pragma omp section
        {
        // Find time used for Vanilla Euler method
        maxErrorVanilla = euler_vanilla_error(xi, tf, h_ratio * h, L0_array, A0_array, k1, k2, mu);
        duration = chrono::system_clock::now() - start;
        cout << "CPU time for Vanilla Euler method : " << duration.count() << endl;
        }

        #pragma omp section
        {
        // Find time used for Feedback Euler method
        start = chrono::system_clock::now();
        maxErrorFeedback = euler_feedback_error(xi, tf, h_ratio * h, mu, alpha, k1, k2, L0_array, A0_array);
        duration = chrono::system_clock::now() - start;
        cout << "CPU time for Feedback Euler method : " << duration.count() << endl;
        }

        #pragma omp section
        {
        // Find time used for Adaptive Feedback Euler method
        start = chrono::system_clock::now();
        maxErrorAdaptive = euler_feedback_adaptive_error(xi, tf, h, mu, c, k1, k2, L0_array, A0_array, lambda);
        duration = chrono::system_clock::now() - start;
        cout << "CPU time for Adaptive Feedback Euler method : " << duration.count() << endl;
        }

        #pragma omp section
        {
        // Find time used for Adaptive Feedback Euler method
        start = chrono::system_clock::now();
        maxErrorAdaptiveLight = euler_feedback_adaptive_error_light(xi, tf, h, mu, k1, k2, L0_array, A0_array, m, lambda, Hmin);
        duration = chrono::system_clock::now() - start;
        cout << "CPU time for Adaptive Feedback Euler method (light) : " << duration.count() << endl;
        }

        #pragma omp section
        {
        // Find time used for Stormer-Verlet method (B)
        start = chrono::system_clock::now();
        maxErrorSVB = stormer_verlet_B_error(xi, tf, h_ratio * h, L0_array, A0_array, k1, k2, mu);
        duration = chrono::system_clock::now() - start;
        cout << "CPU time for Stormer-Verlet method (B) : " << duration.count() << endl;
        }
    }
    // Print Results
    cout << "----------------------------------------------------------------------------" << endl;
    cout << "Max Error Vanilla: " << maxErrorVanilla << "\n";
    cout << "Max Error Feedback: " << maxErrorFeedback << "\n";
    cout << "Max Error Adaptive: " << maxErrorAdaptive <<"\n";
    cout << "Max Error Adaptive (Light): " << maxErrorAdaptiveLight <<"\n";
    cout << "Max Error SVB: " << maxErrorSVB << "\n";
    return 0;
}
