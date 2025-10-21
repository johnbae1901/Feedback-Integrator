#include <iostream>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <math.h>
#include <chrono>
#include "include/basic_operations.h"
#include "include/getError.h"
#include "src/euler_vanilla.h"
#include "src/euler_feedback.h"
#include "src/euler_adaptive.h"

using namespace std;

int main()
{
    // Set final time and time step.
    const double tf = 1e3;
    const double h = 1e-4;
    const double L_update_period = 30;
    const double lambda = 1.1;
    const int m = static_cast<int>(ceil(L_update_period/h));

    // Feedback parameters.
    const double k0 = 50, k1 = 100, k2 = 50;
    const double L0 = 1986;
    double alpha = 1/(h*L0);
    // double alpha = 1;
    
    // Define initial rotation matrix R0 as the identity matrix.
    const double R0[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    // Define the initial angular velocity Omega0.
    const double Omega0[3] = {1.0, 1.0, 1.0};

    // Define a sample inertia matrix I (here, chosen diagonal for simplicity).
    const double I[3][3] = {
        {3.0, 0.0, 0.0},
        {0.0, 2.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    double E0 = E(Omega0, I);
    double Pi0[3];
    Pi(R0, Omega0, I, Pi0);

    // Pointers for the outputs.
    double *R_out_vanilla = nullptr;    // Will store 9*N doubles (3x3 matrices for each step).
    double *Omega_out_vanilla = nullptr; // Will store 3*N doubles.
    double *t_out_vanilla = nullptr;     // Will store N time values.
    int N_vanilla = 0;                   // Number of integration steps.

    double *R_out_feedback_vanilla = nullptr;    // Will store 9*N doubles (3x3 matrices for each step).
    double *Omega_out_feedback_vanilla = nullptr; // Will store 3*N doubles.
    double *t_out_feedback_vanilla = nullptr;     // Will store N time values.
    int N_feedback_vanilla = 0;                   // Number of integration steps.

    double *R_out_feedback = nullptr;    // Will store 9*N doubles (3x3 matrices for each step).
    double *Omega_out_feedback = nullptr; // Will store 3*N doubles.
    double *t_out_feedback = nullptr;     // Will store N time values.
    int N_feedback = 0;                   // Number of integration steps.

    double *R_out_adaptive = nullptr;    // Will store 9*N doubles (3x3 matrices for each step).
    double *Omega_out_adaptive = nullptr; // Will store 3*N doubles.
    double *t_out_adaptive = nullptr;     // Will store N time values.
    int N_adaptive = 0;                   // Number of integration steps.

    double error_vanilla, error_feedback, error_feedback_vanilla, error_adaptive;

    // Run the Euler integration.
    chrono::system_clock::time_point start = chrono::system_clock::now();
    // euler_vanilla(R0, Omega0, I, tf, h, R_out_vanilla, Omega_out_vanilla, t_out_vanilla, N_vanilla);
    error_vanilla = euler_vanilla(R0, Omega0, I, E0, Pi0, k0, k1, k2, tf, h);
    chrono::duration<double> duration = chrono::system_clock::now() - start;
    cout << "CPU time for Vanilla Euler method : " << duration.count() << endl;
    cout << "Max error for Vanilla Euler method : " << error_vanilla << endl;


    // Run Euler integration with feedback.
    start = chrono::system_clock::now();
    error_feedback_vanilla = euler_feedback(R0, Omega0, I,
                                    k0, k1, k2,
                                    E0, Pi0,
                                    1.0, tf, h);
    duration = chrono::system_clock::now() - start;
    cout << "CPU time for Feedback Euler method (Vanilla) : " << duration.count() << endl;
    cout << "Max error for Feedback Euler method (Vanilla) : " << error_feedback_vanilla << endl;


    // Run Euler integration with feedback.
    start = chrono::system_clock::now();
    error_feedback = euler_feedback(R0, Omega0, I,
                                    k0, k1, k2,
                                    E0, Pi0,
                                    1.0, tf, h);
    duration = chrono::system_clock::now() - start;
    cout << "CPU time for Feedback Euler method : " << duration.count() << endl;
    cout << "Max error for Feedback Euler method : " << error_feedback << endl;

    
    // Run Euler integration with adaptive feedback.
    start = chrono::system_clock::now();
    error_adaptive = euler_adaptive(R0, Omega0, I,
                                    k0, k1, k2,
                                    E0, Pi0,
                                    tf, h,
                                    m, lambda);
    duration = chrono::system_clock::now() - start;
    cout << "CPU time for Adaptive Feedback Euler method : " << duration.count() << endl;
    cout << "Max error for Adaptive Euler method : " << error_adaptive << endl;

    // Print final state for verification.
    // std::cout << "Final time: " << t_out_vanilla[N_vanilla-1] << std::endl;
    // std::cout << "Final Rotation Matrix (R) at step " << N_vanilla << ":\n";
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++)
    //         std::cout << R_out_vanilla[(N_vanilla-1)*9 + i*3 + j] << " ";
    //     std::cout << std::endl;
    // }
    // std::cout << "Final Omega: ";
    // for (int i = 0; i < 3; i++)
    //     std::cout << Omega_out_vanilla[(N_vanilla-1)*3 + i] << " ";
    // std::cout << std::endl;

    // std::cout << "Final time: " << t_out_feedback[N_feedback-1] << std::endl;
    // std::cout << "Final Rotation Matrix (R) at step " << N_feedback << ":\n";
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++)
    //         std::cout << R_out_feedback[(N_feedback-1)*9 + i*3 + j] << " ";
    //     std::cout << std::endl;
    // }
    // std::cout << "Final Omega: ";
    // for (int i = 0; i < 3; i++)
    //     std::cout << Omega_out_feedback[(N_feedback-1)*3 + i] << " ";
    // std::cout << std::endl;

    // std::cout << "Final time: " << t_out_adaptive[N_adaptive-1] << std::endl;
    // std::cout << "Final Rotation Matrix (R) at step " << N_adaptive << ":\n";
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++)
    //         std::cout << R_out_adaptive[(N_adaptive-1)*9 + i*3 + j] << " ";
    //     std::cout << std::endl;
    // }
    // std::cout << "Final Omega: ";
    // for (int i = 0; i < 3; i++)
    //     std::cout << Omega_out_adaptive[(N_adaptive-1)*3 + i] << " ";
    // std::cout << std::endl;

    // Clean up allocated memory.
    delete[] R_out_vanilla;
    delete[] Omega_out_vanilla;
    delete[] t_out_vanilla;
    delete[] R_out_feedback;
    delete[] Omega_out_feedback;
    delete[] t_out_feedback;
    delete[] R_out_adaptive;
    delete[] Omega_out_adaptive;
    delete[] t_out_adaptive;

    return 0;
}