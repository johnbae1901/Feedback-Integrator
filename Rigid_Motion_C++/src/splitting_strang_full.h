#pragma once
// History-output version of Strang splitting.
// Produces arrays (caller must delete[]).

void splitting_strang_full(
    const double R0[3][3],
    const double Omega0[3],
    const double I[3][3],
    double tf,
    double h,
    double*& R_hist,    // 9*(N) row-major; we follow Euler code that stored N steps (no +1)
    double*& Omega_hist,// 3*(N)
    double*& t_hist,    // N
    int& N_out
);
