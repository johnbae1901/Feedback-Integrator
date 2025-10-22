#pragma once
// include/splitting_strang.h
// Strang (symmetric) three-rotations splitting for the free rigid body:
//   R3(h/2) -> R2(h/2) -> R1(h) -> R2(h/2) -> R3(h/2)
// Matches the project's "single-scalar error per run" API by returning
// max_t getError(state(t), ...).  'state' follows the 3x4 layout used in Euler.

std::tuple<
    double, 
    double, 
    double, 
    double> 
    splitting_strang(
    const double R0[3][3],
    const double Omega0[3],
    const double I[3][3],
    const double E0,
    const double Pi0[3],
    const double k0,
    const double k1,
    const double k2,
    double tf,
    double h
);
