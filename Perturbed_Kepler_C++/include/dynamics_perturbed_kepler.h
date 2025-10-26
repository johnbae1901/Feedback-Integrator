#pragma once
#include <cmath>
#include "basic_operations.h"  // norm3, cross3, etc. (uses Eigen under the hood)

// ---------- Types ----------
struct PKParams {            // physical params
    double mu   = 1.0;       // gravitational parameter
    double delta= 0.0025;    // perturbation strength in U = -mu/r - delta/r^3
};

struct FIRHSParams {         // FI gains
    double kE = 2.0;         // energy penalty
    double kL = 3.0;         // angular-momentum penalty
};

// ---------- small helpers ----------
inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline void axpy3(double y[3], const double a, const double x[3]) {
    y[0] += a*x[0]; y[1] += a*x[1]; y[2] += a*x[2];
}
inline void scal3(double a, double x[3]) {
    x[0]*=a; x[1]*=a; x[2]*=a;
}
inline void copy3(const double src[3], double dst[3]) {
    dst[0]=src[0]; dst[1]=src[1]; dst[2]=src[2];
}
inline void sub3(const double a[3], const double b[3], double out[3]) {
    out[0]=a[0]-b[0]; out[1]=a[1]-b[1]; out[2]=a[2]-b[2];
}

// ---------- Potential U(r) and derivatives ----------
inline double U_scalar(double r, const PKParams& P) {
    // U(r) = -mu/r - delta/r^3
    return -P.mu/r - P.delta/(r*r*r);
}
inline double Ur_scalar(double r, const PKParams& P) {
    // U'(r) = mu/r^2 + 3*delta/r^4
    return  P.mu/(r*r) + 3.0*P.delta/(r*r*r*r);
}

inline double energy_pk(const double x[3], const double v[3], const PKParams& P){
    return 0.5*dot3(v,v) + U_scalar(norm3(x), P);
}

// a(x) = vdot = -U'(r) * x / r = -(mu/r^3 + 3*delta/r^5) * x
inline void accel(const double x[3], const PKParams& P, double a[3]) {
    const double r = norm3(x);
    const double coef = - Ur_scalar(r, P) / r;
    a[0] = coef * x[0];
    a[1] = coef * x[1];
    a[2] = coef * x[2];
}

// invariants
inline double energy(const double x[3], const double v[3], const PKParams& P){
    return 0.5*dot3(v,v) + U_scalar(norm3(x), P);
}
inline void angmom(const double x[3], const double v[3], double L[3]){
    cross3(x, v, L);
}

// ---------- Feedback Integrator RHS (paper eq. (46) form) ----------
// V = (kE/2)(E-E0)^2 + (kL/2)||L-L0||^2
// xdot = v - kE*dE*Ur(r)*x/r - kL*(v × dL)
// vdot = a(x) - kE*dE*v       - kL*(dL × x)
inline void fi_rhs(const double x[3], const double v[3],
                   const PKParams& P, const FIRHSParams& K,
                   const double E0, const double L0[3],
                   double xdot[3], double vdot[3]) {
    const double r  = norm3(x);
    const double upr= Ur_scalar(r, P);

    // dL = L - L0
    double L[3]; angmom(x, v, L);
    double dL[3]; sub3(L, L0, dL);

    // dE = E - E0
    const double dE = energy(x, v, P) - E0;

    // xdot = v - kE*dE*(upr/r)*x - (v × dL)
    // tmp1 = (upr/r)*x
    double tmp1[3] = { x[0], x[1], x[2] };
    scal3(upr / r, tmp1);

    // v × dL
    double vXdL[3]; cross3(v, dL, vXdL);

    xdot[0] = v[0] - K.kE*dE*tmp1[0] - vXdL[0];
    xdot[1] = v[1] - K.kE*dE*tmp1[1] - vXdL[1];
    xdot[2] = v[2] - K.kE*dE*tmp1[2] - vXdL[2];

    // vdot = a(x) - kE*dE*v - (dL × x)
    double ax[3]; accel(x, P, ax);
    double dLxX[3]; cross3(dL, x, dLxX);

    vdot[0] = ax[0] - K.kE*dE*v[0] - dLxX[0];
    vdot[1] = ax[1] - K.kE*dE*v[1] - dLxX[1];
    vdot[2] = ax[2] - K.kE*dE*v[2] - dLxX[2];
}

// ---------- One-step integrators (same signatures you already use) ----------
inline void step_sv(double x[3], double v[3], double h, const PKParams& P) {
    // velocity-Verlet
    double a0[3]; accel(x, P, a0);
    // v += 0.5*h*a0
    axpy3(v, 0.5*h, a0);
    // x += h*v
    axpy3(x, h, v);
    // a1 at new x
    double a1[3]; accel(x, P, a1);
    // v += 0.5*h*a1
    axpy3(v, 0.5*h, a1);
}

inline void step_fi_euler(double x[3], double v[3], double h,
                          const PKParams& P, const FIRHSParams& K,
                          const double E0, const double L0[3]) {
    double xdot[3], vdot[3];
    fi_rhs(x, v, P, K, E0, L0, xdot, vdot);
    axpy3(x, h, xdot);
    axpy3(v, h, vdot);
}
