// Spherical Pendulum Bench: Euler vs Feedback Integrator vs Stormer–Verlet (RATTLE-like)
// Self-contained C++17. No dependencies beyond <cmath> / <vector> / <iostream> / <fstream>
// Model: x in R^3 with |x|=1, dynamics: x' = v, v' = -lambda x - grad U(x)
// Potential: U(x) = g * e3^T x + eps * x1 * x2
// lambda = v·v - x·gradU (from enforcing x·x=1 and x·v=0)
// Feedback Integrator uses V(x,v) = k0/2 (|x|^2-1)^2 + k1/2 (E-E0)^2 + k2/2 (Lz-Lz0)^2
// and modifies dynamics by subtracting grad_x V, grad_v V from RHS.

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>

using namespace std;

// define PI portable
static constexpr double PI = 3.141592653589793238462643383279502884L;

struct Vec3 {
    double x,y,z;
    Vec3(double x_=0,double y_=0,double z_=0):x(x_),y(y_),z(z_){}
    Vec3 operator+(const Vec3& o) const { return {x+o.x,y+o.y,z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x,y-o.y,z-o.z}; }
    Vec3 operator*(double s) const { return {x*s,y*s,z*s}; }
    Vec3 operator/(double s) const { return {x/s,y/s,z/s}; }
    Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
    Vec3& operator-=(const Vec3& o){ x-=o.x; y-=o.y; z-=o.z; return *this; }
};
static inline double dot(const Vec3& a,const Vec3& b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 cross(const Vec3& a,const Vec3& b){ return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x }; }
static inline double norm(const Vec3& a){ return sqrt(dot(a,a)); }
static inline Vec3 normalize(const Vec3& a){ double n = norm(a); return (n>0)? a/n : a; }

struct State { Vec3 x; Vec3 v; };

struct Params {
    double g = 1.0;       // gravity strength
    double eps = 0.0;    // weak asymmetric perturbation
};

struct FIParams {
    double k0 = 50.0;  // constraint penalty
    double k1 = 20.0;   // energy penalty
    double k2 = 2.0;   // Lz penalty
};

struct Mat3 { double a[3][3]{{0,0,0},{0,0,0},{0,0,0}}; };

static inline Mat3 I3(double s=1.0){
    Mat3 M; M.a[0][0]=M.a[1][1]=M.a[2][2]=s; return M;
}
static inline Mat3 outer(const Vec3& p, const Vec3& q){
    Mat3 M; 
    M.a[0][0]=p.x*q.x; M.a[0][1]=p.x*q.y; M.a[0][2]=p.x*q.z;
    M.a[1][0]=p.y*q.x; M.a[1][1]=p.y*q.y; M.a[1][2]=p.y*q.z;
    M.a[2][0]=p.z*q.x; M.a[2][1]=p.z*q.y; M.a[2][2]=p.z*q.z;
    return M;
}
static inline Mat3 add(const Mat3& A, const Mat3& B){
    Mat3 C; 
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) C.a[i][j]=A.a[i][j]+B.a[i][j];
    return C;
}
static inline Mat3 scal(const Mat3& A, double s){
    Mat3 C; 
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) C.a[i][j]=s*A.a[i][j];
    return C;
}
static inline Mat3 mul(const Mat3& A, const Mat3& B){
    Mat3 C;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++){
        C.a[i][j]=0.0;
        for(int k=0;k<3;k++) C.a[i][j]+=A.a[i][k]*B.a[k][j];
    }
    return C;
}
static inline Vec3 mul(const Mat3& A, const Vec3& x){
    return {
        A.a[0][0]*x.x + A.a[0][1]*x.y + A.a[0][2]*x.z,
        A.a[1][0]*x.x + A.a[1][1]*x.y + A.a[1][2]*x.z,
        A.a[2][0]*x.x + A.a[2][1]*x.y + A.a[2][2]*x.z
    };
}

// Potential and its gradient
static inline double U(const Vec3& x, const Params& P){
    return P.g * x.z + P.eps * x.x * x.y;
}
static inline Vec3 gradU(const Vec3& x, const Params& P){
    return { P.eps * x.y, P.eps * x.x, P.g };
}

// Hessian of U(x) = g z + eps x y
static inline Mat3 hessU(const Params& P){
    Mat3 H;
    H.a[0][1]=H.a[1][0]=P.eps; // d2U/dxdy = d2U/dydx = eps
    // others are zero
    return H;
}

// Build H = ∇²_x V(x,v)
static inline Mat3 hessV_x(const State& s, const Params& P, const FIParams& K, double E0, double Lz0){
    double r2 = dot(s.x,s.x);
    Mat3 H0 = add( scal(I3(1.0), 2.0*K.k0*(r2-1.0) ), scal( outer(s.x,s.x), 4.0*K.k0 ) );

    double E  = 0.5*dot(s.v,s.v) + U(s.x,P);
    double dE = E - E0;
    Vec3 gU = gradU(s.x,P);
    Mat3 H1 = add( scal( outer(gU,gU), K.k1 ), scal( hessU(P), K.k1*dE ) );

    Vec3 dLdx = { s.v.y, -s.v.x, 0.0 };
    Mat3 H2 = scal( outer(dLdx,dLdx), K.k2 );

    return add( add(H0,H1), H2 );
}

// spectral norm ||H||_2 via power iteration on H^T H (here H is symmetric, so H^T H = H^2)
static inline double spectral_norm_2(const Mat3& H){
    Mat3 S = mul(H,H); // symmetric PSD
    Vec3 y{1.0,0.3,0.2};
    double lam = 0.0;
    for(int it=0; it<20; ++it){
        Vec3 z = mul(S,y);
        double nz = norm(z);
        if(nz < 1e-18) { lam = 0.0; break; }
        y = z / nz;
        lam = dot(y, mul(S,y)); // Rayleigh quotient
    }
    if(lam < 0.0) lam = 0.0;
    return sqrt(lam);
}

// Invariants
static inline double energy(const State& s, const Params& P){ return 0.5*dot(s.v,s.v) + U(s.x,P); }
static inline double Lz(const State& s){ return (s.x.x*s.v.y - s.x.y*s.v.x); }

// Lagrange multiplier for exact constrained dynamics (continuous-time formula)
static inline double lambda_ct(const State& s, const Params& P){
    Vec3 gU = gradU(s.x,P);
    return dot(s.v,s.v) - dot(s.x,gU);
}

// ========== Euler (explicit) ==========
void step_euler(State& s, double h, const Params& P){
    double lam = lambda_ct(s,P);
    Vec3 ax = s.v;                   // x' = v
    Vec3 av = (s.x * (-lam)) - gradU(s.x,P); // v' = -lambda x - gradU
    s.x += ax * h;
    s.v += av * h;
}

// Optional: project back to manifold and tangent (used by some variants)
void project_to_sphere_and_tangent(State& s){
    s.x = normalize(s.x);
    // remove radial component of velocity to enforce x·v = 0
    double radial = dot(s.x, s.v);
    s.v -= s.x * radial;
}

// ========== Feedback Integrator ==========
struct FIState { State s; double E0; double Lz0; };

static inline void fi_gradients(const FIState& fis, const Params& P, const FIParams& K, Vec3& dVdx, Vec3& dVdv){
    const State& s = fis.s;
    // V0: constraint penalty
    double r2 = dot(s.x,s.x);
    Vec3 dV0dx = s.x * (2.0 * K.k0 * (r2 - 1.0)); // 2 k0 (|x|^2-1) x
    // V1: energy deviation
    double E = energy(s,P); double dE = E - fis.E0;
    Vec3 dV1dx = gradU(s.x,P) * (K.k1 * dE);
    Vec3 dV1dv = s.v * (K.k1 * dE);
    // V2: Lz deviation
    double lz = Lz(s); double dlz = lz - fis.Lz0;
    Vec3 dLzdx = { s.v.y, -s.v.x, 0.0 };
    Vec3 dLzdv = { -s.x.y, s.x.x, 0.0 };
    Vec3 dV2dx = dLzdx * (K.k2 * dlz);
    Vec3 dV2dv = dLzdv * (K.k2 * dlz);
    dVdx = dV0dx + dV1dx + dV2dx;
    dVdv = dV1dv + dV2dv;
}

void step_fi(State& s, double h, const Params& P, const FIParams& K, const double E0, const double Lz0){
    FIState fis{ s, E0, Lz0 };
    Vec3 dVdx, dVdv; 
    fi_gradients(fis,P,K,dVdx,dVdv);

    // compute scaling alpha = 1 / (h * || ∇_x^2 V(x) ||_2)
    Mat3 H = hessV_x(s, P, K, E0, Lz0);
    double Hn = spectral_norm_2(H);
    double alpha = (Hn > 1e-12) ? (1.0 / (h * Hn)) : (1.0 / h);

    double lam = lambda_ct(s,P);
    Vec3 ax = s.v - dVdx; // modified x'
    Vec3 av = (s.x * (-lam)) - gradU(s.x,P) - (dVdv * alpha); // scaled feedback on v'

    s.x += ax * h;
    s.v += av * h;
}


// ========== Stormer–Verlet (RATTLE-like pragmatic constraint handling) ==========
// This is a practical SV for constrained sphere: half-kick with current lambda, drift, renormalize, second half-kick, re-tangent velocity.
void step_sv_rattle_like(State& s, double h, const Params& P){
    // half kick
    double lam_n = lambda_ct(s,P);
    Vec3 a_half = gradU(s.x,P) + s.x * lam_n; // note: dynamics v' = -gradU - lambda x -> kick uses minus
    Vec3 v_half = s.v - a_half * (0.5 * h);
    // drift
    Vec3 x_new = s.x + v_half * h;
    x_new = normalize(x_new); // position constraint
    // second half kick
    State tmp{ x_new, v_half };
    double lam_np1 = lambda_ct(tmp,P); // use updated x, v_half to estimate lambda
    Vec3 a2 = gradU(x_new,P) + x_new * lam_np1;
    Vec3 v_new = v_half - a2 * (0.5 * h);
    // enforce tangency: remove any radial component
    double vrad = dot(x_new, v_new);
    v_new -= x_new * vrad;
    s.x = x_new; s.v = v_new;
}

struct Metrics { double max_constr=0, max_dE=0, max_dLz=0; };

// Simulate and write CSV
Metrics run_sim(const string& csv_path, State s0, double h, int steps, const Params& P, 
                int method, const FIParams& K = FIParams()){
    ofstream f(csv_path);
    f << "t,x,y,z,vx,vy,vz,|x|-1,E-E0,Lz-Lz0\n";
    State s = s0;
    double E0 = energy(s0,P); double Lz0 = Lz(s0);
    Metrics M;
    for(int n=0;n<=steps;n++){
        double t = n*h;
        double constr = fabs(norm(s.x)-1.0);
        double dE = energy(s,P) - E0;
        double dL = Lz(s) - Lz0;
        M.max_constr = max(M.max_constr, constr);
        M.max_dE = max(M.max_dE, fabs(dE));
        M.max_dLz = max(M.max_dLz, fabs(dL));
        f << fixed << setprecision(10)
          << t << "," << s.x.x << "," << s.x.y << "," << s.x.z << ","
          << s.v.x << "," << s.v.y << "," << s.v.z << ","
          << constr << "," << dE << "," << dL << "\n";
        if(n==steps) break;
        if(method==0){
            step_euler(s,h,P);
        } else if(method==1){
            step_fi(s,h,P,K,E0,Lz0);
        } else if(method==2){
            step_sv_rattle_like(s,h,P);
        }
    }
    f.close();
    return M;
}

int main(){
    // Initial condition: tilt and spin (nutation + precession region)
    Params P; P.g = 1.0; P.eps = 0.0;
    double theta0 = PI/3.0; // 60 degrees
    State s0;
    s0.x = { sin(theta0), 0.0, cos(theta0) };
    s0.v = { 0.0, 1.2, 0.0 }; // tangential initial velocity
    // Ensure constraints (in case of numeric edits)
    s0.x = normalize(s0.x);
    s0.v -= s0.x * dot(s0.x, s0.v);

    double h = 1e-4; // shared step size
    double T = 500.0; // total time
    int steps = (int)llround(T / h);

    // Feedback gains
    FIParams K; K.k0 = 80.0; K.k1 = 3.0; K.k2 = 3.0;

    cerr << "Running with h=" << h << ", steps=" << steps << "\n";

    auto M_eu = run_sim("traj_euler.csv", s0, h, steps, P, /*method=*/0);
    auto M_fi = run_sim("traj_fi.csv",    s0, h, steps, P, /*method=*/1, K);
    auto M_sv = run_sim("traj_sv.csv",    s0, h, steps, P, /*method=*/2);

    cout.setf(std::ios::fixed); cout<<setprecision(6);
    cout << "== Max errors (over [0,"<<T<<"]) ==\n";
    cout << "Euler:       max| |x|-1 |=" << M_eu.max_constr 
         << ", max|E-E0|=" << M_eu.max_dE 
         << ", max|Lz-Lz0|=" << M_eu.max_dLz << "\n";
    cout << "FeedbackInt: max| |x|-1 |=" << M_fi.max_constr 
         << ", max|E-E0|=" << M_fi.max_dE 
         << ", max|Lz-Lz0|=" << M_fi.max_dLz << "\n";
    cout << "SV (RATTLE): max| |x|-1 |=" << M_sv.max_constr 
         << ", max|E-E0|=" << M_sv.max_dE 
         << ", max|Lz-Lz0|=" << M_sv.max_dLz << "\n";

    cerr << "CSV written: traj_euler.csv, traj_fi.csv, traj_sv.csv\n";
    cerr << "Columns: t,x,y,z,vx,vy,vz,|x|-1,E-E0,Lz-Lz0\n";
    return 0;
}
