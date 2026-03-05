#ifndef INTRUDER_H
#define INTRUDER_H
#include "EDMD.h"

#define vn_tol 1e-8
#define INTRUDER_TOL 1e-10

extern double t;
typedef struct polygon polygon;
struct polygon{
    double vertices[20][2];
    double lengths[20];
    double normals[20][2];
    double tangents[20][2];
    int n_vertices;
    double max_extent;
};

typedef struct intruder intruder;
struct intruder{
    double t;
    double x, y, vx, vy, angle, omega;
    double M, Im;
    polygon body_frame_shape;
    polygon real_shape;
    long unsigned int coll;
};

enum intruderShape{
    SQUARE,
    TRIANGLE,
    WHEEL,
    CHIRAL_WHEEL,
    GEN_WHEEL
};

typedef struct collisionInfo collisionInfo;
struct collisionInfo{
    double time;
    double normal[2];
};

typedef struct {
    int enabled;
    double burn_time;
    int has_prev;
    double last_t;
    double phi0, vx0, vy0, om0;   // start-of-flight state (lab) after previous intruder collision
    unsigned long long nsamples;

    // time accumulators (measurement window)
    double sum_dt;                // total physical time used
    double sum_Etrans_dt;         // ∫ (0.5*M*(Vx^2+Vy^2)) dt
    double sum_dU[3];             // sum of body inertial increments: [ΔVb_x, ΔVb_y, ΔΩ]

    // exact free-flight integrals for time averages
    double sum_intVb[2];          // ∫ V_body dt
    double sum_intOm_dt;          // ∫ Ω dt
    double sum_intOmVb[2];        // ∫ Ω V_body dt  (Ω constant per flight)

    // Weighted normal equations for L fit (weights = 1/dt): (X/√dt)^T (X/√dt), etc.
    double XtWX[4][4];
    double XtWY[4][3];

    // Unweighted sums for diffusion (computed from residual increments at the end)
    double XTX[4][4];
    double XTY[4][3];
    double YTY[3][3];
} IntruderKM;

void intruderInit(intruder* intr, double x, double y, double vx, double vy, double angle, double omega, double M, double Im, int shape, double* info, int km_enabled, IntruderKM* km);
polygon square(double side_length);
polygon triangle(double base_length, double height);
void rotateTranslatePolygon(intruder* intr, double angle, double X, double Y);
void printIntruderShapeFromSmallParticles(FILE* fichier, intruder* intr, double rad, int N);
int countFakeParticlesForIntruder(intruder* intr, double rad);
double collisionTimeIntruder(intruder* intr, particle* p);

// Helper functions for collisionTimeIntruder
double signedDistanceParticlePolygon(particle* p, intruder* intru, double t);
void computeCollisionNormal(particle* p, intruder* intru, double t, double* nx, double* ny);
double pointSegmentDistance(double px, double py, double ax, double ay, double bx, double by);


void freeFlyIntruder(intruder* intru);
int particleOutsideIntruder(particle* p, intruder* intru);




void km_finalize_and_dump(IntruderKM km, const char* intruderKMName, intruder intru);
void km_integrate_vbody(double vx, double vy, double phi0, double om, double dt,
                                      double* int_vbx, double* int_vby);

int solve3x3(const double A_in[3][3], const double b_in[3], double x_out[3]);
int solve4x4_3rhs(const double A_in[4][4], const double B_in[4][3], double X_out[4][3]);

#endif