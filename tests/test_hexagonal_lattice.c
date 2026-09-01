#include "EDMD.h"

#include <math.h>
#include <stdio.h>

extern int N, Nbig, Nsmall, Nx, Ny;
#if THREE_D
extern int Nz;
#endif
extern int load, S1, Hex, hexRadius;
extern int controlLenght;
extern particle *particles;
extern double phi, aspectRatio, Einit, Lx, Ly;
#if THREE_D
extern double Lz;
#endif

static int failures;

static void expect_true(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        failures++;
    }
}

static void configure_lattice(void)
{
    N = 1000;
    Nbig = Nsmall = 0;
    load = 0;
    S1 = 0;
    Hex = 1;
    hexRadius = 0;
    controlLenght = 0;
    phi = 0.50;
    aspectRatio = 1.0;
    Einit = 1.0;
    constantInit(0, NULL);
    particlesInit();
}

static void test_dimensions(void)
{
    int nx = 0, ny = 0, nz = 0;
    expect_true(hexagonalLatticeDimensions(1000, 1.0, &nx, &ny, &nz),
                "could not derive ordinary hexagonal-lattice dimensions");
    expect_true(nx > 0 && ny > 0, "hexagonal lattice has a zero dimension");
    expect_true((ny & 1) == 0,
                "triangular rows do not close periodically");
#if THREE_D
    expect_true(nz > 0 && (nz & 1) == 0,
                "HCP AB layers do not close periodically");
    expect_true(!hexagonalLatticeDimensions(3, 1.0, &nx, &ny, &nz),
                "undersized HCP configuration was accepted");
#else
    expect_true(nz == 1, "2D lattice unexpectedly created z layers");
#endif
    expect_true(!hexagonalLatticeDimensions(1000, -1.0, &nx, &ny, &nz),
                "negative lattice aspect ratio was accepted");
}

static void test_realized_state(void)
{
    configure_lattice();
#if THREE_D
    expect_true(N == Nx*Ny*Nz, "realized HCP count is inconsistent");
    expect_true(Lz > 0.0, "HCP box has no z extent");
    expect_true(fabs(Ly/Lx - 1.0) < 0.20 &&
                fabs(Lz/Lx - 1.0) < 0.20,
                "HCP dimension search produced a badly distorted box");
#else
    expect_true(N == Nx*Ny, "realized triangular-lattice count is inconsistent");
#endif
    expect_true(Nbig == N && Nsmall == 0,
                "hexagonal lattice species counts are inconsistent");

    double minimum_distance_squared = INFINITY;
#if THREE_D
    double minimum_z = INFINITY;
    double maximum_z = -INFINITY;
#endif
    for (int i = 0; i < N; ++i) {
        const particle *p = particles + i;
        expect_true(p->num == i, "hexagonal lattice particle index is invalid");
        expect_true(isfinite(p->x) && isfinite(p->y) &&
                    isfinite(p->vx) && isfinite(p->vy),
                    "hexagonal lattice contains a non-finite planar state");
        expect_true(p->x >= 0.0 && p->x < Lx &&
                    p->y >= 0.0 && p->y < Ly,
                    "hexagonal lattice site lies outside its box");
#if THREE_D
        expect_true(isfinite(p->z) && isfinite(p->vz),
                    "HCP lattice did not initialize z/vz");
        expect_true(p->z >= 0.0 && p->z < Lz,
                    "HCP lattice site lies outside z bounds");
        minimum_z = fmin(minimum_z, p->z);
        maximum_z = fmax(maximum_z, p->z);
#endif
        for (int j = i + 1; j < N; ++j) {
            double dx = particles[j].x - p->x;
            double dy = particles[j].y - p->y;
#if THREE_D
            double dz = particles[j].z - p->z;
            PBC(&dx, &dy, &dz);
            const double distance_squared = dx*dx + dy*dy + dz*dz;
#else
            PBC(&dx, &dy);
            const double distance_squared = dx*dx + dy*dy;
#endif
            minimum_distance_squared = fmin(minimum_distance_squared,
                                             distance_squared);
        }
    }
#if THREE_D
    expect_true(maximum_z > minimum_z,
                "HCP initializer collapsed every particle onto one z layer");
#endif
    expect_true(minimum_distance_squared >= 4.0*(1.0 - 1e-11),
                "hexagonal lattice starts with overlapping unit particles");
    freeArrays();
}

int main(void)
{
    test_dimensions();
    test_realized_state();
    if (failures != 0) {
        fprintf(stderr, "%d hexagonal-lattice regression(s) failed\n", failures);
        return 1;
    }
    printf("hexagonal-lattice regressions passed (%dD)\n", THREE_D ? 3 : 2);
    return 0;
}
