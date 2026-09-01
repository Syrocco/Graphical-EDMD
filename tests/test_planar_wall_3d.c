#include "EDMD.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !THREE_D
#error "test_planar_wall_3d requires THREE_D=1"
#endif

extern int N, S1, Hex;
extern int periodicDumps, ther, Nxcells, Nycells, Nzcells;
extern int addWallz;
extern particle *particles;
extern particle **cellList;
extern node **eventList;
extern node *nextEvent;
extern double t, Lx, Ly, Lz, resW, resT;
extern double dp, cumulativeWallImpulse;
#if STRESS
extern double sumKineticZ;
#endif
#if HEAT
extern double sumHeatZ;
#endif

static int failures;

static void expect_true(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        failures++;
    }
}

static void expect_close(double actual, double expected, const char *message)
{
    const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
    if (!isfinite(actual) || fabs(actual - expected) > 1e-12*scale) {
        fprintf(stderr, "%s: got %.17g, expected %.17g\n",
                message, actual, expected);
        failures++;
    }
}

int main(void)
{
    addWallz = 1;
    N = 1;
    S1 = Hex = 0;
    periodicDumps = ther = 0;
    t = 0.0;
    Lx = Ly = Lz = 20.0;
    resW = resT = 1.0;
    dp = cumulativeWallImpulse = 0.0;
    particles = calloc(1, sizeof(*particles));
    expect_true(particles != NULL, "unable to allocate z-wall test particle");
    if (particles == NULL) return 1;
    particles[0].num = 0;
    particles[0].m = 1.0;
    particles[0].rad = 1.0;
    particles[0].x = 10.0;
    particles[0].y = 10.0;
    particles[0].z = 18.0;
    particles[0].vz = 2.0;
    boxConstantHelper();
    cellListInit();
#if (STRESS || HEAT)
    initializeTransportAccumulators();
#endif
    eventListInit();

    node *wall = eventList[N];
    expect_true(wall->type == WALL && wall->j == 3,
                "upper z-wall event was not selected");
    expect_close(wall->t, 0.5, "upper z-wall collision time");
    nextEvent = wall;
    t = wall->t;
    removeEventFromQueue(wall);
    doTheWallNormal();
    expect_close(particles[0].z, 19.0, "z-wall contact projection");
    expect_close(particles[0].vz, -2.0, "z-wall reflection");
    expect_close(dp, 4.0, "z-wall interval impulse");
    expect_close(cumulativeWallImpulse, 4.0,
                 "z-wall cumulative impulse");
#if STRESS
    expect_true(isfinite(sumKineticZ),
                "z-wall produced non-finite stress transport");
#endif
#if HEAT
    expect_true(isfinite(sumHeatZ),
                "z-wall produced non-finite heat transport");
#endif

    double dx = 0.0, dy = 0.0, dz = 18.0;
    PBC(&dx, &dy, &dz);
    expect_close(dz, 18.0, "z wall incorrectly retained periodic z imaging");
    expect_true(particles[0].cell[2] == Nzcells - 1,
                "z-wall particle has an invalid cell membership");
    expect_true(cellList[(Nzcells - 1)*Nycells*Nxcells +
                         particles[0].cell[1]*Nxcells +
                         particles[0].cell[0]] == particles,
                "z-wall particle disappeared from the 3D cell list");

    freeArrays();
    if (failures != 0) {
        fprintf(stderr, "%d full-box regression(s) failed\n", failures);
        return 1;
    }
    printf("full 3D box regressions passed\n");
    return 0;
}
