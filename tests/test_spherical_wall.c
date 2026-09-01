#include "EDMD.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef EDMD_ADD_FIELD
#define EDMD_ADD_FIELD 0
#endif

extern int N, load, S1, Hex;
extern int periodicDumps, ther;
extern int Nxcells, Nycells;
#if THREE_D
extern int Nzcells;
#endif
extern unsigned long int ncol;
extern particle *particles;
extern particle **cellList;
extern node *nextEvent;
extern double t, Lx, Ly, halfLx, halfLy;
#if THREE_D
extern double Lz, halfLz;
#endif
extern double resW, resT;
extern double field;
extern double dp, cumulativeWallImpulse;

#if (STRESS || HEAT)
extern void initializeTransportAccumulators(void);
#endif

static int failures;

static void fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    failures++;
}

static void expect_true(int condition, const char *message)
{
    if (!condition) fail(message);
}

static void expect_close(double actual, double expected, const char *message)
{
    const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
    if (!isfinite(actual) || fabs(actual-expected) > 2e-12*scale) {
        fprintf(stderr, "%s: got %.17g, expected %.17g\n",
                message, actual, expected);
        failures++;
    }
}

static particle *configure_single_particle(void)
{
    N = 1;
    load = 2;
    S1 = 0;
    Hex = 0;
    periodicDumps = 0;
    ther = 0;
    t = 0.0;
    ncol = 0;
    dp = 0.0;
    cumulativeWallImpulse = 0.0;
    Lx = Ly = 20.0;
#if THREE_D
    Lz = 20.0;
#endif
    particles = calloc(1, sizeof(*particles));
    if (particles == NULL) {
        fail("unable to allocate spherical-wall test particle");
        return NULL;
    }
    particle *value = particles;
    value->num = 0;
    value->m = 1.0;
    value->rad = 1.0;
    value->t = 0.0;
    boxConstantHelper();
    value->x = halfLx;
    value->y = halfLy;
#if THREE_D
    value->z = halfLz;
#endif
    return value;
}

static void start_single_particle_engine(void)
{
    cellListInit();
    eventListInit();
#if (STRESS || HEAT)
    initializeTransportAccumulators();
#endif
}

static void process_immediate_wall(void)
{
    nextEvent = findNextEvent();
    expect_true(nextEvent != NULL, "spherical wall event was not scheduled");
    if (nextEvent == NULL) return;
    expect_true(nextEvent->type == WALL,
                "first single-particle event was not the spherical wall");
    expect_close(nextEvent->t, t, "immediate spherical wall event time");
    t = nextEvent->t;
    removeEventFromQueue(nextEvent);
    doTheWallNormal();
}

#if !EDMD_ADD_FIELD
static void test_fully_dissipative_contact_stops_cleanly(void)
{
    particle *value = configure_single_particle();
    if (value == NULL) return;
    value->x = halfLx + (halfLx-value->rad);
    value->vx = 3.0;
    resW = 0.0;
    resT = 0.0;
    start_single_particle_engine();
    process_immediate_wall();

    expect_close(value->vx, 0.0, "fully dissipative radial velocity");
    expect_close(value->vy, 0.0, "fully dissipative tangent velocity");
#if THREE_D
    expect_close(value->vz, 0.0, "fully dissipative second tangent velocity");
#endif
    expect_close(circularWallTimeNormal(value), never,
                 "stationary contact rescheduled a zero-time wall event");
    expect_true(findNextEvent() == NULL,
                "stationary contact retained a bogus crossing/wall event");
    expect_true(ncol == 1, "stationary contact processed more than one impact");
    expect_close(cumulativeWallImpulse, 3.0,
                 "fully dissipative wall impulse accounting");
    freeArrays();
}

static void test_zero_tangent_restitution_stops_cleanly(void)
{
    particle *value = configure_single_particle();
    if (value == NULL) return;
    value->x = halfLx + (halfLx-value->rad);
    value->vy = 2.0;
    resW = 1.0;
    resT = 0.0;
    start_single_particle_engine();
    process_immediate_wall();

    expect_close(value->vx, 0.0, "zero-resT tangent x velocity");
    expect_close(value->vy, 0.0, "zero-resT tangent y velocity");
#if THREE_D
    expect_close(value->vz, 0.0, "zero-resT tangent z velocity");
#endif
    expect_close(circularWallTimeNormal(value), never,
                 "zero-resT tangent contact rescheduled at zero time");
    expect_true(findNextEvent() == NULL,
                "zero-resT tangent contact retained a bogus event");
    expect_close(cumulativeWallImpulse, 0.0,
                 "pure tangent damping recorded a normal wall impulse");
    freeArrays();
}
#endif

#if EDMD_ADD_FIELD
static void test_tangent_field_contact_advances_in_time(void)
{
    particle *value = configure_single_particle();
    if (value == NULL) return;
    value->x = halfLx + (halfLx-value->rad);
    field = -0.02;
    resW = 0.0;
    resT = 0.0;
    start_single_particle_engine();
    process_immediate_wall();

    expect_true(value->vx < 0.0,
                "tangent field contact did not receive an inward regularization");
    expect_close(value->vy, 0.0,
                 "tangent field regularization changed y velocity");
#if THREE_D
    expect_close(value->vz, 0.0,
                 "tangent field regularization changed z velocity");
#endif
    expect_true(circularWallTimeNormal(value) > 0.0,
                "tangent field contact remained a zero-time wall event");
    node *following = findNextEvent();
    /* The equator lies on an axis-aligned cell face in this test.  Downward
       acceleration therefore commits one legitimate zero-time CELLCROSS to
       move the particle from the upper cell into the lower one.  It must not
       turn into a second zero-time WALL event or an unbounded crossing loop. */
    int immediateCrossings = 0;
    while (following != NULL && following->t == t &&
           following->type == CELLCROSS && immediateCrossings < 4) {
        nextEvent = following;
        removeEventFromQueue(nextEvent);
        doTheCrossing();
        immediateCrossings++;
        following = findNextEvent();
    }
    expect_true(immediateCrossings <= 1,
                "tangent field contact produced repeated zero-time cell crossings");
    expect_true(following != NULL && following->t > t,
                "tangent field contact failed to advance the event queue");

    field = 0.0;
    value->vx = value->vy = 0.0;
#if THREE_D
    value->vz = 0.0;
#endif
    expect_close(circularWallTimeNormal(value), never,
                 "zero constant field scheduled a stationary wall contact");
    freeArrays();
}
#endif

static void test_outside_repair_updates_3d_cell_membership(void)
{
    particle *value = configure_single_particle();
    if (value == NULL) return;
#if THREE_D
    const double direction = 1.0/sqrt(3.0);
#else
    const double direction = 1.0/sqrt(2.0);
#endif
    value->x = halfLx + 15.0*direction;
    value->y = halfLy + 15.0*direction;
    value->vx = -direction;
    value->vy = -direction;
#if THREE_D
    value->z = halfLz + 15.0*direction;
    value->vz = -direction;
#endif
    resW = 1.0;
    resT = 1.0;
    start_single_particle_engine();
    const int oldX = value->cell[0];
    const int oldY = value->cell[1];
#if THREE_D
    const int oldZ = value->cell[2];
#endif
    process_immediate_wall();

    const int expectedX = coordToCell(value->x, 1);
    const int expectedY = coordToCell(value->y, 0);
    expect_true(value->cell[0] == expectedX && value->cell[1] == expectedY,
                "projected particle retained stale x/y cell coordinates");
#if THREE_D
    const int expectedZ = coordToCell(value->z, 2);
    expect_true(value->cell[2] == expectedZ,
                "projected particle retained a stale z cell coordinate");
    expect_true(oldX != expectedX || oldY != expectedY || oldZ != expectedZ,
                "outside-repair regression did not cross a 3D cell face");
    const int cellIndex = expectedZ*Nycells*Nxcells +
                          expectedY*Nxcells + expectedX;
#else
    expect_true(oldX != expectedX || oldY != expectedY,
                "outside-repair regression did not cross a cell face");
    const int cellIndex = expectedY*Nxcells + expectedX;
#endif
    expect_true(cellList[cellIndex] == value,
                "projected particle was not relinked in its corrected cell");
    const double dx = value->x-halfLx;
    const double dy = value->y-halfLy;
    double distanceSquared = dx*dx + dy*dy;
#if THREE_D
    const double dz = value->z-halfLz;
    distanceSquared += dz*dz;
#endif
    const double accessible = halfLx-value->rad;
    expect_true(distanceSquared <= accessible*accessible*(1.0+1e-13),
                "outside repair left the particle beyond the spherical wall");
    expect_close(value->vx, -direction,
                 "inward outside repair spuriously reflected x velocity");
    expect_close(value->vy, -direction,
                 "inward outside repair spuriously reflected y velocity");
#if THREE_D
    expect_close(value->vz, -direction,
                 "inward outside repair spuriously reflected z velocity");
#endif
    expect_close(cumulativeWallImpulse, 0.0,
                 "inward outside repair recorded a false wall impulse");
    freeArrays();
}

static void test_off_axis_specular_response(void)
{
    particle *value = configure_single_particle();
    if (value == NULL) return;
#if THREE_D
    const double nx = 1.0/3.0;
    const double ny = 2.0/3.0;
    const double nz = 2.0/3.0;
    const double tx = 2.0;
    const double ty = -1.0;
    const double tz = 0.0;
#else
    const double nx = 3.0/5.0;
    const double ny = 4.0/5.0;
    const double tx = 4.0/5.0;
    const double ty = -3.0/5.0;
#endif
    const double accessible = halfLx-value->rad;
    value->x = halfLx + accessible*nx;
    value->y = halfLy + accessible*ny;
    value->vx = 2.0*nx + tx;
    value->vy = 2.0*ny + ty;
#if THREE_D
    value->z = halfLz + accessible*nz;
    value->vz = 2.0*nz + tz;
#endif
    const double oldEnergy = value->vx*value->vx + value->vy*value->vy
#if THREE_D
        + value->vz*value->vz
#endif
        ;
    resW = 1.0;
    resT = 1.0;
    start_single_particle_engine();
    process_immediate_wall();

    expect_close(value->vx, -2.0*nx + tx,
                 "off-axis specular x response");
    expect_close(value->vy, -2.0*ny + ty,
                 "off-axis specular y response");
#if THREE_D
    expect_close(value->vz, -2.0*nz + tz,
                 "off-axis specular z response");
#endif
    const double newEnergy = value->vx*value->vx + value->vy*value->vy
#if THREE_D
        + value->vz*value->vz
#endif
        ;
    expect_close(newEnergy, oldEnergy,
                 "off-axis specular wall changed kinetic energy");
    expect_close(cumulativeWallImpulse, 4.0,
                 "off-axis specular wall impulse accounting");
    freeArrays();
}

int main(void)
{
#if EDMD_ADD_FIELD
    test_tangent_field_contact_advances_in_time();
#else
    test_fully_dissipative_contact_stops_cleanly();
    test_zero_tangent_restitution_stops_cleanly();
#endif
    test_outside_repair_updates_3d_cell_membership();
    test_off_axis_specular_response();
    if (failures != 0) {
        fprintf(stderr, "%d spherical-wall regression(s) failed\n", failures);
        return 1;
    }
    printf("spherical-wall regressions passed (%dD, stress=%d, heat=%d)\n",
           THREE_D ? 3 : 2, STRESS, HEAT);
    return 0;
}
