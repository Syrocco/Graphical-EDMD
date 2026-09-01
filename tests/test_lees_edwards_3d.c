#include "EDMD.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !THREE_D
#error "test_lees_edwards_3d requires THREE_D=1"
#endif

extern int N, S1, Hex;
extern int periodicDumps, ther, Nxcells, Nycells, Nzcells;
extern particle *particles;
extern node **eventList;
extern node *nextEvent;
extern double t, paulTime, Lx, Ly, Lz, halfLy, strainRate, res;
extern double collTermZ, collTermYZ, collTermZY;
#if STRESS
extern double sumKineticX, sumKineticXY, sumKineticXZ, sumKineticZX;
extern double runningKineticX, runningKineticXY;
extern double runningKineticXZ, runningKineticZX;
#endif
#if HEAT
extern double sumHeatZ, runningHeatZ, collHeatZ, collisionEnergySource;
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
    if (!isfinite(actual) || fabs(actual - expected) > 2e-11*scale) {
        fprintf(stderr, "%s: got %.17g, expected %.17g\n",
                message, actual, expected);
        failures++;
    }
}

static void test_shear_transport_update(void)
{
#if (STRESS || HEAT)
    particle value;
    memset(&value, 0, sizeof(value));
    value.m = 2.0;
    value.x = 5.0;
    value.y = 10.0;
    value.z = 7.0;
    value.vx = 1.0;
    value.vy = 2.0;
    value.vz = 3.0;
    particles = &value;
    N = 1;
    S1 = Hex = 0;
    t = paulTime = 5.0;
    Lx = Ly = Lz = 20.0;
    halfLy = 10.0;
    strainRate = 0.1;
    initializeTransportAccumulators();
#if STRESS
    const double oldX = sumKineticX;
    const double oldXY = sumKineticXY;
    const double oldXZ = sumKineticXZ;
    const double oldZX = sumKineticZX;
#endif
#if HEAT
    const double oldHeatZ = sumHeatZ;
#endif
    /* Advance to an upper-y seam.  Just before remapping c_x has evolved
       continuously from 1 to 2; wrapping y from Ly to zero and adding
       strainRate*Ly to stored vx must leave that peculiar velocity intact. */
    t = paulTime = 10.0;
    value.y = 0.0;
    shearCorrection(&value.x, &value.vx, 1, &value);
    expect_close(value.vx, 3.0, "Lees-Edwards velocity jump");
    expect_close(value.vz, 3.0, "Lees-Edwards changed vz");
#if STRESS
    expect_close(sumKineticX - oldX, 6.0,
                 "Lees-Edwards continuous xx stress state");
    expect_close(sumKineticXY - oldXY, 4.0,
                 "Lees-Edwards continuous xy stress state");
    expect_close(sumKineticXZ - oldXZ, 6.0,
                 "Lees-Edwards continuous xz stress state");
    expect_close(sumKineticZX - oldZX, 6.0,
                 "Lees-Edwards continuous zx stress state");
    expect_close(runningKineticX, 70.0/3.0,
                 "Lees-Edwards integrated xx stress");
    expect_close(runningKineticXY, 30.0,
                 "Lees-Edwards integrated xy stress");
    expect_close(runningKineticXZ, 45.0,
                 "Lees-Edwards integrated xz stress");
    expect_close(runningKineticZX, 45.0,
                 "Lees-Edwards integrated zx stress");
#endif
#if HEAT
    expect_close(sumHeatZ - oldHeatZ, 9.0,
                 "Lees-Edwards continuous z heat-current state");
    expect_close(runningHeatZ, 230.0,
                 "Lees-Edwards integrated z heat current");
#endif
    const double upper_peculiar_vx = value.vx +
        strainRate*(value.y - halfLy);
    value.y = Ly;
    shearCorrection(&value.x, &value.vx, -1, &value);
    expect_close(value.vx + strainRate*(value.y - halfLy),
                 upper_peculiar_vx,
                 "lower Lees-Edwards seam changed peculiar vx");
    particles = NULL;
#endif
}

static void test_lower_seam_crossing(void)
{
    N = 1;
    S1 = Hex = 0;
    periodicDumps = ther = 0;
    t = paulTime = 0.0;
    Lx = Ly = Lz = 20.0;
    strainRate = 0.1;
    particles = calloc(1, sizeof(*particles));
    expect_true(particles != NULL,
                "unable to allocate lower-seam particle");
    if (particles == NULL) return;

    particles[0].num = 0;
    particles[0].m = particles[0].rad = 1.0;
    particles[0].x = 5.0;
    particles[0].y = 0.5;
    particles[0].z = 5.0;
    particles[0].vy = -1.0;
    boxConstantHelper();
    cellListInit();
#if (STRESS || HEAT)
    initializeTransportAccumulators();
#endif
    eventListInit();

    node *crossing = eventList[0];
    expect_true(crossing->type == CELLCROSS && crossing->j == 3,
                "lower-y Lees crossing was not scheduled first");
    const unsigned long old_collision_generation = particles[0].coll;
    const double peculiar_before = particles[0].vx +
        strainRate*(0.0 - halfLy);
    t = crossing->t;
    removeEventFromQueue(crossing);
    nextEvent = crossing;
    doTheCrossing();
    expect_close(particles[0].y, Ly,
                 "lower Lees crossing did not select the upper image");
    expect_close(particles[0].vx +
                 strainRate*(particles[0].y - halfLy), peculiar_before,
                 "lower Lees crossing changed peculiar vx");
    expect_true(particles[0].coll == old_collision_generation + 1,
                "Lees image jump did not invalidate inbound collision events");
    expect_true(particles[0].cell[1] == Nycells - 1,
                "lower Lees crossing entered the wrong y cell");
    freeArrays();
}

static void test_shifted_seam_collision(void)
{
    N = 2;
    S1 = Hex = 0;
    periodicDumps = ther = 0;
    t = paulTime = 5.0;
    Lx = Ly = Lz = 20.0;
    strainRate = 0.1;
    res = 1.0;
    particles = calloc((size_t)N, sizeof(*particles));
    expect_true(particles != NULL,
                "unable to allocate Lees-Edwards seam particles");
    if (particles == NULL) return;

    particles[0].num = 0;
    particles[0].m = particles[0].rad = 1.0;
    particles[0].x = 1.0;
    particles[0].y = 1.5;
    particles[0].z = 0.5;
    particles[0].vy = -0.5;
    particles[0].t = t;

    particles[1].num = 1;
    particles[1].m = particles[1].rad = 1.0;
    particles[1].x = 11.0;
    particles[1].y = 18.5;
    particles[1].z = 19.5;
    particles[1].vx = -2.0;
    particles[1].vy = 0.5;
    particles[1].t = t;

    boxConstantHelper();
    cellListInit();
#if (STRESS || HEAT)
    initializeTransportAccumulators();
#endif
    eventListInit();

    expect_true(particles[0].cell[0] != particles[1].cell[0],
                "seam regression did not require a shifted x-cell scan");
    expect_true(particles[0].cell[2] == 0 &&
                particles[1].cell[2] == Nzcells - 1,
                "seam regression did not require periodic z neighbours");
    node *collision = eventList[N];
    expect_true(collision->type == COLLISION && collision->j == 1,
                "3D Lees-Edwards seam collision was not scheduled");
    expect_close(collision->t - t, 3.0 - sqrt(3.0),
                 "3D Lees-Edwards seam collision time");

    nextEvent = collision;
    t = collision->t;
    removeEventFromQueue(collision);
    doTheCollisionNormal();

    double dx = particles[1].x - particles[0].x;
    double dy = particles[1].y - particles[0].y;
    double dz = particles[1].z - particles[0].z;
    double dvx = particles[1].vx - particles[0].vx;
    correctDistances(particles, particles + 1, particles[1].x,
                     particles[1].vx, &dx, &dvx);
    PBC(&dx, &dy, &dz);
    const double dvy = particles[1].vy - particles[0].vy;
    const double dvz = particles[1].vz - particles[0].vz;
    expect_true(dx*dvx + dy*dvy + dz*dvz > 0.0,
                "Lees-Edwards collision did not reverse the image-normal velocity");
    expect_true(isfinite(collTermZ) && isfinite(collTermYZ) &&
                isfinite(collTermZY) && fabs(collTermZ) > 0.0,
                "Lees-Edwards collision omitted finite z virial accounting");
#if HEAT
    expect_true(isfinite(collHeatZ) && isfinite(collisionEnergySource),
                "Lees-Edwards collision produced invalid heat accounting");
    expect_close(collisionEnergySource, 0.0,
                 "elastic Lees-Edwards collision was marked as an energy source");
#endif
    freeArrays();
}

int main(void)
{
    test_shear_transport_update();
    test_lower_seam_crossing();
    test_shifted_seam_collision();
    if (failures != 0) {
        fprintf(stderr, "%d Lees-Edwards regression(s) failed\n", failures);
        return 1;
    }
    printf("3D Lees-Edwards regressions passed\n");
    return 0;
}
