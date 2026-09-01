#include "EDMD.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

extern double t;
extern double Lx, Ly, halfLx, halfLy;
#if THREE_D
extern double Lz, halfLz;
#endif

static int failures = 0;

static void expect_close(const char *name, double actual, double expected,
                         double tolerance)
{
    const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
    if (!isfinite(actual) || fabs(actual - expected) > tolerance*scale) {
        fprintf(stderr, "%s: got %.17g, expected %.17g\n",
                name, actual, expected);
        failures++;
    }
}

static particle blank_particle(void)
{
    particle value;
    memset(&value, 0, sizeof(value));
    value.rad = 1.0;
    value.m = 1.0;
    return value;
}

static void configure_box(void)
{
    t = 0.0;
    Lx = Ly = 100.0;
    halfLx = halfLy = 50.0;
#if THREE_D
    Lz = 100.0;
    halfLz = 50.0;
#endif
}

static void test_direct_collision(void)
{
    particle left = blank_particle();
    particle right = blank_particle();

#if THREE_D
    left.x = left.y = left.z = 10.0;
    right.x = 13.0;
    right.y = 14.0;
    right.z = 10.0;
    right.vx = -0.6;
    right.vy = -0.8;
    expect_close("3D direct collision", collisionTimeNormal(&left, &right),
                 3.0, 1e-12);
#else
    left.x = left.y = 10.0;
    right.x = 15.0;
    right.y = 10.0;
    left.vx = 1.0;
    right.vx = -1.0;
    expect_close("2D direct collision", collisionTimeNormal(&left, &right),
                 1.5, 1e-12);
#endif
}

static void test_periodic_collision(void)
{
    particle first = blank_particle();
    particle second = blank_particle();
    first.x = 2.0;
    second.x = 98.0;
    first.y = second.y = 20.0;
    first.vx = -1.0;
    second.vx = 1.0;
#if THREE_D
    first.z = second.z = 20.0;
#endif
    expect_close("periodic collision", collisionTimeNormal(&first, &second),
                 1.0, 1e-12);
}

static void test_growth_collision(void)
{
    particle first = blank_particle();
    particle second = blank_particle();
    first.x = first.y = 10.0;
    second.x = 13.0;
    second.y = 10.0;
    first.vr = second.vr = 0.1;
#if THREE_D
    first.z = second.z = 10.0;
#endif
    expect_close("growth collision", collisionTimeGrow(&first, &second),
                 5.0, 1e-11);
}

static void test_circular_wall_time(void)
{
    particle value = blank_particle();
    value.x = halfLx;
    value.y = halfLy;
#if THREE_D
    value.z = halfLz;
    value.vx = 1.0;
    value.vy = 2.0;
    value.vz = 2.0;
    expect_close("3D spherical wall time", circularWallTimeNormal(&value),
                 (halfLx-value.rad)/3.0, 1e-12);
#else
    value.vx = 1.0;
    expect_close("2D circular wall time", circularWallTimeNormal(&value),
                 halfLx-value.rad, 1e-12);
#endif

    /* Near-surface exits must not be rounded down to the zero contact root or
       skipped in favour of the far side of the container. */
    value = blank_particle();
    const double accessible = halfLx - value.rad;
    const double near_gap = 4e-13;
    value.x = halfLx + accessible - near_gap;
    value.y = halfLy;
    value.vx = 1.0;
#if THREE_D
    value.z = halfLz;
#endif
    expect_close("near-surface circular exit", circularWallTimeNormal(&value),
                 near_gap, 2e-14);

    value.x = halfLx + accessible;
    value.vx = -1.0;
    expect_close("surface inward circular return",
                 circularWallTimeNormal(&value), 2.0*accessible, 1e-12);

    value.vx = 0.0;
    value.vy = 1.0;
    expect_close("surface tangent circular contact",
                 circularWallTimeNormal(&value), 0.0, 1e-13);

    value.x = halfLx + accessible + 1e-9;
    value.vx = -1.0;
    value.vy = 0.0;
    expect_close("outside circular state repair",
                 circularWallTimeNormal(&value), 0.0, 1e-13);

    /* Velocity must not enter the geometric surface tolerance.  The previous
       mixed-units tolerance grew as v^2: it treated a genuinely interior fast
       particle as touching the wall and discarded the short far-side root of
       an inward-moving surface particle. */
    value.x = halfLx + accessible - 1.0;
    value.vx = 1e8;
    expect_close("high-speed interior spherical exit",
                 circularWallTimeNormal(&value), 1e-8, 1e-12);

    value.x = halfLx + accessible;
    value.vx = -1e8;
    expect_close("high-speed inward spherical return",
                 circularWallTimeNormal(&value), 2.0*accessible/1e8, 1e-12);

    value.vx = 0.0;
    value.vy = 0.0;
#if THREE_D
    value.vz = 0.0;
#endif
    expect_close("stationary spherical contact has no event",
                 circularWallTimeNormal(&value), never, 1e-15);
}

static void test_growing_wall_time_and_response(void)
{
    particle value = blank_particle();
    value.x = halfLx;
    value.y = halfLy;
    value.vx = 0.25;
    value.vr = 0.5;
#if THREE_D
    value.z = halfLz;
#endif
    expect_close("contracting circular wall time", circularWallTimeGrow(&value),
                 (halfLx-value.rad)/(value.vx+value.vr), 1e-12);

    value.x = halfLx + halfLx - value.rad;
    value.vx = -0.25;
    expect_close("growth overtakes inward particle",
                 circularWallTimeGrow(&value), 0.0, 1e-13);

    value.vr = 0.0;
    value.vx = 0.0;
    expect_close("stationary grown particle has no wall event",
                 circularWallTimeGrow(&value), never, 1e-15);

    value.vy = 1.0;
    expect_close("growing solver retains tangent contact",
                 circularWallTimeGrow(&value), 0.0, 1e-13);
    value.vy = 0.0;

    /* At the left wall vx=+0.1 still approaches the growing radius because
       vr=0.2. In outward-normal coordinates vn=-0.1; the reflected vx must
       be +0.3 rather than the old erroneous -0.5. */
    expect_close("growing left-wall normal response",
                 confiningWallNormalVelocityAfter(-0.1, 0.2, 1.0),
                 -0.3, 1e-13);
    expect_close("growing right-wall normal response",
                 confiningWallNormalVelocityAfter(0.1, 0.2, 1.0),
                 -0.5, 1e-13);
}

static void test_explicit_prediction(void)
{
    particle value = blank_particle();
    value.x = 1.0;
    value.y = 2.0;
    value.vx = 3.0;
    value.vy = -2.0;
#if THREE_D
    value.z = 4.0;
    value.vz = 1.5;
#endif
    const double gamma = 0.5;
    const double decay = exp(-gamma*2.0);
    particlePrediction prediction = predictParticle(
        &value, 2.0, LINEAR_DAMPING_FLIGHT, gamma);
    expect_close("damped x", prediction.x,
                 value.x + value.vx/gamma*(1.0 - decay), 1e-13);
    expect_close("damped y", prediction.y,
                 value.y + value.vy/gamma*(1.0 - decay), 1e-13);
    expect_close("damped vx", prediction.vx, value.vx*decay, 1e-13);
    expect_close("damped vy", prediction.vy, value.vy*decay, 1e-13);
#if THREE_D
    expect_close("damped z", prediction.z,
                 value.z + value.vz/gamma*(1.0 - decay), 1e-13);
    expect_close("damped vz", prediction.vz, value.vz*decay, 1e-13);
#endif
}

static void test_potential_step_solver(void)
{
    double first_after = 0.0;
    double second_after = 0.0;
    bool crossed = potentialStepNormalVelocities(
        1.0, 1.0, 1.0, -1.0, 0.5, &first_after, &second_after);
    if (!crossed) {
        fprintf(stderr, "above-barrier repulsive step was reflected\n");
        failures++;
    }
    expect_close("repulsive step first direction", first_after,
                 sqrt(0.5), 1e-13);
    expect_close("repulsive step second direction", second_after,
                 -sqrt(0.5), 1e-13);
    expect_close("repulsive step momentum", first_after+second_after,
                 0.0, 1e-13);
    expect_close("repulsive step total energy",
                 0.5*(first_after*first_after+second_after*second_after)+0.5,
                 1.0, 1e-13);

    crossed = potentialStepNormalVelocities(
        1.0, 1.0, 0.5, -0.5, 1.0, &first_after, &second_after);
    if (crossed) {
        fprintf(stderr, "below-barrier repulsive step was transmitted\n");
        failures++;
    }
    expect_close("reflected step first", first_after, -0.5, 1e-13);
    expect_close("reflected step second", second_after, 0.5, 1e-13);

    crossed = potentialStepNormalVelocities(
        1.0, 1.0, -0.5, 0.5, -0.5, &first_after, &second_after);
    if (!crossed) {
        fprintf(stderr, "downhill potential exit was reflected\n");
        failures++;
    }
    expect_close("downhill exit first direction", first_after,
                 -sqrt(0.75), 1e-13);
    expect_close("downhill exit second direction", second_after,
                 sqrt(0.75), 1e-13);
    expect_close("downhill exit total energy",
                 0.5*(first_after*first_after+second_after*second_after)-0.5,
                 0.25, 1e-13);
}

int main(void)
{
    configure_box();
    test_direct_collision();
    test_periodic_collision();
    test_growth_collision();
    test_circular_wall_time();
    test_growing_wall_time_and_response();
    test_explicit_prediction();
    test_potential_step_solver();

    if (failures != 0) {
        fprintf(stderr, "%d prediction regression(s) failed\n", failures);
        return 1;
    }
    printf("prediction regressions passed (%dD)\n", THREE_D ? 3 : 2);
    return 0;
}
