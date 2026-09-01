#include "EDMD.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

extern int N;
extern int Nthermo;
extern int addWallx;
extern unsigned long int ncol;
extern particle *particles;
extern FILE *thermo;
extern double t, lastScreen;
extern double Lx, Ly, halfLx, halfLy;
#if THREE_D
extern double Lz, halfLz;
#endif
extern double collTermX, collTermY, collTermXY, collTermYX;
#if THREE_D
extern double collTermZ, collTermXZ, collTermZX, collTermYZ, collTermZY;
#endif
extern double dp, dpLeft, dpRight, dpMid;

static int failures;

static void expect_close(const char *name, double actual, double expected)
{
    const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
    if (!isfinite(actual) || fabs(actual-expected) > 1e-12*scale) {
        fprintf(stderr, "%s: got %.17g, expected %.17g\n",
                name, actual, expected);
        failures++;
    }
}

static void expect_true(const char *name, int condition)
{
    if (!condition) {
        fprintf(stderr, "%s\n", name);
        failures++;
    }
}

static double configured_wall_measure(void)
{
#if THREE_D
    return 2.0*Ly*Lz;
#else
    return 2.0*Ly;
#endif
}

static void configure_state(void)
{
    static particle value;
    memset(&value, 0, sizeof(value));
    value.m = 1.0;
    value.rad = 1.0;
    particles = &value;
    N = 1;
    addWallx = 1;
    ncol = 0;
    Lx = 20.0;
    Ly = 10.0;
    halfLx = 10.0;
    halfLy = 5.0;
#if THREE_D
    Lz = 5.0;
    halfLz = 2.5;
#endif
    t = 0.0;
    resetOutputSamplingState();
}

static void test_output_reset_is_an_interval_reset(void)
{
    cumulativeWallImpulse = 123.0;
    dp = 7.0;
    dpLeft = 1.0;
    dpRight = 2.0;
    dpMid = 3.0;
    collTermX = 4.0;
    collTermY = 5.0;
    collTermXY = 6.0;
    collTermYX = 7.0;
#if THREE_D
    collTermZ = 8.0;
    collTermXZ = 9.0;
    collTermZX = 10.0;
    collTermYZ = 11.0;
    collTermZY = 12.0;
#endif
    t = 2.5;
    resetOutputSamplingState();
    expect_close("output reset time", lastScreen, t);
    expect_close("interval wall impulse reset", dp, 0.0);
    expect_close("left impulse reset", dpLeft, 0.0);
    expect_close("right impulse reset", dpRight, 0.0);
    expect_close("membrane impulse reset", dpMid, 0.0);
    expect_close("x virial reset", collTermX, 0.0);
    expect_close("y virial reset", collTermY, 0.0);
#if THREE_D
    expect_close("z virial reset", collTermZ, 0.0);
#endif
    expect_close("monotonic wall impulse preserved",
                 cumulativeWallImpulse, 123.0);
}

static void test_nthermo_wall_pressure_average(void)
{
    t = 0.0;
    resetOutputSamplingState();
    cumulativeWallImpulse = 17.0;
    Nthermo = 2;
    thermo = tmpfile();
    expect_true("unable to create temporary thermo stream", thermo != NULL);
    if (thermo == NULL) return;

    const double measure = configured_wall_measure();
    t = 1.0;
    dp = measure;
    saveThermo();
    expect_close("first thermo sample clears dp", dp, 0.0);

    /* The second interval is twice as long and contains twice the impulse, so
     * both interval pressures and their Nthermo average are exactly one. */
    t = 3.0;
    dp = 2.0*measure;
    saveThermo();
    expect_close("second thermo sample clears dp", dp, 0.0);
    expect_close("thermo sampling preserves monotonic impulse",
                 cumulativeWallImpulse, 17.0);

    fflush(thermo);
    rewind(thermo);
    double columns[32];
    int count = 0;
    while (count < (int)(sizeof(columns)/sizeof(columns[0])) &&
           fscanf(thermo, "%lf", columns + count) == 1) {
        count++;
    }
#if THREE_D
    const int wall_column = 13
#if HEAT
        + 4
#endif
        ;
#else
    const int wall_column = 8
#if HEAT
        + 2
#endif
        ;
#endif
    expect_true("thermo row did not contain the wall-pressure column",
                count > wall_column);
    if (count > wall_column)
        expect_close("Nthermo wall pressure average",
                     columns[wall_column], 1.0);

    fclose(thermo);
    thermo = NULL;
}

int main(void)
{
    configure_state();
    test_output_reset_is_an_interval_reset();
    test_nthermo_wall_pressure_average();
    if (failures != 0) {
        fprintf(stderr, "%d wall-pressure regression(s) failed\n", failures);
        return 1;
    }
    printf("wall-pressure regressions passed (%dD)\n", THREE_D ? 3 : 2);
    return 0;
}
