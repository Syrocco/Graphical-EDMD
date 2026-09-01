#include "EDMD.h"
#include "live_analysis.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int N = 2;
static particle storage[2];
particle *particles = storage;
double t = 0.0;
double Lx = 100.0, Ly = 100.0, halfLx = 50.0, halfLy = 50.0;
#if THREE_D
double Lz = 100.0;
#endif
bool addWell = false;
bool addField = false;
bool charged = false;
int liquidliquid = 0;
double U = -2.0;
double field = -0.1;
int addWallx = 0, addWally = 0, addCircularWall = 0, addMidWall = 0;
#if THREE_D
int addWallz = 0;
#endif
bool addShear = false;
double strainRate = 0.0;
double collTermX = 0.0, collTermY = 0.0;
#if THREE_D
double collTermZ = 0.0;
#endif
double lastScreen = 0.0;
double dp = 0.0;
double cumulativeWallImpulse = 0.0;

static int failures = 0;

static void expect_true(const char *name, bool value)
{
    if (!value) {
        fprintf(stderr, "%s: expected true\n", name);
        failures++;
    }
}

static void expect_close(const char *name, double actual, double expected)
{
    const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
    if (!isfinite(actual) || fabs(actual - expected) > 1e-10*scale) {
        fprintf(stderr, "%s: got %.17g, expected %.17g\n",
                name, actual, expected);
        failures++;
    }
}

static double latest(const LiveAnalysisSeries *series)
{
    if (series->count == 0) return NAN;
    const size_t index = (series->first + series->count - 1) %
                         LIVE_ANALYSIS_HISTORY_CAPACITY;
    return series->value[index];
}

void freeFlyNormal(particle *value)
{
    (void)value;
}

void freeFlyGrow(particle *value)
{
    (void)value;
}

void (*freeFly)(particle *) = &freeFlyNormal;

static particlePrediction predict(const particle *value, double target,
                                  bool growing)
{
    const double dt = target - value->t;
    particlePrediction result;
    memset(&result, 0, sizeof(result));
    result.x = value->x + dt*value->vx;
    result.y = value->y + dt*value->vy;
    result.vx = value->vx;
    result.vy = value->vy;
    result.rad = value->rad + (growing ? dt*value->vr : 0.0);
    if (addField && !growing) {
        result.y += 0.5*field*dt*dt;
        result.vy += field*dt;
    }
#if THREE_D
    result.z = value->z + dt*value->vz;
    result.vz = value->vz;
#endif
    return result;
}

particlePrediction predictFreeFlyNormal(const particle *value, double target)
{
    return predict(value, target, false);
}

particlePrediction predictFreeFlyGrow(const particle *value, double target)
{
    return predict(value, target, true);
}

#if THREE_D
void PBC(double *dx, double *dy, double *dz)
{
    if (!addWallx && !addCircularWall) {
        if (*dx >= 0.5*Lx) *dx -= Lx;
        else if (*dx < -0.5*Lx) *dx += Lx;
    }
    if (!addWally && !addCircularWall) {
        if (*dy >= 0.5*Ly) *dy -= Ly;
        else if (*dy < -0.5*Ly) *dy += Ly;
    }
    if (!addCircularWall) {
        if (*dz >= 0.5*Lz) *dz -= Lz;
        else if (*dz < -0.5*Lz) *dz += Lz;
    }
}
#else
void PBC(double *dx, double *dy)
{
    if (!addWallx && !addCircularWall) {
        if (*dx >= 0.5*Lx) *dx -= Lx;
        else if (*dx < -0.5*Lx) *dx += Lx;
    }
    if (!addWally && !addCircularWall) {
        if (*dy >= 0.5*Ly) *dy -= Ly;
        else if (*dy < -0.5*Ly) *dy += Ly;
    }
}
#endif

static void configure_particles(void)
{
    memset(storage, 0, sizeof(storage));
    for (int index = 0; index < N; ++index) {
        storage[index].num = index;
        storage[index].m = 1.0;
        storage[index].rad = 1.0;
        storage[index].y = 10.0;
#if THREE_D
        storage[index].z = 10.0;
#endif
    }
    storage[0].x = 25.0;
    storage[1].x = 35.0;
    storage[0].vx = 1.0;
    storage[1].vx = -1.0;
}

int main(void)
{
    configure_particles();
    live_analysis_initialize();
    live_analysis_update(0.0);
    const LiveAnalysisData *data = live_analysis_get_data();
    expect_true("initial energy sample", data->energy_per_particle.count == 1);
    expect_close("kinetic energy per particle",
                 latest(&data->energy_per_particle), 0.5);
    expect_close("translational temperature", latest(&data->temperature),
                 2.0/(THREE_D ? 6.0 : 4.0));
    expect_true("pressure warmup", data->pressure.count == 0 &&
                                      data->pressure_warming_up);

    t = 1.0;
    collTermX = -2.0;
    live_analysis_update(1.0);
    data = live_analysis_get_data();
    const double volume = Lx*Ly
#if THREE_D
        *Lz
#endif
        ;
    const double expected_pressure = 2.0/(THREE_D ? 3.0*volume
                                                  : 2.0*volume) +
                                     2.0/((THREE_D ? 3.0 : 2.0)*volume);
    expect_close("periodic virial pressure", latest(&data->pressure),
                 expected_pressure);

    int first_well[] = {1};
    int second_well[] = {0};
    storage[0].particlesInWell = first_well;
    storage[1].particlesInWell = second_well;
    storage[0].numberOfParticlesInWell = 1;
    storage[1].numberOfParticlesInWell = 1;
    addWell = true;
    t = 2.0;
    live_analysis_update(2.0);
    data = live_analysis_get_data();
    expect_close("square-step total energy", latest(&data->energy_per_particle),
                 -0.5);

    addField = true;
    for (int index = 0; index < N; ++index) storage[index].t = 3.0;
    t = 3.0;
    live_analysis_update(3.0);
    data = live_analysis_get_data();
    expect_close("field potential energy", latest(&data->energy_per_particle),
                 0.5);
    addField = false;
    addWell = false;
    storage[0].numberOfParticlesInWell = 0;
    storage[1].numberOfParticlesInWell = 0;

    LiveAnalysisSettings settings;
    live_analysis_get_settings(&settings);
    settings.mean_squared_displacement = true;
    live_analysis_apply_settings(&settings);
    storage[0].x += 1.0;
    storage[1].x -= 1.0;
    for (int index = 0; index < N; ++index) storage[index].t = 4.0;
    t = 4.0;
    live_analysis_update(4.0);
    data = live_analysis_get_data();
    expect_true("MSD reference ready", data->msd_reference_ready);
    expect_close("COM-corrected MSD",
                 latest(&data->mean_squared_displacement), 1.0);
    live_analysis_clear_histories();
    data = live_analysis_get_data();
    expect_true("clearing MSD history does not fake a new zero baseline",
                data->mean_squared_displacement.count == 0 &&
                data->msd_reference_ready);

    live_analysis_get_settings(&settings);
    settings.pair_correlation = true;
    settings.pair_bins = 16;
    live_analysis_apply_settings(&settings);
    t = 5.0;
    live_analysis_update(5.0);
    data = live_analysis_get_data();
    expect_true("pair correlation sample",
                data->pair_correlation.valid &&
                data->pair_correlation.averages == 1);

    addCircularWall = 1;
    dp = 0.0;
    cumulativeWallImpulse = 0.0;
    live_analysis_reset_generation();
    live_analysis_update(6.0);
    t = 6.0;
    const double wall_measure =
#if THREE_D
        4.0*M_PI*halfLx*halfLx;
#else
        M_PI*Lx;
#endif
    cumulativeWallImpulse = 3.0*wall_measure;
    dp = 0.0; /* Saved-output interval resets must not affect the live plot. */
    live_analysis_update(7.0);
    data = live_analysis_get_data();
    expect_close("confining-wall pressure", latest(&data->pressure), 3.0);

#if THREE_D
    addCircularWall = 0;
    addWallx = addWally = addWallz = 1;
    cumulativeWallImpulse = 0.0;
    live_analysis_reset_generation();
    live_analysis_update(7.1);
    t = 7.0;
    const double xyz_wall_measure = 2.0*(Ly*Lz + Lx*Lz + Lx*Ly);
    cumulativeWallImpulse = 2.5*xyz_wall_measure;
    live_analysis_update(7.5);
    data = live_analysis_get_data();
    expect_close("full XYZ wall pressure", latest(&data->pressure), 2.5);
    addWallx = addWally = addWallz = 0;
#endif

    const uint64_t generation = data->generation;
    freeFly = &freeFlyGrow;
    t = 7.0;
    live_analysis_update(8.0);
    data = live_analysis_get_data();
    expect_true("flight phase resets generation",
                data->generation == generation + 1);

    freeFly = &freeFlyNormal;
    addShear = true;
    live_analysis_reset_generation();
    data = live_analysis_get_data();
    expect_true("shear disables MSD", !data->msd_supported &&
                                      !data->msd_reference_ready);
    expect_true("shear disables orthorhombic g(r)",
                !data->pair_correlation_supported);

    live_analysis_shutdown();
    if (failures != 0) {
        fprintf(stderr, "%d live-analysis regression(s) failed\n", failures);
        return 1;
    }
    printf("live-analysis regressions passed (%s)\n",
           THREE_D ? "3D" : "2D");
    return 0;
}
