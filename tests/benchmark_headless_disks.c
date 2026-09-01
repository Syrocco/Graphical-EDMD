#include "EDMD.h"
#include "graphics.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if G != 1
#error "benchmark_headless_disks needs mutable G=1 dynamics globals"
#endif
#if THREE_D
#error "benchmark_headless_disks is intentionally a 2D benchmark"
#endif
#if POLYGON_PARTICLES
#error "benchmark_headless_disks measures disks only"
#endif

/* EDMD.c is compiled with the mutable graphical configuration for this
 * benchmark, but no window or renderer is linked.  These functions satisfy
 * the unused graphical host while leaving the measured event loop headless. */
window graphicalInit(void)
{
    window result;
    memset(&result, 0, sizeof(result));
    return result;
}

state GUIinit(void)
{
    state result;
    memset(&result, 0, sizeof(result));
    return result;
}

void getInput(window *screen_window, state *screen_state)
{
    (void)screen_window;
    (void)screen_state;
}

void draw(int argc, char **argv, window *screen_window, state *screen_state)
{
    (void)argc;
    (void)argv;
    (void)screen_window;
    (void)screen_state;
}

void graphicsFree(state *screen_state) { (void)screen_state; }
void graphicsPauseSimulation(state *screen_state) { (void)screen_state; }
double graphicsSimulationTimeLimit(double current_time) { return current_time; }
void graphicsRebaseSimulationPacing(double current_time) { (void)current_time; }
void graphicsReportRuntimeError(const char *message) { (void)message; }
void graphicsReportSimulationFault(const char *message) { (void)message; }
void graphicsReportTransientNotice(const char *message) { (void)message; }
int graphicsTargetFrameRate(void) { return 60; }

void CloseWindow(void) {}
bool WindowShouldClose(void) { return true; }
int GetScreenHeight(void) { return 900; }
double GetTime(void) { return 0.0; }
float GetFrameTime(void) { return 1.0f/60.0f; }
void TraceLog(int log_level, const char *text, ...)
{
    (void)log_level;
    (void)text;
}

extern int N, Nbig, Nsmall, load, Hex, S1, version;
extern int periodicDumps, ther, automaticGrowthRate;
extern int Nxcells, Nycells, controlLenght, actualPaulList;
extern int addWallx, addWally, addWallz, addMidWall, addCircularWall, noise;
extern bool addWell, addField, addShear, damping;
extern bool addDelta, addDeltaTangent, addDoubleDelta;
extern bool addDoubleDeltaTangent, polydispersity;
extern int addEvolvingDelta, addExpo, addEnergy, addActiveDelta;
extern double phi, Einit, t, paulTime, E, res, resW, resT;
extern double gamm, delta, aspectRatio, fractionSmallN, Lx, Ly;
extern unsigned long int ncol, ncross;
extern particle *particles;
extern node *nextEvent;
extern void (*doTheCollision)(void);
extern void (*doTheWall)(void);
extern void init_genrand(unsigned long seed);

typedef enum BenchmarkMode {
    MODE_BALLISTIC,
    MODE_DAMPING,
    MODE_DELTA,
    MODE_DAMPING_DELTA
} BenchmarkMode;

typedef struct EventCounters {
    unsigned long events;
    unsigned long physical_collisions;
    unsigned long stale_collisions;
    unsigned long crossings;
} EventCounters;

static double monotonic_seconds(void)
{
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return NAN;
    return (double)now.tv_sec + 1.0e-9*(double)now.tv_nsec;
}

static void fail(const char *message)
{
    fprintf(stderr,
            "headless disk benchmark failure: %s "
            "(N=%d phi=%.9g t=%.17g ncol=%lu)\n",
            message, N, phi, t, ncol);
    if (simulationHasFault())
        fprintf(stderr, "simulation fault: %s\n", simulationFaultMessage());
    exit(EXIT_FAILURE);
}

static bool parse_mode(const char *name, BenchmarkMode *mode)
{
    if (strcmp(name, "ballistic") == 0) *mode = MODE_BALLISTIC;
    else if (strcmp(name, "damping") == 0) *mode = MODE_DAMPING;
    else if (strcmp(name, "delta") == 0) *mode = MODE_DELTA;
    else if (strcmp(name, "damping-delta") == 0)
        *mode = MODE_DAMPING_DELTA;
    else return false;
    return true;
}

static const char *mode_name(BenchmarkMode mode)
{
    switch (mode) {
        case MODE_BALLISTIC: return "ballistic";
        case MODE_DAMPING: return "damping";
        case MODE_DELTA: return "delta";
        case MODE_DAMPING_DELTA: return "damping-delta";
    }
    return "unknown";
}

static void select_mode(BenchmarkMode mode)
{
    damping = mode == MODE_DAMPING || mode == MODE_DAMPING_DELTA;
    addDelta = mode == MODE_DELTA || mode == MODE_DAMPING_DELTA;
    res = addDelta ? 0.90 : 1.0;
    gamm = 0.02;
    delta = 0.50;
}

static enum event dispatch_one(EventCounters *counters)
{
    nextEvent = findNextEvent();
    if (nextEvent == NULL || !isfinite(nextEvent->t) ||
        nextEvent->t < t || nextEvent->t >= never)
        fail("calendar has no finite forward event");

    const enum event type = nextEvent->type;
    const unsigned long collisions_before = ncol;
    t = nextEvent->t;
    removeEventFromQueue(nextEvent);
    switch (type) {
        case COLLISION:
            doTheCollision();
            break;
        case CELLCROSS:
            ++ncross;
            doTheCrossing();
            break;
        case WALL:
            doTheWall();
            break;
        case GROWSTOP:
            stopGrow();
            break;
        default:
            fail("unexpected auxiliary event in output-free benchmark");
    }
    if (simulationHasFault()) fail("event processing raised a fault");

    if (counters != NULL) {
        counters->events++;
        if (type == COLLISION) {
            if (ncol != collisions_before)
                counters->physical_collisions += ncol-collisions_before;
            else
                counters->stale_collisions++;
        } else if (type == CELLCROSS) {
            counters->crossings++;
        }
    }
    return type;
}

static void reset_configuration(int particle_count, double packing_fraction,
                                unsigned long seed, BenchmarkMode mode)
{
    N = particle_count;
    phi = packing_fraction;
    Einit = 1.0;
    version = (int)seed;
    load = Hex = S1 = 0;
    Nbig = Nsmall = 0;
    controlLenght = 0;
    aspectRatio = 1.0;
    fractionSmallN = 0.0;
    automaticGrowthRate = 1;
    periodicDumps = ther = 0;

    addWallx = addWally = addWallz = 0;
    addMidWall = addCircularWall = 0;
    addWell = addField = addShear = false;
    noise = 0;
    addDeltaTangent = addDoubleDelta = addDoubleDeltaTangent = false;
    addEvolvingDelta = addExpo = addEnergy = addActiveDelta = 0;
    polydispersity = false;
    resW = resT = 1.0;
    select_mode(mode);

    t = paulTime = 0.0;
    actualPaulList = 0;
    ncol = ncross = 0;
    simulationClearFault();
    constantInit(0, NULL);
    init_genrand(seed);
    srand((unsigned int)seed);
}

static void run_once(BenchmarkMode mode, int particle_count,
                     double packing_fraction, unsigned long timed_events,
                     unsigned long warmup_events, unsigned long seed,
                     int repetition)
{
    reset_configuration(particle_count, packing_fraction, seed, mode);

    const double setup_started = monotonic_seconds();
    particlesInit();
    adaptCellGridToParticles();
    cellListInit();
    physicalQ();
    eventListInit();
    const double setup_seconds = monotonic_seconds()-setup_started;

    EventCounters growth = {0};
    const double growth_started = monotonic_seconds();
    while (1) {
        nextEvent = findNextEvent();
        if (nextEvent == NULL) fail("growth calendar became empty");
        if (nextEvent->type == GROWSTOP) {
            (void)dispatch_one(NULL);
            break;
        }
        (void)dispatch_one(&growth);
    }
    const double growth_seconds = monotonic_seconds()-growth_started;

    EventCounters warmup = {0};
    while (warmup.events < warmup_events) (void)dispatch_one(&warmup);

    EventCounters measured = {0};
    const double simulation_start = t;
    const double measured_started = monotonic_seconds();
    while (measured.events < timed_events) (void)dispatch_one(&measured);
    const double measured_seconds = monotonic_seconds()-measured_started;
    const double simulation_span = t-simulation_start;
    physicalQ();

    printf("BENCH_HEADLESS mode=%s repetition=%d seed=%lu N=%d phi=%.9g "
           "cells=%dx%d setup_s=%.9g growth_s=%.9g growth_events=%lu "
           "warmup_events=%lu warmup_collisions=%lu events=%lu "
           "wall_s=%.9g events_per_s=%.9g collisions=%lu "
           "collisions_per_s=%.9g stale_collisions=%lu crossings=%lu "
           "sim_dt=%.9g sim_time_per_s=%.9g E_per_N=%.9g\n",
           mode_name(mode), repetition, seed, N, phi, Nxcells, Nycells,
           setup_seconds, growth_seconds, growth.events,
           warmup.events, warmup.physical_collisions, measured.events,
           measured_seconds, measured.events/measured_seconds,
           measured.physical_collisions,
           measured.physical_collisions/measured_seconds,
           measured.stale_collisions, measured.crossings,
           simulation_span, simulation_span/measured_seconds, E/N);
    fflush(stdout);
    freeArrays();
}

static void usage(const char *program)
{
    fprintf(stderr,
            "usage: %s MODE N PHI EVENTS WARMUP_EVENTS FIRST_SEED [REPETITIONS]\n"
            "MODE is ballistic, damping, delta, or damping-delta\n",
            program);
}

int main(int argc, char **argv)
{
    if (argc < 7 || argc > 8) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    BenchmarkMode mode;
    if (!parse_mode(argv[1], &mode)) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    const int particle_count = atoi(argv[2]);
    const double packing_fraction = atof(argv[3]);
    const unsigned long timed_events = strtoul(argv[4], NULL, 10);
    const unsigned long warmup_events = strtoul(argv[5], NULL, 10);
    const unsigned long first_seed = strtoul(argv[6], NULL, 10);
    const int repetitions = argc == 8 ? atoi(argv[7]) : 1;
    if (particle_count < 4 || !isfinite(packing_fraction) ||
        packing_fraction <= 0.0 || packing_fraction >= 0.90 ||
        timed_events == 0 || repetitions < 1) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    for (int repetition = 0; repetition < repetitions; ++repetition) {
        run_once(mode, particle_count, packing_fraction, timed_events,
                 warmup_events, first_seed+(unsigned long)repetition,
                 repetition+1);
    }
    return EXIT_SUCCESS;
}
