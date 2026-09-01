#include "EDMD.h"
#include "polygon_particles.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

extern int N, load, Hex, S1, version, periodicDumps, ther;
extern int automaticGrowthRate, Nxcells, Nycells;
extern int polygonCollisionCellNeighborRange;
extern double phi, Einit, t, vr, E, res, resW, resT, Lx, Ly;
extern unsigned long int ncol, ncross;
extern particle *particles;
extern node *nextEvent;
extern void (*doTheCollision)(void);
extern void (*doTheWall)(void);
extern void init_genrand(unsigned long seed);

typedef struct PhaseCounters {
    unsigned long events;
    unsigned long event_type[16];
    unsigned long physical_collisions;
    unsigned long maximum_same_time_streak;
    unsigned long same_time_streak;
    double last_deadline;
} PhaseCounters;

static double monotonic_seconds(void)
{
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return NAN;
    return (double)now.tv_sec + 1.0e-9*(double)now.tv_nsec;
}

static void benchmark_failure(const char *message)
{
    fprintf(stderr,
        "polygon benchmark failure: %s (t=%.17g phi=%.17g E/N=%.17g ncol=%lu)\n",
        message, t, simulationCurrentPackingFraction(t),
        N > 0 ? E/N : NAN, ncol);
    if (simulationHasFault())
        fprintf(stderr, "simulation fault: %s\n", simulationFaultMessage());
    exit(EXIT_FAILURE);
}

static void configure_triangle(const char *name)
{
    SimParticleShapeConfig shape;
    memset(&shape, 0, sizeof(shape));
    shape.mode = SIM_PARTICLE_CONVEX_POLYGONS;
    shape.vertex_count = 3;
    if (strcmp(name, "skinny") == 0) {
        shape.vertices[0][0] = -1.0;
        shape.vertices[0][1] = -0.3;
        shape.vertices[1][0] = 1.0;
        shape.vertices[1][1] = -0.3;
        shape.vertices[2][0] = 0.0;
        shape.vertices[2][1] = 0.76;
    } else {
        for (int vertex = 0; vertex < 3; ++vertex) {
            const double angle = 0.5*M_PI + 2.0*M_PI*vertex/3.0;
            shape.vertices[vertex][0] = cos(angle);
            shape.vertices[vertex][1] = sin(angle);
        }
    }
    shape.random_orientation = true;
    shape.angular_temperature = 1.0;
    polygon_particles_set_config(&shape);
    if (!polygon_particles_requested())
        benchmark_failure("triangle configuration was rejected");
}

static void count_deadline(PhaseCounters *counters, double deadline)
{
    if (deadline == counters->last_deadline) {
        counters->same_time_streak++;
    } else {
        counters->last_deadline = deadline;
        counters->same_time_streak = 1;
    }
    if (counters->same_time_streak > counters->maximum_same_time_streak)
        counters->maximum_same_time_streak = counters->same_time_streak;
}

static enum event dispatch_event(PhaseCounters *counters)
{
    nextEvent = findNextEvent();
    if (nextEvent == NULL || !isfinite(nextEvent->t) || nextEvent->t < t)
        benchmark_failure("calendar returned an invalid event");
    const enum event type = nextEvent->type;
    const unsigned long collisions_before = ncol;
    count_deadline(counters, nextEvent->t);
    if ((unsigned)type < sizeof(counters->event_type)/
                         sizeof(counters->event_type[0]))
        counters->event_type[type]++;
    counters->events++;
    t = nextEvent->t;
    removeEventFromQueue(nextEvent);
    switch (type) {
        case COLLISION:
            doTheCollision();
            break;
        case WALL:
            doTheWall();
            break;
        case CELLCROSS:
            ++ncross;
            doTheCrossing();
            break;
        case POLYGON_CELLCROSS:
            ++ncross;
            doTheCrossingPolygon();
            break;
        case POLYGON_RECHECK:
            polygon_particles_do_recheck();
            break;
        default:
            benchmark_failure("unexpected event in timed loop");
    }
    if (simulationHasFault())
        benchmark_failure("event processing raised a fault");
    counters->physical_collisions += ncol-collisions_before;
    return type;
}

static void print_phase(const char *name, const PhaseCounters *counters,
                        double seconds, double simulation_span)
{
    const double rate = seconds > 0.0 ? counters->events/seconds : NAN;
    const double collision_rate = seconds > 0.0
        ? counters->physical_collisions/seconds : NAN;
    const double simulation_rate = seconds > 0.0
        ? simulation_span/seconds : NAN;
    const unsigned long collision_events =
        counters->event_type[COLLISION];
    printf("BENCH phase=%s wall_s=%.9g sim_dt=%.9g events=%lu "
           "events_per_s=%.9g sim_time_per_s=%.9g "
           "collision_events=%lu physical_collisions=%lu "
           "physical_collisions_per_s=%.9g stale_collision_events=%lu "
           "crossings=%lu rechecks=%lu walls=%lu max_same_time=%lu\n",
           name, seconds, simulation_span, counters->events, rate,
           simulation_rate, collision_events, counters->physical_collisions,
           collision_rate,
           collision_events >= counters->physical_collisions
               ? collision_events-counters->physical_collisions : 0UL,
           counters->event_type[CELLCROSS] +
               counters->event_type[POLYGON_CELLCROSS],
           counters->event_type[POLYGON_RECHECK],
           counters->event_type[WALL],
           counters->maximum_same_time_streak);
}

static void print_scheduler_phase(
    const char *name, const PolygonParticleSchedulerStats *stats)
{
    printf("BENCH scheduler phase=%s calls=%llu candidates=%llu "
           "broadphase_rejections=%llu ccd_calls=%llu ccd_iterations=%llu "
           "wall_candidates=%llu wall_broadphase_rejections=%llu "
           "wall_manifolds=%llu rechecks_scheduled=%llu\n",
           name,
           stats->collision_event_calls,
           stats->pair_candidates,
           stats->pair_broadphase_rejections,
           stats->pair_ccd_calls,
           stats->pair_ccd_iterations,
           stats->wall_candidates,
           stats->wall_broadphase_rejections,
           stats->wall_manifold_evaluations,
           stats->rechecks_scheduled);
}

static const char *event_name(enum event type)
{
    switch (type) {
        case COLLISION: return "collision";
        case WALL: return "wall";
        case CELLCROSS: return "cell-crossing";
        case POLYGON_CELLCROSS: return "polygon-cell-crossing";
        case POLYGON_RECHECK: return "polygon-recheck";
        case GROWSTOP: return "growth-stop";
        default: return "other";
    }
}

static double verify_at(double query_time)
{
    char message[512];
    const double started = monotonic_seconds();
    const double tolerance = 2.0e-9*(1.0 + fabs(query_time));
    if (!polygon_particles_verify_no_overlaps(
            query_time, tolerance, message, sizeof(message))) {
        fprintf(stderr, "overlap verification: %s\n", message);
        benchmark_failure("one-shot SAT verification failed");
    }
    return monotonic_seconds()-started;
}

int main(int argc, char **argv)
{
    const char *shape_name = argc > 1 ? argv[1] : "regular";
    const int particle_count = argc > 2 ? atoi(argv[2]) : 1000;
    const double requested_phi = argc > 3 ? atof(argv[3]) : 0.80;
    const unsigned long production_event_limit = argc > 4
        ? strtoul(argv[4], NULL, 10) : 250000UL;
    const unsigned long seed = argc > 5
        ? strtoul(argv[5], NULL, 10) : 1UL;
    const int verify = argc > 6 ? atoi(argv[6]) != 0 : 1;
    const double manual_growth_rate = argc > 7 ? atof(argv[7]) : 0.0;
    const unsigned long checkpoint_interval = argc > 8
        ? strtoul(argv[8], NULL, 10) : 0UL;
    const unsigned long checkpoint_start = argc > 9
        ? strtoul(argv[9], NULL, 10) : checkpoint_interval;
    if (particle_count < 4 || !(requested_phi > 0.0) ||
        requested_phi >= 1.0 || production_event_limit == 0)
        benchmark_failure("invalid command-line arguments");

    configure_triangle(shape_name);
    N = particle_count;
    phi = requested_phi;
    Einit = 1.0;
    version = (int)seed;
    load = Hex = S1 = 0;
    periodicDumps = ther = 0;
    res = resW = resT = 1.0;
    if (manual_growth_rate > 0.0) {
        automaticGrowthRate = 0;
        vr = manual_growth_rate;
    }
    simulationClearFault();
    constantInit(0, NULL);
    init_genrand(seed);

    const double setup_started = monotonic_seconds();
    particlesInit();
    adaptCellGridToParticles();
    cellListInit();
    physicalQ();
    eventListInit();
    const double setup_seconds = monotonic_seconds()-setup_started;
    if (!polygon_particles_growing())
        benchmark_failure("initializer did not enter polygon-native growth");

    PhaseCounters growth;
    memset(&growth, 0, sizeof(growth));
    growth.last_deadline = -DBL_MAX;
    polygon_particles_reset_scheduler_stats();
    const double growth_start_time = t;
    const double growth_started = monotonic_seconds();
    while (1) {
        nextEvent = findNextEvent();
        if (nextEvent == NULL || !isfinite(nextEvent->t) || nextEvent->t < t)
            benchmark_failure("invalid event while awaiting GROWSTOP");
        if (nextEvent->type == GROWSTOP) break;
        (void)dispatch_event(&growth);
    }
    const double growth_seconds = monotonic_seconds()-growth_started;
    PolygonParticleSchedulerStats growth_scheduler;
    polygon_particles_get_scheduler_stats(&growth_scheduler);
    const double growth_stop_time = nextEvent->t;
    physicalQ();
    const double temporary_energy_per_particle = E/N;
    const double realized_before_handoff =
        simulationCurrentPackingFraction(growth_stop_time);
    const double growth_verify_seconds = verify
        ? verify_at(growth_stop_time) : 0.0;

    const double handoff_started = monotonic_seconds();
    t = nextEvent->t;
    removeEventFromQueue(nextEvent);
    stopGrow();
    const double handoff_seconds = monotonic_seconds()-handoff_started;
    if (simulationHasFault() || polygon_particles_growing())
        benchmark_failure("GROWSTOP handoff failed");
    physicalQ();
    const double production_initial_energy = E;
    const double production_initial_rotation =
        polygon_particles_rotational_energy();
    const double production_translation = E-production_initial_rotation;

    PhaseCounters production;
    memset(&production, 0, sizeof(production));
    production.last_deadline = -DBL_MAX;
    polygon_particles_reset_scheduler_stats();
    const double production_start_time = t;
    double production_seconds = 0.0;
    double diagnostic_seconds = 0.0;
    unsigned long next_checkpoint = checkpoint_interval > 0
        ? checkpoint_start : production_event_limit;
    if (next_checkpoint == 0) next_checkpoint = checkpoint_interval;
    enum event last_event_type = GROWSTOP;
    while (production.events < production_event_limit) {
        unsigned long chunk_end = production_event_limit;
        if (checkpoint_interval > 0 && next_checkpoint < chunk_end)
            chunk_end = next_checkpoint;
        const double chunk_started = monotonic_seconds();
        while (production.events < chunk_end)
            last_event_type = dispatch_event(&production);
        production_seconds += monotonic_seconds()-chunk_started;
        if (checkpoint_interval > 0 &&
            production.events == next_checkpoint &&
            production.events < production_event_limit) {
            const double diagnostic_started = monotonic_seconds();
            physicalQ();
            printf("BENCH checkpoint events=%lu sim_t=%.12g E_per_N=%.12g "
                   "energy_rel_drift=%.12g last_event=%s "
                   "physical_collisions=%lu stale_collision_events=%lu\n",
                   production.events, t, E/N,
                   fabs(E-production_initial_energy)/
                       fmax(1.0, production_initial_energy),
                   event_name(last_event_type),
                   production.physical_collisions,
                   production.event_type[COLLISION] >=
                           production.physical_collisions
                       ? production.event_type[COLLISION]-
                           production.physical_collisions : 0UL);
            diagnostic_seconds += monotonic_seconds()-diagnostic_started;
            if (production_event_limit-next_checkpoint < checkpoint_interval)
                next_checkpoint = production_event_limit;
            else
                next_checkpoint += checkpoint_interval;
        }
    }
    const double production_end_time = t;
    PolygonParticleSchedulerStats production_scheduler;
    polygon_particles_get_scheduler_stats(&production_scheduler);

    physicalQ();
    const double production_final_energy = E;
    const double production_verify_seconds = verify ? verify_at(t) : 0.0;
    const double relative_energy_drift = fabs(
        production_final_energy-production_initial_energy)/
        fmax(1.0, production_initial_energy);
    const double final_phi = simulationCurrentPackingFraction(t);
    if (!isfinite(relative_energy_drift) || relative_energy_drift > 1.0e-9)
        benchmark_failure("production energy conservation regressed");
    if (!isfinite(final_phi) ||
        fabs(final_phi-requested_phi) > 1.0e-10)
        benchmark_failure("realized packing fraction changed in production");

    printf("BENCH config shape=%s N=%d target_phi=%.9g seed=%lu "
           "vr=%.9g Lx=%.9g Ly=%.9g cells=%dx%d neighbor_range=%d "
           "setup_s=%.9g verification=%d checkpoint_interval=%lu "
           "checkpoint_start=%lu\n",
           shape_name, N, requested_phi, seed, vr, Lx, Ly,
           Nxcells, Nycells, polygonCollisionCellNeighborRange,
           setup_seconds, verify, checkpoint_interval, checkpoint_start);
    print_phase("growth", &growth, growth_seconds,
                growth_stop_time-growth_start_time);
    print_scheduler_phase("growth", &growth_scheduler);
    printf("BENCH growth_state realized_phi=%.12g temporary_E_per_N=%.9g "
           "verify_s=%.9g\n",
           realized_before_handoff, temporary_energy_per_particle,
           growth_verify_seconds);
    printf("BENCH phase=handoff wall_s=%.9g translational_E_per_N=%.9g "
           "rotational_E_per_N=%.9g\n",
           handoff_seconds, production_translation/N,
           production_initial_rotation/N);
    print_phase("production", &production, production_seconds,
                production_end_time-production_start_time);
    print_scheduler_phase("production", &production_scheduler);
    printf("BENCH production_state final_phi=%.12g E_initial_per_N=%.9g "
           "E_final_per_N=%.9g energy_rel_drift=%.9g verify_s=%.9g "
           "diagnostics_s=%.9g\n",
           final_phi, production_initial_energy/N,
           production_final_energy/N, relative_energy_drift,
           production_verify_seconds, diagnostic_seconds);

    freeArrays();
    return EXIT_SUCCESS;
}
