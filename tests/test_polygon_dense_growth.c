#include "EDMD.h"
#include "polygon_particles.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

extern int N, load, Hex, S1, version, periodicDumps, ther;
extern int automaticGrowthRate;
extern int polygonCollisionCellNeighborRange;
extern double phi, Einit, t, vr, E, res, resW, resT;
extern unsigned long int ncol, ncross;
extern particle *particles;
extern node *nextEvent;
extern void (*doTheCollision)(void);
extern void (*doTheWall)(void);
extern void init_genrand(unsigned long seed);

static void fail(const char *message)
{
    fprintf(stderr, "polygon dense-growth failure: %s\n", message);
    fprintf(stderr, "state: t=%.17g packing=%.17g E/N=%.17g collisions=%lu\n",
            t, simulationCurrentPackingFraction(t),
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
        /* Deliberately lower area/circumradius than the regular triangle. Its
           old circumscribed-disk initializer cap was close to the value in
           the reported GUI failure. */
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
    if (!polygon_particles_requested()) fail("triangle configuration rejected");
}

static double translational_energy(void)
{
    double energy = 0.0;
    for (int index = 0; index < N; ++index)
        energy += 0.5*particles[index].m*
            (particles[index].vx*particles[index].vx +
             particles[index].vy*particles[index].vy);
    return energy;
}

static void verify_no_overlap(const char *phase)
{
    char message[512];
    const double tolerance = 2.0e-9*(1.0 + fabs(t));
    if (!polygon_particles_verify_no_overlaps(
            t, tolerance, message, sizeof(message))) {
        fprintf(stderr, "%s: %s\n", phase, message);
        fail("SAT overlap invariant failed");
    }
}

static enum event dispatch_next_event(int *growth_stopped)
{
    nextEvent = findNextEvent();
    if (nextEvent == NULL || !isfinite(nextEvent->t) || nextEvent->t < t) {
        if (nextEvent != NULL)
            fprintf(stderr, "invalid deadline: now=%.17g next=%.17g type=%d i=%d j=%d\n",
                    t, nextEvent->t, (int)nextEvent->type,
                    nextEvent->i, nextEvent->j);
        fail("calendar returned an invalid event deadline");
    }
    const enum event type = nextEvent->type;
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
        case GROWSTOP:
            stopGrow();
            if (growth_stopped != NULL) *growth_stopped = 1;
            break;
        default:
            fail("unexpected auxiliary event in isolated growth test");
    }
    if (simulationHasFault()) fail("event processing raised a fault");
    return type;
}

int main(int argc, char **argv)
{
    const char *shape_name = argc > 1 ? argv[1] : "regular";
    const int count = argc > 2 ? atoi(argv[2]) : 64;
    const double requested_phi = argc > 3 ? atof(argv[3]) : 0.60;
    const unsigned long growth_event_budget = argc > 4
        ? strtoul(argv[4], NULL, 10) : 1000000UL;
    const unsigned long production_event_budget = argc > 5
        ? strtoul(argv[5], NULL, 10) : 20000UL;
    const double manual_growth_rate = argc > 6 ? atof(argv[6]) : 0.0;
    if (count < 4 || !(requested_phi > 0.0) || requested_phi >= 1.0)
        fail("invalid test arguments");

    configure_triangle(shape_name);
    N = count;
    phi = requested_phi;
    Einit = 1.0;
    version = 1;
    load = Hex = S1 = 0;
    periodicDumps = ther = 0;
    res = resW = resT = 1.0;
    if (manual_growth_rate > 0.0) {
        automaticGrowthRate = 0;
        vr = manual_growth_rate;
    }
    simulationClearFault();
    constantInit(0, NULL);
    init_genrand(1);
    particlesInit();
    adaptCellGridToParticles();
    int expected_neighbor_range = (int)floor(16.0*requested_phi + 0.5);
    if (expected_neighbor_range < 2) expected_neighbor_range = 2;
    if (expected_neighbor_range >
        POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE)
        expected_neighbor_range =
            POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE;
    if (polygonCollisionCellNeighborRange != expected_neighbor_range)
        fail("adaptive polygon cell range differs from the density rule");
    cellListInit();
    physicalQ();
    eventListInit();
    if (!polygon_particles_growing())
        fail("random polygon initializer did not enter native growth");

    const double start_energy = E;
    double maximum_growth_energy = E;
    unsigned long growth_events = 0;
    int growth_stopped = 0;
    while (!growth_stopped && growth_events < growth_event_budget) {
        (void)dispatch_next_event(&growth_stopped);
        ++growth_events;
        physicalQ();
        if (!isfinite(E)) fail("non-finite preparation energy");
        if (E > maximum_growth_energy) maximum_growth_energy = E;
        /* The runtime aggregate guard is deliberately generous enough to
           finish dense preparations. This tighter regression bound still has
           ample headroom above the observed native-growth trajectory while
           making the historical 1e9 circumscribed-disk blow-up fail. */
        if (E/N > 140.0)
            fail("preparation energy escaped its configured bound");
    }
    if (!growth_stopped)
        fail("target packing fraction was not reached within the event budget");
    if (polygon_particles_growing())
        fail("growth flag survived GROWSTOP");

    verify_no_overlap("handoff");
    const double realized_phi = simulationCurrentPackingFraction(t);
    if (!isfinite(realized_phi) ||
        fabs(realized_phi-requested_phi) >
            2.0e-10*fmax(1.0, requested_phi))
        fail("realized area fraction differs from the requested fraction");
    for (int index = 0; index < N; ++index) {
        if (fabs(particles[index].rad-1.0) > 2.0e-12 ||
            particles[index].vr != 0.0)
            fail("GROWSTOP did not land exactly on the final particle scale");
    }

    physicalQ();
    const double translation = translational_energy();
    const double rotation = polygon_particles_rotational_energy();
    if (!isfinite(translation) || !isfinite(rotation) ||
        fabs(translation/N-Einit) > 2.0e-11)
        fail("post-growth translational energy does not match Einit");
    if (!(rotation/N > 0.05) || rotation/N > 2.5)
        fail("post-growth angular thermalization is invalid");
    if (fabs(E-translation-rotation) > 2.0e-12*fmax(1.0, E))
        fail("reported total energy omits a rigid-body component");

    const double production_energy = E;
    const double production_start_time = t;
    double maximum_relative_drift = 0.0;
    unsigned long production_events = 0;
    unsigned long production_collisions = 0;
    double last_deadline = -DBL_MAX;
    unsigned long same_time_streak = 0;
    while (production_events < production_event_budget) {
        nextEvent = findNextEvent();
        if (nextEvent == NULL || !isfinite(nextEvent->t) || nextEvent->t < t)
            fail("production calendar returned an invalid deadline");
        if (nextEvent->t == last_deadline) {
            ++same_time_streak;
        } else {
            last_deadline = nextEvent->t;
            same_time_streak = 1;
        }
        if (same_time_streak > 4096)
            fail("production entered a same-time event livelock");
        const unsigned long collisions_before = ncol;
        (void)dispatch_next_event(NULL);
        production_collisions += ncol-collisions_before;
        ++production_events;
        physicalQ();
        if (!isfinite(E)) fail("non-finite production energy");
        const double relative_drift = fabs(E-production_energy)/
            fmax(1.0, production_energy);
        if (relative_drift > maximum_relative_drift)
            maximum_relative_drift = relative_drift;
        if (relative_drift > 2.0e-8)
            fail("elastic polygon dynamics changed total energy");
        if ((production_events & 255UL) == 0)
            verify_no_overlap("production");
    }
    verify_no_overlap("final");
    if (!(t > production_start_time) || production_collisions == 0)
        fail("production run did not make physical progress");

    printf("polygon dense growth ok shape=%s N=%d phi=%.9g vr=%.9g "
           "growth_events=%lu t_stop=%.9g E_start/N=%.9g "
           "E_growth_max/N=%.9g E_translation/N=%.9g E_rotation/N=%.9g "
           "production_events=%lu production_collisions=%lu t_final=%.9g "
           "energy_max_rel_drift=%.3g\n",
           shape_name, N, requested_phi, vr, growth_events,
           production_start_time, start_energy/N, maximum_growth_energy/N,
           translation/N, rotation/N, production_events,
           production_collisions, t, maximum_relative_drift);
    freeArrays();
    return EXIT_SUCCESS;
}
