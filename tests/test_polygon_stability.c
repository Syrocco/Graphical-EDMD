#include "EDMD.h"
#include "convex_polygon.h"
#include "polygon_particles.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#define EDMD_TEST_BOUNDARY_PERIODIC 0
#define EDMD_TEST_BOUNDARY_PLANAR 1
#define EDMD_TEST_BOUNDARY_CIRCULAR 2

#ifndef EDMD_TEST_BOUNDARY_MODE
#error "define EDMD_TEST_BOUNDARY_MODE for the polygon stability fixture"
#endif
#if EDMD_TEST_BOUNDARY_MODE < EDMD_TEST_BOUNDARY_PERIODIC || \
    EDMD_TEST_BOUNDARY_MODE > EDMD_TEST_BOUNDARY_CIRCULAR
#error "invalid EDMD_TEST_BOUNDARY_MODE"
#endif

extern int N, Nbig, Nsmall;
extern int periodicDumps, ther;
extern int addWallx, addWally, addCircularWall;
extern double Lx, Ly, halfLx, halfLy, t, vr;
extern double res, resW, resT;
extern particle *particles;
extern node *nextEvent;
extern unsigned long int ncol, ncross;
extern void init_genrand(unsigned long seed);

typedef struct StabilityRng {
    uint64_t state;
} StabilityRng;

static int diagnostic_event_i = -1;
static int diagnostic_event_j = -1;
static const char *diagnostic_phase = "initial";
static const char *diagnostic_shape = "regular-pentagon";

typedef struct StabilityMinimum {
    double value;
    double time;
    unsigned long seed;
    int dissipative;
    int event_number;
    int first;
    int second;
} StabilityMinimum;

static StabilityMinimum worst_pair;
static StabilityMinimum worst_wall;
static unsigned long stability_event_counts[16];
static unsigned long stability_recheck_owners[64];
static double stability_last_event_time;
static unsigned long stability_same_time_events;
static unsigned long stability_max_same_time_events;
static double stability_initial_energy;
static double stability_minimum_energy;
static double stability_maximum_energy;
static double stability_maximum_relative_energy_error;
static double stability_previous_energy;
static double stability_maximum_step_energy_increase;

static double exact_system_energy(void)
{
    double energy = polygon_particles_rotational_energy();
    for (int index = 0; index < N; ++index) {
        const particlePrediction prediction =
            predictFreeFlyNormal(particles + index, t);
        energy += 0.5*particles[index].m*
            (prediction.vx*prediction.vx + prediction.vy*prediction.vy);
    }
    return energy;
}

static void reset_stability_minima(void)
{
    memset(&worst_pair, 0, sizeof(worst_pair));
    memset(&worst_wall, 0, sizeof(worst_wall));
    worst_pair.first = worst_pair.second = -1;
    worst_wall.first = worst_wall.second = -1;
    memset(stability_event_counts, 0, sizeof(stability_event_counts));
    memset(stability_recheck_owners, 0, sizeof(stability_recheck_owners));
    stability_last_event_time = -DBL_MAX;
    stability_same_time_events = 0;
    stability_max_same_time_events = 0;
    stability_initial_energy = exact_system_energy();
    stability_minimum_energy = stability_initial_energy;
    stability_maximum_energy = stability_initial_energy;
    stability_maximum_relative_energy_error = 0.0;
    stability_previous_energy = stability_initial_energy;
    stability_maximum_step_energy_increase = 0.0;
}

static int check_energy_invariant(unsigned long seed, int dissipative,
                                  int event_number, enum event event_type,
                                  const char *phase)
{
    const double energy = exact_system_energy();
    if (!isfinite(energy) || !(stability_initial_energy > 0.0)) {
        fprintf(stderr,
            "polygon stability invalid energy shape=%s seed=%lu dissipative=%d event=%d type=%d phase=%s t=%.17g energy=%.17g initial=%.17g\n",
            diagnostic_shape, seed, dissipative, event_number,
            (int)event_type, phase, t, energy, stability_initial_energy);
        return 0;
    }
    if (energy < stability_minimum_energy) stability_minimum_energy = energy;
    if (energy > stability_maximum_energy) stability_maximum_energy = energy;
    const double relative_error =
        fabs(energy-stability_initial_energy)/stability_initial_energy;
    const double step_increase = energy-stability_previous_energy;
    if (step_increase > stability_maximum_step_energy_increase)
        stability_maximum_step_energy_increase = step_increase;
    if (relative_error > stability_maximum_relative_energy_error)
        stability_maximum_relative_energy_error = relative_error;

    /* An elastic event must conserve the complete rigid-body energy. A
       dissipative event may remove it, but can never exceed its initial
       value. The tolerance covers only accumulated double rounding. */
    const double upper_tolerance = 5.0e-9*
        (1.0 + (double)event_number/10000.0);
    const double step_tolerance = 5.0e-10*
        fmax(1.0, stability_initial_energy);
    const bool failed = dissipative
        ? step_increase > step_tolerance
        : relative_error > upper_tolerance;
    if (failed) {
        fprintf(stderr,
            "polygon stability energy failure shape=%s seed=%lu dissipative=%d event=%d type=%d owner=%d,%d phase=%s t=%.17g energy=%.17g initial=%.17g relative_error=%.17g min=%.17g max=%.17g tolerance=%.17g\n",
            diagnostic_shape, seed, dissipative, event_number,
            (int)event_type, diagnostic_event_i, diagnostic_event_j,
            phase, t, energy, stability_initial_energy, relative_error,
            stability_minimum_energy, stability_maximum_energy,
            upper_tolerance);
        return 0;
    }
    stability_previous_energy = energy;
    return 1;
}

static void update_stability_minimum(StabilityMinimum *minimum,
                                     double value, unsigned long seed,
                                     int dissipative, int event_number,
                                     int first, int second)
{
    if (minimum == NULL || !(value < minimum->value)) return;
    minimum->value = value;
    minimum->time = t;
    minimum->seed = seed;
    minimum->dissipative = dissipative;
    minimum->event_number = event_number;
    minimum->first = first;
    minimum->second = second;
}

static uint64_t stability_random_u64(StabilityRng *rng)
{
    uint64_t value = (rng->state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30))*UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27))*UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static double stability_uniform(StabilityRng *rng, double minimum,
                                double maximum)
{
    const double unit = (double)(stability_random_u64(rng) >> 11)
        * (1.0/9007199254740992.0);
    return minimum + (maximum-minimum)*unit;
}

static int configure_shape(int shape_case, double angular_temperature)
{
    SimParticleShapeConfig shape;
    memset(&shape, 0, sizeof(shape));
    shape.mode = SIM_PARTICLE_CONVEX_POLYGONS;
    if (shape_case == 0) {
        diagnostic_shape = "regular-pentagon";
        shape.vertex_count = 5;
        for (int vertex = 0; vertex < shape.vertex_count; ++vertex) {
            const double angle = 0.37 + 2.0*M_PI*vertex/shape.vertex_count;
            shape.vertices[vertex][0] = cos(angle);
            shape.vertices[vertex][1] = sin(angle);
        }
    } else if (shape_case == 1) {
        diagnostic_shape = "asymmetric-triangle";
        shape.vertex_count = 3;
        shape.vertices[0][0] = -0.95; shape.vertices[0][1] = -0.45;
        shape.vertices[1][0] =  1.00; shape.vertices[1][1] = -0.30;
        shape.vertices[2][0] = -0.15; shape.vertices[2][1] =  1.00;
    } else if (shape_case == 2) {
        /* Three nearby points on a strictly convex affine circle exercise
           short, almost-parallel adjacent support edges without violating
           the editor/runtime strict-convexity limits. */
        static const double angles[6] = {
            0.0, 0.003, 0.006, 1.6, 3.1, 4.7
        };
        diagnostic_shape = "skew-near-parallel-hexagon";
        shape.vertex_count = 6;
        for (int vertex = 0; vertex < shape.vertex_count; ++vertex) {
            const double cosine = cos(angles[vertex]);
            const double sine = sin(angles[vertex]);
            shape.vertices[vertex][0] = cosine + 0.35*sine;
            shape.vertices[vertex][1] = 0.65*sine;
        }
    } else if (shape_case == 3) {
        diagnostic_shape = "dense-aligned-square";
        shape.vertex_count = 4;
        shape.vertices[0][0] = -1.0; shape.vertices[0][1] = -1.0;
        shape.vertices[1][0] =  1.0; shape.vertices[1][1] = -1.0;
        shape.vertices[2][0] =  1.0; shape.vertices[2][1] =  1.0;
        shape.vertices[3][0] = -1.0; shape.vertices[3][1] =  1.0;
        shape.random_orientation = false;
    } else {
        return 0;
    }
    if (shape_case != 3) shape.random_orientation = true;
    shape.angular_temperature = angular_temperature;
    polygon_particles_set_config(&shape);
    SimParticleShapeConfig accepted;
    polygon_particles_get_config(&accepted);
    return accepted.mode == SIM_PARTICLE_CONVEX_POLYGONS &&
           accepted.vertex_count == shape.vertex_count;
}

static int append_position(double *positions, int count, double x, double y)
{
    positions[2*count] = x;
    positions[2*count+1] = y;
    return count + 1;
}

static int make_initial_positions(double *positions, int capacity)
{
    int count = 0;
#if EDMD_TEST_BOUNDARY_MODE == EDMD_TEST_BOUNDARY_CIRCULAR
    const double center = 5.0;
    for (int row = -3; row <= 3; ++row) {
        for (int column = -3; column <= 3; ++column) {
            const double x = 1.5*column + ((row & 1) ? 0.75 : 0.0);
            const double y = 1.30*row;
            if (hypot(x, y) <= 4.18 && count < capacity)
                count = append_position(positions, count,
                                        center+x, center+y);
        }
    }
#else
    for (int row = 0; row < 6; ++row) {
        for (int column = 0; column < 6; ++column) {
            if (count >= capacity) break;
            count = append_position(positions, count,
                0.82 + 1.66*column + ((row & 1) ? 0.05 : 0.0),
                0.82 + 1.66*row);
        }
    }
#endif
    return count;
}

static int initialize_stability_case(unsigned long seed, int dissipative,
                                     int shape_case)
{
    double positions[2*64];
    StabilityRng rng = {UINT64_C(0xd1b54a32d192ed03) ^ seed};
    N = make_initial_positions(positions, 64);
    if (N < 8) return 0;
    Nbig = N;
    Nsmall = 0;
    Lx = 10.0;
    Ly = 10.0;
    t = 0.0;
    vr = 1.0;
    ncol = 0;
    ncross = 0;
    periodicDumps = 0;
    ther = 0;
    res = dissipative ? 0.72 : 1.0;
    resW = dissipative ? 0.68 : 1.0;
    resT = dissipative ? 0.83 : 1.0;
    simulationClearFault();
    init_genrand(seed);
    if (!configure_shape(shape_case, 3.0)) return 0;

    particles = calloc((size_t)N, sizeof(*particles));
    if (particles == NULL) return 0;
    double mean_vx = 0.0;
    double mean_vy = 0.0;
    for (int index = 0; index < N; ++index) {
        particle *value = particles + index;
        value->num = index;
        value->type = 1;
        value->m = 0.7 + 0.3*(index%3);
#if EDMD_TEST_BOUNDARY_MODE != EDMD_TEST_BOUNDARY_CIRCULAR
        value->rad = shape_case == 3 ? 1.0
            : (index%3 == 0 ? 0.68 : (index%3 == 1 ? 0.56 : 0.62));
#else
        value->rad = index%3 == 0 ? 0.68 : (index%3 == 1 ? 0.56 : 0.62);
#endif
        value->x = positions[2*index];
        value->y = positions[2*index+1];
        value->vx = stability_uniform(&rng, -1.8, 1.8);
        value->vy = stability_uniform(&rng, -1.8, 1.8);
        value->t = t;
        value->lastColl = t;
        mean_vx += value->vx;
        mean_vy += value->vy;
    }
#if EDMD_TEST_BOUNDARY_MODE == EDMD_TEST_BOUNDARY_PERIODIC
    mean_vx /= N;
    mean_vy /= N;
    for (int index = 0; index < N; ++index) {
        particles[index].vx -= mean_vx;
        particles[index].vy -= mean_vy;
    }
#else
    (void)mean_vx;
    (void)mean_vy;
#endif

    adaptCellGridToParticles();
    for (int index = 0; index < N; ++index) {
        particles[index].cell[0] = coordToCell(particles[index].x, 1);
        particles[index].cell[1] = coordToCell(particles[index].y, 0);
    }
    cellListInit();
    eventListInit();
    return polygon_particles_active() && !simulationHasFault();
}

static int current_shape(ConvexPolygon *shape)
{
    memset(shape, 0, sizeof(*shape));
    shape->count = polygon_particles_vertex_count();
    if (shape->count < 3 || shape->count > CONVEX_POLYGON_MAX_VERTICES)
        return 0;
    for (int vertex = 0; vertex < shape->count; ++vertex) {
        if (!polygon_particles_vertex(vertex, &shape->vertices[vertex].x,
                                      &shape->vertices[vertex].y))
            return 0;
    }
    return convex_polygon_refresh(shape);
}

static ConvexPolygonPose particle_pose_now(int index)
{
    const particlePrediction prediction =
        predictFreeFlyNormal(particles + index, t);
    ConvexPolygonPose pose;
    pose.position.x = prediction.x;
    pose.position.y = prediction.y;
    pose.angle = polygon_particles_angle_at(index, t);
    pose.scale = prediction.rad;
    return pose;
}

static double vertex_wall_gap(const ConvexPolygonPose *pose,
                              ConvexPolygonVec2 local)
{
    const double cosine = cos(pose->angle);
    const double sine = sin(pose->angle);
    const double x = pose->position.x + pose->scale*
        (cosine*local.x-sine*local.y);
    const double y = pose->position.y + pose->scale*
        (sine*local.x+cosine*local.y);
    (void)x;
    (void)y;
    double gap = DBL_MAX;
#if EDMD_TEST_BOUNDARY_MODE == EDMD_TEST_BOUNDARY_PLANAR
    gap = fmin(gap, fmin(x, Lx-x));
    gap = fmin(gap, fmin(y, Ly-y));
#endif
#if EDMD_TEST_BOUNDARY_MODE == EDMD_TEST_BOUNDARY_CIRCULAR
    gap = fmin(gap, halfLx-hypot(x-halfLx, y-halfLy));
#endif
    return gap;
}

static int exact_invariants(unsigned long seed, int dissipative,
                            int event_number, enum event event_type)
{
    ConvexPolygon shape;
    if (!current_shape(&shape)) return 0;
    for (int first = 0; first < N; ++first) {
        ConvexPolygonPose pose_a = particle_pose_now(first);
        if (!isfinite(pose_a.position.x) || !isfinite(pose_a.position.y) ||
            !isfinite(pose_a.angle) || !isfinite(pose_a.scale)) {
            fprintf(stderr,
                "polygon stability nonfinite seed=%lu dissipative=%d event=%d type=%d t=%.17g particle=%d\n",
                seed, dissipative, event_number, (int)event_type, t, first);
            return 0;
        }
        for (int vertex = 0; vertex < shape.count; ++vertex) {
            const double gap = vertex_wall_gap(&pose_a,
                                               shape.vertices[vertex]);
            const double tolerance = 2.0e-8*
                (1.0 + pose_a.scale + Lx + Ly);
            const double diagnostic_tolerance = 1.0e-10*
                (1.0 + pose_a.scale + Lx + Ly);
            update_stability_minimum(&worst_wall, gap, seed, dissipative,
                                     event_number, first, vertex);
            if (gap < -diagnostic_tolerance) {
                fprintf(stderr,
                    "polygon stability wall overlap shape=%s seed=%lu dissipative=%d event=%d type=%d owner=%d,%d phase=%s t=%.17g particle=%d vertex=%d gap=%.17g diagnostic_tolerance=%.17g production_tolerance=%.17g production_failure=%d\n",
                    diagnostic_shape, seed, dissipative, event_number, (int)event_type,
                    diagnostic_event_i, diagnostic_event_j, diagnostic_phase, t,
                    first, vertex, gap, diagnostic_tolerance, tolerance,
                    gap < -tolerance);
                return 0;
            }
        }
        for (int second = first+1; second < N; ++second) {
            ConvexPolygonPose pose_b = particle_pose_now(second);
            double dx = pose_b.position.x-pose_a.position.x;
            double dy = pose_b.position.y-pose_a.position.y;
            PBC(&dx, &dy);
            pose_b.position.x = pose_a.position.x+dx;
            pose_b.position.y = pose_a.position.y+dy;
            ConvexPolygonContact contact;
            if (!convex_polygon_sat_contact(&shape, &pose_a, &shape, &pose_b,
                                            &contact)) {
                fprintf(stderr,
                    "polygon stability invalid SAT seed=%lu dissipative=%d event=%d type=%d t=%.17g pair=%d,%d\n",
                    seed, dissipative, event_number, (int)event_type, t,
                    first, second);
                return 0;
            }
            const double tolerance = 2.0e-8*
                (1.0 + pose_a.scale + pose_b.scale);
            const double diagnostic_tolerance = 1.0e-10*
                (1.0 + pose_a.scale + pose_b.scale);
            update_stability_minimum(&worst_pair, contact.separation, seed,
                                     dissipative, event_number, first, second);
            if (contact.separation < -diagnostic_tolerance) {
                fprintf(stderr,
                    "polygon stability pair overlap shape=%s seed=%lu dissipative=%d event=%d type=%d owner=%d,%d phase=%s t=%.17g pair=%d,%d separation=%.17g diagnostic_tolerance=%.17g production_tolerance=%.17g production_failure=%d\n",
                    diagnostic_shape, seed, dissipative, event_number, (int)event_type,
                    diagnostic_event_i, diagnostic_event_j, diagnostic_phase, t,
                    first, second, contact.separation, diagnostic_tolerance,
                    tolerance, contact.separation < -tolerance);
                return 0;
            }
        }
    }
    return 1;
}

static int process_stability_event(unsigned long seed, int dissipative,
                                   int event_number)
{
    nextEvent = findNextEvent();
    if (nextEvent == NULL) return 0;
    if (!isfinite(nextEvent->t) || nextEvent->t < t) {
        fprintf(stderr,
            "polygon stability invalid deadline seed=%lu dissipative=%d event=%d now=%.17g deadline=%.17g\n",
            seed, dissipative, event_number, t, nextEvent->t);
        return -1;
    }
    const enum event event_type = nextEvent->type;
    diagnostic_event_i = nextEvent->i;
    diagnostic_event_j = nextEvent->j;
    const double event_time = nextEvent->t;
    if ((unsigned)event_type < sizeof(stability_event_counts)/
                               sizeof(stability_event_counts[0]))
        stability_event_counts[event_type]++;
    if (event_type == POLYGON_RECHECK && diagnostic_event_i >= 0 &&
        diagnostic_event_i < 64)
        stability_recheck_owners[diagnostic_event_i]++;
    if (event_time == stability_last_event_time) {
        stability_same_time_events++;
    } else {
        stability_same_time_events = 1;
        stability_last_event_time = event_time;
    }
    if (stability_same_time_events > stability_max_same_time_events)
        stability_max_same_time_events = stability_same_time_events;
    if (stability_same_time_events > 4096) {
        fprintf(stderr,
            "polygon stability same-time livelock shape=%s seed=%lu dissipative=%d event=%d type=%d owner=%d,%d t=%.17g streak=%lu\n",
            diagnostic_shape, seed, dissipative, event_number,
            (int)event_type, diagnostic_event_i, diagnostic_event_j,
            event_time, stability_same_time_events);
        return -1;
    }
    const double interval_start = t;
    static const double fractions[] = {0.25, 0.5, 0.75};
    static const char *phases[] = {
        "interval-quarter", "interval-midpoint", "interval-three-quarter"
    };
    for (size_t sample = 0; sample < sizeof(fractions)/sizeof(fractions[0]);
         ++sample) {
        const double sample_time = interval_start +
            fractions[sample]*(event_time-interval_start);
        if (!(sample_time > interval_start) || !(sample_time < event_time))
            continue;
        t = sample_time;
        diagnostic_phase = phases[sample];
        if (!exact_invariants(seed, dissipative, event_number, event_type))
            return -1;
    }
    t = event_time;
    diagnostic_phase = "before";
    if (!exact_invariants(seed, dissipative, event_number, event_type) ||
        !check_energy_invariant(seed, dissipative, event_number, event_type,
                                diagnostic_phase))
        return -1;
    removeEventFromQueue(nextEvent);
    switch (event_type) {
        case CELLCROSS:
            ++ncross;
            doTheCrossing();
            break;
        case POLYGON_CELLCROSS:
            ++ncross;
            doTheCrossingPolygon();
            break;
        case COLLISION:
            polygon_particles_do_collision();
            break;
        case WALL:
            polygon_particles_do_wall();
            break;
        case POLYGON_RECHECK:
            polygon_particles_do_recheck();
            break;
        default:
            fprintf(stderr,
                "polygon stability unexpected event seed=%lu dissipative=%d event=%d type=%d t=%.17g\n",
                seed, dissipative, event_number, (int)event_type, t);
            return -1;
    }
    if (simulationHasFault()) {
        fprintf(stderr,
            "polygon stability runtime fault seed=%lu dissipative=%d event=%d type=%d t=%.17g message=%s\n",
            seed, dissipative, event_number, (int)event_type, t,
            simulationFaultMessage());
        return -1;
    }
    diagnostic_phase = "after";
    return exact_invariants(seed, dissipative, event_number, event_type) &&
           check_energy_invariant(seed, dissipative, event_number, event_type,
                                  diagnostic_phase) ? 1 : -1;
}

int main(int argc, char **argv)
{
    addWallx = addWally =
        EDMD_TEST_BOUNDARY_MODE == EDMD_TEST_BOUNDARY_PLANAR;
    addCircularWall =
        EDMD_TEST_BOUNDARY_MODE == EDMD_TEST_BOUNDARY_CIRCULAR;
    const int seed_count = argc > 1 ? atoi(argv[1]) : 12;
    const int event_limit = argc > 2 ? atoi(argv[2]) : 12000;
    const int requested_shape_count = argc > 3 ? atoi(argv[3]) : 4;
    const int shape_count = requested_shape_count > 0 &&
        requested_shape_count <= 4 ? requested_shape_count : 4;
    for (int shape_case = 0; shape_case < shape_count; ++shape_case) {
      for (int dissipative = 0; dissipative <= 1; ++dissipative) {
        for (int seed_index = 0; seed_index < seed_count; ++seed_index) {
            const unsigned long seed = 1009UL + 7919UL*(unsigned long)seed_index;
            if (!initialize_stability_case(seed, dissipative, shape_case)) {
                fprintf(stderr,
                    "polygon stability initialization failed shape=%s seed=%lu dissipative=%d message=%s\n",
                    diagnostic_shape, seed, dissipative,
                    simulationFaultMessage());
                return EXIT_FAILURE;
            }
            diagnostic_event_i = -1;
            diagnostic_event_j = -1;
            diagnostic_phase = "initial";
            reset_stability_minima();
            if (!exact_invariants(seed, dissipative, 0, COLLISION) ||
                !check_energy_invariant(seed, dissipative, 0, COLLISION,
                                        diagnostic_phase))
                return EXIT_FAILURE;
            int event_number;
            for (event_number = 1; event_number <= event_limit; ++event_number) {
                const int status = process_stability_event(
                    seed, dissipative, event_number);
                if (status < 0) return EXIT_FAILURE;
                if (status == 0) break;
            }
            fprintf(stdout,
                "polygon stability case ok shape=%s seed=%lu dissipative=%d events=%d t=%.9g physical_collisions=%lu collision_events=%lu crossing=%lu wall=%lu recheck=%lu max_same_time=%lu energy_initial=%.9g energy_min=%.9g energy_max=%.9g energy_max_rel_error=%.3g energy_max_step_increase=%.3g worst_pair=%.3g pair=%d,%d event=%d worst_wall=%.3g particle=%d vertex=%d event=%d\n",
                diagnostic_shape, seed, dissipative, event_number-1, t, ncol,
                stability_event_counts[COLLISION],
                stability_event_counts[CELLCROSS] +
                    stability_event_counts[POLYGON_CELLCROSS],
                stability_event_counts[WALL],
                stability_event_counts[POLYGON_RECHECK],
                stability_max_same_time_events,
                stability_initial_energy, stability_minimum_energy,
                stability_maximum_energy,
                stability_maximum_relative_energy_error,
                stability_maximum_step_energy_increase,
                worst_pair.value, worst_pair.first, worst_pair.second,
                worst_pair.event_number, worst_wall.value, worst_wall.first,
                worst_wall.second, worst_wall.event_number);
            if (stability_event_counts[POLYGON_RECHECK] > 0) {
                int maximum_owner = 0;
                for (int index = 1; index < N; ++index)
                    if (stability_recheck_owners[index] >
                        stability_recheck_owners[maximum_owner])
                        maximum_owner = index;
                fprintf(stdout,
                    "polygon stability recheck hotspot shape=%s seed=%lu dissipative=%d owner=%d count=%lu\n",
                    diagnostic_shape, seed, dissipative, maximum_owner,
                    stability_recheck_owners[maximum_owner]);
            }
            freeArrays();
        }
      }
    }
    puts("polygon stability matrix: ok");
    return EXIT_SUCCESS;
}
