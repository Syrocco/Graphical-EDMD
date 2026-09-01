#include "polygon_particles.h"

#include "convex_polygon.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef EDMD_POLYGON_SCHEDULER_STATS
#define EDMD_POLYGON_SCHEDULER_STATS 0
#endif

#if !POLYGON_PARTICLES || THREE_D
#error "polygon_particles.c is a two-dimensional optional source"
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

enum {
    POLYGON_WALL_LEFT = 0,
    POLYGON_WALL_RIGHT,
    POLYGON_WALL_BOTTOM,
    POLYGON_WALL_TOP,
    POLYGON_WALL_CIRCULAR,
    POLYGON_WALL_MID_LEFT,
    POLYGON_WALL_MID_RIGHT,
};

typedef struct PolygonParticleState {
    double angle;
    double angle_time;
    double omega;
    double inertia;
    double target_radius;
} PolygonParticleState;

extern int N;
extern int Nxcells, Nycells;
extern int polygonCollisionCellNeighborRange;
#if G || defined(EDMD_TEST_MUTABLE_CONFIG)
extern int addWallx, addWally, addMidWall, addCircularWall, addExpo;
extern bool addField, damping;
#else
extern const int addWallx, addWally, addMidWall, addCircularWall, addExpo;
extern const bool addField, damping;
#endif
extern double t, Lx, Ly, halfLx, halfLy, vr, Einit;
extern double res, resW, resT;
extern double collTermX, collTermY, collTermXY, collTermYX;
extern double dp, dpLeft, dpRight, dpMid, cumulativeWallImpulse;
extern unsigned long int ncol;
extern particle *particles;
extern particle **cellList;
extern node **eventList;
extern node *nextEvent;
extern void (*collisionEvent)(int);
extern void (*crossingEvent)(int);

static SimParticleShapeConfig configured_shape;
static bool configured_shape_initialized = false;
static ConvexPolygon canonical_shape;
static PolygonParticleState *polygon_state = NULL;
static int polygon_state_count = 0;
static bool polygon_active = false;
static bool polygon_growing = false;
static double recovery_precision_until = -DBL_MAX;
static unsigned long growth_progress_collision_origin = 0;
static double growth_progress_time = 0.0;
static double growth_stop_time = never;
static double growth_kinetic_energy = 0.0;
static bool growth_kinetic_energy_valid = false;

#if EDMD_POLYGON_SCHEDULER_STATS
static PolygonParticleSchedulerStats scheduler_stats;
#define SCHEDULER_INCREMENT(field) (++scheduler_stats.field)
#define SCHEDULER_ADD(field, value) \
    (scheduler_stats.field += (unsigned long long)(value))
#else
#define SCHEDULER_INCREMENT(field) ((void)0)
#define SCHEDULER_ADD(field, value) ((void)0)
#endif

static double maximum(double a, double b) { return a > b ? a : b; }
static double minimum(double a, double b) { return a < b ? a : b; }

static ConvexPolygonVec2 vec(double x, double y) {
    const ConvexPolygonVec2 result = {x, y};
    return result;
}

static ConvexPolygonVec2 add(ConvexPolygonVec2 a, ConvexPolygonVec2 b) {
    return vec(a.x + b.x, a.y + b.y);
}

static ConvexPolygonVec2 subtract(ConvexPolygonVec2 a, ConvexPolygonVec2 b) {
    return vec(a.x - b.x, a.y - b.y);
}

static ConvexPolygonVec2 scale(ConvexPolygonVec2 a, double factor) {
    return vec(a.x*factor, a.y*factor);
}

static double dot(ConvexPolygonVec2 a, ConvexPolygonVec2 b) {
    return a.x*b.x + a.y*b.y;
}

static double cross(ConvexPolygonVec2 a, ConvexPolygonVec2 b) {
    return a.x*b.y - a.y*b.x;
}

static double length(ConvexPolygonVec2 a) {
    return sqrt(dot(a, a));
}

void polygon_particles_reset_scheduler_stats(void) {
#if EDMD_POLYGON_SCHEDULER_STATS
    memset(&scheduler_stats, 0, sizeof(scheduler_stats));
#endif
}

void polygon_particles_get_scheduler_stats(
    PolygonParticleSchedulerStats *stats) {
    if (stats == NULL) return;
#if EDMD_POLYGON_SCHEDULER_STATS
    *stats = scheduler_stats;
#else
    memset(stats, 0, sizeof(*stats));
#endif
}

/* Exact active-set enumeration for the one- or two-contact frictionless LCP
   K*j >= b, j >= 0, j_i*(K*j-b)_i = 0.  Contact matrices here are symmetric
   positive semidefinite rigid-body effective-mass matrices, so at most four
   active sets are required and no iterative solver/tuning is needed. */
static bool solve_contact_lcp(const double matrix[2][2],
                              const double bound[2], int count,
                              double velocity_tolerance,
                              double impulse[2]) {
    static const int masks_two[4] = {3, 1, 2, 0};
    static const int masks_one[2] = {1, 0};
    const int *masks = count == 2 ? masks_two : masks_one;
    const int mask_count = count == 2 ? 4 : 2;
    if (count < 1 || count > 2 || impulse == NULL) return false;
    for (int row = 0; row < count; ++row) {
        if (!(matrix[row][row] > DBL_MIN) ||
            !isfinite(matrix[row][row]) || !isfinite(bound[row])) return false;
        for (int column = 0; column < count; ++column)
            if (!isfinite(matrix[row][column])) return false;
    }

    const double tolerance = maximum(velocity_tolerance,
        2048.0*DBL_EPSILON*(1.0 + fabs(bound[0]) +
                            (count == 2 ? fabs(bound[1]) : 0.0)));
    for (int candidate_index = 0; candidate_index < mask_count;
         ++candidate_index) {
        const int mask = masks[candidate_index];
        double candidate[2] = {0, 0};
        bool solvable = true;
        if (mask == 1) {
            candidate[0] = bound[0]/matrix[0][0];
        } else if (mask == 2) {
            candidate[1] = bound[1]/matrix[1][1];
        } else if (mask == 3) {
            const double determinant = matrix[0][0]*matrix[1][1] -
                matrix[0][1]*matrix[1][0];
            const double determinant_scale =
                fabs(matrix[0][0]*matrix[1][1]) +
                fabs(matrix[0][1]*matrix[1][0]);
            solvable = determinant > 1024.0*DBL_EPSILON*
                maximum(DBL_MIN, determinant_scale) && isfinite(determinant);
            if (solvable) {
                candidate[0] = (bound[0]*matrix[1][1] -
                    bound[1]*matrix[0][1])/determinant;
                candidate[1] = (bound[1]*matrix[0][0] -
                    bound[0]*matrix[1][0])/determinant;
            }
        }
        if (!solvable) continue;
        for (int contact = 0; contact < count; ++contact) {
            const double impulse_tolerance = tolerance/matrix[contact][contact];
            if (!isfinite(candidate[contact]) ||
                candidate[contact] < -impulse_tolerance)
                solvable = false;
            if (candidate[contact] < 0) candidate[contact] = 0;
        }
        if (!solvable) continue;

        for (int row = 0; row < count; ++row) {
            double residual = -bound[row];
            for (int column = 0; column < count; ++column)
                residual += matrix[row][column]*candidate[column];
            if (!isfinite(residual) || residual < -tolerance ||
                ((mask & (1 << row)) && fabs(residual) > 8.0*tolerance)) {
                solvable = false;
                break;
            }
        }
        if (solvable) {
            impulse[0] = candidate[0];
            impulse[1] = candidate[1];
            return true;
        }
    }
    return false;
}

static bool evaluate_contact_impulse(const double matrix[2][2],
                                     const double velocity[2], int count,
                                     const double impulse[2],
                                     double response[2],
                                     double *linear_work,
                                     double *quadratic_work) {
    if (count < 1 || count > 2 || matrix == NULL || velocity == NULL ||
        impulse == NULL || response == NULL || linear_work == NULL ||
        quadratic_work == NULL) return false;
    response[0] = response[1] = 0;
    *linear_work = *quadratic_work = 0;
    for (int row = 0; row < count; ++row) {
        if (!isfinite(velocity[row]) || !isfinite(impulse[row]) ||
            impulse[row] < 0) return false;
        for (int column = 0; column < count; ++column) {
            if (!isfinite(matrix[row][column])) return false;
            response[row] += matrix[row][column]*impulse[column];
        }
        if (!isfinite(response[row])) return false;
        *linear_work += impulse[row]*velocity[row];
        *quadratic_work += impulse[row]*response[row];
    }
    return isfinite(*linear_work) && isfinite(*quadratic_work);
}

static double contact_work_tolerance(double linear_work,
                                     double quadratic_work) {
    return 4096.0*DBL_EPSILON*
        (1.0 + fabs(linear_work) + fabs(quadratic_work));
}

/* A simultaneous frictionless LCP enforces nonpenetration, but a mixed
   two-point manifold (one closing endpoint and one separating endpoint) does
   not, by itself, enforce the kinetic-energy equality for an elastic impact.
   Preserve the LCP impulse direction and project its magnitude onto the
   nonzero elastic root of

       dK(lambda) = lambda*j.u + lambda^2*j.K.j/2 = 0.

   Positive scaling preserves repulsiveness.  The caller additionally gets a
   guarantee that every represented contact is nonclosing after projection;
   if roundoff or a degenerate effective-mass matrix prevents both guarantees,
   the caller can try a propagative resolution before faulting. */
static bool project_elastic_contact_energy(const double matrix[2][2],
                                           const double velocity[2],
                                           int count,
                                           double velocity_tolerance,
                                           double impulse[2]) {
    if (count < 1 || count > 2 || matrix == NULL || velocity == NULL ||
        impulse == NULL || !(velocity_tolerance >= 0) ||
        !isfinite(velocity_tolerance)) return false;

    double response[2] = {0, 0};
    double linear_work = 0;
    double quadratic_work = 0;
    if (!evaluate_contact_impulse(matrix, velocity, count, impulse,
                                  response, &linear_work, &quadratic_work))
        return false;

    const double work_tolerance =
        contact_work_tolerance(linear_work, quadratic_work);
    const double energy_change = linear_work + 0.5*quadratic_work;
    if (!isfinite(work_tolerance) || !isfinite(energy_change)) return false;
    double multiplier = 1.0;
    if (fabs(energy_change) > work_tolerance) {
        if (!(linear_work < 0) || !(quadratic_work > DBL_MIN)) return false;
        multiplier = -2.0*linear_work/quadratic_work;
        if (!(multiplier > 0) || !isfinite(multiplier)) return false;
    }

    for (int row = 0; row < count; ++row) {
        const double outgoing_velocity =
            velocity[row] + multiplier*response[row];
        if (!isfinite(outgoing_velocity) ||
            outgoing_velocity < -velocity_tolerance) return false;
    }
    for (int contact = 0; contact < count; ++contact) {
        impulse[contact] *= multiplier;
        if (!isfinite(impulse[contact])) return false;
    }
    return true;
}

/* If a simultaneous impulse direction cannot be energy-projected without
   making another support point close, resolve a bounded sequence of exact
   one-point elastic reflections.  Each step selects the most-closing point
   and is individually energy-neutral; the final projection removes only
   accumulated floating-point work. */
static bool solve_propagative_elastic_contact(const double matrix[2][2],
                                              const double velocity[2],
                                              int count,
                                              double velocity_tolerance,
                                              double impulse[2]) {
    enum { MAX_PROPAGATIVE_IMPULSES = 512 };
    if (count < 1 || count > 2 || matrix == NULL || velocity == NULL ||
        impulse == NULL || !(velocity_tolerance >= 0) ||
        !isfinite(velocity_tolerance)) return false;
    double current_velocity[2] = {0, 0};
    double candidate[2] = {0, 0};
    for (int contact = 0; contact < count; ++contact) {
        if (!(matrix[contact][contact] > DBL_MIN) ||
            !isfinite(matrix[contact][contact]) ||
            !isfinite(velocity[contact])) return false;
        current_velocity[contact] = velocity[contact];
        for (int column = 0; column < count; ++column)
            if (!isfinite(matrix[contact][column])) return false;
    }

    for (int step = 0; step <= MAX_PROPAGATIVE_IMPULSES; ++step) {
        int closing_contact = -1;
        double closing_velocity = -velocity_tolerance;
        for (int contact = 0; contact < count; ++contact) {
            if (current_velocity[contact] < closing_velocity) {
                closing_velocity = current_velocity[contact];
                closing_contact = contact;
            }
        }
        if (closing_contact < 0) {
            if (!project_elastic_contact_energy(
                    matrix, velocity, count, velocity_tolerance, candidate))
                return false;
            impulse[0] = candidate[0];
            impulse[1] = candidate[1];
            return true;
        }
        if (step == MAX_PROPAGATIVE_IMPULSES) return false;

        const double increment = -2.0*closing_velocity/
            matrix[closing_contact][closing_contact];
        if (!(increment > 0) || !isfinite(increment) ||
            !isfinite(candidate[closing_contact] + increment)) return false;
        candidate[closing_contact] += increment;
        for (int row = 0; row < count; ++row) {
            current_velocity[row] += matrix[row][closing_contact]*increment;
            if (!isfinite(current_velocity[row])) return false;
        }
    }
    return false;
}

/* Production elastic contacts must conserve energy.  Dissipative contacts
   may lose energy, but a mixed manifold is never allowed to add it.  Growth
   callers deliberately bypass this policy because prescribed scale motion
   performs physical work that is absent from the translational/rotational
   effective-mass energy represented here. */
static bool enforce_contact_energy_policy(const double matrix[2][2],
                                          const double velocity[2], int count,
                                          double velocity_tolerance,
                                          bool elastic,
                                          double impulse[2]) {
    double response[2] = {0, 0};
    double linear_work = 0;
    double quadratic_work = 0;
    if (!evaluate_contact_impulse(matrix, velocity, count, impulse,
                                  response, &linear_work, &quadratic_work))
        return false;
    const double energy_change = linear_work + 0.5*quadratic_work;
    const double work_tolerance =
        contact_work_tolerance(linear_work, quadratic_work);
    if (!isfinite(work_tolerance) || !isfinite(energy_change)) return false;
    if (!elastic && energy_change <= work_tolerance) return true;
    if (project_elastic_contact_energy(
            matrix, velocity, count, velocity_tolerance, impulse)) return true;
    return solve_propagative_elastic_contact(
        matrix, velocity, count, velocity_tolerance, impulse);
}

#if EDMD_TEST_PREPARED_STATE
bool polygon_particles_test_solve_contact_lcp(
        const double matrix[2][2], const double bound[2], int count,
        double velocity_tolerance, double impulse[2]) {
    return solve_contact_lcp(matrix, bound, count, velocity_tolerance, impulse);
}

bool polygon_particles_test_project_elastic_contact_energy(
        const double matrix[2][2], const double velocity[2], int count,
        double velocity_tolerance, double impulse[2]) {
    return project_elastic_contact_energy(matrix, velocity, count,
                                          velocity_tolerance, impulse);
}

bool polygon_particles_test_enforce_contact_energy_policy(
        const double matrix[2][2], const double velocity[2], int count,
        double velocity_tolerance, bool elastic, double impulse[2]) {
    return enforce_contact_energy_policy(matrix, velocity, count,
                                         velocity_tolerance, elastic, impulse);
}
#endif

static ConvexPolygonVec2 rotate(ConvexPolygonVec2 value, double angle) {
    const double cosine = cos(angle);
    const double sine = sin(angle);
    return vec(cosine*value.x - sine*value.y,
               sine*value.x + cosine*value.y);
}

static void ensure_default_config(void) {
    if (configured_shape_initialized) return;
    memset(&configured_shape, 0, sizeof(configured_shape));
    configured_shape.mode = SIM_PARTICLE_CIRCLES;
    configured_shape.vertex_count = 6;
    for (int vertex = 0; vertex < configured_shape.vertex_count; ++vertex) {
        const double angle = 2.0*M_PI*vertex/configured_shape.vertex_count;
        configured_shape.vertices[vertex][0] = cos(angle);
        configured_shape.vertices[vertex][1] = sin(angle);
    }
    configured_shape.random_orientation = true;
    configured_shape.angular_temperature = 1.0;
    configured_shape_initialized = true;
}

static bool canonicalize_config(SimParticleShapeConfig *destination,
                                const SimParticleShapeConfig *source) {
    ConvexPolygon candidate;
    if (destination == NULL || source == NULL || source->vertex_count < 3 ||
        source->vertex_count > SIM_POLYGON_MAX_VERTICES) return false;
    memset(&candidate, 0, sizeof(candidate));
    candidate.count = source->vertex_count;
    for (int index = 0; index < candidate.count; ++index) {
        candidate.vertices[index].x = source->vertices[index][0];
        candidate.vertices[index].y = source->vertices[index][1];
    }
    if (!convex_polygon_recenter_and_normalize(&candidate)) return false;
    memset(destination, 0, sizeof(*destination));
    destination->mode = source->mode;
    destination->vertex_count = candidate.count;
    for (int index = 0; index < candidate.count; ++index) {
        destination->vertices[index][0] = candidate.vertices[index].x;
        destination->vertices[index][1] = candidate.vertices[index].y;
    }
    destination->random_orientation = source->random_orientation;
    destination->angular_temperature = source->angular_temperature;
    if (!isfinite(destination->angular_temperature) ||
        destination->angular_temperature < 0) {
        destination->angular_temperature = 0;
    }
    return true;
}

void polygon_particles_get_config(SimParticleShapeConfig *config) {
    ensure_default_config();
    if (config != NULL) *config = configured_shape;
}

void polygon_particles_set_config(const SimParticleShapeConfig *config) {
    SimParticleShapeConfig normalized;
    ensure_default_config();
    if (config == NULL) return;
    if (config->mode == SIM_PARTICLE_CIRCLES) {
        /* Preserve the editable template while the disk model is selected. */
        if (canonicalize_config(&normalized, config)) configured_shape = normalized;
        configured_shape.mode = SIM_PARTICLE_CIRCLES;
        return;
    }
    if (config->mode == SIM_PARTICLE_CONVEX_POLYGONS &&
        canonicalize_config(&normalized, config)) {
        configured_shape = normalized;
    }
}

bool polygon_particles_requested(void) {
    ensure_default_config();
    return configured_shape.mode == SIM_PARTICLE_CONVEX_POLYGONS;
}

bool polygon_particles_active(void) {
    return polygon_active && polygon_state != NULL &&
           polygon_state_count == N;
}

bool polygon_particles_growing(void) {
    return polygon_particles_active() && polygon_growing;
}

static bool rebuild_canonical_shape(void) {
    ensure_default_config();
    memset(&canonical_shape, 0, sizeof(canonical_shape));
    canonical_shape.count = configured_shape.vertex_count;
    for (int index = 0; index < canonical_shape.count; ++index) {
        canonical_shape.vertices[index].x =
            configured_shape.vertices[index][0];
        canonical_shape.vertices[index].y =
            configured_shape.vertices[index][1];
    }
    return convex_polygon_recenter_and_normalize(&canonical_shape) != 0;
}

bool polygon_particles_prepare(void) {
    polygon_particles_release();
    if (!polygon_particles_requested()) return true;
    if (N <= 0 || particles == NULL || !rebuild_canonical_shape()) return false;
    polygon_state = calloc((size_t)N, sizeof(*polygon_state));
    if (polygon_state == NULL) return false;
    polygon_state_count = N;
    for (int index = 0; index < N; ++index) {
        PolygonParticleState *state = polygon_state + index;
        state->angle = configured_shape.random_orientation ?
            drand(0.0, 2.0*M_PI) : 0.0;
        state->angle_time = t;
    }
    return true;
}

void polygon_particles_release(void) {
    free(polygon_state);
    polygon_state = NULL;
    polygon_state_count = 0;
    polygon_active = false;
    polygon_growing = false;
    growth_stop_time = never;
    growth_kinetic_energy = 0.0;
    growth_kinetic_energy_valid = false;
    recovery_precision_until = -DBL_MAX;
}

size_t polygon_particles_checkpoint_bytes(void) {
    return polygon_particles_active()
        ? (size_t)polygon_state_count*sizeof(*polygon_state) : 0;
}

bool polygon_particles_capture_checkpoint(void *destination, size_t bytes) {
    const size_t required = polygon_particles_checkpoint_bytes();
    if (required == 0 || destination == NULL || bytes != required)
        return false;
    memcpy(destination, polygon_state, required);
    return true;
}

bool polygon_particles_restore_checkpoint(const void *source, size_t bytes) {
    const size_t required = polygon_particles_checkpoint_bytes();
    if (required == 0 || source == NULL || bytes != required)
        return false;
    const PolygonParticleState *saved = source;
    for (int index = 0; index < polygon_state_count; ++index) {
        if (!isfinite(saved[index].angle) ||
            !isfinite(saved[index].angle_time) ||
            !isfinite(saved[index].omega) ||
            !isfinite(saved[index].inertia) ||
            !(saved[index].inertia > DBL_MIN) ||
            !isfinite(saved[index].target_radius) ||
            !(saved[index].target_radius > DBL_MIN))
            return false;
    }
    memcpy(polygon_state, source, required);
    return true;
}

void polygon_particles_enable_recovery_precision(double duration) {
    if (!polygon_particles_active() || !isfinite(duration) || duration <= 0.0)
        return;
    const double requested_until = t + duration;
    if (isfinite(requested_until) && requested_until > recovery_precision_until)
        recovery_precision_until = requested_until;
}

static bool recovery_precision_active(void) {
    return polygon_particles_active() && t <= recovery_precision_until;
}

#if EDMD_TEST_PREPARED_STATE
bool polygon_particles_test_recovery_precision_active(void) {
    return recovery_precision_active();
}
#endif

static double gaussian_unit(void) {
    const double first = maximum(drand(0.0, 1.0), DBL_MIN);
    const double second = drand(0.0, 1.0);
    return sqrt(-2.0*log(first))*cos(2.0*M_PI*second);
}

void polygon_particles_activate(void) {
    if (!polygon_particles_requested() || polygon_state == NULL ||
        polygon_state_count != N) return;
    if (!polygon_active) {
        for (int index = 0; index < N; ++index) {
            PolygonParticleState *state = polygon_state + index;
            const particle *value = particles + index;
            state->angle = polygon_particles_angle_at(index, t);
            state->angle_time = t;
            state->inertia = value->m*value->rad*value->rad*
                canonical_shape.inertia_per_mass;
            state->target_radius = value->rad;
            if (state->inertia > DBL_MIN &&
                configured_shape.angular_temperature > 0) {
                state->omega = gaussian_unit()*sqrt(
                    configured_shape.angular_temperature/state->inertia);
            } else {
                state->omega = 0;
            }
        }
    }
    polygon_active = true;
}

void polygon_particles_begin_growth(void) {
    if (!polygon_particles_requested() || polygon_state == NULL ||
        polygon_state_count != N) return;
    polygon_active = true;
    polygon_growing = true;
    growth_stop_time = eventList != NULL && eventList[2*N + 5] != NULL &&
        isfinite(eventList[2*N + 5]->t) ? eventList[2*N + 5]->t :
        (vr > DBL_MIN ? t + 1.0/vr : never);
    growth_progress_collision_origin = ncol;
    growth_progress_time = t;
    for (int index = 0; index < N; ++index) {
        PolygonParticleState *state = polygon_state + index;
        /* Use the final-radius inertia throughout preparation. This stays
           nonsingular near zero scale while allowing collisions to rotate and
           unjam anisotropic bodies. Preparation spin is discarded at the
           production handoff. */
        state->angle_time = t;
        state->omega = 0.0;
        const double target_radius = particles[index].rad +
            (vr > DBL_MIN ? particles[index].vr/vr : 0.0);
        state->target_radius = target_radius;
        state->inertia = particles[index].m*target_radius*target_radius*
            canonical_shape.inertia_per_mass;
    }
    polygon_particles_refresh_growth_energy();
}

static double body_kinetic_energy(const particle *value,
                                  const PolygonParticleState *state) {
    if (value == NULL || state == NULL || !(value->m > 0.0) ||
        !(state->inertia > 0.0) || !isfinite(value->m) ||
        !isfinite(state->inertia)) return NAN;
    return 0.5*value->m*
        (value->vx*value->vx + value->vy*value->vy) +
        0.5*state->inertia*state->omega*state->omega;
}

static double growth_energy_ceiling_per_particle(void) {
    return 128.0*maximum(1.0,
        Einit + configured_shape.angular_temperature);
}

static bool validate_growth_kinetic_energy(void) {
    if (!polygon_growing) return true;
    const double ceiling = (double)N*growth_energy_ceiling_per_particle();
    if (!growth_kinetic_energy_valid ||
        !isfinite(growth_kinetic_energy) || growth_kinetic_energy < 0.0) {
        simulationReportFault(
            "Polygon preparation produced non-finite kinetic energy before reaching the target packing fraction. Lower the target or growth rate and restart.");
        return false;
    }
    if (!isfinite(ceiling) || !(ceiling > 0.0) ||
        growth_kinetic_energy > ceiling) {
        simulationReportFault(
            "Polygon preparation reached its aggregate kinetic-energy safety ceiling before the target packing fraction. Lower the target or growth rate and restart.");
        return false;
    }
    return true;
}

bool polygon_particles_refresh_growth_energy(void) {
    if (!polygon_growing) return true;
    double sum = 0.0;
    double correction = 0.0;
    for (int index = 0; index < N; ++index) {
        const double kinetic = body_kinetic_energy(
            particles + index, polygon_state + index);
        if (!isfinite(kinetic) || kinetic < 0.0) {
            growth_kinetic_energy = NAN;
            growth_kinetic_energy_valid = false;
            return validate_growth_kinetic_energy();
        }
        /* Kahan summation keeps the cache reproducible when a hot outlier is
           present among many ordinary particles. */
        const double adjusted = kinetic - correction;
        const double updated = sum + adjusted;
        correction = (updated - sum) - adjusted;
        sum = updated;
    }
    growth_kinetic_energy = sum;
    growth_kinetic_energy_valid = true;
    return validate_growth_kinetic_energy();
}

bool polygon_particles_set_mass(int index, double mass) {
    if (!polygon_particles_active() || index < 0 ||
        index >= polygon_state_count || !isfinite(mass) || !(mass > 0.0))
        return false;

    particle *value = particles + index;
    PolygonParticleState *state = polygon_state + index;
    /* Preparation intentionally uses final-size inertia even when the current
       homothetic scale is near zero.  Production uses the physical radius
       stored by the ordinary particle engine. */
    const double inertia_radius = polygon_growing
        ? state->target_radius : value->rad;
    const double inertia = mass*inertia_radius*inertia_radius*
        canonical_shape.inertia_per_mass;
    if (!isfinite(inertia_radius) || !(inertia_radius > DBL_MIN) ||
        !isfinite(inertia) || !(inertia > DBL_MIN))
        return false;

    value->m = mass;
    state->inertia = inertia;
    /* Do not rescale omega: a paused live edit must not introduce an angular
       velocity discontinuity.  A full refresh also accounts for velocity
       fields that the same GUI transaction may have changed immediately
       before the mass. */
    return polygon_particles_refresh_growth_energy();
}

static bool update_growth_kinetic_energy(double before, double after) {
    if (!polygon_growing) return true;
    if (!growth_kinetic_energy_valid || !isfinite(before) || before < 0.0 ||
        !isfinite(after) || after < 0.0)
        return polygon_particles_refresh_growth_energy();
    const double updated = growth_kinetic_energy + (after - before);
    if (!isfinite(updated) || updated < 0.0)
        return polygon_particles_refresh_growth_energy();
    growth_kinetic_energy = updated;
    return validate_growth_kinetic_energy();
}

#if EDMD_TEST_PREPARED_STATE
double polygon_particles_test_growth_energy_ceiling(void) {
    return growth_energy_ceiling_per_particle();
}
#endif

static void check_growth_progress(void) {
    if (!polygon_growing || simulationHasFault()) return;
    const unsigned long particle_count = N > 0 ? (unsigned long)N : 1UL;
    const unsigned long progress_interval = 256UL*particle_count;
    if (ncol - growth_progress_collision_origin >= progress_interval) {
        if (!(t > growth_progress_time)) {
            simulationReportFault(
                "Polygon preparation entered a same-time collision jam before reaching the target packing fraction. Lower the target or growth rate and restart.");
            return;
        }
        growth_progress_collision_origin = ncol;
        growth_progress_time = t;
    }
}

void polygon_particles_finish_growth(void) {
    if (!polygon_particles_growing()) return;
    polygon_growing = false;
    growth_stop_time = never;
    growth_kinetic_energy = 0.0;
    growth_kinetic_energy_valid = false;
    for (int index = 0; index < N; ++index) {
        PolygonParticleState *state = polygon_state + index;
        particle *value = particles + index;
        state->angle = polygon_particles_angle_at(index, t);
        state->angle_time = t;
        value->rad = state->target_radius;
        state->inertia = value->m*value->rad*value->rad*
            canonical_shape.inertia_per_mass;
        if (state->inertia > DBL_MIN &&
            configured_shape.angular_temperature > 0.0) {
            state->omega = gaussian_unit()*sqrt(
                configured_shape.angular_temperature/state->inertia);
        } else {
            state->omega = 0.0;
        }
        value->vr = 0.0;
    }
}

double polygon_particles_angle_at(int index, double time) {
    if (polygon_state == NULL || index < 0 || index >= polygon_state_count)
        return 0;
    const PolygonParticleState *state = polygon_state + index;
    double result = state->angle + state->omega*(time - state->angle_time);
    if (isfinite(result) && fabs(result) > 8.0*M_PI) result = remainder(result, 2.0*M_PI);
    return result;
}

double polygon_particles_omega(int index) {
    return polygon_state != NULL && index >= 0 && index < polygon_state_count
        ? polygon_state[index].omega : 0;
}

double polygon_particles_inertia(int index) {
    return polygon_state != NULL && index >= 0 && index < polygon_state_count
        ? polygon_state[index].inertia : 0;
}

static void synchronize_angle(int index) {
    if (polygon_state == NULL || index < 0 || index >= polygon_state_count) return;
    PolygonParticleState *state = polygon_state + index;
    state->angle = polygon_particles_angle_at(index, t);
    state->angle_time = t;
}

void polygon_particles_free_fly(particle *value) {
    if (value == NULL) return;
    synchronize_angle(value->num);
    if (polygon_growing)
        freeFlyGrow(value);
    else
        freeFlyNormal(value);
}

int polygon_particles_vertex_count(void) {
    return polygon_particles_requested() ? canonical_shape.count : 0;
}

bool polygon_particles_vertex(int vertex_index, double *x, double *y) {
    if (!polygon_particles_requested() || vertex_index < 0 ||
        vertex_index >= canonical_shape.count || x == NULL || y == NULL) return false;
    *x = canonical_shape.vertices[vertex_index].x;
    *y = canonical_shape.vertices[vertex_index].y;
    return true;
}

double polygon_particles_unit_area(void) {
    ensure_default_config();
    return simulation_particle_shape_unit_area(&configured_shape);
}

double polygon_particles_rotational_energy(void) {
    if (!polygon_particles_active()) return 0;
    double energy = 0;
    for (int index = 0; index < polygon_state_count; ++index) {
        const PolygonParticleState *state = polygon_state + index;
        energy += 0.5*state->inertia*state->omega*state->omega;
    }
    return energy;
}

bool polygon_particles_apply_angular_langevin(int index,
                                              double temperature,
                                              double gamma,
                                              double interval,
                                              bool euler) {
    if (!polygon_particles_active() || index < 0 || index >= polygon_state_count ||
        !isfinite(temperature) || temperature < 0.0 ||
        !isfinite(gamma) || !(gamma > 0.0) ||
        !isfinite(interval) || !(interval > 0.0) || euler) {
        return false;
    }
    PolygonParticleState *state = polygon_state + index;
    if (!isfinite(state->omega) || !isfinite(state->inertia) ||
        !(state->inertia > DBL_MIN)) {
        return false;
    }

    /* Match the analytical Ornstein-Uhlenbeck update used for translation.
       Gamma is the engine's inverse-time relaxation rate.  The angular
       variance T/I is the equipartition value for one rotational degree of
       freedom, so this remains exact for any kick interval. */
    synchronize_angle(index);
    const double attenuation = exp(-gamma*interval);
    const double variance = temperature*(1.0-attenuation*attenuation)/
        state->inertia;
    if (!isfinite(attenuation) || attenuation < 0.0 ||
        attenuation > 1.0 || !isfinite(variance) || variance < 0.0) {
        return false;
    }
    state->omega = attenuation*state->omega;
    if (variance > 0.0)
        state->omega += sqrt(variance)*gaussian_unit();
    return isfinite(state->omega);
}

static particlePrediction prediction_at(const particle *value, double time) {
    return polygon_growing
        ? predictFreeFlyGrow(value, time)
        : predictFreeFlyNormal(value, time);
}

static ConvexPolygonPose pose_at(int index, double time) {
    const particlePrediction prediction = prediction_at(particles + index, time);
    ConvexPolygonPose pose;
    pose.position = vec(prediction.x, prediction.y);
    pose.angle = polygon_particles_angle_at(index, time);
    pose.scale = prediction.rad;
    return pose;
}

static ConvexPolygonMotion motion_at(int index, double time) {
    const particlePrediction prediction = prediction_at(particles + index, time);
    ConvexPolygonMotion motion;
    motion.pose.position = vec(prediction.x, prediction.y);
    motion.pose.angle = polygon_particles_angle_at(index, time);
    motion.pose.scale = prediction.rad;
    motion.velocity = vec(prediction.vx, prediction.vy);
    motion.angular_velocity = polygon_particles_omega(index);
    motion.scale_velocity = polygon_growing ? particles[index].vr : 0.0;
    return motion;
}

static void minimum_image_shift(ConvexPolygonVec2 first,
                                ConvexPolygonVec2 *second) {
    double dx = second->x - first.x;
    double dy = second->y - first.y;
    PBC(&dx, &dy);
    second->x = first.x + dx;
    second->y = first.y + dy;
}

bool polygon_particles_contains_point(int index, double x, double y,
                                      double time) {
    if (!polygon_particles_requested() || index < 0 || index >= N) return false;
    const ConvexPolygonPose pose = pose_at(index, time);
    double dx = x - pose.position.x;
    double dy = y - pose.position.y;
    PBC(&dx, &dy);
    if (!(pose.scale > 0)) return false;
    const ConvexPolygonVec2 local = scale(rotate(vec(dx, dy), -pose.angle),
                                          1.0/pose.scale);
    const double tolerance = 64.0*DBL_EPSILON;
    for (int edge = 0; edge < canonical_shape.count; ++edge) {
        const ConvexPolygonVec2 a = canonical_shape.vertices[edge];
        const ConvexPolygonVec2 b = canonical_shape.vertices[
            (edge + 1) % canonical_shape.count];
        if (cross(subtract(b, a), subtract(local, a)) < -tolerance) return false;
    }
    return true;
}

static int neighboring_coordinates(
    int center, int count,
    int output[POLYGON_PARTICLES_MAX_CELL_COORDINATES]) {
    if (output == NULL || count <= 0 || center < 0 || center >= count)
        return 0;
    const int range = polygonCollisionCellNeighborRange < 0 ? 0 :
        (polygonCollisionCellNeighborRange >
             POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE
         ? POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE
         : polygonCollisionCellNeighborRange);
    const int stencil_width = 2*range + 1;
    if (count <= stencil_width) {
        /* Every periodic coordinate is a neighbour. Listing it once avoids
           wrap/dedup work in the very small boxes where it matters most. */
        for (int coordinate = 0; coordinate < count; ++coordinate)
            output[coordinate] = coordinate;
        return count;
    }
    for (int offset = -range; offset <= range; ++offset) {
        int coordinate = center + offset;
        if (coordinate < 0) coordinate += count;
        else if (coordinate >= count) coordinate -= count;
        output[offset + range] = coordinate;
    }
    return stencil_width;
}

#if EDMD_TEST_PREPARED_STATE
int polygon_particles_test_neighboring_coordinates(
    int center, int count,
    int output[POLYGON_PARTICLES_MAX_CELL_COORDINATES]) {
    return neighboring_coordinates(center, count, output);
}
#endif

static double scheduling_horizon(int index,
                                 const ConvexPolygonMotion *owner_motion) {
    double horizon = never;
    const double crossing_time = eventList[index]->t;
    if (isfinite(crossing_time) && crossing_time < never)
        horizon = maximum(0.0, crossing_time - t);
    const double angular_speed = owner_motion != NULL
        ? fabs(owner_motion->angular_velocity) :
          fabs(polygon_particles_omega(index));
    if (horizon >= never) {
        if (angular_speed > DBL_EPSILON) {
            /* Rotation does not invalidate linked-cell neighbors because the
               stored radius is the circumscribed radius.  Predict one complete
               orientation period instead of manufacturing a recheck every few
               degrees; the numerical CCD already bounds angular motion. */
            horizon = 2.0*M_PI/angular_speed;
        } else {
            /* A truly static owner needs only the t=now overlap check. Any
               moving partner owns its future prediction, so periodic refresh
               events would prevent a quiescent simulation from pausing. */
            const double linear_speed = owner_motion != NULL
                ? hypot(owner_motion->velocity.x, owner_motion->velocity.y) :
                  hypot(particles[index].vx, particles[index].vy);
            horizon = linear_speed > DBL_EPSILON
                ? hypot(Lx, Ly)/linear_speed : 0.0;
            if (!isfinite(horizon) || horizon <= 0.0) {
                /* Scale-aware CCD can search the whole remaining growth
                   interval for a stationary center. Avoid one-unit rechecks
                   and their repeated full neighbor scans. */
                horizon = polygon_growing &&
                    particles[index].vr > DBL_MIN ? never : 0.0;
            }
        }
    }
    if (polygon_growing && particles[index].vr > DBL_MIN) {
        /* Use the authoritative fixed GROWSTOP deadline, not a radius-derived
           reconstruction. Near the target, independent floating operations
           can otherwise schedule a contact a few ulps after GROWSTOP and let
           the bucketed calendar expose that later event first. */
        const double remaining_growth = isfinite(growth_stop_time) ?
            maximum(0.0, growth_stop_time - t) :
            maximum(0.0,
                (polygon_state[index].target_radius -
                 (owner_motion != NULL ? owner_motion->pose.scale :
                  particles[index].rad))/
                    particles[index].vr);
        if (isfinite(remaining_growth))
            horizon = minimum(horizon, remaining_growth);
    }
    return horizon;
}

static ConvexPolygonCcdStatus pair_time_of_impact(int first, int second,
                                                  double horizon,
                                                  const ConvexPolygonMotion *owner_motion,
                                                  ConvexPolygonCcdResult *result) {
    SCHEDULER_INCREMENT(pair_candidates);
    ConvexPolygonMotion a = owner_motion != NULL
        ? *owner_motion : motion_at(first, t);
    ConvexPolygonMotion b = motion_at(second, t);
    minimum_image_shift(a.pose.position, &b.pose.position);
    ConvexPolygonCcdSettings settings;
    convex_polygon_ccd_default_settings(&settings);
    if (recovery_precision_active()) {
        settings.max_iterations = 4096U;
    } else if (polygon_growing) {
        /* Preparation velocities are discarded at GROWSTOP, so resolving a
           contact to machine epsilon only burns CPU without improving the
           generated configuration.  Resolving contact a few parts per billion
           before exact touching is conservative: it cannot create an overlap,
           and it removes roughly a dozen SAT bisection rounds from the dense
           part of growth.  Production and recovery retain full precision. */
        settings.distance_tolerance = 1.0e-9;
        settings.time_tolerance = 1.0e-9*fmax(1.0, horizon);
        settings.max_iterations = 256U;
    } else {
        settings.max_iterations = 512U;
    }
    ConvexPolygonCcdResult local_result;
    ConvexPolygonCcdResult *stored_result = result != NULL
        ? result : &local_result;
    ConvexPolygonCcdStatus status = convex_polygon_time_of_impact(
        &canonical_shape, &a, &canonical_shape, &b,
        horizon, &settings, stored_result);
    if (status == CONVEX_POLYGON_CCD_NO_HIT &&
        stored_result->iterations == 0)
        SCHEDULER_INCREMENT(pair_broadphase_rejections);
    else
        SCHEDULER_INCREMENT(pair_ccd_calls);
    SCHEDULER_ADD(pair_ccd_iterations, stored_result->iterations);
    return status;
}

static ConvexPolygonVec2 world_vertex(const ConvexPolygonPose *pose,
                                      int vertex_index) {
    return add(pose->position, scale(rotate(
        canonical_shape.vertices[vertex_index], pose->angle), pose->scale));
}

typedef struct PolygonWallManifold {
    double gap;
    int count;
    int tied_count;
    ConvexPolygonVec2 deepest_normal;
    ConvexPolygonVec2 point[2];
    ConvexPolygonVec2 normal[2];
    double outward_speed[2];
} PolygonWallManifold;

static bool wall_manifold_at(int index, double time, int wall,
                             PolygonWallManifold *manifold) {
    if (manifold == NULL) return false;
    memset(manifold, 0, sizeof(*manifold));
    const ConvexPolygonPose pose = pose_at(index, time);
    const particlePrediction prediction = prediction_at(particles + index, time);
    const double omega = polygon_particles_omega(index);
    ConvexPolygonVec2 points[CONVEX_POLYGON_MAX_VERTICES] = {{0}};
    ConvexPolygonVec2 normals[CONVEX_POLYGON_MAX_VERTICES] = {{0}};
    double gaps[CONVEX_POLYGON_MAX_VERTICES] = {0};
    double speeds[CONVEX_POLYGON_MAX_VERTICES] = {0};
    double best = DBL_MAX;
    int best_index = -1;
    for (int vertex_index = 0; vertex_index < canonical_shape.count;
         ++vertex_index) {
        const ConvexPolygonVec2 point = world_vertex(&pose, vertex_index);
        double gap = DBL_MAX;
        ConvexPolygonVec2 normal = vec(0, 0);
        switch (wall) {
            case POLYGON_WALL_LEFT:
                gap = point.x; normal = vec(-1, 0); break;
            case POLYGON_WALL_RIGHT:
                gap = Lx - point.x; normal = vec(1, 0); break;
            case POLYGON_WALL_BOTTOM:
                gap = point.y; normal = vec(0, -1); break;
            case POLYGON_WALL_TOP:
                gap = Ly - point.y; normal = vec(0, 1); break;
            case POLYGON_WALL_MID_LEFT:
                gap = halfLx - point.x; normal = vec(1, 0); break;
            case POLYGON_WALL_MID_RIGHT:
                gap = point.x - halfLx; normal = vec(-1, 0); break;
            case POLYGON_WALL_CIRCULAR: {
                const ConvexPolygonVec2 radial = subtract(point,
                    vec(halfLx, halfLy));
                const double radius = length(radial);
                gap = halfLx - radius;
                normal = radius > DBL_MIN ? scale(radial, 1.0/radius) : vec(1, 0);
                break;
            }
            default: return false;
        }
        points[vertex_index] = point;
        normals[vertex_index] = normal;
        gaps[vertex_index] = gap;
        const ConvexPolygonVec2 arm = subtract(point, pose.position);
        speeds[vertex_index] = prediction.vx*normal.x +
            prediction.vy*normal.y + omega*cross(arm, normal);
        if (polygon_growing && pose.scale > DBL_MIN)
            speeds[vertex_index] += prediction.rad > 0.0
                ? particles[index].vr*dot(arm, normal)/pose.scale : 0.0;
        if (!isfinite(gap) || !isfinite(speeds[vertex_index]) ||
            !isfinite(normal.x) || !isfinite(normal.y)) return false;
        if (gap < best) {
            best = gap;
            best_index = vertex_index;
        }
    }

    /* Strict convexity gives at most two tied support vertices against a
       plane. Circular confinement can tie more only in a zero-clearance
       degenerate state, so count every tie while retaining the two fastest
       ordinary response contacts. Keep the actual points: averaging them
       hides V +/- omega*a and can treat a spinning, approaching endpoint as
       a stationary midpoint. */
    const double support_tolerance = 4096.0*DBL_EPSILON*
        (1.0 + pose.scale + fabs(pose.position.x) +
         fabs(pose.position.y) + fabs(Lx) + fabs(Ly));
    int first = -1;
    int second = -1;
    for (int vertex_index = 0; vertex_index < canonical_shape.count;
         ++vertex_index) {
        if (gaps[vertex_index] <= best + support_tolerance) {
            manifold->tied_count++;
            if (first < 0 || speeds[vertex_index] > speeds[first]) {
                second = first;
                first = vertex_index;
            } else if (second < 0 || speeds[vertex_index] > speeds[second]) {
                second = vertex_index;
            }
        }
    }
    if (first < 0 || best_index < 0 || !isfinite(best)) return false;
    manifold->gap = best;
    manifold->deepest_normal = normals[best_index];
    manifold->count = second >= 0 ? 2 : 1;
    manifold->point[0] = points[first];
    manifold->normal[0] = normals[first];
    manifold->outward_speed[0] = speeds[first];
    if (second >= 0) {
        manifold->point[1] = points[second];
        manifold->normal[1] = normals[second];
        manifold->outward_speed[1] = speeds[second];
    }
    return true;
}

static void overlap_verification_message(char *message, size_t capacity,
                                         const char *format, int first,
                                         int second, double time,
                                         double separation,
                                         double tolerance) {
    if (message == NULL || capacity == 0) return;
    (void)snprintf(message, capacity, format, first, second, time,
                   separation, tolerance);
}

bool polygon_particles_verify_no_overlaps(double time,
                                          double penetration_tolerance,
                                          char *message,
                                          size_t message_capacity) {
    if (message != NULL && message_capacity > 0) message[0] = '\0';
    if (!polygon_particles_active() || !isfinite(time) ||
        !isfinite(penetration_tolerance) || penetration_tolerance < 0.0) {
        if (message != NULL && message_capacity > 0)
            (void)snprintf(message, message_capacity,
                "Convex-particle overlap verification received invalid state or tolerance.");
        return false;
    }

    for (int first = 0; first < N; ++first) {
        int walls[7];
        int wall_count = 0;
        if (addWallx) {
            walls[wall_count++] = POLYGON_WALL_LEFT;
            walls[wall_count++] = POLYGON_WALL_RIGHT;
        }
        if (addWally) {
            walls[wall_count++] = POLYGON_WALL_BOTTOM;
            walls[wall_count++] = POLYGON_WALL_TOP;
        }
        if (addCircularWall) walls[wall_count++] = POLYGON_WALL_CIRCULAR;
        if (addMidWall && particles[first].type == 1) {
            const ConvexPolygonPose pose = pose_at(first, time);
            walls[wall_count++] = pose.position.x <= halfLx ?
                POLYGON_WALL_MID_LEFT : POLYGON_WALL_MID_RIGHT;
        }
        for (int wall_index = 0; wall_index < wall_count; ++wall_index) {
            PolygonWallManifold manifold;
            if (!wall_manifold_at(first, time, walls[wall_index], &manifold)) {
                if (message != NULL && message_capacity > 0)
                    (void)snprintf(message, message_capacity,
                        "Convex-particle overlap verification could not evaluate particle %d against wall %d at t=%.9g.",
                        first, walls[wall_index], time);
                return false;
            }
            if (manifold.gap < -penetration_tolerance) {
                overlap_verification_message(message, message_capacity,
                    "Convex particle %d intersects wall %d at t=%.9g (gap %.17g, tolerance %.17g).",
                    first, walls[wall_index], time, manifold.gap,
                    penetration_tolerance);
                return false;
            }
        }

        const ConvexPolygonPose pose_a = pose_at(first, time);
        for (int second = first + 1; second < N; ++second) {
            ConvexPolygonPose pose_b = pose_at(second, time);
            minimum_image_shift(pose_a.position, &pose_b.position);
            ConvexPolygonContact contact;
            if (!convex_polygon_sat_contact(&canonical_shape, &pose_a,
                    &canonical_shape, &pose_b, &contact)) {
                if (message != NULL && message_capacity > 0)
                    (void)snprintf(message, message_capacity,
                        "Convex-particle overlap verification could not evaluate particles %d and %d at t=%.9g.",
                        first, second, time);
                return false;
            }
            if (contact.separation < -penetration_tolerance) {
                overlap_verification_message(message, message_capacity,
                    "Convex particles %d and %d overlap at t=%.9g (SAT separation %.17g, tolerance %.17g).",
                    first, second, time, contact.separation,
                    penetration_tolerance);
                return false;
            }
        }
    }
    return true;
}

bool polygon_particles_verify_local_no_overlaps(
    double time, double penetration_tolerance,
    char *message, size_t message_capacity) {
    if (message != NULL && message_capacity > 0) message[0] = '\0';
    if (!polygon_particles_active() || !isfinite(time) || time != t ||
        !isfinite(penetration_tolerance) || penetration_tolerance < 0.0 ||
        cellList == NULL || Nxcells <= 0 || Nycells <= 0) {
        if (message != NULL && message_capacity > 0)
            (void)snprintf(message, message_capacity,
                "Local convex-particle overlap verification requires the active linked-cell index at the current simulation time.");
        return false;
    }

    for (int first = 0; first < N; ++first) {
        int walls[7];
        int wall_count = 0;
        if (addWallx) {
            walls[wall_count++] = POLYGON_WALL_LEFT;
            walls[wall_count++] = POLYGON_WALL_RIGHT;
        }
        if (addWally) {
            walls[wall_count++] = POLYGON_WALL_BOTTOM;
            walls[wall_count++] = POLYGON_WALL_TOP;
        }
        if (addCircularWall) walls[wall_count++] = POLYGON_WALL_CIRCULAR;
        if (addMidWall && particles[first].type == 1) {
            const ConvexPolygonPose pose = pose_at(first, time);
            walls[wall_count++] = pose.position.x <= halfLx ?
                POLYGON_WALL_MID_LEFT : POLYGON_WALL_MID_RIGHT;
        }
        for (int wall_index = 0; wall_index < wall_count; ++wall_index) {
            PolygonWallManifold manifold;
            if (!wall_manifold_at(first, time, walls[wall_index], &manifold)) {
                if (message != NULL && message_capacity > 0)
                    (void)snprintf(message, message_capacity,
                        "Local convex-particle overlap verification could not evaluate particle %d against wall %d at t=%.9g.",
                        first, walls[wall_index], time);
                return false;
            }
            if (manifold.gap < -penetration_tolerance) {
                overlap_verification_message(message, message_capacity,
                    "Convex particle %d intersects wall %d at t=%.9g (gap %.17g, tolerance %.17g).",
                    first, walls[wall_index], time, manifold.gap,
                    penetration_tolerance);
                return false;
            }
        }

        const int cell_x = particles[first].cell[0];
        const int cell_y = particles[first].cell[1];
        int neighboring_x[POLYGON_PARTICLES_MAX_CELL_COORDINATES];
        int neighboring_y[POLYGON_PARTICLES_MAX_CELL_COORDINATES];
        const int count_x = neighboring_coordinates(
            cell_x, Nxcells, neighboring_x);
        const int count_y = neighboring_coordinates(
            cell_y, Nycells, neighboring_y);
        if (count_x <= 0 || count_y <= 0) {
            if (message != NULL && message_capacity > 0)
                (void)snprintf(message, message_capacity,
                    "Local convex-particle overlap verification found an invalid cell owner for particle %d.",
                    first);
            return false;
        }

        const ConvexPolygonPose pose_a = pose_at(first, time);
        for (int yi = 0; yi < count_y; ++yi) {
            for (int xi = 0; xi < count_x; ++xi) {
                particle *other = cellList[
                    neighboring_y[yi]*Nxcells + neighboring_x[xi]];
                while (other != NULL) {
                    const int second = other->num;
                    if (second > first) {
                        ConvexPolygonPose pose_b = pose_at(second, time);
                        minimum_image_shift(pose_a.position,
                                            &pose_b.position);
                        const ConvexPolygonVec2 delta = subtract(
                            pose_b.position, pose_a.position);
                        const double envelope = canonical_shape.radius*
                            (pose_a.scale + pose_b.scale);
                        const double filter_tolerance =
                            penetration_tolerance + 2048.0*DBL_EPSILON*
                            (1.0 + fabs(delta.x) + fabs(delta.y) +
                             fabs(envelope));
                        if (!isfinite(envelope) || envelope < 0.0 ||
                            !isfinite(filter_tolerance)) {
                            if (message != NULL && message_capacity > 0)
                                (void)snprintf(message, message_capacity,
                                    "Local convex-particle overlap verification encountered invalid geometry for particles %d and %d.",
                                    first, second);
                            return false;
                        }
                        if (hypot(delta.x, delta.y) <=
                            envelope + filter_tolerance) {
                            ConvexPolygonContact contact;
                            if (!convex_polygon_sat_contact(
                                    &canonical_shape, &pose_a,
                                    &canonical_shape, &pose_b, &contact)) {
                                if (message != NULL && message_capacity > 0)
                                    (void)snprintf(message, message_capacity,
                                        "Local convex-particle overlap verification could not evaluate particles %d and %d at t=%.9g.",
                                        first, second, time);
                                return false;
                            }
                            if (contact.separation <
                                -penetration_tolerance) {
                                overlap_verification_message(message,
                                    message_capacity,
                                    "Convex particles %d and %d overlap at t=%.9g (SAT separation %.17g, tolerance %.17g).",
                                    first, second, time,
                                    contact.separation,
                                    penetration_tolerance);
                                return false;
                            }
                        }
                    }
                    other = other->nxt;
                }
            }
        }
    }
    return true;
}

static bool wall_conservative_advance(double current, double horizon,
                                      double certified_step, double *next) {
    if (next == NULL || !(horizon > current) ||
        !(certified_step > 0) || !isfinite(certified_step)) return false;
    const double candidate = current +
        minimum(certified_step, horizon - current);
    if (!(candidate > current) ||
        !(t + candidate > t + current)) {
        simulationReportFault(polygon_growing
            ? "Polygon preparation jammed against a wall before reaching the target packing fraction. Lower the target or growth rate and restart."
            : "Convex-wall collision prediction reached the resolution limit of the simulation clock; restart from a nearer time origin.");
        return false;
    }
    *next = candidate;
    return true;
}

/* The polygon is contained in its circumscribed disk at every orientation.
   With ballistic center motion and linear scale, each planar envelope gap is
   linear.  For the circular container, |center(t)-boxCenter|+radius(t) is
   convex, so its maximum on an interval is attained at an endpoint. */
static bool wall_envelope_may_contact(int index, int wall, double horizon) {
    if (index < 0 || index >= N || !isfinite(horizon) || horizon < 0.0 ||
        addField || damping)
        return true;
    const particlePrediction start = prediction_at(particles + index, t);
    const particlePrediction finish = prediction_at(
        particles + index, t + horizon);
    const double radius_start = canonical_shape.radius*start.rad;
    const double radius_finish = canonical_shape.radius*finish.rad;
    if (!isfinite(start.x) || !isfinite(start.y) ||
        !isfinite(finish.x) || !isfinite(finish.y) ||
        !isfinite(radius_start) || !isfinite(radius_finish) ||
        radius_start < 0.0 || radius_finish < 0.0)
        return true;
    const double tolerance = 2048.0*DBL_EPSILON*
        (1.0 + fabs(Lx) + fabs(Ly) + radius_start + radius_finish);
    double gap_start;
    double gap_finish;
    switch (wall) {
        case POLYGON_WALL_LEFT:
            gap_start = start.x - radius_start;
            gap_finish = finish.x - radius_finish;
            break;
        case POLYGON_WALL_RIGHT:
            gap_start = Lx - start.x - radius_start;
            gap_finish = Lx - finish.x - radius_finish;
            break;
        case POLYGON_WALL_BOTTOM:
            gap_start = start.y - radius_start;
            gap_finish = finish.y - radius_finish;
            break;
        case POLYGON_WALL_TOP:
            gap_start = Ly - start.y - radius_start;
            gap_finish = Ly - finish.y - radius_finish;
            break;
        case POLYGON_WALL_MID_LEFT:
            gap_start = halfLx - start.x - radius_start;
            gap_finish = halfLx - finish.x - radius_finish;
            break;
        case POLYGON_WALL_MID_RIGHT:
            gap_start = start.x - halfLx - radius_start;
            gap_finish = finish.x - halfLx - radius_finish;
            break;
        case POLYGON_WALL_CIRCULAR:
            gap_start = halfLx - hypot(start.x-halfLx, start.y-halfLy) -
                radius_start;
            gap_finish = halfLx - hypot(finish.x-halfLx,
                                        finish.y-halfLy) - radius_finish;
            break;
        default:
            return true;
    }
    if (!isfinite(gap_start) || !isfinite(gap_finish)) return true;
    return minimum(gap_start, gap_finish) <= tolerance;
}

static bool wall_time_of_impact(int index, int wall, double horizon,
                                double *hit_time, bool *overlap,
                                bool *quiescent) {
    const particlePrediction prediction = prediction_at(particles + index, t);
    const double radius = prediction.rad;
    const double final_radius = polygon_growing
        ? radius + particles[index].vr*horizon : radius;
    const double maximum_radius = maximum(radius, final_radius);
    const double speed_bound = hypot(prediction.vx, prediction.vy) +
        fabs(polygon_particles_omega(index))*maximum_radius +
        (polygon_growing ? fabs(particles[index].vr) : 0.0);
    const double distance_tolerance = 256.0*DBL_EPSILON*
        (1.0 + fabs(Lx) + fabs(Ly) + maximum_radius);
    double local_time = 0;
    if (overlap != NULL) *overlap = false;
    if (quiescent != NULL) *quiescent = false;
    const int maximum_iterations = recovery_precision_active() ? 4096 : 512;
    for (int iteration = 0; iteration < maximum_iterations; ++iteration) {
        PolygonWallManifold manifold;
        SCHEDULER_INCREMENT(wall_manifold_evaluations);
        if (!wall_manifold_at(index, t + local_time, wall, &manifold))
            return false;
        const double gap = manifold.gap;
        if (local_time == 0 && gap < -8.0*distance_tolerance) {
            if (overlap != NULL) *overlap = true;
            return false;
        }
        if (gap <= distance_tolerance) {
            const double velocity_tolerance =
                convex_polygon_ccd_velocity_tolerance(speed_bound);
            if (manifold.tied_count > 2 &&
                manifold.outward_speed[0] > velocity_tolerance) {
                simulationReportFault(
                    "A convex particle touches the circular wall at more than two approaching vertices and has no uniquely resolvable free motion. Reduce the packing fraction and restart.");
                return false;
            }
            if (manifold.tied_count > 2) {
                /* Three or more circular-wall normals span the translational
                   plane.  With no appreciable outward speed this is a
                   quiescent constrained state (not a reason to manufacture
                   epsilon-sized probes forever). */
                if (quiescent != NULL) *quiescent = true;
                return false;
            }
            for (int contact = 0; contact < manifold.count; ++contact)
                if (manifold.outward_speed[contact] > velocity_tolerance) {
                    *hit_time = local_time;
                    return true;
                }
            if (!(speed_bound > 0)) return false;
            const double safe_probe = 0.5*distance_tolerance/speed_bound;
            double next_time;
            if (local_time >= horizon) return false;
            if (!wall_conservative_advance(local_time, horizon, safe_probe,
                                           &next_time)) return false;
            local_time = next_time;
            continue;
        }
        if (!(speed_bound > 0)) return false;
        const double step = 0.8*gap/speed_bound;
        if (!(step > 0) || !isfinite(step)) return false;
        if (step > horizon - local_time) return false;
        double next_time;
        if (!wall_conservative_advance(local_time, horizon, step, &next_time))
            return false;
        local_time = next_time;
    }
    /* Iteration exhaustion is not proof of safety. Wake the owner at the last
       certified time and let the scheduler continue from there. */
    *hit_time = local_time;
    return false;
}

static bool event_deadline_after(double delay, double *deadline) {
    if (deadline == NULL || !isfinite(delay) || delay < 0) {
        simulationReportFault(
            "Convex-particle collision prediction produced an invalid event delay.");
        return false;
    }
    const double candidate = t + delay;
    if (!isfinite(candidate) || (delay > 0 && !(candidate > t))) {
        simulationReportFault(
            "Convex-particle collision prediction reached the resolution limit of the simulation clock; restart from a nearer time origin.");
        return false;
    }
    *deadline = candidate;
    return true;
}

#if EDMD_TEST_PREPARED_STATE
bool polygon_particles_test_event_deadline_after(double delay,
                                                  double *deadline) {
    return event_deadline_after(delay, deadline);
}
#endif

static void schedule_recheck(int index, double delay) {
    double at_time;
    if (!event_deadline_after(delay, &at_time)) return;
    if (polygon_growing && at_time >= growth_stop_time) return;
    if (!(at_time > t) || at_time >= never || at_time > eventList[index]->t)
        return;
    node *event = eventList[N + index];
    event->i = index;
    event->j = -1;
    event->type = POLYGON_RECHECK;
    event->t = at_time;
    event->collActual = particles[index].coll;
    addEventToQueue(event);
    SCHEDULER_INCREMENT(rechecks_scheduled);
}

static void schedule_wall(int index, int wall, double at_time) {
    if (polygon_growing && at_time >= growth_stop_time) return;
    if (at_time > eventList[index]->t) return;
    removeEventFromQueue(eventList[index]);
    node *event = eventList[N + index];
    event->i = index;
    event->j = wall;
    event->type = WALL;
    event->t = at_time;
    event->collActual = particles[index].coll;
    addEventToQueue(event);
}

void polygon_particles_crossing_event(int index) {
    if (polygon_growing) {
        crossingEventGrow(index);
        /* GROWSTOP is the sole owner of the preparation boundary. Suppress a
           cell crossing at or after it so Paul-bucket rounding cannot expose
           a later particle event first. stopGrow() synchronizes every body
           exactly and rebuilds ordinary polygon predictions. */
        if (index >= 0 && index < N && eventList != NULL &&
            eventList[index] != NULL && eventList[index]->q >= 0 &&
            eventList[index]->t >= growth_stop_time)
            removeEventFromQueue(eventList[index]);
    } else {
        crossingEventNormal(index);
    }
    if (index >= 0 && index < N && eventList != NULL &&
        eventList[index] != NULL) {
        eventList[index]->type = POLYGON_CELLCROSS;
    }
}

static void report_overlap(int first, int second, const char *surface) {
    char message[256];
    if (second >= 0) {
        (void)snprintf(message, sizeof(message),
            "Convex particles %d and %d overlap at t=%.9g. The simulation is paused; change the shape or density and restart.",
            first, second, t);
    } else {
        (void)snprintf(message, sizeof(message),
            "Convex particle %d intersects the %s at t=%.9g. The simulation is paused; change the boundary or shape and restart.",
            first, surface != NULL ? surface : "wall", t);
    }
    simulationReportRecoverableOverlapFault(message);
}

void polygon_particles_collision_event(int index) {
    if (!polygon_particles_active() || index < 0 || index >= N) return;
    SCHEDULER_INCREMENT(collision_event_calls);
    const ConvexPolygonMotion owner_motion = motion_at(index, t);
    const double horizon = scheduling_horizon(index, &owner_motion);
    double best_time = horizon;
    double refresh_time = horizon;
    int best_partner = -1;
    int best_wall = -1;

    const particle *value = particles + index;
    int walls[7];
    int wall_count = 0;
    if (addWallx) {
        walls[wall_count++] = POLYGON_WALL_LEFT;
        walls[wall_count++] = POLYGON_WALL_RIGHT;
    }
    if (addWally) {
        walls[wall_count++] = POLYGON_WALL_BOTTOM;
        walls[wall_count++] = POLYGON_WALL_TOP;
    }
    if (addCircularWall) walls[wall_count++] = POLYGON_WALL_CIRCULAR;
    if (addMidWall && value->type == 1) {
        walls[wall_count++] = owner_motion.pose.position.x <= halfLx ?
            POLYGON_WALL_MID_LEFT : POLYGON_WALL_MID_RIGHT;
    }
    for (int wall_index = 0; wall_index < wall_count; ++wall_index) {
        SCHEDULER_INCREMENT(wall_candidates);
        if (!wall_envelope_may_contact(index, walls[wall_index], best_time)) {
            SCHEDULER_INCREMENT(wall_broadphase_rejections);
            continue;
        }
        double wall_time = 0;
        bool overlap = false;
        bool quiescent = false;
        if (wall_time_of_impact(index, walls[wall_index], best_time,
                                &wall_time, &overlap, &quiescent)) {
            if (wall_time <= best_time && best_partner < 0) {
                best_time = wall_time;
                best_wall = walls[wall_index];
                best_partner = -1;
            }
        } else if (overlap) {
            report_overlap(index, -1,
                walls[wall_index] == POLYGON_WALL_CIRCULAR ?
                    "circular wall" : "planar wall");
            return;
        } else if (quiescent) {
            /* The complete multi-support circular constraint has already
               ruled out future wall motion for this stationary center. Pair
               CCD may still install an earlier explicit refresh below. */
            refresh_time = never;
        } else if (wall_time > 0 && wall_time < refresh_time) {
            refresh_time = wall_time;
        }
    }

    int neighboring_x[POLYGON_PARTICLES_MAX_CELL_COORDINATES];
    int neighboring_y[POLYGON_PARTICLES_MAX_CELL_COORDINATES];
    const int count_x = neighboring_coordinates(value->cell[0], Nxcells,
                                                  neighboring_x);
    const int count_y = neighboring_coordinates(value->cell[1], Nycells,
                                                  neighboring_y);
    for (int yi = 0; yi < count_y; ++yi) {
        for (int xi = 0; xi < count_x; ++xi) {
            particle *other = cellList[neighboring_y[yi]*Nxcells +
                                       neighboring_x[xi]];
            while (other != NULL) {
                if (other->num != index) {
                    ConvexPolygonCcdResult result;
                    const ConvexPolygonCcdStatus status = pair_time_of_impact(
                        index, other->num, best_time, &owner_motion, &result);
                    if (status == CONVEX_POLYGON_CCD_INITIAL_OVERLAP) {
                        report_overlap(index, other->num, NULL);
                        return;
                    }
                    if (status == CONVEX_POLYGON_CCD_INVALID) {
                        simulationReportFault(
                            "Convex-particle collision prediction received invalid geometry or motion; restart with a valid shape and finite state.");
                        return;
                    }
                    if (status == CONVEX_POLYGON_CCD_ITERATION_LIMIT) {
                        if (!(result.time > 0)) {
                            simulationReportFault(polygon_growing
                                ? "Polygon preparation jammed before reaching the target packing fraction: collision search cannot make representable forward progress. Lower the target or growth rate and restart."
                                : "Convex-particle collision prediction could not make representable forward progress; restart from a nearer time origin.");
                            return;
                        }
                        if (result.time < refresh_time)
                            refresh_time = result.time;
                    }
                    if (status == CONVEX_POLYGON_CCD_HIT &&
                        result.time <= best_time) {
                        best_time = result.time;
                        best_partner = other->num;
                        best_wall = -1;
                    }
                }
                other = other->nxt;
            }
        }
    }

    if (simulationHasFault()) return;
    /* An iteration-limit result certifies safety only up to refresh_time.
       Never discard that earlier deadline merely because another pair or
       wall has a finite hit later in the current horizon: doing so can leave
       a rotating near-contact asleep until after it overlaps. */
    if (refresh_time > 0 && refresh_time < best_time) {
        schedule_recheck(index, refresh_time);
    } else if (best_partner >= 0) {
        double deadline;
        if (!event_deadline_after(best_time, &deadline)) return;
        if (polygon_growing && deadline >= growth_stop_time) return;
        addCollisionEvent(index, best_partner, deadline);
    } else if (best_wall >= 0) {
        double deadline;
        if (!event_deadline_after(best_time, &deadline)) return;
        schedule_wall(index, best_wall, deadline);
    } else if (refresh_time > 0 && refresh_time < never) {
        schedule_recheck(index, refresh_time);
    }
}

static ConvexPolygonContact contact_now(int first, int second) {
    ConvexPolygonPose a = pose_at(first, t);
    ConvexPolygonPose b = pose_at(second, t);
    minimum_image_shift(a.position, &b.position);
    ConvexPolygonContact contact;
    memset(&contact, 0, sizeof(contact));
    (void)convex_polygon_sat_contact(&canonical_shape, &a,
        &canonical_shape, &b, &contact);
    return contact;
}

static int contact_manifold_now(int first, int second,
                                ConvexPolygonVec2 normal,
                                ConvexPolygonVec2 points[2]) {
    ConvexPolygonPose pose_a = pose_at(first, t);
    ConvexPolygonPose pose_b = pose_at(second, t);
    minimum_image_shift(pose_a.position, &pose_b.position);
    const ConvexPolygonVec2 tangent = vec(-normal.y, normal.x);
    double maximum_a = -DBL_MAX;
    double minimum_b = DBL_MAX;
    ConvexPolygonVec2 vertices_a[CONVEX_POLYGON_MAX_VERTICES];
    ConvexPolygonVec2 vertices_b[CONVEX_POLYGON_MAX_VERTICES];
    for (int vertex_index = 0; vertex_index < canonical_shape.count;
         ++vertex_index) {
        vertices_a[vertex_index] = world_vertex(&pose_a, vertex_index);
        vertices_b[vertex_index] = world_vertex(&pose_b, vertex_index);
        maximum_a = maximum(maximum_a,
                            dot(vertices_a[vertex_index], normal));
        minimum_b = minimum(minimum_b,
                            dot(vertices_b[vertex_index], normal));
    }
    const double feature_tolerance = 2.0e-7*
        (1.0 + pose_a.scale + pose_b.scale) + 512.0*DBL_EPSILON*
        (1.0 + fabs(pose_a.position.x) + fabs(pose_a.position.y) +
         fabs(pose_b.position.x) + fabs(pose_b.position.y));
    double minimum_tangent_a = DBL_MAX;
    double maximum_tangent_a = -DBL_MAX;
    double minimum_tangent_b = DBL_MAX;
    double maximum_tangent_b = -DBL_MAX;
    for (int vertex_index = 0; vertex_index < canonical_shape.count;
         ++vertex_index) {
        if (maximum_a - dot(vertices_a[vertex_index], normal) <=
            feature_tolerance) {
            const double coordinate = dot(vertices_a[vertex_index], tangent);
            minimum_tangent_a = minimum(minimum_tangent_a, coordinate);
            maximum_tangent_a = maximum(maximum_tangent_a, coordinate);
        }
        if (dot(vertices_b[vertex_index], normal) - minimum_b <=
            feature_tolerance) {
            const double coordinate = dot(vertices_b[vertex_index], tangent);
            minimum_tangent_b = minimum(minimum_tangent_b, coordinate);
            maximum_tangent_b = maximum(maximum_tangent_b, coordinate);
        }
    }
    if (minimum_tangent_a == DBL_MAX || minimum_tangent_b == DBL_MAX)
        return 0;
    double lower = maximum(minimum_tangent_a, minimum_tangent_b);
    double upper = minimum(maximum_tangent_a, maximum_tangent_b);
    if (lower > upper) {
        /* A tolerance-level vertex/vertex contact can leave the projected
           support intervals microscopically disjoint. Their midpoint remains
           the unique contact point. */
        const double middle = 0.5*(lower + upper);
        lower = middle;
        upper = middle;
    }
    const double normal_coordinate = 0.5*(maximum_a + minimum_b);
    points[0] = add(scale(normal, normal_coordinate),
                    scale(tangent, lower));
    if (upper - lower > feature_tolerance) {
        points[1] = add(scale(normal, normal_coordinate),
                        scale(tangent, upper));
        return 2;
    }
    return 1;
}

static ConvexPolygonVec2 velocity_at_point(const particle *value,
                                           const PolygonParticleState *state,
                                           ConvexPolygonVec2 arm) {
    ConvexPolygonVec2 velocity = vec(
        value->vx - state->omega*arm.y,
        value->vy + state->omega*arm.x);
    if (polygon_growing && value->rad > DBL_MIN)
        velocity = add(velocity, scale(arm, value->vr/value->rad));
    return velocity;
}

void polygon_particles_do_collision(void) {
    const int first = nextEvent->i;
    const int second = nextEvent->j;
    if (!polygon_particles_active() || first < 0 || second < 0 ||
        first >= N || second >= N) return;
    particle *a = particles + first;
    particle *b = particles + second;
    /* Most discarded pair events are invalidated by a newer collision of the
       partner.  Their owner remains a valid lazy trajectory, so check the
       generation before paying to synchronize its position and angle. */
    if (nextEvent->collActual != b->coll) {
        polygon_particles_collision_event(first);
        return;
    }
    polygon_particles_free_fly(a);
    polygon_particles_free_fly(b);
    ConvexPolygonContact contact = contact_now(first, second);
    const double contact_tolerance = 2.0e-8*
        (1.0 + a->rad + b->rad);
    if (contact.separation < -contact_tolerance) {
        report_overlap(first, second, NULL);
        return;
    }

    PolygonParticleState *state_a = polygon_state + first;
    PolygonParticleState *state_b = polygon_state + second;
    ConvexPolygonVec2 center_delta = vec(b->x - a->x, b->y - a->y);
    PBC(&center_delta.x, &center_delta.y);
    const ConvexPolygonVec2 center_a = vec(a->x, a->y);
    const ConvexPolygonVec2 center_b = add(center_a, center_delta);

    ConvexPolygonVec2 contact_points[2];
    int contact_count = contact_manifold_now(first, second, contact.normal,
                                              contact_points);
    if (contact_count < 1) {
        contact_points[0] = contact.point;
        contact_count = 1;
    }
    ConvexPolygonVec2 arms_a[2], arms_b[2];
    double normal_speed[2] = {0, 0};
    double maximum_closing_speed = 0;
    for (int contact_index = 0; contact_index < contact_count;
         ++contact_index) {
        arms_a[contact_index] = subtract(contact_points[contact_index], center_a);
        arms_b[contact_index] = subtract(contact_points[contact_index], center_b);
        const ConvexPolygonVec2 velocity_a = velocity_at_point(
            a, state_a, arms_a[contact_index]);
        const ConvexPolygonVec2 velocity_b = velocity_at_point(
            b, state_b, arms_b[contact_index]);
        normal_speed[contact_index] = dot(subtract(velocity_b, velocity_a),
                                          contact.normal);
        maximum_closing_speed = maximum(maximum_closing_speed,
                                        -normal_speed[contact_index]);
    }
    const double pair_speed_bound = hypot(b->vx-a->vx, b->vy-a->vy) +
        fabs(state_a->omega)*a->rad + fabs(state_b->omega)*b->rad +
        (polygon_growing ? fabs(a->vr) + fabs(b->vr) : 0.0);
    const double speed_tolerance =
        convex_polygon_ccd_velocity_tolerance(pair_speed_bound);
    bool approaching = false;
    for (int contact_index = 0; contact_index < contact_count; ++contact_index)
        approaching |= normal_speed[contact_index] < -speed_tolerance;
    if (!approaching) {
        polygon_particles_collision_event(first);
        return;
    }
    double restitution = polygon_growing ? 1.0 : res;
    if (!polygon_growing && addExpo)
        restitution = resCoeff(maximum_closing_speed);
    restitution = minimum(1.0, maximum(0.0, restitution));
    double arm_normal_a[2] = {0, 0};
    double arm_normal_b[2] = {0, 0};
    double effective_mass[2][2] = {{0}};
    for (int row = 0; row < contact_count; ++row) {
        arm_normal_a[row] = cross(arms_a[row], contact.normal);
        arm_normal_b[row] = cross(arms_b[row], contact.normal);
    }
    for (int row = 0; row < contact_count; ++row) {
        for (int column = 0; column < contact_count; ++column) {
            effective_mass[row][column] = 1.0/a->m + 1.0/b->m +
                arm_normal_a[row]*arm_normal_a[column]/state_a->inertia +
                arm_normal_b[row]*arm_normal_b[column]/state_b->inertia;
        }
    }
    double impulses[2] = {0, 0};
    double constraints[2] = {0, 0};
    for (int contact_index = 0; contact_index < contact_count;
         ++contact_index) {
        if (normal_speed[contact_index] < -speed_tolerance) {
            constraints[contact_index] =
                -(1.0 + restitution)*normal_speed[contact_index];
        } else if (normal_speed[contact_index] > speed_tolerance) {
            /* A currently separating point may still be driven closed by a
               negatively coupled impulse at the other endpoint.  Its inactive
               constraint is u' >= 0, rather than simply being omitted. */
            constraints[contact_index] = -normal_speed[contact_index];
        }
    }
    if (!solve_contact_lcp(effective_mass, constraints, contact_count,
                           speed_tolerance, impulses)) {
        simulationReportFault(
            "Unable to solve the convex-particle contact constraints with finite repulsive impulses.");
        return;
    }
    if (!polygon_growing && !enforce_contact_energy_policy(
            effective_mass, normal_speed, contact_count, speed_tolerance,
            restitution == 1.0, impulses)) {
        simulationReportFault(
            "Unable to resolve the coupled convex-particle contact constraints without injecting energy.");
        return;
    }

    const double kinetic_before = polygon_growing
        ? body_kinetic_energy(a, state_a) +
          body_kinetic_energy(b, state_b) : 0.0;
    ConvexPolygonVec2 total_impulse_on_a = vec(0, 0);
    for (int contact_index = 0; contact_index < contact_count;
         ++contact_index) {
        if (!(impulses[contact_index] > 0)) continue;
        const ConvexPolygonVec2 impulse_on_a = scale(
            contact.normal, -impulses[contact_index]);
        const ConvexPolygonVec2 impulse_on_b = scale(
            contact.normal, impulses[contact_index]);
        total_impulse_on_a = add(total_impulse_on_a, impulse_on_a);
        a->vx += impulse_on_a.x/a->m;
        a->vy += impulse_on_a.y/a->m;
        b->vx += impulse_on_b.x/b->m;
        b->vy += impulse_on_b.y/b->m;
        state_a->omega += cross(arms_a[contact_index], impulse_on_a)/
            state_a->inertia;
        state_b->omega += cross(arms_b[contact_index], impulse_on_b)/
            state_b->inertia;
    }

    if (!update_growth_kinetic_energy(kinetic_before,
            body_kinetic_energy(a, state_a) +
            body_kinetic_energy(b, state_b))) return;

    if (!(length(total_impulse_on_a) > 0) ||
        !isfinite(total_impulse_on_a.x) ||
        !isfinite(total_impulse_on_a.y)) {
        simulationReportFault(
            "Unable to resolve a convex-particle contact with a finite repulsive impulse.");
        return;
    }

    collTermX += total_impulse_on_a.x*center_delta.x;
    collTermY += total_impulse_on_a.y*center_delta.y;
    collTermXY += total_impulse_on_a.x*center_delta.y;
    collTermYX += total_impulse_on_a.y*center_delta.x;
    a->lastColl = t;
    b->lastColl = t;
    a->coll++;
    b->coll++;
    ncol++;
    check_growth_progress();

    removeEventFromQueue(eventList[first]);
    removeEventFromQueue(eventList[second]);
    crossingEvent(first);
    crossingEvent(second);
    removeEventFromQueue(eventList[N + second]);
    collisionEvent(first);
    collisionEvent(second);
}

typedef enum WallImpulseStatus {
    WALL_IMPULSE_NONE = 0,
    WALL_IMPULSE_APPLIED,
    WALL_IMPULSE_OVERCONSTRAINED,
    WALL_IMPULSE_INVALID
} WallImpulseStatus;

static WallImpulseStatus apply_wall_impulse(
                               int index,
                               const PolygonWallManifold *manifold,
                               double *normal_impulse) {
    particle *value = particles + index;
    PolygonParticleState *state = polygon_state + index;
    const ConvexPolygonVec2 center = vec(value->x, value->y);
    if (manifold == NULL || manifold->count < 1 || manifold->count > 2 ||
        !(value->m > DBL_MIN) || !(state->inertia > DBL_MIN) ||
        !isfinite(value->m) || !isfinite(state->inertia))
        return WALL_IMPULSE_INVALID;

    ConvexPolygonVec2 arm[2];
    bool approaching = false;
    for (int contact = 0; contact < manifold->count; ++contact) {
        arm[contact] = subtract(manifold->point[contact], center);
    }
    const double wall_speed_bound = hypot(value->vx, value->vy) +
        fabs(state->omega)*value->rad +
        (polygon_growing ? fabs(value->vr) : 0.0);
    const double velocity_tolerance =
        convex_polygon_ccd_velocity_tolerance(wall_speed_bound);
    for (int contact = 0; contact < manifold->count; ++contact)
        approaching |= manifold->outward_speed[contact] > velocity_tolerance;
    if (!approaching) return WALL_IMPULSE_NONE;
    if (manifold->tied_count > 2)
        return WALL_IMPULSE_OVERCONSTRAINED;

    /* Preserve the established wall-tangent restitution as one resultant
       support-feature impulse.  For a planar face the averaged point gives the
       common tangential contact velocity; normal endpoint constraints are
       solved separately below. */
    ConvexPolygonVec2 tangent_point = manifold->point[0];
    ConvexPolygonVec2 tangent_normal = manifold->normal[0];
    if (manifold->count == 2) {
        tangent_point = scale(add(manifold->point[0], manifold->point[1]), 0.5);
        const ConvexPolygonVec2 normal_sum = add(manifold->normal[0],
                                                manifold->normal[1]);
        if (length(normal_sum) > DBL_MIN)
            tangent_normal = scale(normal_sum, 1.0/length(normal_sum));
    }
    const ConvexPolygonVec2 tangent = vec(-tangent_normal.y,
                                          tangent_normal.x);
    const ConvexPolygonVec2 tangent_arm = subtract(tangent_point, center);
    const ConvexPolygonVec2 tangent_velocity =
        velocity_at_point(value, state, tangent_arm);
    const double vt = dot(tangent_velocity, tangent);
    const double rt = cross(tangent_arm, tangent);
    const double ktt = 1.0/value->m + rt*rt/state->inertia;
    if (!(ktt > DBL_MIN) || !isfinite(ktt)) return WALL_IMPULSE_INVALID;
    const double tangent_restitution = polygon_growing ? 1.0 :
        minimum(1.0, maximum(-1.0, resT));
    const double jt = (tangent_restitution - 1.0)*vt/ktt;
    if (!isfinite(jt)) return WALL_IMPULSE_INVALID;
    const ConvexPolygonVec2 tangent_impulse = scale(tangent, jt);
    const double trial_vx = value->vx + tangent_impulse.x/value->m;
    const double trial_vy = value->vy + tangent_impulse.y/value->m;
    const double trial_omega = state->omega +
        cross(tangent_arm, tangent_impulse)/state->inertia;

    double normal_matrix[2][2] = {{0}};
    double normal_speed[2] = {0, 0};
    double constraints[2] = {0, 0};
    double impulses[2] = {0, 0};
    double arm_normal[2] = {0, 0};
    const double restitution = polygon_growing ? 1.0 :
        minimum(1.0, maximum(0.0, resW));
    for (int row = 0; row < manifold->count; ++row) {
        arm_normal[row] = cross(arm[row], manifold->normal[row]);
        ConvexPolygonVec2 trial_point_velocity = vec(
            trial_vx - trial_omega*arm[row].y,
            trial_vy + trial_omega*arm[row].x);
        if (polygon_growing && value->rad > DBL_MIN)
            trial_point_velocity = add(trial_point_velocity,
                scale(arm[row], value->vr/value->rad));
        normal_speed[row] = dot(trial_point_velocity, manifold->normal[row]);
        if (normal_speed[row] > velocity_tolerance)
            constraints[row] = (1.0 + restitution)*normal_speed[row];
        else if (normal_speed[row] < -velocity_tolerance)
            constraints[row] = normal_speed[row];
    }
    for (int row = 0; row < manifold->count; ++row)
        for (int column = 0; column < manifold->count; ++column)
            normal_matrix[row][column] =
                dot(manifold->normal[row], manifold->normal[column])/value->m +
                arm_normal[row]*arm_normal[column]/state->inertia;
    if (!solve_contact_lcp(normal_matrix, constraints, manifold->count,
                           velocity_tolerance, impulses))
        return WALL_IMPULSE_INVALID;
    if (!polygon_growing) {
        double impulse_velocity[2] = {
            -normal_speed[0], -normal_speed[1]
        };
        if (!enforce_contact_energy_policy(
                normal_matrix, impulse_velocity, manifold->count,
                velocity_tolerance, restitution == 1.0, impulses))
            return WALL_IMPULSE_INVALID;
    }

    ConvexPolygonVec2 total_impulse = tangent_impulse;
    double angular_impulse = cross(tangent_arm, tangent_impulse);
    double total_normal_impulse = 0;
    for (int contact = 0; contact < manifold->count; ++contact) {
        if (!(impulses[contact] > 0)) continue;
        const ConvexPolygonVec2 contact_impulse =
            scale(manifold->normal[contact], -impulses[contact]);
        total_impulse = add(total_impulse, contact_impulse);
        angular_impulse += cross(arm[contact], contact_impulse);
        total_normal_impulse += impulses[contact];
    }
    if (!isfinite(total_impulse.x) || !isfinite(total_impulse.y) ||
        !isfinite(angular_impulse) || !isfinite(total_normal_impulse))
        return WALL_IMPULSE_INVALID;
    if (!(fabs(jt) + total_normal_impulse > 0))
        return WALL_IMPULSE_INVALID;
    value->vx += total_impulse.x/value->m;
    value->vy += total_impulse.y/value->m;
    state->omega += angular_impulse/state->inertia;
    if (normal_impulse != NULL) *normal_impulse = total_normal_impulse;
    return WALL_IMPULSE_APPLIED;
}

void polygon_particles_do_wall(void) {
    const int index = nextEvent->i;
    const int wall = nextEvent->j;
    if (!polygon_particles_active() || index < 0 || index >= N) return;
    particle *value = particles + index;
    polygon_particles_free_fly(value);
    PolygonWallManifold manifold;
    if (!wall_manifold_at(index, t, wall, &manifold)) {
        simulationReportFault("Unable to evaluate a convex-particle wall contact.");
        return;
    }
    const double penetration_tolerance = 2.0e-8*(1.0 + value->rad + Lx + Ly);
    if (manifold.gap < -penetration_tolerance) {
        report_overlap(index, -1,
            wall == POLYGON_WALL_CIRCULAR ? "circular wall" : "planar wall");
        return;
    }
    if (manifold.gap < 0) {
        value->x += manifold.gap*manifold.deepest_normal.x;
        value->y += manifold.gap*manifold.deepest_normal.y;
        if (!wall_manifold_at(index, t, wall, &manifold)) {
            simulationReportFault("Unable to reevaluate a corrected convex-particle wall contact.");
            return;
        }
        if (manifold.gap < -penetration_tolerance) {
            report_overlap(index, -1,
                wall == POLYGON_WALL_CIRCULAR ? "circular wall" : "planar wall");
            return;
        }
    }
    const double kinetic_before = polygon_growing
        ? body_kinetic_energy(value, polygon_state + index) : 0.0;
    double impulse = 0;
    const WallImpulseStatus impulse_status =
        apply_wall_impulse(index, &manifold, &impulse);
    if (impulse_status == WALL_IMPULSE_OVERCONSTRAINED) {
        simulationReportFault(
            "A convex particle touches the circular wall at more than two approaching vertices and has no uniquely resolvable free motion. Reduce the packing fraction and restart.");
        return;
    }
    if (impulse_status == WALL_IMPULSE_INVALID) {
        simulationReportFault(
            "Unable to resolve a convex-particle wall contact with finite repulsive impulses.");
        return;
    }
    if (impulse_status == WALL_IMPULSE_APPLIED) {
        if (!update_growth_kinetic_energy(kinetic_before,
                body_kinetic_energy(value, polygon_state + index))) return;
        dp += impulse;
        cumulativeWallImpulse += impulse;
        if (wall == POLYGON_WALL_MID_LEFT || wall == POLYGON_WALL_MID_RIGHT)
            dpMid += impulse;
        else if (wall == POLYGON_WALL_LEFT) dpLeft += impulse;
        else if (wall == POLYGON_WALL_RIGHT) dpRight += impulse;
        ncol++;
        check_growth_progress();
    }
    if (impulse_status == WALL_IMPULSE_NONE) {
        /* A stale same-time event may become separating after another contact.
           Rebuild it without changing the collision counter or invalidating
           unrelated predictions owned by other particles. */
        collisionEvent(index);
        return;
    }
    value->lastColl = t;
    value->coll++;
    removeEventFromQueue(eventList[index]);
    crossingEvent(index);
    collisionEvent(index);
}

void polygon_particles_do_recheck(void) {
    const int index = nextEvent->i;
    if (!polygon_particles_active() || index < 0 || index >= N) return;
    polygon_particles_free_fly(particles + index);
    polygon_particles_collision_event(index);
}
