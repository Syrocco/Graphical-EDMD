#ifndef EDMD_POLYGON_PARTICLES_H
#define EDMD_POLYGON_PARTICLES_H

#include "EDMD.h"
#include "sim_config.h"

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Polygon CCD is expensive enough that a finer linked-cell grid is useful.
 * This bound keeps the periodic coordinate lists fixed-size and allocation
 * free; the core disk/potential stencil remains independently capped at 2. */
#define POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE 12
#define POLYGON_PARTICLES_MAX_CELL_COORDINATES \
    (2*POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE + 1)

/*
 * Optional two-dimensional rigid convex-particle adapter.
 *
 * The ordinary particle structure deliberately contains no polygon data.
 * Orientation, spin, inertia and the shared canonical shape are allocated
 * only when this model is selected, leaving the analytical disk hot path and
 * its cache footprint unchanged.
 */
void polygon_particles_get_config(SimParticleShapeConfig *config);
void polygon_particles_set_config(const SimParticleShapeConfig *config);

bool polygon_particles_requested(void);
bool polygon_particles_active(void);
bool polygon_particles_growing(void);

/* Construct/release the per-particle sidecar for the current particle array. */
bool polygon_particles_prepare(void);
void polygon_particles_release(void);

/* Select the rigid-polygon normal-dynamics handlers after initialization or
 * after homothetic growth of the actual convex bodies has completed. */
void polygon_particles_activate(void);
void polygon_particles_begin_growth(void);
void polygon_particles_finish_growth(void);
/* Recompute and validate the preparation kinetic energy after a non-contact
 * velocity update (for example a thermostat kick or interactive drag).
 * Contact impulses update the cached aggregate incrementally. */
bool polygon_particles_refresh_growth_energy(void);

/* Change one particle's mass without leaving the rigid-body sidecar stale.
 * Angular velocity is intentionally preserved, so the instantaneous motion
 * remains continuous; inertia and any active growth-energy cache are updated. */
bool polygon_particles_set_mass(int index, double mass);

/* EDMD function-pointer adapters. */
void polygon_particles_crossing_event(int index);
void polygon_particles_collision_event(int index);
void polygon_particles_do_collision(void);
void polygon_particles_do_wall(void);
void polygon_particles_do_recheck(void);
void polygon_particles_free_fly(particle *value);

/* Rendering/analysis accessors. Angles are predicted without synchronizing
 * the simulation record. The shared vertices are centroid-zero and have
 * circumradius one. */
int polygon_particles_vertex_count(void);
bool polygon_particles_vertex(int vertex, double *x, double *y);
double polygon_particles_unit_area(void);
double polygon_particles_angle_at(int index, double time);
double polygon_particles_omega(int index);
double polygon_particles_inertia(int index);
double polygon_particles_rotational_energy(void);
bool polygon_particles_apply_angular_langevin(int index,
                                              double temperature,
                                              double gamma,
                                              double interval,
                                              bool euler);
bool polygon_particles_contains_point(int index, double x, double y,
                                      double time);

/* Value-only sidecar checkpoint used by the graphical numerical-recovery
 * ring. The particle array itself remains owned by the core engine. */
size_t polygon_particles_checkpoint_bytes(void);
bool polygon_particles_capture_checkpoint(void *destination, size_t bytes);
bool polygon_particles_restore_checkpoint(const void *source, size_t bytes);

/* Temporarily spend more CCD iterations while replaying a GUI checkpoint.
 * Headless runs never call this and retain the ordinary fast path. */
void polygon_particles_enable_recovery_precision(double duration);

/* Exact, read-only invariant check at a lazy simulation time. Every active
 * planar/circular wall vertex gap and every minimum-image pair SAT separation
 * is tested against the caller's absolute penetration tolerance. No event,
 * particle, or fault state is changed. On failure, message describes either
 * the overlap or invalid input/geometry (and is empty after success). */
bool polygon_particles_verify_no_overlaps(double time,
                                          double penetration_tolerance,
                                          char *message,
                                          size_t message_capacity);

/* Current-time safety check using the active linked-cell stencil.  This is
 * O(N * local-neighbours) and is intended for frequent GUI checkpoint/
 * rollback validation.  The exact all-pairs routine above remains available
 * for tests and one-shot diagnostics at arbitrary lazy times. */
bool polygon_particles_verify_local_no_overlaps(
    double time,
    double penetration_tolerance,
    char *message,
    size_t message_capacity);

/* Optional profiling counters.  Define EDMD_POLYGON_SCHEDULER_STATS=1 for a
 * benchmark build; ordinary builds return an all-zero snapshot and compile
 * the hot-path increments out completely. */
typedef struct PolygonParticleSchedulerStats {
    unsigned long long collision_event_calls;
    unsigned long long pair_candidates;
    unsigned long long pair_broadphase_rejections;
    unsigned long long pair_ccd_calls;
    unsigned long long pair_ccd_iterations;
    unsigned long long wall_candidates;
    unsigned long long wall_broadphase_rejections;
    unsigned long long wall_manifold_evaluations;
    unsigned long long rechecks_scheduled;
} PolygonParticleSchedulerStats;

void polygon_particles_reset_scheduler_stats(void);
void polygon_particles_get_scheduler_stats(
    PolygonParticleSchedulerStats *stats);

#if EDMD_TEST_PREPARED_STATE
/* Focused regression hook for the exact two-contact active-set solver. */
bool polygon_particles_test_solve_contact_lcp(
    const double matrix[2][2], const double bound[2], int count,
    double velocity_tolerance, double impulse[2]);
/* Elastic energy projection used after a successful manifold LCP. */
bool polygon_particles_test_project_elastic_contact_energy(
    const double matrix[2][2], const double velocity[2], int count,
    double velocity_tolerance, double impulse[2]);
bool polygon_particles_test_enforce_contact_energy_policy(
    const double matrix[2][2], const double velocity[2], int count,
    double velocity_tolerance, bool elastic, double impulse[2]);
/* Convert a local CCD delay to the global calendar clock with the same
 * representability guard used by collision and recheck scheduling. */
bool polygon_particles_test_event_deadline_after(double delay,
                                                 double *deadline);
bool polygon_particles_test_recovery_precision_active(void);
double polygon_particles_test_growth_energy_ceiling(void);
int polygon_particles_test_neighboring_coordinates(
    int center, int count,
    int output[POLYGON_PARTICLES_MAX_CELL_COORDINATES]);
#endif

#ifdef __cplusplus
}
#endif

#endif
