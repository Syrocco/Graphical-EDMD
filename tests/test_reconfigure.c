#include "EDMD.h"
#include "graphics.h"
#include "sim_config.h"
#include "sim_lifecycle.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if G != 1
#error "test_reconfigure requires mutable graphical globals (G=1)"
#endif

/* EDMD.c is compiled as a graphical core so its configuration globals remain
 * mutable.  The regression itself is headless; these small host/raylib stubs
 * keep window construction and rendering outside the lifecycle test. */
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

extern int N, Nbig, Nsmall;
extern int Nxcells, Nycells;
#if THREE_D
extern int Nzcells;
#endif
extern int addWallx, addWally, addWallz, addMidWall, addCircularWall, noise;
extern bool addWell, addField, addShear;
extern double t, paulTime, Lx, U, sig, res, dtnoise, delta;
extern unsigned long int ncol, ncross;
extern particle *particles;
extern particle **cellList;
extern node **eventList;
extern node *nextEvent;

enum { TEST_PARTICLES = 2 };

static int failures;

static void fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    failures++;
}

static void expect_true(int condition, const char *message)
{
    if (!condition) fail(message);
}

static void expect_close(double actual, double expected, const char *message)
{
    const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
    if (!isfinite(actual) || fabs(actual - expected) > 1e-12*scale) {
        fprintf(stderr, "%s: got %.17g, expected %.17g\n",
                message, actual, expected);
        failures++;
    }
}

typedef struct PreservedState {
    particle *allocation;
    int count;
    double time;
    double paul_time;
    unsigned long collisions;
    unsigned long crossings;
    particle values[TEST_PARTICLES];
} PreservedState;

static PreservedState capture_state(void)
{
    PreservedState result;
    memset(&result, 0, sizeof(result));
    result.allocation = particles;
    result.count = N;
    result.time = t;
    result.paul_time = paulTime;
    result.collisions = ncol;
    result.crossings = ncross;
    for (int index = 0; index < TEST_PARTICLES; ++index)
        result.values[index] = particles[index];
    return result;
}

static void expect_particle_state_preserved(const PreservedState *before,
                                            const char *context)
{
    if (particles != before->allocation || N != before->count ||
        t != before->time || paulTime != before->paul_time ||
        ncol != before->collisions || ncross != before->crossings) {
        fprintf(stderr, "%s replaced/reset the running state\n", context);
        failures++;
    }
    for (int index = 0; index < TEST_PARTICLES; ++index) {
        const particle *actual = particles + index;
        const particle *expected = before->values + index;
        if (actual->x != expected->x || actual->y != expected->y ||
            actual->vx != expected->vx || actual->vy != expected->vy ||
            actual->rad != expected->rad || actual->m != expected->m ||
            actual->t != expected->t || actual->coll != expected->coll ||
            actual->crossX != expected->crossX ||
            actual->crossY != expected->crossY
#if THREE_D
            || actual->z != expected->z || actual->vz != expected->vz ||
            actual->crossZ != expected->crossZ
#endif
        ) {
            fprintf(stderr, "%s changed particle %d\n", context, index);
            failures++;
        }
    }
}

static SimBoundaryMode closed_boundary(void)
{
#if THREE_D
    return SIM_BOUNDARY_WALL_XYZ;
#else
    return SIM_BOUNDARY_WALL_XY;
#endif
}

static SimulationConfig prepare_state_with_membrane(
    SimBoundaryMode boundary, bool potential_enabled, bool center_membrane)
{
    BuildCapabilities capabilities = simulation_build_capabilities();
    SimulationConfig config;

    if (particles != NULL || cellList != NULL || eventList != NULL)
        freeArrays();
    simulation_config_set_defaults(&config, &capabilities);
    config.initialization.particle_count = TEST_PARTICLES;
    config.initialization.packing_fraction = 0.1;
    config.initialization.box_lx = 20.0;
    config.initialization.box_ly = 20.0;
    config.initialization.box_lz = 20.0;
    config.initialization.automatic_growth_rate = false;
    config.initialization.growth_rate = 0.01;
    config.geometry.boundary = boundary;
    config.geometry.center_membrane = center_membrane;
    config.potential.enabled = potential_enabled;
    config.potential.energy = 0.15;
    config.potential.range = 1.5;
    config.output.periodic_dumps = false;
    config.observables.thermo_enabled = false;

    simulation_config_apply_globals(&config);
    particles = calloc(TEST_PARTICLES, sizeof(*particles));
    if (particles == NULL) {
        fail("unable to allocate reconfiguration particles");
        return config;
    }

    t = 7.25;
    paulTime = t;
    ncol = 17;
    ncross = 23;
    Nbig = TEST_PARTICLES;
    Nsmall = 0;
    for (int index = 0; index < TEST_PARTICLES; ++index) {
        particle *value = particles + index;
        value->num = index;
        value->type = 1;
        value->m = 1.0;
        value->rad = 1.0;
        value->t = t;
        value->x = index == 0 ? 5.0 : 15.0;
        value->y = index == 0 ? 5.0 : 15.0;
        value->vx = 0.0;
        value->vy = index == 0 ? 0.4 : -0.3;
#if THREE_D
        value->z = index == 0 ? 5.0 : 15.0;
        value->vz = index == 0 ? 0.2 : -0.15;
#endif
    }
    adaptCellGridToParticles();
    cellListInit();
    if (potential_enabled) pairsInit();
#if (STRESS || HEAT)
    initializeTransportAccumulators();
#endif
    eventListInit();
    return config;
}

static SimulationConfig prepare_state(SimBoundaryMode boundary,
                                      bool potential_enabled)
{
    return prepare_state_with_membrane(boundary, potential_enabled, false);
}

static SimulationReconfigureAction reconfigure(
    const SimulationConfig *current, const SimulationConfig *next,
    SimulationReconfigureAction expected, const char *context,
    char message[512])
{
    message[0] = '\0';
    const SimulationReconfigureAction actual = simulation_reconfigure_core(
        current, next, message, 512);
    if (actual != expected) {
        fprintf(stderr, "%s: action %d, expected %d (%s)\n",
                context, (int)actual, (int)expected,
                message[0] != '\0' ? message : "no diagnostic");
        failures++;
    }
    return actual;
}

static void test_response_only_restitution(void)
{
    SimulationConfig current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    SimulationConfig next = current;
    char message[512];
    const PreservedState before = capture_state();
    node *crossing = eventList[0];
    node *collision = eventList[N];
    const double crossing_time = crossing->t;
    const double collision_time = collision->t;

    next.dynamics.particle_restitution = 0.72;
    reconfigure(&current, &next, SIM_RECONFIGURE_GLOBALS_ONLY,
                "restitution response update", message);
    expect_close(res, 0.72, "restitution global was not updated");
    expect_true(eventList[0] == crossing && eventList[N] == collision &&
                    crossing->t == crossing_time &&
                    collision->t == collision_time,
                "restitution update regenerated event predictions");
    expect_particle_state_preserved(&before, "restitution update");
}

static double execute_head_on_collision(double elapsed_since_last_collision,
                                        double *flight_duration)
{
    const double start_time = t;
    particle *left = particles;
    particle *right = particles + 1;

    /* Start in adjacent linked cells so the pair prediction is immediate, but
       leave a small positive flight interval before contact. */
    left->x = 8.9;
    right->x = 11.1;
    left->y = right->y = 10.0;
    left->vx = 1.0;
    right->vx = -1.0;
    left->vy = right->vy = 0.0;
    left->t = right->t = start_time;
    left->lastColl = right->lastColl =
        start_time - elapsed_since_last_collision;
#if THREE_D
    left->z = right->z = 10.0;
    left->vz = right->vz = 0.0;
#endif

    rebuildSpatialIndexAtCurrentTime();
    rescheduleParticleEventsAtCurrentTime();
#if (STRESS || HEAT)
    initializeTransportAccumulators();
#endif

    node *collision = eventList[N];
    if (collision == NULL || collision->type != COLLISION ||
        collision->j != 1 || !isfinite(collision->t)) {
        fail("head-on injection regression did not schedule pair 0/1");
        return NAN;
    }
    if (flight_duration != NULL)
        *flight_duration = collision->t - start_time;
    nextEvent = collision;
    t = collision->t;
    removeEventFromQueue(collision);
    doTheCollisionNormal();
    nextEvent = NULL;
    return left->vx;
}

static void expect_declarative_delta(double expected, const char *context)
{
    SimulationConfig captured;
    simulation_config_capture_globals(&captured);
    expect_close(delta, expected, context);
    expect_close(captured.dynamics.delta, expected,
                 "captured configuration adopted a per-collision delta");
}

static void test_collision_injection_keeps_configuration_declarative(void)
{
    char message[512];
    SimulationConfig current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    SimulationConfig next = current;

    next.dynamics.injection = SIM_INJECTION_DOUBLE_DELTA_NORMAL;
    next.dynamics.delta = 0.37;
    next.dynamics.delta_small = 0.14;
    next.dynamics.delta_large = 0.62;
    next.dynamics.collision_memory_time = 10.0;
    reconfigure(&current, &next, SIM_RECONFIGURE_GLOBALS_ONLY,
                "double-delta response update", message);
    const double double_velocity = execute_head_on_collision(0.0, NULL);
    expect_close(double_velocity,
                 -next.dynamics.particle_restitution -
                    next.dynamics.delta_large,
                 "double-delta collision used the wrong local amplitude");
    expect_declarative_delta(0.37,
        "double-delta collision overwrote configured delta");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    next = current;
    next.dynamics.injection = SIM_INJECTION_EVOLVING_DELTA;
    next.dynamics.delta = 0.41;
    next.dynamics.evolving_delta_tau = 2.5;
    next.dynamics.evolving_delta_max = 0.80;
    reconfigure(&current, &next, SIM_RECONFIGURE_GLOBALS_ONLY,
                "evolving-delta response update", message);
    double flight_duration = 0.0;
    const double evolving_velocity = execute_head_on_collision(
        3.0, &flight_duration);
    const double elapsed_at_collision = 3.0 + flight_duration;
    const double evolving_amplitude = next.dynamics.evolving_delta_max*
        pow(1.0 - exp(-elapsed_at_collision/
                      next.dynamics.evolving_delta_tau), 3.0);
    expect_close(evolving_velocity,
                 -next.dynamics.particle_restitution - evolving_amplitude,
                 "evolving-delta collision used the wrong local amplitude");
    expect_declarative_delta(0.41,
        "evolving-delta collision overwrote configured delta");
}

static void test_wall_topology(void)
{
    char message[512];

    SimulationConfig current = prepare_state(closed_boundary(), false);
    SimulationConfig next = current;
    const PreservedState removal_before = capture_state();
    next.geometry.boundary = SIM_BOUNDARY_PERIODIC;
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "wall removal", message);
    expect_true(!addWallx && !addWally && !addWallz,
                "wall removal did not restore periodic axes");
    expect_particle_state_preserved(&removal_before, "wall removal");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    next = current;
    const PreservedState addition_before = capture_state();
    next.geometry.boundary = closed_boundary();
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "contained wall addition", message);
    expect_true(addWallx && addWally
#if THREE_D
                    && addWallz
#endif
                , "wall addition did not activate every closed-box face");
    expect_particle_state_preserved(&addition_before, "wall addition");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    next = current;
#if THREE_D
    particles[0].z = 0.5;
#else
    particles[0].x = 0.5;
#endif
    const PreservedState rejected_before = capture_state();
    next.geometry.boundary = closed_boundary();
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "uncontained wall addition", message);
    expect_true(message[0] != '\0',
                "wall rebuild request gave no diagnostic");
    expect_true(!addWallx && !addWally && !addWallz,
                "wall rebuild preflight partially changed topology");
    expect_particle_state_preserved(&rejected_before,
                                    "wall addition awaiting rebuild");
    expect_true(simulation_restart_core(&next, message, sizeof(message)),
                "wall addition could not rebuild the requested system");
    expect_true(addWallx && addWally
#if THREE_D
                    && addWallz
#endif
                , "wall fallback rebuild did not activate the requested faces");
    expect_true(t == 0.0 && !simulationHasFault(),
                "wall fallback rebuild did not produce a clean t=0 state");

    current = prepare_state(closed_boundary(), true);
    next = current;
    const PreservedState potential_wall_before = capture_state();
    particle **potential_wall_cells = cellList;
    next.geometry.boundary = SIM_BOUNDARY_PERIODIC;
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "potential wall removal", message);
    expect_true(cellList == potential_wall_cells,
                "removing walls rebuilt an unchanged potential cell grid");
    expect_particle_state_preserved(&potential_wall_before,
                                    "potential wall removal");
}

static void test_spherical_wall_fallback(void)
{
    char message[512];
    SimulationConfig current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    SimulationConfig next = current;

    /* Put one particle outside the requested curved confinement.  The live
       trajectory cannot be retained, but the configuration itself is valid
       and must therefore request a clean reconstruction instead of rollback. */
    particles[0].x = 0.5;
    particles[0].y = 0.5;
#if THREE_D
    particles[0].z = 0.5;
#endif
    const PreservedState before = capture_state();
    next.geometry.boundary = SIM_BOUNDARY_CIRCULAR_SPHERICAL;
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "uncontained spherical-wall addition", message);
    expect_true(message[0] != '\0',
                "spherical-wall rebuild request gave no diagnostic");
    expect_true(!addCircularWall,
                "spherical-wall preflight partially changed the topology");
    expect_particle_state_preserved(&before,
                                    "spherical-wall addition awaiting rebuild");

    expect_true(simulation_restart_core(&next, message, sizeof(message)),
                "spherical-wall addition could not rebuild the requested system");
    expect_true(addCircularWall && t == 0.0 && !simulationHasFault(),
                "spherical-wall fallback did not produce a clean confined state");
}

static void test_center_membrane_topology(void)
{
#if THREE_D
    /* The planar membrane and its pressure diagnostics are deliberately a
       2D-only capability.  Its validator rejects enabling it in 3D. */
    return;
#else
    char message[512];

    SimulationConfig current = prepare_state_with_membrane(
        SIM_BOUNDARY_PERIODIC, false, true);
    SimulationConfig next = current;
    const PreservedState removal_before = capture_state();
    particle **removal_cells = cellList;
    const int removal_nx = Nxcells;
    const int removal_ny = Nycells;
    next.geometry.center_membrane = false;
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "center-membrane removal", message);
    expect_true(!addMidWall,
                "center-membrane removal left the separator active");
    expect_true(cellList == removal_cells && Nxcells == removal_nx &&
                    Nycells == removal_ny,
                "center-membrane removal rebuilt the spatial index");
    expect_particle_state_preserved(&removal_before,
                                    "center-membrane removal");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    next = current;
    const PreservedState addition_before = capture_state();
    particle **addition_cells = cellList;
    next.geometry.center_membrane = true;
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "contained center-membrane addition", message);
    expect_true(addMidWall,
                "contained center-membrane addition was not activated");
    expect_true(cellList == addition_cells,
                "center-membrane addition rebuilt the spatial index");
    expect_particle_state_preserved(&addition_before,
                                    "center-membrane addition");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    next = current;
    particles[0].x = 0.5*Lx;
    particles[0].type = 1;
    const PreservedState rejected_before = capture_state();
    particle **rejected_cells = cellList;
    node *rejected_crossing = eventList[0];
    node *rejected_collision = eventList[N];
    const double rejected_crossing_time = rejected_crossing->t;
    const double rejected_collision_time = rejected_collision->t;
    next.geometry.center_membrane = true;
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "intersecting center-membrane addition", message);
    expect_true(message[0] != '\0',
                "center-membrane rebuild request gave no diagnostic");
    expect_true(!addMidWall && cellList == rejected_cells &&
                    eventList[0] == rejected_crossing &&
                    eventList[N] == rejected_collision &&
                    rejected_crossing->t == rejected_crossing_time &&
                    rejected_collision->t == rejected_collision_time,
                "center-membrane rebuild preflight changed runtime topology");
    expect_particle_state_preserved(&rejected_before,
                                    "center-membrane addition awaiting rebuild");
    expect_true(simulation_restart_core(&next, message, sizeof(message)),
                "center-membrane addition could not rebuild the requested system");
    expect_true(addMidWall && t == 0.0 && !simulationHasFault(),
                "center-membrane fallback rebuild did not produce the requested clean state");
#endif
}

static void test_potential_tiers(void)
{
    char message[512];

    SimulationConfig current = prepare_state(SIM_BOUNDARY_PERIODIC, true);
    SimulationConfig next = current;
    const PreservedState scalar_before = capture_state();
    node *crossing = eventList[0];
    node *collision = eventList[N];
    const double crossing_time = crossing->t;
    const double collision_time = collision->t;
    const int scalar_cells = Nxcells;
    next.potential.energy = 0.30;
    reconfigure(&current, &next, SIM_RECONFIGURE_GLOBALS_ONLY,
                "nonzero potential-energy scalar", message);
    expect_close(U, 0.30, "potential energy global was not updated");
    expect_true(Nxcells == scalar_cells && eventList[0] == crossing &&
                    eventList[N] == collision &&
                    crossing->t == crossing_time &&
                    collision->t == collision_time,
                "potential-energy scalar rebuilt structural state");
    expect_particle_state_preserved(&scalar_before,
                                    "potential-energy scalar update");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, true);
    next = current;
    const PreservedState same_grid_range_before = capture_state();
    particle **same_grid_cells = cellList;
    const int same_grid_nx = Nxcells;
    next.potential.range = 1.51;
    reconfigure(&current, &next, SIM_RECONFIGURE_SPATIAL_REBUILT,
                "same-grid potential range update", message);
    expect_true(cellList == same_grid_cells && Nxcells == same_grid_nx,
                "small range edit reallocated an unchanged cell grid");
    expect_particle_state_preserved(&same_grid_range_before,
                                    "same-grid potential range update");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, true);
    next = current;
    const PreservedState range_before = capture_state();
    const int old_range_cells = Nxcells;
    next.potential.range = 1.8;
    reconfigure(&current, &next, SIM_RECONFIGURE_SPATIAL_REBUILT,
                "potential range update", message);
    expect_true(Nxcells != old_range_cells,
                "potential range did not rebuild the cutoff grid");
    expect_particle_state_preserved(&range_before, "potential range update");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    next = current;
    const PreservedState enable_before = capture_state();
    const int hard_cells = Nxcells;
    next.potential.enabled = true;
    reconfigure(&current, &next, SIM_RECONFIGURE_SPATIAL_REBUILT,
                "potential enable", message);
    expect_true(addWell && Nxcells != hard_cells,
                "enabling the potential did not rebuild spatial state");
    expect_particle_state_preserved(&enable_before, "potential enable");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    current.potential.energy = 0.0;
    simulation_config_apply_globals(&current);
    next = current;
    const PreservedState zero_enable_before = capture_state();
    particle **zero_enable_cells = cellList;
    node *zero_enable_crossing = eventList[0];
    node *zero_enable_collision = eventList[N];
    const double zero_enable_crossing_time = zero_enable_crossing->t;
    const double zero_enable_collision_time = zero_enable_collision->t;
    const int zero_enable_nx = Nxcells;
    const int zero_enable_ny = Nycells;
#if THREE_D
    const int zero_enable_nz = Nzcells;
#endif
    next.potential.enabled = true;
    reconfigure(&current, &next, SIM_RECONFIGURE_GLOBALS_ONLY,
                "zero-strength potential enable", message);
    expect_true(addWell && U == 0.0 && cellList == zero_enable_cells &&
                    Nxcells == zero_enable_nx && Nycells == zero_enable_ny
#if THREE_D
                    && Nzcells == zero_enable_nz
#endif
                    && eventList[0] == zero_enable_crossing &&
                    eventList[N] == zero_enable_collision &&
                    zero_enable_crossing->t == zero_enable_crossing_time &&
                    zero_enable_collision->t == zero_enable_collision_time,
                "zero-strength potential enable performed structural work");
    expect_particle_state_preserved(&zero_enable_before,
                                    "zero-strength potential enable");

    current = next;
    next = current;
    const PreservedState zero_range_before = capture_state();
    particle **zero_range_cells = cellList;
    node *zero_range_crossing = eventList[0];
    node *zero_range_collision = eventList[N];
    const double zero_range_crossing_time = zero_range_crossing->t;
    const double zero_range_collision_time = zero_range_collision->t;
    const int zero_range_nx = Nxcells;
    const int zero_range_ny = Nycells;
#if THREE_D
    const int zero_range_nz = Nzcells;
#endif
    next.potential.range = 1.8;
    reconfigure(&current, &next, SIM_RECONFIGURE_GLOBALS_ONLY,
                "zero-strength potential range", message);
    expect_close(sig, 1.8, "zero-strength potential range was not retained");
    expect_true(cellList == zero_range_cells &&
                    Nxcells == zero_range_nx && Nycells == zero_range_ny
#if THREE_D
                    && Nzcells == zero_range_nz
#endif
                    && eventList[0] == zero_range_crossing &&
                    eventList[N] == zero_range_collision &&
                    zero_range_crossing->t == zero_range_crossing_time &&
                    zero_range_collision->t == zero_range_collision_time,
                "zero-strength potential range performed structural work");
    expect_particle_state_preserved(&zero_range_before,
                                    "zero-strength potential range");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, true);
    particles[0].x = 8.9;
    particles[1].x = 11.1;
    particles[0].y = particles[1].y = 10.0;
    particles[0].vx = particles[1].vx = 0.0;
    particles[0].vy = particles[1].vy = 0.0;
    particles[0].type = 0;
    particles[1].type = 1;
    current.initialization.composition = SIM_COMPOSITION_BINARY;
    current.initialization.small_fraction = 0.5;
#if THREE_D
    particles[0].z = particles[1].z = 10.0;
    particles[0].vz = particles[1].vz = 0.0;
#endif
    rebuildSpatialIndexAtCurrentTime();
    rescheduleParticleEventsAtCurrentTime();
    expect_true(isParticleInWellList(particles, 1),
                "all-pairs potential omitted an inside membership");
    next = current;
    const PreservedState species_before = capture_state();
    particle **species_cells = cellList;
    const int species_nx = Nxcells;
    const int species_ny = Nycells;
#if THREE_D
    const int species_nz = Nzcells;
#endif
    next.potential.species_mode = SIM_POTENTIAL_SAME_SPECIES;
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "potential species update", message);
    expect_true(cellList == species_cells && Nxcells == species_nx &&
                    Nycells == species_ny
#if THREE_D
                    && Nzcells == species_nz
#endif
                    , "species update unnecessarily rebuilt the cell grid");
    expect_true(!isParticleInWellList(particles, 1),
                "inactive cross-species pair retained a well membership");
    expect_true(eventList[N]->type == COLLISION,
                "inactive cross-species pair retained square-well events");
    expect_particle_state_preserved(&species_before,
                                    "potential species update");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    next = current;
    const PreservedState inactive_charge_before = capture_state();
    particle **inactive_charge_cells = cellList;
    node *inactive_charge_event = eventList[N];
    next.potential.species_mode = SIM_POTENTIAL_CHARGED;
    next.potential.positive_charge_fraction = 0.5;
    next.potential.neutral_charge_fraction = 0.5;
    reconfigure(&current, &next, SIM_RECONFIGURE_GLOBALS_ONLY,
                "inactive charged-potential metadata", message);
    expect_true(cellList == inactive_charge_cells &&
                    eventList[N] == inactive_charge_event &&
                    particles[0].charge == 0 && particles[1].charge == 0,
                "inactive charge fractions touched runtime structures");
    expect_particle_state_preserved(&inactive_charge_before,
                                    "inactive charge metadata");

    current = next;
    next.potential.enabled = true;
    reconfigure(&current, &next, SIM_RECONFIGURE_SPATIAL_REBUILT,
                "charged-potential enable", message);
    expect_true((particles[0].charge == 1 && particles[1].charge == 0) ||
                    (particles[0].charge == 0 && particles[1].charge == 1),
                "enabling a remembered charged potential did not realize its fractions");
}

static void configure_changed_flight(SimulationConfig *next)
{
#if (STRESS || HEAT)
    /* Transport builds intentionally expose only ballistic and Lees-Edwards
       flight because those are the models with complete moment accounting. */
    next->dynamics.flight = SIM_FLIGHT_LEES_EDWARDS_SHEAR;
    next->dynamics.shear_strain_rate = 0.02;
#else
    next->dynamics.flight = SIM_FLIGHT_CONSTANT_FIELD;
    next->dynamics.constant_field = -0.01;
#endif
}

static void test_lazy_particle_flight_synchronization(void)
{
    char message[512];
    SimulationConfig current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    SimulationConfig next = current;
    particle *lazy = particles;
    const double lag = 0.75;
    const double expected_x = lazy->x;
    const double expected_y = lazy->y;
    const double expected_vx = lazy->vx;
    const double expected_vy = lazy->vy;
#if THREE_D
    const double expected_z = lazy->z;
    const double expected_vz = lazy->vz;
#endif
    const particle other_before = particles[1];
    particle *allocation_before = particles;
    const unsigned long collisions_before = ncol;
    const unsigned long crossings_before = ncross;

    /* Store particle zero at an earlier lazy timestamp while leaving its cell
       membership at the position it reaches at the current global time. */
    lazy->t = t - lag;
    lazy->x = expected_x - lag*expected_vx;
    lazy->y = expected_y - lag*expected_vy;
#if THREE_D
    lazy->z = expected_z - lag*expected_vz;
#endif
    configure_changed_flight(&next);
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "lazy-particle flight update", message);
    expect_close(lazy->t, t,
                 "flight update did not synchronize the lazy timestamp");
    expect_close(lazy->x, expected_x,
                 "flight update advanced lazy x with the new law");
    expect_close(lazy->y, expected_y,
                 "flight update advanced lazy y with the new law");
    expect_close(lazy->vx, expected_vx,
                 "flight update changed lazy vx while synchronizing");
    expect_close(lazy->vy, expected_vy,
                 "flight update changed lazy vy while synchronizing");
#if THREE_D
    expect_close(lazy->z, expected_z,
                 "flight update advanced lazy z with the new law");
    expect_close(lazy->vz, expected_vz,
                 "flight update changed lazy vz while synchronizing");
#endif
    expect_true(particles == allocation_before && ncol == collisions_before &&
                    ncross == crossings_before &&
                    memcmp(particles + 1, &other_before,
                           sizeof(other_before)) == 0,
                "lazy-particle flight update reset unrelated runtime state");
}

static void test_flight_and_thermostat_rescheduling(void)
{
    char message[512];

    SimulationConfig current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    SimulationConfig next = current;
    const PreservedState flight_before = capture_state();
    const double old_crossing_time = eventList[0]->t;
    configure_changed_flight(&next);
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "flight-law update", message);
#if (STRESS || HEAT)
    expect_true(addShear, "Lees-Edwards flight global was not activated");
#else
    expect_true(addField, "field flight global was not activated");
    expect_true(eventList[0]->t != old_crossing_time,
                "flight-law update did not replace particle predictions");
#endif
    (void)old_crossing_time;
    expect_particle_state_preserved(&flight_before, "flight-law update");

    current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    next = current;
    const PreservedState thermostat_before = capture_state();
    const double particle_crossing_time = eventList[0]->t;
    const double particle_collision_time = eventList[N]->t;
    next.driving.thermostat = SIM_THERMOSTAT_LANGEVIN;
    next.driving.thermostat_interval = 0.25;
    reconfigure(&current, &next, SIM_RECONFIGURE_EVENTS_REBUILT,
                "thermostat update", message);
    expect_true(noise == SIM_THERMOSTAT_LANGEVIN,
                "thermostat global was not activated");
    expect_true(eventList[2*N + 2]->type == ADDINGNOISE,
                "thermostat auxiliary event was not scheduled");
    expect_close(eventList[2*N + 2]->t, t + dtnoise,
                 "thermostat auxiliary deadline is wrong");
    expect_true(eventList[0]->t == particle_crossing_time &&
                    eventList[N]->t == particle_collision_time,
                "thermostat update unnecessarily rescheduled particle events");
    expect_particle_state_preserved(&thermostat_before, "thermostat update");
}

static void test_restart_required_and_rendering_only(void)
{
    char message[512];
    SimulationConfig current = prepare_state(SIM_BOUNDARY_PERIODIC, false);
    SimulationConfig next = current;
    const PreservedState before = capture_state();

    next.initialization.random_seed++;
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "initializer update", message);
    expect_particle_state_preserved(&before, "initializer restart request");

    next = current;
    next.particle_shape.vertices[0][0] -= 0.01;
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "particle-shape template update", message);
    expect_particle_state_preserved(&before,
                                    "particle-shape restart request");

#if !THREE_D
    /* Use a binary random-growth setup so both initializer identities are
       valid while every field other than the mode remains identical.  The
       membrane initializer itself is intentionally unavailable in 3D. */
    current.initialization.composition = SIM_COMPOSITION_BINARY;
    current.initialization.small_fraction = 0.5;
    simulation_config_apply_globals(&current);
    next = current;
    next.initialization.mode = SIM_INIT_MEMBRANE_OSMOSIS;
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "random-to-membrane initializer update", message);
    expect_particle_state_preserved(
        &before, "random-to-membrane initializer restart request");

    current.initialization.mode = SIM_INIT_MEMBRANE_OSMOSIS;
    simulation_config_apply_globals(&current);
    next = current;
    next.initialization.mode = SIM_INIT_RANDOM_GROWTH;
    const PreservedState reverse_before = capture_state();
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "membrane-to-random initializer update", message);
    expect_particle_state_preserved(
        &reverse_before, "membrane-to-random initializer restart request");

    current = next;
    current.initialization.composition = SIM_COMPOSITION_MONODISPERSE;
    current.initialization.small_fraction = 0.0;
    simulation_config_apply_globals(&current);
#endif

    next = current;
    next.initialization.aspect_ratio = 1.25;
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "box update", message);
    expect_particle_state_preserved(&before, "box restart request");

    next = current;
    next.initialization.composition = SIM_COMPOSITION_BINARY;
    next.initialization.small_fraction = 0.5;
    reconfigure(&current, &next, SIM_RECONFIGURE_RESTART_REQUIRED,
                "composition update", message);
    expect_particle_state_preserved(&before, "composition restart request");

    next = current;
    next.rendering.target_fps++;
    next.rendering.velocity_arrows = !next.rendering.velocity_arrows;
    node *crossing = eventList[0];
    node *collision = eventList[N];
    reconfigure(&current, &next, SIM_RECONFIGURE_UNCHANGED,
                "rendering-only update", message);
    expect_true(eventList[0] == crossing && eventList[N] == collision,
                "rendering-only update performed physics work");
    expect_particle_state_preserved(&before, "rendering-only update");
}

int main(void)
{
    test_response_only_restitution();
    test_collision_injection_keeps_configuration_declarative();
    test_wall_topology();
    test_spherical_wall_fallback();
    test_center_membrane_topology();
    test_potential_tiers();
    test_lazy_particle_flight_synchronization();
    test_flight_and_thermostat_rescheduling();
    test_restart_required_and_rendering_only();
    if (particles != NULL || cellList != NULL || eventList != NULL)
        freeArrays();

    if (failures != 0) {
        fprintf(stderr, "%d incremental reconfiguration regression(s) failed\n",
                failures);
        return 1;
    }
    printf("%dD incremental reconfiguration regressions passed\n",
           THREE_D ? 3 : 2);
    return 0;
}
