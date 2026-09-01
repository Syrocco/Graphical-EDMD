#include "EDMD.h"
#include "convex_polygon.h"
#include "polygon_particles.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

extern int N, Nbig, Nsmall;
extern int periodicDumps, ther;
extern int addCircularWall;
extern int polygonCollisionCellNeighborRange;
extern double Lx, Ly, t, vr, Einit;
extern double cumulativeWallImpulse;
extern particle *particles;
extern node *nextEvent;
extern unsigned long int ncol, ncross;
extern void init_genrand(unsigned long seed);

static int validation_has_error(const SimulationConfig *config,
                                const BuildCapabilities *capabilities,
                                const char *field)
{
    ConfigIssue issues[64];
    const size_t count = simulation_config_validate(
        config, capabilities, issues, sizeof(issues)/sizeof(issues[0]));
    const size_t stored = count < sizeof(issues)/sizeof(issues[0])
        ? count : sizeof(issues)/sizeof(issues[0]);

    for (size_t index = 0; index < stored; ++index) {
        if (issues[index].severity == CONFIG_ISSUE_ERROR &&
            strcmp(issues[index].field, field) == 0) {
            return 1;
        }
    }
    return 0;
}

static ConvexPolygon polygon_from_config(const SimParticleShapeConfig *shape)
{
    ConvexPolygon polygon;
    memset(&polygon, 0, sizeof(polygon));
    polygon.count = shape->vertex_count;
    for (int index = 0; index < polygon.count; ++index) {
        polygon.vertices[index].x = shape->vertices[index][0];
        polygon.vertices[index].y = shape->vertices[index][1];
    }
    return polygon;
}

static void test_config_uses_runtime_geometry_tolerances(void)
{
    SimulationConfig config;
    BuildCapabilities capabilities = simulation_build_capabilities();
    ConvexPolygon polygon;

    simulation_config_set_defaults(&config, &capabilities);
    capabilities.dimensions = 2;
    capabilities.graphics = true;
    capabilities.convex_polygon_particles = true;
    capabilities.stress = false;
    capabilities.heat = false;
    capabilities.rough_particles = false;
    capabilities.intruder = false;
    capabilities.additive_contacts = true;
    /* The debug verifier now evaluates exact SAT and polygon wall vertices,
       so enabling it must not make polygon mode an invalid configuration. */
    capabilities.overlap_verification = true;
    config.particle_shape.mode = SIM_PARTICLE_CONVEX_POLYGONS;
    config.dynamics.injection = SIM_INJECTION_NONE;
    config.observables.pressure = false;
    config.output.periodic_dumps = false;

    /* This centrally symmetric unit-circumradius polygon is canonical and
       strictly convex, but its two shortest edges are below the geometry
       module's scale-relative runtime threshold. */
    {
        const double angle = 1.0e-7;
        const double cosine = cos(angle);
        const double sine = sin(angle);
        const double vertices[6][2] = {
            {1.0, 0.0}, {cosine, sine}, {0.0, 1.0},
            {-1.0, 0.0}, {-cosine, -sine}, {0.0, -1.0},
        };
        config.particle_shape.vertex_count = 6;
        memset(config.particle_shape.vertices, 0,
               sizeof(config.particle_shape.vertices));
        memcpy(config.particle_shape.vertices, vertices, sizeof(vertices));
    }
    polygon = polygon_from_config(&config.particle_shape);
    assert(convex_polygon_validate(&polygon) ==
           CONVEX_POLYGON_EDGE_TOO_SHORT);
    assert(!validation_has_error(&config, &capabilities,
                                 "particle_shape.mode"));
    assert(validation_has_error(&config, &capabilities,
                                "particle_shape.vertices"));

    /* No edge is short here, but the two nearly flat vertices have a turn
       sine below the runtime threshold. Central symmetry keeps the area
       centroid exactly at the origin and the corner radius is one. */
    {
        const double scale = sqrt(0.5);
        const double bulge = 1.0e-8;
        const double vertices[6][2] = {
            {-scale, -scale}, {0.0, -scale*(1.0 + bulge)},
            {scale, -scale}, {scale, scale},
            {0.0, scale*(1.0 + bulge)}, {-scale, scale},
        };
        memset(config.particle_shape.vertices, 0,
               sizeof(config.particle_shape.vertices));
        memcpy(config.particle_shape.vertices, vertices, sizeof(vertices));
    }
    polygon = polygon_from_config(&config.particle_shape);
    assert(convex_polygon_validate(&polygon) ==
           CONVEX_POLYGON_NOT_STRICTLY_CCW_CONVEX);
    assert(validation_has_error(&config, &capabilities,
                                "particle_shape.vertices"));
}

static void configure_square_particles(void)
{
    SimParticleShapeConfig shape;
    memset(&shape, 0, sizeof(shape));
    shape.mode = SIM_PARTICLE_CONVEX_POLYGONS;
    shape.vertex_count = 4;
    shape.vertices[0][0] = -1.0; shape.vertices[0][1] = -1.0;
    shape.vertices[1][0] =  1.0; shape.vertices[1][1] = -1.0;
    shape.vertices[2][0] =  1.0; shape.vertices[2][1] =  1.0;
    shape.vertices[3][0] = -1.0; shape.vertices[3][1] =  1.0;
    shape.random_orientation = false;
    shape.angular_temperature = 0.0;
    polygon_particles_set_config(&shape);
}

static void configure_spinning_square_particles(void)
{
    SimParticleShapeConfig shape;
    memset(&shape, 0, sizeof(shape));
    shape.mode = SIM_PARTICLE_CONVEX_POLYGONS;
    shape.vertex_count = 4;
    shape.vertices[0][0] = -1.0; shape.vertices[0][1] = -1.0;
    shape.vertices[1][0] =  1.0; shape.vertices[1][1] = -1.0;
    shape.vertices[2][0] =  1.0; shape.vertices[2][1] =  1.0;
    shape.vertices[3][0] = -1.0; shape.vertices[3][1] =  1.0;
    shape.random_orientation = false;
    shape.angular_temperature = 2.0;
    polygon_particles_set_config(&shape);
}

static void initialize_pair(void)
{
    N = 2;
    Nbig = 2;
    Nsmall = 0;
    Lx = 20.0;
    Ly = 20.0;
    t = 0.0;
    vr = 1.0;
    ncol = 0;
    ncross = 0;
    periodicDumps = 0;
    ther = 0;
    simulationClearFault();

    particles = calloc((size_t)N, sizeof(*particles));
    assert(particles != NULL);
    for (int index = 0; index < N; ++index) {
        particle *value = particles + index;
        value->num = index;
        value->type = 1;
        value->m = 1.0;
        value->rad = 1.0;
        value->t = 0.0;
        value->lastColl = 0.0;
    }
    particles[0].x = 6.0;
    particles[0].y = 10.0;
    particles[0].vx = 1.0;
    particles[1].x = 13.5;
    particles[1].y = 10.0;
    particles[1].vx = -1.0;

    adaptCellGridToParticles();
    /* Diagnostics must use the particles that actually exist, including an
       arbitrary loaded/binary radius, rather than Nbig/Nsmall or target phi.
       The compiled fixture has a circular container of radius Lx/2 and the
       canonical square has unit-radius area two. */
    particles[1].rad = 0.5;
    assert(fabs(simulationCurrentPackingFraction(t) -
                2.0*(1.0 + 0.25)/(M_PI*100.0)) < 1.0e-14);
    particles[1].rad = 1.0;
    assert(fabs(simulationCurrentPackingFraction(t) -
                4.0/(M_PI*100.0)) < 1.0e-14);
    for (int index = 0; index < N; ++index) {
        particles[index].cell[0] = coordToCell(particles[index].x, 1);
        particles[index].cell[1] = coordToCell(particles[index].y, 0);
    }
    cellListInit();
    eventListInit();
    assert(polygon_particles_active());
    assert(!simulationHasFault());

    /* Polygon-aware preparation keeps its own handler while the radius is
       growing. The public packing metric must still evaluate each lazy
       circumscale at the requested time. */
    particles[0].rad = particles[1].rad = 0.5;
    particles[0].vr = particles[1].vr = 0.1;
    polygon_particles_begin_growth();
    assert(polygon_particles_growing());
    assert(fabs(simulationCurrentPackingFraction(1.0) -
                2.0*(0.6*0.6 + 0.6*0.6)/(M_PI*100.0)) < 1.0e-14);
    particles[0].rad = particles[1].rad = 1.0;
    polygon_particles_finish_growth();
    assert(!polygon_particles_growing());
    assert(fabs(simulationCurrentPackingFraction(t) -
                2.0*(0.6*0.6 + 0.6*0.6)/(M_PI*100.0)) < 1.0e-14);
    particles[0].rad = particles[1].rad = 1.0;
}

static void test_local_overlap_verifier(void)
{
    char message[256];
    assert(polygon_particles_verify_local_no_overlaps(
        t, 1.0e-10, message, sizeof(message)));
    assert(message[0] == '\0');
    assert(!polygon_particles_verify_local_no_overlaps(
        t + 0.25, 1.0e-10, message, sizeof(message)));
    assert(strstr(message, "current simulation time") != NULL);

    const double saved_x = particles[1].x;
    const double saved_y = particles[1].y;
    const double saved_particle_time = particles[1].t;
    removeFromCell(1);
    particles[1].x = particles[0].x + 0.25;
    particles[1].y = particles[0].y;
    particles[1].t = t;
    addToCell(1);
    assert(!polygon_particles_verify_local_no_overlaps(
        t, 1.0e-10, message, sizeof(message)));
    assert(strstr(message, "particles 0 and 1") != NULL);

    removeFromCell(1);
    particles[1].x = saved_x;
    particles[1].y = saved_y;
    particles[1].t = saved_particle_time;
    addToCell(1);
    assert(polygon_particles_verify_local_no_overlaps(
        t, 1.0e-10, message, sizeof(message)));
}

static void test_periodic_neighbor_enumeration(void)
{
    const int cell_counts[] = {1, 2, 5, 24, 25, 26};
    const int saved_range = polygonCollisionCellNeighborRange;
    polygonCollisionCellNeighborRange =
        POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE;
    for (size_t fixture = 0;
         fixture < sizeof(cell_counts)/sizeof(cell_counts[0]); ++fixture) {
        const int count = cell_counts[fixture];
        for (int center = 0; center < count; ++center) {
            int coordinates[POLYGON_PARTICLES_MAX_CELL_COORDINATES];
            bool expected[26] = {false};
            bool observed[26] = {false};
            for (int offset = -POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE;
                 offset <= POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE;
                 ++offset) {
                int coordinate = (center + offset) % count;
                if (coordinate < 0) coordinate += count;
                expected[coordinate] = true;
            }
            const int found = polygon_particles_test_neighboring_coordinates(
                center, count, coordinates);
            const int expected_count = count <
                    POLYGON_PARTICLES_MAX_CELL_COORDINATES
                ? count : POLYGON_PARTICLES_MAX_CELL_COORDINATES;
            assert(found == expected_count);
            for (int index = 0; index < found; ++index) {
                assert(coordinates[index] >= 0 && coordinates[index] < count);
                assert(!observed[coordinates[index]]);
                observed[coordinates[index]] = true;
            }
            for (int coordinate = 0; coordinate < count; ++coordinate)
                assert(observed[coordinate] == expected[coordinate]);
        }
    }
    polygonCollisionCellNeighborRange = saved_range;
}

static void process_until_first_collision(void)
{
    for (int event_count = 0; event_count < 64 && ncol == 0; ++event_count) {
        nextEvent = findNextEvent();
        assert(nextEvent != NULL);
        assert(isfinite(nextEvent->t));
        assert(nextEvent->t >= t);
        t = nextEvent->t;
        removeEventFromQueue(nextEvent);
        switch (nextEvent->type) {
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
                assert(!"unexpected auxiliary event in polygon pair test");
        }
        assert(!simulationHasFault());
    }
    assert(ncol == 1);
}

static void process_next_polygon_event(void)
{
    const double previous_time = t;
    nextEvent = findNextEvent();
    assert(nextEvent != NULL);
    assert(isfinite(nextEvent->t));
    assert(nextEvent->t >= previous_time);
    assert(nextEvent->t < never);
    t = nextEvent->t;
    removeEventFromQueue(nextEvent);
    switch (nextEvent->type) {
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
            assert(!"unexpected event in polygon stress test");
    }
    assert(!simulationHasFault());
}

static void initialize_single_for_circular_wall(void)
{
    N = 1;
    Nbig = 1;
    Nsmall = 0;
    Lx = 20.0;
    Ly = 20.0;
    t = 0.0;
    vr = 1.0;
    ncol = 0;
    ncross = 0;
    periodicDumps = 0;
    ther = 0;
    simulationClearFault();

    particles = calloc((size_t)N, sizeof(*particles));
    assert(particles != NULL);
    particles[0].num = 0;
    particles[0].type = 1;
    particles[0].m = 1.0;
    particles[0].rad = 1.0;
    particles[0].x = 10.0;
    particles[0].y = 10.0;
    particles[0].vx = 1.0;

    adaptCellGridToParticles();
    particles[0].cell[0] = coordToCell(particles[0].x, 1);
    particles[0].cell[1] = coordToCell(particles[0].y, 0);
    cellListInit();
    eventListInit();
    assert(polygon_particles_active());
    assert(!simulationHasFault());
}

static void initialize_centered_zero_clearance_square(double vx)
{
    N = 1;
    Nbig = 1;
    Nsmall = 0;
    Lx = 2.0;
    Ly = 2.0;
    t = 0.0;
    vr = 1.0;
    ncol = 0;
    ncross = 0;
    periodicDumps = 0;
    ther = 0;
    simulationClearFault();

    particles = calloc((size_t)N, sizeof(*particles));
    assert(particles != NULL);
    particles[0].num = 0;
    particles[0].type = 1;
    particles[0].m = 1.0;
    particles[0].rad = 1.0;
    particles[0].x = 1.0;
    particles[0].y = 1.0;
    particles[0].vx = vx;

    adaptCellGridToParticles();
    particles[0].cell[0] = coordToCell(particles[0].x, 1);
    particles[0].cell[1] = coordToCell(particles[0].y, 0);
    cellListInit();
    eventListInit();
    assert(polygon_particles_active());
}

static void initialize_spinning_face_at_circular_wall(void)
{
    const double half_vertex = M_SQRT1_2;
    const double wall_radius = 10.0;
    const double radial_x = sqrt(wall_radius*wall_radius -
                                 half_vertex*half_vertex);
    N = 1;
    Nbig = 1;
    Nsmall = 0;
    Lx = 20.0;
    Ly = 20.0;
    t = 0.0;
    vr = 1.0;
    ncol = 0;
    ncross = 0;
    periodicDumps = 0;
    ther = 0;
    simulationClearFault();
    init_genrand(8128UL);

    particles = calloc((size_t)N, sizeof(*particles));
    assert(particles != NULL);
    particles[0].num = 0;
    particles[0].type = 1;
    particles[0].m = 1.0;
    particles[0].rad = 1.0;
    /* The two right-face vertices are exactly tied on the circle.  Their
       normal speeds are opposite under pure spin, so an averaged midpoint
       incorrectly appears stationary. */
    particles[0].x = 10.0 + radial_x - half_vertex;
    particles[0].y = 10.0;

    adaptCellGridToParticles();
    particles[0].cell[0] = coordToCell(particles[0].x, 1);
    particles[0].cell[1] = coordToCell(particles[0].y, 0);
    cellListInit();
    eventListInit();
    assert(polygon_particles_active());
    assert(fabs(polygon_particles_omega(0)) > 1.0e-6);
    assert(!simulationHasFault());
}

static double spinning_face_maximum_outward_speed(void)
{
    const double half_vertex = M_SQRT1_2;
    const double wall_radius = 10.0;
    const double radial_x = sqrt(wall_radius*wall_radius -
                                 half_vertex*half_vertex);
    const double omega = polygon_particles_omega(0);
    double maximum_speed = -INFINITY;
    for (int sign = -1; sign <= 1; sign += 2) {
        const double arm_x = half_vertex;
        const double arm_y = sign*half_vertex;
        const double normal_x = radial_x/wall_radius;
        const double normal_y = arm_y/wall_radius;
        const double point_vx = particles[0].vx - omega*arm_y;
        const double point_vy = particles[0].vy + omega*arm_x;
        const double speed = point_vx*normal_x + point_vy*normal_y;
        if (speed > maximum_speed) maximum_speed = speed;
    }
    return maximum_speed;
}

static double polygon_system_energy(void)
{
    double energy = polygon_particles_rotational_energy();
    for (int index = 0; index < N; ++index) {
        energy += 0.5*particles[index].m*
            (particles[index].vx*particles[index].vx +
             particles[index].vy*particles[index].vy);
    }
    return energy;
}

static double polygon_translational_energy(void)
{
    double energy = 0.0;
    for (int index = 0; index < N; ++index) {
        energy += 0.5*particles[index].m*
            (particles[index].vx*particles[index].vx +
             particles[index].vy*particles[index].vy);
    }
    return energy;
}

static void test_live_mass_edit_updates_polygon_inertia(void)
{
    const double original_mass = particles[0].m;
    const double original_omega = polygon_particles_omega(0);
    const double edited_mass = 2.75*original_mass;

    assert(original_mass > 0.0);
    /* Normalize the deliberately hand-edited radius used by initialize_pair:
       the public transaction must derive inertia from the live radius rather
       than merely multiplying a potentially stale sidecar value. */
    assert(polygon_particles_set_mass(0, original_mass));
    const double original_inertia = polygon_particles_inertia(0);
    assert(original_inertia > 0.0);
    assert(polygon_particles_omega(0) == original_omega);
    assert(polygon_particles_set_mass(0, edited_mass));
    assert(particles[0].m == edited_mass);
    assert(fabs(polygon_particles_inertia(0) -
                2.75*original_inertia) <
           1.0e-13*(1.0 + 2.75*original_inertia));
    /* A mass edit changes angular momentum and rotational energy, not the
       instantaneous angular velocity selected by the running state. */
    assert(polygon_particles_omega(0) == original_omega);

    const double edited_inertia = polygon_particles_inertia(0);
    assert(!polygon_particles_set_mass(0, 0.0));
    assert(particles[0].m == edited_mass);
    assert(polygon_particles_inertia(0) == edited_inertia);
    assert(polygon_particles_omega(0) == original_omega);

    assert(polygon_particles_set_mass(0, original_mass));
    assert(fabs(polygon_particles_inertia(0) - original_inertia) <
           1.0e-13*(1.0 + original_inertia));
    assert(polygon_particles_omega(0) == original_omega);
}

static void test_growth_energy_guard(void)
{
    particles[0].rad = particles[1].rad = 1.0;
    particles[0].vr = particles[1].vr = 0.1;
    polygon_particles_begin_growth();
    assert(polygon_particles_growing());
    particles[0].vx = 1.0e9;
    particles[0].vy = -1.0e9;
    /* A live mass transaction refreshes the aggregate preparation energy,
       including velocity fields changed earlier in the same transaction. */
    assert(!polygon_particles_set_mass(0, particles[0].m));
    assert(simulationHasFault());
    assert(strstr(simulationFaultMessage(), "aggregate kinetic-energy safety ceiling") != NULL);
    simulationClearFault();
    particles[0].vx = 0.0;
    particles[0].vy = 0.0;
    polygon_particles_finish_growth();
    assert(!polygon_particles_growing());
}

static void test_graphical_recovery_primitives(void)
{
    const size_t checkpoint_bytes = polygon_particles_checkpoint_bytes();
    assert(checkpoint_bytes > 0);
    void *checkpoint = malloc(checkpoint_bytes);
    assert(checkpoint != NULL);

    const double saved_time = t;
    const double saved_angle = polygon_particles_angle_at(0, saved_time);
    const double saved_omega = polygon_particles_omega(0);
    assert(isfinite(saved_angle));
    assert(isfinite(saved_omega));
    assert(polygon_particles_capture_checkpoint(checkpoint, checkpoint_bytes));
    assert(!polygon_particles_capture_checkpoint(NULL, checkpoint_bytes));
    assert(!polygon_particles_capture_checkpoint(
        checkpoint, checkpoint_bytes - 1));

    /* A thermostat kick changes both the synchronized pose record and spin.
       Restoring the value-only sidecar must recover both exactly even though
       the global clock remains at the later replay time. */
    t = saved_time + 0.125;
    assert(polygon_particles_apply_angular_langevin(
        0, 0.0, 1.0, 0.5, false));
    assert(polygon_particles_omega(0) != saved_omega);
    assert(polygon_particles_restore_checkpoint(checkpoint, checkpoint_bytes));
    assert(polygon_particles_omega(0) == saved_omega);
    assert(fabs(polygon_particles_angle_at(0, t) -
                (saved_angle + saved_omega*(t - saved_time))) < 1.0e-13);

    /* Validation precedes the memcpy, so a torn/corrupt ring entry cannot
       partially damage the live sidecar. */
    ((double *)checkpoint)[0] = NAN;
    const double angle_before_rejected_restore =
        polygon_particles_angle_at(0, t);
    assert(!polygon_particles_restore_checkpoint(
        checkpoint, checkpoint_bytes));
    assert(polygon_particles_omega(0) == saved_omega);
    assert(polygon_particles_angle_at(0, t) ==
           angle_before_rejected_restore);
    free(checkpoint);

    /* Recovery precision is opt-in and expires on simulation time, so an
       ordinary (including headless) run never remains on the expensive CCD
       path merely because wall-clock rendering paused. */
    assert(!polygon_particles_test_recovery_precision_active());
    polygon_particles_enable_recovery_precision(0.25);
    assert(polygon_particles_test_recovery_precision_active());
    t = nextafter(t + 0.25, INFINITY);
    assert(!polygon_particles_test_recovery_precision_active());
    polygon_particles_enable_recovery_precision(-1.0);
    assert(!polygon_particles_test_recovery_precision_active());

    simulationClearFault();
    simulationReportFault("ordinary fault");
    assert(simulationHasFault());
    assert(!simulationFaultIsRecoverableOverlap());
    simulationClearFault();
    simulationReportRecoverableOverlapFault("small polygon overlap");
    assert(simulationHasFault());
    assert(simulationFaultIsRecoverableOverlap());
    simulationClearFault();
}

static void test_angular_langevin_equilibrium(void)
{
    const double bath_temperature = 2.75;
    const double gamma = 1.3;
    const double independent_interval = 20.0/gamma;
    const int samples = 50000;
    const double inertia = polygon_particles_inertia(0);
    assert(inertia > 0.0);
    init_genrand(7331UL);

    double squared_sum = 0.0;
    double mean_sum = 0.0;
    for (int sample = 0; sample < samples; ++sample) {
        assert(polygon_particles_apply_angular_langevin(
            0, bath_temperature, gamma, independent_interval, false));
        const double omega = polygon_particles_omega(0);
        mean_sum += omega;
        squared_sum += inertia*omega*omega;
    }
    const double measured_temperature = squared_sum/samples;
    assert(fabs(measured_temperature-bath_temperature) <
           0.035*bath_temperature);
    assert(fabs(mean_sum/samples) <
           0.035*sqrt(bath_temperature/inertia));

    /* The analytical update is continuous in orientation even when invoked
       after a lazily represented free flight. */
    t += 0.375;
    const double angle_before = polygon_particles_angle_at(0, t);
    const double omega_before_cooling = polygon_particles_omega(0);
    assert(polygon_particles_apply_angular_langevin(
        0, 0.0, gamma, 0.25, false));
    assert(fabs(polygon_particles_angle_at(0, t)-angle_before) < 1.0e-13);
    assert(fabs(polygon_particles_omega(0) -
                omega_before_cooling*exp(-gamma*0.25)) < 1.0e-13);

    const double omega_before_rejected_euler = polygon_particles_omega(0);
    assert(!polygon_particles_apply_angular_langevin(
        0, bath_temperature, gamma, 0.25, true));
    assert(polygon_particles_omega(0) == omega_before_rejected_euler);
}

static void test_two_contact_lcp(void)
{
    double impulse[2] = {0, 0};
    const double pair_matrix[2][2] = {{5.0, -1.0}, {-1.0, 5.0}};
    const double mixed_pair_bound[2] = {2.0, -0.1};
    assert(polygon_particles_test_solve_contact_lcp(
        pair_matrix, mixed_pair_bound, 2, 1.0e-12, impulse));
    assert(fabs(impulse[0] - 0.4125) < 1.0e-12);
    assert(fabs(impulse[1] - 0.0625) < 1.0e-12);

    const double symmetric_bound[2] = {2.0, 2.0};
    assert(polygon_particles_test_solve_contact_lcp(
        pair_matrix, symmetric_bound, 2, 1.0e-12, impulse));
    assert(fabs(impulse[0] - 0.5) < 1.0e-12);
    assert(fabs(impulse[1] - 0.5) < 1.0e-12);

    const double wall_matrix[2][2] = {{2.5, -0.5}, {-0.5, 2.5}};
    const double mixed_wall_bound[2] = {2.0, -0.1};
    assert(polygon_particles_test_solve_contact_lcp(
        wall_matrix, mixed_wall_bound, 2, 1.0e-12, impulse));
    assert(fabs(impulse[0] - 0.825) < 1.0e-12);
    assert(fabs(impulse[1] - 0.125) < 1.0e-12);

    /* Recorded dense-triangle event: one face endpoint was separating while
       the other was closing.  The raw LCP was nonpenetrating but injected
       2.26e-2 kinetic-energy units.  Elastic projection must retain a
       repulsive impulse, leave both endpoints nonclosing, and conserve
       energy to roundoff. */
    const double recorded_matrix[2][2] = {
        {6.9347787326198072, -2.8184370167088848},
        {-2.8184370167088848, 6.9347767594809575}
    };
    const double recorded_velocity[2] = {
        0.11340094959327045, -2.981961793236839
    };
    const double recorded_bound[2] = {
        -0.11340094959327045, 5.9639235864736779
    };
    assert(polygon_particles_test_solve_contact_lcp(
        recorded_matrix, recorded_bound, 2, 1.0e-12, impulse));
    double raw_response[2] = {0, 0};
    double raw_energy_change = 0;
    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column)
            raw_response[row] +=
                recorded_matrix[row][column]*impulse[column];
        raw_energy_change += impulse[row]*recorded_velocity[row] +
            0.5*impulse[row]*raw_response[row];
    }
    assert(raw_energy_change > 2.2e-2);
    assert(polygon_particles_test_project_elastic_contact_energy(
        recorded_matrix, recorded_velocity, 2, 1.0e-12, impulse));
    double projected_response[2] = {0, 0};
    double projected_energy_change = 0;
    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column)
            projected_response[row] +=
                recorded_matrix[row][column]*impulse[column];
        projected_energy_change += impulse[row]*recorded_velocity[row] +
            0.5*impulse[row]*projected_response[row];
        assert(recorded_velocity[row] + projected_response[row] >=
               -1.0e-12);
        assert(impulse[row] > 0);
    }
    assert(fabs(projected_energy_change) < 1.0e-13);

    /* Strong negative coupling can make the simultaneous LCP direction's
       elastic root leave the originally closing endpoint still closing.
       The energy policy must then use bounded propagative one-point
       reflections instead of accepting the injecting LCP. */
    const double fallback_matrix[2][2] = {
        {1.0, -0.9}, {-0.9, 1.0}
    };
    const double fallback_velocity[2] = {1.0, -1.0};
    const double fallback_bound[2] = {-1.0, 2.0};
    assert(polygon_particles_test_solve_contact_lcp(
        fallback_matrix, fallback_bound, 2, 1.0e-12, impulse));
    double unprojectable_impulse[2] = {impulse[0], impulse[1]};
    assert(!polygon_particles_test_project_elastic_contact_energy(
        fallback_matrix, fallback_velocity, 2, 1.0e-12,
        unprojectable_impulse));
    assert(polygon_particles_test_enforce_contact_energy_policy(
        fallback_matrix, fallback_velocity, 2, 1.0e-12, true, impulse));
    double fallback_response[2] = {0, 0};
    double fallback_energy_change = 0;
    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column)
            fallback_response[row] +=
                fallback_matrix[row][column]*impulse[column];
        fallback_energy_change += impulse[row]*fallback_velocity[row] +
            0.5*impulse[row]*fallback_response[row];
        assert(fallback_velocity[row] + fallback_response[row] >=
               -1.0e-12);
        assert(impulse[row] >= 0);
    }
    assert(fabs(fallback_energy_change) < 1.0e-12);

    /* The same guard caps a nominally dissipative mixed-manifold solution
       when the LCP would nevertheless add energy. */
    assert(polygon_particles_test_solve_contact_lcp(
        recorded_matrix, recorded_bound, 2, 1.0e-12, impulse));
    assert(polygon_particles_test_enforce_contact_energy_policy(
        recorded_matrix, recorded_velocity, 2, 1.0e-12, false, impulse));
    double capped_energy_change = 0;
    for (int row = 0; row < 2; ++row) {
        double response = 0;
        for (int column = 0; column < 2; ++column)
            response += recorded_matrix[row][column]*impulse[column];
        capped_energy_change += impulse[row]*recorded_velocity[row] +
            0.5*impulse[row]*response;
    }
    assert(capped_energy_change < 1.0e-13);
}

static void test_absolute_clock_deadline_guard(void)
{
    const double saved_time = t;
    double deadline = 0.0;
    t = 1.0e20;
    simulationClearFault();
    assert(!polygon_particles_test_event_deadline_after(1.0, &deadline));
    assert(simulationHasFault());
    assert(strstr(simulationFaultMessage(), "resolution limit") != NULL);

    simulationClearFault();
    const double representable_delay = nextafter(t, INFINITY) - t;
    assert(representable_delay > 0.0);
    assert(polygon_particles_test_event_deadline_after(
        representable_delay, &deadline));
    assert(deadline > t);
    assert(!simulationHasFault());
    t = saved_time;
}

static void initialize_many_for_calendar_stress(void)
{
    N = 12;
    Nbig = N;
    Nsmall = 0;
    Lx = 20.0;
    Ly = 20.0;
    t = 0.0;
    vr = 1.0;
    ncol = 0;
    ncross = 0;
    periodicDumps = 0;
    ther = 0;
    simulationClearFault();

    particles = calloc((size_t)N, sizeof(*particles));
    assert(particles != NULL);
    double mean_vx = 0.0;
    double mean_vy = 0.0;
    for (int index = 0; index < N; ++index) {
        particle *value = particles + index;
        value->num = index;
        value->type = 1;
        value->m = 1.0;
        value->rad = 0.8;
        value->x = 5.5 + 3.0*(index % 4);
        value->y = 7.0 + 3.0*(index / 4);
        value->vx = 0.8*cos(0.73*(index + 1));
        value->vy = 0.8*sin(1.13*(index + 1));
        mean_vx += value->vx;
        mean_vy += value->vy;
    }
    mean_vx /= N;
    mean_vy /= N;
    for (int index = 0; index < N; ++index) {
        particles[index].vx -= mean_vx;
        particles[index].vy -= mean_vy;
    }

    adaptCellGridToParticles();
    for (int index = 0; index < N; ++index) {
        particles[index].cell[0] = coordToCell(particles[index].x, 1);
        particles[index].cell[1] = coordToCell(particles[index].y, 0);
    }
    cellListInit();
    eventListInit();
    assert(polygon_particles_active());
    assert(!simulationHasFault());
}

int main(void)
{
    addCircularWall = 1;
    test_periodic_neighbor_enumeration();
    test_config_uses_runtime_geometry_tolerances();
    test_two_contact_lcp();
    test_absolute_clock_deadline_guard();

    configure_spinning_square_particles();
    init_genrand(611953UL);
    initialize_pair();
    test_live_mass_edit_updates_polygon_inertia();
    test_graphical_recovery_primitives();
    test_growth_energy_guard();
    freeArrays();

    configure_square_particles();
    assert(fabs(polygon_particles_unit_area() - 2.0) < 1.0e-14);
    initialize_pair();
    test_local_overlap_verifier();
    {
        char overlap_message[256];
        assert(polygon_particles_verify_no_overlaps(
            t, 1.0e-10, overlap_message, sizeof(overlap_message)));
        assert(overlap_message[0] == '\0');
        /* Verify lazy poses directly: no synchronization or calendar event is
           processed, but the future ballistic poses overlap after contact. */
        assert(!polygon_particles_verify_no_overlaps(
            3.2, 1.0e-10, overlap_message, sizeof(overlap_message)));
        assert(strstr(overlap_message, "particles 0 and 1") != NULL);
        assert(t == 0.0);
        assert(particles[0].t == 0.0 && particles[1].t == 0.0);
    }
    const double initial_energy =
        0.5*particles[0].m*particles[0].vx*particles[0].vx +
        0.5*particles[1].m*particles[1].vx*particles[1].vx;

    process_until_first_collision();
    assert(fabs(t - 0.5*(7.5 - sqrt(2.0))) < 1.0e-7);
    assert(fabs(particles[0].vx + 1.0) < 1.0e-8);
    assert(fabs(particles[1].vx - 1.0) < 1.0e-8);
    assert(fabs(polygon_particles_omega(0)) < 1.0e-9);
    assert(fabs(polygon_particles_omega(1)) < 1.0e-9);
    const double final_energy =
        0.5*particles[0].m*particles[0].vx*particles[0].vx +
        0.5*particles[1].m*particles[1].vx*particles[1].vx +
        polygon_particles_rotational_energy();
    assert(fabs(final_energy - initial_energy) < 1.0e-9);

    test_angular_langevin_equilibrium();

    freeArrays();

    /* Linear-energy normalization must not count the independently sampled
       polygon spin energy in a factor that is applied only to vx/vy. */
    configure_spinning_square_particles();
    initialize_pair();
    Einit = 1.0;
    assert(polygon_particles_rotational_energy() > 0.0);
    normalizePhysicalQ();
    assert(fabs(polygon_translational_energy() - N*Einit) < 1.0e-10);
    freeArrays();

    configure_square_particles();
    initialize_single_for_circular_wall();
    const double wall_impulse_before = cumulativeWallImpulse;
    process_until_first_collision();
    assert(fabs(t - (sqrt(100.0 - 0.5) - sqrt(0.5))) < 1.0e-7);
    assert(fabs(particles[0].vx + 1.0) < 1.0e-8);
    assert(fabs(particles[0].vy) < 1.0e-9);
    assert(fabs(polygon_particles_omega(0)) < 1.0e-9);
    assert(fabs(polygon_particles_rotational_energy()) < 1.0e-10);
    /* Both face vertices contact the curved wall.  The wall-pressure signal
       is the sum of their outward-normal impulses, not merely the x component
       of the particle momentum change. */
    const double expected_wall_impulse = 2.0/sqrt(0.995);
    assert(fabs((cumulativeWallImpulse - wall_impulse_before) -
                expected_wall_impulse) < 1.0e-8);
    freeArrays();

    configure_square_particles();
    initialize_centered_zero_clearance_square(1.0);
    assert(simulationHasFault());
    assert(strstr(simulationFaultMessage(), "more than two") != NULL);
    freeArrays();

    /* A regular polygon can rotate in place when its circumcircle exactly
       matches the wall.  All radial contact speeds are zero, so the
       overconstrained-translation guard must not reject this motion. */
    configure_spinning_square_particles();
    init_genrand(9217UL);
    initialize_centered_zero_clearance_square(0.0);
    assert(fabs(polygon_particles_omega(0)) > 1.0e-6);
    assert(!simulationHasFault());
    assert(findNextEvent() == NULL);
    freeArrays();

    configure_spinning_square_particles();
    initialize_spinning_face_at_circular_wall();
    const double spinning_initial_energy = polygon_system_energy();
    assert(spinning_face_maximum_outward_speed() > 1.0e-6);
    process_until_first_collision();
    assert(fabs(t) < 1.0e-12);
    assert(spinning_face_maximum_outward_speed() <= 1.0e-9*
           (1.0 + fabs(polygon_particles_omega(0))));
    assert(fabs(polygon_system_energy() - spinning_initial_energy) <=
           2.0e-9*(1.0 + spinning_initial_energy));
    freeArrays();

    configure_square_particles();
    initialize_many_for_calendar_stress();
    const double many_initial_energy = polygon_system_energy();
    int processed_events = 0;
    while (ncol < 500 && processed_events < 50000) {
        process_next_polygon_event();
        processed_events++;
        if ((processed_events % 127) == 0) {
            assert(fabs(polygon_system_energy() - many_initial_energy) <=
                   2.0e-8*(1.0 + many_initial_energy));
        }
    }
    assert(ncol >= 500);
    assert(fabs(polygon_system_energy() - many_initial_energy) <=
           2.0e-8*(1.0 + many_initial_energy));
    freeArrays();

    puts("polygon particle calendar integration: ok");
    return 0;
}
