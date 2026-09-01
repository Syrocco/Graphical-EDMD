#define _POSIX_C_SOURCE 200809L

#include "EDMD.h"
#include "sim_config.h"
#include "sim_lifecycle.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

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
    if (!isfinite(actual) || fabs(actual - expected) > 1e-13*scale)
        fail(message);
}

static void default_config(SimulationConfig *config,
                           BuildCapabilities *capabilities)
{
    *capabilities = simulation_build_capabilities();
    simulation_config_set_defaults(config, capabilities);
}

static void expect_error_field(const SimulationConfig *config,
                               const BuildCapabilities *capabilities,
                               const char *field, const char *message)
{
    ConfigIssue issues[128];
    const size_t count = simulation_config_validate(
        config, capabilities, issues, sizeof(issues)/sizeof(issues[0]));
    const size_t stored = count < sizeof(issues)/sizeof(issues[0])
        ? count : sizeof(issues)/sizeof(issues[0]);
    for (size_t index = 0; index < stored; ++index) {
        if (issues[index].severity == CONFIG_ISSUE_ERROR &&
            strcmp(issues[index].field, field) == 0) {
            return;
        }
    }
    fail(message);
}

static void expect_no_validation_errors(const SimulationConfig *config,
                                        const BuildCapabilities *capabilities,
                                        const char *message)
{
    ConfigIssue issues[128];
    const size_t count = simulation_config_validate(
        config, capabilities, issues, sizeof(issues)/sizeof(issues[0]));
    const size_t stored = count < sizeof(issues)/sizeof(issues[0])
        ? count : sizeof(issues)/sizeof(issues[0]);
    for (size_t index = 0; index < stored; ++index) {
        if (issues[index].severity == CONFIG_ISSUE_ERROR) {
            fprintf(stderr, "%s (%s: %s)\n", message,
                    issues[index].field, issues[index].message);
            failures++;
            return;
        }
    }
}

static void test_3d_geometry_capabilities(void)
{
    SimulationConfig config;
    BuildCapabilities capabilities;

    default_config(&config, &capabilities);
    expect_true(capabilities.shear_3d,
                "Lees-Edwards capability is not exposed");
    expect_true(capabilities.full_3d_thermostats,
                "complete thermostat capability is not exposed");
    expect_true(capabilities.lattice_initializers_3d,
                "HCP initializer capability is not exposed");
    expect_true(capabilities.wall_z == (THREE_D != 0),
                "z-wall capability does not match dimensionality");

    config.initialization.mode = SIM_INIT_HEX_LATTICE;
    config.initialization.composition = SIM_COMPOSITION_MONODISPERSE;
    config.initialization.particle_count = 256;
    config.initialization.packing_fraction = 0.5;
    expect_no_validation_errors(&config, &capabilities,
                                "dimension-complete hex/HCP lattice was rejected");

    config.initialization.packing_fraction = capabilities.dimensions == 3
        ? 0.745 : 0.91;
    expect_error_field(&config, &capabilities,
                       "initialization.packing_fraction",
                       "overlapping close-packed lattice was accepted");

    default_config(&config, &capabilities);
    config.geometry.boundary = SIM_BOUNDARY_WALL_XYZ;
    if (THREE_D) {
        expect_no_validation_errors(&config, &capabilities,
                                    "full x/y/z wall box was rejected in 3D");
    } else {
        expect_error_field(&config, &capabilities, "geometry.boundary",
                           "full x/y/z wall box was accepted in 2D");
    }

    default_config(&config, &capabilities);
    capabilities.stress = true;
    capabilities.heat = true;
    config.dynamics.flight = SIM_FLIGHT_LEES_EDWARDS_SHEAR;
    config.dynamics.shear_strain_rate = 0.1;
    expect_no_validation_errors(&config, &capabilities,
                                "Lees-Edwards was rejected with STRESS/HEAT bookkeeping");
}

static void test_polygon_particle_configuration(void)
{
    SimulationConfig config;
    SimulationConfig valid_polygon;
    BuildCapabilities capabilities;

    default_config(&config, &capabilities);
    expect_true(config.particle_shape.mode == SIM_PARTICLE_CIRCLES,
                "circles are not the runtime default particle shape");
    expect_true(config.particle_shape.vertex_count == 6,
                "the default editable polygon is not a hexagon");
    expect_true(config.particle_shape.random_orientation,
                "polygon orientations are not randomized by default");
    expect_close(config.particle_shape.angular_temperature, 1.0,
                 "default polygon angular temperature is incorrect");
    expect_close(simulation_particle_shape_unit_area(&config.particle_shape),
                 M_PI,
                 "a unit-radius circle does not report area pi");
    for (int vertex = 6; vertex < SIM_POLYGON_MAX_VERTICES; ++vertex) {
        expect_close(config.particle_shape.vertices[vertex][0], 0.0,
                     "unused polygon X coordinate was not zero initialized");
        expect_close(config.particle_shape.vertices[vertex][1], 0.0,
                     "unused polygon Y coordinate was not zero initialized");
    }
    expect_true(capabilities.convex_polygon_particles ==
                    (G != 0 && POLYGON_PARTICLES != 0 && THREE_D == 0),
                "polygon capability does not match the build macros");

    /* Exercise the validator independently of this test binary's optional
       feature flag; graphical frontends use the same capability contract. */
    capabilities.dimensions = 2;
    capabilities.convex_polygon_particles = true;
    capabilities.additive_contacts = true;
    capabilities.stress = false;
    capabilities.heat = false;
    capabilities.rough_particles = false;
    config.particle_shape.mode = SIM_PARTICLE_CONVEX_POLYGONS;
    config.dynamics.injection = SIM_INJECTION_NONE;
    config.observables.pressure = false;
    config.output.periodic_dumps = false;
    valid_polygon = config;
    expect_no_validation_errors(&valid_polygon, &capabilities,
                                "a normalized regular convex polygon was rejected");
    expect_close(simulation_particle_shape_unit_area(
                     &valid_polygon.particle_shape),
                 3.0*sqrt(3.0)/2.0,
                 "the unit-circumradius regular hexagon area is incorrect");

    /* Reconstruct the exact case shown in the dense-system report.  A regular
       unit-circumradius triangle has area 3*sqrt(3)/4, so N=1000 in the box
       derived at the reported phi=3/8 must have L=sqrt(2000*sqrt(3)) and its
       occupied-area fraction must be exactly 0.375. */
    {
        SimParticleShapeConfig triangle = valid_polygon.particle_shape;
        const double triangle_area = 3.0*sqrt(3.0)/4.0;
        const double requested_fraction = 3.0/8.0;
        const double box_length = sqrt(1000.0*triangle_area/
                                       requested_fraction);
        triangle.vertex_count = 3;
        memset(triangle.vertices, 0, sizeof(triangle.vertices));
        for (int vertex = 0; vertex < triangle.vertex_count; ++vertex) {
            const double angle = M_PI/2.0 +
                2.0*M_PI*vertex/triangle.vertex_count;
            triangle.vertices[vertex][0] = cos(angle);
            triangle.vertices[vertex][1] = sin(angle);
        }
        expect_close(simulation_particle_shape_unit_area(&triangle),
                     triangle_area,
                     "the unit-circumradius triangle area is incorrect");
        expect_close(box_length, sqrt(2000.0*sqrt(3.0)),
                     "the triangle box length does not follow occupied area / box area");
        expect_close(simulation_particle_shape_packing_fraction_2d(
                         &triangle, 1000.0, box_length*box_length),
                     requested_fraction,
                     "monodisperse triangle packing is not sum(area)/box area");

        /* The same invariant must use the exact squared radii for mixtures,
           rather than disk area or a nominal species count. */
        expect_close(simulation_particle_shape_packing_fraction_2d(
                         &triangle, 2.0 + 2.0*0.25, 10.0),
                     triangle_area*2.5/10.0,
                     "binary triangle packing ignored the realized radii");
    }

    /* A custom canonical diamond has unit area two.  This guards against
       silently substituting pi or a regular-polygon formula for editor-made
       bodies. */
    {
        SimParticleShapeConfig diamond = valid_polygon.particle_shape;
        const double vertices[4][2] = {
            {1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}, {0.0, -1.0}
        };
        diamond.vertex_count = 4;
        memset(diamond.vertices, 0, sizeof(diamond.vertices));
        memcpy(diamond.vertices, vertices, sizeof(vertices));
        expect_close(simulation_particle_shape_packing_fraction_2d(
                         &diamond, 3.25, 13.0),
                     0.5,
                     "custom polygon packing is not true scaled body area / box area");
        expect_true(isnan(simulation_particle_shape_packing_fraction_2d(
                        &diamond, 3.25, 0.0)),
                    "zero container area produced a finite packing fraction");
    }

    config = valid_polygon;
    config.initialization.packing_fraction = 0.9;
    expect_no_validation_errors(&config, &capabilities,
                                "polygon-native growth retained the obsolete circumscribed-disk packing cap");

    config.initialization.packing_fraction = 1.000001;
    expect_error_field(&config, &capabilities,
                       "initialization.packing_fraction",
                       "packing fraction above the physical area bound was accepted");

    config = valid_polygon;
    config.initialization.mode = SIM_INIT_HEX_LATTICE;
    config.initialization.composition = SIM_COMPOSITION_MONODISPERSE;
    config.initialization.packing_fraction =
        simulation_particle_shape_unit_area(&config.particle_shape)/
        (2.0*sqrt(3.0)) + 1.0e-6;
    expect_error_field(&config, &capabilities,
                       "initialization.packing_fraction",
                       "polygon lattice accepted overlapping circumdisks");

    config = valid_polygon;
    config.particle_shape.vertex_count = 2;
    expect_error_field(&config, &capabilities,
                       "particle_shape.vertex_count",
                       "a two-vertex particle shape was accepted");

    config = valid_polygon;
    config.particle_shape.vertices[1][0] = 0.0;
    config.particle_shape.vertices[1][1] = 0.0;
    expect_error_field(&config, &capabilities, "particle_shape.vertices",
                       "a non-strictly-convex particle shape was accepted");

    config = valid_polygon;
    for (int vertex = 0; vertex < config.particle_shape.vertex_count;
         ++vertex) {
        config.particle_shape.vertices[vertex][0] *= 0.5;
        config.particle_shape.vertices[vertex][1] *= 0.5;
    }
    expect_error_field(&config, &capabilities, "particle_shape.vertices",
                       "a noncanonical polygon circumradius was accepted");

    config = valid_polygon;
    config.particle_shape.angular_temperature = -0.01;
    expect_error_field(&config, &capabilities,
                       "particle_shape.angular_temperature",
                       "a negative angular temperature was accepted");

    config = valid_polygon;
    capabilities.convex_polygon_particles = false;
    expect_error_field(&config, &capabilities, "particle_shape.mode",
                       "polygon mode was accepted without compiled support");
    capabilities.convex_polygon_particles = true;

    capabilities.dimensions = 3;
    expect_error_field(&config, &capabilities, "particle_shape.mode",
                       "polygon mode was accepted in 3D");
    capabilities.dimensions = 2;

    config = valid_polygon;
    config.initialization.mode = SIM_INIT_DUMP_FILE;
    expect_error_field(&config, &capabilities, "initialization.mode",
                       "polygon mode accepted a legacy dump restart");

    config = valid_polygon;
    config.dynamics.flight = SIM_FLIGHT_CONSTANT_FIELD;
    expect_error_field(&config, &capabilities, "dynamics.flight",
                       "polygon mode accepted non-ballistic flight");

    config = valid_polygon;
    config.dynamics.injection = SIM_INJECTION_DELTA_NORMAL;
    expect_error_field(&config, &capabilities, "dynamics.injection",
                       "polygon mode accepted collision injection");

    config = valid_polygon;
    config.driving.thermostat = SIM_THERMOSTAT_LANGEVIN;
    config.driving.euler_noise = false;
    expect_no_validation_errors(
        &config, &capabilities,
        "polygon mode rejected the translation-and-rotation Langevin thermostat");

    config.driving.euler_noise = true;
    expect_error_field(&config, &capabilities, "driving.euler_noise",
                       "polygon mode accepted the inexact Euler Langevin update");

    config = valid_polygon;
    config.driving.thermostat = SIM_THERMOSTAT_REMOVE_COM;
    expect_error_field(&config, &capabilities, "driving.thermostat",
                       "polygon mode accepted a translation-only thermostat");

    config = valid_polygon;
    config.potential.enabled = true;
    expect_error_field(&config, &capabilities, "potential.enabled",
                       "polygon mode accepted a square potential");

    config = valid_polygon;
    config.intruder.enabled = true;
    expect_error_field(&config, &capabilities, "intruder.enabled",
                       "polygon mode accepted the separate intruder");

    config = valid_polygon;
    config.geometry.boundary = SIM_BOUNDARY_WALL_X;
    config.geometry.thermal_wall = true;
    expect_error_field(&config, &capabilities, "geometry.thermal_wall",
                       "polygon mode accepted a thermal wall");

    config = valid_polygon;
    config.geometry.boundary = SIM_BOUNDARY_WALL_Y;
    config.geometry.shear_wall = true;
    expect_error_field(&config, &capabilities, "geometry.shear_wall",
                       "polygon mode accepted a moving shear wall");

    config = valid_polygon;
    capabilities.additive_contacts = false;
    expect_error_field(&config, &capabilities, "particle_shape.mode",
                       "polygon mode accepted NONADDITIVE=0");
    capabilities.additive_contacts = true;

    capabilities.stress = true;
    expect_error_field(&config, &capabilities, "particle_shape.mode",
                       "polygon mode accepted incomplete STRESS bookkeeping");
    capabilities.stress = false;

    config = valid_polygon;
    config.observables.pressure = true;
    expect_error_field(&config, &capabilities, "observables.pressure",
                       "polygon mode accepted the disk-oriented saved pressure column");

    config = valid_polygon;
    config.geometry.center_membrane = true;
    config.observables.thermo_enabled = true;
    expect_error_field(&config, &capabilities,
                       "observables.thermo_enabled",
                       "polygon mode accepted disk-oriented membrane thermodynamics");

    config = valid_polygon;
    config.observables.pressure_field = true;
    expect_error_field(&config, &capabilities,
                       "observables.pressure_field",
                       "polygon mode accepted an unconnected local pressure field");

    config = valid_polygon;
    config.observables.temperature_profile = true;
    expect_error_field(&config, &capabilities,
                       "observables.temperature_profile",
                       "polygon mode accepted an incomplete rotational heat profile");

    config = valid_polygon;
    config.observables.interface_density = true;
    expect_error_field(&config, &capabilities,
                       "observables.interface_density",
                       "polygon mode accepted disk-area interface density");

    config = valid_polygon;
    config.observables.interface_temperature = true;
    expect_error_field(&config, &capabilities,
                       "observables.interface_temperature",
                       "polygon mode accepted the disk-area interface field");

    config = valid_polygon;
    config.observables.osmosis = true;
    expect_error_field(&config, &capabilities, "observables.osmosis",
                       "polygon mode accepted disk-area osmosis profiles");

    config = valid_polygon;
    config.observables.voronoi_area = true;
    expect_error_field(&config, &capabilities,
                       "observables.voronoi_area",
                       "polygon mode accepted disk-area Voronoi packing");

    config = valid_polygon;
    config.observables.critical_regions = true;
    expect_error_field(&config, &capabilities,
                       "observables.critical_regions",
                       "polygon mode accepted incomplete local energy fields");

    config = valid_polygon;
    config.observables.bubble_stop = true;
    expect_error_field(&config, &capabilities, "observables.bubble_stop",
                       "polygon mode accepted hard-disk bubble thresholds");

    config = valid_polygon;
    config.observables.velocity_profile = true;
    config.observables.pair_correlation = true;
    config.observables.structure_factor = SIM_STRUCTURE_FACTOR_DENSITY;
    config.observables.clusters = true;
    config.observables.bond_order = SIM_BOOP_CUTOFF;
    config.rendering.graphical_structure_factor = true;
    expect_no_validation_errors(
        &config, &capabilities,
        "polygon-safe center/velocity observables were rejected");
}

static void test_thermostat_validation(void)
{
    SimulationConfig config;
    BuildCapabilities capabilities;

    default_config(&config, &capabilities);
    config.driving.thermostat = (SimThermostatMode)3;
    expect_error_field(&config, &capabilities, "driving.thermostat",
                       "reserved thermostat id 3 was accepted");

    default_config(&config, &capabilities);
    config.driving.thermostat = (SimThermostatMode)5;
    expect_error_field(&config, &capabilities, "driving.thermostat",
                       "reserved thermostat id 5 was accepted");

    default_config(&config, &capabilities);
    config.driving.thermostat = (SimThermostatMode)10;
    expect_error_field(&config, &capabilities, "driving.thermostat",
                       "reserved thermostat id 10 was accepted");

    default_config(&config, &capabilities);
    config.driving.thermostat = SIM_THERMOSTAT_VELOCITY_RESCALE;
    config.driving.temperature = 0.0;
    expect_error_field(&config, &capabilities, "driving.temperature",
                       "velocity rescaling accepted zero temperature");

    default_config(&config, &capabilities);
    config.geometry.boundary = SIM_BOUNDARY_WALL_X;
    config.geometry.thermal_wall = true;
    config.driving.temperature = 0.0;
    expect_error_field(&config, &capabilities, "driving.temperature",
                       "thermal wall accepted zero temperature");
}

static void test_3d_lighting_configuration(void)
{
    SimulationConfig config;
    BuildCapabilities capabilities;
    const SimLightingMode modes[] = {
        SIM_LIGHTING_SIMPLE,
        SIM_LIGHTING_STUDIO,
        SIM_LIGHTING_DEPTH_ENHANCED,
    };

    default_config(&config, &capabilities);
    expect_true(config.rendering.lighting_mode == SIM_LIGHTING_STUDIO,
                "Studio lighting is not the default");
    for (size_t index = 0; index < sizeof(modes)/sizeof(modes[0]); ++index) {
        config.rendering.lighting_mode = modes[index];
        expect_no_validation_errors(&config, &capabilities,
                                    "a supported lighting mode was rejected");
    }
    config.rendering.lighting_mode = (SimLightingMode)99;
    expect_error_field(&config, &capabilities,
                       "rendering.lighting_mode",
                       "an unknown lighting mode was accepted");
}

static void test_binary_mixture_round_trip(void)
{
    SimulationConfig config;
    SimulationConfig captured;
    BuildCapabilities capabilities;

    default_config(&config, &capabilities);
    config.initialization.composition = SIM_COMPOSITION_BINARY;
    config.initialization.small_fraction = 0.35;
    simulation_config_apply_globals(&config);
    simulation_config_capture_globals(&captured);
    expect_true(captured.initialization.composition == SIM_COMPOSITION_BINARY,
                "binary composition was lost after apply/capture");
    expect_close(captured.initialization.small_fraction, 0.35,
                 "binary fraction was lost after apply/capture");

    config.initialization.mode = SIM_INIT_SQUARE_LATTICE;
    config.initialization.small_fraction = 0.0;
    simulation_config_apply_globals(&config);
    simulation_config_capture_globals(&captured);
    expect_true(captured.initialization.composition == SIM_COMPOSITION_BINARY,
                "forced binary lattice was captured as monodisperse");

    config.initialization.mode = SIM_INIT_RANDOM_GROWTH;
    config.initialization.small_fraction = 0.0;
    expect_error_field(&config, &capabilities,
                       "initialization.small_fraction",
                       "degenerate zero-fraction binary mixture was accepted");
    config.initialization.small_fraction = 1.0;
    expect_error_field(&config, &capabilities,
                       "initialization.small_fraction",
                       "degenerate one-fraction binary mixture was accepted");
}

static void test_default_time_horizon(void)
{
    SimulationConfig config;
    BuildCapabilities capabilities = simulation_build_capabilities();

    capabilities.graphics = false;
    simulation_config_set_defaults(&config, &capabilities);
    expect_close(config.output.maximum_time, 30000.0,
                 "headless default time horizon is not the historical 30000");

    capabilities.graphics = true;
    simulation_config_set_defaults(&config, &capabilities);
    expect_true(config.output.maximum_time == never,
                "graphical default time horizon is not unbounded");
}

static void test_exposed_numeric_ranges(void)
{
    SimulationConfig config;
    BuildCapabilities capabilities;

    default_config(&config, &capabilities);
    config.initialization.packing_fraction = 1.0;
    config.rendering.target_fps = 1;
    config.rendering.graphical_structure_factor = true;
    config.rendering.structure_factor_sample_cap = 1;
    expect_no_validation_errors(&config, &capabilities,
                                "lower exposed rendering bounds or packing fraction were rejected");

    config.initialization.packing_fraction = 0.5;
    config.rendering.target_fps = 1000;
    config.rendering.structure_factor_sample_cap = 8192;
    expect_no_validation_errors(&config, &capabilities,
                                "upper exposed rendering bounds were rejected");

    config.rendering.target_fps = 1001;
    expect_error_field(&config, &capabilities, "rendering.target_fps",
                       "target FPS above the exposed maximum was accepted");
    config.rendering.target_fps = 1000;
    config.rendering.structure_factor_sample_cap = 8193;
    expect_error_field(&config, &capabilities,
                       "rendering.structure_factor_sample_cap",
                       "structure-factor cap above the exposed maximum was accepted");
}

static void test_velocity_rescale_factor(void)
{
    double factor = -1.0;
    expect_true(velocityRescaleFactor(4.0, 2, 8.0, &factor),
                "valid velocity rescale inputs were rejected");
    expect_close(factor, 2.0, "velocity rescale factor is incorrect");
    expect_true(!velocityRescaleFactor(0.0, 2, 1.0, &factor),
                "zero kinetic energy was accepted for velocity rescaling");
    expect_true(!velocityRescaleFactor(NAN, 2, 1.0, &factor),
                "NaN kinetic energy was accepted for velocity rescaling");
    expect_true(!velocityRescaleFactor(1.0, 2, 0.0, &factor),
                "zero target temperature was accepted for velocity rescaling");
    expect_true(!velocityRescaleFactor(1.0, 2, INFINITY, &factor),
                "infinite target temperature was accepted for velocity rescaling");
    expect_true(!velocityRescaleFactor(1.0, 0, 1.0, &factor),
                "zero particle count was accepted for velocity rescaling");
}

static int open_temporary_path(char path[SIM_CONFIG_PATH_CAPACITY])
{
    char pattern[] = "/tmp/edmd-config-test-XXXXXX";
    const int descriptor = mkstemp(pattern);
    if (descriptor < 0) return -1;
    (void)snprintf(path, SIM_CONFIG_PATH_CAPACITY, "%s", pattern);
    return descriptor;
}

static bool write_dump(char path[SIM_CONFIG_PATH_CAPACITY],
                       double x_lower, double x_upper,
                       double y_lower, double y_upper,
                       double particle_x, double particle_y,
                       double radius)
{
    const int descriptor = open_temporary_path(path);
    if (descriptor < 0) return false;
    FILE *stream = fdopen(descriptor, "w");
    if (stream == NULL) {
        (void)close(descriptor);
        (void)unlink(path);
        return false;
    }
    fprintf(stream,
            "ITEM: TIMESTEP\n0\n"
            "ITEM: NUMBER OF ATOMS\n1\n"
            "ITEM: BOX BOUNDS pp pp pp\n"
            "%.17g %.17g\n%.17g %.17g\n",
            x_lower, x_upper, y_lower, y_upper);
#if THREE_D
    fprintf(stream,
            "0 10\n"
            "ITEM: ATOMS type x y z vx vy vz radius m\n"
            "1 %.17g %.17g 5 0 0 0 %.17g 1\n",
            particle_x, particle_y, radius);
#else
    fprintf(stream,
            "0 1\n"
            "ITEM: ATOMS type x y vx vy radius m\n"
            "1 %.17g %.17g 0 0 %.17g 1\n",
            particle_x, particle_y, radius);
#endif
    if (fclose(stream) != 0) {
        (void)unlink(path);
        return false;
    }
    return true;
}

#if THREE_D
static bool write_dump_3d_z(char path[SIM_CONFIG_PATH_CAPACITY],
                            double particle_z, double radius)
{
    const int descriptor = open_temporary_path(path);
    if (descriptor < 0) return false;
    FILE *stream = fdopen(descriptor, "w");
    if (stream == NULL) {
        (void)close(descriptor);
        (void)unlink(path);
        return false;
    }
    fprintf(stream,
            "ITEM: TIMESTEP\n0\n"
            "ITEM: NUMBER OF ATOMS\n1\n"
            "ITEM: BOX BOUNDS pp pp pp\n"
            "0 10\n0 10\n0 10\n"
            "ITEM: ATOMS type x y z vx vy vz radius m\n"
            "1 5 5 %.17g 0 0 0 %.17g 1\n",
            particle_z, radius);
    if (fclose(stream) != 0) {
        (void)unlink(path);
        return false;
    }
    return true;
}
#endif

static void configure_file_restart(SimulationConfig *config,
                                   BuildCapabilities *capabilities,
                                   SimInitMode mode, const char *path,
                                   SimBoundaryMode boundary)
{
    default_config(config, capabilities);
    config->initialization.mode = mode;
    config->geometry.boundary = boundary;
    (void)snprintf(config->initialization.load_path,
                   sizeof(config->initialization.load_path), "%s", path);
}

static void expect_restart_error(const SimulationConfig *config,
                                 const char *needle, const char *message)
{
    char error[1024] = {0};
    const bool restarted = simulation_restart_core(config, error,
                                                    sizeof(error));
    if (restarted || strstr(error, needle) == NULL) {
        fprintf(stderr, "%s (error was: %s)\n", message, error);
        failures++;
    }
}

static void test_dump_origin_and_wall_containment(void)
{
    char path[SIM_CONFIG_PATH_CAPACITY];
    SimulationConfig config;
    BuildCapabilities capabilities;

    expect_true(write_dump(path, -5.0, 5.0, 0.0, 10.0,
                           0.0, 5.0, 0.5),
                "could not create nonzero-origin dump fixture");
    configure_file_restart(&config, &capabilities, SIM_INIT_DUMP_FILE,
                           path, SIM_BOUNDARY_PERIODIC);
    expect_restart_error(&config, "unsupported nonzero x-box origin",
                         "nonzero-origin dump was not rejected clearly");
    (void)unlink(path);

    expect_true(write_dump(path, 0.0, 10.0, 0.0, 10.0,
                           0.25, 5.0, 0.5),
                "could not create X-wall containment dump fixture");
    configure_file_restart(&config, &capabilities, SIM_INIT_DUMP_FILE,
                           path, SIM_BOUNDARY_WALL_X);
    expect_restart_error(&config, "outside the selected X walls",
                         "dump particle outside X walls was accepted");
    configure_file_restart(&config, &capabilities, SIM_INIT_DUMP_FILE,
                           path, SIM_BOUNDARY_PERIODIC);
    expect_restart_error(&config, "only available in a graphical build",
                         "periodic dump was rejected by wall containment");
    (void)unlink(path);

    expect_true(write_dump(path, 0.0, 10.0, 0.0, 10.0,
                           5.0, 9.75, 0.5),
                "could not create Y-wall containment dump fixture");
    configure_file_restart(&config, &capabilities, SIM_INIT_DUMP_FILE,
                           path, SIM_BOUNDARY_WALL_Y);
    expect_restart_error(&config, "outside the selected Y walls",
                         "dump particle outside Y walls was accepted");
    (void)unlink(path);

#if THREE_D
    expect_true(write_dump_3d_z(path, 9.75, 0.5),
                "could not create Z-wall containment dump fixture");
    configure_file_restart(&config, &capabilities, SIM_INIT_DUMP_FILE,
                           path, SIM_BOUNDARY_WALL_XYZ);
    expect_restart_error(&config, "outside the selected Z walls",
                         "dump particle outside Z walls was accepted");
    (void)unlink(path);
#endif
}

int main(void)
{
    test_3d_geometry_capabilities();
    test_polygon_particle_configuration();
    test_thermostat_validation();
    test_3d_lighting_configuration();
    test_binary_mixture_round_trip();
    test_default_time_horizon();
    test_exposed_numeric_ranges();
    test_velocity_rescale_factor();
    test_dump_origin_and_wall_containment();
    if (failures != 0) {
        fprintf(stderr, "%d configuration regression(s) failed\n", failures);
        return 1;
    }
    printf("configuration regressions passed (%dD)\n", THREE_D ? 3 : 2);
    return 0;
}
