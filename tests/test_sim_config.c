#define _POSIX_C_SOURCE 200809L

#include "EDMD.h"
#include "sim_config.h"
#include "sim_lifecycle.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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
    test_thermostat_validation();
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
