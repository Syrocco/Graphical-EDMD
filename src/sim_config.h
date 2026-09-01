#ifndef EDMD_SIM_CONFIG_H
#define EDMD_SIM_CONFIG_H

/*
 * Declarative simulation configuration shared by the command line and GUI.
 *
 * This header deliberately has no raylib, ImGui, or EDMD-internal includes so
 * it remains usable from both C and C++ front ends.  Options which are not yet
 * mutable in the legacy engine are still represented: BuildCapabilities and
 * simulation_config_validate() let a front end display them as unavailable
 * instead of silently accepting an ineffective combination.
 */

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SIM_CONFIG_SCHEMA_VERSION 3
#define SIM_CONFIG_PATH_CAPACITY 512
#define SIM_CONFIG_FIELD_CAPACITY 64
#define SIM_CONFIG_MESSAGE_CAPACITY 256
#define SIM_POLYGON_MAX_VERTICES 16

typedef enum SimParticleShapeMode {
    SIM_PARTICLE_CIRCLES = 0,
    SIM_PARTICLE_CONVEX_POLYGONS
} SimParticleShapeMode;

/* A normalized, counter-clockwise convex particle template.  The runtime
   scales this common shape by each particle's existing radius. */
typedef struct SimParticleShapeConfig {
    SimParticleShapeMode mode;
    int vertex_count;
    double vertices[SIM_POLYGON_MAX_VERTICES][2];
    bool random_orientation;
    double angular_temperature;
} SimParticleShapeConfig;

typedef enum SimInitMode {
    SIM_INIT_RANDOM_GROWTH = 0,
    SIM_INIT_DUMP_FILE,
    SIM_INIT_SQUARE_LATTICE,
    SIM_INIT_HEX_LATTICE,
    SIM_INIT_MEMBRANE_OSMOSIS
} SimInitMode;

typedef enum SimCompositionMode {
    SIM_COMPOSITION_MONODISPERSE = 0,
    SIM_COMPOSITION_BINARY,
    SIM_COMPOSITION_POLYDISPERSE
} SimCompositionMode;

typedef enum SimBoxLengthMode {
    SIM_BOX_FROM_PACKING_FRACTION = 0,
    SIM_BOX_FIXED_LX,
    SIM_BOX_FIXED_LY
} SimBoxLengthMode;

typedef enum SimBoundaryMode {
    SIM_BOUNDARY_PERIODIC = 0,
    SIM_BOUNDARY_WALL_X,
    SIM_BOUNDARY_WALL_Y,
    SIM_BOUNDARY_WALL_XY,
    SIM_BOUNDARY_CIRCULAR_SPHERICAL,
    /* A genuinely closed cuboid: x, y, and z each have two planar faces.
       This mode is rejected by 2D builds. */
    SIM_BOUNDARY_WALL_XYZ
} SimBoundaryMode;

typedef enum SimFlightMode {
    SIM_FLIGHT_BALLISTIC = 0,
    SIM_FLIGHT_CONSTANT_FIELD,
    SIM_FLIGHT_LINEAR_DAMPING,
    SIM_FLIGHT_LEES_EDWARDS_SHEAR
} SimFlightMode;

typedef enum SimThermostatMode {
    SIM_THERMOSTAT_NONE = 0,
    SIM_THERMOSTAT_LANGEVIN = 1,
    SIM_THERMOSTAT_VELOCITY_RESCALE = 2,
    SIM_THERMOSTAT_REMOVE_COM = 4,
    SIM_THERMOSTAT_MAXWELL_REDRAW = 8
} SimThermostatMode;

typedef enum SimInjectionMode {
    SIM_INJECTION_NONE = 0,
    SIM_INJECTION_DELTA_NORMAL,
    SIM_INJECTION_DELTA_TANGENT,
    SIM_INJECTION_DOUBLE_DELTA_NORMAL,
    SIM_INJECTION_DOUBLE_DELTA_TANGENT,
    SIM_INJECTION_EVOLVING_DELTA,
    SIM_INJECTION_COLLISION_ENERGY,
    SIM_INJECTION_ACTIVE_TANGENTIAL,
    SIM_INJECTION_ROUGH_ROTATIONAL
} SimInjectionMode;

typedef enum SimRestitutionMode {
    SIM_RESTITUTION_CONSTANT = 0,
    SIM_RESTITUTION_EXPONENTIAL
} SimRestitutionMode;

typedef enum SimPotentialSpeciesMode {
    SIM_POTENTIAL_ALL_PAIRS = 0,
    SIM_POTENTIAL_CHARGED,
    SIM_POTENTIAL_SAME_SPECIES
} SimPotentialSpeciesMode;

typedef enum SimStructureFactorMode {
    SIM_STRUCTURE_FACTOR_OFF = 0,
    SIM_STRUCTURE_FACTOR_DENSITY,
    SIM_STRUCTURE_FACTOR_VELOCITY
} SimStructureFactorMode;

typedef enum SimBoopMode {
    SIM_BOOP_OFF = 0,
    SIM_BOOP_VORONOI = 1,
    SIM_BOOP_CUTOFF = 2
} SimBoopMode;

typedef enum SimColorMode {
    SIM_COLOR_BY_TYPE = 0,
    SIM_COLOR_BY_COLLISIONS,
    SIM_COLOR_BY_SPEED,
    SIM_COLOR_BY_HEX_ORDER,
    SIM_COLOR_BY_SQUARE_ORDER,
    SIM_COLOR_BY_RADIUS,
    SIM_COLOR_BY_CHARGE
} SimColorMode;

typedef enum SimColorPalette {
    SIM_PALETTE_PLASMA = 0,
    SIM_PALETTE_VIRIDIS,
    SIM_PALETTE_COPPER
} SimColorPalette;

typedef enum SimSphereQuality {
    SIM_SPHERE_QUALITY_AUTO = 0,
    SIM_SPHERE_QUALITY_LOW,
    SIM_SPHERE_QUALITY_MEDIUM,
    SIM_SPHERE_QUALITY_HIGH,
    SIM_SPHERE_QUALITY_ULTRA
} SimSphereQuality;

/* Portable, single-pass lighting profiles for the instanced 3D spheres.
   The depth-enhanced profile adds view-dependent depth cues, not screen-space
   ambient occlusion or ray tracing. */
typedef enum SimLightingMode {
    SIM_LIGHTING_SIMPLE = 0,
    SIM_LIGHTING_STUDIO,
    SIM_LIGHTING_DEPTH_ENHANCED
} SimLightingMode;

typedef enum ConfigIssueSeverity {
    CONFIG_ISSUE_INFO = 0,
    CONFIG_ISSUE_WARNING,
    CONFIG_ISSUE_ERROR
} ConfigIssueSeverity;

typedef struct ConfigIssue {
    ConfigIssueSeverity severity;
    char field[SIM_CONFIG_FIELD_CAPACITY];
    char message[SIM_CONFIG_MESSAGE_CAPACITY];
} ConfigIssue;

typedef struct BuildCapabilities {
    int dimensions;
    bool graphics;
    bool stress;
    bool heat;
    bool rough_particles;
    bool intruder;
    /* True when hard-contact distance is r_i + r_j.  The legacy build flag
       selecting this rule is unfortunately named NONADDITIVE. */
    bool additive_contacts;
    bool runtime_physics_toggles;
    bool runtime_basic_injection_toggles;
    bool runtime_advanced_injection_toggles;
    bool runtime_exponential_restitution;
    bool runtime_charged_potential;
    bool runtime_same_species_potential;
    bool runtime_output_lifecycle;
    bool runtime_output_column_toggles;
    bool runtime_intruder_parameters;
    bool convex_polygon_particles;

    bool spherical_confinement;
    bool tangent_delta_3d;
    bool wall_z;
    /* True only when every exposed thermostat has a complete 3D definition. */
    bool full_3d_thermostats;
    /* Refers to Lees-Edwards shear, not dimension-safe moving Y walls. */
    bool shear_3d;
    bool lattice_initializers_3d;
    bool intruder_3d;
    bool pressure_field_3d;
    bool projected_xy_observables;

    bool queue_verification;
    bool overlap_verification;
    bool collision_accounting_verification;
} BuildCapabilities;

typedef struct SimInitializationConfig {
    SimInitMode mode;
    char load_path[SIM_CONFIG_PATH_CAPACITY];
    int particle_count;
    double packing_fraction;
    SimCompositionMode composition;
    double size_ratio;
    double small_fraction;
    double polydisperse_a;
    double polydisperse_b;
    double initial_energy;
    int random_seed;

    SimBoxLengthMode box_length_mode;
    double aspect_ratio;
    double box_lx;
    double box_ly;
    double box_lz;

    double liquid_packing_fraction;
    double gas_packing_fraction;
    double hex_hole_radius;

    bool automatic_growth_rate;
    double growth_rate;
} SimInitializationConfig;

typedef struct SimGeometryConfig {
    SimBoundaryMode boundary;
    bool center_membrane;
    double wall_normal_restitution;
    double wall_tangent_restitution;
    bool thermal_wall;
    bool shear_wall;
    double shear_wall_speed;
} SimGeometryConfig;

typedef struct SimDynamicsConfig {
    SimFlightMode flight;
    double constant_field;
    double damping_gamma;
    double shear_strain_rate;

    SimRestitutionMode restitution_model;
    double particle_restitution;
    double exponential_velocity_scale;
    double exponential_amplitude;

    SimInjectionMode injection;
    double delta;
    double delta_small;
    double delta_large;
    double collision_memory_time;
    double evolving_delta_tau;
    double evolving_delta_max;

    double injected_energy;
    double background_injected_energy;
    double injection_recovery_time;
    double injection_exponent;
    bool randomize_injected_energy;
    bool simplified_energy_injection;

    /* The legacy TANGENTIAL build always uses rough rotational contacts, but
       may additionally apply its older translational tangential-delta kick. */
    bool rough_tangential_delta;
    double rough_tangential_restitution;
    double rough_angular_kick;
} SimDynamicsConfig;

typedef struct SimDrivingConfig {
    SimThermostatMode thermostat;
    double thermostat_interval;
    double temperature;
    double gamma;
    bool euler_noise;
} SimDrivingConfig;

typedef struct SimPotentialConfig {
    bool enabled;
    double energy;
    double range;
    SimPotentialSpeciesMode species_mode;
    double positive_charge_fraction;
    double neutral_charge_fraction;
} SimPotentialConfig;

typedef struct SimObservablesConfig {
    bool thermo_enabled;
    bool collision_count;
    bool energy;
    bool pressure;
    bool center_of_mass_energy;
    bool kurtosis;
    bool sonine;
    bool mean_squared_displacement;

    bool pressure_field;
    double pressure_field_spacing;
    bool temperature_profile;
    int temperature_profile_slabs;
    bool velocity_profile;
    int velocity_profile_bins;

    bool interface_density;
    bool interface_temperature;
    bool osmosis;
    double osmosis_slice_width;

    SimStructureFactorMode structure_factor;
    bool intermediate_scattering;
    double structure_factor_qmax;
    bool transverse_velocity_correlation;
    int tvacf_mode_m;
    int tvacf_mode_n;
    int tvacf_samples;
    double tvacf_time_span;
    int tvacf_checkpoint_interval;

    bool clusters;
    bool cluster_histogram;
    int cluster_stop_size;
    int cluster_mode;
    int cluster_neighbor_threshold;
    double cluster_cutoff;

    SimBoopMode bond_order;
    bool voronoi_area;
    bool pair_correlation;
    bool bond_order_pair_correlation;
    bool g6_pair_correlation;
    bool velocity_distribution;

    bool critical_regions;
    bool critical_snapshot_labels;
    double critical_fraction;
    bool bubble_stop;
} SimObservablesConfig;

typedef struct SimOutputConfig {
    double maximum_time;
    bool periodic_dumps;
    double first_dump_time;
    double dump_interval;
    double next_dump_time_override;
    double first_thermo_time;
    double thermo_interval;
    int thermo_average_samples;
    bool logarithmic_thermo_spacing;
    double logarithmic_spacing_base;
    bool reduced_dump;
    bool unwrap_coordinates;
    int version;
    char output_path[SIM_CONFIG_PATH_CAPACITY];
} SimOutputConfig;

typedef struct SimRenderingConfig {
    bool start_running;
    int target_fps;
    bool velocity_arrows;
    SimColorMode color_mode;
    SimColorPalette palette;
    bool graphical_structure_factor;
    bool show_fps;
    bool show_axes;
    bool show_boundary;
    bool show_selection_outline;
    SimLightingMode lighting_mode;
    SimSphereQuality sphere_quality;
    double scene_supersampling;
    double camera_fov_degrees;
    int structure_factor_sample_cap;
    double structure_factor_refresh_seconds;
    double color_refresh_seconds;
    uint32_t type0_rgba;
    uint32_t type1_rgba;
    uint32_t background_rgba;
    uint32_t boundary_rgba;
} SimRenderingConfig;

typedef struct SimIntruderConfig {
    bool enabled;
    int shape;
    double shape_parameters[4];
    double mass;
    double moment_of_inertia;
    double initial_x;
    double initial_y;
    double initial_vx;
    double initial_vy;
    double initial_angle;
    double initial_omega;
    double restitution;
    double delta;
    double bath_particle_radius;
    bool rethermalize_after_collision;
    double display_particle_radius;
    bool kinetic_model_measurement;
    double kinetic_model_burn_time;
} SimIntruderConfig;

typedef struct SimulationConfig {
    int schema_version;
    SimParticleShapeConfig particle_shape;
    SimInitializationConfig initialization;
    SimGeometryConfig geometry;
    SimDynamicsConfig dynamics;
    SimDrivingConfig driving;
    SimPotentialConfig potential;
    SimObservablesConfig observables;
    SimOutputConfig output;
    SimRenderingConfig rendering;
    SimIntruderConfig intruder;
} SimulationConfig;

BuildCapabilities simulation_build_capabilities(void);
void simulation_config_set_defaults(SimulationConfig *config,
                                    const BuildCapabilities *capabilities);
void simulation_config_capture_globals(SimulationConfig *config);
void simulation_config_apply_globals(const SimulationConfig *config);

/* Physical area of a unit-radius 2D particle.  Circle radius is the usual
 * geometric radius; polygon radius is the canonical template circumradius.
 * Returns NaN for an invalid/unknown shape description. */
double simulation_particle_shape_unit_area(
    const SimParticleShapeConfig *particle_shape);

/* Occupied-area fraction for similarly shaped 2D particles whose individual
 * circumradii have squared sum `sum_radius_squared`.  Keeping this calculation
 * in the shared configuration layer lets box sizing, diagnostics, and tests
 * use the exact same definition: sum(particle areas) / container area. */
double simulation_particle_shape_packing_fraction_2d(
    const SimParticleShapeConfig *particle_shape,
    double sum_radius_squared,
    double container_area);

/* Returns the total number of issues, which may exceed issue_capacity. */
size_t simulation_config_validate(const SimulationConfig *config,
                                  const BuildCapabilities *capabilities,
                                  ConfigIssue *issues,
                                  size_t issue_capacity);

bool simulation_config_equal(const SimulationConfig *left,
                             const SimulationConfig *right);

#ifdef __cplusplus
}
#endif

#endif
