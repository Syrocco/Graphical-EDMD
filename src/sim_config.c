#include "sim_config.h"

#include "EDMD.h"

#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#ifndef G
#define G 0
#endif
#ifndef STRESS
#define STRESS 0
#endif
#ifndef HEAT
#define HEAT 0
#endif
#ifndef NONADDITIVE
#define NONADDITIVE 1
#endif
#ifndef EDMD_VERIFY_QUEUE
#define EDMD_VERIFY_QUEUE 0
#endif
#ifndef EDMD_VERIFY_OVERLAPS
#define EDMD_VERIFY_OVERLAPS 0
#endif
#ifndef EDMD_VERIFY_COLLISION_ACCOUNTING
#define EDMD_VERIFY_COLLISION_ACCOUNTING 0
#endif

enum {
    /* eventListInit() indexes 2*N particle events plus ten auxiliary slots. */
    EVENT_LIST_AUXILIARY_SLOTS = 10,
    UNUSUALLY_LARGE_PARTICLE_COUNT = 1000000,
    RENDER_TARGET_FPS_MAX = 1000,
    RENDER_STRUCTURE_FACTOR_SAMPLE_CAP_MAX = 8192
};

static const double RENDER_SUPERSAMPLING_MIN = 0.5;
static const double RENDER_SUPERSAMPLING_MAX = 3.0;
static const double RENDER_CAMERA_FOV_MIN = 10.0;
static const double RENDER_CAMERA_FOV_MAX = 120.0;
static const double RENDER_REFRESH_SECONDS_MIN = 1.0/1000.0;
static const double RENDER_REFRESH_SECONDS_MAX = 3600.0;

/* Legacy globals.  This translation unit is intentionally the only bridge
 * between the declarative configuration and EDMD.c's current global state. */
extern int N;
extern double phi, sizeratio, fractionSmallN, aspectRatio;
extern int controlLenght;
extern double Lx, Ly, Lz;
extern double phig, phil;

extern int tempProfileThermo, tempProfileNSlabs;

extern int load, S1, Hex, MembraneInit, hexRadius;
extern char initialFilePath[350];
extern int osmosisThermo;
extern double sigma;
extern int interfaceThermo, interfaceThermoTemp;
extern int velocityDistributionThermo, profileThermo, nProfileBins;
extern int clusterThermo, dumpCluster, killCluster, clusterMode;
extern int clusterNeighborThreshold;
extern double clusterCutoff;
extern int strucThermo, doFQT;
extern double qmax;
extern int bubbleStop, areaThermo, boopThermo;
extern int pcfThermo, pcfBondOrderThermo, pcfg6Thermo;
extern int critical, snapshotCritical;
extern double fracCritical;

extern double tmax, firstScreen, dtime, dtimeThermo, firstThermo;
extern int periodicDumps;
extern int Nthermo, logSpacing;
extern double base, dr_pressure, nextScreen;

extern double delta, deltaM, deltam, ts, tau, deltaMax;
extern double vo, ao, sig, U;
extern double deltaE, beta, taur, additionalEnergy;
extern double field, Einit, resW, resT, res, resBeta, deltaOmega;
extern double gamm, T, dtnoise, strainRate;
extern double a, b;
extern double proportionPositivelyCharged, proportionNeutralyCharged;
extern double vr;
extern int automaticGrowthRate;
extern int version, outputVersion;

#if G
extern int ther, reduce, msd;
extern int addEnergyThermo, addPressureThermo;
extern int thermoPressureField, addCollisionNumberThermo;
extern int addEnergyCenterOfMassThermo, thermoKurt, thermoSonine;
extern int unwrap;
#else
extern int ther, msd;
extern const int reduce;
extern const int addEnergyThermo, addPressureThermo;
extern const int thermoPressureField, addCollisionNumberThermo;
extern const int addEnergyCenterOfMassThermo, thermoKurt, thermoSonine;
extern const int unwrap;
#endif

#if G
extern bool addWell, addField, damping;
extern bool addDelta, addDeltaTangent, addDoubleDelta;
extern bool addDoubleDeltaTangent, polydispersity, thermoWall;
extern bool charged, addShear, addShearWall;
extern int noise, addWally, addWallx, addWallz, addMidWall, addCircularWall;
extern int addEvolvingDelta, addExpo, addEnergy, addActiveDelta;
extern int liquidliquid;
extern int euler;
extern double randomInjection;
extern int simplifyInjection;
extern double shearWallSpeed;
#else
extern const int addWell, addField, damping;
extern const int addDelta, addDeltaTangent, addDoubleDelta;
extern const int addDoubleDeltaTangent, polydispersity, thermoWall;
extern const int charged, addShear, addShearWall;
extern const int noise, addWally, addWallx, addWallz, addMidWall, addCircularWall;
extern const int addEvolvingDelta, addExpo, addEnergy, addActiveDelta;
extern const int liquidliquid;
extern const int euler;
extern const double randomInjection;
extern const int simplifyInjection;
extern const double shearWallSpeed;
#endif

#if INTRUDER
#include "intruder.h"
extern double res_intruder, delta_intruder, bathRad;
extern double Im, M, angle, omega, VX, VY;
extern enum intruderShape shape;
extern double paramIntruder[];
extern const int rethermalize_after_intruder_collision;
extern double displayRad;
extern int km_enabled;
#endif

static bool is_finite_positive(double value)
{
    return isfinite(value) && value > 0.0;
}

static bool init_mode_uses_generic_particle_builder(SimInitMode mode)
{
    return mode == SIM_INIT_RANDOM_GROWTH ||
           mode == SIM_INIT_MEMBRANE_OSMOSIS;
}

static bool init_mode_supports_fixed_box_length(SimInitMode mode)
{
    return mode == SIM_INIT_RANDOM_GROWTH;
}

static bool init_mode_uses_phase_fractions(SimInitMode mode)
{
    return mode == SIM_INIT_MEMBRANE_OSMOSIS;
}

static SimInitMode capture_init_mode(void)
{
    if (load == 2) return SIM_INIT_DUMP_FILE;
    if (S1) return SIM_INIT_SQUARE_LATTICE;
    if (Hex) return SIM_INIT_HEX_LATTICE;
    if (MembraneInit) return SIM_INIT_MEMBRANE_OSMOSIS;
    return SIM_INIT_RANDOM_GROWTH;
}

static SimBoundaryMode capture_boundary_mode(void)
{
    if (addCircularWall) return SIM_BOUNDARY_CIRCULAR_SPHERICAL;
    if (addWallz) return SIM_BOUNDARY_WALL_XYZ;
    if (addWallx && addWally) return SIM_BOUNDARY_WALL_XY;
    if (addWallx) return SIM_BOUNDARY_WALL_X;
    if (addWally) return SIM_BOUNDARY_WALL_Y;
    return SIM_BOUNDARY_PERIODIC;
}

static SimInjectionMode capture_injection_mode(void)
{
#if TANGENTIAL
    return SIM_INJECTION_ROUGH_ROTATIONAL;
#else
    if (addDoubleDeltaTangent) return SIM_INJECTION_DOUBLE_DELTA_TANGENT;
    if (addDeltaTangent) return SIM_INJECTION_DELTA_TANGENT;
    if (addDoubleDelta) return SIM_INJECTION_DOUBLE_DELTA_NORMAL;
    if (addEvolvingDelta) return SIM_INJECTION_EVOLVING_DELTA;
    if (addDelta) return SIM_INJECTION_DELTA_NORMAL;
    if (addEnergy) return SIM_INJECTION_COLLISION_ENERGY;
    if (addActiveDelta) return SIM_INJECTION_ACTIVE_TANGENTIAL;
    return SIM_INJECTION_NONE;
#endif
}

BuildCapabilities simulation_build_capabilities(void)
{
    BuildCapabilities result;
    memset(&result, 0, sizeof(result));
    result.dimensions = THREE_D ? 3 : 2;
    result.graphics = G != 0;
    result.stress = STRESS != 0;
    result.heat = HEAT != 0;
    result.rough_particles = TANGENTIAL != 0;
    result.intruder = INTRUDER != 0;
    result.additive_contacts = NONADDITIVE != 0;
    result.runtime_physics_toggles = G != 0;
    result.runtime_basic_injection_toggles = G != 0;
    result.runtime_advanced_injection_toggles = G != 0;
    result.runtime_exponential_restitution = G != 0;
    result.runtime_charged_potential = G != 0;
    result.runtime_same_species_potential = G != 0;
    result.runtime_output_lifecycle = G != 0;
    result.runtime_output_column_toggles = G != 0;
    result.runtime_intruder_parameters = INTRUDER != 0;

    result.spherical_confinement = true;
    result.tangent_delta_3d = true;
    result.wall_z = THREE_D != 0;
    result.full_3d_thermostats = true;
    result.shear_3d = true;
    /* The HCP/hexagonal initializer is dimension-complete. Other legacy
       lattice layouts remain explicitly planar in the validator. */
    result.lattice_initializers_3d = true;
    result.intruder_3d = false;
    result.pressure_field_3d = false;
    result.projected_xy_observables = THREE_D != 0;
    result.queue_verification = EDMD_VERIFY_QUEUE != 0;
    result.overlap_verification = EDMD_VERIFY_OVERLAPS != 0;
    result.collision_accounting_verification =
        EDMD_VERIFY_COLLISION_ACCOUNTING != 0;
    return result;
}

void simulation_config_set_defaults(SimulationConfig *config,
                                    const BuildCapabilities *capabilities)
{
    BuildCapabilities detected;
    bool graphical;
    int dimensions;
    if (config == NULL) return;
    if (capabilities == NULL) {
        detected = simulation_build_capabilities();
        capabilities = &detected;
    }
    graphical = capabilities->graphics;
    dimensions = capabilities->dimensions;
    memset(config, 0, sizeof(*config));
    config->schema_version = SIM_CONFIG_SCHEMA_VERSION;

    config->initialization.mode = SIM_INIT_RANDOM_GROWTH;
    config->initialization.particle_count = graphical ? 1000 : 500;
    config->initialization.packing_fraction = graphical ? 0.5 : 0.38;
    config->initialization.composition = SIM_COMPOSITION_MONODISPERSE;
    config->initialization.size_ratio = 0.5;
    config->initialization.small_fraction = 0.0;
    config->initialization.polydisperse_a = 2.0;
    config->initialization.polydisperse_b = 5.0;
    config->initialization.initial_energy = 1.0;
    config->initialization.random_seed = 1;
    config->initialization.box_length_mode = SIM_BOX_FROM_PACKING_FRACTION;
    config->initialization.aspect_ratio = 1.0;
    config->initialization.liquid_packing_fraction = 0.2;
    config->initialization.gas_packing_fraction = 0.6;
    config->initialization.automatic_growth_rate = true;

    config->geometry.boundary = SIM_BOUNDARY_PERIODIC;
    config->geometry.wall_normal_restitution = 1.0;
    config->geometry.wall_tangent_restitution = 1.0;
    config->geometry.shear_wall_speed = 0.1;

    config->dynamics.flight = SIM_FLIGHT_BALLISTIC;
    config->dynamics.constant_field = graphical ? -0.001 : -0.02;
    config->dynamics.damping_gamma = graphical ? 0.025 : 0.02;
    config->dynamics.shear_strain_rate = -0.01;
    config->dynamics.restitution_model = SIM_RESTITUTION_CONSTANT;
    config->dynamics.particle_restitution = 1.0;
    config->dynamics.exponential_velocity_scale = 3.5;
    config->dynamics.exponential_amplitude = 1.5;
    config->dynamics.injection = capabilities->rough_particles ?
        SIM_INJECTION_ROUGH_ROTATIONAL : SIM_INJECTION_NONE;
    config->dynamics.delta = graphical ? 1.0 : 0.5;
    config->dynamics.delta_small = 0.1;
    config->dynamics.delta_large = 2.0;
    config->dynamics.collision_memory_time = 10.0;
    config->dynamics.evolving_delta_tau = 5.0;
    config->dynamics.evolving_delta_max = 0.1;
    config->dynamics.injected_energy = 50.0;
    config->dynamics.background_injected_energy = 0.01;
    config->dynamics.injection_recovery_time = 3.0;
    config->dynamics.injection_exponent = 10.0;
    config->dynamics.rough_tangential_delta = false;
    config->dynamics.rough_tangential_restitution = 0.0;
    config->dynamics.rough_angular_kick = 100.0;

    config->driving.thermostat = SIM_THERMOSTAT_NONE;
    config->driving.thermostat_interval = graphical ? 1.0 : 0.5;
    config->driving.temperature = graphical ? 0.01 : 0.05;
    config->driving.gamma = graphical ? 0.025 : 0.02;

    config->potential.energy = graphical ? 0.15 : -0.8;
    config->potential.range = graphical ? 1.5 : 1.6;
    config->potential.species_mode = SIM_POTENTIAL_ALL_PAIRS;
    config->potential.positive_charge_fraction = 0.2;
    config->potential.neutral_charge_fraction = 0.1;

    /* The desktop workspace computes live plots in memory.  Scheduled file
       streams remain a headless concern; manual snapshots are independent. */
    config->observables.thermo_enabled = !graphical;
    config->observables.collision_count = true;
    config->observables.energy = true;
    config->observables.pressure = true;
    config->observables.sonine = true;
    config->observables.pressure_field_spacing = 1.0;
    config->observables.temperature_profile_slabs = 40;
    config->observables.velocity_profile_bins = 100;
    config->observables.osmosis_slice_width = 2.0;
    config->observables.structure_factor_qmax = 0.3;
    config->observables.tvacf_mode_m = 3;
    config->observables.tvacf_mode_n = 0;
    config->observables.tvacf_samples = 101;
    config->observables.tvacf_time_span = 3000.0;
    config->observables.tvacf_checkpoint_interval = 100;
    config->observables.cluster_neighbor_threshold = 3;
    config->observables.cluster_cutoff = 3.55;
    config->observables.critical_fraction = -1.0/6.0;
    /* Interactive runs are controlled by the user.  Batch runs retain the
       historical finite horizon unless --time overrides it. */
    config->output.maximum_time = graphical ? never : 30000.0;
    config->output.periodic_dumps = !graphical;
    config->output.first_dump_time = graphical ? 0.0 : 100.0;
    config->output.dump_interval = graphical ? 1.0 : 100000.0;
    config->output.next_dump_time_override = -1.0;
    config->output.first_thermo_time = 1000.0;
    config->output.thermo_interval = 1000.0;
    config->output.thermo_average_samples = 1;
    config->output.logarithmic_spacing_base = 1.01;
    config->output.version = 1;

    config->rendering.start_running = true;
    config->rendering.target_fps = 144;
    config->rendering.color_mode = SIM_COLOR_BY_TYPE;
    config->rendering.palette = SIM_PALETTE_PLASMA;
    config->rendering.show_fps = true;
    config->rendering.show_axes = dimensions == 3;
    config->rendering.show_boundary = true;
    config->rendering.show_selection_outline = true;
    config->rendering.sphere_quality = SIM_SPHERE_QUALITY_AUTO;
    config->rendering.scene_supersampling = 1.0;
    config->rendering.camera_fov_degrees = 45.0;
    config->rendering.structure_factor_sample_cap = 4096;
    config->rendering.structure_factor_refresh_seconds = 0.5;
    config->rendering.color_refresh_seconds = 1.0/30.0;
    config->rendering.type0_rgba = UINT32_C(0x4296faff);
    config->rendering.type1_rgba = UINT32_C(0xf29d49ff);
    config->rendering.background_rgba = UINT32_C(0x181a1fff);
    config->rendering.boundary_rgba = UINT32_C(0x8b949eff);

    config->intruder.enabled = false;
    config->intruder.shape = 7; /* POLYGON */
    config->intruder.shape_parameters[0] = 14.1421356237;
    config->intruder.shape_parameters[1] = 4.0;
    config->intruder.mass = 1.0e16;
    config->intruder.moment_of_inertia = 1.0e15;
    config->intruder.restitution = 1.0;
    config->intruder.bath_particle_radius = 1.0;
    config->intruder.kinetic_model_measurement = true;
}

void simulation_config_capture_globals(SimulationConfig *config)
{
    BuildCapabilities capabilities = simulation_build_capabilities();
    if (config == NULL) return;
    simulation_config_set_defaults(config, &capabilities);

    config->initialization.mode = capture_init_mode();
    (void)strncpy(config->initialization.load_path, initialFilePath,
                  sizeof(config->initialization.load_path) - 1);
    config->initialization.load_path[
        sizeof(config->initialization.load_path) - 1] = '\0';
    config->initialization.particle_count = N;
    config->initialization.packing_fraction = phi;
    config->initialization.composition = polydispersity ?
        SIM_COMPOSITION_POLYDISPERSE :
        ((config->initialization.mode == SIM_INIT_SQUARE_LATTICE ||
          config->initialization.mode == SIM_INIT_MEMBRANE_OSMOSIS ||
          fractionSmallN > 0.0) ? SIM_COMPOSITION_BINARY :
                                  SIM_COMPOSITION_MONODISPERSE);
    config->initialization.size_ratio = sizeratio;
    config->initialization.small_fraction = fractionSmallN;
    config->initialization.polydisperse_a = a;
    config->initialization.polydisperse_b = b;
    config->initialization.initial_energy = Einit;
    config->initialization.random_seed = version;
    config->initialization.box_length_mode = controlLenght == 2 ?
        SIM_BOX_FIXED_LX : (controlLenght == 1 ? SIM_BOX_FIXED_LY :
                                                   SIM_BOX_FROM_PACKING_FRACTION);
    config->initialization.aspect_ratio = aspectRatio;
    config->initialization.box_lx = Lx;
    config->initialization.box_ly = Ly;
    config->initialization.box_lz = Lz;
    config->initialization.liquid_packing_fraction = phil;
    config->initialization.gas_packing_fraction = phig;
    config->initialization.hex_hole_radius = hexRadius;
    config->initialization.automatic_growth_rate = automaticGrowthRate != 0;
    config->initialization.growth_rate = vr;

    config->geometry.boundary = capture_boundary_mode();
    config->geometry.center_membrane = addMidWall != 0;
    config->geometry.wall_normal_restitution = resW;
    config->geometry.wall_tangent_restitution = resT;
    config->geometry.thermal_wall = thermoWall != 0;
    config->geometry.shear_wall = addShearWall != 0;
    config->geometry.shear_wall_speed = shearWallSpeed;

    config->dynamics.flight = addShear ? SIM_FLIGHT_LEES_EDWARDS_SHEAR :
        (addField ? SIM_FLIGHT_CONSTANT_FIELD :
         (damping ? SIM_FLIGHT_LINEAR_DAMPING : SIM_FLIGHT_BALLISTIC));
    config->dynamics.constant_field = field;
    config->dynamics.damping_gamma = gamm;
    config->dynamics.shear_strain_rate = strainRate;
    config->dynamics.restitution_model = addExpo ?
        SIM_RESTITUTION_EXPONENTIAL : SIM_RESTITUTION_CONSTANT;
    config->dynamics.particle_restitution = res;
    config->dynamics.exponential_velocity_scale = vo;
    config->dynamics.exponential_amplitude = ao;
    config->dynamics.injection = capture_injection_mode();
    config->dynamics.delta = delta;
    config->dynamics.delta_small = deltam;
    config->dynamics.delta_large = deltaM;
    config->dynamics.collision_memory_time = ts;
    config->dynamics.evolving_delta_tau = tau;
    config->dynamics.evolving_delta_max = deltaMax;
    config->dynamics.injected_energy = deltaE;
    config->dynamics.background_injected_energy = additionalEnergy;
    config->dynamics.injection_recovery_time = taur;
    config->dynamics.injection_exponent = beta;
    config->dynamics.randomize_injected_energy = randomInjection != 0.0;
    config->dynamics.simplified_energy_injection = simplifyInjection != 0;
    config->dynamics.rough_tangential_delta = addDeltaTangent != 0;
    config->dynamics.rough_tangential_restitution = resBeta;
    config->dynamics.rough_angular_kick = deltaOmega;

    config->driving.thermostat = (SimThermostatMode)noise;
    config->driving.thermostat_interval = dtnoise;
    config->driving.temperature = T;
    config->driving.gamma = gamm;
    config->driving.euler_noise = euler != 0;

    config->potential.enabled = addWell != 0;
    config->potential.energy = U;
    config->potential.range = sig;
    config->potential.species_mode = charged ? SIM_POTENTIAL_CHARGED :
        (liquidliquid ? SIM_POTENTIAL_SAME_SPECIES :
                       SIM_POTENTIAL_ALL_PAIRS);
    config->potential.positive_charge_fraction = proportionPositivelyCharged;
    config->potential.neutral_charge_fraction = proportionNeutralyCharged;

    config->observables.thermo_enabled = ther != 0;
    config->observables.collision_count = addCollisionNumberThermo != 0;
    config->observables.energy = addEnergyThermo != 0;
    config->observables.pressure = addPressureThermo != 0;
    config->observables.center_of_mass_energy =
        addEnergyCenterOfMassThermo != 0;
    config->observables.kurtosis = thermoKurt != 0;
    config->observables.sonine = thermoSonine != 0;
    config->observables.mean_squared_displacement = msd != 0;
    config->observables.pressure_field = thermoPressureField != 0;
    config->observables.pressure_field_spacing = dr_pressure;
    config->observables.temperature_profile = tempProfileThermo != 0;
    config->observables.temperature_profile_slabs = tempProfileNSlabs;
    config->observables.velocity_profile = profileThermo != 0;
    config->observables.velocity_profile_bins = nProfileBins;
    config->observables.interface_density = interfaceThermo != 0;
    config->observables.interface_temperature = interfaceThermoTemp != 0;
    config->observables.osmosis = osmosisThermo != 0;
    config->observables.osmosis_slice_width = sigma;
    config->observables.structure_factor =
        (strucThermo >= SIM_STRUCTURE_FACTOR_DENSITY &&
         strucThermo <= SIM_STRUCTURE_FACTOR_VELOCITY)
        ? (SimStructureFactorMode)strucThermo : SIM_STRUCTURE_FACTOR_OFF;
    config->observables.intermediate_scattering = doFQT != 0;
    config->observables.structure_factor_qmax = qmax;
    {
        int tvacf_enabled = 0;
        getTransverseVelocityCorrelationConfig(
            &tvacf_enabled, &config->observables.tvacf_mode_m,
            &config->observables.tvacf_mode_n,
            &config->observables.tvacf_samples,
            &config->observables.tvacf_time_span,
            &config->observables.tvacf_checkpoint_interval);
        config->observables.transverse_velocity_correlation =
            tvacf_enabled != 0;
    }
    config->observables.clusters = clusterThermo != 0;
    config->observables.cluster_histogram = dumpCluster != 0;
    config->observables.cluster_stop_size = killCluster;
    config->observables.cluster_mode = clusterMode;
    config->observables.cluster_neighbor_threshold = clusterNeighborThreshold;
    config->observables.cluster_cutoff = clusterCutoff;
    config->observables.bond_order = (SimBoopMode)boopThermo;
    config->observables.voronoi_area = areaThermo != 0;
    config->observables.pair_correlation = pcfThermo != 0;
    config->observables.bond_order_pair_correlation =
        pcfBondOrderThermo != 0;
    config->observables.g6_pair_correlation = pcfg6Thermo != 0;
    config->observables.velocity_distribution =
        velocityDistributionThermo != 0;
    config->observables.critical_regions = critical != 0;
    config->observables.critical_snapshot_labels = snapshotCritical != 0;
    config->observables.critical_fraction = fracCritical;
    config->observables.bubble_stop = bubbleStop != 0;

    config->output.maximum_time = tmax;
    config->output.periodic_dumps = periodicDumps != 0;
    config->output.first_dump_time = firstScreen;
    config->output.dump_interval = dtime;
    config->output.next_dump_time_override = nextScreen;
    config->output.first_thermo_time = firstThermo;
    config->output.thermo_interval = dtimeThermo;
    config->output.thermo_average_samples = Nthermo;
    config->output.logarithmic_thermo_spacing = logSpacing != 0;
    config->output.logarithmic_spacing_base = base;
    config->output.reduced_dump = reduce != 0;
    config->output.unwrap_coordinates = unwrap != 0;
    config->output.version = outputVersion;

#if INTRUDER
    config->intruder.enabled = true;
    config->intruder.shape = (int)shape;
    config->intruder.shape_parameters[0] = paramIntruder[0];
    config->intruder.shape_parameters[1] = paramIntruder[1];
    config->intruder.shape_parameters[2] = paramIntruder[2];
    config->intruder.shape_parameters[3] = paramIntruder[3];
    config->intruder.mass = M;
    config->intruder.moment_of_inertia = Im;
    config->intruder.initial_vx = VX;
    config->intruder.initial_vy = VY;
    config->intruder.initial_angle = angle;
    config->intruder.initial_omega = omega;
    config->intruder.restitution = res_intruder;
    config->intruder.delta = delta_intruder;
    config->intruder.bath_particle_radius = bathRad;
    config->intruder.rethermalize_after_collision =
        rethermalize_after_intruder_collision != 0;
    config->intruder.display_particle_radius = displayRad;
    config->intruder.kinetic_model_measurement = km_enabled != 0;
    config->intruder.kinetic_model_burn_time =
        getIntruderKineticModelBurnTime();
#endif
}

static void apply_initialization_mode(SimInitMode mode)
{
    load = 0;
    S1 = 0;
    Hex = 0;
    MembraneInit = 0;
    switch (mode) {
        case SIM_INIT_DUMP_FILE: load = 2; break;
        case SIM_INIT_SQUARE_LATTICE: S1 = 1; break;
        case SIM_INIT_HEX_LATTICE: Hex = 1; break;
        case SIM_INIT_MEMBRANE_OSMOSIS: MembraneInit = 1; break;
        case SIM_INIT_RANDOM_GROWTH:
        default:
            break;
    }
}

void simulation_config_apply_globals(const SimulationConfig *config)
{
    if (config == NULL) return;
    apply_initialization_mode(config->initialization.mode);
    if (config->initialization.load_path[0] != '\0') {
        (void)strncpy(initialFilePath, config->initialization.load_path, 349);
        initialFilePath[349] = '\0';
    }
    N = config->initialization.particle_count;
    phi = config->initialization.packing_fraction;
    sizeratio = config->initialization.size_ratio;
    fractionSmallN = config->initialization.composition == SIM_COMPOSITION_BINARY ?
        config->initialization.small_fraction : 0.0;
    a = config->initialization.polydisperse_a;
    b = config->initialization.polydisperse_b;
    Einit = config->initialization.initial_energy;
    version = config->initialization.random_seed;
    controlLenght = config->geometry.boundary ==
        SIM_BOUNDARY_CIRCULAR_SPHERICAL ? 0 :
        (config->initialization.box_length_mode == SIM_BOX_FIXED_LX ? 2 :
         (config->initialization.box_length_mode == SIM_BOX_FIXED_LY ? 1 : 0));
    aspectRatio = config->geometry.boundary ==
        SIM_BOUNDARY_CIRCULAR_SPHERICAL ? 1.0 :
        config->initialization.aspect_ratio;
    Lx = config->initialization.box_lx;
    Ly = config->initialization.box_ly;
    Lz = config->initialization.box_lz;
    phil = config->initialization.liquid_packing_fraction;
    phig = config->initialization.gas_packing_fraction;
    hexRadius = (int)config->initialization.hex_hole_radius;
    automaticGrowthRate = config->initialization.automatic_growth_rate;
    if (!automaticGrowthRate) vr = config->initialization.growth_rate;

    resW = config->geometry.wall_normal_restitution;
    resT = config->geometry.wall_tangent_restitution;
    field = config->dynamics.constant_field;
    gamm = config->dynamics.damping_gamma;
    strainRate = config->dynamics.shear_strain_rate;
    res = config->dynamics.particle_restitution;
    vo = config->dynamics.exponential_velocity_scale;
    ao = config->dynamics.exponential_amplitude;
    delta = config->dynamics.delta;
    deltam = config->dynamics.delta_small;
    deltaM = config->dynamics.delta_large;
    ts = config->dynamics.collision_memory_time;
    tau = config->dynamics.evolving_delta_tau;
    deltaMax = config->dynamics.evolving_delta_max;
    deltaE = config->dynamics.injected_energy;
    additionalEnergy = config->dynamics.background_injected_energy;
    taur = config->dynamics.injection_recovery_time;
    beta = config->dynamics.injection_exponent;
    resBeta = config->dynamics.rough_tangential_restitution;
    deltaOmega = config->dynamics.rough_angular_kick;

    dtnoise = config->driving.thermostat_interval;
    T = config->driving.temperature;
    /* Gamma is shared by the legacy damping and thermostat implementations. */
    gamm = config->dynamics.flight == SIM_FLIGHT_LINEAR_DAMPING ?
        config->dynamics.damping_gamma : config->driving.gamma;

    U = config->potential.energy;
    sig = config->potential.range;
    proportionPositivelyCharged = config->potential.positive_charge_fraction;
    proportionNeutralyCharged = config->potential.neutral_charge_fraction;

    tempProfileThermo = config->observables.temperature_profile;
    tempProfileNSlabs = config->observables.temperature_profile_slabs;
    profileThermo = config->observables.velocity_profile;
    nProfileBins = config->observables.velocity_profile_bins;
    interfaceThermo = config->observables.interface_density;
    interfaceThermoTemp = config->observables.interface_temperature;
    osmosisThermo = config->observables.osmosis;
    sigma = config->observables.osmosis_slice_width;
    strucThermo = (int)config->observables.structure_factor;
    doFQT = config->observables.intermediate_scattering;
    qmax = config->observables.structure_factor_qmax;
    setTransverseVelocityCorrelationConfig(
        config->observables.transverse_velocity_correlation,
        config->observables.tvacf_mode_m,
        config->observables.tvacf_mode_n,
        config->observables.tvacf_samples,
        config->observables.tvacf_time_span,
        config->observables.tvacf_checkpoint_interval);
    clusterThermo = config->observables.clusters;
    dumpCluster = config->observables.cluster_histogram;
    killCluster = config->observables.cluster_stop_size;
    clusterMode = config->observables.cluster_mode;
    clusterNeighborThreshold = config->observables.cluster_neighbor_threshold;
    clusterCutoff = config->observables.cluster_cutoff;
    boopThermo = (int)config->observables.bond_order;
    areaThermo = config->observables.voronoi_area;
    pcfThermo = config->observables.pair_correlation;
    pcfBondOrderThermo = config->observables.bond_order_pair_correlation;
    pcfg6Thermo = config->observables.g6_pair_correlation;
    velocityDistributionThermo = config->observables.velocity_distribution;
    critical = config->observables.critical_regions;
    snapshotCritical = config->observables.critical_snapshot_labels;
    fracCritical = config->observables.critical_fraction;
    bubbleStop = config->observables.bubble_stop;

    tmax = config->output.maximum_time;
    periodicDumps = config->output.periodic_dumps;
    firstScreen = config->output.first_dump_time;
    dtime = config->output.dump_interval;
    nextScreen = config->output.next_dump_time_override;
    firstThermo = config->output.first_thermo_time;
    dtimeThermo = config->output.thermo_interval;
    Nthermo = config->output.thermo_average_samples;
    logSpacing = config->output.logarithmic_thermo_spacing;
    base = config->output.logarithmic_spacing_base;
    outputVersion = config->output.version;
    dr_pressure = config->observables.pressure_field_spacing;

#if G
    ther = config->observables.thermo_enabled;
    addCollisionNumberThermo = config->observables.collision_count;
    addEnergyThermo = config->observables.energy;
    addPressureThermo = config->observables.pressure;
    addEnergyCenterOfMassThermo = config->observables.center_of_mass_energy;
    thermoKurt = config->observables.kurtosis;
    thermoSonine = config->observables.sonine;
    msd = config->observables.mean_squared_displacement;
    thermoPressureField = config->observables.pressure_field;
    reduce = config->output.reduced_dump;
    unwrap = config->output.unwrap_coordinates;
#endif

#if G
    polydispersity = config->initialization.composition ==
        SIM_COMPOSITION_POLYDISPERSE;
    addWallx = config->geometry.boundary == SIM_BOUNDARY_WALL_X ||
        config->geometry.boundary == SIM_BOUNDARY_WALL_XY ||
        config->geometry.boundary == SIM_BOUNDARY_WALL_XYZ;
    addWally = config->geometry.boundary == SIM_BOUNDARY_WALL_Y ||
        config->geometry.boundary == SIM_BOUNDARY_WALL_XY ||
        config->geometry.boundary == SIM_BOUNDARY_WALL_XYZ;
    addWallz = config->geometry.boundary == SIM_BOUNDARY_WALL_XYZ;
    addCircularWall = config->geometry.boundary ==
        SIM_BOUNDARY_CIRCULAR_SPHERICAL;
    addMidWall = config->geometry.center_membrane;
    thermoWall = config->geometry.thermal_wall;
    addShearWall = config->geometry.shear_wall;
    shearWallSpeed = config->geometry.shear_wall_speed;

    addField = config->dynamics.flight == SIM_FLIGHT_CONSTANT_FIELD;
    damping = config->dynamics.flight == SIM_FLIGHT_LINEAR_DAMPING;
    addShear = config->dynamics.flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR;
    addExpo = config->dynamics.restitution_model ==
        SIM_RESTITUTION_EXPONENTIAL;
    #if TANGENTIAL
    /* Rough contacts are fixed by the build.  Their optional translational
       tangential kick is an independent, previously hidden runtime switch. */
    addDelta = false;
    addDeltaTangent = config->dynamics.rough_tangential_delta;
    addDoubleDelta = false;
    addDoubleDeltaTangent = false;
    addEvolvingDelta = false;
    addEnergy = false;
    addActiveDelta = false;
    #else
    addDelta = config->dynamics.injection == SIM_INJECTION_DELTA_NORMAL;
    addDeltaTangent = config->dynamics.injection == SIM_INJECTION_DELTA_TANGENT;
    addDoubleDelta = config->dynamics.injection ==
        SIM_INJECTION_DOUBLE_DELTA_NORMAL;
    addDoubleDeltaTangent = config->dynamics.injection ==
        SIM_INJECTION_DOUBLE_DELTA_TANGENT;
    addEvolvingDelta = config->dynamics.injection ==
        SIM_INJECTION_EVOLVING_DELTA;
    addEnergy = config->dynamics.injection ==
        SIM_INJECTION_COLLISION_ENERGY;
    addActiveDelta = config->dynamics.injection ==
        SIM_INJECTION_ACTIVE_TANGENTIAL;
    #endif
    randomInjection = config->dynamics.randomize_injected_energy ? 1.0 : 0.0;
    simplifyInjection = config->dynamics.simplified_energy_injection;
    noise = (int)config->driving.thermostat;
    euler = config->driving.euler_noise;
    addWell = config->potential.enabled;
    /* Species rules select which square-well events are constructed.  Keep a
       remembered GUI choice inert while the potential itself is disabled. */
    charged = config->potential.enabled &&
        config->potential.species_mode == SIM_POTENTIAL_CHARGED;
    liquidliquid = config->potential.enabled &&
        config->potential.species_mode == SIM_POTENTIAL_SAME_SPECIES;
#endif

#if INTRUDER
    M = config->intruder.mass;
    Im = config->intruder.moment_of_inertia;
    VX = config->intruder.initial_vx;
    VY = config->intruder.initial_vy;
    angle = config->intruder.initial_angle;
    omega = config->intruder.initial_omega;
    res_intruder = config->intruder.restitution;
    delta_intruder = config->intruder.delta;
    bathRad = config->intruder.bath_particle_radius;
    displayRad = config->intruder.display_particle_radius;
    km_enabled = config->intruder.kinetic_model_measurement;
    setIntruderKineticModelBurnTime(config->intruder.kinetic_model_burn_time);
    shape = (enum intruderShape)config->intruder.shape;
    paramIntruder[0] = config->intruder.shape_parameters[0];
    paramIntruder[1] = config->intruder.shape_parameters[1];
    paramIntruder[2] = config->intruder.shape_parameters[2];
    paramIntruder[3] = config->intruder.shape_parameters[3];
#endif
}

typedef struct IssueWriter {
    ConfigIssue *issues;
    size_t capacity;
    size_t count;
} IssueWriter;

static void add_issue(IssueWriter *writer, ConfigIssueSeverity severity,
                      const char *field, const char *format, ...)
{
    if (writer->count < writer->capacity && writer->issues != NULL) {
        ConfigIssue *issue = writer->issues + writer->count;
        va_list arguments;
        issue->severity = severity;
        snprintf(issue->field, sizeof(issue->field), "%s",
                 field != NULL ? field : "");
        va_start(arguments, format);
        vsnprintf(issue->message, sizeof(issue->message), format, arguments);
        va_end(arguments);
    }
    writer->count++;
}

static bool thermostat_uses_temperature(SimThermostatMode mode)
{
    return mode == SIM_THERMOSTAT_LANGEVIN ||
           mode == SIM_THERMOSTAT_VELOCITY_RESCALE ||
           mode == SIM_THERMOSTAT_MAXWELL_REDRAW;
}

static bool thermostat_uses_gamma(SimThermostatMode mode)
{
    return mode == SIM_THERMOSTAT_LANGEVIN;
}

static bool thermostat_uses_euler_choice(SimThermostatMode mode)
{
    return mode == SIM_THERMOSTAT_LANGEVIN;
}

static bool init_mode_is_planar(SimInitMode mode)
{
    /* Dump restarts carry z and vz columns and particlesInit() restores them
       in a THREE_D build. The HCP constructor is also dimension-complete;
       the remaining explicitly planar constructors stay disabled. */
    return mode == SIM_INIT_SQUARE_LATTICE ||
           mode == SIM_INIT_MEMBRANE_OSMOSIS;
}

static bool injection_supplies_energy(SimInjectionMode mode)
{
    return mode != SIM_INJECTION_NONE &&
           mode != SIM_INJECTION_ROUGH_ROTATIONAL;
}

size_t simulation_config_validate(const SimulationConfig *config,
                                  const BuildCapabilities *capabilities,
                                  ConfigIssue *issues,
                                  size_t issue_capacity)
{
    BuildCapabilities detected;
    IssueWriter writer = { issues, issue_capacity, 0 };
    const SimInitializationConfig *initialization;
    const SimGeometryConfig *geometry;
    const SimDynamicsConfig *dynamics;
    const SimDrivingConfig *driving;
    const SimPotentialConfig *potential;
    const SimObservablesConfig *observables;
    const SimOutputConfig *output;
    const SimRenderingConfig *rendering;
    const SimIntruderConfig *intruder;
    bool init_mode_valid;
    bool file_initializer;
    bool generic_composition_initializer;
    bool packing_fraction_active;
    bool size_ratio_active;
    bool small_fraction_active;
    bool polydispersity_active;
    bool aspect_ratio_active;
    bool particle_count_valid;
    bool aspect_ratio_valid;
    bool wall_enabled;
    bool wall_normal_restitution_active;
    bool wall_tangent_restitution_active;
    bool collision_injection_active;
    bool dissipative;
    bool driven;

    if (config == NULL) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "config",
                  "The simulation configuration is null.");
        return writer.count;
    }
    if (capabilities == NULL) {
        detected = simulation_build_capabilities();
        capabilities = &detected;
    }
    initialization = &config->initialization;
    geometry = &config->geometry;
    dynamics = &config->dynamics;
    driving = &config->driving;
    potential = &config->potential;
    observables = &config->observables;
    output = &config->output;
    rendering = &config->rendering;
    intruder = &config->intruder;

    if (config->schema_version != SIM_CONFIG_SCHEMA_VERSION) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "schema_version",
                  "Unsupported configuration schema %d (expected %d).",
                  config->schema_version, SIM_CONFIG_SCHEMA_VERSION);
    }
    if (capabilities->dimensions != 2 && capabilities->dimensions != 3) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "capabilities.dimensions",
                  "The build must have either two or three dimensions.");
    }

    init_mode_valid = initialization->mode >= SIM_INIT_RANDOM_GROWTH &&
        initialization->mode <= SIM_INIT_MEMBRANE_OSMOSIS;
    file_initializer = initialization->mode == SIM_INIT_DUMP_FILE;
    generic_composition_initializer =
        initialization->mode == SIM_INIT_RANDOM_GROWTH;
    packing_fraction_active =
        initialization->mode == SIM_INIT_RANDOM_GROWTH ||
        initialization->mode == SIM_INIT_SQUARE_LATTICE ||
        initialization->mode == SIM_INIT_HEX_LATTICE;
    size_ratio_active =
        initialization->mode == SIM_INIT_SQUARE_LATTICE ||
        initialization->mode == SIM_INIT_MEMBRANE_OSMOSIS ||
        (generic_composition_initializer &&
         initialization->composition == SIM_COMPOSITION_BINARY);
    small_fraction_active = generic_composition_initializer &&
        initialization->composition == SIM_COMPOSITION_BINARY;
    polydispersity_active = generic_composition_initializer &&
        initialization->composition == SIM_COMPOSITION_POLYDISPERSE;
    aspect_ratio_active = !file_initializer &&
        geometry->boundary != SIM_BOUNDARY_CIRCULAR_SPHERICAL;
    if (!init_mode_valid) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "initialization.mode",
                  "Unknown initialization mode.");
    }
    if (initialization->mode == SIM_INIT_DUMP_FILE &&
        initialization->load_path[0] == '\0') {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "initialization.load_path",
                  "A configuration path is required for file initialization.");
    }
    if (capabilities->dimensions == 3 &&
        init_mode_is_planar(initialization->mode)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "initialization.mode",
                  "This initializer does not create complete 3D positions and velocities.");
    }
    particle_count_valid = !file_initializer &&
        initialization->particle_count > 0 &&
        initialization->particle_count <=
            (INT_MAX - EVENT_LIST_AUXILIARY_SLOTS) / 2;
    if (!file_initializer && initialization->particle_count <= 0) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.particle_count",
                  "Particle count must be greater than zero.");
    } else if (!file_initializer && initialization->particle_count >
               (INT_MAX - EVENT_LIST_AUXILIARY_SLOTS) / 2) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.particle_count",
                  "Particle count is too large for the engine's signed 2*N event-array indexing.");
    } else if (!file_initializer && initialization->particle_count >
               UNUSUALLY_LARGE_PARTICLE_COUNT) {
        add_issue(&writer, CONFIG_ISSUE_WARNING,
                  "initialization.particle_count",
                  "This unusually large particle count may require many gigabytes of memory.");
    }
    if (packing_fraction_active &&
        !is_finite_positive(initialization->packing_fraction)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.packing_fraction",
                  "Packing fraction must be finite and greater than zero.");
    } else if (packing_fraction_active &&
               ((capabilities->dimensions == 2 &&
                initialization->packing_fraction > 0.91) ||
               (capabilities->dimensions == 3 &&
                initialization->packing_fraction > 0.76))) {
        add_issue(&writer, CONFIG_ISSUE_WARNING,
                  "initialization.packing_fraction",
                  "This packing fraction is above the usual close-packing range and initialization may fail.");
    }
    if (initialization->mode == SIM_INIT_HEX_LATTICE &&
        is_finite_positive(initialization->packing_fraction)) {
        const double close_packing = capabilities->dimensions == 3
            ? M_PI/(3.0*sqrt(2.0))
            : M_PI/(2.0*sqrt(3.0));
        if (initialization->packing_fraction > close_packing) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "initialization.packing_fraction",
                      "The requested lattice packing fraction exceeds exact %s close packing (%.12g) and would start with overlapping particles.",
                      capabilities->dimensions == 3 ? "HCP" : "triangular",
                      close_packing);
        }
    }
    if (!file_initializer &&
        (initialization->composition < SIM_COMPOSITION_MONODISPERSE ||
         initialization->composition > SIM_COMPOSITION_POLYDISPERSE)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "initialization.composition",
                  "Unknown composition mode.");
    }
    if (size_ratio_active &&
        (!is_finite_positive(initialization->size_ratio) ||
         initialization->size_ratio > 1.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "initialization.size_ratio",
                  "Small/large radius ratio must lie in (0, 1]; the cell-list stencil assumes the large-particle radius is one.");
    }
    if (small_fraction_active &&
        (!isfinite(initialization->small_fraction) ||
         initialization->small_fraction <= 0.0 ||
         initialization->small_fraction >= 1.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.small_fraction",
                  "A binary mixture requires both species: the small-particle fraction must lie in (0, 1).");
    }
    if (generic_composition_initializer &&
        initialization->composition == SIM_COMPOSITION_MONODISPERSE &&
        initialization->small_fraction != 0.0) {
        add_issue(&writer, CONFIG_ISSUE_WARNING,
                  "initialization.small_fraction",
                  "The small-particle fraction is ignored for a monodisperse composition.");
    }
    if (polydispersity_active) {
        if (!is_finite_positive(initialization->polydisperse_a)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "initialization.polydisperse_a",
                      "Polydispersity parameter a must be greater than zero.");
        }
        if (!is_finite_positive(initialization->polydisperse_b)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "initialization.polydisperse_b",
                      "Polydispersity parameter b must be greater than zero.");
        }
        if (initialization->small_fraction != 0.0) {
            add_issue(&writer, CONFIG_ISSUE_WARNING,
                      "initialization.small_fraction",
                      "The binary small fraction is ignored for a polydisperse composition.");
        }
    }
    if (init_mode_valid) {
        if ((initialization->mode == SIM_INIT_SQUARE_LATTICE ||
             initialization->mode == SIM_INIT_MEMBRANE_OSMOSIS) &&
            initialization->composition != SIM_COMPOSITION_BINARY) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "initialization.composition",
                      "This initializer constructs a binary mixture and requires binary composition.");
        }
        if (initialization->mode == SIM_INIT_HEX_LATTICE &&
            initialization->composition != SIM_COMPOSITION_MONODISPERSE) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "initialization.composition",
                      "This lattice initializer constructs equal-radius particles and requires monodisperse composition.");
        }
        if (initialization->mode == SIM_INIT_SQUARE_LATTICE) {
            add_issue(&writer, CONFIG_ISSUE_INFO,
                      "initialization.small_fraction",
                      "The square initializer fixes the lattice mixture close to 50:50; the requested small-particle fraction is ignored.");
        }
        if (initialization->mode == SIM_INIT_MEMBRANE_OSMOSIS) {
            add_issue(&writer, CONFIG_ISSUE_INFO,
                      "initialization.small_fraction",
                      "Membrane initialization derives the species fraction from phil, phig, and the size ratio.");
        }
        if (initialization->mode == SIM_INIT_DUMP_FILE) {
            add_issue(&writer, CONFIG_ISSUE_INFO,
                      "initialization.composition",
                      "Particle count, box size, radii, masses, and types are read from the selected file.");
        }
        if (initialization->mode == SIM_INIT_DUMP_FILE) {
            add_issue(&writer, CONFIG_ISSUE_INFO,
                      "initialization.initial_energy",
                      "Dump initialization restores velocities; initial energy is ignored.");
        }
    }
    if (initialization->mode != SIM_INIT_DUMP_FILE &&
        !is_finite_positive(initialization->initial_energy)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.initial_energy",
                  "Initial kinetic energy per particle must be greater than zero.");
    }
    aspect_ratio_valid = aspect_ratio_active &&
        is_finite_positive(initialization->aspect_ratio);
    if (aspect_ratio_active && !aspect_ratio_valid) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.aspect_ratio",
                  "Aspect ratio must be greater than zero.");
    }
    if (initialization->box_length_mode < SIM_BOX_FROM_PACKING_FRACTION ||
        initialization->box_length_mode > SIM_BOX_FIXED_LY) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.box_length_mode",
                  "Unknown box-size control mode.");
    }
    if (capabilities->dimensions == 2 &&
        geometry->boundary != SIM_BOUNDARY_CIRCULAR_SPHERICAL &&
        init_mode_supports_fixed_box_length(initialization->mode) &&
        initialization->box_length_mode == SIM_BOX_FIXED_LX &&
        !is_finite_positive(initialization->box_lx)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "initialization.box_lx",
                  "Fixed Lx must be greater than zero.");
    }
    if (capabilities->dimensions == 2 &&
        geometry->boundary != SIM_BOUNDARY_CIRCULAR_SPHERICAL &&
        init_mode_supports_fixed_box_length(initialization->mode) &&
        initialization->box_length_mode == SIM_BOX_FIXED_LY &&
        !is_finite_positive(initialization->box_ly)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "initialization.box_ly",
                  "Fixed Ly must be greater than zero.");
    }
    if (capabilities->dimensions == 3 &&
        initialization->box_length_mode != SIM_BOX_FROM_PACKING_FRACTION) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.box_length_mode",
                  "Legacy fixed-length sizing uses 2D packing formulas and is unavailable in 3D.");
    }
    if (init_mode_valid &&
        initialization->box_length_mode != SIM_BOX_FROM_PACKING_FRACTION &&
        !init_mode_supports_fixed_box_length(initialization->mode)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.box_length_mode",
                  "The selected initializer computes or loads its own box dimensions; fixed Lx/Ly is unavailable.");
    }
    if (init_mode_valid &&
        init_mode_uses_phase_fractions(initialization->mode)) {
        if (!is_finite_positive(initialization->liquid_packing_fraction) ||
            !is_finite_positive(initialization->gas_packing_fraction)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "initialization.phase_fractions",
                      "Both membrane-region packing fractions must be finite and greater than zero.");
        }
    }
    if (init_mode_valid && initialization->mode == SIM_INIT_HEX_LATTICE &&
        (!isfinite(initialization->hex_hole_radius) ||
         initialization->hex_hole_radius < 0.0 ||
         floor(initialization->hex_hole_radius) !=
             initialization->hex_hole_radius ||
         initialization->hex_hole_radius > INT_MAX)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.hex_hole_radius",
                  "The legacy hexagonal-hole radius must be a nonnegative integer representable by the engine.");
    }
    if (init_mode_valid &&
        init_mode_uses_generic_particle_builder(initialization->mode) &&
        !initialization->automatic_growth_rate &&
        !is_finite_positive(initialization->growth_rate)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.growth_rate",
                  "A manual particle growth rate must be greater than zero.");
    }

    /* The legacy lattice builders truncate floating dimensions to int and
       immediately divide by them.  Check those derived dimensions before a
       restart can reach the unsafe divisions. */
    if (init_mode_valid && particle_count_valid && aspect_ratio_valid &&
        initialization->mode == SIM_INIT_SQUARE_LATTICE) {
        const int half_count = initialization->particle_count/2;
        const double nx_value = sqrt(half_count/initialization->aspect_ratio);
        const double ny_value = sqrt(half_count*initialization->aspect_ratio);
        if (half_count <= 0 || !isfinite(nx_value) || !isfinite(ny_value) ||
            nx_value < 1.0 || ny_value < 1.0 ||
            nx_value > INT_MAX || ny_value > INT_MAX) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "initialization.particle_count",
                      "Particle count and aspect ratio produce a zero or unrepresentable square-lattice dimension.");
        }
    }
    if (init_mode_valid && particle_count_valid && aspect_ratio_valid &&
        initialization->mode == SIM_INIT_HEX_LATTICE) {
        int nx = 0;
        int ny = 0;
        int nz = 1;
        if (!hexagonalLatticeDimensions(initialization->particle_count,
                                        initialization->aspect_ratio,
                                        &nx, &ny, &nz)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "initialization.particle_count",
                      "Particle count and aspect ratio do not produce a complete representable periodic hexagonal lattice.");
        } else {
            const long long lattice_count = (long long)nx*ny*nz;
            if (lattice_count <= 0 || lattice_count >
                    (INT_MAX - EVENT_LIST_AUXILIARY_SLOTS)/2) {
                add_issue(&writer, CONFIG_ISSUE_ERROR,
                          "initialization.particle_count",
                          "The nearest complete periodic hexagonal lattice is too large for the engine's event indexing.");
            } else if (is_finite_positive(initialization->packing_fraction)) {
                double lx = 0.0;
                double ly = 0.0;
                double lz = 0.0;
#if THREE_D
                const double spacing = cbrt((4.0/3.0)*M_PI/
                    (initialization->packing_fraction*(sqrt(2.0)/2.0)));
                lx = nx*spacing;
                ly = ny*(sqrt(3.0)/2.0)*spacing;
                lz = nz*sqrt(2.0/3.0)*spacing;
#else
                ly = sqrt((sqrt(3.0)/2.0)*
                    ((double)ny/nx)*M_PI*lattice_count/
                    initialization->packing_fraction);
                lx = (2.0/sqrt(3.0))*((double)nx/ny)*ly;
#endif
                if (!is_finite_positive(lx) || !is_finite_positive(ly)
#if THREE_D
                    || !is_finite_positive(lz)
#endif
                    ) {
                    add_issue(&writer, CONFIG_ISSUE_ERROR,
                              "initialization.packing_fraction",
                              "Hex-lattice parameters produce non-finite box dimensions.");
                } else if (isfinite(initialization->hex_hole_radius) &&
                           initialization->hex_hole_radius > 0.0 &&
                           floor(initialization->hex_hole_radius) ==
                               initialization->hex_hole_radius) {
                    const double dx = lx/nx;
                    const double dy = ly/ny;
#if THREE_D
                    const double dz = lz/nz;
#endif
                    const int candidate_x[] = {0, 1, nx - 2, nx - 1};
                    const int candidate_y[] = {0, 1, ny - 2, ny - 1};
                    const int candidate_z[] = {0, 1, nz - 2, nz - 1};
                    double farthest_squared = 0.0;
                    for (size_t z_choice = 0;
                         z_choice < sizeof(candidate_z)/sizeof(candidate_z[0]);
                         ++z_choice) {
                        const int z_index = candidate_z[z_choice];
                        if (z_index < 0 || z_index >= nz) continue;
#if THREE_D
                        const int layer = z_index & 1;
#endif
                        for (size_t y_choice = 0;
                             y_choice < sizeof(candidate_y)/sizeof(candidate_y[0]);
                             ++y_choice) {
                            const int y_index = candidate_y[y_choice];
                            if (y_index < 0 || y_index >= ny) continue;
                            for (size_t x_choice = 0;
                                 x_choice < sizeof(candidate_x)/sizeof(candidate_x[0]);
                                 ++x_choice) {
                                const int x_index = candidate_x[x_choice];
                                if (x_index < 0 || x_index >= nx) continue;
#if THREE_D
                                const double x = fmod((x_index +
                                    0.5*(y_index & 1) + 0.5*layer)*dx, lx);
                                const double y = fmod(y_index*dy +
                                    layer*dy/3.0, ly);
                                const double z = z_index*dz;
#else
                                const double x = x_index*dx +
                                    (y_index & 1)*(dx/2.0);
                                const double y = y_index*dy;
                                const double z = 0.0;
#endif
                                const double rx = x - lx/2.0;
                                const double ry = y - ly/2.0;
                                const double rz = z - lz/2.0;
                                const double distance_squared =
                                    rx*rx + ry*ry + rz*rz;
                                if (distance_squared > farthest_squared)
                                    farthest_squared = distance_squared;
                            }
                        }
                    }
                    if (initialization->hex_hole_radius*
                        initialization->hex_hole_radius > farthest_squared) {
                        add_issue(&writer, CONFIG_ISSUE_ERROR,
                                  "initialization.hex_hole_radius",
                                  "The central hole removes every site from the derived hexagonal lattice.");
                    }
                }
            }
        }
    }
    if (geometry->boundary < SIM_BOUNDARY_PERIODIC ||
        geometry->boundary > SIM_BOUNDARY_WALL_XYZ) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.boundary",
                  "Unknown outer-boundary mode.");
    }
    wall_enabled = geometry->boundary != SIM_BOUNDARY_PERIODIC ||
        geometry->center_membrane ||
        initialization->mode == SIM_INIT_MEMBRANE_OSMOSIS;
    wall_normal_restitution_active = wall_enabled &&
        !geometry->thermal_wall;
    wall_tangent_restitution_active = !geometry->thermal_wall &&
        (geometry->boundary == SIM_BOUNDARY_WALL_X ||
         geometry->boundary == SIM_BOUNDARY_WALL_XY ||
         geometry->boundary == SIM_BOUNDARY_WALL_XYZ ||
         geometry->boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL ||
         geometry->center_membrane ||
         initialization->mode == SIM_INIT_MEMBRANE_OSMOSIS ||
         ((geometry->boundary == SIM_BOUNDARY_WALL_Y ||
           geometry->boundary == SIM_BOUNDARY_WALL_XY ||
           geometry->boundary == SIM_BOUNDARY_WALL_XYZ) &&
          (capabilities->dimensions == 3 || !geometry->shear_wall)));
    if (geometry->boundary == SIM_BOUNDARY_WALL_XYZ &&
        (capabilities->dimensions != 3 || !capabilities->wall_z)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.boundary",
                  "The fully closed x/y/z box requires a 3D build with z-wall support.");
    }
    if (geometry->boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL &&
        !capabilities->spherical_confinement) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.boundary",
                  "Circular/spherical confinement is not available in this build.");
    }
    if (geometry->boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL &&
        (geometry->center_membrane ||
         initialization->mode == SIM_INIT_MEMBRANE_OSMOSIS)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.center_membrane",
                  "Circular/spherical confinement cannot be combined with the planar center membrane.");
    }
    if (geometry->boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL &&
        init_mode_valid &&
        initialization->mode != SIM_INIT_RANDOM_GROWTH) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "initialization.mode",
                  "Circular/spherical confinement currently requires the generic random-growth initializer.");
    }
    if (geometry->boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL &&
        initialization->box_length_mode !=
            SIM_BOX_FROM_PACKING_FRACTION) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "initialization.box_length_mode",
                  "Circular/spherical confinement derives its diameter from packing fraction; fixed box lengths are unavailable.");
    }
    /* The center separator is runtime geometry, not ownership of an initial
     * state.  Keeping it independent lets a running membrane be removed (or
     * safely added after an extent preflight) without regenerating particles.
     * Edge-anchored lattice recipes remain incompatible with enabling it at
     * t=0 because their unit-radius sites can intersect x=Lx/2. */
    if (init_mode_valid && geometry->center_membrane &&
        (initialization->mode == SIM_INIT_SQUARE_LATTICE ||
         initialization->mode == SIM_INIT_HEX_LATTICE ||
         initialization->mode == SIM_INIT_DUMP_FILE)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.center_membrane",
                  "This initial-state source cannot safely enable the center membrane without a runtime particle-extent preflight.");
    }
    if (init_mode_valid &&
        (geometry->boundary == SIM_BOUNDARY_WALL_X ||
         geometry->boundary == SIM_BOUNDARY_WALL_Y ||
         geometry->boundary == SIM_BOUNDARY_WALL_XY ||
         geometry->boundary == SIM_BOUNDARY_WALL_XYZ) &&
        (initialization->mode == SIM_INIT_SQUARE_LATTICE ||
         initialization->mode == SIM_INIT_HEX_LATTICE)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.boundary",
                  "This lattice places radius-one particle centers on x=0 or y=0 and cannot start against planar outer walls.");
    }
    if (capabilities->dimensions == 3 &&
        (geometry->boundary == SIM_BOUNDARY_WALL_X ||
         geometry->boundary == SIM_BOUNDARY_WALL_Y ||
         geometry->boundary == SIM_BOUNDARY_WALL_XY)) {
        add_issue(&writer, CONFIG_ISSUE_INFO, "geometry.boundary",
                  "Planar X/Y walls leave the Z direction periodic; this is not a fully closed 3D box.");
    }
    if (geometry->center_membrane && capabilities->dimensions == 3) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.center_membrane",
                  "The membrane initialization and pressure diagnostics are not complete in 3D.");
    }
    if (wall_normal_restitution_active &&
        (!isfinite(geometry->wall_normal_restitution) ||
         geometry->wall_normal_restitution < 0.0 ||
         geometry->wall_normal_restitution > 1.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "geometry.wall_normal_restitution",
                  "Normal wall restitution must lie in [0, 1].");
    }
    if (wall_tangent_restitution_active &&
        (!isfinite(geometry->wall_tangent_restitution) ||
         geometry->wall_tangent_restitution < -1.0 ||
         geometry->wall_tangent_restitution > 1.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "geometry.wall_tangent_restitution",
                  "Tangential wall restitution must lie in [-1, 1].");
    }
    if (geometry->thermal_wall && !wall_enabled) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.thermal_wall",
                  "A thermal wall requires a non-periodic boundary.");
    }
    if (geometry->shear_wall &&
        geometry->boundary != SIM_BOUNDARY_WALL_Y &&
        geometry->boundary != SIM_BOUNDARY_WALL_XY &&
        geometry->boundary != SIM_BOUNDARY_WALL_XYZ) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "geometry.shear_wall",
                  "A shear wall requires Y walls.");
    }
    if (geometry->shear_wall &&
        !isfinite(geometry->shear_wall_speed)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "geometry.shear_wall_speed",
                  "Shear-wall speed must be finite.");
    }

    if (dynamics->flight < SIM_FLIGHT_BALLISTIC ||
        dynamics->flight > SIM_FLIGHT_LEES_EDWARDS_SHEAR) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.flight",
                  "Unknown free-flight model.");
    }
    if ((capabilities->stress || capabilities->heat) &&
        dynamics->flight != SIM_FLIGHT_BALLISTIC &&
        dynamics->flight != SIM_FLIGHT_LEES_EDWARDS_SHEAR) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.flight",
                  "STRESS/HEAT transport bookkeeping supports ballistic and Lees-Edwards free flight.");
    }
    if (dynamics->flight == SIM_FLIGHT_CONSTANT_FIELD &&
        (!isfinite(dynamics->constant_field) ||
         dynamics->constant_field >= 0.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "dynamics.constant_field",
                  "The current constant-field solver requires a strictly negative Y acceleration.");
    }
    if (dynamics->flight == SIM_FLIGHT_LINEAR_DAMPING &&
        !is_finite_positive(dynamics->damping_gamma)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "dynamics.damping_gamma",
                  "Linear damping gamma must be greater than zero.");
    }
    if (dynamics->flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR &&
        !isfinite(dynamics->shear_strain_rate)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "dynamics.shear_strain_rate",
                  "The active shear strain rate must be finite.");
    }
    if (dynamics->flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR) {
        if (capabilities->dimensions == 3 && !capabilities->shear_3d) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.flight",
                      "Lees-Edwards shear accounting is not complete in 3D.");
        }
        if (wall_enabled || geometry->center_membrane || potential->enabled) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.flight",
                      "Lees-Edwards shear requires periodic boundaries and no square potential.");
        }
    }
    if (dynamics->restitution_model < SIM_RESTITUTION_CONSTANT ||
        dynamics->restitution_model > SIM_RESTITUTION_EXPONENTIAL) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "dynamics.restitution_model",
                  "Unknown restitution model.");
    }
    if (dynamics->restitution_model == SIM_RESTITUTION_CONSTANT &&
        (!isfinite(dynamics->particle_restitution) ||
         dynamics->particle_restitution < 0.0 ||
         dynamics->particle_restitution > 1.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "dynamics.particle_restitution",
                  "Constant particle restitution must lie in [0, 1].");
    }
    if (dynamics->restitution_model == SIM_RESTITUTION_EXPONENTIAL) {
        if (!is_finite_positive(dynamics->exponential_velocity_scale)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.exponential_velocity_scale",
                      "Exponential velocity scale must be greater than zero.");
        }
        if (!isfinite(dynamics->exponential_amplitude) ||
            dynamics->exponential_amplitude <= 1.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.exponential_amplitude",
                      "The exponential collision model requires amplitude greater than one.");
        }
    }
    if (dynamics->injection < SIM_INJECTION_NONE ||
        dynamics->injection > SIM_INJECTION_ROUGH_ROTATIONAL) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.injection",
                  "Unknown collision energy-input model.");
    }
    if (dynamics->injection == SIM_INJECTION_ROUGH_ROTATIONAL &&
        !capabilities->rough_particles) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.injection",
                  "Rough rotational particles were not compiled into this build.");
    }
    if (capabilities->rough_particles &&
        dynamics->injection != SIM_INJECTION_ROUGH_ROTATIONAL) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.injection",
                  "This TANGENTIAL build always uses rough rotational contacts; rebuild without TANGENTIAL to select another contact law.");
    }
    if (dynamics->injection == SIM_INJECTION_ROUGH_ROTATIONAL &&
        capabilities->dimensions == 3) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.injection",
                  "The rough rotational collision model is currently 2D-only.");
    }
    if ((dynamics->injection == SIM_INJECTION_DELTA_TANGENT ||
         dynamics->injection == SIM_INJECTION_DOUBLE_DELTA_TANGENT) &&
        capabilities->dimensions == 3 && !capabilities->tangent_delta_3d) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.injection",
                  "Fixed-axis tangent delta is unavailable in this 3D build.");
    }
    if ((dynamics->injection == SIM_INJECTION_DELTA_TANGENT ||
         dynamics->injection == SIM_INJECTION_DOUBLE_DELTA_TANGENT ||
         (dynamics->injection == SIM_INJECTION_ROUGH_ROTATIONAL &&
          dynamics->rough_tangential_delta)) &&
        geometry->boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL &&
        !geometry->thermal_wall &&
        geometry->wall_tangent_restitution >= 1.0 - 1e-12) {
        add_issue(&writer, CONFIG_ISSUE_WARNING,
                  "dynamics.injection",
                  "Fixed-axis tangential delta injects one-signed angular momentum, while a specular circular/spherical wall provides no tangential torque sink. Strong driving can spin up without a steady state; use tangential wall dissipation, a thermal wall, or damping.");
    }
    if ((dynamics->injection == SIM_INJECTION_DELTA_NORMAL ||
         dynamics->injection == SIM_INJECTION_DELTA_TANGENT ||
         dynamics->injection == SIM_INJECTION_ACTIVE_TANGENTIAL ||
         (dynamics->injection == SIM_INJECTION_ROUGH_ROTATIONAL &&
          dynamics->rough_tangential_delta)) &&
        (!isfinite(dynamics->delta) || dynamics->delta < 0.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.delta",
                  "The active delta amplitude must be finite and nonnegative.");
    }
    if ((dynamics->injection == SIM_INJECTION_DOUBLE_DELTA_NORMAL ||
         dynamics->injection == SIM_INJECTION_DOUBLE_DELTA_TANGENT) &&
        (!isfinite(dynamics->delta_small) || dynamics->delta_small < 0.0 ||
         !isfinite(dynamics->delta_large) || dynamics->delta_large < 0.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.delta",
                  "Both active two-scale delta amplitudes must be finite and nonnegative.");
    }
    if ((dynamics->injection == SIM_INJECTION_DOUBLE_DELTA_NORMAL ||
         dynamics->injection == SIM_INJECTION_DOUBLE_DELTA_TANGENT) &&
        (!isfinite(dynamics->collision_memory_time) ||
         dynamics->collision_memory_time < 0.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "dynamics.collision_memory_time",
                  "Collision memory time must be finite and nonnegative.");
    }
    if (dynamics->injection == SIM_INJECTION_EVOLVING_DELTA) {
        if (!is_finite_positive(dynamics->evolving_delta_tau)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.evolving_delta_tau",
                      "Evolving-delta recovery time must be greater than zero.");
        }
        if (!isfinite(dynamics->evolving_delta_max) ||
            dynamics->evolving_delta_max < 0.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.evolving_delta_max",
                      "Maximum evolving delta must be finite and nonnegative.");
        }
    }
    if (dynamics->injection == SIM_INJECTION_COLLISION_ENERGY) {
        if (!isfinite(dynamics->injected_energy) ||
            dynamics->injected_energy < 0.0 ||
            !isfinite(dynamics->background_injected_energy) ||
            dynamics->background_injected_energy < 0.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.injected_energy",
                      "Collision and background injected energies must be finite and nonnegative.");
        }
        if (!is_finite_positive(dynamics->injection_recovery_time)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.injection_recovery_time",
                      "Energy-injection recovery time must be greater than zero.");
        }
        if (!is_finite_positive(dynamics->injection_exponent)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.injection_exponent",
                      "Energy-injection exponent must be greater than zero.");
        }
    }
    if (dynamics->injection == SIM_INJECTION_ROUGH_ROTATIONAL) {
        if (!isfinite(dynamics->rough_tangential_restitution) ||
            dynamics->rough_tangential_restitution < -1.0 ||
            dynamics->rough_tangential_restitution > 1.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.rough_tangential_restitution",
                      "Rough tangential restitution must lie in [-1, 1].");
        }
        if (!isfinite(dynamics->rough_angular_kick)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "dynamics.rough_angular_kick",
                      "Rough angular kick must be finite.");
        }
    }

    switch (driving->thermostat) {
        case SIM_THERMOSTAT_NONE:
        case SIM_THERMOSTAT_LANGEVIN:
        case SIM_THERMOSTAT_VELOCITY_RESCALE:
        case SIM_THERMOSTAT_REMOVE_COM:
        case SIM_THERMOSTAT_MAXWELL_REDRAW:
            break;
        default:
            add_issue(&writer, CONFIG_ISSUE_ERROR, "driving.thermostat",
                      "Unknown thermostat/forcing mode.");
            break;
    }
    if (driving->thermostat != SIM_THERMOSTAT_NONE &&
        !is_finite_positive(driving->thermostat_interval)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "driving.thermostat_interval",
                  "Thermostat interval must be greater than zero.");
    }
    if ((geometry->thermal_wall ||
         thermostat_uses_temperature(driving->thermostat)) &&
        (!isfinite(driving->temperature) || driving->temperature < 0.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "driving.temperature",
                  "The active bath temperature must be finite and nonnegative.");
    } else if ((geometry->thermal_wall ||
                driving->thermostat == SIM_THERMOSTAT_VELOCITY_RESCALE) &&
               driving->temperature == 0.0) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "driving.temperature",
                  geometry->thermal_wall
                      ? "A thermal wall requires a strictly positive temperature."
                      : "Global velocity rescaling requires a strictly positive target temperature.");
    }
    if (thermostat_uses_gamma(driving->thermostat) &&
        !is_finite_positive(driving->gamma)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "driving.gamma",
                  "The active thermostat gamma must be greater than zero.");
    }
    if (thermostat_uses_euler_choice(driving->thermostat) &&
        driving->euler_noise &&
        dynamics->flight == SIM_FLIGHT_LINEAR_DAMPING) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "driving.euler_noise",
                  "Euler noise cannot be combined with exact linear damping.");
    }
    if (dynamics->flight == SIM_FLIGHT_LINEAR_DAMPING &&
        thermostat_uses_gamma(driving->thermostat) &&
        fabs(dynamics->damping_gamma - driving->gamma) >
            1e-12*fmax(1.0, fabs(driving->gamma))) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "driving.gamma",
                  "The legacy engine shares one gamma between damping and thermostat; both values must match.");
    }
    if (potential->enabled) {
        if (!isfinite(potential->energy)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "potential.energy",
                      "Potential energy must be finite.");
        }
        if (!isfinite(potential->range) || potential->range <= 1.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "potential.range",
                      "Square-potential range must be greater than the hard-core diameter (sigma > 1). ");
        }
    }
    if (potential->enabled &&
        (potential->species_mode < SIM_POTENTIAL_ALL_PAIRS ||
         potential->species_mode > SIM_POTENTIAL_SAME_SPECIES)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "potential.species_mode",
                  "Unknown potential species mode.");
    }
    if (potential->enabled &&
        potential->species_mode == SIM_POTENTIAL_CHARGED) {
        double total = potential->positive_charge_fraction +
                       potential->neutral_charge_fraction;
        if (!isfinite(potential->positive_charge_fraction) ||
            potential->positive_charge_fraction < 0.0 ||
            !isfinite(potential->neutral_charge_fraction) ||
            potential->neutral_charge_fraction < 0.0 || total > 1.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "potential.charge_fractions",
                      "Positive and neutral fractions must be nonnegative and sum to at most one.");
        }
        if (init_mode_valid &&
            !init_mode_uses_generic_particle_builder(initialization->mode)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "potential.species_mode",
                      "Charge assignment is implemented only by random-growth and membrane initializers.");
        }
    }
    if (potential->enabled &&
        potential->species_mode == SIM_POTENTIAL_SAME_SPECIES &&
        init_mode_valid &&
        initialization->mode != SIM_INIT_DUMP_FILE &&
        initialization->mode != SIM_INIT_HEX_LATTICE &&
        initialization->composition != SIM_COMPOSITION_BINARY) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "potential.species_mode",
                  "Same-species potential requires binary particle types (the hex initializer creates them specifically for this rule).");
    }
    if (potential->enabled &&
        dynamics->flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR &&
        potential->species_mode == SIM_POTENTIAL_CHARGED) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "potential.species_mode",
                  "Charged interactions are incompatible with Lees-Edwards shear.");
    }

    collision_injection_active = injection_supplies_energy(dynamics->injection) ||
        (dynamics->injection == SIM_INJECTION_ROUGH_ROTATIONAL &&
         ((dynamics->rough_tangential_delta && dynamics->delta > 0.0) ||
          dynamics->rough_angular_kick != 0.0));
    dissipative =
        (dynamics->restitution_model == SIM_RESTITUTION_CONSTANT &&
         dynamics->particle_restitution < 1.0) ||
        dynamics->flight == SIM_FLIGHT_LINEAR_DAMPING ||
        (wall_normal_restitution_active &&
         geometry->wall_normal_restitution < 1.0) ||
        (wall_tangent_restitution_active &&
         fabs(geometry->wall_tangent_restitution) < 1.0);
    driven = collision_injection_active ||
        dynamics->restitution_model == SIM_RESTITUTION_EXPONENTIAL ||
        driving->thermostat != SIM_THERMOSTAT_NONE ||
        geometry->thermal_wall;
    if (dissipative && !driven) {
        add_issue(&writer, CONFIG_ISSUE_WARNING, "dynamics",
                  "The system is dissipative but has no energy input or thermostat.");
    }
    if (!dissipative && collision_injection_active) {
        add_issue(&writer, CONFIG_ISSUE_WARNING, "dynamics.injection",
                  "Collision injection in a fully elastic system causes unbounded heating.");
    }
    if (driving->thermostat != SIM_THERMOSTAT_NONE &&
        collision_injection_active) {
        add_issue(&writer, CONFIG_ISSUE_WARNING, "driving.thermostat",
                  "A thermostat and collision-based energy injection are both active.");
    }

    if (observables->pressure_field) {
        if (!is_finite_positive(observables->pressure_field_spacing)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.pressure_field_spacing",
                      "Pressure-field spacing must be greater than zero.");
        }
        if (capabilities->dimensions == 3 &&
            !capabilities->pressure_field_3d) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.pressure_field",
                      "The pressure-field implementation bins XY area only and is unavailable in 3D.");
        }
    }
    if (observables->temperature_profile &&
        (observables->temperature_profile_slabs < 1 ||
         observables->temperature_profile_slabs > 40)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "observables.temperature_profile_slabs",
                  "Temperature-profile slab count must lie in [1, 40].");
    }
    if (observables->velocity_profile &&
        (observables->velocity_profile_bins <= 0 ||
         observables->velocity_profile_bins > 4096)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "observables.velocity_profile_bins",
                  "Velocity-profile bin count must lie in [1, 4096].");
    }
    if (observables->interface_temperature &&
        !observables->interface_density) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "observables.interface_temperature",
                  "Interface temperature is a variant of the interface-density stream; enable Interface density as well.");
    }
    if (observables->osmosis) {
        if (!is_finite_positive(observables->osmosis_slice_width)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.osmosis_slice_width",
                      "Osmosis slice width must be greater than zero.");
        }
        if (capabilities->dimensions == 3) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "observables.osmosis",
                      "Osmosis density and energy normalization is currently 2D-only.");
        }
    }
    if (observables->structure_factor < SIM_STRUCTURE_FACTOR_OFF ||
        observables->structure_factor > SIM_STRUCTURE_FACTOR_VELOCITY) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "observables.structure_factor",
                  "Unknown saved structure-factor mode.");
    } else if (observables->structure_factor != SIM_STRUCTURE_FACTOR_OFF) {
        if (!is_finite_positive(observables->structure_factor_qmax)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.structure_factor_qmax",
                      "Structure-factor q maximum must be greater than zero.");
        }
        if (capabilities->dimensions == 3) {
            add_issue(&writer, CONFIG_ISSUE_INFO,
                      "observables.structure_factor",
                      "The structure factor is an XY projection in a 3D build.");
        }
    }
    if (observables->transverse_velocity_correlation) {
        if (observables->tvacf_mode_m == 0 &&
            observables->tvacf_mode_n == 0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.tvacf_mode_m",
                      "TVACF wavevector cannot be (0, 0).");
        }
        if (observables->tvacf_samples < 2 ||
            !is_finite_positive(observables->tvacf_time_span) ||
            observables->tvacf_checkpoint_interval <= 0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.tvacf_samples",
                      "TVACF requires at least two samples, a positive time span, and a positive checkpoint interval.");
        }
        if (output->logarithmic_thermo_spacing) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "output.logarithmic_thermo_spacing",
                      "TVACF requires uniform thermodynamic sampling.");
        }
        if (capabilities->dimensions == 3) {
            add_issue(&writer, CONFIG_ISSUE_INFO,
                      "observables.transverse_velocity_correlation",
                      "TVACF uses XY wavevectors and velocities in 3D.");
        }
    }
    if (observables->clusters) {
        if (observables->cluster_mode != 0 &&
            observables->cluster_mode != 1) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.cluster_mode",
                      "Cluster mode must be 0 (dense) or 1 (dilute).");
        }
        if (observables->cluster_neighbor_threshold < 0 ||
            !is_finite_positive(observables->cluster_cutoff)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.cluster_neighbor_threshold",
                      "Cluster threshold must be nonnegative and cutoff positive.");
        }
    }
    if (observables->bond_order < SIM_BOOP_OFF ||
        observables->bond_order > SIM_BOOP_CUTOFF) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "observables.bond_order",
                  "Unknown bond-orientational-order mode.");
    }
    if (capabilities->dimensions == 3 &&
        (observables->bond_order != SIM_BOOP_OFF ||
         observables->voronoi_area)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "observables.bond_order",
                  "Manual-snapshot BOOP and Voronoi columns are implemented only in 2D.");
    }
    if (capabilities->dimensions == 3 &&
        (observables->bond_order_pair_correlation ||
         observables->g6_pair_correlation)) {
        add_issue(&writer, CONFIG_ISSUE_INFO, "observables",
                  "Angular bond-order correlations are XY projections in 3D.");
    }
    if (observables->critical_regions ||
        observables->critical_snapshot_labels) {
        if (!isfinite(observables->critical_fraction) ||
            observables->critical_fraction == 0.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "observables.critical_fraction",
                      "Critical-region fraction must be finite and nonzero; a negative value requests automatic sizing.");
        }
    }
    if (observables->critical_snapshot_labels &&
        !observables->critical_regions) {
        add_issue(&writer, CONFIG_ISSUE_WARNING,
                  "observables.critical_snapshot_labels",
                  "Critical snapshot labels are enabled without critical-region thermodynamics.");
    }
    if (!is_finite_positive(output->maximum_time)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "output.maximum_time",
                  "Maximum simulation time must be greater than zero.");
    }
    if (output->periodic_dumps &&
        (!isfinite(output->first_dump_time) ||
         output->first_dump_time < -1.0 ||
         !is_finite_positive(output->dump_interval))) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "output.dump_interval",
                  "Dump interval must be positive and first dump time at least -1.");
    }
    if (output->periodic_dumps &&
        (!isfinite(output->next_dump_time_override) ||
         output->next_dump_time_override < -1.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "output.next_dump_time_override",
                  "Next-dump override must be -1 (automatic) or a nonnegative time.");
    }
    if (!isfinite(output->first_thermo_time) ||
        output->first_thermo_time < 0.0 ||
        !is_finite_positive(output->thermo_interval)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "output.thermo_interval",
                  "Thermodynamic sampling times are invalid.");
    }
    if (output->thermo_average_samples <= 0) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "output.thermo_average_samples",
                  "Thermodynamic averaging count must be greater than zero.");
    }
    if (output->logarithmic_thermo_spacing &&
        (!isfinite(output->logarithmic_spacing_base) ||
         output->logarithmic_spacing_base <= 1.0)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "output.logarithmic_spacing_base",
                  "Logarithmic spacing base must be greater than one.");
    }

    if (rendering->target_fps <= 0 ||
        rendering->target_fps > RENDER_TARGET_FPS_MAX) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "rendering.target_fps",
                  "Target frame rate must lie in [1, 1000].");
    }
    if (!isfinite(rendering->scene_supersampling) ||
        rendering->scene_supersampling < RENDER_SUPERSAMPLING_MIN ||
        rendering->scene_supersampling > RENDER_SUPERSAMPLING_MAX) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "rendering.scene_supersampling",
                  "Scene supersampling factor must lie in [0.5, 3].");
    }
    if (capabilities->dimensions == 3 &&
        (!isfinite(rendering->camera_fov_degrees) ||
         rendering->camera_fov_degrees < RENDER_CAMERA_FOV_MIN ||
         rendering->camera_fov_degrees > RENDER_CAMERA_FOV_MAX)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "rendering.camera_fov_degrees",
                  "Perspective field of view must lie in [10, 120] degrees.");
    }
    if (rendering->graphical_structure_factor &&
        (rendering->structure_factor_sample_cap <= 0 ||
         rendering->structure_factor_sample_cap >
            RENDER_STRUCTURE_FACTOR_SAMPLE_CAP_MAX)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "rendering.structure_factor_sample_cap",
                  "Structure-factor sample cap must lie in [1, 8192].");
    }
    if (rendering->graphical_structure_factor &&
        (!isfinite(rendering->structure_factor_refresh_seconds) ||
         rendering->structure_factor_refresh_seconds <
            RENDER_REFRESH_SECONDS_MIN ||
         rendering->structure_factor_refresh_seconds >
            RENDER_REFRESH_SECONDS_MAX)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "rendering.structure_factor_refresh_seconds",
                  "Structure-factor refresh interval must lie in [0.001, 3600] seconds.");
    }
    if (!isfinite(rendering->color_refresh_seconds) ||
        rendering->color_refresh_seconds < RENDER_REFRESH_SECONDS_MIN ||
        rendering->color_refresh_seconds > RENDER_REFRESH_SECONDS_MAX) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "rendering.color_refresh_seconds",
                  "Color refresh interval must lie in [0.001, 3600] seconds.");
    }
    if (rendering->color_mode < SIM_COLOR_BY_TYPE ||
        rendering->color_mode > SIM_COLOR_BY_CHARGE) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "rendering.color_mode",
                  "Unknown particle color mode.");
    }
    if (rendering->palette < SIM_PALETTE_PLASMA ||
        rendering->palette > SIM_PALETTE_COPPER) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "rendering.palette",
                  "Unknown continuous color palette.");
    }
    if (capabilities->dimensions == 3 &&
        (rendering->sphere_quality < SIM_SPHERE_QUALITY_AUTO ||
         rendering->sphere_quality > SIM_SPHERE_QUALITY_ULTRA)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR,
                  "rendering.sphere_quality",
                  "Unknown sphere tessellation quality.");
    }
    if (capabilities->dimensions == 3 &&
        (rendering->color_mode == SIM_COLOR_BY_HEX_ORDER ||
         rendering->color_mode == SIM_COLOR_BY_SQUARE_ORDER)) {
        add_issue(&writer, CONFIG_ISSUE_ERROR, "rendering.color_mode",
                  "Planar bond-order coloring is disabled in 3D.");
    }

    if (intruder->enabled) {
        if (dynamics->flight != SIM_FLIGHT_BALLISTIC) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "dynamics.flight",
                      "The polygonal intruder predictor currently requires ballistic particle flight.");
        }
        if (capabilities->stress || capabilities->heat) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "intruder.enabled",
                      "Intruder collisions do not yet maintain complete STRESS/HEAT transport bookkeeping.");
        }
        if (!capabilities->intruder) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "intruder.enabled",
                      "Intruder support was not compiled into this build.");
        }
        if (capabilities->dimensions == 3 && !capabilities->intruder_3d) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "intruder.enabled",
                      "The intruder collision model is currently 2D-only.");
        }
        if (!is_finite_positive(intruder->mass)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "intruder.mass",
                      "Intruder mass must be finite and greater than zero.");
        }
        if (!isfinite(intruder->moment_of_inertia)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "intruder.moment_of_inertia",
                      "Intruder moment of inertia must be finite; values at or below zero request automatic shape-based inertia.");
        }
        if (!isfinite(intruder->initial_vx) ||
            !isfinite(intruder->initial_vy) ||
            !isfinite(intruder->initial_angle) ||
            !isfinite(intruder->initial_omega)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "intruder.initial_state",
                      "Every applied intruder velocity, angle, and angular velocity must be finite.");
        }
        if (!isfinite(intruder->restitution) ||
            intruder->restitution < 0.0 || intruder->restitution > 1.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "intruder.restitution",
                      "Intruder restitution must lie in [0, 1].");
        }
        if (!isfinite(intruder->delta) || intruder->delta < 0.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "intruder.delta",
                      "Intruder tangential delta must be finite and nonnegative.");
        }
        if (!is_finite_positive(intruder->bath_particle_radius)) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "intruder.bath_particle_radius",
                      "Intruder bath-particle radius must be finite and greater than zero.");
        }
        if (!isfinite(intruder->display_particle_radius) ||
            intruder->display_particle_radius < 0.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "intruder.display_particle_radius",
                      "Intruder display-particle radius must be finite and nonnegative (zero disables the outline particles).");
        }
        if (!isfinite(intruder->kinetic_model_burn_time) ||
            intruder->kinetic_model_burn_time < 0.0) {
            add_issue(&writer, CONFIG_ISSUE_ERROR,
                      "intruder.kinetic_model_burn_time",
                      "Kinetic-model burn time must be finite and nonnegative.");
        }
        if (intruder->shape < 0 || intruder->shape > 7) {
            add_issue(&writer, CONFIG_ISSUE_ERROR, "intruder.shape",
                      "Unknown intruder shape.");
        } else {
            const double *parameter = intruder->shape_parameters;
            int parameter_count = 0;

            switch (intruder->shape) {
                case 0: /* SQUARE: side length */
                    parameter_count = 1;
                    break;
                case 1: /* TRIANGLE: base and height */
                case 6: /* CROSS: extension and width */
                    parameter_count = 2;
                    break;
                case 2: /* WHEEL: inner/outer radii and spike count */
                case 3: /* CHIRAL_WHEEL */
                    parameter_count = 3;
                    break;
                case 4: /* GEN_WHEEL: plus chirality */
                case 5: /* COMPLEX_WHEEL: four continuous parameters */
                    parameter_count = 4;
                    break;
                case 7: /* POLYGON: radius and side count */
                    parameter_count = 2;
                    break;
            }

            for (int index = 0; index < parameter_count; ++index) {
                if (!isfinite(parameter[index])) {
                    add_issue(&writer, CONFIG_ISSUE_ERROR,
                              "intruder.shape_parameters",
                              "Every parameter required by the selected intruder shape must be finite.");
                    break;
                }
            }

            if ((intruder->shape == 0 || intruder->shape == 1 ||
                 intruder->shape == 2 || intruder->shape == 3 ||
                 intruder->shape == 4 || intruder->shape == 6 ||
                 intruder->shape == 7) &&
                !is_finite_positive(parameter[0])) {
                add_issue(&writer, CONFIG_ISSUE_ERROR,
                          "intruder.shape_parameters[0]",
                          "The selected shape requires a positive size or inner radius.");
            }
            if ((intruder->shape == 1 || intruder->shape == 2 ||
                 intruder->shape == 3 || intruder->shape == 4 ||
                 intruder->shape == 6) &&
                !is_finite_positive(parameter[1])) {
                add_issue(&writer, CONFIG_ISSUE_ERROR,
                          "intruder.shape_parameters[1]",
                          "The selected shape requires a positive height, width, or outer radius.");
            }
            if (intruder->shape == 5 &&
                !is_finite_positive(parameter[0])) {
                add_issue(&writer, CONFIG_ISSUE_ERROR,
                          "intruder.shape_parameters[0]",
                          "The complex wheel requires a positive base radius.");
            }
            if (intruder->shape >= 2 && intruder->shape <= 4) {
                const double spikes = parameter[2];
                if (!isfinite(spikes) || floor(spikes) != spikes ||
                    spikes < 2.0 || spikes > 32.0) {
                    add_issue(&writer, CONFIG_ISSUE_ERROR,
                              "intruder.shape_parameters[2]",
                              "Wheel spike count must be an integer from 2 to 32 (at most 64 polygon vertices).");
                }
                if (is_finite_positive(parameter[0]) &&
                    is_finite_positive(parameter[1]) &&
                    parameter[1] < parameter[0]) {
                    add_issue(&writer, CONFIG_ISSUE_WARNING,
                              "intruder.shape_parameters[1]",
                              "The nominal outer wheel radius is smaller than the inner radius.");
                }
            }
            if (intruder->shape == 7) {
                const double sides = parameter[1];
                if (!isfinite(sides) || floor(sides) != sides ||
                    sides < 3.0 || sides > 64.0) {
                    add_issue(&writer, CONFIG_ISSUE_ERROR,
                              "intruder.shape_parameters[1]",
                              "Regular-polygon side count must be an integer from 3 to 64.");
                }
            }
        }
    }

    return writer.count;
}

bool simulation_config_equal(const SimulationConfig *left,
                             const SimulationConfig *right)
{
    if (left == right) return true;
    if (left == NULL || right == NULL) return false;
    /* All constructors in this module zero-initialize the object first.  This
     * makes the representation stable (including padding) and keeps equality
     * exhaustive as fields are added to the schema. */
    return memcmp(left, right, sizeof(*left)) == 0;
}
