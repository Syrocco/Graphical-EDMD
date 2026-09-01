#include "live_analysis.h"

#include "EDMD.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern int N;
extern particle *particles;
extern double t;
extern double Lx, Ly, halfLx, halfLy;
#if THREE_D
extern double Lz;
#endif
extern bool addWell;
extern bool addField;
extern bool charged;
extern int liquidliquid;
extern double U;
extern double field;
extern int addWallx, addWally, addCircularWall, addMidWall;
#if THREE_D
extern int addWallz;
#endif
extern bool addShear;
extern double strainRate;
extern double collTermX, collTermY;
#if THREE_D
extern double collTermZ;
#endif
extern double lastScreen;
extern double cumulativeWallImpulse;
extern void (*freeFly)(particle *);

typedef struct LiveAnalysisState {
    LiveAnalysisSettings settings;
    LiveAnalysisData data;
    bool initialized;
    int observed_particles;
    int observed_dimensions;
    double observed_lx, observed_ly, observed_lz;
    double observed_time;
    void (*observed_free_fly)(particle *);
    double last_history_wall_time;
    double last_history_simulation_time;
    double last_pair_wall_time;

    double *msd_reference;
    double *msd_current;
    size_t msd_capacity;
    double msd_reference_com[3];

    double *pair_positions;
    size_t pair_position_capacity;

    bool pressure_baseline_ready;
    double pressure_time;
    double pressure_last_engine_reset;
    double pressure_virial[3];
    double pressure_wall_impulse;
} LiveAnalysisState;

static LiveAnalysisState state;

static int dimensions(void)
{
    return THREE_D ? 3 : 2;
}

static void clear_series(LiveAnalysisSeries *series)
{
    if (series == NULL) return;
    series->count = 0;
    series->first = 0;
}

static void append_series(LiveAnalysisSeries *series, double time,
                          double value)
{
    if (series == NULL || !isfinite(time) || !isfinite(value)) return;
    size_t index;
    if (series->count < LIVE_ANALYSIS_HISTORY_CAPACITY) {
        index = (series->first + series->count) %
                LIVE_ANALYSIS_HISTORY_CAPACITY;
        series->count++;
    } else {
        index = series->first;
        series->first = (series->first + 1) %
                        LIVE_ANALYSIS_HISTORY_CAPACITY;
    }
    series->simulation_time[index] = time;
    series->value[index] = value;
}

static particlePrediction predicted_particle(const particle *value)
{
    if (freeFly == &freeFlyGrow) return predictFreeFlyGrow(value, t);
    return predictFreeFlyNormal(value, t);
}

static double unwrapped_coordinate(double coordinate, int crossings,
                                   double length)
{
    return coordinate + crossings*length;
}

static bool ensure_msd_storage(void)
{
    const size_t needed = (size_t)N*(size_t)dimensions();
    if (needed == 0) return false;
    if (needed <= state.msd_capacity && state.msd_reference != NULL &&
        state.msd_current != NULL) {
        return true;
    }
    double *reference = realloc(state.msd_reference, needed*sizeof(double));
    if (reference == NULL) {
        state.data.msd_allocation_failed = true;
        return false;
    }
    state.msd_reference = reference;
    double *current = realloc(state.msd_current, needed*sizeof(double));
    if (current == NULL) {
        state.data.msd_allocation_failed = true;
        return false;
    }
    state.msd_current = current;
    state.msd_capacity = needed;
    state.data.msd_allocation_failed = false;
    return true;
}

static void fill_unwrapped_positions(double *destination, double com[3])
{
    const int dims = dimensions();
    com[0] = com[1] = com[2] = 0.0;
    for (int index = 0; index < N; ++index) {
        const particle *value = particles + index;
        const particlePrediction prediction = predicted_particle(value);
        const double x = unwrapped_coordinate(prediction.x, value->crossX, Lx);
        const double y = unwrapped_coordinate(prediction.y, value->crossY, Ly);
        destination[(size_t)index*dims] = x;
        destination[(size_t)index*dims + 1] = y;
        com[0] += x;
        com[1] += y;
#if THREE_D
        const double z = unwrapped_coordinate(prediction.z, value->crossZ, Lz);
        destination[(size_t)index*dims + 2] = z;
        com[2] += z;
#endif
    }
    if (N > 0) {
        com[0] /= N;
        com[1] /= N;
        com[2] /= N;
    }
}

static void capture_msd_reference(void)
{
    state.data.msd_supported = !addShear;
    state.data.msd_reference_ready = false;
    state.data.msd_allocation_failed = false;
    clear_series(&state.data.mean_squared_displacement);
    if (!state.settings.mean_squared_displacement ||
        !state.data.msd_supported || N <= 0 ||
        particles == NULL || !ensure_msd_storage()) {
        return;
    }
    fill_unwrapped_positions(state.msd_reference,
                             state.msd_reference_com);
    state.data.msd_reference_time = t;
    state.data.msd_reference_ready = true;
    append_series(&state.data.mean_squared_displacement, t, 0.0);
}

static bool sample_msd(double *result)
{
    if (result == NULL || !state.data.msd_reference_ready ||
        !ensure_msd_storage() || N <= 0) {
        return false;
    }
    double current_com[3];
    fill_unwrapped_positions(state.msd_current, current_com);
    const int dims = dimensions();
    double sum = 0.0;
    for (int index = 0; index < N; ++index) {
        for (int axis = 0; axis < dims; ++axis) {
            const size_t offset = (size_t)index*dims + axis;
            const double displacement = state.msd_current[offset] -
                state.msd_reference[offset] -
                (current_com[axis] - state.msd_reference_com[axis]);
            sum += displacement*displacement;
        }
    }
    *result = sum/N;
    return isfinite(*result);
}

static double pair_step_energy(const particle *first,
                               const particle *second)
{
    if (charged) {
        const int product = first->charge*second->charge;
        if (product < 0) return -fabs(U);
        if (product > 0) return fabs(U);
        return 0.0;
    }
    if (liquidliquid && first->type != second->type) return 0.0;
    return U;
}

static void sample_energetics(double *total_energy,
                              double *translational_temperature,
                              double *kinetic_energy)
{
    double translational = 0.0;
    double rotational = 0.0;
    double potential = 0.0;
    for (int index = 0; index < N; ++index) {
        const particle *value = particles + index;
        const particlePrediction prediction = predicted_particle(value);
        double vx = prediction.vx;
        if (addShear) vx += strainRate*(prediction.y - halfLy);
        double speed2 = vx*vx +
                        prediction.vy*prediction.vy;
#if THREE_D
        speed2 += prediction.vz*prediction.vz;
#endif
        translational += 0.5*value->m*speed2;
        if (addField) potential -= value->m*field*prediction.y;
#if TANGENTIAL
        rotational += 0.5*value->J*value->omega*value->omega;
#endif
    }

    if (addWell && U != 0.0) {
        for (int index = 0; index < N; ++index) {
            const particle *first = particles + index;
            for (int member = 0;
                 member < first->numberOfParticlesInWell; ++member) {
                const int partner = first->particlesInWell[member];
                if (partner > index && partner < N) {
                    potential += pair_step_energy(first,
                                                  particles + partner);
                }
            }
        }
    }
    *kinetic_energy = translational;
    *total_energy = translational + rotational + potential;
    *translational_temperature = N > 0
        ? 2.0*translational/((double)dimensions()*N) : 0.0;
}

static double container_volume(void)
{
    if (addCircularWall) {
#if THREE_D
        return (4.0/3.0)*M_PI*halfLx*halfLx*halfLx;
#else
        return M_PI*halfLx*halfLx;
#endif
    }
    double volume = Lx*Ly;
#if THREE_D
    volume *= Lz;
#endif
    return volume;
}

static double current_virial_sum(void)
{
    double result = collTermX + collTermY;
#if THREE_D
    result += collTermZ;
#endif
    return result;
}

static void capture_pressure_baseline(void)
{
    state.pressure_baseline_ready = true;
    state.pressure_time = t;
    state.pressure_last_engine_reset = lastScreen;
    state.pressure_virial[0] = collTermX;
    state.pressure_virial[1] = collTermY;
#if THREE_D
    state.pressure_virial[2] = collTermZ;
#else
    state.pressure_virial[2] = 0.0;
#endif
    state.pressure_wall_impulse = cumulativeWallImpulse;
    state.data.pressure_warming_up = true;
}

static double active_wall_measure(void)
{
    double measure = 0.0;
#if THREE_D
    if (addCircularWall) measure += 4.0*M_PI*halfLx*halfLx;
    if (addWallx) measure += 2.0*Ly*Lz;
    if (addWally) measure += 2.0*Lx*Lz;
    if (addWallz) measure += 2.0*Lx*Ly;
#else
    if (addCircularWall) measure += M_PI*Lx;
    if (addWallx) measure += 2.0*Ly;
    if (addWally) measure += 2.0*Lx;
#endif
    return measure;
}

static bool sample_pressure(double kinetic_energy, double *result)
{
    const bool physical_confinement = addWallx || addWally ||
                                      addCircularWall
#if THREE_D
                                      || addWallz
#endif
                                      ;
    state.data.pressure_supported = !addMidWall && !addShear &&
                                    freeFly != &freeFlyGrow;
    if (!state.data.pressure_supported || result == NULL) {
        capture_pressure_baseline();
        return false;
    }
    if (!state.pressure_baseline_ready) {
        capture_pressure_baseline();
        return false;
    }

    double start_time = state.pressure_time;
    if (physical_confinement) {
        const double elapsed = t - start_time;
        const double wall_measure = active_wall_measure();
        const double impulse_delta = cumulativeWallImpulse -
                                     state.pressure_wall_impulse;

        state.pressure_time = t;
        state.pressure_last_engine_reset = lastScreen;
        state.pressure_wall_impulse = cumulativeWallImpulse;
        state.pressure_virial[0] = collTermX;
        state.pressure_virial[1] = collTermY;
#if THREE_D
        state.pressure_virial[2] = collTermZ;
#else
        state.pressure_virial[2] = 0.0;
#endif
        if (!(elapsed > DBL_EPSILON) || !(wall_measure > 0.0) ||
            impulse_delta < 0.0) {
            state.data.pressure_warming_up = true;
            return false;
        }
        *result = impulse_delta/(elapsed*wall_measure);
        state.data.pressure_warming_up = !isfinite(*result);
        return isfinite(*result);
    }

    double previous_sum = state.pressure_virial[0] +
                          state.pressure_virial[1] +
                          state.pressure_virial[2];
    if (lastScreen != state.pressure_last_engine_reset) {
        /* saveThermo() resets the engine virial at lastScreen. Use only the
         * post-reset interval instead of mixing two accumulation epochs. */
        start_time = lastScreen;
        previous_sum = 0.0;
    }
    const double elapsed = t - start_time;
    const double volume = container_volume();
    const double virial_delta = current_virial_sum() - previous_sum;

    state.pressure_time = t;
    state.pressure_last_engine_reset = lastScreen;
    state.pressure_virial[0] = collTermX;
    state.pressure_virial[1] = collTermY;
#if THREE_D
    state.pressure_virial[2] = collTermZ;
#else
    state.pressure_virial[2] = 0.0;
#endif

    if (!(elapsed > DBL_EPSILON) || !(volume > 0.0)) {
        state.data.pressure_warming_up = true;
        return false;
    }
    *result = 2.0*kinetic_energy/(dimensions()*volume) -
              virial_delta/(dimensions()*volume*elapsed);
    state.data.pressure_warming_up = !isfinite(*result);
    return isfinite(*result);
}

static bool ensure_pair_storage(int particles_to_sample)
{
    const size_t needed = (size_t)particles_to_sample*(size_t)dimensions();
    if (needed == 0) return false;
    if (needed <= state.pair_position_capacity &&
        state.pair_positions != NULL) return true;
    double *resized = realloc(state.pair_positions,
                              needed*sizeof(double));
    if (resized == NULL) return false;
    state.pair_positions = resized;
    state.pair_position_capacity = needed;
    return true;
}

static double automatic_pair_radius(void)
{
    double shortest = fmin(Lx, Ly);
#if THREE_D
    shortest = fmin(shortest, Lz);
#endif
    return 0.5*shortest;
}

static void sample_pair_correlation(void)
{
    LivePairCorrelation *pair = &state.data.pair_correlation;
    state.data.pair_correlation_supported = !addShear;
    if (!state.data.pair_correlation_supported) {
        pair->valid = false;
        return;
    }
    const int bins = state.settings.pair_bins;
    const int sample_count = N < state.settings.pair_particle_cap
        ? N : state.settings.pair_particle_cap;
    const double automatic_radius = automatic_pair_radius();
    double maximum_radius = state.settings.pair_max_radius > 0.0
        ? state.settings.pair_max_radius : automatic_radius;
    maximum_radius = fmin(maximum_radius, automatic_radius);
    const double volume = container_volume();
    if (bins <= 0 || sample_count < 2 || !(maximum_radius > 0.0) ||
        !(volume > 0.0) || !ensure_pair_storage(sample_count)) {
        pair->valid = false;
        return;
    }

    const int dims = dimensions();
    for (int sample = 0; sample < sample_count; ++sample) {
        const int particle_index = (int)(((size_t)sample*(size_t)N)/
                                         (size_t)sample_count);
        const particlePrediction prediction =
            predicted_particle(particles + particle_index);
        state.pair_positions[(size_t)sample*dims] = prediction.x;
        state.pair_positions[(size_t)sample*dims + 1] = prediction.y;
#if THREE_D
        state.pair_positions[(size_t)sample*dims + 2] = prediction.z;
#endif
    }

    double histogram[LIVE_ANALYSIS_PAIR_MAX_BINS] = {0.0};
    const double inverse_bin_width = bins/maximum_radius;
    for (int first = 0; first < sample_count; ++first) {
        for (int second = first + 1; second < sample_count; ++second) {
            double dx = state.pair_positions[(size_t)second*dims] -
                        state.pair_positions[(size_t)first*dims];
            double dy = state.pair_positions[(size_t)second*dims + 1] -
                        state.pair_positions[(size_t)first*dims + 1];
#if THREE_D
            double dz = state.pair_positions[(size_t)second*dims + 2] -
                        state.pair_positions[(size_t)first*dims + 2];
            PBC(&dx, &dy, &dz);
            const double radius = sqrt(dx*dx + dy*dy + dz*dz);
#else
            PBC(&dx, &dy);
            const double radius = sqrt(dx*dx + dy*dy);
#endif
            const int bin = (int)(radius*inverse_bin_width);
            if (bin >= 0 && bin < bins) histogram[bin] += 1.0;
        }
    }

    const uint64_t next_average = pair->averages + 1;
    const double pair_count = 0.5*sample_count*(sample_count - 1.0);
    const double bin_width = maximum_radius/bins;
    for (int bin = 0; bin < bins; ++bin) {
        const double lower = bin*bin_width;
        const double upper = lower + bin_width;
#if THREE_D
        const double shell = (4.0/3.0)*M_PI*
            (upper*upper*upper - lower*lower*lower);
#else
        const double shell = M_PI*(upper*upper - lower*lower);
#endif
        const double expected = pair_count*shell/volume;
        const double sample_value = expected > 0.0
            ? histogram[bin]/expected : 0.0;
        if (pair->averages == 0) {
            pair->value[bin] = sample_value;
        } else {
            pair->value[bin] +=
                (sample_value - pair->value[bin])/(double)next_average;
        }
        pair->radius[bin] = lower + 0.5*bin_width;
    }
    pair->bins = bins;
    pair->particles_sampled = sample_count;
    pair->averages = next_average;
    pair->maximum_radius = maximum_radius;
    pair->simulation_time = t;
    pair->valid = true;
}

void live_analysis_initialize(void)
{
    if (state.initialized) return;
    memset(&state, 0, sizeof(state));
    state.settings.energy = true;
    state.settings.temperature = true;
    state.settings.pressure = true;
    state.settings.history_refresh_seconds = 0.25;
    state.settings.pair_refresh_seconds = 1.0;
    state.settings.pair_bins = 96;
    state.settings.pair_particle_cap = 512;
    state.last_history_wall_time = -DBL_MAX;
    state.last_history_simulation_time = -DBL_MAX;
    state.last_pair_wall_time = -DBL_MAX;
    state.observed_time = t;
    state.observed_free_fly = freeFly;
    state.initialized = true;
    live_analysis_reset_generation();
}

void live_analysis_shutdown(void)
{
    free(state.msd_reference);
    free(state.msd_current);
    free(state.pair_positions);
    memset(&state, 0, sizeof(state));
}

void live_analysis_clear_histories(void)
{
    clear_series(&state.data.energy_per_particle);
    clear_series(&state.data.temperature);
    clear_series(&state.data.pressure);
    clear_series(&state.data.mean_squared_displacement);
}

void live_analysis_reset_pair_average(void)
{
    memset(&state.data.pair_correlation, 0,
           sizeof(state.data.pair_correlation));
    state.last_pair_wall_time = -DBL_MAX;
}

void live_analysis_reset_msd_reference(void)
{
    capture_msd_reference();
}

void live_analysis_reset_generation(void)
{
    state.data.generation++;
    live_analysis_clear_histories();
    live_analysis_reset_pair_average();
    state.pressure_baseline_ready = false;
    state.data.pressure_warming_up = true;
    state.observed_particles = N;
    state.observed_dimensions = dimensions();
    state.observed_lx = Lx;
    state.observed_ly = Ly;
#if THREE_D
    state.observed_lz = Lz;
#else
    state.observed_lz = 0.0;
#endif
    state.observed_time = t;
    state.observed_free_fly = freeFly;
    state.data.msd_supported = !addShear;
    state.data.pair_correlation_supported = !addShear;
    state.last_history_wall_time = -DBL_MAX;
    state.last_history_simulation_time = -DBL_MAX;
    if (state.settings.mean_squared_displacement) capture_msd_reference();
}

void live_analysis_get_settings(LiveAnalysisSettings *settings)
{
    if (settings != NULL) *settings = state.settings;
}

static double clamp_double(double value, double minimum, double maximum,
                           double fallback)
{
    if (!isfinite(value)) return fallback;
    if (value < minimum) return minimum;
    if (value > maximum) return maximum;
    return value;
}

void live_analysis_apply_settings(const LiveAnalysisSettings *settings)
{
    if (settings == NULL) return;
    LiveAnalysisSettings sanitized = *settings;
    sanitized.history_refresh_seconds = clamp_double(
        sanitized.history_refresh_seconds, 0.02, 60.0, 0.25);
    sanitized.pair_refresh_seconds = clamp_double(
        sanitized.pair_refresh_seconds, 0.1, 3600.0, 1.0);
    if (sanitized.pair_bins < 16) sanitized.pair_bins = 16;
    if (sanitized.pair_bins > LIVE_ANALYSIS_PAIR_MAX_BINS)
        sanitized.pair_bins = LIVE_ANALYSIS_PAIR_MAX_BINS;
    if (sanitized.pair_particle_cap < 32)
        sanitized.pair_particle_cap = 32;
    if (sanitized.pair_particle_cap > 4096)
        sanitized.pair_particle_cap = 4096;
    sanitized.pair_max_radius = clamp_double(
        sanitized.pair_max_radius, 0.0, 1e12, 0.0);

    const bool enable_msd = !state.settings.mean_squared_displacement &&
                            sanitized.mean_squared_displacement;
    const bool enable_energy = !state.settings.energy && sanitized.energy;
    const bool enable_temperature = !state.settings.temperature &&
                                    sanitized.temperature;
    const bool enable_pressure = !state.settings.pressure &&
                                 sanitized.pressure;
    const bool pair_changed =
        state.settings.pair_correlation != sanitized.pair_correlation ||
        state.settings.pair_bins != sanitized.pair_bins ||
        state.settings.pair_particle_cap != sanitized.pair_particle_cap ||
        state.settings.pair_max_radius != sanitized.pair_max_radius;
    state.settings = sanitized;

    if (enable_energy) clear_series(&state.data.energy_per_particle);
    if (enable_temperature) clear_series(&state.data.temperature);
    if (enable_pressure) {
        clear_series(&state.data.pressure);
        state.pressure_baseline_ready = false;
    }
    if (!state.settings.mean_squared_displacement) {
        state.data.msd_reference_ready = false;
        clear_series(&state.data.mean_squared_displacement);
    } else if (enable_msd) {
        capture_msd_reference();
    }
    if (pair_changed) live_analysis_reset_pair_average();
}

static bool observed_generation_changed(void)
{
#if THREE_D
    const double current_lz = Lz;
#else
    const double current_lz = 0.0;
#endif
    return N != state.observed_particles ||
        dimensions() != state.observed_dimensions ||
        Lx != state.observed_lx || Ly != state.observed_ly ||
        current_lz != state.observed_lz ||
        freeFly != state.observed_free_fly ||
        t + 1e-12 < state.observed_time;
}

void live_analysis_update(double wall_time)
{
    if (!state.initialized) live_analysis_initialize();
    if (particles == NULL || N <= 0) return;
    if (observed_generation_changed()) live_analysis_reset_generation();
    state.observed_time = t;

    const bool simulation_advanced =
        t != state.last_history_simulation_time;
    const bool energetics_requested = state.settings.energy ||
        state.settings.temperature || state.settings.pressure;
    const bool history_requested = energetics_requested ||
        state.settings.mean_squared_displacement;
    if (simulation_advanced && history_requested &&
        wall_time - state.last_history_wall_time >=
            state.settings.history_refresh_seconds) {
        double total_energy = 0.0;
        double temperature = 0.0;
        double kinetic_energy = 0.0;
        /* Dense square-well systems may have many active memberships. Avoid
           walking them merely because an independent MSD curve is enabled. */
        if (energetics_requested)
            sample_energetics(&total_energy, &temperature, &kinetic_energy);
        if (state.settings.energy)
            append_series(&state.data.energy_per_particle, t,
                          total_energy/N);
        if (state.settings.temperature)
            append_series(&state.data.temperature, t, temperature);
        if (state.settings.pressure) {
            double pressure = 0.0;
            if (sample_pressure(kinetic_energy, &pressure))
                append_series(&state.data.pressure, t, pressure);
        }
        if (state.settings.mean_squared_displacement) {
            double msd = 0.0;
            if (sample_msd(&msd))
                append_series(&state.data.mean_squared_displacement, t, msd);
        }
        state.last_history_wall_time = wall_time;
        state.last_history_simulation_time = t;
    }

    if (state.settings.pair_correlation &&
        wall_time - state.last_pair_wall_time >=
            state.settings.pair_refresh_seconds &&
        (state.data.pair_correlation.averages == 0 ||
         t != state.data.pair_correlation.simulation_time)) {
        sample_pair_correlation();
        state.last_pair_wall_time = wall_time;
    }
}

const LiveAnalysisData *live_analysis_get_data(void)
{
    return &state.data;
}
