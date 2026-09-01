#ifndef EDMD_LIVE_ANALYSIS_H
#define EDMD_LIVE_ANALYSIS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
    LIVE_ANALYSIS_HISTORY_CAPACITY = 512,
    LIVE_ANALYSIS_PAIR_MAX_BINS = 256
};

typedef struct LiveAnalysisSettings {
    bool energy;
    bool temperature;
    bool pressure;
    bool mean_squared_displacement;
    bool pair_correlation;
    double history_refresh_seconds;
    double pair_refresh_seconds;
    int pair_bins;
    int pair_particle_cap;
    /* Zero selects half the shortest box length. */
    double pair_max_radius;
} LiveAnalysisSettings;

typedef struct LiveAnalysisSeries {
    double simulation_time[LIVE_ANALYSIS_HISTORY_CAPACITY];
    double value[LIVE_ANALYSIS_HISTORY_CAPACITY];
    size_t count;
    size_t first;
} LiveAnalysisSeries;

typedef struct LivePairCorrelation {
    double radius[LIVE_ANALYSIS_PAIR_MAX_BINS];
    double value[LIVE_ANALYSIS_PAIR_MAX_BINS];
    int bins;
    int particles_sampled;
    uint64_t averages;
    double maximum_radius;
    double simulation_time;
    bool valid;
} LivePairCorrelation;

typedef struct LiveAnalysisData {
    uint64_t generation;
    LiveAnalysisSeries energy_per_particle;
    LiveAnalysisSeries temperature;
    LiveAnalysisSeries pressure;
    LiveAnalysisSeries mean_squared_displacement;
    LivePairCorrelation pair_correlation;
    bool pressure_supported;
    bool pressure_warming_up;
    bool msd_supported;
    bool msd_reference_ready;
    bool msd_allocation_failed;
    bool pair_correlation_supported;
    double msd_reference_time;
} LiveAnalysisData;

void live_analysis_initialize(void);
void live_analysis_shutdown(void);

/* Starts a fresh history generation while preserving the user's live-only
 * display settings. If MSD is enabled, its reference is captured immediately. */
void live_analysis_reset_generation(void);
void live_analysis_clear_histories(void);
void live_analysis_reset_msd_reference(void);
void live_analysis_reset_pair_average(void);

void live_analysis_get_settings(LiveAnalysisSettings *settings);
void live_analysis_apply_settings(const LiveAnalysisSettings *settings);

/* Called from the graphical frame loop. wall_time is monotonic wall-clock
 * time; history x coordinates always use simulation time. */
void live_analysis_update(double wall_time);
const LiveAnalysisData *live_analysis_get_data(void);

#ifdef __cplusplus
}
#endif

#endif
