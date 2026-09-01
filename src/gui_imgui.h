#ifndef EDMD_GUI_IMGUI_H
#define EDMD_GUI_IMGUI_H

#include "raylib.h"
#include "sim_config.h"

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum GuiCommandType {
    GUI_COMMAND_NONE = 0,
    /* Apply a validated dirty configuration as a live update.  The host may
     * update globals/events/cells incrementally; this is distinct from an
     * explicit restart at t = 0. */
    GUI_COMMAND_APPLY_CONFIG,
    GUI_COMMAND_RESTART_APPLIED,
    GUI_COMMAND_SAVE_SNAPSHOT,
    GUI_COMMAND_STEP_ONE_EVENT,
    GUI_COMMAND_UPDATE_PARTICLE,
    GUI_COMMAND_CLEAR_SELECTION,
    GUI_COMMAND_RESET_VIEW,
    GUI_COMMAND_TOGGLE_FULLSCREEN
} GuiCommandType;

typedef enum GuiConfigApplyOutcome {
    /* A live update succeeded.  Fault/selection/analysis state is preserved. */
    GUI_CONFIG_APPLY_SUCCEEDED = 0,
    /* A live update was rejected and the previous applied state was retained. */
    GUI_CONFIG_APPLY_REJECTED,
    /* The host completed an explicit or required full restart at t = 0. */
    GUI_CONFIG_APPLY_RESTARTED
} GuiConfigApplyOutcome;

typedef struct GuiParticleUpdate {
    int index;
    double vx;
    double vy;
    double vz;
    double mass;
} GuiParticleUpdate;

typedef struct GuiCommand {
    GuiCommandType type;
    /* Populated for APPLY_CONFIG and RESTART_APPLIED. */
    SimulationConfig config;
    /* Populated for UPDATE_PARTICLE. */
    GuiParticleUpdate particle;
} GuiCommand;

typedef struct GuiViewport {
    float x;
    float y;
    float width;
    float height;
    bool hovered;
    bool focused;
} GuiViewport;

typedef struct GuiMetrics {
    int dimensions;
    int particle_count;
    int visible_particle_count;
    int frustum_culled_particle_count;
    int cells_x;
    int cells_y;
    int cells_z;
    int queued_events;
    int event_capacity;
    double simulation_time;
    double box_lx;
    double box_ly;
    double box_lz;
    double packing_fraction;
    double total_energy;
    double energy_per_particle;
    double pressure;
    double bath_temperature;
    double temperature;
    double event_rate;
    double simulation_time_rate;
    double target_simulation_time_rate;
    double simulation_time_lag;
    uint64_t collisions;
    uint64_t crossings;
    float frames_per_second;
    bool growing;
    bool pressure_valid;
    bool simulation_time_pacing;
} GuiMetrics;

typedef struct GuiParticleSnapshot {
    bool valid;
    int index;
    int type;
    int charge;
    uint64_t collisions;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double radius;
    double mass;
    double last_collision_time;
} GuiParticleSnapshot;

/*
 * Rendering and interaction settings which may change immediately.  Physics
 * panels never write simulation globals; they edit a validated draft. Pointer
 * sliders/drags emit rate-limited live updates while held, whereas exact
 * numeric entry emits a debounced update only after editing has finished.
 */
typedef struct GuiLiveState {
    bool running;
    bool pace_simulation_time;
    double target_simulation_time_rate;
    bool show_scene;
    bool velocity_arrows;
    bool graphical_structure_factor;
    bool show_fps;
    bool show_help;
    bool show_axes;
    bool show_boundary;
    bool show_selection_outline;
    bool energy_log_log;
    bool temperature_log_log;
    bool pressure_log_log;
    bool msd_log_log;
    bool pair_correlation_log_log;
    SimColorMode color_mode;
    SimColorPalette palette;
    SimLightingMode lighting_mode;
    SimSphereQuality sphere_quality;
    int target_fps;
    double scene_supersampling;
    double camera_fov_degrees;
    int structure_factor_sample_cap;
    double structure_factor_refresh_seconds;
    double color_refresh_seconds;
    uint32_t type0_rgba;
    uint32_t type1_rgba;
    uint32_t background_rgba;
    uint32_t boundary_rgba;
} GuiLiveState;

typedef struct GuiFrameResult {
    GuiCommand command;
    GuiViewport viewport;
    bool wants_mouse;
    bool wants_keyboard;
} GuiFrameResult;

/* Must be called after InitWindow(). */
bool gui_imgui_init(const SimulationConfig *applied_config);

/* Must be called before CloseWindow(). */
void gui_imgui_shutdown(void);

/* Acknowledge a live configuration update or full restart.  For SUCCEEDED and
 * RESTARTED, applied_config must be non-null.  It may be null for REJECTED; a
 * non-null value records the effective configuration retained by the host.
 *
 * SUCCEEDED deliberately does not clear an engine fault, particle-selection
 * state, or live-analysis history.  RESTARTED retains the historical full
 * reset behavior.  For RESTARTED, a non-empty message requests a short,
 * automatically fading acknowledgement; rejection diagnostics remain visible
 * until the user dismisses them or a subsequent update succeeds. */
void gui_imgui_acknowledge_applied_config(
    const SimulationConfig *applied_config,
    GuiConfigApplyOutcome outcome,
    const char *message);

/* Backwards-compatible full-restart acknowledgement. */
void gui_imgui_set_applied_config(const SimulationConfig *applied_config);

/* Display a host/lifecycle failure in the workspace instead of only logging
 * it to the launching terminal. */
void gui_imgui_report_error(const char *message);

/* Report an engine state which requires a restart.  Unlike an ordinary
 * runtime warning, a fault keeps Run and Step disabled until a successful
 * restart clears it.  Both kinds of messages are non-modal. */
void gui_imgui_report_fault(const char *message);

/* True while the GUI is deliberately preventing event processing after an
 * engine fault. */
bool gui_imgui_simulation_faulted(void);

/*
 * Call between BeginDrawing() and EndDrawing().  The function owns the
 * matching rlImGuiBegin()/rlImGuiEnd() pair.
 */
GuiFrameResult gui_imgui_frame(const RenderTexture2D *scene,
                               const Texture2D *structure_factor,
                               const GuiMetrics *metrics,
                               const GuiParticleSnapshot *selected,
                               GuiLiveState *live);

#ifdef __cplusplus
}
#endif

#endif
