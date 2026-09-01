#ifndef EDMD_GUI_IMGUI_INTERNAL_HPP
#define EDMD_GUI_IMGUI_INTERNAL_HPP

#include "gui_imgui.h"

#include <cstddef>

namespace edmd_gui {

constexpr std::size_t kIssueCapacity = 192;
constexpr std::size_t kSearchCapacity = 128;
constexpr std::size_t kRuntimeErrorCapacity = 512;
constexpr std::size_t kToastMessageCapacity = 128;

struct Context {
    SimulationConfig applied{};
    SimulationConfig draft{};
    BuildCapabilities capabilities{};
    ConfigIssue issues[kIssueCapacity]{};
    std::size_t stored_issue_count = 0;
    std::size_t total_issue_count = 0;
    char search[kSearchCapacity]{};
    char runtime_error[kRuntimeErrorCapacity]{};
    char toast_message[kToastMessageCapacity]{};
    bool initialized = false;
    bool live_initialized = false;
    bool reset_dock_layout = true;
    bool runtime_message_visible = false;
    bool runtime_message_from_apply = false;
    bool toast_visible = false;
    bool simulation_faulted = false;
    bool automatic_apply_pending = false;
    bool automatic_apply_blocked = false;
    bool apply_command_outstanding = false;
    double last_physics_edit_time = 0;
    double toast_started_at = 0;
    /* Set only by pointer-driven numeric controls in the configuration
       panel.  Exact text entry deliberately remains on the ordinary debounce
       path so partially typed values can never update the simulation. */
    bool reactive_numeric_active = false;
    bool reactive_numeric_deactivated_after_edit = false;
    double last_reactive_apply_time = -1.0e30;
    bool particle_edit_loaded = false;
    int particle_edit_index = -1;
    double particle_vx = 0;
    double particle_vy = 0;
    double particle_vz = 0;
    double particle_mass = 1;
};

bool matches_filter(const Context &context, const char *label,
                    const char *keywords = nullptr);
void help_marker(const char *text);
void disabled_reason(const char *reason);
void refresh_validation(Context &context);
bool has_validation_errors(const Context &context);
void draw_issue_summary(Context &context, bool detailed);
void draw_configuration(Context &context, GuiLiveState &live);
void draw_live_analysis(Context &context, GuiLiveState &live);
void copy_live_to_rendering(const GuiLiveState &live,
                            SimRenderingConfig &rendering);
void copy_rendering_to_live(const SimRenderingConfig &rendering,
                            GuiLiveState &live);

} // namespace edmd_gui

#endif
