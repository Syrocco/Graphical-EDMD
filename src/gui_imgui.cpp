#include "gui_imgui_internal.hpp"

#include "imgui.h"
#include "imgui_internal.h"
#include "live_analysis.h"
#include "rlImGui.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <new>

namespace {

edmd_gui::Context *g_context = nullptr;

/* Applying cells/events can be much more expensive than drawing a GUI frame.
   Fifteen updates per second feels continuous while coalescing the many
   mouse-motion samples produced by high-refresh-rate displays. */
constexpr double kReactiveApplyIntervalSeconds = 1.0/15.0;

ImVec4 severity_color(ConfigIssueSeverity severity) {
    switch (severity) {
        case CONFIG_ISSUE_ERROR:
            return ImVec4(0.95f, 0.30f, 0.28f, 1.0f);
        case CONFIG_ISSUE_WARNING:
            return ImVec4(0.95f, 0.68f, 0.20f, 1.0f);
        default:
            return ImVec4(0.35f, 0.70f, 0.95f, 1.0f);
    }
}

const char *severity_name(ConfigIssueSeverity severity) {
    switch (severity) {
        case CONFIG_ISSUE_ERROR: return "Error";
        case CONFIG_ISSUE_WARNING: return "Warning";
        default: return "Info";
    }
}

void setup_default_dock_layout(ImGuiID dockspace_id, const ImVec2 &size) {
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id,
                              ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(dockspace_id, size);

    ImGuiID center = dockspace_id;
    ImGuiID sidebar = 0;
    ImGui::DockBuilderSplitNode(center, ImGuiDir_Right, 0.34f, &sidebar,
                                &center);

    ImGuiID lower_sidebar = 0;
    ImGuiID configuration = sidebar;
    ImGui::DockBuilderSplitNode(configuration, ImGuiDir_Down, 0.40f,
                                &lower_sidebar, &configuration);

    ImGui::DockBuilderDockWindow("Simulation", center);
    ImGui::DockBuilderDockWindow("Configuration", configuration);
    ImGui::DockBuilderDockWindow("Inspector", lower_sidebar);
    ImGui::DockBuilderDockWindow("Diagnostics", lower_sidebar);
    ImGui::DockBuilderDockWindow("Live analysis", lower_sidebar);
    ImGui::DockBuilderDockWindow("Structure factor", lower_sidebar);
    ImGui::DockBuilderDockWindow("EDMD Help", lower_sidebar);
    ImGui::DockBuilderFinish(dockspace_id);
}

bool config_is_dirty(edmd_gui::Context &context) {
    return !simulation_config_equal(&context.draft, &context.applied);
}

void assign_command(GuiCommand &command, GuiCommandType type,
                    const SimulationConfig &config) {
    if (command.type != GUI_COMMAND_NONE) return;
    command.type = type;
    command.config = config;
}

void assign_restart_command(GuiCommand &command, GuiCommandType type,
                            const SimulationConfig &config, bool running) {
    SimulationConfig restart = config;
    restart.rendering.start_running = running;
    assign_command(command, type, restart);
}

void assign_restart(GuiCommand &command, const SimulationConfig &config,
                    bool running) {
    assign_restart_command(command, GUI_COMMAND_RESTART_APPLIED,
                           config, running);
}

void assign_config_apply(GuiCommand &command, const SimulationConfig &config,
                         bool running) {
    SimulationConfig update = config;
    /* Space may toggle the live run state later in the same frame than the
       panel copy.  Carry that state in the acknowledgement snapshot without
       giving the update restart semantics. */
    update.rendering.start_running = running;
    assign_command(command, GUI_COMMAND_APPLY_CONFIG, update);
}

void draw_toolbar(edmd_gui::Context &context, GuiLiveState &live,
                  const GuiMetrics *metrics, GuiCommand &command) {
    const bool invalid = edmd_gui::has_validation_errors(context);
    const bool dirty = config_is_dirty(context);

    ImGui::BeginDisabled(context.simulation_faulted);
    if (ImGui::Button(live.running ? "Pause" : "Run")) {
        live.running = !live.running;
    }
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip(context.simulation_faulted
            ? "Restart the simulation before resuming event processing."
            : "Pause or resume event processing (Space). Physics settings are unaffected.");
    }
    ImGui::SameLine();

    ImGui::Checkbox("Pace", &live.pace_simulation_time);
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip(
            "Tie simulated-time advancement to wall time. If event processing "
            "cannot keep up, the accumulated simulated-time lag is shown here.");
    }
    ImGui::SameLine();
    const double minimum_exact_rate = 1e-6;
    const double maximum_exact_rate = 1e9;
    const double minimum_slider_rate = 1e-3;
    const double maximum_slider_rate = 1e3;
    if (!std::isfinite(live.target_simulation_time_rate) ||
        live.target_simulation_time_rate <= 0.0) {
        live.target_simulation_time_rate = 1.0;
    }
    live.target_simulation_time_rate = std::clamp(
        live.target_simulation_time_rate,
        minimum_exact_rate, maximum_exact_rate);

    /* Six logarithmic decades cover interactive slow motion through fast
       observation without making the slider hypersensitive.  The adjacent
       exact entry deliberately retains the wider host-supported range. */
    double slider_rate = std::clamp(live.target_simulation_time_rate,
                                    minimum_slider_rate,
                                    maximum_slider_rate);
    ImGui::BeginDisabled(!live.pace_simulation_time);
    ImGui::SetNextItemWidth(138.0f);
    if (ImGui::SliderScalar("##pace-rate-slider", ImGuiDataType_Double,
                            &slider_rate, &minimum_slider_rate,
                            &maximum_slider_rate, "%.3g",
                            ImGuiSliderFlags_AlwaysClamp |
                            ImGuiSliderFlags_Logarithmic)) {
        live.target_simulation_time_rate = slider_rate;
    }
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        ImGui::SetTooltip(
            "Logarithmic simulated-time rate slider (0.001 to 1000 sim s/s)."
        );
    }
    ImGui::SameLine();
    ImGui::SetNextItemWidth(88.0f);
    if (ImGui::InputDouble("##pace-rate-input",
                           &live.target_simulation_time_rate,
                           0.0, 0.0, "%.7g")) {
        if (!std::isfinite(live.target_simulation_time_rate) ||
            live.target_simulation_time_rate <= 0.0) {
            live.target_simulation_time_rate = minimum_exact_rate;
        } else {
            live.target_simulation_time_rate = std::clamp(
                live.target_simulation_time_rate,
                minimum_exact_rate, maximum_exact_rate);
        }
    }
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        ImGui::SetTooltip(
            "Exact requested simulated-time units per real second "
            "(1e-6 to 1e9). This is a live pacing control and does not restart "
            "or rebuild the simulation.");
    }
    ImGui::SameLine();
    ImGui::TextDisabled("sim s/s");
    if (live.pace_simulation_time && metrics != nullptr) {
        ImGui::SameLine();
        const double lag = std::max(0.0, metrics->simulation_time_lag);
        const bool visibly_behind = lag >
            std::max(1e-9, live.target_simulation_time_rate*0.10);
        if (visibly_behind) {
            ImGui::TextColored(ImVec4(0.95f, 0.58f, 0.22f, 1.0f),
                               "lag %.4g sim s", lag);
        } else {
            ImGui::TextDisabled("lag %.3g sim s", lag);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Achieved %.4g sim s/s; requested %.4g.",
                              metrics->simulation_time_rate,
                              live.target_simulation_time_rate);
        }
    }
    ImGui::SameLine();

    ImGui::BeginDisabled(live.running || context.simulation_faulted);
    if (ImGui::Button("Step event")) {
        assign_command(command, GUI_COMMAND_STEP_ONE_EVENT, context.applied);
    }
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        ImGui::SetTooltip(context.simulation_faulted
            ? "Restart the simulation before stepping."
            : live.running
            ? "Pause the simulation before stepping one event."
            : "Commit exactly one scheduled event.");
    }
    ImGui::SameLine();

    if (ImGui::Button("Restart")) {
        assign_restart(command, context.applied, live.running);
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Rebuild at t = 0 and preserve the current run/pause state (R).");
    }
    ImGui::SameLine();

    if (ImGui::Button("Save Snapshot")) {
        assign_command(command, GUI_COMMAND_SAVE_SNAPSHOT, context.applied);
    }
    ImGui::SameLine();

    if (context.simulation_faulted) {
        ImGui::TextColored(ImVec4(0.95f, 0.30f, 0.28f, 1.0f),
                           "Restart required");
    } else if (invalid) {
        ImGui::TextColored(ImVec4(0.95f, 0.30f, 0.28f, 1.0f),
                           "Invalid edit - running previous configuration");
    } else if (context.automatic_apply_blocked && dirty) {
        if (ImGui::SmallButton("Retry live update")) {
            context.automatic_apply_blocked = false;
            assign_config_apply(command, context.draft, live.running);
            context.apply_command_outstanding = true;
        }
    } else if (dirty) {
        ImGui::TextColored(ImVec4(0.95f, 0.68f, 0.20f, 1.0f),
                           "Live update pending...");
    } else {
        ImGui::TextDisabled("Configuration live");
    }
}

void draw_live_analysis_window(edmd_gui::Context &context,
                               GuiLiveState &live) {
    if (ImGui::Begin("Live analysis"))
        edmd_gui::draw_live_analysis(context, live);
    ImGui::End();
}

ImGuiID draw_dock_host(edmd_gui::Context &context, GuiLiveState &live,
                       const GuiMetrics *metrics, GuiCommand &command) {
    const ImGuiViewport *viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->WorkPos);
    ImGui::SetNextWindowSize(viewport->WorkSize);
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));

    const ImGuiWindowFlags flags = ImGuiWindowFlags_MenuBar |
        ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
        ImGuiWindowFlags_NoNavFocus;
    ImGui::Begin("EDMD Workspace", nullptr, flags);
    ImGui::PopStyleVar(3);

    if (ImGui::BeginMenuBar()) {
        draw_toolbar(context, live, metrics, command);
        const float layout_width = ImGui::CalcTextSize("Layout").x + 28.0f;
        if (ImGui::GetContentRegionAvail().x > layout_width) {
            ImGui::SetCursorPosX(ImGui::GetWindowWidth() - layout_width);
        }
        if (ImGui::BeginMenu("Layout")) {
            if (ImGui::MenuItem("Reset default workspace")) {
                context.reset_dock_layout = true;
            }
            ImGui::MenuItem("Help", nullptr, &live.show_help);
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    const ImGuiID dockspace_id = ImGui::GetID("EDMD.MainDockspace");
    if (context.reset_dock_layout ||
        ImGui::DockBuilderGetNode(dockspace_id) == nullptr) {
        setup_default_dock_layout(dockspace_id, viewport->WorkSize);
        context.reset_dock_layout = false;
    }
    ImGui::DockSpace(dockspace_id, ImVec2(0, 0));
    ImGui::End();
    return dockspace_id;
}

void draw_simulation_window(const RenderTexture2D *scene,
                            const GuiLiveState &live,
                            const GuiMetrics *metrics,
                            GuiViewport &result) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("Simulation");

    const ImVec2 available = ImGui::GetContentRegionAvail();
    const ImVec2 content_pos = ImGui::GetCursorScreenPos();
    /* Always publish the current content rectangle.  In particular, this
     * clears the host's previous hover state when the scene is hidden or its
     * render target is temporarily unavailable. */
    result.x = content_pos.x;
    result.y = content_pos.y;
    result.width = std::max(0.0f, available.x);
    result.height = std::max(0.0f, available.y);
    result.hovered = false;
    result.focused = false;
    if (!live.show_scene) {
        const char *message = "Scene hidden in Appearance";
        const ImVec2 size = ImGui::CalcTextSize(message);
        ImGui::SetCursorPos(ImVec2(std::max(0.0f, (available.x-size.x)*0.5f),
                                  std::max(0.0f, (available.y-size.y)*0.5f)));
        ImGui::TextDisabled("%s", message);
    } else if (scene == nullptr || scene->texture.id == 0 ||
               scene->texture.width <= 0 || scene->texture.height <= 0) {
        ImGui::TextDisabled("The raylib scene target is not ready.");
    } else {
        /* The renderer rebuilds its target from this content rectangle on the
           next frame. Fill the dock node now instead of preserving the old
           target's aspect ratio: preserving it creates a feedback loop that
           permanently letterboxes the viewport after a dock resize. */
        const ImVec2 image_size(std::max(0.0f, available.x),
                                std::max(0.0f, available.y));
        const ImVec2 image_pos = ImGui::GetCursorScreenPos();
        ImGui::Image(static_cast<ImTextureID>(scene->texture.id), image_size,
                     ImVec2(0, 1), ImVec2(1, 0));

        result.x = image_pos.x;
        result.y = image_pos.y;
        result.width = image_size.x;
        result.height = image_size.y;
        result.hovered = ImGui::IsItemHovered();
        result.focused = ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows);

        if (live.show_fps && metrics != nullptr) {
            char fps[48];
            std::snprintf(fps, sizeof(fps), "%.1f FPS", metrics->frames_per_second);
            const ImVec2 padding(8.0f, 5.0f);
            const ImVec2 text_size = ImGui::CalcTextSize(fps);
            const ImVec2 top_left(image_pos.x + 10.0f, image_pos.y + 10.0f);
            const ImVec2 bottom_right(top_left.x + text_size.x + 2.0f*padding.x,
                                      top_left.y + text_size.y + 2.0f*padding.y);
            ImDrawList *draw_list = ImGui::GetWindowDrawList();
            draw_list->AddRectFilled(top_left, bottom_right,
                                     IM_COL32(15, 18, 22, 185), 5.0f);
            draw_list->AddText(ImVec2(top_left.x + padding.x,
                                     top_left.y + padding.y),
                               IM_COL32(235, 240, 245, 255), fps);
        }
    }
    ImGui::End();
    ImGui::PopStyleVar();
}

void table_double_row(const char *label, double value, const char *format) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::TextDisabled("%s", label);
    ImGui::TableSetColumnIndex(1);
    ImGui::Text(format, value);
}

void table_int_row(const char *label, long long value) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::TextDisabled("%s", label);
    ImGui::TableSetColumnIndex(1);
    ImGui::Text("%lld", value);
}

void draw_inspector(edmd_gui::Context &context,
                    const GuiParticleSnapshot *selected, int dimensions,
                    const GuiLiveState &live, GuiCommand &command) {
    const bool transport_build = context.capabilities.stress ||
        context.capabilities.heat;
    ImGui::Begin("Inspector");
    if (selected == nullptr || !selected->valid) {
        context.particle_edit_loaded = false;
        context.particle_edit_index = -1;
        ImGui::TextWrapped("Select a particle in the simulation viewport to inspect it.");
        ImGui::Separator();
        ImGui::TextDisabled("Pause the simulation after selecting a particle to edit its "
                            "velocity. Mass is also editable when transport observables "
                            "are not compiled in.");
        ImGui::End();
        return;
    }

    ImGui::Text("Particle %d", selected->index);
    ImGui::SameLine();
    ImGui::TextDisabled("type %d", selected->type);
    if (ImGui::BeginTable("particle-properties", 2,
                          ImGuiTableFlags_RowBg |
                          ImGuiTableFlags_BordersInnerV |
                          ImGuiTableFlags_SizingStretchProp)) {
        table_double_row("x", selected->x, "%.8g");
        table_double_row("y", selected->y, "%.8g");
        if (dimensions == 3) table_double_row("z", selected->z, "%.8g");
        table_double_row("vx", selected->vx, "%.8g");
        table_double_row("vy", selected->vy, "%.8g");
        if (dimensions == 3) table_double_row("vz", selected->vz, "%.8g");
        table_double_row("radius", selected->radius, "%.8g");
        table_double_row("mass", selected->mass, "%.8g");
        table_double_row("last collision", selected->last_collision_time, "%.8g");
        table_int_row("charge", selected->charge);
        table_int_row("collisions", static_cast<long long>(selected->collisions));
        ImGui::EndTable();
    }

    if (live.running || context.simulation_faulted ||
        !context.particle_edit_loaded ||
        context.particle_edit_index != selected->index) {
        context.particle_edit_loaded = true;
        context.particle_edit_index = selected->index;
        context.particle_vx = selected->vx;
        context.particle_vy = selected->vy;
        context.particle_vz = selected->vz;
        context.particle_mass = selected->mass;
    }

    ImGui::SeparatorText("Paused edit");
    if (context.simulation_faulted) {
        ImGui::TextColored(severity_color(CONFIG_ISSUE_ERROR),
                           "Particle edits are locked until restart.");
    }
    const bool particle_edit_disabled = live.running ||
                                        context.simulation_faulted;
    ImGui::BeginDisabled(particle_edit_disabled);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputDouble("vx##particle-edit", &context.particle_vx,
                       0.0, 0.0, "%.10g");
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputDouble("vy##particle-edit", &context.particle_vy,
                       0.0, 0.0, "%.10g");
    if (dimensions == 3) {
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputDouble("vz##particle-edit", &context.particle_vz,
                           0.0, 0.0, "%.10g");
    }
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::BeginDisabled(transport_build);
    ImGui::InputDouble("mass##particle-edit", &context.particle_mass,
                       0.0, 0.0, "%.10g");
    if (!std::isfinite(context.particle_mass) || context.particle_mass < 0.0)
        context.particle_mass = 0.0;
    ImGui::EndDisabled();
    if (transport_build &&
        ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        ImGui::SetTooltip("Mass edits are disabled because this build accumulates "
                          "stress or heat using the current particle mass.");
    }
    if (transport_build) context.particle_mass = selected->mass;
    const bool finite = std::isfinite(context.particle_vx) &&
        std::isfinite(context.particle_vy) &&
        std::isfinite(context.particle_vz) &&
        std::isfinite(context.particle_mass) && context.particle_mass > 0;
    ImGui::BeginDisabled(!finite);
    if (ImGui::Button("Apply particle")) {
        command.type = GUI_COMMAND_UPDATE_PARTICLE;
        command.particle.index = selected->index;
        command.particle.vx = context.particle_vx;
        command.particle.vy = context.particle_vy;
        command.particle.vz = dimensions == 3 ? context.particle_vz : 0.0;
        command.particle.mass = context.particle_mass;
    }
    ImGui::EndDisabled();
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        if (context.simulation_faulted)
            ImGui::SetTooltip("Restart the faulted simulation before changing a particle.");
        else if (live.running) ImGui::SetTooltip("Pause before editing a particle.");
        else if (!finite) ImGui::SetTooltip("Velocities must be finite and mass must be positive.");
        else ImGui::SetTooltip("The simulation host will synchronize and reschedule this particle once.");
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear selection")) {
        command.type = GUI_COMMAND_CLEAR_SELECTION;
        context.particle_edit_loaded = false;
        context.particle_edit_index = -1;
    }
    ImGui::End();
}

void draw_structure_factor(const Texture2D *texture, int dimensions,
                           GuiLiveState &live) {
    if (!live.graphical_structure_factor) return;
    if (!ImGui::Begin("Structure factor", &live.graphical_structure_factor)) {
        ImGui::End();
        return;
    }
    if (dimensions == 3) {
        ImGui::TextDisabled("XY projection of the 3D structure factor");
    }
    if (texture == nullptr || texture->id == 0 ||
        texture->width <= 0 || texture->height <= 0) {
        ImGui::TextDisabled("Structure-factor texture is not ready.");
    } else {
        const ImVec2 available = ImGui::GetContentRegionAvail();
        const float scale = std::max(0.0f, std::min(
            available.x/static_cast<float>(texture->width),
            available.y/static_cast<float>(texture->height)));
        const ImVec2 size(texture->width*scale, texture->height*scale);
        ImVec2 pos = ImGui::GetCursorScreenPos();
        pos.x += std::max(0.0f, (available.x-size.x)*0.5f);
        pos.y += std::max(0.0f, (available.y-size.y)*0.5f);
        ImGui::SetCursorScreenPos(pos);
        ImGui::Image(static_cast<ImTextureID>(texture->id), size);
    }
    ImGui::End();
}

void draw_diagnostics(edmd_gui::Context &context, const GuiMetrics *metrics) {
    ImGui::Begin("Diagnostics");
    if (metrics == nullptr) {
        ImGui::TextDisabled("Metrics are not available.");
    } else {
        ImGui::Text("t = %.8g", metrics->simulation_time);
        ImGui::SameLine();
        ImGui::TextDisabled("%dD  N=%d", metrics->dimensions,
                            metrics->particle_count);
        if (metrics->growing) {
            ImGui::SameLine();
            ImGui::TextColored(ImVec4(0.40f, 0.80f, 0.45f, 1), "growing");
        }
        if (ImGui::BeginTable("simulation-metrics", 2,
                              ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_BordersInnerV |
                              ImGuiTableFlags_SizingStretchProp)) {
            table_double_row("Lx", metrics->box_lx, "%.7g");
            table_double_row("Ly", metrics->box_ly, "%.7g");
            if (metrics->dimensions == 3)
                table_double_row("Lz", metrics->box_lz, "%.7g");
            table_double_row("packing fraction", metrics->packing_fraction, "%.6g");
            table_double_row("total energy", metrics->total_energy, "%.8g");
            table_double_row("energy / particle", metrics->energy_per_particle, "%.8g");
            table_double_row("kinetic temperature", metrics->temperature, "%.7g");
            table_double_row("bath / target temperature",
                             metrics->bath_temperature, "%.7g");
            if (metrics->pressure_valid) {
                table_double_row("pressure", metrics->pressure, "%.7g");
            } else {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::TextDisabled("pressure");
                ImGui::TableSetColumnIndex(1);
                ImGui::TextDisabled("not sampled in graphics mode");
            }
            table_double_row("events / second", metrics->event_rate, "%.4g");
            table_double_row("simulation time / second", metrics->simulation_time_rate, "%.4g");
            if (metrics->simulation_time_pacing) {
                table_double_row("requested simulation time / second",
                                 metrics->target_simulation_time_rate, "%.4g");
                table_double_row("simulation-time lag",
                                 metrics->simulation_time_lag, "%.6g");
            }
            table_double_row("render FPS", metrics->frames_per_second, "%.1f");
            if (metrics->dimensions == 3) {
                table_int_row("visible particles",
                              metrics->visible_particle_count);
                table_int_row("frustum-culled particles",
                              metrics->frustum_culled_particle_count);
            }
            table_int_row("collision / interaction events",
                          static_cast<long long>(metrics->collisions));
            table_int_row("cell crossings", static_cast<long long>(metrics->crossings));
            table_int_row("scheduled events", metrics->queued_events);
            table_int_row("fixed event slots", metrics->event_capacity);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextDisabled("cell grid");
            ImGui::TableSetColumnIndex(1);
            if (metrics->dimensions == 3)
                ImGui::Text("%d x %d x %d", metrics->cells_x,
                            metrics->cells_y, metrics->cells_z);
            else
                ImGui::Text("%d x %d", metrics->cells_x, metrics->cells_y);
            ImGui::EndTable();
        }
    }

    ImGui::SeparatorText("Configuration checks");
    edmd_gui::draw_issue_summary(context, true);
    ImGui::End();
}

void draw_help_window(GuiLiveState &live, int dimensions) {
    if (!live.show_help) return;
    ImGui::SetNextWindowSize(ImVec2(560, 500), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("EDMD Help", &live.show_help)) {
        ImGui::End();
        return;
    }
    ImGui::TextWrapped(
        "Configuration changes are validated and sent as live updates. "
        "Appearance settings are immediate; compatible geometry and dynamics "
        "changes preserve the current simulation state. If a new constraint "
        "cannot contain the current particles, the requested configuration is "
        "applied by rebuilding automatically at t = 0.");
    ImGui::SeparatorText("Workflow");
    ImGui::BulletText("Choose the initial state, dynamics, and boundaries.");
    ImGui::BulletText("Resolve red validation errors; the previous valid run remains available meanwhile.");
    ImGui::BulletText("Sliders and drags send live updates while held; exact numeric entry updates after editing finishes.");
    ImGui::BulletText("A small fading notice acknowledges an edit that rebuilt the particles at t = 0.");
    ImGui::BulletText("Restart explicitly rebuilds the last successfully applied configuration at t = 0.");
    ImGui::SeparatorText("Viewport");
    if (dimensions == 3) {
        ImGui::BulletText("Wheel zooms; right-drag orbits; Home or C resets the camera.");
    } else {
        ImGui::BulletText("Wheel zooms; right-drag pans; Home or C resets the view.");
    }
    ImGui::BulletText("R restarts; F toggles fullscreen; S saves a snapshot; Space toggles run/pause.");
    ImGui::BulletText("Pace ties simulated time to the requested simulated-time units per real second and reports lag.");
    ImGui::BulletText("Left-click inspects; left-drag pulls a hard particle with a damped force without teleporting or pausing.");
    ImGui::BulletText("Particle edits are locked after an engine fault until the simulation is restarted.");
    ImGui::SeparatorText("Unavailable controls");
    ImGui::TextWrapped(
        "Controls that do not exist in this compiled dimension are hidden. "
        "Build-dependent and temporarily incompatible choices remain visible; "
        "hover them to see the exact reason they cannot currently be selected.");
    ImGui::End();
}

void draw_rebuild_toast(edmd_gui::Context &context) {
    if (!context.toast_visible) return;

    constexpr double hold_seconds = 2.5;
    constexpr double fade_seconds = 3.5;
    const double elapsed = std::max(0.0,
        ImGui::GetTime() - context.toast_started_at);
    if (elapsed >= hold_seconds + fade_seconds) {
        context.toast_visible = false;
        return;
    }

    float alpha = 1.0f;
    if (elapsed > hold_seconds) {
        const float linear = static_cast<float>(
            1.0 - (elapsed - hold_seconds)/fade_seconds);
        /* Smoothstep avoids an abrupt change of opacity when fading starts. */
        alpha = linear*linear*(3.0f - 2.0f*linear);
    }

    const ImGuiViewport *viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(
        ImVec2(viewport->WorkPos.x + 16.0f,
               viewport->WorkPos.y + ImGui::GetFrameHeight() + 16.0f),
        ImGuiCond_Always);
    ImGui::SetNextWindowSizeConstraints(ImVec2(220.0f, 0.0f),
                                        ImVec2(380.0f, FLT_MAX));
    ImGui::SetNextWindowBgAlpha(0.90f*alpha);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, alpha);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 7.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(12.0f, 9.0f));
    const ImGuiWindowFlags flags = ImGuiWindowFlags_AlwaysAutoResize |
        ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoInputs |
        ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoSavedSettings;
    if (ImGui::Begin("##simulation-rebuild-toast", nullptr, flags)) {
        ImGui::TextColored(ImVec4(0.35f, 0.72f, 0.95f, 1.0f), "%s",
                           context.toast_message);
    }
    ImGui::End();
    ImGui::PopStyleVar(3);
}

void draw_runtime_status(edmd_gui::Context &context, GuiLiveState &live,
                         GuiCommand &command) {
    if (!context.runtime_message_visible) return;
    const ImGuiViewport *viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(
        ImVec2(viewport->WorkPos.x + viewport->WorkSize.x*0.5f,
               viewport->WorkPos.y + ImGui::GetFrameHeight() + 12.0f),
        ImGuiCond_Appearing, ImVec2(0.5f, 0.0f));
    ImGui::SetNextWindowSizeConstraints(ImVec2(360, 0), ImVec2(680, FLT_MAX));
    bool open = true;
    const ImGuiWindowFlags flags = ImGuiWindowFlags_AlwaysAutoResize |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings;
    if (ImGui::Begin("Simulation status", &open, flags)) {
        const ImVec4 color = context.simulation_faulted
            ? severity_color(CONFIG_ISSUE_ERROR)
            : severity_color(CONFIG_ISSUE_WARNING);
        ImGui::TextColored(color, "%s",
            context.simulation_faulted ? "Simulation paused - restart required"
                                       : "Simulation notice");
        ImGui::PushTextWrapPos(ImGui::GetFontSize()*42.0f);
        ImGui::TextUnformatted(context.runtime_error);
        ImGui::PopTextWrapPos();
        ImGui::Spacing();
        if (context.simulation_faulted) {
            if (ImGui::Button("Restart paused")) {
                assign_restart(command, context.applied, live.running);
            }
            ImGui::SameLine();
            ImGui::TextDisabled("Recovery stays paused; inspection and configuration remain available.");
        } else if (ImGui::Button("Dismiss")) {
            context.runtime_message_visible = false;
        }
    }
    ImGui::End();
    /* Fault notices remain available until a restart succeeds. */
    if (!open && !context.simulation_faulted) {
        context.runtime_message_visible = false;
    }
}

} // namespace

namespace edmd_gui {

void copy_live_to_rendering(const GuiLiveState &live,
                            SimRenderingConfig &rendering) {
    rendering.start_running = live.running;
    rendering.target_fps = live.target_fps;
    rendering.velocity_arrows = live.velocity_arrows;
    rendering.color_mode = live.color_mode;
    rendering.palette = live.palette;
    rendering.graphical_structure_factor = live.graphical_structure_factor;
    rendering.show_fps = live.show_fps;
    rendering.show_axes = live.show_axes;
    rendering.show_boundary = live.show_boundary;
    rendering.show_selection_outline = live.show_selection_outline;
    rendering.lighting_mode = live.lighting_mode;
    rendering.sphere_quality = live.sphere_quality;
    rendering.scene_supersampling = live.scene_supersampling;
    rendering.camera_fov_degrees = live.camera_fov_degrees;
    rendering.structure_factor_sample_cap = live.structure_factor_sample_cap;
    rendering.structure_factor_refresh_seconds = live.structure_factor_refresh_seconds;
    rendering.color_refresh_seconds = live.color_refresh_seconds;
    rendering.type0_rgba = live.type0_rgba;
    rendering.type1_rgba = live.type1_rgba;
    rendering.background_rgba = live.background_rgba;
    rendering.boundary_rgba = live.boundary_rgba;
}

void copy_rendering_to_live(const SimRenderingConfig &rendering,
                            GuiLiveState &live) {
    live.running = rendering.start_running;
    live.velocity_arrows = rendering.velocity_arrows;
    live.graphical_structure_factor = rendering.graphical_structure_factor;
    live.show_fps = rendering.show_fps;
    live.show_axes = rendering.show_axes;
    live.show_boundary = rendering.show_boundary;
    live.show_selection_outline = rendering.show_selection_outline;
    live.lighting_mode = rendering.lighting_mode;
    live.color_mode = rendering.color_mode;
    live.palette = rendering.palette;
    live.sphere_quality = rendering.sphere_quality;
    live.target_fps = rendering.target_fps;
    live.scene_supersampling = rendering.scene_supersampling;
    live.camera_fov_degrees = rendering.camera_fov_degrees;
    live.structure_factor_sample_cap = rendering.structure_factor_sample_cap;
    live.structure_factor_refresh_seconds = rendering.structure_factor_refresh_seconds;
    live.color_refresh_seconds = rendering.color_refresh_seconds;
    live.type0_rgba = rendering.type0_rgba;
    live.type1_rgba = rendering.type1_rgba;
    live.background_rgba = rendering.background_rgba;
    live.boundary_rgba = rendering.boundary_rgba;
}

void refresh_validation(Context &context) {
    std::memset(context.issues, 0, sizeof(context.issues));
    context.total_issue_count = simulation_config_validate(
        &context.draft, &context.capabilities, context.issues, kIssueCapacity);
    context.stored_issue_count = std::min(context.total_issue_count,
                                          kIssueCapacity);
}

bool has_validation_errors(const Context &context) {
    for (std::size_t i = 0; i < context.stored_issue_count; ++i) {
        if (context.issues[i].severity == CONFIG_ISSUE_ERROR) return true;
    }
    return context.total_issue_count > context.stored_issue_count;
}

void draw_issue_summary(Context &context, bool detailed) {
    int errors = 0;
    int warnings = 0;
    int infos = 0;
    for (std::size_t i = 0; i < context.stored_issue_count; ++i) {
        switch (context.issues[i].severity) {
            case CONFIG_ISSUE_ERROR: ++errors; break;
            case CONFIG_ISSUE_WARNING: ++warnings; break;
            default: ++infos; break;
        }
    }
    if (context.total_issue_count == 0) {
        ImGui::TextColored(ImVec4(0.40f, 0.80f, 0.45f, 1.0f),
                           "Configuration is valid");
        return;
    }
    ImGui::TextColored(errors ? severity_color(CONFIG_ISSUE_ERROR)
                              : severity_color(CONFIG_ISSUE_WARNING),
                       "%d error%s, %d warning%s, %d info",
                       errors, errors == 1 ? "" : "s",
                       warnings, warnings == 1 ? "" : "s", infos);
    if (!detailed) return;

    for (std::size_t i = 0; i < context.stored_issue_count; ++i) {
        const ConfigIssue &issue = context.issues[i];
        ImGui::PushID(static_cast<int>(i));
        ImGui::TextColored(severity_color(issue.severity), "%s",
                           severity_name(issue.severity));
        ImGui::SameLine();
        ImGui::TextWrapped("%s", issue.message);
        if (issue.field[0] != '\0') {
            ImGui::Indent();
            ImGui::TextDisabled("%s", issue.field);
            ImGui::Unindent();
        }
        ImGui::PopID();
    }
    if (context.total_issue_count > context.stored_issue_count) {
        ImGui::TextDisabled("%zu additional issues were omitted.",
                            context.total_issue_count-context.stored_issue_count);
    }
}

} // namespace edmd_gui

extern "C" bool gui_imgui_init(const SimulationConfig *applied_config) {
    if (g_context != nullptr) return true;
    g_context = new (std::nothrow) edmd_gui::Context();
    if (g_context == nullptr) return false;

    /* Build the font atlas at backing-store resolution instead of magnifying
       the 13 px bitmap font on Retina displays. */
    rlImGuiBeginInitImGui();
    ImGui::StyleColorsDark();
    ImGuiIO &io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags &= ~ImGuiConfigFlags_ViewportsEnable;

    ImGuiStyle &style = ImGui::GetStyle();
    style.WindowRounding = 5.0f;
    style.ChildRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.GrabRounding = 3.0f;
    style.TabRounding = 4.0f;
    style.WindowMenuButtonPosition = ImGuiDir_Right;
#if defined(__APPLE__) && !defined(__EMSCRIPTEN__)
    const Vector2 dpi = GetWindowScaleDPI();
    const float dpi_scale = std::isfinite(dpi.x) && std::isfinite(dpi.y)
        ? std::clamp(std::max(dpi.x, dpi.y), 1.0f, 2.5f) : 1.0f;
    ImFont *native_font = io.Fonts->AddFontFromFileTTF(
        "/System/Library/Fonts/SFNS.ttf", 15.0f*dpi_scale);
    if (native_font != nullptr) {
        io.FontDefault = native_font;
        /* ImGui 1.92 moved the application-wide font scale into the style.
           Rasterize at backing-store density, then retain 15 logical pixels. */
        style.FontScaleMain = 1.0f/dpi_scale;
    }
#endif
    rlImGuiEndInitImGui();

    g_context->capabilities = simulation_build_capabilities();
    if (applied_config != nullptr) {
        g_context->applied = *applied_config;
    } else {
        simulation_config_set_defaults(&g_context->applied,
                                       &g_context->capabilities);
    }
    g_context->draft = g_context->applied;
    g_context->initialized = true;
    edmd_gui::refresh_validation(*g_context);
    live_analysis_initialize();
    return true;
}

extern "C" void gui_imgui_shutdown(void) {
    if (g_context == nullptr) return;
    live_analysis_shutdown();
    rlImGuiShutdown();
    delete g_context;
    g_context = nullptr;
}

extern "C" void gui_imgui_acknowledge_applied_config(
    const SimulationConfig *applied_config,
    GuiConfigApplyOutcome outcome,
    const char *message) {
    if (g_context == nullptr) return;
    if (outcome != GUI_CONFIG_APPLY_REJECTED && applied_config == nullptr)
        return;

    switch (outcome) {
        case GUI_CONFIG_APPLY_SUCCEEDED:
            g_context->applied = *applied_config;
            g_context->draft = *applied_config;
            if (g_context->runtime_message_from_apply) {
                g_context->runtime_message_visible = false;
                g_context->runtime_message_from_apply = false;
            }
            g_context->automatic_apply_pending = false;
            g_context->automatic_apply_blocked = false;
            g_context->apply_command_outstanding = false;
            /* A compatible update is deliberately non-destructive: engine
               faults, selection-side editor state, and analysis history have
               independent lifetimes and are not reset here. */
            break;
        case GUI_CONFIG_APPLY_REJECTED:
            if (applied_config != nullptr)
                g_context->applied = *applied_config;
            g_context->automatic_apply_pending = false;
            g_context->automatic_apply_blocked = true;
            g_context->apply_command_outstanding = false;
            g_context->runtime_message_from_apply = true;
            g_context->toast_visible = false;
            break;
        case GUI_CONFIG_APPLY_RESTARTED:
            g_context->applied = *applied_config;
            g_context->draft = *applied_config;
            g_context->simulation_faulted = false;
            g_context->automatic_apply_pending = false;
            g_context->automatic_apply_blocked = false;
            g_context->apply_command_outstanding = false;
            g_context->runtime_message_visible = false;
            g_context->runtime_message_from_apply = false;
            g_context->toast_visible = false;
            live_analysis_reset_generation();
            break;
        default:
            return;
    }

    if (outcome == GUI_CONFIG_APPLY_RESTARTED &&
        message != nullptr && message[0] != '\0') {
        /* The detailed lifecycle diagnostic is useful in logs, but a
           successful automatic rebuild needs only a quiet acknowledgement. */
        std::snprintf(g_context->toast_message,
                      sizeof(g_context->toast_message), "%s",
                      "Simulation rebuilt for the new configuration");
        g_context->toast_started_at = ImGui::GetTime();
        g_context->toast_visible = true;
    } else if (message != nullptr && message[0] != '\0') {
        std::snprintf(g_context->runtime_error,
                      sizeof(g_context->runtime_error), "%s", message);
        g_context->runtime_message_visible = true;
        g_context->runtime_message_from_apply =
            outcome == GUI_CONFIG_APPLY_REJECTED;
    } else if (outcome == GUI_CONFIG_APPLY_REJECTED) {
        std::snprintf(g_context->runtime_error,
                      sizeof(g_context->runtime_error), "%s",
                      "Live update rolled back; the previous valid "
                      "configuration remains active.");
        g_context->runtime_message_visible = true;
        g_context->runtime_message_from_apply = true;
    }
    edmd_gui::refresh_validation(*g_context);
}

extern "C" void gui_imgui_set_applied_config(
    const SimulationConfig *applied_config) {
    gui_imgui_acknowledge_applied_config(
        applied_config, GUI_CONFIG_APPLY_RESTARTED, nullptr);
}

extern "C" void gui_imgui_report_error(const char *message) {
    if (g_context == nullptr) return;
    std::snprintf(g_context->runtime_error,
                  sizeof(g_context->runtime_error), "%s",
                  message != nullptr ? message : "Unknown simulation error");
    g_context->runtime_message_visible = true;
    g_context->runtime_message_from_apply = false;
    g_context->toast_visible = false;
    if (g_context->apply_command_outstanding) {
        g_context->automatic_apply_blocked = true;
        g_context->apply_command_outstanding = false;
    }
}

extern "C" void gui_imgui_report_fault(const char *message) {
    if (g_context == nullptr) return;
    std::snprintf(g_context->runtime_error,
                  sizeof(g_context->runtime_error), "%s",
                  message != nullptr ? message : "The simulation requires a restart");
    g_context->runtime_message_visible = true;
    g_context->runtime_message_from_apply = false;
    g_context->toast_visible = false;
    g_context->simulation_faulted = true;
}

extern "C" bool gui_imgui_simulation_faulted(void) {
    return g_context != nullptr && g_context->simulation_faulted;
}

extern "C" GuiFrameResult gui_imgui_frame(
    const RenderTexture2D *scene, const Texture2D *structure_factor,
    const GuiMetrics *metrics,
    const GuiParticleSnapshot *selected, GuiLiveState *live) {
    GuiFrameResult result{};
    result.command.type = GUI_COMMAND_NONE;
    if (g_context == nullptr || live == nullptr) return result;

    if (!g_context->live_initialized) {
        edmd_gui::copy_rendering_to_live(g_context->applied.rendering, *live);
        g_context->live_initialized = true;
    }
    edmd_gui::copy_live_to_rendering(*live, g_context->applied.rendering);
    edmd_gui::copy_live_to_rendering(*live, g_context->draft.rendering);
    edmd_gui::refresh_validation(*g_context);
    const SimulationConfig draft_before_frame = g_context->draft;
    g_context->reactive_numeric_active = false;
    g_context->reactive_numeric_deactivated_after_edit = false;

    live_analysis_update(ImGui::GetTime());
    rlImGuiBegin();
    draw_dock_host(*g_context, *live, metrics, result.command);
    draw_simulation_window(scene, *live, metrics, result.viewport);
    edmd_gui::draw_configuration(*g_context, *live);
    const int dimensions = metrics ? metrics->dimensions
                                   : g_context->capabilities.dimensions;
    draw_inspector(*g_context, selected, dimensions, *live, result.command);
    draw_diagnostics(*g_context, metrics);
    draw_live_analysis_window(*g_context, *live);
    draw_structure_factor(structure_factor, dimensions, *live);
    draw_help_window(*live, dimensions);

    edmd_gui::copy_live_to_rendering(*live, g_context->applied.rendering);
    edmd_gui::copy_live_to_rendering(*live, g_context->draft.rendering);
    edmd_gui::refresh_validation(*g_context);
    const bool physics_changed =
        !simulation_config_equal(&draft_before_frame, &g_context->draft);
    const bool dirty = config_is_dirty(*g_context);
    const bool invalid = edmd_gui::has_validation_errors(*g_context);
    const double now = ImGui::GetTime();
    if (physics_changed) {
        g_context->last_physics_edit_time = now;
        g_context->automatic_apply_pending = true;
        g_context->automatic_apply_blocked = false;
    }
    if (!dirty || invalid) {
        g_context->automatic_apply_pending = false;
    } else if (!g_context->automatic_apply_blocked &&
               !g_context->apply_command_outstanding) {
        g_context->automatic_apply_pending = true;
    }
    /* Read global shortcuts from ImGui's event stream.  rlImGui may consume
       the corresponding raylib edge before getInput() runs, which was why the
       documented shortcuts appeared inert. */
    const bool shortcut_available = !ImGui::GetIO().WantTextInput &&
                                    !ImGui::IsAnyItemActive();
    if (shortcut_available) {
        if (ImGui::IsKeyPressed(ImGuiKey_Space, false) &&
            !g_context->simulation_faulted) {
            live->running = !live->running;
        }
        if (ImGui::IsKeyPressed(ImGuiKey_R, false)) {
            assign_restart(result.command, g_context->applied, live->running);
        } else if (ImGui::IsKeyPressed(ImGuiKey_Home, false) ||
                   ImGui::IsKeyPressed(ImGuiKey_C, false)) {
            assign_command(result.command, GUI_COMMAND_RESET_VIEW,
                           g_context->applied);
        } else if (ImGui::IsKeyPressed(ImGuiKey_F, false)) {
            assign_command(result.command, GUI_COMMAND_TOGGLE_FULLSCREEN,
                           g_context->applied);
        } else if (ImGui::IsKeyPressed(ImGuiKey_S, false)) {
            assign_command(result.command, GUI_COMMAND_SAVE_SNAPSHOT,
                           g_context->applied);
        }
    }
    /* Pointer-driven numeric controls update at a bounded cadence while held
       and always commit their last value on release.  Exact text entry and
       discrete controls retain the settled debounce, so an incomplete number
       can never trigger a live update and no Apply step is required anywhere. */
    const bool reactive_numeric_apply = dirty &&
        ((g_context->reactive_numeric_active &&
          now - g_context->last_reactive_apply_time >=
              kReactiveApplyIntervalSeconds) ||
         g_context->reactive_numeric_deactivated_after_edit);
    const bool settled_apply = dirty && !ImGui::IsAnyItemActive() &&
        now - g_context->last_physics_edit_time >= 0.30;
    if (g_context->automatic_apply_pending && !invalid &&
        !g_context->automatic_apply_blocked &&
        !g_context->apply_command_outstanding &&
        (reactive_numeric_apply || settled_apply) &&
        result.command.type == GUI_COMMAND_NONE) {
        /* Space is handled above in this same frame.  Snapshot the current
           live state instead of the pre-input rendering field in draft. */
        assign_config_apply(result.command, g_context->draft, live->running);
        if (reactive_numeric_apply) {
            g_context->last_reactive_apply_time = now;
        }
        g_context->automatic_apply_pending = false;
        g_context->apply_command_outstanding = true;
    }
    draw_runtime_status(*g_context, *live, result.command);
    draw_rebuild_toast(*g_context);
    result.wants_mouse = ImGui::GetIO().WantCaptureMouse;
    /* Docking/navigation may set WantCaptureKeyboard for an entire window,
       which previously swallowed R and Space.  Reserve global shortcuts only
       while the user is actually editing a widget or entering text. */
    result.wants_keyboard = ImGui::GetIO().WantTextInput ||
                            ImGui::IsAnyItemActive();
    rlImGuiEnd();
    return result;
}
