#include "gui_imgui_internal.hpp"
#if POLYGON_PARTICLES && !THREE_D
#include "gui_polygon_editor.hpp"
#endif

#include "imgui.h"

#if defined(PLATFORM_WEB)
#include <emscripten/emscripten.h>
#endif

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <string>

#if defined(PLATFORM_WEB)
extern "C" EMSCRIPTEN_KEEPALIVE void edmd_gui_web_restart_file_loaded(
    edmd_gui::Context *context, const char *virtual_path) {
    if (context == nullptr || virtual_path == nullptr || virtual_path[0] == '\0') {
        return;
    }
    std::snprintf(context->draft.initialization.load_path,
                  sizeof(context->draft.initialization.load_path), "%s",
                  virtual_path);
    context->last_physics_edit_time = ImGui::GetTime();
    context->automatic_apply_pending = true;
    context->automatic_apply_blocked = false;
    edmd_gui::refresh_validation(*context);
}

EM_JS(void, choose_browser_restart_file, (void *context), {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = '.dump,.lammpstrj,.txt,text/plain';
    input.style.display = 'none';
    input.addEventListener('change', async () => {
        try {
            const file = input.files && input.files[0];
            if (!file) return;
            const safeName = file.name.replace(/[^A-Za-z0-9._-]/g, '_') ||
                'restart.dump';
            const directory = '/uploads';
            try { FS.mkdir(directory); } catch (error) {
                if (!FS.analyzePath(directory).exists) throw error;
            }
            if (Module.edmdRestartUploadPath) {
                try { FS.unlink(Module.edmdRestartUploadPath); } catch (_) {}
            }
            const virtualPath = directory + '/' + Date.now() + '-' + safeName;
            const bytes = new Uint8Array(await file.arrayBuffer());
            FS.writeFile(virtualPath, bytes);
            Module.edmdRestartUploadPath = virtualPath;
            Module.ccall('edmd_gui_web_restart_file_loaded', null,
                ['number', 'string'], [context, virtualPath]);
        } catch (error) {
            console.error('Unable to import EDMD restart file', error);
        } finally {
            input.remove();
        }
    }, {once: true});
    document.body.appendChild(input);
    input.click();
});
#endif

namespace {

using edmd_gui::Context;

struct EnumOption {
    int value;
    const char *label;
    const char *disabled_reason;
};

bool contains_case_insensitive(const char *haystack, const std::string &needle) {
    if (needle.empty()) return true;
    if (haystack == nullptr) return false;
    std::string lower;
    lower.reserve(std::strlen(haystack));
    for (const unsigned char c : std::string(haystack)) {
        lower.push_back(static_cast<char>(std::tolower(c)));
    }
    return lower.find(needle) != std::string::npos;
}

bool begin_properties(const char *id) {
    return ImGui::BeginTable(id, 2,
        ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV |
        ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_PadOuterX);
}

void property_label(const char *label, const char *help) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(label);
    if (help != nullptr && help[0] != '\0') {
        ImGui::SameLine();
        edmd_gui::help_marker(help);
    }
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
}

void finish_disabled(bool disabled, const char *reason) {
    ImGui::EndDisabled();
    if (disabled) edmd_gui::disabled_reason(reason);
}

/* ImGui reports activity for the most recently submitted item, so capture it
   immediately after the slider/drag and before submitting the exact input in
   the same row.  The frame scheduler uses these flags to rebuild at a bounded
   cadence while the pointer is held and once more on release. */
void track_reactive_numeric_edit(Context &context, bool reactive) {
    if (!reactive) return;
    context.reactive_numeric_active |= ImGui::IsItemActive();
    context.reactive_numeric_deactivated_after_edit |=
        ImGui::IsItemDeactivatedAfterEdit();
}

bool checkbox_row(Context &context, const char *label, const char *keywords,
                  bool *value, const char *help = nullptr,
                  const char *reason = nullptr) {
    if (!edmd_gui::matches_filter(context, label, keywords)) return false;
    property_label(label, help);
    const bool disabled = reason != nullptr;
    ImGui::BeginDisabled(disabled);
    ImGui::PushID(label);
    const bool changed = ImGui::Checkbox("##value", value);
    ImGui::PopID();
    finish_disabled(disabled, reason);
    return changed;
}

bool input_int_row(Context &context, const char *label, const char *keywords,
                   int *value, const char *help = nullptr,
                   const char *reason = nullptr, int step = 1,
                   bool reactive = true) {
    if (!edmd_gui::matches_filter(context, label, keywords)) return false;
    property_label(label, help);
    const bool disabled = reason != nullptr;
    ImGui::BeginDisabled(disabled);
    ImGui::PushID(label);
    ImGui::BeginGroup();
    const float available = ImGui::GetContentRegionAvail().x;
    const float gap = ImGui::GetStyle().ItemInnerSpacing.x;
    ImGui::SetNextItemWidth(std::max(70.0f, available*0.56f));
    bool changed = ImGui::DragInt("##drag", value, (float)std::max(step, 1));
    track_reactive_numeric_edit(context, reactive);
    ImGui::SameLine(0.0f, gap);
    ImGui::SetNextItemWidth(-FLT_MIN);
    changed |= ImGui::InputInt("##input", value, step, step*10);
    ImGui::EndGroup();
    ImGui::PopID();
    finish_disabled(disabled, reason);
    return changed;
}

bool input_double_row(Context &context, const char *label, const char *keywords,
                      double *value, const char *help = nullptr,
                      const char *reason = nullptr,
                      const char *format = "%.10g",
                      bool reactive = true,
                      double drag_minimum = -DBL_MAX,
                      double drag_maximum = DBL_MAX) {
    if (!edmd_gui::matches_filter(context, label, keywords)) return false;
    property_label(label, help);
    const bool disabled = reason != nullptr;
    ImGui::BeginDisabled(disabled);
    ImGui::PushID(label);
    ImGui::BeginGroup();
    const float available = ImGui::GetContentRegionAvail().x;
    const float gap = ImGui::GetStyle().ItemInnerSpacing.x;
    const double magnitude = std::isfinite(*value) ? std::fabs(*value) : 0.0;
    const float speed = (float)std::clamp(magnitude*0.01, 1e-3, 1e6);
    ImGui::SetNextItemWidth(std::max(70.0f, available*0.56f));
    bool changed = ImGui::DragScalar("##drag", ImGuiDataType_Double, value,
                                     speed, &drag_minimum, &drag_maximum,
                                     format, ImGuiSliderFlags_AlwaysClamp);
    track_reactive_numeric_edit(context, reactive);
    ImGui::SameLine(0.0f, gap);
    ImGui::SetNextItemWidth(-FLT_MIN);
    changed |= ImGui::InputDouble("##input", value, 0.0, 0.0, format);
    const double sanitized = std::isfinite(*value)
        ? std::clamp(*value, drag_minimum, drag_maximum)
        : std::clamp(0.0, drag_minimum, drag_maximum);
    if (*value != sanitized) {
        *value = sanitized;
        changed = true;
    }
    ImGui::EndGroup();
    ImGui::PopID();
    finish_disabled(disabled, reason);
    return changed;
}

bool input_nonnegative_double_row(Context &context, const char *label,
                                  const char *keywords, double *value,
                                  const char *help = nullptr,
                                  const char *reason = nullptr,
                                  const char *format = "%.10g",
                                  bool reactive = true) {
    return input_double_row(context, label, keywords, value, help, reason,
                            format, reactive, 0.0, DBL_MAX);
}

bool input_bounded_double_row(Context &context, const char *label,
                              const char *keywords, double *value,
                              double minimum, double maximum,
                              const char *help = nullptr,
                              const char *reason = nullptr,
                              const char *format = "%.10g",
                              bool reactive = true) {
    return input_double_row(context, label, keywords, value, help, reason,
                            format, reactive, minimum, maximum);
}

/* A bounded slider and an exact text entry deliberately share the same row.
   The slider is convenient for exploration; the input keeps uncommon or
   high-precision values possible.  Structural slider values are live-reloaded
   at a bounded cadence while held and committed once more on release. */
bool slider_int_row(Context &context, const char *label, const char *keywords,
                    int *value, int minimum, int maximum,
                    const char *help = nullptr,
                    const char *reason = nullptr, int step = 1,
                    bool reactive = true) {
    if (!edmd_gui::matches_filter(context, label, keywords)) return false;
    property_label(label, help);
    const bool disabled = reason != nullptr;
    bool changed = false;
    ImGui::BeginDisabled(disabled);
    ImGui::PushID(label);
    ImGui::BeginGroup();
    const float available = ImGui::GetContentRegionAvail().x;
    const float gap = ImGui::GetStyle().ItemInnerSpacing.x;
    ImGui::SetNextItemWidth(std::max(70.0f, available*0.56f));
    changed |= ImGui::SliderInt("##slider", value, minimum, maximum);
    track_reactive_numeric_edit(context, reactive);
    ImGui::SameLine(0.0f, gap);
    ImGui::SetNextItemWidth(-FLT_MIN);
    changed |= ImGui::InputInt("##input", value, step, step*10);
    const int sanitized = std::clamp(*value, minimum, maximum);
    if (*value != sanitized) {
        *value = sanitized;
        changed = true;
    }
    ImGui::EndGroup();
    ImGui::PopID();
    finish_disabled(disabled, reason);
    return changed;
}

bool slider_double_row(Context &context, const char *label,
                       const char *keywords, double *value,
                       double minimum, double maximum,
                       const char *help = nullptr,
                       const char *reason = nullptr,
                       const char *format = "%.6g",
                       bool logarithmic = false,
                       bool reactive = true) {
    if (!edmd_gui::matches_filter(context, label, keywords)) return false;
    property_label(label, help);
    const bool disabled = reason != nullptr;
    bool changed = false;
    ImGui::BeginDisabled(disabled);
    ImGui::PushID(label);
    ImGui::BeginGroup();
    const float available = ImGui::GetContentRegionAvail().x;
    const float gap = ImGui::GetStyle().ItemInnerSpacing.x;
    ImGui::SetNextItemWidth(std::max(70.0f, available*0.56f));
    const ImGuiSliderFlags flags = logarithmic
        ? ImGuiSliderFlags_Logarithmic : ImGuiSliderFlags_None;
    changed |= ImGui::SliderScalar("##slider", ImGuiDataType_Double, value,
                                   &minimum, &maximum, format, flags);
    track_reactive_numeric_edit(context, reactive);
    ImGui::SameLine(0.0f, gap);
    ImGui::SetNextItemWidth(-FLT_MIN);
    changed |= ImGui::InputDouble("##input", value, 0.0, 0.0, "%.10g");
    const double sanitized = std::isfinite(*value)
        ? std::clamp(*value, minimum, maximum)
        : std::clamp(0.0, minimum, maximum);
    if (*value != sanitized) {
        *value = sanitized;
        changed = true;
    }
    ImGui::EndGroup();
    ImGui::PopID();
    finish_disabled(disabled, reason);
    return changed;
}

#if !defined(PLATFORM_WEB)
bool input_text_row(Context &context, const char *label, const char *keywords,
                    char *value, std::size_t capacity,
                    const char *help = nullptr,
                    const char *reason = nullptr) {
    if (!edmd_gui::matches_filter(context, label, keywords)) return false;
    property_label(label, help);
    const bool disabled = reason != nullptr;
    ImGui::BeginDisabled(disabled);
    ImGui::PushID(label);
    const bool changed = ImGui::InputText("##value", value, capacity);
    ImGui::PopID();
    finish_disabled(disabled, reason);
    return changed;
}
#endif

#if defined(PLATFORM_WEB)
void browser_restart_file_row(Context &context, char *value,
                              std::size_t capacity) {
    if (!edmd_gui::matches_filter(context, "Restart file", "load path dump")) {
        return;
    }
    property_label("Restart file",
                   "Choose a local LAMMPS-style dump. It is copied into the browser's temporary filesystem and the last frame is used.");
    ImGui::PushID("Restart file");
    const float available = ImGui::GetContentRegionAvail().x;
    const float button_width = std::min(150.0f, available*0.42f);
    if (ImGui::Button("Choose local file", ImVec2(button_width, 0))) {
        choose_browser_restart_file(&context);
    }
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputText("##uploaded-path", value, capacity,
                     ImGuiInputTextFlags_ReadOnly);
    ImGui::PopID();
}
#endif

bool enum_row(Context &context, const char *label, const char *keywords,
              int *value, const EnumOption *options, std::size_t count,
              const char *help = nullptr, const char *reason = nullptr) {
    if (!edmd_gui::matches_filter(context, label, keywords)) return false;
    property_label(label, help);
    const char *preview = "Unknown";
    for (std::size_t i = 0; i < count; ++i) {
        if (options[i].value == *value) preview = options[i].label;
    }
    const bool disabled = reason != nullptr;
    bool changed = false;
    ImGui::BeginDisabled(disabled);
    ImGui::PushID(label);
    if (ImGui::BeginCombo("##value", preview)) {
        for (std::size_t i = 0; i < count; ++i) {
            const bool option_disabled = options[i].disabled_reason != nullptr;
            ImGui::BeginDisabled(option_disabled);
            const bool selected = options[i].value == *value;
            if (ImGui::Selectable(options[i].label, selected) && !option_disabled) {
                *value = options[i].value;
                changed = true;
            }
            if (selected) ImGui::SetItemDefaultFocus();
            if (option_disabled) {
                edmd_gui::disabled_reason(options[i].disabled_reason);
            }
            ImGui::EndDisabled();
        }
        ImGui::EndCombo();
    }
    ImGui::PopID();
    finish_disabled(disabled, reason);
    return changed;
}

void section_note(const char *text) {
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s", text);
    ImGui::PopTextWrapPos();
}

bool collapsing_header(const Context &context, const char *label,
                       ImGuiTreeNodeFlags flags = 0) {
    /* A text filter should reveal matches below nested sections instead of
       requiring the user to guess which closed heading contains them. */
    if (context.search[0] != '\0') {
        ImGui::SetNextItemOpen(true, ImGuiCond_Always);
    }
    return ImGui::CollapsingHeader(label, flags);
}

void draw_first_error(const Context &context) {
    for (std::size_t i = 0; i < context.stored_issue_count; ++i) {
        if (context.issues[i].severity != CONFIG_ISSUE_ERROR) continue;
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(ImVec4(0.95f, 0.35f, 0.30f, 1.0f),
                           "First error: %s", context.issues[i].message);
        ImGui::PopTextWrapPos();
        return;
    }
}

void clamp_live_double(double &value, double minimum, double maximum,
                       double fallback) {
    if (!std::isfinite(value)) value = fallback;
    value = std::clamp(value, minimum, maximum);
}

void disable_polygon_incompatible_saved_analysis(SimulationConfig &config) {
    SimObservablesConfig &observables = config.observables;
    /* These legacy file-oriented measurements either receive impulses only
       from the disk handler or normalize local occupancy with pi*r^2. Live
       center-based plots use a separate pipeline and remain available. */
    observables.pressure = false;
    observables.pressure_field = false;
    observables.temperature_profile = false;
    observables.interface_density = false;
    observables.interface_temperature = false;
    observables.osmosis = false;
    observables.voronoi_area = false;
    observables.critical_regions = false;
    observables.bubble_stop = false;
    if (config.geometry.center_membrane ||
        config.initialization.mode == SIM_INIT_MEMBRANE_OSMOSIS) {
        observables.thermo_enabled = false;
    }
    config.output.periodic_dumps = false;
}

void draw_setup(Context &context) {
    SimInitializationConfig &init = context.draft.initialization;
    const bool is3d = context.capabilities.dimensions == 3;
    if (is3d) {
        /* Recover old files/CLI-derived drafts before hiding their planar-only
           controls.  The GUI must never leave an invalid value with no route
           back to a supported configuration. */
        if (init.mode == SIM_INIT_SQUARE_LATTICE ||
            init.mode == SIM_INIT_MEMBRANE_OSMOSIS) {
            init.mode = SIM_INIT_RANDOM_GROWTH;
            context.draft.geometry.center_membrane = false;
        }
        if (init.box_length_mode != SIM_BOX_FROM_PACKING_FRACTION) {
            init.box_length_mode = SIM_BOX_FROM_PACKING_FRACTION;
        }
    }
    constexpr double triangular_close_packing =
        0.90689968211710892530;
    constexpr double hcp_close_packing =
        0.74048048969306104117;
    constexpr double pi = 3.141592653589793238462643383279502884;
    const bool curved_boundary = context.draft.geometry.boundary ==
                                 SIM_BOUNDARY_CIRCULAR_SPHERICAL;
    const bool file_initializer = init.mode == SIM_INIT_DUMP_FILE;
    const bool generic_initializer = init.mode == SIM_INIT_RANDOM_GROWTH;
    const bool polygon_mode = !is3d &&
        context.draft.particle_shape.mode == SIM_PARTICLE_CONVEX_POLYGONS;
    const double polygon_unit_area = polygon_mode
        ? simulation_particle_shape_unit_area(&context.draft.particle_shape)
        : pi;
    /* The explicit hex initializer still uses a disk-safe center lattice.
       Random preparation grows the actual polygons and has no envelope cap. */
    const double polygon_lattice_limit = polygon_mode &&
        std::isfinite(polygon_unit_area) && polygon_unit_area > 0.0
        ? polygon_unit_area/(2.0*std::sqrt(3.0))
        : triangular_close_packing;
    const double hex_close_packing = polygon_mode
        ? polygon_lattice_limit
        : (is3d ? hcp_close_packing : triangular_close_packing);
    if (polygon_mode)
        disable_polygon_incompatible_saved_analysis(context.draft);
    const char *hex_initializer_reason =
        is3d && !context.capabilities.lattice_initializers_3d
        ? "This build does not provide the complete 3D HCP initializer."
        : nullptr;
    const char *curved_initializer_reason = curved_boundary
        ? "Circular/spherical confinement currently supports only Random placement + growth."
        : nullptr;

    #if POLYGON_PARTICLES && !THREE_D
    if (!is3d && context.capabilities.convex_polygon_particles) {
        SimParticleShapeConfig &particle_shape = context.draft.particle_shape;
        ImGui::SeparatorText("Particle body");
        if (begin_properties("setup-particle-body")) {
            const char *polygon_build_reason =
                (context.capabilities.stress || context.capabilities.heat)
                    ? "This build accumulates transport terms that do not yet include polygon rotation."
                : (context.capabilities.intruder
                    ? "The legacy intruder collision path is compiled into this build."
                : (context.capabilities.rough_particles
                    ? "The legacy TANGENTIAL rough-disk model is compiled into this build."
                : (!context.capabilities.additive_contacts
                    ? "Convex polygons require additive radius scaling (NONADDITIVE=1)."
                    : nullptr)));
            const EnumOption body_modes[] = {
                {SIM_PARTICLE_CIRCLES, "Circular disks", nullptr},
                {SIM_PARTICLE_CONVEX_POLYGONS,
                    "Rigid convex polygons", polygon_build_reason},
            };
            int body_mode = particle_shape.mode;
            if (enum_row(context, "Shape model", "circle disk polygon rigid body",
                         &body_mode, body_modes, std::size(body_modes),
                         "Polygon particles use the same editable convex body shape, scaled by each particle radius.")) {
                particle_shape.mode =
                    static_cast<SimParticleShapeMode>(body_mode);
                if (particle_shape.mode == SIM_PARTICLE_CONVEX_POLYGONS) {
                    /* Enter the model through a complete, valid combination.
                       The disabled controls in their own panels explain why
                       these choices cannot be re-enabled for polygons. */
                    if (init.mode == SIM_INIT_DUMP_FILE)
                        init.mode = SIM_INIT_RANDOM_GROWTH;
                    context.draft.dynamics.flight = SIM_FLIGHT_BALLISTIC;
                    context.draft.dynamics.injection = SIM_INJECTION_NONE;
                    if (context.draft.driving.thermostat !=
                            SIM_THERMOSTAT_LANGEVIN)
                        context.draft.driving.thermostat =
                            SIM_THERMOSTAT_NONE;
                    context.draft.driving.euler_noise = false;
                    context.draft.potential.enabled = false;
                    context.draft.geometry.thermal_wall = false;
                    context.draft.geometry.shear_wall = false;
                    context.draft.intruder.enabled = false;
                    disable_polygon_incompatible_saved_analysis(context.draft);
                }
            }
            if (particle_shape.mode == SIM_PARTICLE_CONVEX_POLYGONS) {
                checkbox_row(context, "Random initial orientation",
                             "polygon angle orientation",
                             &particle_shape.random_orientation,
                             "Otherwise every body starts with the editor's orientation.");
                slider_double_row(context, "Initial angular temperature",
                                  "polygon spin omega rotation",
                                  &particle_shape.angular_temperature,
                                  0.0, 10.0,
                                  "Sets the variance of initial angular velocity; zero starts without spin.");
            }
            ImGui::EndTable();
        }
        if (particle_shape.mode == SIM_PARTICLE_CONVEX_POLYGONS) {
            section_note(
                "Drag a numbered vertex or an edge handle. The editor stops at the convexity boundary; releasing commits one restart.");
            (void)edmd_gui::draw_polygon_editor(particle_shape);
            section_note(
                "Radius is the circumradius used only by the collision grid. Packing fraction is the polygon area occupied by all particles divided by the container area; random preparation grows and collides the actual polygons.");
            section_note(
                "Live energy, translational temperature, pressure, MSD, g(r), and structure factor remain available. Legacy saved local pressure, packing, and heat profiles are disabled because they still assume disks.");
        }
    }
    #endif

    ImGui::SeparatorText("Starting point");
    if (begin_properties("setup-start")) {
        const SimInitMode previous_mode = init.mode;
        int mode = init.mode;
        bool mode_changed = false;
        if (is3d) {
            /* Square-lattice and membrane initializers are intrinsically 2D;
               omitting them is clearer than presenting dead menu entries. */
            const EnumOption modes_3d[] = {
                {SIM_INIT_RANDOM_GROWTH, "Random placement + growth", nullptr},
                {SIM_INIT_DUMP_FILE, "LAMMPS-style dump restart",
                    curved_initializer_reason},
                {SIM_INIT_HEX_LATTICE,
                    "Hexagonal close-packed (HCP) lattice",
                    curved_initializer_reason != nullptr
                        ? curved_initializer_reason : hex_initializer_reason},
            };
            mode_changed = enum_row(
                context, "Initial layout", "start source load lattice HCP",
                &mode, modes_3d, std::size(modes_3d),
                "Changing the initializer automatically constructs a new t = 0 state.");
        } else {
            const EnumOption modes_2d[] = {
                {SIM_INIT_RANDOM_GROWTH,
                    polygon_mode ? "Random placement + polygon growth"
                                 : "Random placement + growth", nullptr},
                {SIM_INIT_DUMP_FILE, "LAMMPS-style dump restart",
                    polygon_mode
                        ? "Legacy dumps do not contain polygon angles or angular velocities."
                        : curved_initializer_reason},
                {SIM_INIT_SQUARE_LATTICE, "Binary square lattice",
                    curved_initializer_reason},
                {SIM_INIT_HEX_LATTICE, "Hexagonal lattice",
                    curved_initializer_reason},
                {SIM_INIT_MEMBRANE_OSMOSIS, "Membrane / osmosis",
                    curved_initializer_reason},
            };
            mode_changed = enum_row(
                context, "Initial layout",
                "start source load lattice membrane osmosis",
                &mode, modes_2d, std::size(modes_2d),
                "Changing the initializer automatically constructs a new t = 0 state.");
        }
        if (mode_changed) {
            init.mode = static_cast<SimInitMode>(mode);
            if (init.mode == SIM_INIT_SQUARE_LATTICE ||
                init.mode == SIM_INIT_MEMBRANE_OSMOSIS) {
                init.composition = SIM_COMPOSITION_BINARY;
            } else if (init.mode == SIM_INIT_HEX_LATTICE) {
                init.composition = SIM_COMPOSITION_MONODISPERSE;
                init.packing_fraction = std::min(init.packing_fraction,
                                                 hex_close_packing);
            }
            if (init.mode == SIM_INIT_RANDOM_GROWTH &&
                init.composition == SIM_COMPOSITION_BINARY &&
                (init.small_fraction <= 0.0 || init.small_fraction >= 1.0)) {
                init.small_fraction = 0.5;
            }
            if (init.mode == SIM_INIT_MEMBRANE_OSMOSIS) {
                context.draft.geometry.center_membrane = true;
            } else if (previous_mode == SIM_INIT_MEMBRANE_OSMOSIS) {
                /* The membrane initializer is the only supported owner of this
                   legacy separator.  Do not leave a hidden incompatible wall
                   behind when the starting layout changes. */
                context.draft.geometry.center_membrane = false;
            }
        }
        if (init.mode == SIM_INIT_DUMP_FILE) {
			#if defined(PLATFORM_WEB)
			browser_restart_file_row(context, init.load_path,
			                         sizeof(init.load_path));
			#else
            input_text_row(context, "Restart file", "load path dump",
                           init.load_path, sizeof(init.load_path),
                           "The last frame is used for dump files.");
			#endif
        }
        input_int_row(context, "Random seed", "seed reproducibility",
                      &init.random_seed,
                      "Seeds particle placement and stochastic dynamics on the next automatic reload.");
        ImGui::EndTable();
    }

    ImGui::SeparatorText("System size and composition");
    if (begin_properties("setup-system")) {
        const char *file_count_reason = file_initializer
            ? "Particle count is read from the restart file." : nullptr;
        slider_int_row(context, "Particle count", "N population number",
                      &init.particle_count,
                      1, 100000,
                      "The realized count can differ for lattice initializers that require complete rows.",
                      file_count_reason, 10);
        const char *packing_reason = file_initializer
            ? "Packing fraction is derived from the loaded particle and box data."
            : (init.mode == SIM_INIT_MEMBRANE_OSMOSIS
                ? "This initializer derives the overall packing fraction from its liquid and gas inputs."
                : nullptr);
        const double packing_slider_max = init.mode == SIM_INIT_HEX_LATTICE
            ? hex_close_packing : 1.0;
        const char *packing_label = init.mode == SIM_INIT_RANDOM_GROWTH
            ? "Target packing fraction" : "Packing fraction";
        slider_double_row(context, packing_label, "phi density volume fraction",
                         &init.packing_fraction, 0.001, packing_slider_max,
                         polygon_mode && init.mode == SIM_INIT_RANDOM_GROWTH
                            ? "Target occupied polygon area divided by container area. Diagnostics shows the smaller instantaneous value while polygon-native preparation grows toward it; difficult near-unity targets can jam and pause clearly."
                         : polygon_mode
                            ? "Occupied polygon area divided by container area. This explicit lattice remains capped by its conservative center spacing."
                         : init.mode == SIM_INIT_HEX_LATTICE
                            ? (is3d
                                ? "HCP is capped at exact hard-sphere close packing; larger values overlap at t = 0."
                                : "The triangular lattice is capped at exact hard-disk close packing; larger values overlap at t = 0.")
                            : "Values near one may be impossible to initialize without overlaps; validation warnings remain visible.",
                         packing_reason);
        const char *initial_energy_label = polygon_mode
            ? "Initial translational energy / particle"
            : "Initial energy / particle";
        slider_double_row(context, initial_energy_label,
                         "Einit temperature velocity",
                         &init.initial_energy, 0.0, 10.0,
                         polygon_mode
                            ? "Sets translational kinetic energy after preparation. Rotational energy is initialized separately by Initial angular temperature and is included in the displayed total energy."
                            : nullptr,
                         init.mode == SIM_INIT_DUMP_FILE
                            ? "A dump restart restores particle velocities, so this energy is not used."
                            : nullptr);

        const char *composition_reason = nullptr;
        if (file_initializer) {
            composition_reason = "Species, radii, and masses are read from the restart file.";
        } else if (init.mode == SIM_INIT_SQUARE_LATTICE) {
            composition_reason = "This initializer always constructs a binary square lattice.";
        } else if (init.mode == SIM_INIT_HEX_LATTICE) {
            composition_reason = "This initializer constructs monodisperse particles.";
        } else if (init.mode == SIM_INIT_MEMBRANE_OSMOSIS) {
            composition_reason = "The membrane initializer derives its binary composition from the liquid and gas densities.";
        }
        const EnumOption compositions[] = {
            {SIM_COMPOSITION_MONODISPERSE, "Monodisperse", nullptr},
            {SIM_COMPOSITION_BINARY, "Binary mixture", nullptr},
            {SIM_COMPOSITION_POLYDISPERSE, "Continuous polydispersity", nullptr},
        };
        int composition = init.composition;
        if (enum_row(context, "Composition", "mixture bidisperse polydisperse",
                     &composition, compositions, std::size(compositions), nullptr,
                     composition_reason)) {
            init.composition = static_cast<SimCompositionMode>(composition);
            if (init.composition == SIM_COMPOSITION_BINARY &&
                (init.small_fraction <= 0.0 || init.small_fraction >= 1.0)) {
                /* The legacy globals infer binary composition from this
                   fraction.  Seed a real two-species mixture so the applied
                   configuration cannot collapse back to monodisperse. */
                init.small_fraction = 0.5;
            }
        }
        const bool forced_size_ratio =
            init.mode == SIM_INIT_SQUARE_LATTICE ||
            init.mode == SIM_INIT_MEMBRANE_OSMOSIS;
        if (init.composition == SIM_COMPOSITION_BINARY || forced_size_ratio) {
            const char *size_ratio_reason = forced_size_ratio ? nullptr
                : composition_reason;
            slider_double_row(context, "Small / large radius", "size ratio sizeratio",
                             &init.size_ratio, 0.05, 1.0,
                             "Radius of the small species divided by the unit-radius large species; valid range (0, 1].",
                             size_ratio_reason);
            const char *fraction_reason = composition_reason;
            if (init.mode == SIM_INIT_SQUARE_LATTICE) {
                fraction_reason =
                    "This binary lattice constructs equal complete sublattices; the fraction is fixed by realized row counts.";
            } else if (init.mode == SIM_INIT_MEMBRANE_OSMOSIS) {
                fraction_reason =
                    "The membrane initializer derives the small-particle fraction from liquid/gas densities and the size ratio.";
            }
            slider_double_row(context, "Small-particle fraction", "small fraction composition",
                             &init.small_fraction, 0.01, 0.99,
                             "Both species remain present; use Monodisperse when only one radius is desired.",
                             fraction_reason);
        } else if (init.composition == SIM_COMPOSITION_POLYDISPERSE) {
            input_nonnegative_double_row(context, "Distribution shape a", "polydisperse beta distribution",
                                         &init.polydisperse_a, nullptr, composition_reason);
            input_nonnegative_double_row(context, "Distribution shape b", "polydisperse beta distribution",
                                         &init.polydisperse_b, nullptr, composition_reason);
        }
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Box sizing");
    if (begin_properties("setup-box")) {
        const char *derived_box_reason = file_initializer
            ? "Box lengths are read from the restart file."
            : (!generic_initializer
                ? "This initializer derives its box from its own geometry."
                : nullptr);
        if (!is3d) {
            const EnumOption box_modes_2d[] = {
                {SIM_BOX_FROM_PACKING_FRACTION,
                    "Derive from packing fraction", nullptr},
                {SIM_BOX_FIXED_LX, "Fix Lx",
                    curved_boundary
                        ? "Circular confinement derives its box from packing fraction."
                        : derived_box_reason},
                {SIM_BOX_FIXED_LY, "Fix Ly",
                    curved_boundary
                        ? "Circular confinement derives its box from packing fraction."
                        : derived_box_reason},
            };
            int box_mode = init.box_length_mode;
            if (enum_row(context, "Length constraint",
                         "box Lx Ly geometry aspect", &box_mode,
                         box_modes_2d, std::size(box_modes_2d))) {
                init.box_length_mode = static_cast<SimBoxLengthMode>(box_mode);
            }
        }
        const char *aspect_reason = file_initializer
            ? "Aspect ratio is read from the restart box."
            : (curved_boundary
                ? "Circular/spherical confinement forces equal box lengths."
                : nullptr);
        slider_double_row(context, "Aspect ratio Ly / Lx", "aspect box shape",
                         &init.aspect_ratio, 0.1, 10.0, nullptr, aspect_reason,
                         "%.5g", true);
        if (!is3d && init.box_length_mode == SIM_BOX_FIXED_LX) {
            input_nonnegative_double_row(context, "Fixed Lx", "box length", &init.box_lx,
                                         nullptr,
                             curved_boundary ? "Circular/spherical confinement derives its diameter from packing fraction."
                             : derived_box_reason);
        } else if (!is3d && init.box_length_mode == SIM_BOX_FIXED_LY) {
            input_nonnegative_double_row(context, "Fixed Ly", "box length", &init.box_ly,
                                         nullptr,
                             curved_boundary ? "Circular/spherical confinement derives its diameter from packing fraction."
                             : derived_box_reason);
        }
        ImGui::EndTable();
    }

    if (init.mode == SIM_INIT_RANDOM_GROWTH ||
        init.mode == SIM_INIT_MEMBRANE_OSMOSIS) {
        ImGui::SeparatorText("Growth preparation");
        if (begin_properties("setup-growth")) {
            checkbox_row(context, "Automatic growth rate", "growth LS compression",
                         &init.automatic_growth_rate,
                         "Let EDMD select a conservative rate for the requested density.");
            if (!init.automatic_growth_rate) {
                input_nonnegative_double_row(context, "Growth rate", "growth vr compression",
                                             &init.growth_rate);
            }
            ImGui::EndTable();
        }
    }

    const bool phase_initializer = init.mode == SIM_INIT_MEMBRANE_OSMOSIS;
    if (phase_initializer || init.mode == SIM_INIT_HEX_LATTICE) {
        ImGui::SeparatorText("Initializer-specific parameters");
        if (begin_properties("setup-special")) {
            if (phase_initializer) {
                input_nonnegative_double_row(context, "Liquid packing fraction", "membrane dense phil",
                                             &init.liquid_packing_fraction);
                input_nonnegative_double_row(context, "Gas packing fraction", "membrane dilute phig",
                                             &init.gas_packing_fraction);
            }
            if (init.mode == SIM_INIT_HEX_LATTICE) {
                input_nonnegative_double_row(context, "Central hole radius", "hex hole annulus",
                                             &init.hex_hole_radius);
            }
            ImGui::EndTable();
        }
    }
}

void draw_dynamics(Context &context) {
    SimDynamicsConfig &dynamics = context.draft.dynamics;
    SimDrivingConfig &driving = context.draft.driving;
    SimPotentialConfig &potential = context.draft.potential;
    const bool is3d = context.capabilities.dimensions == 3;
    const bool polygon_mode = !is3d &&
        context.draft.particle_shape.mode == SIM_PARTICLE_CONVEX_POLYGONS;
    const char *polygon_ballistic_reason = polygon_mode
        ? "Rigid convex polygons currently use conservative collision prediction with ballistic free flight."
        : nullptr;
    const char *transport_flight_reason =
        context.capabilities.stress || context.capabilities.heat
        ? "STRESS/HEAT transport bookkeeping currently requires ballistic free flight."
        : nullptr;
    const char *nonballistic_flight_reason = polygon_ballistic_reason != nullptr
        ? polygon_ballistic_reason : (context.draft.intruder.enabled
        ? "The polygonal intruder predictor currently requires ballistic particle flight."
        : transport_flight_reason);
    const bool physical_boundary =
        context.draft.geometry.boundary != SIM_BOUNDARY_PERIODIC ||
        context.draft.geometry.center_membrane ||
        context.draft.initialization.mode == SIM_INIT_MEMBRANE_OSMOSIS;
    /* STRESS/HEAT has complete Lees-Edwards bookkeeping.  The transport
       restriction applies only to field and damping trajectories. */
    const char *lees_reason = polygon_ballistic_reason != nullptr
        ? polygon_ballistic_reason : (context.draft.intruder.enabled
        ? "The polygonal intruder predictor currently requires ballistic particle flight."
        : nullptr);
    if (lees_reason == nullptr && is3d && !context.capabilities.shear_3d) {
        lees_reason = "Lees-Edwards shear has not been validated for the 3D build.";
    }
    if (lees_reason == nullptr && physical_boundary) {
        lees_reason = "Lees-Edwards shear requires a fully periodic domain.";
    }
    if (lees_reason == nullptr && potential.enabled) {
        lees_reason = "Lees-Edwards shear cannot be combined with square-well events.";
    }

    ImGui::SeparatorText("Free flight");
    if (begin_properties("dynamics-flight")) {
        const EnumOption modes[] = {
            {SIM_FLIGHT_BALLISTIC, "Ballistic", nullptr},
            {SIM_FLIGHT_CONSTANT_FIELD, "Constant field", nonballistic_flight_reason},
            {SIM_FLIGHT_LINEAR_DAMPING, "Linear damping", nonballistic_flight_reason},
            {SIM_FLIGHT_LEES_EDWARDS_SHEAR, "Lees-Edwards shear", lees_reason},
        };
        int flight = dynamics.flight;
        if (enum_row(context, "Trajectory law", "flight damping field shear ballistic",
                     &flight, modes, std::size(modes),
                     "One analytical free-flight law is active at a time, preventing field and damping from silently shadowing each other.")) {
            dynamics.flight = static_cast<SimFlightMode>(flight);
            if (dynamics.flight == SIM_FLIGHT_LINEAR_DAMPING) {
                driving.gamma = dynamics.damping_gamma;
            }
        }
        if (dynamics.flight == SIM_FLIGHT_CONSTANT_FIELD) {
            input_double_row(context, "Field acceleration", "gravity field g",
                             &dynamics.constant_field,
                             "The current analytical solver requires a strictly negative y acceleration.");
        } else if (dynamics.flight == SIM_FLIGHT_LINEAR_DAMPING) {
            if (input_nonnegative_double_row(context, "Damping gamma", "drag flight gamma",
                                             &dynamics.damping_gamma,
                                             "The legacy engine shares gamma with Langevin thermostats; both schema values are kept synchronized.")) {
                driving.gamma = dynamics.damping_gamma;
            }
        } else if (dynamics.flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR) {
            input_double_row(context, "Shear strain rate", "Lees Edwards rate",
                             &dynamics.shear_strain_rate);
        }
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Contact restitution");
    if (begin_properties("dynamics-restitution")) {
        const char *restitution_toggle_reason =
            context.capabilities.runtime_exponential_restitution ? nullptr
            : "The constant/exponential restitution choice is fixed at compile time in this build.";
        const EnumOption restitution[] = {
            {SIM_RESTITUTION_CONSTANT, "Constant coefficient", nullptr},
            {SIM_RESTITUTION_EXPONENTIAL, "Speed-dependent exponential", nullptr},
        };
        int restitution_mode = dynamics.restitution_model;
        if (enum_row(context, "Restitution model", "collision dissipation exponential",
                     &restitution_mode, restitution, std::size(restitution),
                     nullptr, restitution_toggle_reason)) {
            dynamics.restitution_model =
                static_cast<SimRestitutionMode>(restitution_mode);
        }
        slider_double_row(context, "Particle restitution", "res alpha collision",
                         &dynamics.particle_restitution, 0.0, 1.0, nullptr,
                         dynamics.restitution_model == SIM_RESTITUTION_EXPONENTIAL
                            ? "The exponential restitution law uses its velocity scale and amplitude instead of a constant coefficient."
                            : nullptr);
        if (dynamics.restitution_model == SIM_RESTITUTION_EXPONENTIAL) {
            input_nonnegative_double_row(context, "Velocity scale", "exponential vo restitution",
                                         &dynamics.exponential_velocity_scale);
            input_nonnegative_double_row(context, "Exponential amplitude", "ao restitution",
                                         &dynamics.exponential_amplitude);
        }
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Collision energy input");
    if (begin_properties("dynamics-injection")) {
        const char *basic_injection_reason =
            context.capabilities.runtime_basic_injection_toggles ? nullptr
            : "Basic collision-injection laws are fixed at compile time in this build.";
        const char *advanced_injection_reason =
            context.capabilities.runtime_advanced_injection_toggles ? nullptr
            : "This advanced injection law is fixed at compile time and cannot be selected at restart.";
        const char *rough_reason = !context.capabilities.rough_particles
            ? "Requires a build with the rotational rough-particle model."
            : nullptr;
        const char *rough_fixed_reason = context.capabilities.rough_particles
            ? "Rough-particle collision dynamics are selected by the TANGENTIAL build option and cannot be changed at runtime."
            : nullptr;
        const bool current_advanced_injection =
            dynamics.injection == SIM_INJECTION_EVOLVING_DELTA ||
            dynamics.injection == SIM_INJECTION_COLLISION_ENERGY ||
            dynamics.injection == SIM_INJECTION_ACTIVE_TANGENTIAL;
        const char *injection_selector_reason = polygon_mode
            ? "Collision-injection laws are not defined for rigid convex contacts."
            : (rough_fixed_reason != nullptr
            ? rough_fixed_reason
            : (current_advanced_injection &&
               !context.capabilities.runtime_advanced_injection_toggles
                ? advanced_injection_reason
                : (!current_advanced_injection &&
                   !context.capabilities.runtime_basic_injection_toggles
                    ? basic_injection_reason : nullptr)));
        int injection = dynamics.injection;
        bool injection_changed = false;
        if (is3d) {
            const char *tangent_reason =
                !context.capabilities.tangent_delta_3d
                    ? "Fixed-axis tangential delta is unavailable in this 3D build."
                    : basic_injection_reason;
            const EnumOption injections_3d[] = {
                {SIM_INJECTION_NONE, "None", basic_injection_reason},
                {SIM_INJECTION_DELTA_NORMAL, "Normal delta", basic_injection_reason},
                {SIM_INJECTION_DELTA_TANGENT,
                    "Tangential delta (u_z x n)", tangent_reason},
                {SIM_INJECTION_DOUBLE_DELTA_NORMAL,
                    "Two-scale normal delta", basic_injection_reason},
                {SIM_INJECTION_DOUBLE_DELTA_TANGENT,
                    "Two-scale tangential delta (u_z x n)", tangent_reason},
                {SIM_INJECTION_EVOLVING_DELTA,
                    "Collision-age evolving delta", advanced_injection_reason},
                {SIM_INJECTION_COLLISION_ENERGY,
                    "Prescribed collision energy", advanced_injection_reason},
                {SIM_INJECTION_ACTIVE_TANGENTIAL,
                    "Relative tangential activation", advanced_injection_reason},
            };
            injection_changed = enum_row(
                context, "Injection law", "delta drive active tangent energy",
                &injection, injections_3d, std::size(injections_3d),
                "Restitution and collision injection are independent, so exponential restitution may be combined with a delta law.",
                injection_selector_reason);
        } else {
            const EnumOption injections_2d[] = {
                {SIM_INJECTION_NONE, "None", basic_injection_reason},
                {SIM_INJECTION_DELTA_NORMAL, "Normal delta", basic_injection_reason},
                {SIM_INJECTION_DELTA_TANGENT,
                    "Tangential delta", basic_injection_reason},
                {SIM_INJECTION_DOUBLE_DELTA_NORMAL,
                    "Two-scale normal delta", basic_injection_reason},
                {SIM_INJECTION_DOUBLE_DELTA_TANGENT,
                    "Two-scale tangential delta", basic_injection_reason},
                {SIM_INJECTION_EVOLVING_DELTA,
                    "Collision-age evolving delta", advanced_injection_reason},
                {SIM_INJECTION_COLLISION_ENERGY,
                    "Prescribed collision energy", advanced_injection_reason},
                {SIM_INJECTION_ACTIVE_TANGENTIAL,
                    "Relative tangential activation", advanced_injection_reason},
                {SIM_INJECTION_ROUGH_ROTATIONAL,
                    "Rough rotational particles", rough_reason},
            };
            injection_changed = enum_row(
                context, "Injection law", "delta drive active tangent energy",
                &injection, injections_2d, std::size(injections_2d),
                "Restitution and collision injection are independent, so exponential restitution may be combined with a delta law.",
                injection_selector_reason);
        }
        if (injection_changed) {
            dynamics.injection = static_cast<SimInjectionMode>(injection);
        }
        switch (dynamics.injection) {
            case SIM_INJECTION_DELTA_NORMAL:
            case SIM_INJECTION_DELTA_TANGENT:
            case SIM_INJECTION_ACTIVE_TANGENTIAL:
                input_nonnegative_double_row(context, "Delta", "collision kick injection",
                                             &dynamics.delta);
                break;
            case SIM_INJECTION_DOUBLE_DELTA_NORMAL:
            case SIM_INJECTION_DOUBLE_DELTA_TANGENT:
                input_nonnegative_double_row(context, "Small delta", "deltam two scale",
                                             &dynamics.delta_small);
                input_nonnegative_double_row(context, "Large delta", "deltaM two scale",
                                             &dynamics.delta_large);
                input_nonnegative_double_row(context, "Collision memory time", "ts memory",
                                             &dynamics.collision_memory_time);
                break;
            case SIM_INJECTION_EVOLVING_DELTA:
                input_nonnegative_double_row(context, "Evolution time tau", "delta age tau",
                                             &dynamics.evolving_delta_tau);
                input_nonnegative_double_row(context, "Maximum delta", "delta max",
                                             &dynamics.evolving_delta_max);
                break;
            case SIM_INJECTION_COLLISION_ENERGY:
                input_nonnegative_double_row(context, "Recovered collision energy", "deltaE injected",
                                             &dynamics.injected_energy);
                input_nonnegative_double_row(context, "Background injected energy", "additional energy",
                                             &dynamics.background_injected_energy);
                input_nonnegative_double_row(context, "Recovery time", "taur injection",
                                             &dynamics.injection_recovery_time);
                input_nonnegative_double_row(context, "Recovery exponent", "beta injection",
                                             &dynamics.injection_exponent);
                checkbox_row(context, "Randomize injected energy", "stochastic injection",
                             &dynamics.randomize_injected_energy);
                checkbox_row(context, "Simplified recovery", "piecewise injection",
                             &dynamics.simplified_energy_injection);
                break;
            case SIM_INJECTION_ROUGH_ROTATIONAL:
                input_bounded_double_row(context, "Tangential restitution", "rough beta spin",
                                         &dynamics.rough_tangential_restitution,
                                         -1.0, 1.0);
                input_double_row(context, "Angular kick", "delta omega spin",
                                 &dynamics.rough_angular_kick);
                checkbox_row(context, "Additional tangential delta",
                             "rough translational delta kick",
                             &dynamics.rough_tangential_delta,
                             "Adds the legacy translational tangential kick on top of rough rotational contact dynamics.");
                input_nonnegative_double_row(context, "Tangential delta", "rough delta kick",
                                             &dynamics.delta, nullptr,
                                 dynamics.rough_tangential_delta ? nullptr
                                    : "Enable Additional tangential delta first.");
                break;
            default: break;
        }
        ImGui::EndTable();
        const bool fixed_axis_tangent =
            dynamics.injection == SIM_INJECTION_DELTA_TANGENT ||
            dynamics.injection == SIM_INJECTION_DOUBLE_DELTA_TANGENT ||
            (dynamics.injection == SIM_INJECTION_ROUGH_ROTATIONAL &&
             dynamics.rough_tangential_delta);
        if (fixed_axis_tangent &&
            context.draft.geometry.boundary ==
                SIM_BOUNDARY_CIRCULAR_SPHERICAL &&
            !context.draft.geometry.thermal_wall &&
            context.draft.geometry.wall_tangent_restitution >= 1.0 - 1e-12) {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextColored(ImVec4(0.95f, 0.68f, 0.20f, 1.0f),
                "A specular curved wall cannot absorb the one-signed angular momentum injected by u_z x n. Strong drive can spin up indefinitely; lower wall tangential restitution or use a thermal wall/damping.");
            ImGui::PopTextWrapPos();
        }
    }

    ImGui::SeparatorText("Square-well potential");
    if (begin_properties("dynamics-potential")) {
        const char *physics_toggle_reason =
            context.capabilities.runtime_physics_toggles ? nullptr
            : "Runtime physics selection is unavailable in this build.";
        const bool lees_edwards =
            dynamics.flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR;
        const char *potential_toggle_reason = physics_toggle_reason;
        if (polygon_mode) {
            potential_toggle_reason =
                "Square-well events are not implemented for rotating convex particles.";
        }
        if (potential_toggle_reason == nullptr && lees_edwards && !potential.enabled) {
            potential_toggle_reason =
                "Square-well events cannot be enabled with Lees-Edwards shear.";
        }
        checkbox_row(context, "Enable square well", "potential attraction well",
                     &potential.enabled,
                     "Adds exact well-entry and well-exit events around the hard core.",
                     potential_toggle_reason);
        const char *well_disabled = potential.enabled ? nullptr
            : "Enable the square-well potential first.";
        slider_double_row(context, "Well energy U", "potential depth attraction repulsion",
                         &potential.energy, -10.0, 10.0, nullptr, well_disabled);
        slider_double_row(context, "Well range / contact", "sigma potential radius",
                         &potential.range, 1.0, 5.0,
                         "Multiplicative range relative to contact distance.", well_disabled);
        const bool charged_initializer =
            context.draft.initialization.mode == SIM_INIT_RANDOM_GROWTH ||
            context.draft.initialization.mode == SIM_INIT_MEMBRANE_OSMOSIS;
        const char *charged_reason = !context.capabilities.runtime_charged_potential
            ? "Charge-dependent potential selection is fixed at compile time."
            : (!charged_initializer
                ? "This initializer does not assign particle charges; use random growth or the membrane initializer."
                : nullptr);
        const bool file_species =
            context.draft.initialization.mode == SIM_INIT_DUMP_FILE;
        const bool forced_binary_species =
            context.draft.initialization.mode == SIM_INIT_SQUARE_LATTICE ||
            context.draft.initialization.mode == SIM_INIT_MEMBRANE_OSMOSIS;
        const bool effectively_monodisperse = !file_species &&
            !forced_binary_species &&
            context.draft.initialization.mode != SIM_INIT_HEX_LATTICE &&
            context.draft.initialization.composition ==
                SIM_COMPOSITION_MONODISPERSE;
        const char *same_species_reason =
            !context.capabilities.runtime_same_species_potential
                ? "Same-species potential selection is fixed at compile time."
                : (effectively_monodisperse
                    ? "All particles have the same species, so this rule is indistinguishable from All pairs."
                    : nullptr);
        const EnumOption species[] = {
            {SIM_POTENTIAL_ALL_PAIRS, "All pairs", nullptr},
            {SIM_POTENTIAL_CHARGED, "Charge-dependent", charged_reason},
            {SIM_POTENTIAL_SAME_SPECIES, "Same species only", same_species_reason},
        };
        int species_mode = potential.species_mode;
        const char *species_fixed =
            potential.species_mode == SIM_POTENTIAL_CHARGED &&
            !context.capabilities.runtime_charged_potential
                ? "The compiled charge-dependent interaction rule cannot be changed at runtime."
                : (potential.species_mode == SIM_POTENTIAL_SAME_SPECIES &&
                   !context.capabilities.runtime_same_species_potential
                    ? "The compiled same-species interaction rule cannot be changed at runtime."
                    : well_disabled);
        if (enum_row(context, "Species rule", "charged liquid liquid pairs",
                     &species_mode, species, std::size(species), nullptr,
                     species_fixed)) {
            potential.species_mode =
                static_cast<SimPotentialSpeciesMode>(species_mode);
        }
        if (potential.species_mode == SIM_POTENTIAL_CHARGED) {
            input_bounded_double_row(context, "Positive charge fraction", "charge plus composition",
                                     &potential.positive_charge_fraction, 0.0, 1.0,
                                     nullptr, well_disabled);
            input_bounded_double_row(context, "Neutral charge fraction", "charge neutral composition",
                                     &potential.neutral_charge_fraction, 0.0, 1.0,
                                     nullptr, well_disabled);
        }
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Thermostat and global driving");
    if (begin_properties("dynamics-thermostat")) {
        const char *polygon_translation_only_reason = polygon_mode
            ? "This thermostat changes translation without consistently thermalizing polygon rotation."
            : nullptr;
        const EnumOption thermostats[] = {
            {SIM_THERMOSTAT_NONE, "None", nullptr},
            {SIM_THERMOSTAT_LANGEVIN, "Langevin / OU kicks", nullptr},
            {SIM_THERMOSTAT_VELOCITY_RESCALE, "Global velocity rescale",
                polygon_translation_only_reason},
            {SIM_THERMOSTAT_REMOVE_COM, "Remove center-of-mass drift",
                polygon_translation_only_reason},
            {SIM_THERMOSTAT_MAXWELL_REDRAW, "Maxwell velocity redraw",
                polygon_translation_only_reason},
        };
        int thermostat = driving.thermostat;
        if (enum_row(context, "Thermostat", "noise temperature rescale active OU",
                     &thermostat, thermostats, std::size(thermostats),
                     polygon_mode
                        ? "Analytical Langevin kicks thermalize both polygon translation and rotation."
                        : nullptr)) {
            driving.thermostat = static_cast<SimThermostatMode>(thermostat);
            if (polygon_mode) driving.euler_noise = false;
        }
        const bool thermostat_enabled =
            driving.thermostat != SIM_THERMOSTAT_NONE;
        const char *thermostat_off = thermostat_enabled ? nullptr
            : "Select a thermostat first.";
        const bool thermostat_uses_temperature =
            driving.thermostat == SIM_THERMOSTAT_LANGEVIN ||
            driving.thermostat == SIM_THERMOSTAT_VELOCITY_RESCALE ||
            driving.thermostat == SIM_THERMOSTAT_MAXWELL_REDRAW;
        const bool thermostat_uses_gamma =
            driving.thermostat == SIM_THERMOSTAT_LANGEVIN;
        const bool thermostat_uses_euler =
            driving.thermostat == SIM_THERMOSTAT_LANGEVIN;
        const char *temperature_off =
            thermostat_uses_temperature || context.draft.geometry.thermal_wall
                ? nullptr
                : "This forcing mode does not use a bath temperature; enable a thermal wall or a temperature-based thermostat.";
        const char *gamma_off = thermostat_uses_gamma ? nullptr
            : (dynamics.flight == SIM_FLIGHT_LINEAR_DAMPING
                ? "Adjust the shared gamma with Damping gamma in the Free flight section."
                : "Only the Langevin thermostat uses noise gamma.");
        const char *euler_off = polygon_mode
            ? "Polygon rotation uses the analytical OU update to preserve the exact bath equilibrium."
            : (thermostat_uses_euler ? nullptr
                : "The Euler/OU choice applies only to Langevin thermostat kicks.");
        input_nonnegative_double_row(context, "Thermostat interval", "dt noise kick",
                                     &driving.thermostat_interval, nullptr, thermostat_off);
        input_nonnegative_double_row(context, "Temperature", "T thermostat bath",
                                     &driving.temperature,
                         context.draft.geometry.thermal_wall
                            ? "Shared by the selected thermostat and thermal-wall velocity distribution."
                            : nullptr,
                         temperature_off);
        if (input_nonnegative_double_row(context, "Noise gamma", "gamma thermostat OU",
                                         &driving.gamma,
                                         "The legacy engine shares this value with analytical flight damping.",
                                         gamma_off) &&
            dynamics.flight == SIM_FLIGHT_LINEAR_DAMPING) {
            dynamics.damping_gamma = driving.gamma;
        }
        checkbox_row(context, "Euler noise integrator", "Euler thermostat",
                     &driving.euler_noise,
                     "Use the legacy Euler kick instead of the analytical OU update.",
                     euler_off);
        ImGui::EndTable();
    }

}

void draw_boundaries(Context &context) {
    SimGeometryConfig &geometry = context.draft.geometry;
    const bool is3d = context.capabilities.dimensions == 3;
    const bool polygon_mode = !is3d &&
        context.draft.particle_shape.mode == SIM_PARTICLE_CONVEX_POLYGONS;
    if (!is3d && geometry.boundary == SIM_BOUNDARY_WALL_XYZ) {
        /* XYZ from an imported 3D configuration has an exact 2D analogue. */
        geometry.boundary = SIM_BOUNDARY_WALL_XY;
    }
    if (is3d) {
        geometry.center_membrane = false;
    }
    const bool file_initialization =
        context.draft.initialization.mode == SIM_INIT_DUMP_FILE;
    const bool membrane_initialization =
        context.draft.initialization.mode == SIM_INIT_MEMBRANE_OSMOSIS;
    if (membrane_initialization) {
        /* The initializer places species around the center separator. Keep the
           draft honest instead of offering a checkbox that would be ignored. */
        geometry.center_membrane = true;
    }
    const bool lees_edwards =
        context.draft.dynamics.flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR;
    const bool periodic_initializer_required =
        context.draft.initialization.mode == SIM_INIT_SQUARE_LATTICE ||
        context.draft.initialization.mode == SIM_INIT_HEX_LATTICE;
    const char *planar_wall_reason = lees_edwards
        ? "Lees-Edwards shear requires a fully periodic domain."
        : (periodic_initializer_required
            ? "This lattice initializer places radius-1 particle centers on box edges and therefore requires periodic boundaries."
            : nullptr);
    const char *full_box_reason = !context.capabilities.wall_z
            ? "This 3D build does not provide walls normal to z."
            : planar_wall_reason;
    const char *curved_boundary_reason = nullptr;
    if (!context.capabilities.spherical_confinement) {
        curved_boundary_reason =
            "Circular/spherical confinement is unavailable in this build.";
    } else if (context.draft.initialization.mode != SIM_INIT_RANDOM_GROWTH) {
        curved_boundary_reason =
            "Circular/spherical confinement currently supports only Random placement + growth.";
    } else if (geometry.center_membrane || membrane_initialization) {
        curved_boundary_reason =
            "Circular/spherical confinement cannot be combined with the planar center membrane.";
    } else if (file_initialization) {
        curved_boundary_reason =
            "Circular/spherical confinement cannot load a file because loaded particles are not checked against the curved wall.";
    } else if (lees_edwards) {
        curved_boundary_reason =
            "Lees-Edwards shear requires a fully periodic domain.";
    }
    ImGui::SeparatorText("Domain topology");
    section_note(is3d
        ? "Individual x/y walls remain periodic in z; the full wall box and spherical confinement constrain all three axes."
        : "A wall replaces periodic wrapping only along its normal axis.");

    if (begin_properties("boundary-topology")) {
        int boundary = geometry.boundary;
        bool boundary_changed = false;
        if (is3d) {
            const EnumOption boundaries_3d[] = {
                {SIM_BOUNDARY_PERIODIC, "Periodic", nullptr},
                {SIM_BOUNDARY_WALL_X, "Walls normal to x", planar_wall_reason},
                {SIM_BOUNDARY_WALL_Y, "Walls normal to y", planar_wall_reason},
                {SIM_BOUNDARY_WALL_XY, "x/y wall box", planar_wall_reason},
                {SIM_BOUNDARY_WALL_XYZ, "Full x/y/z wall box", full_box_reason},
                {SIM_BOUNDARY_CIRCULAR_SPHERICAL,
                    "Spherical confinement", curved_boundary_reason},
            };
            boundary_changed = enum_row(
                context, "Boundary", "periodic walls box sphere confinement",
                &boundary, boundaries_3d, std::size(boundaries_3d),
                "Boundary choices are exclusive; contradictory wall toggles cannot be constructed.");
        } else {
            const EnumOption boundaries_2d[] = {
                {SIM_BOUNDARY_PERIODIC, "Periodic", nullptr},
                {SIM_BOUNDARY_WALL_X, "Walls normal to x", planar_wall_reason},
                {SIM_BOUNDARY_WALL_Y, "Walls normal to y", planar_wall_reason},
                {SIM_BOUNDARY_WALL_XY, "Rectangular wall box", planar_wall_reason},
                {SIM_BOUNDARY_CIRCULAR_SPHERICAL,
                    "Circular confinement", curved_boundary_reason},
            };
            boundary_changed = enum_row(
                context, "Boundary", "periodic walls box circle confinement",
                &boundary, boundaries_2d, std::size(boundaries_2d),
                "Boundary choices are exclusive; contradictory wall toggles cannot be constructed.");
        }
        if (boundary_changed) {
            geometry.boundary = static_cast<SimBoundaryMode>(boundary);
            if (geometry.boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL) {
                /* The engine derives an isotropic diameter for curved
                   confinement; canonicalize hidden box fields immediately. */
                context.draft.initialization.aspect_ratio = 1.0;
                context.draft.initialization.box_length_mode =
                    SIM_BOX_FROM_PACKING_FRACTION;
            }
        }
        const char *membrane_reason = membrane_initialization
            ? "The membrane/osmosis initializer requires this center membrane."
            : (geometry.boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL
               && !geometry.center_membrane
                ? "The planar center membrane cannot be combined with circular confinement."
                : (periodic_initializer_required && !geometry.center_membrane
                    ? "This initializer requires a fully periodic domain."
                : (lees_edwards && !geometry.center_membrane
                    ? "Lees-Edwards shear requires a fully periodic domain."
                    : nullptr)));
        if (!is3d) {
            checkbox_row(context, "Center membrane",
                         "mid wall osmosis separator",
                         &geometry.center_membrane,
                         "Adds the species-selective wall at the box center.",
                         membrane_reason);
        }
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Wall collision law");
    if (begin_properties("boundary-law")) {
        const bool has_wall = geometry.boundary != SIM_BOUNDARY_PERIODIC ||
                              geometry.center_membrane ||
                              membrane_initialization;
        const char *no_wall = has_wall ? nullptr
            : "Select a physical wall boundary first.";
        const char *normal_restitution_reason = geometry.thermal_wall
            ? "A thermal wall redraws the outgoing normal velocity, replacing normal restitution."
            : no_wall;
        slider_double_row(context, "Normal restitution", "wall resW bounce",
                         &geometry.wall_normal_restitution, 0.0, 1.0,
                         "Restitution of the velocity component normal to the wall.",
                         normal_restitution_reason);
        const bool has_x_or_curved_wall =
            geometry.boundary == SIM_BOUNDARY_WALL_X ||
            geometry.boundary == SIM_BOUNDARY_WALL_XY ||
            geometry.boundary == SIM_BOUNDARY_WALL_XYZ ||
            geometry.boundary == SIM_BOUNDARY_CIRCULAR_SPHERICAL ||
            geometry.center_membrane || membrane_initialization;
        const bool tangent_restitution_active = !geometry.thermal_wall &&
            (has_x_or_curved_wall ||
             ((geometry.boundary == SIM_BOUNDARY_WALL_Y ||
               geometry.boundary == SIM_BOUNDARY_WALL_XY ||
               geometry.boundary == SIM_BOUNDARY_WALL_XYZ) &&
              (is3d || !geometry.shear_wall)));
        const char *tangent_help = geometry.shear_wall && is3d
            ? "Moving y walls prescribe outgoing vx; this multiplier still applies to the other tangential component vz."
            : (geometry.shear_wall && has_x_or_curved_wall
                ? "Moving y walls prescribe outgoing vx there; this multiplier still applies at x or curved walls."
                : "Multiplier applied to tangential wall velocity components.");
        const char *tangent_restitution_reason = geometry.thermal_wall
            ? "A thermal wall redraws tangential velocity components instead of applying tangential restitution."
            : (!tangent_restitution_active
                ? (no_wall != nullptr ? no_wall
                    : "Moving y walls prescribe the only tangential velocity component in this 2D geometry.")
                : nullptr);
        slider_double_row(context, "Tangential restitution", "wall resT slip",
                         &geometry.wall_tangent_restitution, -1.0, 1.0,
                         tangent_help, tangent_restitution_reason);
        checkbox_row(context, "Thermal wall", "Maxwell diffuse temperature wall",
                     &geometry.thermal_wall,
                     "Draws outgoing velocities from the wall-temperature distribution.",
                     polygon_mode
                        ? "Thermal-wall redraws do not yet sample polygon spin consistently."
                        : (no_wall != nullptr && !geometry.thermal_wall ? no_wall : nullptr));
        if (geometry.thermal_wall) {
            input_nonnegative_double_row(
                context, "Thermal-wall temperature",
                "Maxwell diffuse wall bath T",
                &context.draft.driving.temperature,
                "The legacy engine shares this temperature with thermostat baths.",
                no_wall);
        }

        const bool y_wall = geometry.boundary == SIM_BOUNDARY_WALL_Y ||
                            geometry.boundary == SIM_BOUNDARY_WALL_XY ||
                            geometry.boundary == SIM_BOUNDARY_WALL_XYZ;
        const char *shear_reason = polygon_mode
            ? "Moving shear-wall contact is not yet implemented for rotating convex particles."
            : (!y_wall && !geometry.shear_wall
                ? "Shear walls require walls normal to y." : nullptr);
        checkbox_row(context, "Moving shear walls", "Couette shear wall",
                     &geometry.shear_wall,
                     is3d
                        ? "The y-wall speed prescribes vx; vz remains a separate tangential component."
                        : nullptr,
                     shear_reason);
        input_double_row(context, "Shear wall speed", "Couette wall velocity",
                         &geometry.shear_wall_speed, nullptr,
                         geometry.shear_wall && y_wall ? nullptr
                            : (geometry.shear_wall
                                ? "Moving shear walls require walls normal to y."
                                : "Enable moving shear walls first."));
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Compatibility");
    ImGui::TextWrapped(
        "The validator also checks the selected trajectory law. Lees-Edwards "
        "shear cannot be composed with physical walls, square-well/charged "
        "events, constant field, or damping in the current analytical engine.");
}

void rgba_to_float(uint32_t rgba, float color[4]) {
    color[0] = static_cast<float>((rgba >> 24) & 0xffu)/255.0f;
    color[1] = static_cast<float>((rgba >> 16) & 0xffu)/255.0f;
    color[2] = static_cast<float>((rgba >> 8) & 0xffu)/255.0f;
    color[3] = static_cast<float>(rgba & 0xffu)/255.0f;
}

uint32_t float_to_rgba(const float color[4]) {
    const auto byte = [](float value) -> uint32_t {
        return static_cast<uint32_t>(std::lround(
            std::clamp(value, 0.0f, 1.0f)*255.0f));
    };
    return (byte(color[0]) << 24) | (byte(color[1]) << 16) |
           (byte(color[2]) << 8) | byte(color[3]);
}

bool color_row(Context &context, const char *label, const char *keywords,
               uint32_t *rgba, const char *help = nullptr) {
    if (!edmd_gui::matches_filter(context, label, keywords)) return false;
    property_label(label, help);
    float color[4];
    rgba_to_float(*rgba, color);
    ImGui::PushID(label);
    const bool changed = ImGui::ColorEdit4("##value", color,
        ImGuiColorEditFlags_AlphaBar | ImGuiColorEditFlags_AlphaPreviewHalf);
    ImGui::PopID();
    if (changed) *rgba = float_to_rgba(color);
    return changed;
}

void draw_appearance(Context &context, GuiLiveState &live) {
    const bool is3d = context.capabilities.dimensions == 3;
    if (is3d && (live.color_mode == SIM_COLOR_BY_HEX_ORDER ||
                 live.color_mode == SIM_COLOR_BY_SQUARE_ORDER)) {
        live.color_mode = SIM_COLOR_BY_TYPE;
    }
    ImGui::SeparatorText("Scene");
    section_note("Appearance changes are live. They do not rebuild particles, cells, or events.");
    if (begin_properties("appearance-scene")) {
        checkbox_row(context, "Show simulation scene", "viewport visibility",
                     &live.show_scene);
        checkbox_row(context, "Show boundary", "box walls outline",
                     &live.show_boundary);
        if (is3d) {
            checkbox_row(context, "Show axes", "3D axis orientation",
                         &live.show_axes);
        }
        checkbox_row(context, "Selection outline", "selected particle highlight",
                     &live.show_selection_outline);
        checkbox_row(context, "Velocity arrows", "vector arrows velocity",
                     &live.velocity_arrows);
        checkbox_row(context, "Show render FPS", "frames performance",
                     &live.show_fps);
        checkbox_row(context, "Show help", "instructions shortcuts",
                     &live.show_help);
        if (slider_int_row(context, "Target FPS", "frame rate vsync render",
                          &live.target_fps, 1, 1000,
                          "The host applies this limit without changing simulated time. Range: 1-1000.",
                          nullptr, 1, false)) {
            live.target_fps = std::clamp(live.target_fps, 1, 1000);
        }
        if (slider_double_row(context, "Scene supersampling", "render target quality retina",
                             &live.scene_supersampling, 0.5, 3.0,
                             "Multiplier on the native-DPI render target. 1 is one scene pixel per physical display pixel; range: 0.5-3.",
                             nullptr, "%.6g", false, false)) {
            clamp_live_double(live.scene_supersampling, 0.5, 3.0, 1.0);
        }
        if (is3d) {
            const EnumOption lighting_modes_3d[] = {
                {SIM_LIGHTING_SIMPLE, "Simple", nullptr},
                {SIM_LIGHTING_STUDIO, "Studio", nullptr},
                {SIM_LIGHTING_DEPTH_ENHANCED, "Depth enhanced", nullptr},
            };
            int lighting = live.lighting_mode;
            if (enum_row(
                    context, "Lighting", "3D light shading depth rim specular",
                    &lighting, lighting_modes_3d,
                    std::size(lighting_modes_3d),
                    "Simple is the least expensive Lambert shading. Studio adds hemisphere, key/fill, and soft specular light. Depth enhanced also adds subtle silhouette, rim, and distance cues; it is a single-pass effect, not ambient occlusion or ray tracing.")) {
                live.lighting_mode = static_cast<SimLightingMode>(lighting);
            }

            if (slider_double_row(
                    context, "Camera field of view", "3D FOV perspective",
                    &live.camera_fov_degrees, 10.0, 120.0,
                    "Vertical field of view in degrees. Range: 10-120.",
                    nullptr, "%.6g", false, false)) {
                clamp_live_double(live.camera_fov_degrees,
                                  10.0, 120.0, 45.0);
            }

            const EnumOption sphere_quality_3d[] = {
                {SIM_SPHERE_QUALITY_AUTO, "Automatic", nullptr},
                {SIM_SPHERE_QUALITY_LOW, "Low", nullptr},
                {SIM_SPHERE_QUALITY_MEDIUM, "Medium", nullptr},
                {SIM_SPHERE_QUALITY_HIGH, "High", nullptr},
                {SIM_SPHERE_QUALITY_ULTRA, "Ultra", nullptr},
            };
            int quality = live.sphere_quality;
            if (enum_row(
                    context, "Sphere tessellation",
                    "sphere quality rings slices", &quality,
                    sphere_quality_3d, std::size(sphere_quality_3d),
                    "Automatic quality reduces geometry as particle count grows.")) {
                live.sphere_quality = static_cast<SimSphereQuality>(quality);
            }
        }
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Particle coloring");
    if (begin_properties("appearance-color")) {
        int color_mode = live.color_mode;
        bool color_mode_changed = false;
        if (is3d) {
            const EnumOption color_modes_3d[] = {
                {SIM_COLOR_BY_TYPE, "Particle type", nullptr},
                {SIM_COLOR_BY_COLLISIONS, "Collision count", nullptr},
                {SIM_COLOR_BY_SPEED, "Speed", nullptr},
                {SIM_COLOR_BY_RADIUS, "Radius", nullptr},
                {SIM_COLOR_BY_CHARGE, "Charge", nullptr},
            };
            color_mode_changed = enum_row(
                context, "Color mode", "particle color speed charge radius",
                &color_mode, color_modes_3d, std::size(color_modes_3d));
        } else {
            const EnumOption color_modes_2d[] = {
                {SIM_COLOR_BY_TYPE, "Particle type", nullptr},
                {SIM_COLOR_BY_COLLISIONS, "Collision count", nullptr},
                {SIM_COLOR_BY_SPEED, "Speed", nullptr},
                {SIM_COLOR_BY_HEX_ORDER, "Hexagonal order", nullptr},
                {SIM_COLOR_BY_SQUARE_ORDER, "Square order", nullptr},
                {SIM_COLOR_BY_RADIUS, "Radius", nullptr},
                {SIM_COLOR_BY_CHARGE, "Charge", nullptr},
            };
            color_mode_changed = enum_row(
                context, "Color mode",
                "particle color BOOP speed charge radius", &color_mode,
                color_modes_2d, std::size(color_modes_2d));
        }
        if (color_mode_changed) {
            live.color_mode = static_cast<SimColorMode>(color_mode);
        }
        const EnumOption palettes[] = {
            {SIM_PALETTE_PLASMA, "Plasma", nullptr},
            {SIM_PALETTE_VIRIDIS, "Viridis", nullptr},
            {SIM_PALETTE_COPPER, "Copper", nullptr},
        };
        int palette = live.palette;
        if (enum_row(context, "Continuous palette", "colormap plasma viridis copper",
                     &palette, palettes, std::size(palettes),
                     "Used by collision, speed, radius, charge, and order metrics.")) {
            live.palette = static_cast<SimColorPalette>(palette);
        }
        if (input_double_row(context, "Color refresh interval", "color cache rate",
                             &live.color_refresh_seconds,
                             "Wall-clock seconds between live color-cache updates. Range: 0.001-3600.",
                             nullptr, "%.10g", false, 0.001, 3600.0)) {
            clamp_live_double(live.color_refresh_seconds, 0.001, 3600.0,
                              1.0/30.0);
        }
        color_row(context, "Type 0 color", "small particle rgba",
                  &live.type0_rgba);
        color_row(context, "Type 1 color", "large particle rgba",
                  &live.type1_rgba);
        color_row(context, "Background color", "scene clear rgba",
                  &live.background_rgba);
        color_row(context, "Boundary color", "wall box sphere rgba",
                  &live.boundary_rgba);
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Live structure factor");
    if (begin_properties("appearance-structure")) {
        checkbox_row(context, "Structure-factor window", "graphical S q live",
                     &live.graphical_structure_factor,
                     is3d ? "Displays the current XY projection in 3D." : nullptr);
        if (slider_int_row(context, "Particle sample cap", "structure factor sample performance",
                          &live.structure_factor_sample_cap, 1, 8192,
                          "At most 8192 particles are sampled.",
                          live.graphical_structure_factor ? nullptr
                              : "Enable the structure-factor window first.",
                          1, false)) {
            live.structure_factor_sample_cap =
                std::clamp(live.structure_factor_sample_cap, 1, 8192);
        }
        ImGui::EndTable();
    }

    /* Also sanitize values loaded from older configuration schemas before the
       renderer consumes these live fields. */
    live.target_fps = std::clamp(live.target_fps, 1, 1000);
    if (live.lighting_mode < SIM_LIGHTING_SIMPLE ||
        live.lighting_mode > SIM_LIGHTING_DEPTH_ENHANCED) {
        live.lighting_mode = SIM_LIGHTING_STUDIO;
    }
    clamp_live_double(live.scene_supersampling, 0.5, 3.0, 1.0);
    clamp_live_double(live.camera_fov_degrees, 10.0, 120.0, 45.0);
    live.structure_factor_sample_cap =
        std::clamp(live.structure_factor_sample_cap, 1, 8192);
    clamp_live_double(live.structure_factor_refresh_seconds,
                      0.001, 3600.0, 0.5);
    clamp_live_double(live.color_refresh_seconds,
                      0.001, 3600.0, 1.0/30.0);
}

void draw_capability_row(const char *name, bool available,
                         const char *unavailable_text = "not compiled") {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::TextUnformatted(name);
    ImGui::TableSetColumnIndex(1);
    if (available) {
        ImGui::TextColored(ImVec4(0.40f, 0.80f, 0.45f, 1), "available");
    } else {
        ImGui::TextDisabled("%s", unavailable_text);
    }
}

void draw_advanced(Context &context) {
    SimIntruderConfig &intruder = context.draft.intruder;
    const bool is3d = context.capabilities.dimensions == 3;
    const bool polygon_mode = !is3d &&
        context.draft.particle_shape.mode == SIM_PARTICLE_CONVEX_POLYGONS;
    if (is3d && !context.capabilities.intruder_3d) {
        intruder.enabled = false;
    }

    if (collapsing_header(context, "Build capabilities", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::TextDisabled("Configuration schema %d; %dD graphics build",
                            context.draft.schema_version,
                            context.capabilities.dimensions);
        if (ImGui::BeginTable("advanced-capabilities", 2,
                              ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_BordersInnerV |
                              ImGuiTableFlags_SizingStretchProp)) {
            draw_capability_row("Runtime physics selection",
                                context.capabilities.runtime_physics_toggles);
            draw_capability_row("Runtime basic collision injection",
                                context.capabilities.runtime_basic_injection_toggles);
            draw_capability_row("Runtime advanced collision injection",
                                context.capabilities.runtime_advanced_injection_toggles);
            draw_capability_row("Runtime exponential restitution",
                                context.capabilities.runtime_exponential_restitution);
            draw_capability_row("Runtime charged potential",
                                context.capabilities.runtime_charged_potential);
            draw_capability_row("Runtime same-species potential",
                                context.capabilities.runtime_same_species_potential);
            draw_capability_row("Runtime saved-output lifecycle",
                                context.capabilities.runtime_output_lifecycle,
                                "unavailable");
            draw_capability_row("Runtime output-column selection",
                                context.capabilities.runtime_output_column_toggles,
                                "compile-time fixed");
            draw_capability_row("Runtime intruder parameters",
                                context.capabilities.runtime_intruder_parameters);
            draw_capability_row("Full stress transport", context.capabilities.stress);
            draw_capability_row("Heat-current transport", context.capabilities.heat);
            draw_capability_row("Rotational rough particles",
                                context.capabilities.rough_particles);
            draw_capability_row("Intruder", context.capabilities.intruder);
            if (!is3d) {
                draw_capability_row("Rigid convex particles",
                                    context.capabilities.convex_polygon_particles);
            }
            draw_capability_row("Radius-sum / additive contacts",
                                context.capabilities.additive_contacts);
            draw_capability_row("Circular / spherical confinement",
                                context.capabilities.spherical_confinement);
            draw_capability_row("3D tangential delta",
                                !is3d || context.capabilities.tangent_delta_3d);
            draw_capability_row("3D z-normal walls",
                                !is3d || context.capabilities.wall_z,
                                "full x/y/z box unavailable");
            draw_capability_row("3D Lees-Edwards shear",
                                !is3d || context.capabilities.shear_3d);
            draw_capability_row("Full 3D thermostats",
                                !is3d || context.capabilities.full_3d_thermostats,
                                "some exposed modes remain planar");
            draw_capability_row("3D HCP lattice initializer",
                                !is3d || context.capabilities.lattice_initializers_3d,
                                "complete z/vz state unavailable");
            draw_capability_row("3D intruder collisions",
                                !is3d || context.capabilities.intruder_3d);
            draw_capability_row(
                "Saved local pressure field",
                (!is3d || context.capabilities.pressure_field_3d) &&
                    !polygon_mode,
                polygon_mode
                    ? "disk collision impulses only"
                    : "3D field unavailable");
            draw_capability_row("XY-projected 3D observables",
                                !is3d || context.capabilities.projected_xy_observables);
            draw_capability_row("Queue invariant checks",
                                context.capabilities.queue_verification);
            draw_capability_row("Overlap invariant checks",
                                context.capabilities.overlap_verification);
            draw_capability_row("Collision-accounting checks",
                                context.capabilities.collision_accounting_verification);
            ImGui::EndTable();
        }
    }

    if ((!is3d || context.capabilities.intruder_3d) &&
        collapsing_header(context, "Intruder")) {
        const char *intruder_build_reason = !context.capabilities.intruder
            ? "Requires a build with INTRUDER enabled."
            : ((context.capabilities.stress || context.capabilities.heat)
                    ? "Intruder collisions do not yet maintain complete STRESS/HEAT transport bookkeeping."
                    : nullptr);
        const char *intruder_enable_reason = polygon_mode
            ? "The intruder and rigid convex-particle collision models cannot be active together."
            : (intruder_build_reason != nullptr
            ? intruder_build_reason
            : "Intruder presence is fixed by the INTRUDER build option; live reload cannot add or remove it.");
        if (begin_properties("advanced-intruder")) {
            checkbox_row(context, "Enable intruder", "polygon obstacle body",
                         &intruder.enabled, nullptr, intruder_enable_reason);
            const char *off = intruder_build_reason != nullptr
                ? intruder_build_reason
                : (!context.capabilities.runtime_intruder_parameters
                    ? "Intruder parameters are fixed at compile time in this build."
                    : (!intruder.enabled
                        ? "The intruder is not enabled in this build."
                        : nullptr));
            const EnumOption shapes[] = {
                {0, "Square", nullptr},
                {1, "Triangle", nullptr},
                {2, "Wheel", nullptr},
                {3, "Chiral wheel", nullptr},
                {4, "Generalized wheel", nullptr},
                {5, "Complex wheel", nullptr},
                {6, "Cross", nullptr},
                {7, "Regular polygon", nullptr},
            };
            enum_row(context, "Shape", "intruder polygon wheel cross triangle",
                     &intruder.shape, shapes, std::size(shapes), nullptr, off);
            switch (intruder.shape) {
                case 0:
                    input_nonnegative_double_row(
                        context, "Side length", "intruder square size",
                        &intruder.shape_parameters[0], nullptr, off);
                    break;
                case 1:
                    input_nonnegative_double_row(
                        context, "Base length", "intruder triangle size",
                        &intruder.shape_parameters[0], nullptr, off);
                    input_nonnegative_double_row(
                        context, "Height", "intruder triangle size",
                        &intruder.shape_parameters[1], nullptr, off);
                    break;
                case 2:
                case 3:
                case 4:
                    input_nonnegative_double_row(
                        context, "Inner radius", "intruder wheel valley radius",
                        &intruder.shape_parameters[0], nullptr, off);
                    input_nonnegative_double_row(
                        context, "Outer radius", "intruder wheel tip radius",
                        &intruder.shape_parameters[1], nullptr, off);
                    input_bounded_double_row(
                        context, "Spike count", "intruder wheel vertices integer",
                        &intruder.shape_parameters[2], 2.0, 32.0,
                        "Integer from 2 to 32.", off, "%.0f");
                    if (intruder.shape == 4) {
                        input_double_row(context, "Chirality", "intruder generalized wheel shift",
                                         &intruder.shape_parameters[3],
                                         "Angular offset from the achiral halfway position.", off);
                    }
                    break;
                case 5:
                    input_nonnegative_double_row(
                        context, "Base radius", "intruder complex wheel size",
                        &intruder.shape_parameters[0], nullptr, off);
                    input_double_row(context, "Spike amplitude", "intruder complex wheel ripple",
                                     &intruder.shape_parameters[1], nullptr, off);
                    input_double_row(context, "Triangular amplitude", "intruder complex wheel anisotropy",
                                     &intruder.shape_parameters[2], nullptr, off);
                    input_double_row(context, "Phase", "intruder complex wheel angle",
                                     &intruder.shape_parameters[3], nullptr, off);
                    break;
                case 6:
                    input_nonnegative_double_row(
                        context, "Arm extension", "intruder cross geometry",
                        &intruder.shape_parameters[0], nullptr, off);
                    input_nonnegative_double_row(
                        context, "Arm width", "intruder cross geometry",
                        &intruder.shape_parameters[1], nullptr, off);
                    break;
                case 7:
                    input_nonnegative_double_row(
                        context, "Radius", "intruder regular polygon size",
                        &intruder.shape_parameters[0], nullptr, off);
                    input_bounded_double_row(
                        context, "Side count", "intruder polygon vertices integer",
                        &intruder.shape_parameters[1], 3.0, 64.0,
                        "Integer from 3 to 64.", off, "%.0f");
                    break;
                default:
                    section_note("Select a recognized shape to edit its geometry.");
                    break;
            }
            input_nonnegative_double_row(context, "Mass", "intruder M",
                                         &intruder.mass, nullptr, off);
            input_double_row(context, "Moment of inertia", "intruder Im rotation",
                             &intruder.moment_of_inertia, nullptr, off);
            input_double_row(context, "Initial x", "intruder position",
                             &intruder.initial_x, nullptr,
                             "The legacy intruder initializer places its center from the box geometry; this schema field is not wired.");
            input_double_row(context, "Initial y", "intruder position",
                             &intruder.initial_y, nullptr,
                             "The legacy intruder initializer places its center from the box geometry; this schema field is not wired.");
            input_double_row(context, "Initial vx", "intruder velocity",
                             &intruder.initial_vx, nullptr, off);
            input_double_row(context, "Initial vy", "intruder velocity",
                             &intruder.initial_vy, nullptr, off);
            input_double_row(context, "Initial angle", "intruder orientation",
                             &intruder.initial_angle, nullptr, off);
            input_double_row(context, "Initial angular velocity", "intruder omega",
                             &intruder.initial_omega, nullptr, off);
            slider_double_row(context, "Collision restitution", "intruder res",
                              &intruder.restitution, 0.0, 1.0, nullptr, off);
            input_nonnegative_double_row(
                context, "Collision delta", "intruder kick",
                &intruder.delta, nullptr, off);
            input_nonnegative_double_row(
                context, "Bath particle radius", "intruder bath size",
                &intruder.bath_particle_radius, nullptr, off);
            checkbox_row(context, "Rethermalize after collision", "intruder thermostat",
                         &intruder.rethermalize_after_collision, nullptr,
                         "Intruder-collision rethermalization is a compile-time constant in the current engine.");
            input_nonnegative_double_row(
                context, "Display particle radius", "intruder dump outline",
                &intruder.display_particle_radius, nullptr, off);
            checkbox_row(context, "Kinetic-model measurement", "intruder KM statistics",
                         &intruder.kinetic_model_measurement, nullptr, off);
            input_nonnegative_double_row(
                context, "Kinetic-model burn time", "intruder KM burn",
                &intruder.kinetic_model_burn_time,
                "Measurements accumulated before this simulation time are discarded.",
                off);
            ImGui::EndTable();
        }
    }

    if (collapsing_header(context, "Configuration diagnostics")) {
        edmd_gui::draw_issue_summary(context, true);
    }
}

} // namespace

namespace edmd_gui {

bool matches_filter(const Context &context, const char *label,
                    const char *keywords) {
    if (context.search[0] == '\0') return true;
    std::string haystack = label ? label : "";
    if (keywords != nullptr) {
        haystack.push_back(' ');
        haystack += keywords;
    }

    const char *cursor = context.search;
    while (*cursor != '\0') {
        while (*cursor != '\0' && std::isspace(static_cast<unsigned char>(*cursor)))
            ++cursor;
        std::string token;
        while (*cursor != '\0' && !std::isspace(static_cast<unsigned char>(*cursor))) {
            token.push_back(static_cast<char>(
                std::tolower(static_cast<unsigned char>(*cursor))));
            ++cursor;
        }
        if (!token.empty() && !contains_case_insensitive(haystack.c_str(), token))
            return false;
    }
    return true;
}

void help_marker(const char *text) {
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered()) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize()*32.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

void disabled_reason(const char *reason) {
    if (reason == nullptr) return;
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize()*32.0f);
        ImGui::TextUnformatted(reason);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

void draw_configuration(Context &context, GuiLiveState &live) {
    ImGui::Begin("Configuration");
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputTextWithHint("##configuration-search",
                             "Search every setting (e.g. thermostat, cluster, wall)",
                             context.search, sizeof(context.search));

    draw_issue_summary(context, false);
    draw_first_error(context);
    ImGui::Separator();
    ImGui::TextDisabled(
        "Changes reload automatically. Numeric sliders update the simulation "
        "while held; exact entry reloads after editing finishes.");

    if (context.search[0] != '\0') {
        ImGui::SeparatorText("Filtered settings");
        if (collapsing_header(context, "Setup", ImGuiTreeNodeFlags_DefaultOpen))
            draw_setup(context);
        if (collapsing_header(context, "Dynamics", ImGuiTreeNodeFlags_DefaultOpen))
            draw_dynamics(context);
        if (collapsing_header(context, "Boundaries", ImGuiTreeNodeFlags_DefaultOpen))
            draw_boundaries(context);
        if (collapsing_header(context, "Appearance", ImGuiTreeNodeFlags_DefaultOpen))
            draw_appearance(context, live);
        if (collapsing_header(context, "Advanced", ImGuiTreeNodeFlags_DefaultOpen))
            draw_advanced(context);
    } else if (ImGui::BeginTabBar("configuration-tabs",
                                  ImGuiTabBarFlags_Reorderable |
                                  ImGuiTabBarFlags_FittingPolicyScroll)) {
        if (ImGui::BeginTabItem("Setup")) {
            draw_setup(context);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Dynamics")) {
            draw_dynamics(context);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Boundaries")) {
            draw_boundaries(context);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Appearance")) {
            draw_appearance(context, live);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Advanced")) {
            draw_advanced(context);
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }

    refresh_validation(context);
    ImGui::End();
}

} // namespace edmd_gui
