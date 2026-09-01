#include "gui_imgui_internal.hpp"

#include "imgui.h"
#include "live_analysis.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>

namespace {

size_t series_index(const LiveAnalysisSeries &series, size_t ordinal) {
    return (series.first + ordinal) % LIVE_ANALYSIS_HISTORY_CAPACITY;
}

void draw_empty_plot_message(const ImVec2 &minimum, const ImVec2 &maximum,
                             const char *message) {
    ImDrawList *draw = ImGui::GetWindowDrawList();
    const ImVec2 text = ImGui::CalcTextSize(message);
    draw->AddText(ImVec2((minimum.x + maximum.x - text.x)*0.5f,
                         (minimum.y + maximum.y - text.y)*0.5f),
                  IM_COL32(145, 153, 166, 255), message);
}

bool series_sample_is_plottable(const LiveAnalysisSeries &series,
                                size_t ordinal, bool log_log) {
    const size_t index = series_index(series, ordinal);
    const double sample_time = series.simulation_time[index];
    const double sample_value = series.value[index];
    return std::isfinite(sample_time) && std::isfinite(sample_value) &&
        (!log_log || (sample_time > 0.0 && sample_value > 0.0));
}

double plot_coordinate(double value, bool logarithmic) {
    return logarithmic ? std::log10(value) : value;
}

void draw_series_plot(const char *id, const char *label,
                      const LiveAnalysisSeries &series, ImU32 color,
                      bool *log_log, bool include_zero = false) {
    ImGui::PushID(id);
    ImGui::Text("%s", label);
    if (series.count > 0) {
        const size_t latest_index = series_index(series, series.count - 1);
        ImGui::SameLine();
        ImGui::TextDisabled("latest %.7g", series.value[latest_index]);
    }
    ImGui::SameLine();
    ImGui::Checkbox("log-log", log_log);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Use logarithmic simulation-time and value axes. "
                          "Zero, negative, NaN, and infinite samples are skipped.");
    const bool use_log_log = *log_log;

    const float width = std::max(80.0f, ImGui::GetContentRegionAvail().x);
    const ImVec2 size(width, 112.0f);
    const ImVec2 minimum = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("plot", size);
    const ImVec2 maximum(minimum.x + size.x, minimum.y + size.y);
    ImDrawList *draw = ImGui::GetWindowDrawList();
    draw->AddRectFilled(minimum, maximum, IM_COL32(18, 21, 27, 255), 4.0f);
    draw->AddRect(minimum, maximum, IM_COL32(68, 75, 88, 255), 4.0f);

    size_t valid_count = 0;
    size_t first_valid = 0;
    size_t last_valid = 0;
    for (size_t ordinal = 0; ordinal < series.count; ++ordinal) {
        if (!series_sample_is_plottable(series, ordinal, use_log_log))
            continue;
        if (valid_count == 0) first_valid = ordinal;
        last_valid = ordinal;
        ++valid_count;
    }
    if (valid_count == 0) {
        draw_empty_plot_message(minimum, maximum,
            series.count == 0 ? "Waiting for simulation progress"
            : (use_log_log ? "Log-log view needs positive finite samples"
                           : "No finite samples are available"));
        ImGui::PopID();
        return;
    }

    const size_t first_index = series_index(series, first_valid);
    const size_t last_index = series_index(series, last_valid);
    const double raw_time_min = series.simulation_time[first_index];
    const double raw_time_max = series.simulation_time[last_index];
    double time_min = plot_coordinate(raw_time_min, use_log_log);
    double time_max = plot_coordinate(raw_time_max, use_log_log);
    double value_min = plot_coordinate(series.value[first_index], use_log_log);
    double value_max = value_min;
    for (size_t ordinal = first_valid + 1; ordinal < series.count; ++ordinal) {
        if (!series_sample_is_plottable(series, ordinal, use_log_log))
            continue;
        const double value = plot_coordinate(
            series.value[series_index(series, ordinal)], use_log_log);
        value_min = std::min(value_min, value);
        value_max = std::max(value_max, value);
    }
    if (include_zero && !use_log_log) value_min = std::min(0.0, value_min);
    if (!(time_max > time_min)) time_max = time_min + 1.0;
    if (!(value_max > value_min)) {
        const double padding = std::max(1e-9, std::fabs(value_min)*0.08);
        value_min -= padding;
        value_max += padding;
    } else {
        const double padding = (value_max - value_min)*0.08;
        value_min -= padding;
        value_max += padding;
    }

    const ImVec2 plot_min(minimum.x + 8.0f, minimum.y + 8.0f);
    const ImVec2 plot_max(maximum.x - 8.0f, maximum.y - 18.0f);
    for (int grid = 1; grid < 4; ++grid) {
        const float fraction = grid/4.0f;
        const float y = plot_min.y + fraction*(plot_max.y - plot_min.y);
        draw->AddLine(ImVec2(plot_min.x, y), ImVec2(plot_max.x, y),
                      IM_COL32(50, 55, 65, 150));
    }

    auto point_for = [&](size_t ordinal) {
        const size_t index = series_index(series, ordinal);
        const double sample_time = plot_coordinate(
            series.simulation_time[index], use_log_log);
        const double sample_value = plot_coordinate(
            series.value[index], use_log_log);
        const float x = plot_min.x + static_cast<float>(
            (sample_time - time_min)/(time_max - time_min))*
            (plot_max.x - plot_min.x);
        const float y = plot_max.y - static_cast<float>(
            (sample_value - value_min)/(value_max - value_min))*
            (plot_max.y - plot_min.y);
        return ImVec2(x, y);
    };

    draw->PushClipRect(plot_min, plot_max, true);
    if (valid_count == 1) {
        draw->AddCircleFilled(point_for(first_valid), 2.5f, color);
    } else {
        ImVec2 previous{};
        bool previous_valid = false;
        for (size_t ordinal = 0; ordinal < series.count; ++ordinal) {
            if (!series_sample_is_plottable(series, ordinal, use_log_log)) {
                /* Do not bridge across a zero/negative gap in log space. */
                previous_valid = false;
                continue;
            }
            const ImVec2 current = point_for(ordinal);
            if (previous_valid)
                draw->AddLine(previous, current, color, 1.8f);
            previous = current;
            previous_valid = true;
        }
    }
    draw->PopClipRect();

    char time_range[96];
    if (use_log_log) {
        std::snprintf(time_range, sizeof(time_range),
                      "log-log  t %.4g - %.4g   %zu positive%s",
                      raw_time_min, raw_time_max, valid_count,
                      valid_count == series.count ? "" : " (others skipped)");
    } else {
        std::snprintf(time_range, sizeof(time_range), "t %.4g - %.4g",
                      raw_time_min, raw_time_max);
    }
    draw->AddText(ImVec2(plot_min.x, maximum.y - 15.0f),
                  IM_COL32(135, 143, 155, 255), time_range);

    if (ImGui::IsItemHovered()) {
        const float fraction = std::clamp(
            (ImGui::GetIO().MousePos.x - plot_min.x)/
            std::max(1.0f, plot_max.x - plot_min.x), 0.0f, 1.0f);
        const double target_time_coordinate =
            time_min + fraction*(time_max-time_min);
        size_t nearest = first_valid;
        double distance = DBL_MAX;
        for (size_t ordinal = 0; ordinal < series.count; ++ordinal) {
            if (!series_sample_is_plottable(series, ordinal, use_log_log))
                continue;
            const size_t index = series_index(series, ordinal);
            const double candidate = std::fabs(
                plot_coordinate(series.simulation_time[index], use_log_log) -
                target_time_coordinate);
            if (candidate < distance) {
                distance = candidate;
                nearest = ordinal;
            }
        }
        const size_t index = series_index(series, nearest);
        ImGui::BeginTooltip();
        ImGui::Text("t = %.9g", series.simulation_time[index]);
        ImGui::Text("%s = %.9g", label, series.value[index]);
        ImGui::EndTooltip();
    }
    ImGui::PopID();
}

bool pair_sample_is_plottable(const LivePairCorrelation &pair, int bin,
                              bool log_log) {
    const double radius = pair.radius[bin];
    const double value = pair.value[bin];
    return std::isfinite(radius) && std::isfinite(value) &&
        (!log_log || (radius > 0.0 && value > 0.0));
}

void draw_pair_plot(const LivePairCorrelation &pair, bool *log_log) {
    ImGui::Text("Averaged radial distribution g(r)");
    ImGui::SameLine();
    ImGui::Checkbox("log-log", log_log);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Use logarithmic radius and g(r) axes. Zero, negative, "
                          "NaN, and infinite bins are skipped.");
    const bool use_log_log = *log_log;
    const float width = std::max(80.0f, ImGui::GetContentRegionAvail().x);
    const ImVec2 size(width, 150.0f);
    const ImVec2 minimum = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("pair-correlation-plot", size);
    const ImVec2 maximum(minimum.x + size.x, minimum.y + size.y);
    ImDrawList *draw = ImGui::GetWindowDrawList();
    draw->AddRectFilled(minimum, maximum, IM_COL32(18, 21, 27, 255), 4.0f);
    draw->AddRect(minimum, maximum, IM_COL32(68, 75, 88, 255), 4.0f);
    if (!pair.valid || pair.bins <= 0 || pair.averages == 0) {
        draw_empty_plot_message(minimum, maximum,
                                "Waiting for first pair sample");
        return;
    }

    int first_valid = -1;
    int valid_count = 0;
    double radius_min = 0.0;
    double radius_max = 0.0;
    double value_min = 0.0;
    double value_max = 0.0;
    for (int bin = 0; bin < pair.bins; ++bin) {
        if (!pair_sample_is_plottable(pair, bin, use_log_log))
            continue;
        const double radius = plot_coordinate(pair.radius[bin], use_log_log);
        const double value = plot_coordinate(pair.value[bin], use_log_log);
        if (first_valid < 0) {
            first_valid = bin;
            radius_min = radius_max = radius;
            value_min = value_max = value;
        } else {
            radius_min = std::min(radius_min, radius);
            radius_max = std::max(radius_max, radius);
            value_min = std::min(value_min, value);
            value_max = std::max(value_max, value);
        }
        ++valid_count;
    }
    if (valid_count == 0) {
        draw_empty_plot_message(minimum, maximum,
            use_log_log ? "Log-log view needs positive finite bins"
                        : "No finite pair-correlation bins are available");
        return;
    }
    if (!(radius_max > radius_min)) radius_max = radius_min + 1.0;
    if (!use_log_log) value_min = std::min(0.0, value_min);
    if (!(value_max > value_min)) {
        const double padding = std::max(1e-9, std::fabs(value_min)*0.08);
        value_min -= padding;
        value_max += padding;
    } else {
        const double padding = (value_max - value_min)*0.08;
        value_min -= padding;
        value_max += padding;
    }
    const ImVec2 plot_min(minimum.x + 8.0f, minimum.y + 8.0f);
    const ImVec2 plot_max(maximum.x - 8.0f, maximum.y - 22.0f);
    const double reference_value = plot_coordinate(1.0, use_log_log);
    if (reference_value >= value_min && reference_value <= value_max) {
        const float reference_y = plot_max.y - static_cast<float>(
            (reference_value-value_min)/(value_max-value_min))*
            (plot_max.y-plot_min.y);
        draw->AddLine(ImVec2(plot_min.x, reference_y),
                      ImVec2(plot_max.x, reference_y),
                      IM_COL32(125, 132, 145, 150), 1.0f);
    }

    auto point_for = [&](int bin) {
        const float x = plot_min.x + static_cast<float>(
            (plot_coordinate(pair.radius[bin], use_log_log)-radius_min)/
            (radius_max-radius_min))*(plot_max.x-plot_min.x);
        const float y = plot_max.y - static_cast<float>(
            (plot_coordinate(pair.value[bin], use_log_log)-value_min)/
            (value_max-value_min))*(plot_max.y-plot_min.y);
        return ImVec2(x, y);
    };
    draw->PushClipRect(plot_min, plot_max, true);
    ImVec2 previous{};
    bool previous_valid = false;
    for (int bin = 0; bin < pair.bins; ++bin) {
        if (!pair_sample_is_plottable(pair, bin, use_log_log)) {
            previous_valid = false;
            continue;
        }
        const ImVec2 current = point_for(bin);
        if (previous_valid)
            draw->AddLine(previous, current,
                          IM_COL32(195, 125, 255, 255), 1.8f);
        else
            draw->AddCircleFilled(current, 2.0f,
                                  IM_COL32(195, 125, 255, 255));
        previous = current;
        previous_valid = true;
    }
    draw->PopClipRect();
    char footer[160];
    if (use_log_log) {
        std::snprintf(footer, sizeof(footer),
                      "log-log   %d positive bins   %llu averages   %d particles",
                      valid_count,
                      static_cast<unsigned long long>(pair.averages),
                      pair.particles_sampled);
    } else {
        std::snprintf(footer, sizeof(footer),
                      "r <= %.4g   %llu averages   %d particles",
                      pair.maximum_radius,
                      static_cast<unsigned long long>(pair.averages),
                      pair.particles_sampled);
    }
    draw->AddText(ImVec2(plot_min.x, maximum.y - 17.0f),
                  IM_COL32(135, 143, 155, 255), footer);

    if (ImGui::IsItemHovered()) {
        const float fraction = std::clamp(
            (ImGui::GetIO().MousePos.x - plot_min.x)/
            std::max(1.0f, plot_max.x - plot_min.x), 0.0f, 1.0f);
        const double target_radius = radius_min +
            fraction*(radius_max-radius_min);
        int nearest = first_valid;
        double distance = DBL_MAX;
        for (int bin = 0; bin < pair.bins; ++bin) {
            if (!pair_sample_is_plottable(pair, bin, use_log_log))
                continue;
            const double candidate = std::fabs(
                plot_coordinate(pair.radius[bin], use_log_log)-target_radius);
            if (candidate < distance) {
                distance = candidate;
                nearest = bin;
            }
        }
        ImGui::BeginTooltip();
        ImGui::Text("r = %.7g", pair.radius[nearest]);
        ImGui::Text("g(r) = %.7g", pair.value[nearest]);
        ImGui::Text("averaged at t = %.7g", pair.simulation_time);
        ImGui::EndTooltip();
    }
}

void checkbox_setting(const char *label, bool *value, const char *help) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::Checkbox(label, value);
    ImGui::TableSetColumnIndex(1);
    ImGui::TextDisabled("%s", help);
}

} // namespace

namespace edmd_gui {

void draw_live_analysis(Context &context, GuiLiveState &live) {
    LiveAnalysisSettings settings{};
    live_analysis_get_settings(&settings);

    ImGui::TextWrapped(
        "Live analysis samples the running state in memory. It never opens an "
        "output stream or changes the event schedule; headless saving remains "
        "controlled by the command-line/configuration pipeline.");

    if (ImGui::Button("Clear plot history"))
        live_analysis_clear_histories();
    ImGui::SameLine();
    const LiveAnalysisData *data = live_analysis_get_data();
    ImGui::TextDisabled("generation %llu",
        static_cast<unsigned long long>(data->generation));

    ImGui::SeparatorText("Time series");
    if (ImGui::BeginTable("live-series-settings", 2,
                          ImGuiTableFlags_RowBg |
                          ImGuiTableFlags_BordersInnerV |
                          ImGuiTableFlags_SizingStretchProp)) {
        checkbox_setting("Energy / particle", &settings.energy,
                         "kinetic + active conservative potentials");
        checkbox_setting("Temperature", &settings.temperature,
                         "translational kinetic temperature");
        checkbox_setting("Pressure", &settings.pressure,
                         "bulk virial or confining-wall impulse");
        checkbox_setting("Mean-squared displacement",
                         &settings.mean_squared_displacement,
                         "center-of-mass drift removed");
        ImGui::EndTable();
    }

    live_analysis_apply_settings(&settings);
    data = live_analysis_get_data();

    if (settings.energy)
        draw_series_plot("energy", "Energy / particle",
                         data->energy_per_particle,
                         IM_COL32(255, 184, 92, 255),
                         &live.energy_log_log);
    if (settings.temperature)
        draw_series_plot("temperature", "Temperature", data->temperature,
                         IM_COL32(92, 196, 255, 255),
                         &live.temperature_log_log, true);
    if (settings.pressure) {
        if (!data->pressure_supported) {
            ImGui::TextColored(ImVec4(0.95f, 0.68f, 0.25f, 1.0f),
                "Pressure plot unavailable for growth, shear, or a center membrane.");
        } else if (data->pressure_warming_up && data->pressure.count == 0) {
            ImGui::TextDisabled(
                "Pressure needs two live samples to form a measurement interval.");
        }
        draw_series_plot("pressure", "Pressure", data->pressure,
                         IM_COL32(110, 224, 158, 255),
                         &live.pressure_log_log);
    }
    if (settings.mean_squared_displacement) {
        if (!data->msd_supported) {
            ImGui::TextColored(ImVec4(0.95f, 0.68f, 0.25f, 1.0f),
                               "MSD is unavailable with Lees-Edwards shear.");
        } else if (data->msd_allocation_failed) {
            ImGui::TextColored(ImVec4(0.95f, 0.35f, 0.30f, 1.0f),
                               "Unable to allocate the MSD reference state.");
        } else if (data->msd_reference_ready) {
            ImGui::TextDisabled("MSD reference t = %.8g",
                                data->msd_reference_time);
        }
        ImGui::SameLine();
        ImGui::BeginDisabled(!data->msd_supported);
        if (ImGui::SmallButton("Reset MSD reference"))
            live_analysis_reset_msd_reference();
        ImGui::EndDisabled();
        draw_series_plot("msd", "MSD", data->mean_squared_displacement,
                         IM_COL32(224, 118, 196, 255),
                         &live.msd_log_log, true);
    }

    ImGui::SeparatorText("Averaged pair correlation");
    if (ImGui::Checkbox("Enable live g(r)", &settings.pair_correlation)) {
        live_analysis_apply_settings(&settings);
    }
    ImGui::TextDisabled(
        "Deterministic particle subsampling bounds the O(M^2) pair work; each "
        "curve is an online average across completed snapshots.");
    ImGui::BeginDisabled(!settings.pair_correlation);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::SliderInt("Radial bins", &settings.pair_bins, 16,
                     LIVE_ANALYSIS_PAIR_MAX_BINS);
    settings.pair_bins = std::clamp<int>(settings.pair_bins, 16,
                                         LIVE_ANALYSIS_PAIR_MAX_BINS);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::SliderInt("Particle cap", &settings.pair_particle_cap, 32, 4096);
    settings.pair_particle_cap = std::clamp(
        settings.pair_particle_cap, 32, 4096);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputDouble("Maximum radius (0 = auto)",
                       &settings.pair_max_radius, 0.0, 0.0, "%.6g");
    if (!std::isfinite(settings.pair_max_radius) ||
        settings.pair_max_radius < 0.0)
        settings.pair_max_radius = 0.0;
    if (ImGui::Button("Reset g(r) average"))
        live_analysis_reset_pair_average();
    ImGui::EndDisabled();
    live_analysis_apply_settings(&settings);
    data = live_analysis_get_data();
    if (settings.pair_correlation && !data->pair_correlation_supported) {
        ImGui::TextColored(ImVec4(0.95f, 0.68f, 0.25f, 1.0f),
            "g(r) is unavailable with Lees-Edwards shear because its "
            "time-dependent image displacement requires a sheared metric.");
    } else if (settings.pair_correlation) {
        draw_pair_plot(data->pair_correlation,
                       &live.pair_correlation_log_log);
    }
    if (context.capabilities.dimensions == 3)
        ImGui::TextDisabled("g(r) uses 3D shell-volume normalization.");
    if (context.applied.geometry.boundary != SIM_BOUNDARY_PERIODIC)
        ImGui::TextDisabled(
            "For confined systems, finite-wall shell truncation is not corrected.");

    ImGui::SeparatorText("Structure factor");
    ImGui::Checkbox("Live structure-factor window",
                    &live.graphical_structure_factor);
    ImGui::BeginDisabled(!live.graphical_structure_factor);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::SliderInt("Structure-factor particle cap",
                     &live.structure_factor_sample_cap, 1, 8192);
    live.structure_factor_sample_cap = std::clamp(
        live.structure_factor_sample_cap, 1, 8192);
    ImGui::EndDisabled();
    live.structure_factor_sample_cap = std::clamp(
        live.structure_factor_sample_cap, 1, 8192);
    if (context.capabilities.dimensions == 3)
        ImGui::TextDisabled(
            "The existing live structure-factor renderer uses the XY projection.");
    ImGui::TextDisabled(
        "Structure-factor intensity is already shown on a logarithmic color scale; "
        "its signed q axes remain linear.");
}

} // namespace edmd_gui
