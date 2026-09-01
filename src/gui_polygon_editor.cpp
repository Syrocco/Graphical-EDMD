#include "gui_polygon_editor.hpp"

#include "convex_polygon.h"
#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

namespace {

constexpr int kMinVertices = 3;
constexpr double kPi = 3.141592653589793238462643383279502884;

static_assert(SIM_POLYGON_MAX_VERTICES == CONVEX_POLYGON_MAX_VERTICES,
              "GUI and collision geometry vertex limits must agree");

struct Point {
    double x = 0.0;
    double y = 0.0;
};

struct EditablePolygon {
    int count = kMinVertices;
    Point vertices[SIM_POLYGON_MAX_VERTICES]{};
};

enum class DragKind {
    None,
    Vertex,
    Edge
};

struct EditorState {
    const SimParticleShapeConfig *owner = nullptr;
    SimParticleShapeConfig synchronized_shape{};
    EditablePolygon staged{};
    EditablePolygon drag_origin{};
    Point drag_start{};
    DragKind drag_kind = DragKind::None;
    int drag_index = -1;
    bool drag_changed = false;
    bool side_count_active = false;
    bool initialized = false;
};

Point subtract(const Point &a, const Point &b) {
    return {a.x-b.x, a.y-b.y};
}

Point add(const Point &a, const Point &b) {
    return {a.x+b.x, a.y+b.y};
}

Point multiply(const Point &point, double factor) {
    return {point.x*factor, point.y*factor};
}

double squared_distance(const ImVec2 &a, const ImVec2 &b) {
    const double dx = static_cast<double>(a.x-b.x);
    const double dy = static_cast<double>(a.y-b.y);
    return dx*dx + dy*dy;
}

EditablePolygon regular_polygon(int count) {
    EditablePolygon polygon{};
    polygon.count = std::clamp(count, kMinVertices,
                               static_cast<int>(SIM_POLYGON_MAX_VERTICES));
    for (int i = 0; i < polygon.count; ++i) {
        const double angle = kPi*0.5 + 2.0*kPi*static_cast<double>(i)/
            static_cast<double>(polygon.count);
        polygon.vertices[i] = {std::cos(angle), std::sin(angle)};
    }
    return polygon;
}

bool valid_strictly_convex(const EditablePolygon &polygon) {
    if (polygon.count < kMinVertices ||
        polygon.count > static_cast<int>(SIM_POLYGON_MAX_VERTICES)) {
        return false;
    }

    ConvexPolygon geometry{};
    geometry.count = polygon.count;
    for (int i = 0; i < polygon.count; ++i) {
        const Point &point = polygon.vertices[i];
        if (!std::isfinite(point.x) || !std::isfinite(point.y)) return false;
        geometry.vertices[i] = {point.x, point.y};
    }
    return convex_polygon_validate(&geometry) == CONVEX_POLYGON_VALID;
}

bool normalize_polygon(EditablePolygon &polygon) {
    ConvexPolygon geometry{};
    geometry.count = polygon.count;
    for (int i = 0; i < polygon.count; ++i) {
        geometry.vertices[i] = {polygon.vertices[i].x,
                                polygon.vertices[i].y};
    }
    if (!convex_polygon_recenter_and_normalize(&geometry)) return false;
    for (int i = 0; i < polygon.count; ++i) {
        polygon.vertices[i] = {geometry.vertices[i].x,
                               geometry.vertices[i].y};
    }
    return true;
}

EditablePolygon read_polygon(const SimParticleShapeConfig &shape) {
    EditablePolygon polygon{};
    polygon.count = shape.vertex_count;
    if (polygon.count < kMinVertices ||
        polygon.count > static_cast<int>(SIM_POLYGON_MAX_VERTICES)) {
        return regular_polygon(6);
    }
    for (int i = 0; i < polygon.count; ++i) {
        polygon.vertices[i] = {shape.vertices[i][0], shape.vertices[i][1]};
    }
    if (!normalize_polygon(polygon)) return regular_polygon(polygon.count);
    return polygon;
}

bool same_geometry(const SimParticleShapeConfig &shape,
                   const EditablePolygon &polygon) {
    if (shape.vertex_count != polygon.count) return false;
    for (int i = 0; i < polygon.count; ++i) {
        if (std::abs(shape.vertices[i][0]-polygon.vertices[i].x) > 1.0e-12 ||
            std::abs(shape.vertices[i][1]-polygon.vertices[i].y) > 1.0e-12) {
            return false;
        }
    }
    for (int i = polygon.count; i < SIM_POLYGON_MAX_VERTICES; ++i) {
        if (shape.vertices[i][0] != 0.0 || shape.vertices[i][1] != 0.0)
            return false;
    }
    return true;
}

void write_polygon(SimParticleShapeConfig &shape,
                   const EditablePolygon &polygon) {
    shape.vertex_count = polygon.count;
    std::memset(shape.vertices, 0, sizeof(shape.vertices));
    for (int i = 0; i < polygon.count; ++i) {
        shape.vertices[i][0] = polygon.vertices[i].x;
        shape.vertices[i][1] = polygon.vertices[i].y;
    }
}

bool synchronized_geometry_equal(const SimParticleShapeConfig &left,
                                 const SimParticleShapeConfig &right) {
    if (left.vertex_count != right.vertex_count) return false;
    const int count = std::clamp(left.vertex_count, 0,
                                 static_cast<int>(SIM_POLYGON_MAX_VERTICES));
    for (int i = 0; i < count; ++i) {
        if (left.vertices[i][0] != right.vertices[i][0] ||
            left.vertices[i][1] != right.vertices[i][1]) {
            return false;
        }
    }
    return true;
}

bool commit(EditorState &state, SimParticleShapeConfig &shape) {
    EditablePolygon normalized = state.staged;
    if (!normalize_polygon(normalized)) {
        /* This should be unreachable because every staged pointer edit uses
           the same validator. Recover transactionally if floating-point
           roundoff at the exact convexity boundary still defeats the final
           normalization. */
        state.staged = read_polygon(shape);
        state.synchronized_shape = shape;
        return false;
    }
    state.staged = normalized;
    const bool changed = !same_geometry(shape, normalized);
    if (changed) write_polygon(shape, normalized);
    state.synchronized_shape = shape;
    return changed;
}

EditablePolygon interpolate_valid(const EditablePolygon &origin,
                                  const EditablePolygon &requested) {
    if (valid_strictly_convex(requested)) return requested;

    double low = 0.0;
    double high = 1.0;
    EditablePolygon result = origin;
    for (int iteration = 0; iteration < 40; ++iteration) {
        const double fraction = 0.5*(low+high);
        EditablePolygon candidate = origin;
        for (int i = 0; i < candidate.count; ++i) {
            candidate.vertices[i].x = origin.vertices[i].x +
                fraction*(requested.vertices[i].x-origin.vertices[i].x);
            candidate.vertices[i].y = origin.vertices[i].y +
                fraction*(requested.vertices[i].y-origin.vertices[i].y);
        }
        if (valid_strictly_convex(candidate)) {
            low = fraction;
            result = candidate;
        } else {
            high = fraction;
        }
    }
    return result;
}

bool polygons_differ(const EditablePolygon &left,
                     const EditablePolygon &right) {
    if (left.count != right.count) return true;
    for (int i = 0; i < left.count; ++i) {
        if (std::abs(left.vertices[i].x-right.vertices[i].x) > 1.0e-10 ||
            std::abs(left.vertices[i].y-right.vertices[i].y) > 1.0e-10) {
            return true;
        }
    }
    return false;
}

ImVec2 screen_point(const Point &point, const ImVec2 &center, float scale) {
    return ImVec2(center.x+static_cast<float>(point.x)*scale,
                  center.y-static_cast<float>(point.y)*scale);
}

Point local_point(const ImVec2 &point, const ImVec2 &center, float scale) {
    return {
        static_cast<double>(point.x-center.x)/static_cast<double>(scale),
        static_cast<double>(center.y-point.y)/static_cast<double>(scale)
    };
}

int hovered_handle(const EditablePolygon &polygon, const ImVec2 &mouse,
                   const ImVec2 &center, float scale, DragKind &kind) {
    constexpr float vertex_radius = 9.0f;
    constexpr float edge_radius = 7.0f;
    int best = -1;
    double best_distance = static_cast<double>(vertex_radius*vertex_radius);
    kind = DragKind::None;

    /* Vertices deliberately win over edge handles where a very short edge
       makes their hit regions overlap. */
    for (int i = 0; i < polygon.count; ++i) {
        const double distance = squared_distance(
            mouse, screen_point(polygon.vertices[i], center, scale));
        if (distance <= best_distance) {
            best_distance = distance;
            best = i;
            kind = DragKind::Vertex;
        }
    }
    if (best >= 0) return best;

    best_distance = static_cast<double>(edge_radius*edge_radius);
    for (int i = 0; i < polygon.count; ++i) {
        const Point midpoint = multiply(
            add(polygon.vertices[i], polygon.vertices[(i+1)%polygon.count]),
            0.5);
        const double distance = squared_distance(
            mouse, screen_point(midpoint, center, scale));
        if (distance <= best_distance) {
            best_distance = distance;
            best = i;
            kind = DragKind::Edge;
        }
    }
    return best;
}

} // namespace

namespace edmd_gui {

bool draw_polygon_editor(SimParticleShapeConfig &shape) {
    static EditorState state;
    bool changed = false;

    const bool interacting = state.drag_kind != DragKind::None ||
                             state.side_count_active;
    if (!state.initialized || state.owner != &shape ||
        (!interacting &&
         !synchronized_geometry_equal(shape, state.synchronized_shape))) {
        state = EditorState{};
        state.owner = &shape;
        state.staged = read_polygon(shape);
        state.synchronized_shape = shape;
        state.initialized = true;
    }

    ImGui::PushID("convex-polygon-editor");
    ImGui::TextWrapped(
        "Create one normalized convex shape shared by all particles. "
        "Drag a numbered vertex, or drag a smaller edge handle to move both "
        "endpoints. The editor stops at the convexity limit.");

    int requested_count = state.staged.count;
    ImGui::SetNextItemWidth(std::min(240.0f, ImGui::GetContentRegionAvail().x));
    const bool count_changed = ImGui::SliderInt(
        "Sides", &requested_count, kMinVertices,
        static_cast<int>(SIM_POLYGON_MAX_VERTICES), "%d",
        ImGuiSliderFlags_AlwaysClamp);
    const bool count_is_active = ImGui::IsItemActive();
    const bool count_was_deactivated = ImGui::IsItemDeactivated();
    const bool count_was_committed = ImGui::IsItemDeactivatedAfterEdit();
    if (count_changed) {
        state.staged = regular_polygon(requested_count);
        state.side_count_active = count_is_active;
    }
    if (count_was_committed || (count_changed && !count_is_active)) {
        changed |= commit(state, shape);
        state.side_count_active = false;
    } else if (count_is_active) {
        state.side_count_active = true;
    } else if (count_was_deactivated) {
        /* Clicking the slider without moving it is not an edit, but it still
           has to end the staging session. */
        state.side_count_active = false;
    }

    ImGui::SameLine();
    if (ImGui::Button("Reset regular")) {
        state.staged = regular_polygon(state.staged.count);
        changed |= commit(state, shape);
        state.side_count_active = false;
    }
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Restore an equally spaced regular polygon.");

    const float available_width = ImGui::GetContentRegionAvail().x;
    const float canvas_side = std::max(96.0f, std::min(360.0f, available_width));
    const ImVec2 canvas_min = ImGui::GetCursorScreenPos();
    const ImVec2 canvas_size(canvas_side, canvas_side);
    const ImVec2 canvas_max(canvas_min.x+canvas_side, canvas_min.y+canvas_side);
    const ImVec2 center(canvas_min.x+0.5f*canvas_side,
                        canvas_min.y+0.5f*canvas_side);
    const float polygon_scale = 0.36f*canvas_side;

    ImGui::InvisibleButton("shape-canvas", canvas_size,
                           ImGuiButtonFlags_MouseButtonLeft);
    const bool canvas_hovered = ImGui::IsItemHovered();
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    DragKind hover_kind = DragKind::None;
    const int hover_index = canvas_hovered
        ? hovered_handle(state.staged, mouse, center, polygon_scale, hover_kind)
        : -1;

    if (canvas_hovered && hover_index >= 0)
        ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);

    if (ImGui::IsItemActivated() && hover_index >= 0 &&
        !state.side_count_active) {
        state.drag_kind = hover_kind;
        state.drag_index = hover_index;
        state.drag_origin = state.staged;
        state.drag_start = local_point(mouse, center, polygon_scale);
        state.drag_changed = false;
    }

    if (state.drag_kind != DragKind::None &&
        ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        const Point current = local_point(mouse, center, polygon_scale);
        EditablePolygon requested = state.drag_origin;
        if (state.drag_kind == DragKind::Vertex) {
            /* Preserve the point inside the handle that was grabbed. Merely
               pressing/releasing near a handle's edge must not reshape the
               polygon and trigger an unnecessary restart. */
            const Point movement = subtract(current, state.drag_start);
            requested.vertices[state.drag_index] = add(
                requested.vertices[state.drag_index], movement);
        } else {
            const Point movement = subtract(current, state.drag_start);
            const int first = state.drag_index;
            const int second = (first+1)%requested.count;
            requested.vertices[first] = add(requested.vertices[first], movement);
            requested.vertices[second] = add(requested.vertices[second], movement);
        }
        state.staged = interpolate_valid(state.drag_origin, requested);
        state.drag_changed = polygons_differ(state.staged, state.drag_origin);
    }

    if (state.drag_kind != DragKind::None &&
        !ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        if (state.drag_changed) changed |= commit(state, shape);
        else state.staged = state.drag_origin;
        state.drag_kind = DragKind::None;
        state.drag_index = -1;
        state.drag_changed = false;
    }

    ImDrawList *draw_list = ImGui::GetWindowDrawList();
    const ImU32 background = ImGui::GetColorU32(ImGuiCol_FrameBg);
    const ImU32 border = ImGui::GetColorU32(ImGuiCol_Border);
    const ImU32 grid = ImGui::GetColorU32(ImVec4(0.42f, 0.48f, 0.56f, 0.16f));
    const ImU32 fill = ImGui::GetColorU32(ImVec4(0.25f, 0.62f, 0.91f, 0.27f));
    const ImU32 outline = ImGui::GetColorU32(ImVec4(0.34f, 0.72f, 1.0f, 1.0f));
    const ImU32 vertex_color = ImGui::GetColorU32(ImVec4(0.95f, 0.66f, 0.23f, 1.0f));
    const ImU32 edge_color = ImGui::GetColorU32(ImVec4(0.56f, 0.80f, 0.91f, 1.0f));
    const ImU32 active_color = ImGui::GetColorU32(ImGuiCol_ButtonActive);
    const ImU32 text_color = ImGui::GetColorU32(ImGuiCol_Text);

    draw_list->AddRectFilled(canvas_min, canvas_max, background, 5.0f);
    draw_list->AddRect(canvas_min, canvas_max, border, 5.0f);
    draw_list->PushClipRect(canvas_min, canvas_max, true);
    draw_list->AddLine(ImVec2(center.x, canvas_min.y+8.0f),
                       ImVec2(center.x, canvas_max.y-8.0f), grid);
    draw_list->AddLine(ImVec2(canvas_min.x+8.0f, center.y),
                       ImVec2(canvas_max.x-8.0f, center.y), grid);

    ImVec2 screen_vertices[SIM_POLYGON_MAX_VERTICES];
    ImVec2 fill_vertices[SIM_POLYGON_MAX_VERTICES];
    for (int i = 0; i < state.staged.count; ++i) {
        screen_vertices[i] = screen_point(state.staged.vertices[i], center,
                                          polygon_scale);
        /* The physics template is CCW in y-up coordinates.  screen_point()
           flips y, leaving the visual order counter-clockwise; Dear ImGui's
           anti-aliased convex fill requires clockwise screen-space points. */
        fill_vertices[state.staged.count - 1 - i] = screen_vertices[i];
    }
    draw_list->AddConvexPolyFilled(fill_vertices, state.staged.count, fill);
    draw_list->AddPolyline(screen_vertices, state.staged.count, outline,
                           ImDrawFlags_Closed, 2.25f);

    for (int i = 0; i < state.staged.count; ++i) {
        const int next = (i+1)%state.staged.count;
        const ImVec2 midpoint(0.5f*(screen_vertices[i].x+screen_vertices[next].x),
                              0.5f*(screen_vertices[i].y+screen_vertices[next].y));
        const bool active = state.drag_kind == DragKind::Edge &&
                            state.drag_index == i;
        const bool hovered = hover_kind == DragKind::Edge && hover_index == i;
        draw_list->AddCircleFilled(midpoint, active ? 6.0f : (hovered ? 5.5f : 4.5f),
                                   active || hovered ? active_color : edge_color, 12);
    }

    for (int i = 0; i < state.staged.count; ++i) {
        const bool active = state.drag_kind == DragKind::Vertex &&
                            state.drag_index == i;
        const bool hovered = hover_kind == DragKind::Vertex && hover_index == i;
        const float radius = active ? 8.0f : (hovered ? 7.5f : 6.5f);
        draw_list->AddCircleFilled(screen_vertices[i], radius,
                                   active || hovered ? active_color : vertex_color,
                                   16);
        char label[8];
        std::snprintf(label, sizeof(label), "%d", i+1);
        const ImVec2 label_size = ImGui::CalcTextSize(label);
        draw_list->AddText(ImVec2(screen_vertices[i].x-0.5f*label_size.x,
                                  screen_vertices[i].y-0.5f*label_size.y),
                           text_color, label);
    }
    draw_list->PopClipRect();

    ImGui::TextDisabled(
        "Numbered handles reshape vertices; small handles translate edges. "
        "The centroid and outer radius are normalized on release.");
    ImGui::PopID();
    return changed;
}

} // namespace edmd_gui
