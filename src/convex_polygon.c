#include "convex_polygon.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define CONVEX_POLYGON_EDIT_BISECTIONS 64

static ConvexPolygonVec2 vec2(double x, double y)
{
    ConvexPolygonVec2 result;
    result.x = x;
    result.y = y;
    return result;
}

static ConvexPolygonVec2 add2(ConvexPolygonVec2 a, ConvexPolygonVec2 b)
{
    return vec2(a.x + b.x, a.y + b.y);
}

static ConvexPolygonVec2 sub2(ConvexPolygonVec2 a, ConvexPolygonVec2 b)
{
    return vec2(a.x - b.x, a.y - b.y);
}

static ConvexPolygonVec2 mul2(ConvexPolygonVec2 a, double s)
{
    return vec2(a.x * s, a.y * s);
}

static double dot2(ConvexPolygonVec2 a, ConvexPolygonVec2 b)
{
    return a.x * b.x + a.y * b.y;
}

static double cross2(ConvexPolygonVec2 a, ConvexPolygonVec2 b)
{
    return a.x * b.y - a.y * b.x;
}

static double length2(ConvexPolygonVec2 a)
{
    return hypot(a.x, a.y);
}

static int finite_vec2(ConvexPolygonVec2 a)
{
    return isfinite(a.x) && isfinite(a.y);
}

static ConvexPolygonVec2 rotate2(ConvexPolygonVec2 v, double cosine, double sine)
{
    return vec2(cosine * v.x - sine * v.y,
                sine * v.x + cosine * v.y);
}

static ConvexPolygonVec2 angular_velocity_at_point(
    double angular_velocity,
    ConvexPolygonVec2 lever)
{
    return vec2(-angular_velocity * lever.y,
                 angular_velocity * lever.x);
}

static int basic_shape_ok(const ConvexPolygon *polygon)
{
    return polygon != NULL
        && polygon->count >= 3
        && polygon->count <= CONVEX_POLYGON_MAX_VERTICES;
}

static int pose_ok(const ConvexPolygonPose *pose)
{
    return pose != NULL
        && finite_vec2(pose->position)
        && isfinite(pose->angle)
        && isfinite(pose->scale)
        && pose->scale > 0.0;
}

static int motion_ok(const ConvexPolygonMotion *motion)
{
    return motion != NULL
        && pose_ok(&motion->pose)
        && finite_vec2(motion->velocity)
        && isfinite(motion->angular_velocity)
        && isfinite(motion->scale_velocity);
}

static void projection_interval(
    const ConvexPolygonVec2 *vertices,
    int count,
    ConvexPolygonVec2 axis,
    double *minimum,
    double *maximum);

static void compute_properties_unchecked(
    const ConvexPolygon *polygon,
    double *area,
    ConvexPolygonVec2 *centroid,
    double *radius,
    double *inertia_per_mass)
{
    double twice_area = 0.0;
    double centroid_x_numerator = 0.0;
    double centroid_y_numerator = 0.0;
    double polar_numerator = 0.0;
    int i;

    for (i = 0; i < polygon->count; ++i) {
        const int j = (i + 1) % polygon->count;
        const ConvexPolygonVec2 a = polygon->vertices[i];
        const ConvexPolygonVec2 b = polygon->vertices[j];
        const double edge_cross = cross2(a, b);

        twice_area += edge_cross;
        centroid_x_numerator += (a.x + b.x) * edge_cross;
        centroid_y_numerator += (a.y + b.y) * edge_cross;
        polar_numerator += edge_cross
            * (dot2(a, a) + dot2(a, b) + dot2(b, b));
    }

    *area = 0.5 * twice_area;
    centroid->x = centroid_x_numerator / (3.0 * twice_area);
    centroid->y = centroid_y_numerator / (3.0 * twice_area);

    *radius = 0.0;
    for (i = 0; i < polygon->count; ++i) {
        const double r = length2(sub2(polygon->vertices[i], *centroid));
        if (r > *radius) {
            *radius = r;
        }
    }

    /* Integral r^2 dA about the origin, followed by the parallel-axis shift. */
    {
        double polar_about_centroid = polar_numerator / 12.0
            - (*area) * dot2(*centroid, *centroid);
        const double roundoff_scale = DBL_EPSILON
            * (1.0 + fabs(polar_numerator) + fabs(*area));
        if (polar_about_centroid < 0.0
            && polar_about_centroid >= -64.0 * roundoff_scale) {
            polar_about_centroid = 0.0;
        }
        *inertia_per_mass = polar_about_centroid / (*area);
    }
}

ConvexPolygonValidation convex_polygon_validate_with_tolerances(
    const ConvexPolygon *polygon,
    double min_edge_fraction,
    double min_turn_sine)
{
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double scale;
    double minimum_edge;
    double twice_area = 0.0;
    double geometric_epsilon;
    int i;
    int j;

    if (polygon == NULL
        || !isfinite(min_edge_fraction)
        || !isfinite(min_turn_sine)
        || min_edge_fraction < 0.0
        || min_turn_sine < 0.0
        || min_turn_sine >= 1.0) {
        return CONVEX_POLYGON_INVALID_ARGUMENT;
    }
    if (polygon->count < 3
        || polygon->count > CONVEX_POLYGON_MAX_VERTICES) {
        return CONVEX_POLYGON_INVALID_VERTEX_COUNT;
    }

    for (i = 0; i < polygon->count; ++i) {
        if (!finite_vec2(polygon->vertices[i])) {
            return CONVEX_POLYGON_NONFINITE_VERTEX;
        }
    }

    min_x = max_x = polygon->vertices[0].x;
    min_y = max_y = polygon->vertices[0].y;
    for (i = 1; i < polygon->count; ++i) {
        const ConvexPolygonVec2 v = polygon->vertices[i];
        if (v.x < min_x) min_x = v.x;
        if (v.x > max_x) max_x = v.x;
        if (v.y < min_y) min_y = v.y;
        if (v.y > max_y) max_y = v.y;
    }
    scale = hypot(max_x - min_x, max_y - min_y);
    if (!isfinite(scale) || scale <= DBL_MIN) {
        return CONVEX_POLYGON_DEGENERATE;
    }

    minimum_edge = min_edge_fraction * scale;
    geometric_epsilon = 128.0 * DBL_EPSILON * scale * scale;

    for (i = 0; i < polygon->count; ++i) {
        const int previous = (i + polygon->count - 1) % polygon->count;
        const int next = (i + 1) % polygon->count;
        const ConvexPolygonVec2 incoming = sub2(
            polygon->vertices[i], polygon->vertices[previous]);
        const ConvexPolygonVec2 outgoing = sub2(
            polygon->vertices[next], polygon->vertices[i]);
        const double incoming_length = length2(incoming);
        const double outgoing_length = length2(outgoing);
        const double turn = cross2(incoming, outgoing);

        if (incoming_length <= minimum_edge || outgoing_length <= minimum_edge) {
            return CONVEX_POLYGON_EDGE_TOO_SHORT;
        }
        if (turn <= min_turn_sine * incoming_length * outgoing_length) {
            return CONVEX_POLYGON_NOT_STRICTLY_CCW_CONVEX;
        }
        twice_area += cross2(polygon->vertices[i], polygon->vertices[next]);
    }

    if (!isfinite(twice_area) || twice_area <= geometric_epsilon) {
        return CONVEX_POLYGON_DEGENERATE;
    }

    /*
     * Positive local turns alone can admit a multiply-wound/self-intersecting
     * ordering.  A convex CCW polygon has every other vertex on the left of
     * every directed edge.
     */
    for (i = 0; i < polygon->count; ++i) {
        const int next = (i + 1) % polygon->count;
        const ConvexPolygonVec2 origin = polygon->vertices[i];
        const ConvexPolygonVec2 edge = sub2(polygon->vertices[next], origin);
        for (j = 0; j < polygon->count; ++j) {
            if (j == i || j == next) {
                continue;
            }
            if (cross2(edge, sub2(polygon->vertices[j], origin))
                < -geometric_epsilon) {
                return CONVEX_POLYGON_NOT_STRICTLY_CCW_CONVEX;
            }
        }
    }

    return CONVEX_POLYGON_VALID;
}

ConvexPolygonValidation convex_polygon_validate(const ConvexPolygon *polygon)
{
    return convex_polygon_validate_with_tolerances(
        polygon,
        CONVEX_POLYGON_DEFAULT_MIN_EDGE_FRACTION,
        CONVEX_POLYGON_DEFAULT_MIN_TURN_SINE);
}

const char *convex_polygon_validation_message(ConvexPolygonValidation result)
{
    switch (result) {
    case CONVEX_POLYGON_VALID:
        return "valid";
    case CONVEX_POLYGON_INVALID_ARGUMENT:
        return "invalid polygon or tolerance";
    case CONVEX_POLYGON_INVALID_VERTEX_COUNT:
        return "a polygon needs between 3 and 16 vertices";
    case CONVEX_POLYGON_NONFINITE_VERTEX:
        return "polygon coordinates must be finite";
    case CONVEX_POLYGON_DEGENERATE:
        return "polygon area is degenerate";
    case CONVEX_POLYGON_EDGE_TOO_SHORT:
        return "polygon edge is too short";
    case CONVEX_POLYGON_NOT_STRICTLY_CCW_CONVEX:
        return "vertices must form a strictly convex counter-clockwise polygon";
    default:
        return "unknown polygon validation result";
    }
}

int convex_polygon_refresh(ConvexPolygon *polygon)
{
    int i;
    if (convex_polygon_validate(polygon) != CONVEX_POLYGON_VALID) {
        return 0;
    }
    compute_properties_unchecked(
        polygon,
        &polygon->area,
        &polygon->centroid,
        &polygon->radius,
        &polygon->inertia_per_mass);
    for (i = 0; i < polygon->count; ++i) {
        const int next = i + 1 == polygon->count ? 0 : i + 1;
        const ConvexPolygonVec2 edge = sub2(
            polygon->vertices[next], polygon->vertices[i]);
        const double edge_length = length2(edge);
        double minimum;
        double maximum;
        if (!(edge_length > DBL_MIN) || !isfinite(edge_length)) return 0;
        polygon->edge_inverse_lengths[i] = 1.0/edge_length;
        polygon->edge_normals[i] = vec2(
            edge.y*polygon->edge_inverse_lengths[i],
            -edge.x*polygon->edge_inverse_lengths[i]);
        projection_interval(
            polygon->vertices, polygon->count, polygon->edge_normals[i],
            &minimum, &maximum);
        polygon->edge_projection_minima[i] = minimum;
        polygon->edge_projection_maxima[i] = maximum;
    }
    return isfinite(polygon->area)
        && isfinite(polygon->centroid.x)
        && isfinite(polygon->centroid.y)
        && isfinite(polygon->radius)
        && isfinite(polygon->inertia_per_mass)
        && polygon->area > 0.0
        && polygon->radius > 0.0
        && polygon->inertia_per_mass > 0.0;
}

int convex_polygon_recenter_and_normalize(ConvexPolygon *polygon)
{
    ConvexPolygonVec2 centroid;
    double area;
    double radius;
    double inertia;
    int i;

    if (convex_polygon_validate(polygon) != CONVEX_POLYGON_VALID) {
        return 0;
    }
    compute_properties_unchecked(polygon, &area, &centroid, &radius, &inertia);
    if (!isfinite(radius) || radius <= DBL_MIN) {
        return 0;
    }

    for (i = 0; i < polygon->count; ++i) {
        polygon->vertices[i] = mul2(
            sub2(polygon->vertices[i], centroid), 1.0 / radius);
    }
    if (convex_polygon_validate(polygon) != CONVEX_POLYGON_VALID) {
        return 0;
    }
    if (!convex_polygon_refresh(polygon)) {
        return 0;
    }
    /* Keep the cached broadphase radius conservative despite roundoff. */
    if (polygon->radius < 1.0) {
        polygon->radius = 1.0;
    }
    return 1;
}

int convex_polygon_make_regular(ConvexPolygon *polygon, int vertex_count)
{
    int i;

    if (polygon == NULL
        || vertex_count < 3
        || vertex_count > CONVEX_POLYGON_MAX_VERTICES) {
        return 0;
    }
    memset(polygon, 0, sizeof(*polygon));
    polygon->count = vertex_count;
    for (i = 0; i < vertex_count; ++i) {
        const double angle = 0.5 * M_PI
            + 2.0 * M_PI * (double)i / (double)vertex_count;
        polygon->vertices[i] = vec2(cos(angle), sin(angle));
    }
    return convex_polygon_recenter_and_normalize(polygon);
}

static int finite_edit_target(ConvexPolygonVec2 value)
{
    return isfinite(value.x) && isfinite(value.y);
}

double convex_polygon_move_vertex_clamped(
    ConvexPolygon *polygon,
    int vertex_index,
    ConvexPolygonVec2 requested_position)
{
    ConvexPolygon original;
    ConvexPolygon candidate;
    ConvexPolygonVec2 start;
    double low;
    double high;
    int iteration;

    if (polygon == NULL
        || vertex_index < 0
        || vertex_index >= polygon->count
        || !finite_edit_target(requested_position)
        || convex_polygon_validate(polygon) != CONVEX_POLYGON_VALID) {
        return -1.0;
    }

    original = *polygon;
    start = original.vertices[vertex_index];
    candidate = original;
    candidate.vertices[vertex_index] = requested_position;
    if (convex_polygon_validate(&candidate) == CONVEX_POLYGON_VALID) {
        *polygon = candidate;
        (void)convex_polygon_refresh(polygon);
        return 1.0;
    }

    low = 0.0;
    high = 1.0;
    for (iteration = 0; iteration < CONVEX_POLYGON_EDIT_BISECTIONS; ++iteration) {
        const double middle = 0.5 * (low + high);
        candidate = original;
        candidate.vertices[vertex_index] = add2(
            start, mul2(sub2(requested_position, start), middle));
        if (convex_polygon_validate(&candidate) == CONVEX_POLYGON_VALID) {
            low = middle;
        } else {
            high = middle;
        }
    }

    candidate = original;
    candidate.vertices[vertex_index] = add2(
        start, mul2(sub2(requested_position, start), low));
    *polygon = candidate;
    (void)convex_polygon_refresh(polygon);
    return low;
}

double convex_polygon_translate_edge_clamped(
    ConvexPolygon *polygon,
    int edge_index,
    ConvexPolygonVec2 requested_delta)
{
    ConvexPolygon original;
    ConvexPolygon candidate;
    const int next_index = polygon != NULL && polygon->count > 0
        ? (edge_index + 1) % polygon->count : 0;
    double low;
    double high;
    int iteration;

    if (polygon == NULL
        || edge_index < 0
        || edge_index >= polygon->count
        || !finite_edit_target(requested_delta)
        || convex_polygon_validate(polygon) != CONVEX_POLYGON_VALID) {
        return -1.0;
    }

    original = *polygon;
    candidate = original;
    candidate.vertices[edge_index] = add2(
        original.vertices[edge_index], requested_delta);
    candidate.vertices[next_index] = add2(
        original.vertices[next_index], requested_delta);
    if (convex_polygon_validate(&candidate) == CONVEX_POLYGON_VALID) {
        *polygon = candidate;
        (void)convex_polygon_refresh(polygon);
        return 1.0;
    }

    low = 0.0;
    high = 1.0;
    for (iteration = 0; iteration < CONVEX_POLYGON_EDIT_BISECTIONS; ++iteration) {
        const double middle = 0.5 * (low + high);
        const ConvexPolygonVec2 accepted_delta = mul2(requested_delta, middle);
        candidate = original;
        candidate.vertices[edge_index] = add2(
            original.vertices[edge_index], accepted_delta);
        candidate.vertices[next_index] = add2(
            original.vertices[next_index], accepted_delta);
        if (convex_polygon_validate(&candidate) == CONVEX_POLYGON_VALID) {
            low = middle;
        } else {
            high = middle;
        }
    }

    candidate = original;
    candidate.vertices[edge_index] = add2(
        original.vertices[edge_index], mul2(requested_delta, low));
    candidate.vertices[next_index] = add2(
        original.vertices[next_index], mul2(requested_delta, low));
    *polygon = candidate;
    (void)convex_polygon_refresh(polygon);
    return low;
}

static int transform_vertices(
    const ConvexPolygon *polygon,
    const ConvexPolygonPose *pose,
    ConvexPolygonVec2 transformed[CONVEX_POLYGON_MAX_VERTICES],
    double *cosine_out,
    double *sine_out)
{
    double cosine;
    double sine;
    int i;

    if (!basic_shape_ok(polygon) || !pose_ok(pose)) {
        return 0;
    }
    {
        /* Polygon state keeps ordinary angles within this range. Avoid the
           comparatively expensive general argument reduction on the hot SAT
           path; retain it for long prediction horizons and external callers. */
        const double reduced_angle = fabs(pose->angle) <= 8.0*M_PI
            ? pose->angle : remainder(pose->angle, 2.0*M_PI);
        cosine = cos(reduced_angle);
        sine = sin(reduced_angle);
    }
    if (cosine_out != NULL) *cosine_out = cosine;
    if (sine_out != NULL) *sine_out = sine;
    for (i = 0; i < polygon->count; ++i) {
        const ConvexPolygonVec2 local = mul2(polygon->vertices[i], pose->scale);
        transformed[i] = add2(pose->position, rotate2(local, cosine, sine));
        if (!finite_vec2(transformed[i])) {
            return 0;
        }
    }
    return 1;
}

static void projection_interval(
    const ConvexPolygonVec2 *vertices,
    int count,
    ConvexPolygonVec2 axis,
    double *minimum,
    double *maximum)
{
    int i;
    *minimum = *maximum = dot2(vertices[0], axis);
    for (i = 1; i < count; ++i) {
        const double projection = dot2(vertices[i], axis);
        if (projection < *minimum) *minimum = projection;
        if (projection > *maximum) *maximum = projection;
    }
}

static ConvexPolygonVec2 averaged_support(
    const ConvexPolygonVec2 *vertices,
    int count,
    ConvexPolygonVec2 direction,
    int maximum_support)
{
    double extreme = dot2(vertices[0], direction);
    double tolerance;
    ConvexPolygonVec2 sum = vec2(0.0, 0.0);
    int matches = 0;
    int i;

    for (i = 1; i < count; ++i) {
        const double projection = dot2(vertices[i], direction);
        if ((maximum_support && projection > extreme)
            || (!maximum_support && projection < extreme)) {
            extreme = projection;
        }
    }
    tolerance = 256.0 * DBL_EPSILON * (1.0 + fabs(extreme));
    for (i = 0; i < count; ++i) {
        const double projection = dot2(vertices[i], direction);
        if ((maximum_support && projection >= extreme - tolerance)
            || (!maximum_support && projection <= extreme + tolerance)) {
            sum = add2(sum, vertices[i]);
            ++matches;
        }
    }
    return mul2(sum, 1.0 / (double)matches);
}

static void support_feature_interval(
    const ConvexPolygonVec2 *vertices,
    int count,
    ConvexPolygonVec2 normal,
    int maximum_support,
    double *normal_coordinate,
    double *minimum_tangent,
    double *maximum_tangent)
{
    const ConvexPolygonVec2 tangent = vec2(-normal.y, normal.x);
    double extreme = dot2(vertices[0], normal);
    double coordinate_scale = fabs(extreme);
    int i;

    for (i = 1; i < count; ++i) {
        const double projection = dot2(vertices[i], normal);
        if ((maximum_support && projection > extreme)
            || (!maximum_support && projection < extreme)) {
            extreme = projection;
        }
        coordinate_scale = fmax(coordinate_scale,
                                fabs(vertices[i].x) + fabs(vertices[i].y));
    }

    const double tolerance = 512.0 * DBL_EPSILON
        * (1.0 + coordinate_scale + fabs(extreme));
    *minimum_tangent = DBL_MAX;
    *maximum_tangent = -DBL_MAX;
    for (i = 0; i < count; ++i) {
        const double projection = dot2(vertices[i], normal);
        if ((maximum_support && projection >= extreme - tolerance)
            || (!maximum_support && projection <= extreme + tolerance)) {
            const double tangent_coordinate = dot2(vertices[i], tangent);
            if (tangent_coordinate < *minimum_tangent)
                *minimum_tangent = tangent_coordinate;
            if (tangent_coordinate > *maximum_tangent)
                *maximum_tangent = tangent_coordinate;
        }
    }
    *normal_coordinate = extreme;
}

typedef struct SatAxisResult {
    double separation;
    ConvexPolygonVec2 normal;
    int source;
    int edge;
} SatAxisResult;

/* Shared SAT axis scan. Keeping witness construction out of this routine lets
 * conservative advancement ask only for separation on its overwhelmingly
 * common non-contact iterations. */
static int sat_axis_scan(
    const ConvexPolygon *a,
    const ConvexPolygonPose *pose_a,
    const ConvexPolygonVec2 vertices_a[CONVEX_POLYGON_MAX_VERTICES],
    double cosine_a,
    double sine_a,
    const ConvexPolygon *b,
    const ConvexPolygonPose *pose_b,
    const ConvexPolygonVec2 vertices_b[CONVEX_POLYGON_MAX_VERTICES],
    double cosine_b,
    double sine_b,
    SatAxisResult *result)
{
    const ConvexPolygonVec2 center_delta =
        sub2(pose_b->position, pose_a->position);
    double best_separation = -DBL_MAX;
    double best_alignment = -DBL_MAX;
    ConvexPolygonVec2 best_normal = vec2(1.0, 0.0);
    int best_source = 0;
    int best_edge = 0;
    int source;
    int i;

    for (source = 0; source < 2; ++source) {
        const int source_count = source == 0 ? a->count : b->count;
        const ConvexPolygon *source_polygon = source == 0 ? a : b;
        const double source_cosine = source == 0 ? cosine_a : cosine_b;
        const double source_sine = source == 0 ? sine_a : sine_b;
        for (i = 0; i < source_count; ++i) {
            ConvexPolygonVec2 axis;
            ConvexPolygonVec2 candidate_normal;
            double min_a;
            double max_a;
            double min_b;
            double max_b;
            double gap_a_to_b;
            double gap_b_to_a;
            double candidate_separation;
            double alignment;
            double comparison_tolerance;

            if (!finite_vec2(source_polygon->edge_normals[i])
                || (source_polygon->edge_normals[i].x == 0.0
                    && source_polygon->edge_normals[i].y == 0.0)) {
                return 0;
            }
            axis = rotate2(source_polygon->edge_normals[i],
                           source_cosine, source_sine);
            if (source == 0) {
                const double origin_projection = dot2(pose_a->position, axis);
                min_a = origin_projection + pose_a->scale*
                    source_polygon->edge_projection_minima[i];
                max_a = origin_projection + pose_a->scale*
                    source_polygon->edge_projection_maxima[i];
                projection_interval(vertices_b, b->count, axis,
                                    &min_b, &max_b);
            } else {
                const double origin_projection = dot2(pose_b->position, axis);
                min_b = origin_projection + pose_b->scale*
                    source_polygon->edge_projection_minima[i];
                max_b = origin_projection + pose_b->scale*
                    source_polygon->edge_projection_maxima[i];
                projection_interval(vertices_a, a->count, axis,
                                    &min_a, &max_a);
            }
            gap_a_to_b = min_b - max_a;
            gap_b_to_a = min_a - max_b;

            if (gap_a_to_b > gap_b_to_a) {
                candidate_normal = axis;
                candidate_separation = gap_a_to_b;
            } else if (gap_b_to_a > gap_a_to_b) {
                candidate_normal = mul2(axis, -1.0);
                candidate_separation = gap_b_to_a;
            } else if (dot2(center_delta, axis) >= 0.0) {
                candidate_normal = axis;
                candidate_separation = gap_a_to_b;
            } else {
                candidate_normal = mul2(axis, -1.0);
                candidate_separation = gap_b_to_a;
            }

            alignment = dot2(candidate_normal, center_delta);
            comparison_tolerance = 128.0 * DBL_EPSILON
                * (1.0 + fabs(best_separation) + fabs(candidate_separation));
            if (candidate_separation > best_separation + comparison_tolerance
                || (fabs(candidate_separation - best_separation)
                        <= comparison_tolerance
                    && alignment > best_alignment)) {
                best_separation = candidate_separation;
                best_alignment = alignment;
                best_normal = candidate_normal;
                best_source = source;
                best_edge = i;
            }
        }
    }

    result->separation = best_separation;
    result->normal = best_normal;
    result->source = best_source;
    result->edge = best_edge;
    return isfinite(best_separation) && finite_vec2(best_normal);
}

static int convex_polygon_sat_separation(
    const ConvexPolygon *a,
    const ConvexPolygonPose *pose_a,
    const ConvexPolygon *b,
    const ConvexPolygonPose *pose_b,
    double *separation)
{
    ConvexPolygonVec2 vertices_a[CONVEX_POLYGON_MAX_VERTICES];
    ConvexPolygonVec2 vertices_b[CONVEX_POLYGON_MAX_VERTICES];
    SatAxisResult axes;
    double cosine_a;
    double sine_a;
    double cosine_b;
    double sine_b;
    if (separation == NULL
        || !transform_vertices(a, pose_a, vertices_a, &cosine_a, &sine_a)
        || !transform_vertices(b, pose_b, vertices_b, &cosine_b, &sine_b)
        || !sat_axis_scan(a, pose_a, vertices_a, cosine_a, sine_a,
                          b, pose_b, vertices_b, cosine_b, sine_b, &axes)) {
        return 0;
    }
    *separation = axes.separation;
    return 1;
}

int convex_polygon_sat_contact(
    const ConvexPolygon *a,
    const ConvexPolygonPose *pose_a,
    const ConvexPolygon *b,
    const ConvexPolygonPose *pose_b,
    ConvexPolygonContact *contact)
{
    ConvexPolygonVec2 vertices_a[CONVEX_POLYGON_MAX_VERTICES];
    ConvexPolygonVec2 vertices_b[CONVEX_POLYGON_MAX_VERTICES];
    SatAxisResult axes;
    double cosine_a;
    double sine_a;
    double cosine_b;
    double sine_b;

    if (contact == NULL
        || !transform_vertices(a, pose_a, vertices_a, &cosine_a, &sine_a)
        || !transform_vertices(b, pose_b, vertices_b, &cosine_b, &sine_b)
        || !sat_axis_scan(a, pose_a, vertices_a, cosine_a, sine_a,
                          b, pose_b, vertices_b, cosine_b, sine_b, &axes)) {
        return 0;
    }

    contact->separation = axes.separation;
    contact->normal = axes.normal;
    {
        const ConvexPolygonVec2 tangent = vec2(-axes.normal.y,
                                                axes.normal.x);
        double normal_a;
        double normal_b;
        double tangent_min_a;
        double tangent_max_a;
        double tangent_min_b;
        double tangent_max_b;
        support_feature_interval(vertices_a, a->count, axes.normal, 1,
                                 &normal_a, &tangent_min_a, &tangent_max_a);
        support_feature_interval(vertices_b, b->count, axes.normal, 0,
                                 &normal_b, &tangent_min_b, &tangent_max_b);
        const double lower = fmax(tangent_min_a, tangent_min_b);
        const double upper = fmin(tangent_max_a, tangent_max_b);
        /* Midpoint of the overlapping support features.  For a vertex/face
           contact this selects the vertex coordinate, whereas averaging each
           complete support independently can place the witness away from the
           true contact and give CCD the wrong angular approach velocity. */
        const double tangent_coordinate = 0.5 * (lower + upper);
        contact->point_a = add2(mul2(axes.normal, normal_a),
                                mul2(tangent, tangent_coordinate));
        contact->point_b = add2(mul2(axes.normal, normal_b),
                                mul2(tangent, tangent_coordinate));
        contact->point = mul2(add2(contact->point_a, contact->point_b), 0.5);
    }
    contact->axis_from_a = axes.source == 0;
    contact->edge_index = axes.edge;
    return isfinite(contact->separation)
        && finite_vec2(contact->normal)
        && finite_vec2(contact->point);
}

void convex_polygon_ccd_default_settings(ConvexPolygonCcdSettings *settings)
{
    if (settings != NULL) {
        settings->distance_tolerance = 0.0;
        settings->time_tolerance = 0.0;
        settings->max_iterations = 256U;
    }
}

double convex_polygon_ccd_velocity_tolerance(double speed_bound)
{
    if (!isfinite(speed_bound) || speed_bound < 0.0) return NAN;
    return 256.0 * DBL_EPSILON * (1.0 + speed_bound);
}

static ConvexPolygonPose pose_at_time(
    const ConvexPolygonMotion *motion,
    double time)
{
    ConvexPolygonPose pose = motion->pose;
    pose.position = add2(pose.position, mul2(motion->velocity, time));
    pose.angle += motion->angular_velocity * time;
    pose.scale += motion->scale_velocity * time;
    return pose;
}

static void resolved_ccd_settings(
    const ConvexPolygonCcdSettings *requested,
    double length_scale,
    double horizon,
    double *distance_tolerance,
    double *time_tolerance,
    unsigned int *max_iterations)
{
    *distance_tolerance = requested != NULL
            && isfinite(requested->distance_tolerance)
            && requested->distance_tolerance > 0.0
        ? requested->distance_tolerance
        : 1024.0 * DBL_EPSILON * fmax(1.0, length_scale);
    *time_tolerance = requested != NULL
            && isfinite(requested->time_tolerance)
            && requested->time_tolerance > 0.0
        ? requested->time_tolerance
        : 128.0 * DBL_EPSILON * fmax(1.0, horizon);
    *max_iterations = requested != NULL && requested->max_iterations > 0U
        ? requested->max_iterations : 256U;
}

/* Intersect [0,horizon] with
 *
 *   |relative_position + relative_velocity*t|
 *       <= initial_radius
 *          + (final_radius-initial_radius)*t/horizon.
 *
 * We solve in normalized time x=t/horizon and normalized length, which keeps
 * the quadratic coefficients bounded even for very large horizons.  Both
 * endpoint radii are supplied because they have already been checked by the
 * caller.  The feasible set is an interval: norm(affine)-affine is convex
 * while the radius remains positive. */
static int circumscribed_disk_interval(
    ConvexPolygonVec2 relative_position,
    ConvexPolygonVec2 relative_velocity,
    double initial_radius,
    double final_radius,
    double horizon,
    double *interval_start,
    double *interval_end)
{
    ConvexPolygonVec2 displacement;
    double length_scale;
    double px;
    double py;
    double dx;
    double dy;
    double radius_start;
    double radius_delta;
    double coefficient_a;
    double coefficient_b;
    double coefficient_c;
    double first = 0.0;
    double second = 1.0;
    int has_interval = 0;

    if (interval_start == NULL || interval_end == NULL
        || !finite_vec2(relative_position)
        || !finite_vec2(relative_velocity)
        || !isfinite(initial_radius) || !(initial_radius > 0.0)
        || !isfinite(final_radius) || !(final_radius > 0.0)
        || !isfinite(horizon) || horizon < 0.0) {
        return 0;
    }
    if (horizon == 0.0) {
        const double center_distance = length2(relative_position);
        if (!isfinite(center_distance) || center_distance > initial_radius)
            return 0;
        *interval_start = 0.0;
        *interval_end = 0.0;
        return 1;
    }

    displacement = mul2(relative_velocity, horizon);
    if (!finite_vec2(displacement)) {
        /* An unsafe broadphase calculation must never become a false NO_HIT.
           Let conservative advancement diagnose the extreme motion. */
        *interval_start = 0.0;
        *interval_end = horizon;
        return 1;
    }
    length_scale = fmax(fmax(fabs(relative_position.x),
                             fabs(relative_position.y)),
                        fmax(fabs(displacement.x), fabs(displacement.y)));
    length_scale = fmax(length_scale,
                        fmax(initial_radius, final_radius));
    if (!isfinite(length_scale) || !(length_scale > DBL_MIN)) {
        *interval_start = 0.0;
        *interval_end = horizon;
        return 1;
    }

    px = relative_position.x/length_scale;
    py = relative_position.y/length_scale;
    dx = displacement.x/length_scale;
    dy = displacement.y/length_scale;
    radius_start = initial_radius/length_scale;
    radius_delta = (final_radius - initial_radius)/length_scale;
    {
        /* Inflate both endpoint radii by a constant normalized ulp guard.
           The resulting interval may be microscopically wider, never narrower
           than the true disk-overlap interval. */
        const double roundoff_padding = 256.0*DBL_EPSILON*
            (1.0 + fabs(px) + fabs(py) + fabs(dx) + fabs(dy) +
             fabs(radius_start) + fabs(radius_delta));
        radius_start += roundoff_padding;
    }

    coefficient_a = dx*dx + dy*dy - radius_delta*radius_delta;
    coefficient_b = 2.0*(px*dx + py*dy -
                         radius_start*radius_delta);
    coefficient_c = px*px + py*py - radius_start*radius_start;

    if (coefficient_a == 0.0) {
        if (coefficient_b == 0.0) {
            has_interval = coefficient_c <= 0.0;
        } else {
            const double root = -coefficient_c/coefficient_b;
            if (coefficient_b > 0.0) {
                first = 0.0;
                second = fmin(1.0, root);
                has_interval = second >= 0.0;
            } else {
                first = fmax(0.0, root);
                second = 1.0;
                has_interval = first <= 1.0;
            }
        }
    } else {
        double discriminant = coefficient_b*coefficient_b -
            4.0*coefficient_a*coefficient_c;
        const double discriminant_roundoff = 256.0*DBL_EPSILON*
            (1.0 + coefficient_b*coefficient_b +
             fabs(4.0*coefficient_a*coefficient_c));
        if (discriminant < 0.0 &&
            discriminant >= -discriminant_roundoff)
            discriminant = 0.0;
        if (discriminant < 0.0) {
            has_interval = coefficient_a < 0.0;
        } else {
            const double square_root = sqrt(discriminant);
            double low_root;
            double high_root;
            if (square_root == 0.0) {
                low_root = high_root =
                    -0.5*coefficient_b/coefficient_a;
            } else {
                const double q = -0.5*(coefficient_b +
                    copysign(square_root, coefficient_b));
                low_root = q/coefficient_a;
                high_root = coefficient_c/q;
                if (low_root > high_root) {
                    const double swap = low_root;
                    low_root = high_root;
                    high_root = swap;
                }
            }
            if (coefficient_a > 0.0) {
                first = fmax(0.0, low_root);
                second = fmin(1.0, high_root);
                has_interval = first <= second;
            } else {
                const int left_exists = low_root >= 0.0;
                const int right_exists = high_root <= 1.0;
                if (left_exists && right_exists) {
                    /* This cannot occur for positive linearly varying radii
                       in exact arithmetic (the feasible set is convex). A
                       full interval is the safe response to roundoff. */
                    first = 0.0;
                    second = 1.0;
                    has_interval = 1;
                } else if (left_exists) {
                    first = 0.0;
                    second = fmin(1.0, low_root);
                    has_interval = second >= 0.0;
                } else if (right_exists) {
                    first = fmax(0.0, high_root);
                    second = 1.0;
                    has_interval = first <= 1.0;
                }
            }
        }
    }
    if (!has_interval || !isfinite(first) || !isfinite(second)) return 0;

    *interval_start = first <= 0.0 ? 0.0 :
        fmax(0.0, nextafter(first*horizon, -INFINITY));
    *interval_end = second >= 1.0 ? horizon :
        fmin(horizon, nextafter(second*horizon, INFINITY));
    return *interval_start <= *interval_end;
}

static double motion_contact_normal_velocity(
    const ConvexPolygonMotion *a,
    const ConvexPolygonPose *pose_a,
    const ConvexPolygonMotion *b,
    const ConvexPolygonPose *pose_b,
    const ConvexPolygonContact *contact)
{
    const ConvexPolygonVec2 lever_a = sub2(contact->point, pose_a->position);
    const ConvexPolygonVec2 lever_b = sub2(contact->point, pose_b->position);
    ConvexPolygonVec2 velocity_a = add2(
        a->velocity,
        angular_velocity_at_point(a->angular_velocity, lever_a));
    ConvexPolygonVec2 velocity_b = add2(
        b->velocity,
        angular_velocity_at_point(b->angular_velocity, lever_b));
    if (pose_a->scale > DBL_MIN && a->scale_velocity != 0.0)
        velocity_a = add2(velocity_a, mul2(
            lever_a, a->scale_velocity/pose_a->scale));
    if (pose_b->scale > DBL_MIN && b->scale_velocity != 0.0)
        velocity_b = add2(velocity_b, mul2(
            lever_b, b->scale_velocity/pose_b->scale));
    return dot2(sub2(velocity_b, velocity_a), contact->normal);
}

static void store_ccd_result(
    ConvexPolygonCcdResult *result,
    ConvexPolygonCcdStatus status,
    double time,
    unsigned int iterations,
    const ConvexPolygonContact *contact)
{
    result->status = status;
    result->time = time;
    result->iterations = iterations;
    if (contact != NULL) {
        result->contact = *contact;
    } else {
        memset(&result->contact, 0, sizeof(result->contact));
    }
}

static int bisect_pair_contact(
    const ConvexPolygon *a,
    const ConvexPolygonMotion *motion_a,
    const ConvexPolygon *b,
    const ConvexPolygonMotion *motion_b,
    double low,
    double high,
    double distance_tolerance,
    double time_tolerance,
    unsigned int max_iterations,
    unsigned int *iterations,
    double *hit_time,
    ConvexPolygonContact *hit_contact)
{
    while (high - low > time_tolerance && *iterations < max_iterations) {
        const double middle = low + 0.5 * (high - low);
        const ConvexPolygonPose pose_a = pose_at_time(motion_a, middle);
        const ConvexPolygonPose pose_b = pose_at_time(motion_b, middle);
        double separation;
        ++(*iterations);
        if (!convex_polygon_sat_separation(
                a, &pose_a, b, &pose_b, &separation)) {
            return 0;
        }
        if (separation <= distance_tolerance) {
            high = middle;
        } else {
            low = middle;
        }
    }
    {
        const ConvexPolygonPose pose_a = pose_at_time(motion_a, high);
        const ConvexPolygonPose pose_b = pose_at_time(motion_b, high);
        if (!convex_polygon_sat_contact(
                a, &pose_a, b, &pose_b, hit_contact)) return 0;
    }
    *hit_time = high;
    return 1;
}

/* Advance by no more than a certified conservative step.  In particular, do
   not manufacture a minimum time step when floating-point addition cannot
   represent the requested progress: doing so can jump completely over a thin,
   high-speed contact.  The caller reports ITERATION_LIMIT in that case so it
   can rebase/recheck instead of incorrectly claiming NO_HIT. */
static int conservative_time_advance(
    double current,
    double limit,
    double certified_step,
    double *next)
{
    double accepted_step;
    double candidate;

    if (next == NULL
        || !isfinite(current)
        || !isfinite(limit)
        || !isfinite(certified_step)
        || !(limit > current)
        || !(certified_step > 0.0)) {
        return 0;
    }
    accepted_step = fmin(certified_step, limit - current);
    candidate = current + accepted_step;
    if (!(candidate > current) || candidate > limit) {
        return 0;
    }
    *next = candidate;
    return 1;
}

ConvexPolygonCcdStatus convex_polygon_time_of_impact(
    const ConvexPolygon *a,
    const ConvexPolygonMotion *motion_a,
    const ConvexPolygon *b,
    const ConvexPolygonMotion *motion_b,
    double horizon,
    const ConvexPolygonCcdSettings *settings,
    ConvexPolygonCcdResult *result)
{
    ConvexPolygonCcdResult local_result;
    const double radius_a = a != NULL && motion_a != NULL
        ? a->radius * motion_a->pose.scale : 0.0;
    const double radius_b = b != NULL && motion_b != NULL
        ? b->radius * motion_b->pose.scale : 0.0;
    const double final_radius_a = a != NULL && motion_a != NULL
        ? a->radius * (motion_a->pose.scale +
            motion_a->scale_velocity*horizon) : 0.0;
    const double final_radius_b = b != NULL && motion_b != NULL
        ? b->radius * (motion_b->pose.scale +
            motion_b->scale_velocity*horizon) : 0.0;
    const ConvexPolygonVec2 relative_velocity = motion_a != NULL && motion_b != NULL
        ? sub2(motion_b->velocity, motion_a->velocity) : vec2(0.0, 0.0);
    double distance_tolerance;
    double time_tolerance;
    unsigned int max_iterations;
    double interval_start;
    double interval_end;
    double speed_bound;
    double velocity_epsilon;
    double time;
    double previous_time;
    double previous_separation = DBL_MAX;
    int initial_contact_suppressed = 0;
    unsigned int iterations = 0U;

    if (result == NULL) {
        result = &local_result;
    }
    store_ccd_result(result, CONVEX_POLYGON_CCD_INVALID, 0.0, 0U, NULL);
    if (!basic_shape_ok(a)
        || !basic_shape_ok(b)
        || !motion_ok(motion_a)
        || !motion_ok(motion_b)
        || !isfinite(radius_a)
        || !isfinite(radius_b)
        || !isfinite(final_radius_a)
        || !isfinite(final_radius_b)
        || radius_a <= 0.0
        || radius_b <= 0.0
        || final_radius_a <= 0.0
        || final_radius_b <= 0.0
        || !isfinite(horizon)
        || horizon < 0.0) {
        return CONVEX_POLYGON_CCD_INVALID;
    }

    const double maximum_radius_a = fmax(radius_a, final_radius_a);
    const double maximum_radius_b = fmax(radius_b, final_radius_b);
    resolved_ccd_settings(
        settings, maximum_radius_a + maximum_radius_b, horizon,
        &distance_tolerance, &time_tolerance, &max_iterations);
    speed_bound = length2(relative_velocity)
        + fabs(motion_a->angular_velocity) * maximum_radius_a
        + fabs(motion_b->angular_velocity) * maximum_radius_b
        + a->radius*fabs(motion_a->scale_velocity)
        + b->radius*fabs(motion_b->scale_velocity);
    velocity_epsilon = convex_polygon_ccd_velocity_tolerance(speed_bound);

    if (!circumscribed_disk_interval(
            sub2(motion_b->pose.position, motion_a->pose.position),
            relative_velocity,
            radius_a + radius_b + distance_tolerance,
            final_radius_a + final_radius_b + distance_tolerance,
            horizon,
            &interval_start,
            &interval_end)) {
        store_ccd_result(result, CONVEX_POLYGON_CCD_NO_HIT,
                         horizon, iterations, NULL);
        return CONVEX_POLYGON_CCD_NO_HIT;
    }

    time = interval_start;
    previous_time = time;
    while (iterations < max_iterations) {
        const ConvexPolygonPose pose_a = pose_at_time(motion_a, time);
        const ConvexPolygonPose pose_b = pose_at_time(motion_b, time);
        ConvexPolygonContact contact;
        double separation;
        double next_time;

        ++iterations;
        if (!convex_polygon_sat_separation(
                a, &pose_a, b, &pose_b, &separation)) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_INVALID,
                             time, iterations, NULL);
            return CONVEX_POLYGON_CCD_INVALID;
        }

        if (separation <= distance_tolerance &&
            !convex_polygon_sat_contact(
                a, &pose_a, b, &pose_b, &contact)) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_INVALID,
                             time, iterations, NULL);
            return CONVEX_POLYGON_CCD_INVALID;
        }

        if (time == 0.0 && separation < -distance_tolerance) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_INITIAL_OVERLAP,
                             0.0, iterations, &contact);
            return CONVEX_POLYGON_CCD_INITIAL_OVERLAP;
        }

        if (separation <= distance_tolerance) {
            if (time == 0.0 && !initial_contact_suppressed) {
                const double normal_velocity = motion_contact_normal_velocity(
                    motion_a, &pose_a, motion_b, &pose_b, &contact);
                if (normal_velocity < -velocity_epsilon) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_HIT,
                                     0.0, iterations, &contact);
                    return CONVEX_POLYGON_CCD_HIT;
                }
                initial_contact_suppressed = 1;
            } else if (!initial_contact_suppressed) {
                double hit_time = time;
                ConvexPolygonContact hit_contact = contact;
                if (previous_separation > distance_tolerance
                    && previous_time < time
                    && !bisect_pair_contact(
                        a, motion_a, b, motion_b,
                        previous_time, time,
                        distance_tolerance, time_tolerance,
                        max_iterations, &iterations,
                        &hit_time, &hit_contact)) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_INVALID,
                                     time, iterations, NULL);
                    return CONVEX_POLYGON_CCD_INVALID;
                }
                store_ccd_result(result, CONVEX_POLYGON_CCD_HIT,
                                 hit_time, iterations, &hit_contact);
                return CONVEX_POLYGON_CCD_HIT;
            } else {
                const double normal_velocity = motion_contact_normal_velocity(
                    motion_a, &pose_a, motion_b, &pose_b, &contact);
                if (normal_velocity < -velocity_epsilon && time > 0.0) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_HIT,
                                     time, iterations, &contact);
                    return CONVEX_POLYGON_CCD_HIT;
                }
                if (speed_bound <= velocity_epsilon) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_NO_HIT,
                                     horizon, iterations, NULL);
                    return CONVEX_POLYGON_CCD_NO_HIT;
                }
                if (time >= interval_end) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_NO_HIT,
                                     horizon, iterations, NULL);
                    return CONVEX_POLYGON_CCD_NO_HIT;
                }
                /* Stay inside the numerical contact band until the initially
                   separating/tangent contact is demonstrably clear.  The
                   speed bound limits total relative boundary motion, unlike a
                   normal-velocity extrapolation whose direction can change as
                   a rotating support feature changes. */
                if (!conservative_time_advance(
                        time, interval_end,
                        0.5 * distance_tolerance / speed_bound,
                        &next_time)) {
                    store_ccd_result(result,
                                     CONVEX_POLYGON_CCD_ITERATION_LIMIT,
                                     time, iterations, &contact);
                    return CONVEX_POLYGON_CCD_ITERATION_LIMIT;
                }
                previous_time = time;
                previous_separation = separation;
                time = next_time;
                continue;
            }
        } else if (initial_contact_suppressed) {
            initial_contact_suppressed = 0;
        }

        if (time >= interval_end) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_NO_HIT,
                             horizon, iterations, NULL);
            return CONVEX_POLYGON_CCD_NO_HIT;
        }
        if (speed_bound <= velocity_epsilon) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_NO_HIT,
                             horizon, iterations, NULL);
            return CONVEX_POLYGON_CCD_NO_HIT;
        }

        previous_time = time;
        previous_separation = separation;
        if (!conservative_time_advance(
                time, interval_end,
                0.9 * separation / speed_bound,
                &next_time)) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_ITERATION_LIMIT,
                             time, iterations, NULL);
            return CONVEX_POLYGON_CCD_ITERATION_LIMIT;
        }
        time = next_time;
    }

    store_ccd_result(result, CONVEX_POLYGON_CCD_ITERATION_LIMIT,
                     time, iterations, NULL);
    return CONVEX_POLYGON_CCD_ITERATION_LIMIT;
}

int convex_polygon_plane_contact(
    const ConvexPolygon *polygon,
    const ConvexPolygonPose *pose,
    const ConvexPolygonPlane *plane,
    ConvexPolygonContact *contact)
{
    ConvexPolygonVec2 vertices[CONVEX_POLYGON_MAX_VERTICES];
    ConvexPolygonVec2 normal;
    double normal_length;
    double normalized_offset;
    double maximum;
    int i;

    if (plane == NULL
        || contact == NULL
        || !finite_vec2(plane->normal)
        || !isfinite(plane->offset)
        || !transform_vertices(polygon, pose, vertices, NULL, NULL)) {
        return 0;
    }
    normal_length = length2(plane->normal);
    if (!(normal_length > DBL_MIN) || !isfinite(normal_length)) {
        return 0;
    }
    normal = mul2(plane->normal, 1.0 / normal_length);
    normalized_offset = plane->offset / normal_length;
    maximum = dot2(vertices[0], normal);
    for (i = 1; i < polygon->count; ++i) {
        const double projection = dot2(vertices[i], normal);
        if (projection > maximum) maximum = projection;
    }

    contact->separation = normalized_offset - maximum;
    contact->normal = normal;
    contact->point_a = averaged_support(
        vertices, polygon->count, normal, 1);
    contact->point_b = add2(
        contact->point_a, mul2(normal, contact->separation));
    contact->point = mul2(add2(contact->point_a, contact->point_b), 0.5);
    contact->axis_from_a = 1;
    contact->edge_index = -1;
    return isfinite(contact->separation) && finite_vec2(contact->point);
}

static double plane_contact_normal_velocity(
    const ConvexPolygonMotion *motion,
    const ConvexPolygonPose *pose,
    const ConvexPolygonContact *contact)
{
    const ConvexPolygonVec2 lever = sub2(contact->point, pose->position);
    ConvexPolygonVec2 point_velocity = add2(
        motion->velocity,
        angular_velocity_at_point(motion->angular_velocity, lever));
    if (pose->scale > DBL_MIN && motion->scale_velocity != 0.0)
        point_velocity = add2(point_velocity, mul2(
            lever, motion->scale_velocity/pose->scale));
    /* Relative wall-minus-body velocity along body-to-wall normal. */
    return -dot2(point_velocity, contact->normal);
}

ConvexPolygonCcdStatus convex_polygon_plane_time_of_impact(
    const ConvexPolygon *polygon,
    const ConvexPolygonMotion *motion,
    const ConvexPolygonPlane *plane,
    double horizon,
    const ConvexPolygonCcdSettings *settings,
    ConvexPolygonCcdResult *result)
{
    ConvexPolygonCcdResult local_result;
    const double radius = polygon != NULL && motion != NULL
        ? polygon->radius * motion->pose.scale : 0.0;
    const double final_radius = polygon != NULL && motion != NULL
        ? polygon->radius * (motion->pose.scale +
            motion->scale_velocity*horizon) : 0.0;
    double distance_tolerance;
    double time_tolerance;
    unsigned int max_iterations;
    double speed_bound;
    double velocity_epsilon;
    double time = 0.0;
    double previous_time = 0.0;
    double previous_separation = DBL_MAX;
    int initial_contact_suppressed = 0;
    unsigned int iterations = 0U;

    if (result == NULL) {
        result = &local_result;
    }
    store_ccd_result(result, CONVEX_POLYGON_CCD_INVALID, 0.0, 0U, NULL);
    if (!basic_shape_ok(polygon)
        || !motion_ok(motion)
        || plane == NULL
        || !finite_vec2(plane->normal)
        || !isfinite(plane->offset)
        || !(length2(plane->normal) > DBL_MIN)
        || !isfinite(radius)
        || !isfinite(final_radius)
        || radius <= 0.0
        || final_radius <= 0.0
        || !isfinite(horizon)
        || horizon < 0.0) {
        return CONVEX_POLYGON_CCD_INVALID;
    }

    const double maximum_radius = fmax(radius, final_radius);
    resolved_ccd_settings(
        settings, maximum_radius, horizon,
        &distance_tolerance, &time_tolerance, &max_iterations);
    {
        const ConvexPolygonVec2 normal = mul2(
            plane->normal, 1.0 / length2(plane->normal));
        speed_bound = fabs(dot2(motion->velocity, normal))
            + fabs(motion->angular_velocity) * maximum_radius
            + polygon->radius*fabs(motion->scale_velocity);
    }
    velocity_epsilon = convex_polygon_ccd_velocity_tolerance(speed_bound);

    while (iterations < max_iterations) {
        const ConvexPolygonPose pose = pose_at_time(motion, time);
        ConvexPolygonContact contact;
        double next_time;

        ++iterations;
        if (!convex_polygon_plane_contact(polygon, &pose, plane, &contact)) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_INVALID,
                             time, iterations, NULL);
            return CONVEX_POLYGON_CCD_INVALID;
        }
        if (time == 0.0 && contact.separation < -distance_tolerance) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_INITIAL_OVERLAP,
                             0.0, iterations, &contact);
            return CONVEX_POLYGON_CCD_INITIAL_OVERLAP;
        }

        if (contact.separation <= distance_tolerance) {
            if (time == 0.0 && !initial_contact_suppressed) {
                const double normal_velocity = plane_contact_normal_velocity(
                    motion, &pose, &contact);
                if (normal_velocity < -velocity_epsilon) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_HIT,
                                     0.0, iterations, &contact);
                    return CONVEX_POLYGON_CCD_HIT;
                }
                initial_contact_suppressed = 1;
            } else if (!initial_contact_suppressed) {
                double low = previous_time;
                double high = time;
                ConvexPolygonContact hit_contact = contact;
                while (previous_separation > distance_tolerance
                    && high - low > time_tolerance
                    && iterations < max_iterations) {
                    const double middle = low + 0.5 * (high - low);
                    const ConvexPolygonPose middle_pose = pose_at_time(motion, middle);
                    ConvexPolygonContact middle_contact;
                    ++iterations;
                    if (!convex_polygon_plane_contact(
                            polygon, &middle_pose, plane, &middle_contact)) {
                        store_ccd_result(result, CONVEX_POLYGON_CCD_INVALID,
                                         time, iterations, NULL);
                        return CONVEX_POLYGON_CCD_INVALID;
                    }
                    if (middle_contact.separation <= distance_tolerance) {
                        high = middle;
                        hit_contact = middle_contact;
                    } else {
                        low = middle;
                    }
                }
                store_ccd_result(result, CONVEX_POLYGON_CCD_HIT,
                                 high, iterations, &hit_contact);
                return CONVEX_POLYGON_CCD_HIT;
            } else {
                const double normal_velocity = plane_contact_normal_velocity(
                    motion, &pose, &contact);
                if (normal_velocity < -velocity_epsilon && time > 0.0) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_HIT,
                                     time, iterations, &contact);
                    return CONVEX_POLYGON_CCD_HIT;
                }
                if (speed_bound <= velocity_epsilon) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_NO_HIT,
                                     horizon, iterations, NULL);
                    return CONVEX_POLYGON_CCD_NO_HIT;
                }
                if (time >= horizon) {
                    store_ccd_result(result, CONVEX_POLYGON_CCD_NO_HIT,
                                     horizon, iterations, NULL);
                    return CONVEX_POLYGON_CCD_NO_HIT;
                }
                if (!conservative_time_advance(
                        time, horizon,
                        0.5 * distance_tolerance / speed_bound,
                        &next_time)) {
                    store_ccd_result(result,
                                     CONVEX_POLYGON_CCD_ITERATION_LIMIT,
                                     time, iterations, &contact);
                    return CONVEX_POLYGON_CCD_ITERATION_LIMIT;
                }
                previous_time = time;
                previous_separation = contact.separation;
                time = next_time;
                continue;
            }
        } else if (initial_contact_suppressed) {
            initial_contact_suppressed = 0;
        }

        if (time >= horizon || speed_bound <= velocity_epsilon) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_NO_HIT,
                             horizon, iterations, NULL);
            return CONVEX_POLYGON_CCD_NO_HIT;
        }
        previous_time = time;
        previous_separation = contact.separation;
        if (!conservative_time_advance(
                time, horizon,
                0.9 * contact.separation / speed_bound,
                &next_time)) {
            store_ccd_result(result, CONVEX_POLYGON_CCD_ITERATION_LIMIT,
                             time, iterations, &contact);
            return CONVEX_POLYGON_CCD_ITERATION_LIMIT;
        }
        time = next_time;
    }

    store_ccd_result(result, CONVEX_POLYGON_CCD_ITERATION_LIMIT,
                     time, iterations, NULL);
    return CONVEX_POLYGON_CCD_ITERATION_LIMIT;
}

double convex_polygon_contact_normal_velocity(
    const ConvexPolygonRigidBody *a,
    const ConvexPolygonRigidBody *b,
    const ConvexPolygonContact *contact)
{
    ConvexPolygonVec2 normal;
    double normal_length;
    ConvexPolygonVec2 velocity_a;
    ConvexPolygonVec2 velocity_b;

    if (a == NULL || b == NULL || contact == NULL
        || !finite_vec2(a->position) || !finite_vec2(b->position)
        || !finite_vec2(a->velocity) || !finite_vec2(b->velocity)
        || !isfinite(a->angular_velocity) || !isfinite(b->angular_velocity)
        || !finite_vec2(contact->point) || !finite_vec2(contact->normal)) {
        return NAN;
    }
    normal_length = length2(contact->normal);
    if (!(normal_length > DBL_MIN) || !isfinite(normal_length)) {
        return NAN;
    }
    normal = mul2(contact->normal, 1.0 / normal_length);
    velocity_a = add2(
        a->velocity,
        angular_velocity_at_point(
            a->angular_velocity, sub2(contact->point, a->position)));
    velocity_b = add2(
        b->velocity,
        angular_velocity_at_point(
            b->angular_velocity, sub2(contact->point, b->position)));
    return dot2(sub2(velocity_b, velocity_a), normal);
}

static double inverse_if_dynamic(double value)
{
    return isfinite(value) && value > 0.0 ? 1.0 / value : 0.0;
}

double convex_polygon_apply_pair_impulse(
    ConvexPolygonRigidBody *a,
    ConvexPolygonRigidBody *b,
    const ConvexPolygonContact *contact,
    double restitution)
{
    ConvexPolygonVec2 normal;
    ConvexPolygonVec2 lever_a;
    ConvexPolygonVec2 lever_b;
    double normal_length;
    double inverse_mass_a;
    double inverse_mass_b;
    double inverse_inertia_a;
    double inverse_inertia_b;
    double lever_cross_a;
    double lever_cross_b;
    double normal_velocity;
    double effective_inverse_mass;
    double impulse;

    if (a == NULL || b == NULL || contact == NULL
        || !isfinite(restitution) || restitution < 0.0
        || !finite_vec2(contact->normal) || !finite_vec2(contact->point)) {
        return 0.0;
    }
    normal_length = length2(contact->normal);
    if (!(normal_length > DBL_MIN) || !isfinite(normal_length)) {
        return 0.0;
    }
    normal = mul2(contact->normal, 1.0 / normal_length);
    inverse_mass_a = inverse_if_dynamic(a->mass);
    inverse_mass_b = inverse_if_dynamic(b->mass);
    inverse_inertia_a = inverse_if_dynamic(a->inertia);
    inverse_inertia_b = inverse_if_dynamic(b->inertia);
    lever_a = sub2(contact->point, a->position);
    lever_b = sub2(contact->point, b->position);
    lever_cross_a = cross2(lever_a, normal);
    lever_cross_b = cross2(lever_b, normal);
    normal_velocity = convex_polygon_contact_normal_velocity(a, b, contact);
    effective_inverse_mass = inverse_mass_a + inverse_mass_b
        + lever_cross_a * lever_cross_a * inverse_inertia_a
        + lever_cross_b * lever_cross_b * inverse_inertia_b;

    if (!isfinite(normal_velocity)
        || normal_velocity >= 0.0
        || !(effective_inverse_mass > 0.0)
        || !isfinite(effective_inverse_mass)) {
        return 0.0;
    }
    impulse = -(1.0 + restitution) * normal_velocity
        / effective_inverse_mass;
    a->velocity = sub2(a->velocity, mul2(normal, impulse * inverse_mass_a));
    b->velocity = add2(b->velocity, mul2(normal, impulse * inverse_mass_b));
    a->angular_velocity -= impulse * lever_cross_a * inverse_inertia_a;
    b->angular_velocity += impulse * lever_cross_b * inverse_inertia_b;
    return impulse;
}

double convex_polygon_apply_fixed_wall_impulse(
    ConvexPolygonRigidBody *body,
    ConvexPolygonVec2 contact_point,
    ConvexPolygonVec2 wall_normal,
    double restitution)
{
    ConvexPolygonVec2 normal;
    ConvexPolygonVec2 lever;
    ConvexPolygonVec2 point_velocity;
    double normal_length;
    double inverse_mass;
    double inverse_inertia;
    double lever_cross;
    double normal_velocity;
    double effective_inverse_mass;
    double impulse;

    if (body == NULL
        || !finite_vec2(contact_point)
        || !finite_vec2(wall_normal)
        || !isfinite(restitution)
        || restitution < 0.0) {
        return 0.0;
    }
    normal_length = length2(wall_normal);
    if (!(normal_length > DBL_MIN) || !isfinite(normal_length)) {
        return 0.0;
    }
    normal = mul2(wall_normal, 1.0 / normal_length);
    inverse_mass = inverse_if_dynamic(body->mass);
    inverse_inertia = inverse_if_dynamic(body->inertia);
    lever = sub2(contact_point, body->position);
    lever_cross = cross2(lever, normal);
    point_velocity = add2(
        body->velocity,
        angular_velocity_at_point(body->angular_velocity, lever));
    normal_velocity = -dot2(point_velocity, normal);
    effective_inverse_mass = inverse_mass
        + lever_cross * lever_cross * inverse_inertia;

    if (!isfinite(normal_velocity)
        || normal_velocity >= 0.0
        || !(effective_inverse_mass > 0.0)
        || !isfinite(effective_inverse_mass)) {
        return 0.0;
    }
    impulse = -(1.0 + restitution) * normal_velocity
        / effective_inverse_mass;
    body->velocity = sub2(
        body->velocity, mul2(normal, impulse * inverse_mass));
    body->angular_velocity -= impulse * lever_cross * inverse_inertia;
    return impulse;
}
