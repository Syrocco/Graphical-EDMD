#include "convex_polygon.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440084436210484904
#endif

static int failures = 0;

#define CHECK(condition)                                                       \
    do {                                                                       \
        if (!(condition)) {                                                    \
            fprintf(stderr, "%s:%d: check failed: %s\n",                    \
                    __FILE__, __LINE__, #condition);                           \
            ++failures;                                                        \
        }                                                                      \
    } while (0)

#define CHECK_NEAR(actual, expected, tolerance)                                \
    do {                                                                       \
        const double check_actual_ = (actual);                                 \
        const double check_expected_ = (expected);                             \
        const double check_tolerance_ = (tolerance);                           \
        if (!(fabs(check_actual_ - check_expected_) <= check_tolerance_)) {    \
            fprintf(stderr, "%s:%d: %.17g != %.17g within %.3g\n",          \
                    __FILE__, __LINE__, check_actual_, check_expected_,        \
                    check_tolerance_);                                         \
            ++failures;                                                        \
        }                                                                      \
    } while (0)

static ConvexPolygonVec2 v2(double x, double y)
{
    ConvexPolygonVec2 value;
    value.x = x;
    value.y = y;
    return value;
}

static ConvexPolygon axis_aligned_square(void)
{
    ConvexPolygon polygon;
    memset(&polygon, 0, sizeof(polygon));
    polygon.count = 4;
    polygon.vertices[0] = v2(-1.0, -1.0);
    polygon.vertices[1] = v2( 1.0, -1.0);
    polygon.vertices[2] = v2( 1.0,  1.0);
    polygon.vertices[3] = v2(-1.0,  1.0);
    CHECK(convex_polygon_recenter_and_normalize(&polygon));
    return polygon;
}

static ConvexPolygon thin_rectangle(void)
{
    ConvexPolygon polygon;
    memset(&polygon, 0, sizeof(polygon));
    polygon.count = 4;
    polygon.vertices[0] = v2(-1.0, -0.2);
    polygon.vertices[1] = v2( 1.0, -0.2);
    polygon.vertices[2] = v2( 1.0,  0.2);
    polygon.vertices[3] = v2(-1.0,  0.2);
    CHECK(convex_polygon_recenter_and_normalize(&polygon));
    return polygon;
}

static ConvexPolygonPose pose(double x, double y, double angle, double scale)
{
    ConvexPolygonPose result;
    result.position = v2(x, y);
    result.angle = angle;
    result.scale = scale;
    return result;
}

static ConvexPolygonMotion motion(
    double x, double y, double angle, double scale,
    double vx, double vy, double angular_velocity)
{
    ConvexPolygonMotion result;
    result.pose = pose(x, y, angle, scale);
    result.velocity = v2(vx, vy);
    result.angular_velocity = angular_velocity;
    result.scale_velocity = 0.0;
    return result;
}

static double deterministic_fraction(unsigned int *state)
{
    *state = 1664525U*(*state) + 1013904223U;
    return (double)(*state >> 8)*(1.0/16777216.0);
}

static double brute_sat_separation(
    const ConvexPolygon *a,
    const ConvexPolygonPose *pose_a,
    const ConvexPolygon *b,
    const ConvexPolygonPose *pose_b)
{
    ConvexPolygonVec2 vertices_a[CONVEX_POLYGON_MAX_VERTICES];
    ConvexPolygonVec2 vertices_b[CONVEX_POLYGON_MAX_VERTICES];
    double best = -INFINITY;
    int source;
    int i;

    for (i = 0; i < a->count; ++i) {
        const double cosine = cos(pose_a->angle);
        const double sine = sin(pose_a->angle);
        const ConvexPolygonVec2 local = a->vertices[i];
        vertices_a[i] = v2(
            pose_a->position.x + pose_a->scale*
                (cosine*local.x-sine*local.y),
            pose_a->position.y + pose_a->scale*
                (sine*local.x+cosine*local.y));
    }
    for (i = 0; i < b->count; ++i) {
        const double cosine = cos(pose_b->angle);
        const double sine = sin(pose_b->angle);
        const ConvexPolygonVec2 local = b->vertices[i];
        vertices_b[i] = v2(
            pose_b->position.x + pose_b->scale*
                (cosine*local.x-sine*local.y),
            pose_b->position.y + pose_b->scale*
                (sine*local.x+cosine*local.y));
    }

    for (source = 0; source < 2; ++source) {
        const ConvexPolygonVec2 *vertices = source == 0
            ? vertices_a : vertices_b;
        const int count = source == 0 ? a->count : b->count;
        for (i = 0; i < count; ++i) {
            const int next = i + 1 == count ? 0 : i + 1;
            const double edge_x = vertices[next].x-vertices[i].x;
            const double edge_y = vertices[next].y-vertices[i].y;
            const double inverse_length = 1.0/hypot(edge_x, edge_y);
            const ConvexPolygonVec2 axis = v2(
                edge_y*inverse_length, -edge_x*inverse_length);
            double minimum_a = vertices_a[0].x*axis.x
                + vertices_a[0].y*axis.y;
            double maximum_a = minimum_a;
            double minimum_b = vertices_b[0].x*axis.x
                + vertices_b[0].y*axis.y;
            double maximum_b = minimum_b;
            int vertex;
            for (vertex = 1; vertex < a->count; ++vertex) {
                const double projection = vertices_a[vertex].x*axis.x
                    + vertices_a[vertex].y*axis.y;
                if (projection < minimum_a) minimum_a = projection;
                if (projection > maximum_a) maximum_a = projection;
            }
            for (vertex = 1; vertex < b->count; ++vertex) {
                const double projection = vertices_b[vertex].x*axis.x
                    + vertices_b[vertex].y*axis.y;
                if (projection < minimum_b) minimum_b = projection;
                if (projection > maximum_b) maximum_b = projection;
            }
            {
                const double candidate = fmax(minimum_b-maximum_a,
                                              minimum_a-maximum_b);
                if (candidate > best) best = candidate;
            }
        }
    }
    return best;
}

static double body_energy(const ConvexPolygonRigidBody *body)
{
    return 0.5 * body->mass
        * (body->velocity.x * body->velocity.x
           + body->velocity.y * body->velocity.y)
        + 0.5 * body->inertia
        * body->angular_velocity * body->angular_velocity;
}

static void test_regular_polygons_and_properties(void)
{
    int count;
    for (count = 3; count <= CONVEX_POLYGON_MAX_VERTICES; ++count) {
        ConvexPolygon polygon;
        CHECK(convex_polygon_make_regular(&polygon, count));
        CHECK(polygon.count == count);
        CHECK(convex_polygon_validate(&polygon) == CONVEX_POLYGON_VALID);
        CHECK_NEAR(polygon.centroid.x, 0.0, 2.0e-15);
        CHECK_NEAR(polygon.centroid.y, 0.0, 2.0e-15);
        CHECK_NEAR(polygon.radius, 1.0, 2.0e-15);
        CHECK(polygon.area > 0.0);
        CHECK(polygon.inertia_per_mass > 0.0);
    }

    {
        ConvexPolygon square = axis_aligned_square();
        CHECK_NEAR(square.inertia_per_mass, 1.0 / 3.0, 2.0e-14);
        CHECK_NEAR(square.area, 2.0, 2.0e-14);
    }

    {
        ConvexPolygon translated;
        memset(&translated, 0, sizeof(translated));
        translated.count = 4;
        translated.vertices[0] = v2(7.0, -5.0);
        translated.vertices[1] = v2(11.0, -5.0);
        translated.vertices[2] = v2(11.0, -3.0);
        translated.vertices[3] = v2(7.0, -3.0);
        CHECK(convex_polygon_recenter_and_normalize(&translated));
        CHECK_NEAR(translated.centroid.x, 0.0, 1.0e-14);
        CHECK_NEAR(translated.centroid.y, 0.0, 1.0e-14);
        CHECK_NEAR(translated.radius, 1.0, 1.0e-14);
    }
}

static void test_validation_failures(void)
{
    ConvexPolygon polygon = axis_aligned_square();
    ConvexPolygon bad;

    bad = polygon;
    bad.count = 2;
    CHECK(convex_polygon_validate(&bad)
          == CONVEX_POLYGON_INVALID_VERTEX_COUNT);

    bad = polygon;
    bad.vertices[1].x = NAN;
    CHECK(convex_polygon_validate(&bad)
          == CONVEX_POLYGON_NONFINITE_VERTEX);

    bad = polygon;
    {
        const ConvexPolygonVec2 swap = bad.vertices[1];
        bad.vertices[1] = bad.vertices[3];
        bad.vertices[3] = swap;
    }
    CHECK(convex_polygon_validate(&bad)
          == CONVEX_POLYGON_NOT_STRICTLY_CCW_CONVEX);

    bad = polygon;
    bad.vertices[2] = v2(0.0, 0.0);
    CHECK(convex_polygon_validate(&bad)
          == CONVEX_POLYGON_NOT_STRICTLY_CCW_CONVEX);

    bad = polygon;
    bad.vertices[2] = bad.vertices[1];
    CHECK(convex_polygon_validate(&bad)
          == CONVEX_POLYGON_EDGE_TOO_SHORT);

    CHECK(strcmp(convex_polygon_validation_message(CONVEX_POLYGON_VALID),
                 "valid") == 0);
}

static void test_clamped_interactive_edits(void)
{
    ConvexPolygon polygon = axis_aligned_square();
    double accepted;

    accepted = convex_polygon_move_vertex_clamped(
        &polygon, 2, v2(-1.5, -1.5));
    CHECK(accepted >= 0.0 && accepted < 1.0);
    CHECK(convex_polygon_validate(&polygon) == CONVEX_POLYGON_VALID);

    accepted = convex_polygon_move_vertex_clamped(
        &polygon, 2, v2(polygon.vertices[2].x + 0.1,
                        polygon.vertices[2].y + 0.1));
    CHECK_NEAR(accepted, 1.0, 0.0);
    CHECK(convex_polygon_validate(&polygon) == CONVEX_POLYGON_VALID);

    polygon = axis_aligned_square();
    accepted = convex_polygon_translate_edge_clamped(
        &polygon, 2, v2(0.0, -3.0));
    CHECK(accepted >= 0.0 && accepted < 1.0);
    CHECK(convex_polygon_validate(&polygon) == CONVEX_POLYGON_VALID);

    CHECK(convex_polygon_recenter_and_normalize(&polygon));
    CHECK_NEAR(polygon.radius, 1.0, 2.0e-14);
}

static void test_sat_separation_and_witnesses(void)
{
    const ConvexPolygon square = axis_aligned_square();
    const double half_width = M_SQRT1_2;
    ConvexPolygonPose a = pose(0.0, 0.0, 0.0, 1.0);
    ConvexPolygonPose b = pose(2.0 * half_width + 0.25, 0.0, 0.0, 1.0);
    ConvexPolygonContact contact;

    CHECK(convex_polygon_sat_contact(&square, &a, &square, &b, &contact));
    CHECK_NEAR(contact.separation, 0.25, 2.0e-14);
    CHECK_NEAR(contact.normal.x, 1.0, 2.0e-14);
    CHECK_NEAR(contact.normal.y, 0.0, 2.0e-14);
    CHECK_NEAR(contact.point_b.x - contact.point_a.x, 0.25, 2.0e-14);

    b.position.x = 2.0 * half_width;
    CHECK(convex_polygon_sat_contact(&square, &a, &square, &b, &contact));
    CHECK_NEAR(contact.separation, 0.0, 2.0e-14);

    b.position.x -= 0.20;
    CHECK(convex_polygon_sat_contact(&square, &a, &square, &b, &contact));
    CHECK_NEAR(contact.separation, -0.20, 2.0e-14);

    b = pose(0.0, 2.2, M_PI / 4.0, 1.0);
    CHECK(convex_polygon_sat_contact(&square, &a, &square, &b, &contact));
    CHECK(contact.separation > 0.0);
    CHECK(contact.normal.y > 0.9);
}

static void test_cached_sat_matches_brute_projection(void)
{
    unsigned int state = 0x51a7c0deU;
    int sample;

    for (sample = 0; sample < 1024; ++sample) {
        ConvexPolygon a;
        ConvexPolygon b;
        ConvexPolygonContact contact;
        ConvexPolygonPose pose_a;
        ConvexPolygonPose pose_b;
        double expected;
        double tolerance;
        int vertex;

        CHECK(convex_polygon_make_regular(&a, 3 + sample%14));
        CHECK(convex_polygon_make_regular(&b, 3 + (sample*5)%14));
        for (vertex = 0; vertex < a.count; ++vertex) {
            const ConvexPolygonVec2 original = a.vertices[vertex];
            a.vertices[vertex] = v2(original.x + 0.23*original.y,
                                    0.61*original.y);
        }
        for (vertex = 0; vertex < b.count; ++vertex) {
            const ConvexPolygonVec2 original = b.vertices[vertex];
            b.vertices[vertex] = v2(0.77*original.x - 0.17*original.y,
                                    original.y);
        }
        CHECK(convex_polygon_recenter_and_normalize(&a));
        CHECK(convex_polygon_recenter_and_normalize(&b));

        pose_a = pose(
            40.0*deterministic_fraction(&state)-20.0,
            40.0*deterministic_fraction(&state)-20.0,
            2.0*M_PI*deterministic_fraction(&state)-M_PI,
            0.2 + 1.8*deterministic_fraction(&state));
        pose_b = pose(
            pose_a.position.x + 6.0*deterministic_fraction(&state)-3.0,
            pose_a.position.y + 6.0*deterministic_fraction(&state)-3.0,
            2.0*M_PI*deterministic_fraction(&state)-M_PI,
            0.2 + 1.8*deterministic_fraction(&state));
        expected = brute_sat_separation(&a, &pose_a, &b, &pose_b);
        CHECK(convex_polygon_sat_contact(
            &a, &pose_a, &b, &pose_b, &contact));
        tolerance = 2.0e-13*(1.0 + fabs(expected));
        CHECK_NEAR(contact.separation, expected, tolerance);
        if (fabs(expected) > 4.0*tolerance) {
            CHECK((contact.separation > 0.0) == (expected > 0.0));
        }
    }
}

static ConvexPolygonCcdSettings tight_ccd_settings(void)
{
    ConvexPolygonCcdSettings settings;
    settings.distance_tolerance = 1.0e-12;
    settings.time_tolerance = 1.0e-12;
    settings.max_iterations = 512U;
    return settings;
}

static void test_translational_ccd(void)
{
    const ConvexPolygon square = axis_aligned_square();
    const double half_width = M_SQRT1_2;
    const double expected = (6.0 - 2.0 * half_width) / 2.0;
    const ConvexPolygonCcdSettings settings = tight_ccd_settings();
    ConvexPolygonMotion a = motion(-3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    ConvexPolygonMotion b = motion( 3.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0);
    ConvexPolygonCcdResult result;

    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 4.0, &settings, &result)
        == CONVEX_POLYGON_CCD_HIT);
    CHECK_NEAR(result.time, expected, 2.0e-10);
    CHECK(result.contact.normal.x > 0.99);

    a.velocity.x = -1.0;
    b.velocity.x = 1.0;
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 10.0, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);

    a = motion(-half_width, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0);
    b = motion( half_width, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 1.0, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);

    a.velocity.x = 1.0;
    b.velocity.x = -1.0;
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 1.0, &settings, &result)
        == CONVEX_POLYGON_CCD_HIT);
    CHECK_NEAR(result.time, 0.0, 0.0);

    a = motion(-0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    b = motion( 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 1.0, &settings, &result)
        == CONVEX_POLYGON_CCD_INITIAL_OVERLAP);

    a = motion(-100.0, 0.0, 0.0, 1.0, 1000.0, 0.0, 0.0);
    b = motion( 100.0, 0.0, 0.0, 1.0, -1000.0, 0.0, 0.0);
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 0.2, &settings, &result)
        == CONVEX_POLYGON_CCD_HIT);
    CHECK_NEAR(result.time,
               (200.0 - 2.0 * half_width) / 2000.0,
               2.0e-10);

    /* A default time tolerance scales with the requested horizon.  It is a
       bisection stopping criterion, never permission to inflate a conservative
       advancement step: with the huge horizon below, such inflation jumps
       across the complete high-speed overlap interval. */
    a = motion(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    b = motion(3.0, 0.0, 0.0, 1.0, -1.0e6, 0.0, 0.0);
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 1.0e16, NULL, &result)
        == CONVEX_POLYGON_CCD_HIT);
    CHECK_NEAR(result.time,
               (3.0 - 2.0 * half_width) / 1.0e6,
               2.0e-15);
}

static void test_rotational_ccd(void)
{
    const ConvexPolygon bar = thin_rectangle();
    ConvexPolygon small_square = axis_aligned_square();
    const ConvexPolygonCcdSettings settings = tight_ccd_settings();
    ConvexPolygonMotion a = motion(
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, M_PI / 2.0);
    ConvexPolygonMotion b = motion(
        0.0, 1.08, 0.0, 0.2, 0.0, 0.0, 0.0);
    ConvexPolygonCcdResult result;

    CHECK(convex_polygon_time_of_impact(
        &bar, &a, &small_square, &b, 1.0, &settings, &result)
        == CONVEX_POLYGON_CCD_HIT);
    CHECK(result.time > 0.2 && result.time < 1.0);

    CHECK(convex_polygon_time_of_impact(
        &bar, &a, &small_square, &b, 0.1, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);

    /* Vertex-on-face at t=0: the true vertex is separating even though the
       midpoint between independently averaged support features would appear
       to approach.  CCD must use the clipped support witness and suppress the
       zero-time event. */
    {
        const ConvexPolygon square = axis_aligned_square();
        a = motion(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0);
        b = motion(1.0 + M_SQRT1_2, 0.5, M_PI/4.0, 1.0,
                   -0.375, 0.0, 0.0);
        CHECK(convex_polygon_time_of_impact(
            &square, &a, &square, &b, 1.0e-4, &settings, &result)
            == CONVEX_POLYGON_CCD_NO_HIT);
    }
}

static void test_homothetic_growth_ccd(void)
{
    const ConvexPolygon square = axis_aligned_square();
    const ConvexPolygonCcdSettings settings = tight_ccd_settings();
    ConvexPolygonMotion a = motion(
        -2.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0);
    ConvexPolygonMotion b = motion(
         2.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0);
    ConvexPolygonCcdResult result;
    const double expected_pair =
        (4.0/(2.0*M_SQRT1_2) - 0.25)/0.5;
    a.scale_velocity = 0.5;
    b.scale_velocity = 0.5;
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 6.0, &settings, &result)
        == CONVEX_POLYGON_CCD_HIT);
    CHECK_NEAR(result.time, expected_pair, 3.0e-10);

    {
        const ConvexPolygonPlane plane = { { 1.0, 0.0 }, 2.0 };
        ConvexPolygonMotion body = motion(
            0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0);
        const double expected_plane =
            (2.0/M_SQRT1_2 - 0.25)/0.5;
        body.scale_velocity = 0.5;
        CHECK(convex_polygon_plane_time_of_impact(
            &square, &body, &plane, 6.0, &settings, &result)
            == CONVEX_POLYGON_CCD_HIT);
        CHECK_NEAR(result.time, expected_plane, 3.0e-10);
    }
}

static void test_scale_aware_disk_broadphase_edges(void)
{
    const ConvexPolygon square = axis_aligned_square();
    const ConvexPolygonCcdSettings settings = tight_ccd_settings();
    ConvexPolygonCcdResult result;
    ConvexPolygonMotion a;
    ConvexPolygonMotion b;

    /* Stationary centers whose final expanding envelopes remain disjoint are
       rejected without running even one SAT iteration. */
    a = motion(-5.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0);
    b = motion( 5.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0);
    a.scale_velocity = b.scale_velocity = 0.01;
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 10.0, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);
    CHECK(result.iterations == 0U);

    /* Equal radial and center-separation rates reduce the squared envelope
       equation to a linear one; a fixed positive gap stays a no-hit. */
    a = motion(0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0);
    b = motion(4.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
    a.scale_velocity = b.scale_velocity = 0.5;
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 10.0, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);
    CHECK(result.iterations == 0U);

    /* Expanding envelopes can catch centers that are moving apart. */
    a = motion(0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0);
    b = motion(4.0, 0.0, 0.0, 0.25, 0.25, 0.0, 0.0);
    a.scale_velocity = b.scale_velocity = 0.5;
    {
        const double expected = (4.0 - 0.5*M_SQRT1_2)/
            (M_SQRT1_2 - 0.25);
        CHECK(convex_polygon_time_of_impact(
            &square, &a, &square, &b, 10.0, &settings, &result)
            == CONVEX_POLYGON_CCD_HIT);
        CHECK_NEAR(result.time, expected, 4.0e-10);
    }

    /* If center separation grows faster than the two envelopes, expansion
       must not manufacture a candidate SAT interval. */
    b.velocity.x = 1.0;
    a.scale_velocity = b.scale_velocity = 0.1;
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 20.0, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);
    CHECK(result.iterations == 0U);

    /* Shrinking bodies remain supported: translation can still close the gap
       before either scale reaches zero. */
    a = motion(-3.0, 0.0, 0.0, 2.0, 0.5, 0.0, 0.0);
    b = motion( 3.0, 0.0, 0.0, 2.0, -0.5, 0.0, 0.0);
    a.scale_velocity = b.scale_velocity = -0.1;
    {
        const double expected = (6.0 - 4.0*M_SQRT1_2)/
            (1.0 - 0.2*M_SQRT1_2);
        CHECK(convex_polygon_time_of_impact(
            &square, &a, &square, &b, 6.0, &settings, &result)
            == CONVEX_POLYGON_CCD_HIT);
        CHECK_NEAR(result.time, expected, 4.0e-10);
    }

    /* A double-root disk tangency is retained conservatively, but the SAT
       narrow phase correctly rejects the non-touching squares. */
    a = motion(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    b = motion(-5.0, 2.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 10.0, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);
    CHECK(result.iterations > 0U);

    /* Normalized-time coefficients avoid overflow/cancellation for a huge
       horizon and a correspondingly tiny expansion rate. */
    a = motion(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    b = motion(10.0, 0.0, 0.0, 1.0, 1.0e-16, 0.0, 0.0);
    a.scale_velocity = b.scale_velocity = 5.0e-17;
    CHECK(convex_polygon_time_of_impact(
        &square, &a, &square, &b, 1.0e16, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);
    CHECK(result.iterations == 0U);
}

static void test_plane_contact_and_ccd(void)
{
    const ConvexPolygon square = axis_aligned_square();
    const ConvexPolygonPlane plane = { { 1.0, 0.0 }, 3.0 };
    const ConvexPolygonCcdSettings settings = tight_ccd_settings();
    ConvexPolygonMotion body = motion(
        0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    ConvexPolygonContact contact;
    ConvexPolygonCcdResult result;

    CHECK(convex_polygon_plane_contact(
        &square, &body.pose, &plane, &contact));
    CHECK_NEAR(contact.separation, 3.0 - M_SQRT1_2, 2.0e-14);
    CHECK_NEAR(contact.normal.x, 1.0, 2.0e-14);

    CHECK(convex_polygon_plane_time_of_impact(
        &square, &body, &plane, 4.0, &settings, &result)
        == CONVEX_POLYGON_CCD_HIT);
    CHECK_NEAR(result.time, 3.0 - M_SQRT1_2, 2.0e-10);

    body.velocity.x = -1.0;
    CHECK(convex_polygon_plane_time_of_impact(
        &square, &body, &plane, 4.0, &settings, &result)
        == CONVEX_POLYGON_CCD_NO_HIT);

    body = motion(0.0, 0.0, 0.0, 1.0, 1.0e6, 0.0, 0.0);
    CHECK(convex_polygon_plane_time_of_impact(
        &square, &body, &plane, 1.0e16, NULL, &result)
        == CONVEX_POLYGON_CCD_HIT);
    CHECK_NEAR(result.time,
               (3.0 - M_SQRT1_2) / 1.0e6,
               2.0e-15);
}

static void test_pair_impulses(void)
{
    ConvexPolygonContact contact;
    ConvexPolygonRigidBody a;
    ConvexPolygonRigidBody b;
    double energy_before;
    double momentum_x_before;
    double normal_velocity_before;
    double impulse;

    memset(&contact, 0, sizeof(contact));
    contact.normal = v2(1.0, 0.0);
    contact.point = v2(0.0, 0.0);

    a.position = v2(-1.0, 0.0);
    a.velocity = v2(1.0, 0.0);
    a.angular_velocity = 0.0;
    a.mass = 1.0;
    a.inertia = 1.0;
    b.position = v2(1.0, 0.0);
    b.velocity = v2(-1.0, 0.0);
    b.angular_velocity = 0.0;
    b.mass = 1.0;
    b.inertia = 1.0;
    energy_before = body_energy(&a) + body_energy(&b);
    momentum_x_before = a.mass * a.velocity.x + b.mass * b.velocity.x;
    normal_velocity_before = convex_polygon_contact_normal_velocity(
        &a, &b, &contact);
    impulse = convex_polygon_apply_pair_impulse(&a, &b, &contact, 1.0);
    CHECK(impulse > 0.0);
    CHECK_NEAR(a.velocity.x, -1.0, 2.0e-14);
    CHECK_NEAR(b.velocity.x, 1.0, 2.0e-14);
    CHECK_NEAR(a.mass * a.velocity.x + b.mass * b.velocity.x,
               momentum_x_before, 2.0e-14);
    CHECK_NEAR(body_energy(&a) + body_energy(&b), energy_before, 2.0e-14);
    CHECK_NEAR(convex_polygon_contact_normal_velocity(&a, &b, &contact),
               -normal_velocity_before, 2.0e-14);

    contact.point = v2(0.0, 0.5);
    a.velocity = v2(1.0, 0.0);
    b.velocity = v2(-1.0, 0.0);
    a.angular_velocity = 0.0;
    b.angular_velocity = 0.0;
    energy_before = body_energy(&a) + body_energy(&b);
    momentum_x_before = a.mass * a.velocity.x + b.mass * b.velocity.x;
    normal_velocity_before = convex_polygon_contact_normal_velocity(
        &a, &b, &contact);
    impulse = convex_polygon_apply_pair_impulse(&a, &b, &contact, 1.0);
    CHECK(impulse > 0.0);
    CHECK(a.angular_velocity > 0.0);
    CHECK(b.angular_velocity < 0.0);
    CHECK_NEAR(a.mass * a.velocity.x + b.mass * b.velocity.x,
               momentum_x_before, 2.0e-14);
    CHECK_NEAR(body_energy(&a) + body_energy(&b), energy_before, 2.0e-14);
    CHECK_NEAR(convex_polygon_contact_normal_velocity(&a, &b, &contact),
               -normal_velocity_before, 2.0e-14);

    a.velocity = v2(1.0, 0.0);
    b.velocity = v2(-1.0, 0.0);
    a.angular_velocity = 0.0;
    b.angular_velocity = 0.0;
    normal_velocity_before = convex_polygon_contact_normal_velocity(
        &a, &b, &contact);
    impulse = convex_polygon_apply_pair_impulse(&a, &b, &contact, 0.5);
    CHECK(impulse > 0.0);
    CHECK_NEAR(convex_polygon_contact_normal_velocity(&a, &b, &contact),
               -0.5 * normal_velocity_before, 2.0e-14);
}

static void test_fixed_wall_impulse(void)
{
    ConvexPolygonRigidBody body;
    const ConvexPolygonVec2 point = { 1.0, 0.5 };
    const ConvexPolygonVec2 normal = { 1.0, 0.0 };
    double energy_before;
    double point_normal_velocity_before;
    double point_normal_velocity_after;

    body.position = v2(0.0, 0.0);
    body.velocity = v2(1.0, 0.0);
    body.angular_velocity = 0.0;
    body.mass = 2.0;
    body.inertia = 1.0;
    energy_before = body_energy(&body);
    point_normal_velocity_before = body.velocity.x
        - body.angular_velocity * point.y;
    CHECK(convex_polygon_apply_fixed_wall_impulse(
        &body, point, normal, 1.0) > 0.0);
    point_normal_velocity_after = body.velocity.x
        - body.angular_velocity * point.y;
    CHECK_NEAR(point_normal_velocity_after,
               -point_normal_velocity_before, 2.0e-14);
    CHECK_NEAR(body_energy(&body), energy_before, 2.0e-14);
}

int main(void)
{
    test_regular_polygons_and_properties();
    test_validation_failures();
    test_clamped_interactive_edits();
    test_sat_separation_and_witnesses();
    test_cached_sat_matches_brute_projection();
    test_translational_ccd();
    test_rotational_ccd();
    test_homothetic_growth_ccd();
    test_scale_aware_disk_broadphase_edges();
    test_plane_contact_and_ccd();
    test_pair_impulses();
    test_fixed_wall_impulse();

    if (failures != 0) {
        fprintf(stderr, "convex polygon tests: %d failure(s)\n", failures);
        return EXIT_FAILURE;
    }
    puts("convex polygon tests: ok");
    return EXIT_SUCCESS;
}
