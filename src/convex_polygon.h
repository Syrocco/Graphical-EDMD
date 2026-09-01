#ifndef EDMD_CONVEX_POLYGON_H
#define EDMD_CONVEX_POLYGON_H

/*
 * Small, allocation-free geometry/rigid-body helpers for the optional 2-D
 * convex-polygon particle mode.  This header deliberately has no dependency
 * on EDMD.h: the disk solver need not include it, and the routines can be
 * tested independently.
 */

#ifdef __cplusplus
extern "C" {
#endif

#define CONVEX_POLYGON_MAX_VERTICES 16
#define CONVEX_POLYGON_DEFAULT_MIN_EDGE_FRACTION 1.0e-5
#define CONVEX_POLYGON_DEFAULT_MIN_TURN_SINE 1.0e-6

typedef struct ConvexPolygonVec2 {
    double x;
    double y;
} ConvexPolygonVec2;

/*
 * Vertices are stored counter-clockwise.  area, centroid, radius,
 * inertia_per_mass, edge_inverse_lengths, edge_normals, and edge projection
 * bounds are derived values filled by convex_polygon_refresh() or
 * convex_polygon_recenter_and_normalize().
 * inertia_per_mass is the polar moment about centroid for a lamina of unit
 * mass. Cached inverse edge lengths avoid square roots in repeated SAT scans.
 */
typedef struct ConvexPolygon {
    int count;
    ConvexPolygonVec2 vertices[CONVEX_POLYGON_MAX_VERTICES];
    double area;
    ConvexPolygonVec2 centroid;
    double radius;
    double inertia_per_mass;
    double edge_inverse_lengths[CONVEX_POLYGON_MAX_VERTICES];
    ConvexPolygonVec2 edge_normals[CONVEX_POLYGON_MAX_VERTICES];
    double edge_projection_minima[CONVEX_POLYGON_MAX_VERTICES];
    double edge_projection_maxima[CONVEX_POLYGON_MAX_VERTICES];
} ConvexPolygon;

typedef enum ConvexPolygonValidation {
    CONVEX_POLYGON_VALID = 0,
    CONVEX_POLYGON_INVALID_ARGUMENT,
    CONVEX_POLYGON_INVALID_VERTEX_COUNT,
    CONVEX_POLYGON_NONFINITE_VERTEX,
    CONVEX_POLYGON_DEGENERATE,
    CONVEX_POLYGON_EDGE_TOO_SHORT,
    CONVEX_POLYGON_NOT_STRICTLY_CCW_CONVEX
} ConvexPolygonValidation;

/* Validate using scale-relative edge and turn tolerances. */
ConvexPolygonValidation convex_polygon_validate_with_tolerances(
    const ConvexPolygon *polygon,
    double min_edge_fraction,
    double min_turn_sine);
ConvexPolygonValidation convex_polygon_validate(const ConvexPolygon *polygon);
const char *convex_polygon_validation_message(ConvexPolygonValidation result);

/*
 * Refresh derived values without changing vertices.  Recenter/normalize also
 * moves the area centroid to (0,0) and scales the furthest vertex to radius 1.
 * Both return nonzero on success.
 */
int convex_polygon_refresh(ConvexPolygon *polygon);
int convex_polygon_recenter_and_normalize(ConvexPolygon *polygon);
int convex_polygon_make_regular(ConvexPolygon *polygon, int vertex_count);

/*
 * Interactive-edit helpers.  The original polygon must be valid.  A requested
 * invalid edit is clamped by binary search to the last valid point on the drag
 * segment.  The return value is the accepted fraction in [0,1], or -1 on bad
 * input.  They preserve the current coordinate frame; call recenter/normalize
 * after the edit is committed for use as a rigid-body template.
 */
double convex_polygon_move_vertex_clamped(
    ConvexPolygon *polygon,
    int vertex_index,
    ConvexPolygonVec2 requested_position);
double convex_polygon_translate_edge_clamped(
    ConvexPolygon *polygon,
    int edge_index,
    ConvexPolygonVec2 requested_delta);

typedef struct ConvexPolygonPose {
    ConvexPolygonVec2 position;
    double angle;
    double scale;
} ConvexPolygonPose;

typedef struct ConvexPolygonMotion {
    ConvexPolygonPose pose;
    ConvexPolygonVec2 velocity;
    double angular_velocity;
    /* Homothetic scale rate. Zero is the ordinary rigid-body case. */
    double scale_velocity;
} ConvexPolygonMotion;

/*
 * separation is positive for a gap, zero at touching, and negative for SAT
 * overlap. normal points from A toward B. point_a/point_b are support-feature
 * witnesses and point is their midpoint, suitable for a frictionless impulse.
 */
typedef struct ConvexPolygonContact {
    double separation;
    ConvexPolygonVec2 normal;
    ConvexPolygonVec2 point_a;
    ConvexPolygonVec2 point_b;
    ConvexPolygonVec2 point;
    int axis_from_a;
    int edge_index;
} ConvexPolygonContact;

int convex_polygon_sat_contact(
    const ConvexPolygon *a,
    const ConvexPolygonPose *pose_a,
    const ConvexPolygon *b,
    const ConvexPolygonPose *pose_b,
    ConvexPolygonContact *contact);

typedef struct ConvexPolygonCcdSettings {
    /* Values <= 0 select scale/machine-precision based defaults. */
    double distance_tolerance;
    double time_tolerance;
    unsigned int max_iterations;
} ConvexPolygonCcdSettings;

typedef enum ConvexPolygonCcdStatus {
    CONVEX_POLYGON_CCD_INVALID = -1,
    CONVEX_POLYGON_CCD_ITERATION_LIMIT = -2,
    CONVEX_POLYGON_CCD_NO_HIT = 0,
    CONVEX_POLYGON_CCD_HIT = 1,
    CONVEX_POLYGON_CCD_INITIAL_OVERLAP = 2
} ConvexPolygonCcdStatus;

typedef struct ConvexPolygonCcdResult {
    ConvexPolygonCcdStatus status;
    double time;
    unsigned int iterations;
    ConvexPolygonContact contact;
} ConvexPolygonCcdResult;

void convex_polygon_ccd_default_settings(ConvexPolygonCcdSettings *settings);

/* Shared numerical classification for an initial touching contact.  Runtime
 * response code must use the same threshold as CCD so a contact cannot be
 * predicted as an immediate hit and then rejected as nonapproaching forever.
 */
double convex_polygon_ccd_velocity_tolerance(double speed_bound);

/*
 * Conservative advancement over [0,horizon].  Poses advance with constant
 * linear/angular velocity.  Circumscribed disks are used only as a safe
 * broadphase.  A separating contact at t=0 is ignored; an approaching one is
 * reported.  The returned enum is also stored in result->status.
 */
ConvexPolygonCcdStatus convex_polygon_time_of_impact(
    const ConvexPolygon *a,
    const ConvexPolygonMotion *motion_a,
    const ConvexPolygon *b,
    const ConvexPolygonMotion *motion_b,
    double horizon,
    const ConvexPolygonCcdSettings *settings,
    ConvexPolygonCcdResult *result);

/* A fixed plane bounds dot(x, normal) <= offset; normal points at the wall. */
typedef struct ConvexPolygonPlane {
    ConvexPolygonVec2 normal;
    double offset;
} ConvexPolygonPlane;

int convex_polygon_plane_contact(
    const ConvexPolygon *polygon,
    const ConvexPolygonPose *pose,
    const ConvexPolygonPlane *plane,
    ConvexPolygonContact *contact);

ConvexPolygonCcdStatus convex_polygon_plane_time_of_impact(
    const ConvexPolygon *polygon,
    const ConvexPolygonMotion *motion,
    const ConvexPolygonPlane *plane,
    double horizon,
    const ConvexPolygonCcdSettings *settings,
    ConvexPolygonCcdResult *result);

typedef struct ConvexPolygonRigidBody {
    ConvexPolygonVec2 position;
    ConvexPolygonVec2 velocity;
    double angular_velocity;
    double mass;
    double inertia;
} ConvexPolygonRigidBody;

/* Relative normal contact velocity v_B(point)-v_A(point), dotted with A->B. */
double convex_polygon_contact_normal_velocity(
    const ConvexPolygonRigidBody *a,
    const ConvexPolygonRigidBody *b,
    const ConvexPolygonContact *contact);

/* Return the nonnegative impulse magnitude actually applied, or zero. */
double convex_polygon_apply_pair_impulse(
    ConvexPolygonRigidBody *a,
    ConvexPolygonRigidBody *b,
    const ConvexPolygonContact *contact,
    double restitution);

/*
 * Fixed-wall variant. wall_normal points from the body toward the wall and
 * contact_point is in world coordinates.
 */
double convex_polygon_apply_fixed_wall_impulse(
    ConvexPolygonRigidBody *body,
    ConvexPolygonVec2 contact_point,
    ConvexPolygonVec2 wall_normal,
    double restitution);

#ifdef __cplusplus
}
#endif

#endif /* EDMD_CONVEX_POLYGON_H */
