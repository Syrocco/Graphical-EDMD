#ifndef EDMD_CAMERA_FRUSTUM_H
#define EDMD_CAMERA_FRUSTUM_H

#include "raylib.h"

#include <stdbool.h>

typedef enum CameraFrustumPlane3DIndex {
	CAMERA_FRUSTUM_PLANE_NEAR = 0,
	CAMERA_FRUSTUM_PLANE_FAR,
	CAMERA_FRUSTUM_PLANE_LEFT,
	CAMERA_FRUSTUM_PLANE_RIGHT,
	CAMERA_FRUSTUM_PLANE_BOTTOM,
	CAMERA_FRUSTUM_PLANE_TOP,
	CAMERA_FRUSTUM_PLANE_COUNT
} CameraFrustumPlane3DIndex;

/* Plane normals point into the visible volume.  Each plane is normalized and
 * uses the equation dot(normal, point) + offset >= 0 for points inside. */
typedef struct CameraFrustumPlane3D {
	double normal_x;
	double normal_y;
	double normal_z;
	double offset;
} CameraFrustumPlane3D;

typedef struct CameraFrustum3D {
	CameraFrustumPlane3D planes[CAMERA_FRUSTUM_PLANE_COUNT];
	bool valid;
} CameraFrustum3D;

/* Build the same symmetric perspective frustum used by raylib's BeginMode3D.
 * `aspect_ratio` must be the active framebuffer width divided by its height;
 * for EDMD's offscreen scene this is sceneTargetWidth/sceneTargetHeight.
 * `near_distance` and `far_distance` should match rlGetCullDistanceNear() and
 * rlGetCullDistanceFar().
 *
 * On malformed or degenerate input this returns false and leaves `frustum`
 * invalid.  Intersection functions deliberately fail open for an invalid
 * frustum, ensuring that an optimization can never make geometry disappear. */
bool camera_frustum_3d_build_perspective(CameraFrustum3D* frustum,
	const Camera3D* camera, double aspect_ratio, double near_distance,
	double far_distance);

/* Conservative sphere/frustum intersection.  A sphere tangent to any plane
 * remains visible.  Invalid/non-finite geometry also fails open. */
bool camera_frustum_3d_intersects_sphere(const CameraFrustum3D* frustum,
	Vector3 center, double radius);

/* Conservative bound for an arrow or other swept segment.  The geometry must
 * lie within `radial_extent` of the segment from start to tip.  Internally it
 * uses the bounding sphere
 *
 *   center = (start + tip)/2
 *   radius = hypot(length(start, tip)/2, radial_extent)
 *
 * which includes a cylinder shaft and cone head.  Invalid input fails open. */
bool camera_frustum_3d_intersects_arrow(const CameraFrustum3D* frustum,
	Vector3 start, Vector3 tip, double radial_extent);

#endif
