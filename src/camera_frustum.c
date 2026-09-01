#include "camera_frustum.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

typedef struct FrustumVector3D {
	double x;
	double y;
	double z;
} FrustumVector3D;

static bool vector_is_finite(FrustumVector3D value){
	return isfinite(value.x) && isfinite(value.y) && isfinite(value.z);
}

static FrustumVector3D vector_from_raylib(Vector3 value){
	return (FrustumVector3D){value.x, value.y, value.z};
}

static FrustumVector3D vector_subtract(FrustumVector3D left,
	FrustumVector3D right){
	return (FrustumVector3D){
		left.x - right.x,
		left.y - right.y,
		left.z - right.z,
	};
}

static FrustumVector3D vector_scale(FrustumVector3D value, double scale){
	return (FrustumVector3D){
		value.x*scale,
		value.y*scale,
		value.z*scale,
	};
}

static FrustumVector3D vector_add(FrustumVector3D left,
	FrustumVector3D right){
	return (FrustumVector3D){
		left.x + right.x,
		left.y + right.y,
		left.z + right.z,
	};
}

static double vector_dot(FrustumVector3D left, FrustumVector3D right){
	return left.x*right.x + left.y*right.y + left.z*right.z;
}

static FrustumVector3D vector_cross(FrustumVector3D left,
	FrustumVector3D right){
	return (FrustumVector3D){
		left.y*right.z - left.z*right.y,
		left.z*right.x - left.x*right.z,
		left.x*right.y - left.y*right.x,
	};
}

static bool vector_normalize(FrustumVector3D value,
	FrustumVector3D* normalized){
	if (normalized == NULL || !vector_is_finite(value)) return false;
	const double length_squared = vector_dot(value, value);
	if (!isfinite(length_squared) || length_squared <= DBL_MIN) return false;
	const double inverse_length = 1.0/sqrt(length_squared);
	if (!isfinite(inverse_length)) return false;
	*normalized = vector_scale(value, inverse_length);
	return vector_is_finite(*normalized);
}

static CameraFrustumPlane3D plane_from_normal_and_eye(
	FrustumVector3D normal, FrustumVector3D eye, double distance_offset){
	return (CameraFrustumPlane3D){
		.normal_x = normal.x,
		.normal_y = normal.y,
		.normal_z = normal.z,
		.offset = -vector_dot(normal, eye) + distance_offset,
	};
}

bool camera_frustum_3d_build_perspective(CameraFrustum3D* frustum,
	const Camera3D* camera, double aspect_ratio, double near_distance,
	double far_distance){
	if (frustum == NULL) return false;
	memset(frustum, 0, sizeof(*frustum));
	if (camera == NULL || camera->projection != CAMERA_PERSPECTIVE ||
		!isfinite(aspect_ratio) || aspect_ratio <= 0 ||
		!isfinite(near_distance) || near_distance <= 0 ||
		!isfinite(far_distance) || far_distance <= near_distance ||
		!isfinite(camera->fovy) || camera->fovy <= 0 || camera->fovy >= 180){
		return false;
	}

	const FrustumVector3D eye = vector_from_raylib(camera->position);
	const FrustumVector3D target = vector_from_raylib(camera->target);
	const FrustumVector3D camera_up = vector_from_raylib(camera->up);
	if (!vector_is_finite(eye) || !vector_is_finite(target) ||
		!vector_is_finite(camera_up)){
		return false;
	}

	FrustumVector3D forward;
	FrustumVector3D right;
	FrustumVector3D up;
	if (!vector_normalize(vector_subtract(target, eye), &forward) ||
		!vector_normalize(vector_cross(forward, camera_up), &right) ||
		!vector_normalize(vector_cross(right, forward), &up)){
		return false;
	}

	const double half_vertical = 0.5*(double)camera->fovy*DEG2RAD;
	const double half_horizontal = atan(tan(half_vertical)*aspect_ratio);
	const double sin_vertical = sin(half_vertical);
	const double cos_vertical = cos(half_vertical);
	const double sin_horizontal = sin(half_horizontal);
	const double cos_horizontal = cos(half_horizontal);
	if (!isfinite(half_horizontal) || !isfinite(sin_vertical) ||
		!isfinite(cos_vertical) || !isfinite(sin_horizontal) ||
		!isfinite(cos_horizontal)){
		return false;
	}

	frustum->planes[CAMERA_FRUSTUM_PLANE_NEAR] =
		plane_from_normal_and_eye(forward, eye, -near_distance);
	frustum->planes[CAMERA_FRUSTUM_PLANE_FAR] =
		plane_from_normal_and_eye(vector_scale(forward, -1), eye,
			far_distance);
	frustum->planes[CAMERA_FRUSTUM_PLANE_LEFT] =
		plane_from_normal_and_eye(vector_add(vector_scale(right, cos_horizontal),
			vector_scale(forward, sin_horizontal)), eye, 0);
	frustum->planes[CAMERA_FRUSTUM_PLANE_RIGHT] =
		plane_from_normal_and_eye(vector_add(vector_scale(right, -cos_horizontal),
			vector_scale(forward, sin_horizontal)), eye, 0);
	frustum->planes[CAMERA_FRUSTUM_PLANE_BOTTOM] =
		plane_from_normal_and_eye(vector_add(vector_scale(up, cos_vertical),
			vector_scale(forward, sin_vertical)), eye, 0);
	frustum->planes[CAMERA_FRUSTUM_PLANE_TOP] =
		plane_from_normal_and_eye(vector_add(vector_scale(up, -cos_vertical),
			vector_scale(forward, sin_vertical)), eye, 0);

	for (int i = 0; i < CAMERA_FRUSTUM_PLANE_COUNT; i++){
		const CameraFrustumPlane3D* plane = frustum->planes + i;
		if (!isfinite(plane->normal_x) || !isfinite(plane->normal_y) ||
			!isfinite(plane->normal_z) || !isfinite(plane->offset)){
			memset(frustum, 0, sizeof(*frustum));
			return false;
		}
	}
	frustum->valid = true;
	return true;
}

static bool intersects_sphere_components(const CameraFrustum3D* frustum,
	double center_x, double center_y, double center_z, double radius){
	if (frustum == NULL || !frustum->valid || !isfinite(center_x) ||
		!isfinite(center_y) || !isfinite(center_z) || !isfinite(radius) ||
		radius < 0){
		return true;
	}

	for (int i = 0; i < CAMERA_FRUSTUM_PLANE_COUNT; i++){
		const CameraFrustumPlane3D* plane = frustum->planes + i;
		const double term_x = plane->normal_x*center_x;
		const double term_y = plane->normal_y*center_y;
		const double term_z = plane->normal_z*center_z;
		const double signed_distance =
			term_x + term_y + term_z + plane->offset;
		if (!isfinite(signed_distance)) return true;
		/* Rendering uses float positions and matrices.  Expanding every plane by
		 * a few float ulps prevents a mathematically tangent sphere from
		 * flickering out because CPU and GPU evaluation round differently. */
		const double scale = 1.0 + fabs(term_x) + fabs(term_y) +
			fabs(term_z) + fabs(plane->offset) + radius;
		const double tolerance = 16.0*FLT_EPSILON*scale;
		if (signed_distance < -radius - tolerance) return false;
	}
	return true;
}

bool camera_frustum_3d_intersects_sphere(const CameraFrustum3D* frustum,
	Vector3 center, double radius){
	return intersects_sphere_components(frustum, center.x, center.y, center.z,
		radius);
}

bool camera_frustum_3d_intersects_arrow(const CameraFrustum3D* frustum,
	Vector3 start, Vector3 tip, double radial_extent){
	const FrustumVector3D segment_start = vector_from_raylib(start);
	const FrustumVector3D segment_tip = vector_from_raylib(tip);
	if (!vector_is_finite(segment_start) || !vector_is_finite(segment_tip) ||
		!isfinite(radial_extent) || radial_extent < 0){
		return true;
	}
	const FrustumVector3D difference =
		vector_subtract(segment_tip, segment_start);
	const double length_squared = vector_dot(difference, difference);
	if (!isfinite(length_squared) || length_squared < 0) return true;
	const double half_length = 0.5*sqrt(length_squared);
	const double bound_radius = hypot(half_length, radial_extent);
	const FrustumVector3D midpoint = vector_add(segment_start,
		vector_scale(difference, 0.5));
	return intersects_sphere_components(frustum, midpoint.x, midpoint.y,
		midpoint.z, bound_radius);
}
