#include "camera_frustum.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void expect_true(const char* name, bool value){
	if (value) return;
	fprintf(stderr, "FAIL: %s\n", name);
	failures++;
}

static void expect_false(const char* name, bool value){
	expect_true(name, !value);
}

static Camera3D forward_z_camera(void){
	return (Camera3D){
		.position = {0, 0, 0},
		.target = {0, 0, 1},
		.up = {0, 1, 0},
		.fovy = 90,
		.projection = CAMERA_PERSPECTIVE,
	};
}

static CameraFrustum3D unit_frustum(void){
	CameraFrustum3D frustum;
	Camera3D camera = forward_z_camera();
	expect_true("build unit perspective frustum",
		camera_frustum_3d_build_perspective(&frustum, &camera, 1.0, 1.0,
			10.0));
	return frustum;
}

static void test_center_and_six_planes(void){
	const CameraFrustum3D frustum = unit_frustum();
	expect_true("center point is visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 0, 5}, 0));

	expect_true("left tangent point remains visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){-5, 0, 5}, 0));
	expect_false("point beyond left plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){-5.01f, 0, 5}, 0));
	expect_true("right tangent point remains visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){5, 0, 5}, 0));
	expect_false("point beyond right plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){5.01f, 0, 5}, 0));
	expect_true("bottom tangent point remains visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, -5, 5}, 0));
	expect_false("point beyond bottom plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, -5.01f, 5}, 0));
	expect_true("top tangent point remains visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 5, 5}, 0));
	expect_false("point beyond top plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 5.01f, 5}, 0));

	expect_true("sphere tangent to near plane remains visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 0, 0.5f}, 0.5));
	expect_false("sphere wholly before near plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 0, 0.49f}, 0.5));
	expect_true("sphere tangent to far plane remains visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 0, 10.5f}, 0.5));
	expect_false("sphere wholly beyond far plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 0, 10.51f}, 0.5));
}

static void test_radius_expands_side_planes(void){
	const CameraFrustum3D frustum = unit_frustum();
	expect_true("side-straddling sphere is visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){5.5f, 0, 5}, 0.4));
	expect_false("sphere wholly beyond side plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){5.7f, 0, 5}, 0.4));
	expect_true("sphere crossing near plane from behind eye is visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 0, -0.1f}, 1.1));
	expect_false("sphere wholly behind near plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){0, 0, -0.2f}, 1.1));
}

static void test_aspect_and_orbit(void){
	CameraFrustum3D frustum;
	Camera3D camera = forward_z_camera();
	expect_true("build wide frustum",
		camera_frustum_3d_build_perspective(&frustum, &camera, 2.0, 1.0,
			10.0));
	expect_true("wide horizontal tangent remains visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){10, 0, 5}, 0));
	expect_false("point beyond wide horizontal plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){10.01f, 0, 5}, 0));

	camera.position = (Vector3){10, 0, 0};
	camera.target = (Vector3){5, 0, 0};
	expect_true("build rotated frustum",
		camera_frustum_3d_build_perspective(&frustum, &camera, 1.0, 1.0,
			10.0));
	expect_true("rotated forward point is visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){5, 0, 0}, 0));
	expect_true("rotated right tangent remains visible",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){5, 0, -5}, 0));
	expect_false("point beyond rotated right plane is rejected",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){5, 0, -5.01f}, 0));
}

static void test_arrow_bound(void){
	const CameraFrustum3D frustum = unit_frustum();
	expect_false("offscreen particle sphere is culled",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){6, 0, 5}, 0.1));
	expect_true("arrow entering the view remains visible",
		camera_frustum_3d_intersects_arrow(&frustum,
			(Vector3){6, 0, 5}, (Vector3){4, 0, 5}, 0.1));
	expect_false("wholly offscreen arrow is culled",
		camera_frustum_3d_intersects_arrow(&frustum,
			(Vector3){7, 0, 5}, (Vector3){6, 0, 5}, 0.1));
	expect_true("zero-length arrow uses its radial bound",
		camera_frustum_3d_intersects_arrow(&frustum,
			(Vector3){5.1f, 0, 5}, (Vector3){5.1f, 0, 5}, 0.2));
}

static void test_fail_open(void){
	CameraFrustum3D frustum;
	Camera3D camera = forward_z_camera();
	camera.target = camera.position;
	expect_false("degenerate forward vector rejects build",
		camera_frustum_3d_build_perspective(&frustum, &camera, 1.0, 1.0,
			10.0));
	expect_true("invalid frustum fails open",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){100, 100, 100}, 0));
	expect_true("null frustum fails open",
		camera_frustum_3d_intersects_sphere(NULL, (Vector3){100, 100, 100}, 0));

	camera = forward_z_camera();
	camera.up = (Vector3){0, 0, 1};
	expect_false("parallel up vector rejects build",
		camera_frustum_3d_build_perspective(&frustum, &camera, 1.0, 1.0,
			10.0));
	camera = forward_z_camera();
	expect_false("zero aspect rejects build",
		camera_frustum_3d_build_perspective(&frustum, &camera, 0, 1.0, 10.0));
	expect_false("reversed clip distances reject build",
		camera_frustum_3d_build_perspective(&frustum, &camera, 1.0, 10.0, 1.0));
	camera.projection = CAMERA_ORTHOGRAPHIC;
	expect_false("orthographic camera rejects perspective build",
		camera_frustum_3d_build_perspective(&frustum, &camera, 1.0, 1.0,
			10.0));

	frustum = unit_frustum();
	expect_true("NaN center fails open",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){NAN, 0, 5}, 0));
	expect_true("negative radius fails open",
		camera_frustum_3d_intersects_sphere(&frustum, (Vector3){100, 0, 5}, -1));
	expect_true("invalid arrow radius fails open",
		camera_frustum_3d_intersects_arrow(&frustum,
			(Vector3){100, 0, 5}, (Vector3){101, 0, 5}, NAN));
}

int main(void){
	test_center_and_six_planes();
	test_radius_expands_side_planes();
	test_aspect_and_orbit();
	test_arrow_bound();
	test_fail_open();
	if (failures != 0){
		fprintf(stderr, "%d camera-frustum test(s) failed\n", failures);
		return EXIT_FAILURE;
	}
	puts("camera-frustum tests passed");
	return EXIT_SUCCESS;
}
