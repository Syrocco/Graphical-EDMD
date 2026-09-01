#include "graphics.h"
#include "color.h"
#include "gui_imgui.h"
#include "live_analysis.h"
#include "sim_lifecycle.h"
#include "EDMD.h"
#if POLYGON_PARTICLES && !THREE_D
#include "polygon_particles.h"
#endif
#include "raymath.h"
#if THREE_D
#include "camera_frustum.h"
#include "rlgl.h"
#endif
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#if !defined(PLATFORM_WEB)
#include <pthread.h>
#else
#include <emscripten/emscripten.h>
#include <emscripten/html5.h>
#endif
#include <stdlib.h>
#include <string.h>

extern double E;
extern double pressure;
extern double tmax;
extern const int suppSizeEvent;
#if THREE_D
extern int Nzcells;
#endif

#define qN 180
#define SIZE 400
#define NUM_THREADS 4
#define PARTICLE_DRAG_THRESHOLD_PIXELS 5.0f
#define PARTICLE_DRAG_REFERENCE_FPS 60.0
#define PARTICLE_DRAG_VELOCITY_BLEND 0.30
#define STRUCTURE_FACTOR_DEFAULT_SAMPLE_CAP 4096
#define STRUCTURE_FACTOR_HARD_SAMPLE_CAP 8192
#define STRUCTURE_FACTOR_DEFAULT_REFRESH_SECONDS 0.5
#define COLOR_DEFAULT_REFRESH_SECONDS (1.0/30.0)
#define RENDER_REFRESH_MIN_SECONDS (1.0/1000.0)
#define RENDER_REFRESH_MAX_SECONDS 3600.0
#if !THREE_D
int nSym = 4;
#endif

double qxG[qN] = {0};
float structFactorGraphics[qN][qN] = {0};
bool structFactorActivated = 0;
position* positions;
#if !defined(PLATFORM_WEB)
pthread_t threads[NUM_THREADS];
threadArg threadArgs[NUM_THREADS];
#endif

#if !THREE_D
typedef struct circleBatch circleBatch;
struct circleBatch{
	Mesh mesh;
	int capacity;
	int count;
	bool loaded;
};

static circleBatch particleCircleBatch = {0};
enum { VELOCITY_ARROW_VERTICES = 18 };
typedef struct velocityArrowBatch velocityArrowBatch;
struct velocityArrowBatch{
	Mesh mesh;
	int capacity;
	int count;
	bool loaded;
};
static velocityArrowBatch arrowBatch = {0};
static Material arrowMaterial = {0};
static bool arrowMaterialLoaded = false;
static Material proceduralCircleMaterial = {0};
static bool proceduralCircleMaterialLoaded = false;
#if POLYGON_PARTICLES
typedef struct polygonTriangleBatch polygonTriangleBatch;
struct polygonTriangleBatch{
	Mesh mesh;
	Material material;
	int capacity;
	int shapeVertexCount;
	int verticesPerParticle;
	int count;
	bool meshLoaded;
	bool materialLoaded;
};
static polygonTriangleBatch particlePolygonBatch = {0};
#endif
#endif
#if THREE_D
typedef struct sphereBatch sphereBatch;
struct sphereBatch{
	Mesh mesh;
	Material material;
	Shader shader;
	Matrix* transforms;
	int lightingModeLocation;
	int cameraPositionLocation;
	int backgroundColorLocation;
	int sceneRadiusLocation;
	int capacity;
	int count;
	int rings;
	int slices;
	bool meshLoaded;
	bool materialLoaded;
	bool shaderLoaded;
	bool loaded;
};
static sphereBatch particleSphereBatch = {0};
#endif
static particlePrediction* renderPredictions = NULL;
static int renderPredictionCapacity = 0;
#if THREE_D
static double maximumPredictedParticleRadius3D = 0;
#else
static double maximumPredictedParticleRadius = 0;
static double maximumPredictedParticleRadiusTime = NAN;
static int maximumPredictedParticleRadiusCount = -1;
#endif

typedef struct particleColorCache particleColorCache;
struct particleColorCache{
	bool valid;
	int mode;
	int particleCount;
	double simulationTime;
	double wallTime;
	double min;
	double max;
};

static particleColorCache colorCache = {0};

#if !defined(PLATFORM_WEB)
static pthread_mutex_t structMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t structJobReady = PTHREAD_COND_INITIALIZER;
static pthread_cond_t structJobFinished = PTHREAD_COND_INITIALIZER;
#else
static int structWebNextRow = 0;
#endif
static unsigned long structJobGeneration = 0;
#if !defined(PLATFORM_WEB)
static int structWorkersRemaining = 0;
#endif
static int structSampleCount = 0;
static int structSampleCapacity = 0;
static int structRejectedSampleCount = -1;
static bool structJobRunning = false;
static bool structJobCancelRequested = false;
static bool structThreadPoolStopping = false;
static bool structThreadPoolStarted = false;
static double lastStructLaunchTime = -STRUCTURE_FACTOR_DEFAULT_REFRESH_SECONDS;

Image image;
Texture2D texture;

static RenderTexture2D sceneTarget = {0};
static int sceneTargetWidth = 0;
static int sceneTargetHeight = 0;
static GuiViewport guiViewport = {0};
static bool imguiInitialized = false;
static char pendingRuntimeError[1024] = {0};
static bool pendingRuntimeFault = false;
static char pendingTransientNotice[512] = {0};
/* Keep the declarative configuration separate from mutable engine working
 * state.  Recapturing derived runtime values as the user's applied
 * configuration could make an unrelated slider edit silently reset or
 * misclassify that state. */
static SimulationConfig hostAppliedConfig = {0};
static bool hostAppliedConfigValid = false;
static GuiLiveState guiLive = {
	.running = true,
	.pace_simulation_time = false,
	.target_simulation_time_rate = 1.0,
	.show_scene = true,
	.show_fps = true,
	.show_help = true,
	.show_axes = true,
	.show_boundary = true,
	.show_selection_outline = true,
	.color_mode = SIM_COLOR_BY_TYPE,
	.palette = SIM_PALETTE_PLASMA,
	.lighting_mode = SIM_LIGHTING_STUDIO,
	.sphere_quality = SIM_SPHERE_QUALITY_AUTO,
	.target_fps = 144,
	.scene_supersampling = 1.0,
	.camera_fov_degrees = 45.0,
	.structure_factor_sample_cap = STRUCTURE_FACTOR_DEFAULT_SAMPLE_CAP,
	.structure_factor_refresh_seconds = STRUCTURE_FACTOR_DEFAULT_REFRESH_SECONDS,
	.color_refresh_seconds = COLOR_DEFAULT_REFRESH_SECONDS,
	.type0_rgba = 0x4296faffu,
	.type1_rgba = 0xf29d49ffu,
	.background_rgba = 0x181a1fffu,
	.boundary_rgba = 0x8b949effu,
};
static GuiMetrics guiMetrics = {0};
static double lastMetricsWallTime = -1.0;
static double lastMetricsSimulationTime = 0.0;
static unsigned long lastMetricsEvents = 0;
#if !defined(PLATFORM_WEB)
static int appliedTargetFps = 144;
#endif

typedef struct simulationPacingState simulationPacingState;
struct simulationPacingState{
	bool initialized;
	bool running;
	bool enabled;
	double rate;
	double wallAnchor;
	double simulationAnchor;
	double targetSimulationTime;
	double lag;
	double lastSimulationTime;
};

static simulationPacingState pacingState = {0};

static double sanitizedSimulationPacingRate(void){
	double rate = guiLive.target_simulation_time_rate;
	if (!isfinite(rate) || rate <= 0)
		rate = 1.0;
	return fmin(fmax(rate, 1e-6), 1e9);
}

static void synchronizeSimulationPacing(bool running,
	double currentSimulationTime){
	const double now = GetTime();
	const bool enabled = guiLive.pace_simulation_time;
	const double rate = sanitizedSimulationPacingRate();
	const double timeTolerance = 64*DBL_EPSILON*
		(fabs(currentSimulationTime) + fabs(pacingState.lastSimulationTime) + 1);
	const bool simulationRestarted = pacingState.initialized &&
		currentSimulationTime < pacingState.lastSimulationTime - timeTolerance;
	const bool modeChanged = !pacingState.initialized ||
		pacingState.running != running || pacingState.enabled != enabled ||
		pacingState.rate != rate || simulationRestarted;

	if (modeChanged || !running){
		pacingState.wallAnchor = now;
		pacingState.simulationAnchor = currentSimulationTime;
	}
	pacingState.initialized = true;
	pacingState.running = running;
	pacingState.enabled = enabled;
	pacingState.rate = rate;
	pacingState.lastSimulationTime = currentSimulationTime;

	if (!running || !enabled){
		pacingState.targetSimulationTime = currentSimulationTime;
		pacingState.lag = 0;
		return;
	}
	const double wallElapsed = fmax(0.0, now - pacingState.wallAnchor);
	const double target = pacingState.simulationAnchor + rate*wallElapsed;
	if (!isfinite(target)){
		/* A malformed clock/rate must never lock event processing. */
		pacingState.enabled = false;
		pacingState.targetSimulationTime = currentSimulationTime;
		pacingState.lag = 0;
		return;
	}
	pacingState.targetSimulationTime = target;
	pacingState.lag = fmax(0.0, target - currentSimulationTime);
}

typedef struct particleDrag particleDrag;
struct particleDrag{
	bool active;
	int index;
	double lastWallTime;
	#if THREE_D
	Vector3 planePoint;
	Vector3 planeNormal;
	#endif
};

typedef struct particlePress particlePress;
struct particlePress{
	bool pending;
	int index;
	Vector2 screenPoint;
};

static particleDrag activeParticleDrag = {.index = -1};
static particlePress pendingParticlePress = {.index = -1};
#if THREE_D
static bool camera3DNeedsReset = true;
#else
static bool camera2DNeedsReset = true;
#endif

static particlePrediction graphicsParticlePrediction(const particle* p);
#if !THREE_D
static double graphicsMaximumParticleRadius(void);
#endif
static void invalidateParticleColorCache(void);

#if !THREE_D
/* Return the smallest symmetric cell stencil that can contain a point within
 * cutoff of any point in the center cell. */
static int graphicsCellNeighborRange(double cutoff, double boxLength, int cellCount){
	if (!(cutoff > 0) || !(boxLength > 0) || (cellCount <= 1)){
		return 0;
	}
	const double range = ceil(cutoff*cellCount/boxLength);
	if (!isfinite(range) || (range >= cellCount)){
		return cellCount;
	}
	return (int)range;
}

/* Consecutive periodic offsets are unique until their span covers the whole
 * axis.  Enumerating all coordinates in that case also avoids out-of-bounds
 * one-step wrapping in boxes with only one or two cells. */
static int graphicsNeighborCellCount(int cellCount, int range){
	const long long stencilWidth = 2LL*range + 1;
	return (stencilWidth >= cellCount) ? cellCount : (int)stencilWidth;
}

static int graphicsNeighborCellCoordinate(int center, int cellCount,
	int range, int index){
	if (graphicsNeighborCellCount(cellCount, range) == cellCount){
		return index;
	}
	long long coordinate = (long long)center - range + index;
	coordinate %= cellCount;
	if (coordinate < 0) coordinate += cellCount;
	return (int)coordinate;
}
#endif

#if THREE_D
static void resetCamera3D(window* screenWindow);
static void ensureCamera3D(window* screenWindow);
static void drawParticlesAndBox3D(const Camera3D* camera,
	const CameraFrustum3D* frustum, state* screenState,
	bool predictionsReady);
#else
static void resetCamera2D(window* screenWindow);
#endif
static bool ensureSceneTarget(int width, int height);
static bool screenPointToScene(Vector2 screenPoint, Vector2* scenePoint);

#if !THREE_D
static void graphicsPBCXY(double* dx, double* dy){
	PBC(dx, dy);
}
#endif

static bool screenPointToScene(Vector2 screenPoint, Vector2* scenePoint){
	if (scenePoint == NULL || sceneTargetWidth <= 0 || sceneTargetHeight <= 0 ||
		guiViewport.width <= 0 || guiViewport.height <= 0){
		return false;
	}
	if (screenPoint.x < guiViewport.x || screenPoint.y < guiViewport.y ||
		screenPoint.x >= guiViewport.x + guiViewport.width ||
		screenPoint.y >= guiViewport.y + guiViewport.height){
		return false;
	}
	scenePoint->x = (screenPoint.x - guiViewport.x)*
		(float)sceneTargetWidth/guiViewport.width;
	scenePoint->y = (screenPoint.y - guiViewport.y)*
		(float)sceneTargetHeight/guiViewport.height;
	return true;
}

static Color colorFromRgba(uint32_t rgba){
	return (Color){
		.r = (unsigned char)((rgba >> 24) & 0xffu),
		.g = (unsigned char)((rgba >> 16) & 0xffu),
		.b = (unsigned char)((rgba >> 8) & 0xffu),
		.a = (unsigned char)(rgba & 0xffu),
	};
}

int getWhatsUnderClick(window* screenWindow, state* screenState){
	const Vector2 mouse = GetMousePosition();
	Vector2 sceneMouse = {0};
	if (!screenPointToScene(mouse, &sceneMouse)){
		return -1;
	}

	#if THREE_D
	ensureCamera3D(screenWindow);
	const Ray ray = GetScreenToWorldRayEx(sceneMouse, screenWindow->cam3D,
		sceneTargetWidth, sceneTargetHeight);
	double closestDistance = INFINITY;
	particle* closestParticle = NULL;
	for (int i = 0; i < N; i++){
		const particlePrediction predicted = graphicsParticlePrediction(particles + i);
		const Vector3 center = {
			(float)predicted.x, (float)predicted.y, (float)predicted.z,
		};
		const RayCollision collision = GetRayCollisionSphere(
			ray, center, (float)predicted.rad);
		if (collision.hit && (collision.distance < closestDistance)){
			closestDistance = collision.distance;
			closestParticle = particles + i;
		}
	}
	if (closestParticle != NULL){
		screenState->particleUnderClick = closestParticle;
		return 1;
	}
	return 0;
	#else
	double factor = screenWindow->factor;
	Vector2 pos = GetScreenToWorld2D(sceneMouse, screenWindow->cam);
	double x = pos.x/factor;
	double y = Ly - pos.y/factor; 

	if (t > 1/vr){
		if ((addWallx) && (x > Lx - 5) && (x < Lx + 5)){
			return 2;
		}
		if ((addWally) && (y > Ly - 5) && (y < Ly + 5)){
			return 3;
		}
		}
	if ((x >= 0) && (x < Lx) && (y >= 0) && (y < Ly)){
		int X = coordToCell(x, 1);
		int Y = coordToCell(y, 0);
		const double maximumRadius = graphicsMaximumParticleRadius();
		const int rangeX = graphicsCellNeighborRange(maximumRadius, Lx, Nxcells);
		const int rangeY = graphicsCellNeighborRange(maximumRadius, Ly, Nycells);
		const int numberX = graphicsNeighborCellCount(Nxcells, rangeX);
		const int numberY = graphicsNeighborCellCount(Nycells, rangeY);

		for (int yi = 0; yi < numberY; yi++){
			const int neighborY = graphicsNeighborCellCoordinate(Y, Nycells, rangeY, yi);
			for (int xi = 0; xi < numberX; xi++){
				const int neighborX = graphicsNeighborCellCoordinate(X, Nxcells, rangeX, xi);
				screenState->particleUnderClick =
					cellList[neighborY*Nxcells + neighborX];
				while (screenState->particleUnderClick != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					const particlePrediction predicted = graphicsParticlePrediction(screenState->particleUnderClick);
					double dx = predicted.x - x;
					double dy = predicted.y - y;
					graphicsPBCXY(&dx, &dy);
					const double distanceSquared = dx*dx + dy*dy;
					const double radiusSquared = predicted.rad*predicted.rad;
					const double pickTolerance = 128.0*DBL_EPSILON*
						(1.0 + distanceSquared + radiusSquared);
					if (distanceSquared <= radiusSquared + pickTolerance){
						#if POLYGON_PARTICLES
						if (polygon_particles_requested() &&
							!polygon_particles_contains_point(
								screenState->particleUnderClick->num, x, y, t)){
							screenState->particleUnderClick =
								screenState->particleUnderClick->nxt;
							continue;
						}
						#endif
						return 1;
					}
					screenState->particleUnderClick = screenState->particleUnderClick->nxt;
				}
			}
		}
	}
	return 0;
	#endif
}

static bool ensureSceneTarget(int width, int height){
	if ((width <= 0) || (height <= 0)){
		return false;
	}
	if ((sceneTarget.texture.id != 0) &&
		(sceneTargetWidth == width) && (sceneTargetHeight == height)){
		return true;
	}
	if (sceneTarget.texture.id != 0){
		UnloadRenderTexture(sceneTarget);
		sceneTarget = (RenderTexture2D){0};
	}
	sceneTarget = LoadRenderTexture(width, height);
	if (sceneTarget.texture.id == 0){
		sceneTargetWidth = 0;
		sceneTargetHeight = 0;
		TraceLog(LOG_ERROR, "Unable to allocate the simulation render target");
		return false;
	}
	sceneTargetWidth = width;
	sceneTargetHeight = height;
	SetTextureFilter(sceneTarget.texture, TEXTURE_FILTER_BILINEAR);
	#if THREE_D
	camera3DNeedsReset = true;
	#else
	camera2DNeedsReset = true;
	#endif
	return true;
}

#if THREE_D
static float boxRadius3D(void){
	const double diagonalSquared = Lx*Lx + Ly*Ly + Lz*Lz;
	return (float)(0.5*sqrt(fmax(diagonalSquared, 1.0)));
}

typedef struct cameraClipPlanes3D {
	double nearDistance;
	double farDistance;
} cameraClipPlanes3D;

/* raylib's fixed far plane is 1000, which can cut through a large simulation
 * even though resetCamera3D has correctly moved the camera far enough away.
 * Bound the complete box, particle radii and (when enabled) the longest
 * velocity arrow, while keeping the near plane as far forward as is safely
 * possible for useful depth-buffer precision. */
static cameraClipPlanes3D computeCameraClipPlanes3D(const Camera3D* camera,
	bool includeVelocityArrows, double maximumParticleRadius){
	cameraClipPlanes3D result = {0.01, 1000.0};
	if (camera == NULL) return result;

	const double sceneRadius = boxRadius3D();
	const double centerX = 0.5*Lx;
	const double centerY = 0.5*Ly;
	const double centerZ = 0.5*Lz;
	const double dx = (double)camera->position.x - centerX;
	const double dy = (double)camera->position.y - centerY;
	const double dz = (double)camera->position.z - centerZ;
	const double cameraDistance = sqrt(dx*dx + dy*dy + dz*dz);
	if (!isfinite(sceneRadius) || sceneRadius <= 0 ||
		!isfinite(cameraDistance)){
		return result;
	}

	if (!isfinite(maximumParticleRadius) || maximumParticleRadius < 0){
		/* A failed prediction snapshot must widen, never narrow, the clip volume. */
		maximumParticleRadius = sceneRadius;
	}
	/* Particle centers lie within the box bounding sphere.  Add only the actual
	 * largest predicted radius: 1.1 covers selection wires, while 5.9 covers
	 * the longest 5.5-radius arrow and its 0.3-radius cone with margin. */
	const double particleExtent = (includeVelocityArrows ? 5.9 : 1.1)*
		maximumParticleRadius;
	const double padding = fmax(0.02*sceneRadius, 1e-3);
	const double contentRadius = sceneRadius + particleExtent + padding;
	const double minimumNear = fmax(1e-4, 1e-5*sceneRadius);
	const double nearestContent = cameraDistance - contentRadius;
	result.nearDistance = (nearestContent > 0) ?
		fmax(minimumNear, 0.5*nearestContent) : minimumNear;
	result.farDistance = fmax(result.nearDistance + fmax(sceneRadius, 1.0),
		cameraDistance + contentRadius);
	/* Projection entries are floats inside rlgl, and rlFrustum() multiplies
	 * near*far after casting both operands to float. Keeping each plane below
	 * sqrt(FLT_MAX)/4 keeps that intermediate finite even for a malformed,
	 * astronomically large scene. Ordinary large boxes still receive a far
	 * plane many orders of magnitude beyond raylib's fixed default. */
	result.farDistance = fmin(result.farDistance,
		0.25*sqrt((double)FLT_MAX));
	if (!(result.farDistance > result.nearDistance)){
		result.nearDistance = fmax(1e-6, 1e-6*result.farDistance);
	}
	return result;
}

static void resetCamera3D(window* screenWindow){
	const Vector3 center = {
		(float)(0.5*Lx), (float)(0.5*Ly), (float)(0.5*Lz),
	};
	const Vector3 viewDirection = Vector3Normalize((Vector3){1.25f, 0.9f, 1.4f});
	const float radius = boxRadius3D();
	const float verticalFov = Clamp((float)guiLive.camera_fov_degrees,
		10.0f, 120.0f)*DEG2RAD;
	const float aspect = (sceneTargetWidth > 0 && sceneTargetHeight > 0)
		? (float)sceneTargetWidth/(float)sceneTargetHeight : 1.0f;
	const float horizontalFov = 2.0f*atanf(tanf(0.5f*verticalFov)*aspect);
	const float limitingHalfFov = 0.5f*fminf(verticalFov, horizontalFov);
	/* Fit the complete box bounding sphere in both the horizontal and vertical
	   field of view.  The previous fixed 3-radius distance looked arbitrary in
	   tall/narrow dock layouts. */
	const float distance = fmaxf(
		1.12f*radius/fmaxf(sinf(limitingHalfFov), 0.08f), 5.0f);
	screenWindow->cam3D = (Camera3D){
		.position = Vector3Add(center, Vector3Scale(viewDirection, distance)),
		.target = center,
		.up = (Vector3){0.0f, 1.0f, 0.0f},
		.fovy = 45.0f,
		.projection = CAMERA_PERSPECTIVE,
	};
	screenWindow->cameraBoxLx = Lx;
	screenWindow->cameraBoxLy = Ly;
	screenWindow->cameraBoxLz = Lz;
	camera3DNeedsReset = false;
}

static void ensureCamera3D(window* screenWindow){
	if (camera3DNeedsReset || (screenWindow->cameraBoxLx <= 0) ||
		(fabs(screenWindow->cameraBoxLx - Lx) > 1e-9) ||
		(fabs(screenWindow->cameraBoxLy - Ly) > 1e-9) ||
		(fabs(screenWindow->cameraBoxLz - Lz) > 1e-9)){
		resetCamera3D(screenWindow);
	}
}

static void orbitCamera3D(window* screenWindow, Vector2 mouseDelta){
	const Vector3 worldUp = {0.0f, 1.0f, 0.0f};
	Vector3 offset = Vector3Subtract(screenWindow->cam3D.position,
		screenWindow->cam3D.target);
	offset = Vector3RotateByAxisAngle(offset, worldUp, -0.006f*mouseDelta.x);
	Vector3 right = Vector3CrossProduct(worldUp, offset);
	if (Vector3LengthSqr(right) > 1e-12f){
		right = Vector3Normalize(right);
		const Vector3 pitched = Vector3RotateByAxisAngle(
			offset, right, -0.006f*mouseDelta.y);
		const Vector3 direction = Vector3Normalize(pitched);
		if (fabsf(Vector3DotProduct(direction, worldUp)) < 0.98f){
			offset = pitched;
		}
	}
	screenWindow->cam3D.position = Vector3Add(screenWindow->cam3D.target, offset);
	screenWindow->cam3D.up = worldUp;
}

static void zoomCamera3D(window* screenWindow, float wheel){
	Vector3 offset = Vector3Subtract(screenWindow->cam3D.position,
		screenWindow->cam3D.target);
	const float oldDistance = Vector3Length(offset);
	if (oldDistance <= 1e-6f){
		resetCamera3D(screenWindow);
		return;
	}
	const float radius = boxRadius3D();
	const float minimumDistance = fmaxf(0.2f*radius, 0.1f);
	const float maximumDistance = fmaxf(2.0f*minimumDistance,
		fminf(20.0f*radius, 800.0f));
	const float newDistance = Clamp(oldDistance*powf(0.82f, wheel),
		minimumDistance, maximumDistance);
	offset = Vector3Scale(offset, newDistance/oldDistance);
	screenWindow->cam3D.position = Vector3Add(screenWindow->cam3D.target, offset);
}
#else
static void resetCamera2D(window* screenWindow){
	if (sceneTargetWidth <= 0 || sceneTargetHeight <= 0 ||
		!(Lx > 0) || !(Ly > 0)){
		return;
	}
	const double fitX = 0.90*sceneTargetWidth/Lx;
	const double fitY = 0.90*sceneTargetHeight/Ly;
	screenWindow->factor = fmax(1e-9, fmin(fitX, fitY));
	screenWindow->cam = (Camera2D){
		.offset = {(float)sceneTargetWidth*0.5f,
			(float)sceneTargetHeight*0.5f},
		.target = {(float)(0.5*Lx*screenWindow->factor),
			(float)(0.5*Ly*screenWindow->factor)},
		.rotation = 0.0f,
		.zoom = 1.0f,
	};
	camera2DNeedsReset = false;
}
#endif

static void updateWindowLayout(window* screenWindow, bool rememberWindowedSize){
	const int width = GetScreenWidth();
	const int height = GetScreenHeight();
	if (rememberWindowedSize){
		screenWindow->screenWidth = width;
		screenWindow->screenHeight = height;
	}
	screenWindow->factor = height/Ly;
	screenWindow->start = height;
	screenWindow->xGUI = width/1800.0;
	screenWindow->yGUI = height/900.0;
}

static void toggleFullscreenWindow(window* screenWindow){
	const int display = GetCurrentMonitor();
	if (IsWindowFullscreen()){
		ToggleFullscreen();
		SetWindowSize(screenWindow->screenWidth, screenWindow->screenHeight);
	}
	else{
		screenWindow->screenWidth = GetScreenWidth();
		screenWindow->screenHeight = GetScreenHeight();
		SetWindowSize(GetMonitorWidth(display), GetMonitorHeight(display));
		ToggleFullscreen();
	}
	updateWindowLayout(screenWindow, false);
}

#if THREE_D
static bool pointOnDragPlane(Ray ray, Vector3 point, Vector3 normal,
	Vector3* result){
	const float denominator = Vector3DotProduct(ray.direction, normal);
	if (result == NULL || fabsf(denominator) < 1e-6f){
		return false;
	}
	const float distance = Vector3DotProduct(
		Vector3Subtract(point, ray.position), normal)/denominator;
	if (!isfinite(distance)){
		return false;
	}
	*result = Vector3Add(ray.position, Vector3Scale(ray.direction, distance));
	return true;
}
#endif

static bool beginParticleDrag(window* screenWindow, state* screenState,
	Vector2 sceneMouse){
	particle* p = screenState->particleUnderClick;
	if (p == NULL || p->num < 0 || p->num >= N){
		return false;
	}
	if (gui_imgui_simulation_faulted()){
		return false;
	}

	/* Fix the 3D pull plane at the particle's current depth.  This makes the
	   mouse act as a damped spring in a camera-facing plane without ever
	   teleporting the particle or invalidating pair membership. */
	#if THREE_D
	const particlePrediction predicted = graphicsParticlePrediction(p);
	ensureCamera3D(screenWindow);
	const Vector3 planePoint = {
		(float)predicted.x, (float)predicted.y, (float)predicted.z};
	const Vector3 planeNormal = Vector3Normalize(Vector3Subtract(
		screenWindow->cam3D.position, screenWindow->cam3D.target));
	const Ray ray = GetScreenToWorldRayEx(sceneMouse, screenWindow->cam3D,
		sceneTargetWidth, sceneTargetHeight);
	Vector3 hit = planePoint;
	if (!pointOnDragPlane(ray, planePoint, planeNormal, &hit)){
		return false;
	}
	#else
	if (!(screenWindow->factor > 0) || !isfinite(screenWindow->factor)){
		return false;
	}
	const Vector2 world = GetScreenToWorld2D(sceneMouse, screenWindow->cam);
	const double pointerX = world.x/screenWindow->factor;
	const double pointerY = Ly - world.y/screenWindow->factor;
	if (!isfinite(pointerX) || !isfinite(pointerY)){
		return false;
	}
	#endif

	const double now = GetTime();
	activeParticleDrag = (particleDrag){
		.active = true,
		.index = p->num,
		/* Give the first update the same impulse as one 60 Hz frame. */
		.lastWallTime = now - 1.0/PARTICLE_DRAG_REFERENCE_FPS,
		#if THREE_D
		.planePoint = planePoint,
		.planeNormal = planeNormal,
		#endif
	};
	screenState->leftClicked = true;
	return true;
}

static void updateParticleDrag(window* screenWindow, Vector2 sceneMouse){
	if (!activeParticleDrag.active || activeParticleDrag.index < 0 ||
		activeParticleDrag.index >= N){
		return;
	}
	particle* p = particles + activeParticleDrag.index;
	double targetX = 0;
	double targetY = 0;
	#if THREE_D
	double targetZ = 0;
	#endif
	#if THREE_D
	const Ray ray = GetScreenToWorldRayEx(sceneMouse, screenWindow->cam3D,
		sceneTargetWidth, sceneTargetHeight);
	Vector3 hit;
	if (!pointOnDragPlane(ray, activeParticleDrag.planePoint,
		activeParticleDrag.planeNormal, &hit)) return;
	targetX = hit.x;
	targetY = hit.y;
	targetZ = hit.z;
	#else
	if (!(screenWindow->factor > 0) || !isfinite(screenWindow->factor)) return;
	const Vector2 world = GetScreenToWorld2D(sceneMouse, screenWindow->cam);
	targetX = world.x/screenWindow->factor;
	targetY = Ly - world.y/screenWindow->factor;
	#endif
	if (!isfinite(targetX) || !isfinite(targetY)
		#if THREE_D
		|| !isfinite(targetZ)
		#endif
	) return;

	#if (STRESS || HEAT)
	beginTransportVelocityBatch();
	#endif
	freeFly(p);
	const double oldVx = p->vx;
	const double oldVy = p->vy;
	#if THREE_D
	const double oldVz = p->vz;
	#else
	const double oldVz = 0;
	#endif

	double dx = targetX - p->x;
	double dy = targetY - p->y;
	#if THREE_D
	double dz = targetZ - p->z;
	PBC(&dx, &dy, &dz);
	#else
	PBC(&dx, &dy);
	#endif
	const double now = GetTime();
	double wallDt = now - activeParticleDrag.lastWallTime;
	if (!isfinite(wallDt) || wallDt <= 0){
		wallDt = 1.0/PARTICLE_DRAG_REFERENCE_FPS;
	}
	/* Clamp long UI stalls: resuming the window must not produce an enormous
	   one-frame kick.  The exponential form exactly reproduces the old 0.3
	   spring/damping blend at 60 Hz and is stable at other frame rates. */
	wallDt = fmin(wallDt, 0.05);
	const double response = 1.0 - pow(1.0 - PARTICLE_DRAG_VELOCITY_BLEND,
		wallDt*PARTICLE_DRAG_REFERENCE_FPS);
	p->vx += response*(dx - p->vx);
	p->vy += response*(dy - p->vy);
	#if THREE_D
	p->vz += response*(dz - p->vz);
	#endif
	activeParticleDrag.lastWallTime = now;

	#if (STRESS || HEAT)
	recordTransportVelocityChange(p, oldVx, oldVy, oldVz);
	#else
	(void)oldVx;
	(void)oldVy;
	(void)oldVz;
	#endif
	#if POLYGON_PARTICLES && !THREE_D
	if (!polygon_particles_refresh_growth_energy()) return;
	#endif
	/* A velocity impulse invalidates predictions but leaves coordinates, cell
	   links, potential membership, and MSD continuity intact. */
	p->coll++;
	removeEventFromQueue(eventList[p->num]);
	removeEventFromQueue(eventList[N + p->num]);
	crossingEvent(p->num);
	collisionEvent(p->num);
	physicalQ();
	invalidateParticleColorCache();
}

static void finishParticleDrag(state* screenState){
	if (!activeParticleDrag.active || activeParticleDrag.index < 0 ||
		activeParticleDrag.index >= N){
		activeParticleDrag = (particleDrag){.index = -1};
		return;
	}
	particle* p = particles + activeParticleDrag.index;
	screenState->particleUnderClick = p;
	screenState->selected = true;
	screenState->leftClicked = false;
	activeParticleDrag = (particleDrag){.index = -1};
}

static void clearPendingParticlePress(void){
	pendingParticlePress = (particlePress){.index = -1};
}

static bool pendingParticlePressMoved(Vector2 currentScreen){
	const float dx = currentScreen.x - pendingParticlePress.screenPoint.x;
	const float dy = currentScreen.y - pendingParticlePress.screenPoint.y;
	const float thresholdSquared = PARTICLE_DRAG_THRESHOLD_PIXELS*
		PARTICLE_DRAG_THRESHOLD_PIXELS;
	return dx*dx + dy*dy > thresholdSquared;
}

static bool startPendingParticleDrag(window* screenWindow, state* screenState,
	Vector2 currentScreen, bool finishImmediately){
	const int pendingIndex = pendingParticlePress.index;
	Vector2 dragMouse = {0};
	clearPendingParticlePress();
	if (!screenPointToScene(currentScreen, &dragMouse) ||
		screenState->particleUnderClick == NULL ||
		screenState->particleUnderClick->num != pendingIndex ||
		!beginParticleDrag(screenWindow, screenState, dragMouse)){
		return false;
	}
	updateParticleDrag(screenWindow, dragMouse);
	if (finishImmediately){
		finishParticleDrag(screenState);
	}
	return true;
}


void getInput(window* screenWindow, state* screenState){
	if (IsWindowResized()){
		updateWindowLayout(screenWindow, !IsWindowFullscreen());
	}

	Vector2 sceneMouse = {0};
	const bool mouseInScene = guiViewport.hovered &&
		screenPointToScene(GetMousePosition(), &sceneMouse);

	#if THREE_D
	ensureCamera3D(screenWindow);
	if (mouseInScene){
		const float wheel = GetMouseWheelMove();
		if (wheel != 0){
			zoomCamera3D(screenWindow, wheel);
		}
		if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)){
			const Vector2 mouseDelta = GetMouseDelta();
			const Vector2 sceneDelta = {
				mouseDelta.x*(float)sceneTargetWidth/guiViewport.width,
				mouseDelta.y*(float)sceneTargetHeight/guiViewport.height,
			};
			orbitCamera3D(screenWindow, sceneDelta);
		}
	}
	#else
	if (mouseInScene){
		const float wheel = GetMouseWheelMove();
		if (wheel != 0){
			const Vector2 mouseWorldPos =
				GetScreenToWorld2D(sceneMouse, screenWindow->cam);
			screenWindow->cam.offset = sceneMouse;
			screenWindow->cam.target = mouseWorldPos;
			const double speed = 0.0625*100.0/Ly;
			screenWindow->cam.zoom += (float)(wheel*speed*
				screenWindow->cam.zoom*screenWindow->cam.zoom);
			screenWindow->cam.zoom = Clamp(screenWindow->cam.zoom, 0.2f, 10.0f);
		}
		if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)){
			Vector2 delta = GetMouseDelta();
			delta.x *= (float)sceneTargetWidth/guiViewport.width;
			delta.y *= (float)sceneTargetHeight/guiViewport.height;
			const float scale = -1.0f/screenWindow->cam.zoom;
			screenWindow->cam.target.x += delta.x*scale;
			screenWindow->cam.target.y += delta.y*scale;
		}
	}
	#endif

	if (activeParticleDrag.active){
		Vector2 dragMouse = {0};
		if (gui_imgui_simulation_faulted()){
			finishParticleDrag(screenState);
		}
		else if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) &&
			screenPointToScene(GetMousePosition(), &dragMouse)){
			updateParticleDrag(screenWindow, dragMouse);
		}
		if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)){
			finishParticleDrag(screenState);
		}
	}
	else if (pendingParticlePress.pending){
		const Vector2 currentScreen = GetMousePosition();
		const bool moved = pendingParticlePressMoved(currentScreen);
		if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)){
			if (moved){
				startPendingParticleDrag(screenWindow, screenState,
					currentScreen, true);
			}
			else{
				clearPendingParticlePress();
			}
		}
		else if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && moved){
			startPendingParticleDrag(screenWindow, screenState,
				currentScreen, false);
		}
	}
	else if (mouseInScene && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)){
		const int result = getWhatsUnderClick(screenWindow, screenState);
		screenState->selected = (result == 1);
		screenState->leftClicked = false;
		if (result != 1){
			screenState->particleUnderClick = NULL;
		}
		else if (!gui_imgui_simulation_faulted()){
			pendingParticlePress = (particlePress){
				.pending = true,
				.index = screenState->particleUnderClick->num,
				.screenPoint = GetMousePosition(),
			};
		}
	}

}


static void invalidateParticleColorCache(void){
	colorCache.valid = false;
	#if !THREE_D
	maximumPredictedParticleRadiusCount = -1;
	#endif
}

static particlePrediction graphicsParticlePrediction(const particle* p){
	particlePrediction prediction;
	if (freeFly == freeFlyGrow
		#if POLYGON_PARTICLES && !THREE_D
		|| polygon_particles_growing()
		#endif
	){
		prediction = predictFreeFlyGrow(p, t);
	}
	else{
		prediction = predictFreeFlyNormal(p, t);
	}
	if (!addWallx && !addCircularWall)
		PBCpostX(&prediction.x);
	if (!addWally && !addCircularWall)
		PBCpostY(&prediction.y);
	#if THREE_D
	if (!addCircularWall)
		PBCpostZ(&prediction.z);
	#endif
	return prediction;
}

#if !THREE_D
static double graphicsMaximumParticleRadius(void){
	if ((maximumPredictedParticleRadiusCount == N) &&
		(maximumPredictedParticleRadiusTime == t)){
		return maximumPredictedParticleRadius;
	}

	double maximumRadius = 0;
	for (int i = 0; i < N; i++){
		const particlePrediction predicted = graphicsParticlePrediction(particles + i);
		if (isfinite(predicted.rad) && (predicted.rad > maximumRadius)){
			maximumRadius = predicted.rad;
		}
	}
	maximumPredictedParticleRadius = maximumRadius;
	maximumPredictedParticleRadiusTime = t;
	maximumPredictedParticleRadiusCount = N;
	return maximumRadius;
}
#endif

static bool refreshRenderPredictions(void){
	if (renderPredictionCapacity != N){
		particlePrediction* resized = realloc(renderPredictions,
			(size_t)N*sizeof(particlePrediction));
		if (resized == NULL){
			TraceLog(LOG_WARNING, "Unable to allocate the particle render snapshot");
			return false;
		}
		renderPredictions = resized;
		renderPredictionCapacity = N;
	}
	#if THREE_D
	double maximumRadius3D = 0;
	#else
	double maximumRadius = 0;
	#endif
	for (int i = 0; i < N; i++){
		renderPredictions[i] = graphicsParticlePrediction(particles + i);
		#if THREE_D
		if (isfinite(renderPredictions[i].rad) &&
			(renderPredictions[i].rad > maximumRadius3D)){
			maximumRadius3D = renderPredictions[i].rad;
		}
		#else
		if (isfinite(renderPredictions[i].rad) &&
			(renderPredictions[i].rad > maximumRadius)){
			maximumRadius = renderPredictions[i].rad;
		}
		#endif
	}
	#if THREE_D
	maximumPredictedParticleRadius3D = maximumRadius3D;
	#else
	maximumPredictedParticleRadius = maximumRadius;
	maximumPredictedParticleRadiusTime = t;
	maximumPredictedParticleRadiusCount = N;
	#endif
	return true;
}

static particlePrediction cachedGraphicsParticlePrediction(const particle* p){
	if ((renderPredictions != NULL) && (renderPredictionCapacity == N) &&
		(p->num >= 0) && (p->num < N)){
		return renderPredictions[p->num];
	}
	return graphicsParticlePrediction(p);
}

static bool ensureColorFunctionArray(state* screenState){
	if ((screenState->colorFunctionArray != NULL) &&
		(screenState->colorFunctionCapacity == N)){
		return true;
	}
	double* resized = realloc(screenState->colorFunctionArray,
		(size_t)N*sizeof(double));
	if (resized == NULL){
		TraceLog(LOG_WARNING, "Unable to allocate the particle color cache");
		free(screenState->colorFunctionArray);
		screenState->colorFunctionArray = NULL;
		screenState->colorFunctionCapacity = 0;
		screenState->colorParam2 = 0;
		invalidateParticleColorCache();
		return false;
	}
	screenState->colorFunctionArray = resized;
	screenState->colorFunctionCapacity = N;
	invalidateParticleColorCache();
	return true;
}

static double validatedRenderRefreshSeconds(double requested, double fallback){
	if (!isfinite(requested) || (requested <= 0)){
		requested = fallback;
	}
	if (requested < RENDER_REFRESH_MIN_SECONDS){
		return RENDER_REFRESH_MIN_SECONDS;
	}
	if (requested > RENDER_REFRESH_MAX_SECONDS){
		return RENDER_REFRESH_MAX_SECONDS;
	}
	return requested;
}

static bool refreshParticleColorCache(state* screenState){
	if ((screenState->colorParam2 == 0) || !ensureColorFunctionArray(screenState)){
		return false;
	}

	const double now = GetTime();
	const double refreshInterval = validatedRenderRefreshSeconds(
		guiLive.color_refresh_seconds, COLOR_DEFAULT_REFRESH_SECONDS);
	const bool configurationChanged = !colorCache.valid ||
		(colorCache.mode != screenState->colorParam2) ||
		(colorCache.particleCount != N);
	const bool simulationChanged = (colorCache.simulationTime != t);
	const bool refreshDue = simulationChanged &&
		(!screenState->running || (now - colorCache.wallTime >= refreshInterval));
	if (!configurationChanged && !refreshDue){
		return true;
	}

	double min = INFINITY;
	double max = -INFINITY;
	for (int i = 0; i < N; i++){
		double value = screenState->colorFunction(particles + i);
		if (!isfinite(value)) value = 0;
		screenState->colorFunctionArray[i] = value;
		if (value < min) min = value;
		if (value > max) max = value;
	}
	if (!isfinite(min) || !isfinite(max)){
		min = 0;
		max = 1;
	}
	else if (min == max){
		max = min + 1;
	}
	colorCache = (particleColorCache){
		.valid = true,
		.mode = screenState->colorParam2,
		.particleCount = N,
		.simulationTime = t,
		.wallTime = now,
		.min = min,
		.max = max,
	};
	return true;
}

#if !THREE_D
/* A signed-distance circle keeps the dense 2D path to one mesh draw without
 * shipping or magnifying a bitmap image.  fwidth provides approximately one
 * physical pixel of antialiasing at every camera zoom and framebuffer DPI. */
static const char* proceduralCircleFragmentShader2D =
	#if defined(PLATFORM_WEB)
	"#version 300 es\n"
	"precision highp float;\n"
	#else
	"#version 330\n"
	#endif
	"in vec2 fragTexCoord;\n"
	"in vec4 fragColor;\n"
	"uniform vec4 colDiffuse;\n"
	"uniform float circleEdgeScale;\n"
	"out vec4 finalColor;\n"
	"void main(void)\n"
	"{\n"
	"    vec2 centered = 2.0*fragTexCoord - vec2(1.0);\n"
	"    float distanceSquared = dot(centered, centered);\n"
	"    float edgeWidth = max(circleEdgeScale*fwidth(distanceSquared), 0.002);\n"
	"    float coverage = 1.0 - smoothstep(1.0-edgeWidth, 1.0+edgeWidth,\n"
	"                                      distanceSquared);\n"
	"    if (coverage <= 0.0) discard;\n"
	"    vec4 color = fragColor*colDiffuse;\n"
	"    finalColor = vec4(color.rgb, color.a*coverage);\n"
	"}\n";

static bool loadProceduralCircleMaterial(void){
	Shader shader = LoadShaderFromMemory(NULL,
		proceduralCircleFragmentShader2D);
	/* A failed raylib shader load can yield its default shader.  The unique
	 * uniform distinguishes that fallback from the program requested above. */
	if ((shader.id == 0) || (shader.locs == NULL)){
		if (shader.id != 0) UnloadShader(shader);
		TraceLog(LOG_WARNING,
			"Procedural particle shader unavailable; using immediate circles");
		return false;
	}
	const int edgeScaleLocation = GetShaderLocation(shader, "circleEdgeScale");
	if (edgeScaleLocation < 0){
		/* LoadShaderFromMemory returns the shared default program on a compile
		 * failure, but with a newly allocated location table.  Do not ask
		 * UnloadShader() to unload that shared program; release only its table. */
		if (arrowMaterialLoaded && shader.id == arrowMaterial.shader.id)
			MemFree(shader.locs);
		else
			UnloadShader(shader);
		TraceLog(LOG_WARNING,
			"Procedural particle shader unavailable; using immediate circles");
		return false;
	}

	Material material = LoadMaterialDefault();
	if (material.maps == NULL){
		UnloadShader(shader);
		TraceLog(LOG_WARNING,
			"Procedural particle material unavailable; using immediate circles");
		return false;
	}
	material.shader = shader;
	const float edgeScale = 1.0f;
	SetShaderValue(shader, edgeScaleLocation, &edgeScale, SHADER_UNIFORM_FLOAT);
	proceduralCircleMaterial = material;
	proceduralCircleMaterialLoaded = true;
	return true;
}

static void unloadCircleBatch(circleBatch* batch){
	if (batch->loaded){
		UnloadMesh(batch->mesh);
	}
	*batch = (circleBatch){0};
}

static bool loadCircleBatch(circleBatch* batch, int capacity){
	if (capacity <= 0){
		return false;
	}

	const size_t vertexCount = (size_t)capacity*6;
	Mesh mesh = {0};
	mesh.vertexCount = (int)vertexCount;
	mesh.triangleCount = capacity*2;
	mesh.vertices = calloc(vertexCount*3, sizeof(float));
	mesh.texcoords = malloc(vertexCount*2*sizeof(float));
	mesh.colors = malloc(vertexCount*4*sizeof(unsigned char));
	if ((mesh.vertices == NULL) || (mesh.texcoords == NULL) || (mesh.colors == NULL)){
		free(mesh.vertices);
		free(mesh.texcoords);
		free(mesh.colors);
		TraceLog(LOG_ERROR, "Unable to allocate the particle render batch");
		return false;
	}

	static const float quadUV[12] = {
		0.0f, 0.0f,  0.0f, 1.0f,  1.0f, 1.0f,
		0.0f, 0.0f,  1.0f, 1.0f,  1.0f, 0.0f,
	};
	for (int circle = 0; circle < capacity; circle++){
		memcpy(mesh.texcoords + (size_t)circle*12, quadUV, sizeof(quadUV));
	}
	memset(mesh.colors, 255, vertexCount*4*sizeof(unsigned char));

	UploadMesh(&mesh, true);
	batch->mesh = mesh;
	batch->capacity = capacity;
	batch->count = 0;
	batch->loaded = true;
	return true;
}

static bool ensureCircleBatch(circleBatch* batch, int capacity){
	if (batch->loaded && (batch->capacity == capacity)){
		return true;
	}
	unloadCircleBatch(batch);
	return loadCircleBatch(batch, capacity);
}

static void setCircleBatchVertex(circleBatch* batch, int vertex,
	float x, float y, Color color){
	float* position = batch->mesh.vertices + (size_t)vertex*3;
	position[0] = x;
	position[1] = y;
	position[2] = 0.0f;

	unsigned char* vertexColor = batch->mesh.colors + (size_t)vertex*4;
	vertexColor[0] = color.r;
	vertexColor[1] = color.g;
	vertexColor[2] = color.b;
	vertexColor[3] = color.a;
}

static bool appendCircleQuad(circleBatch* batch, const Vector2 corners[4],
	Color color){
	if (!batch->loaded || (batch->count >= batch->capacity)){
		return false;
	}

	const int firstVertex = batch->count*6;
	static const int cornerOrder[6] = {0, 1, 2, 0, 2, 3};
	for (int i = 0; i < 6; i++){
		const Vector2 corner = corners[cornerOrder[i]];
		setCircleBatchVertex(batch, firstVertex + i, corner.x, corner.y, color);
	}
	batch->count++;
	return true;
}

static void drawCircleBatch(circleBatch* batch){
	if (!batch->loaded || (batch->count == 0)){
		return;
	}

	const int vertexCount = batch->count*6;
	UpdateMeshBuffer(batch->mesh, 0, batch->mesh.vertices,
		vertexCount*3*(int)sizeof(float), 0);
	UpdateMeshBuffer(batch->mesh, 3, batch->mesh.colors,
		vertexCount*4*(int)sizeof(unsigned char), 0);

	Mesh visibleMesh = batch->mesh;
	visibleMesh.vertexCount = vertexCount;
	visibleMesh.triangleCount = batch->count*2;
	DrawMesh(visibleMesh, proceduralCircleMaterial, MatrixIdentity());
}

#if POLYGON_PARTICLES
static void unloadPolygonTriangleBatch(polygonTriangleBatch* batch){
	if (batch->meshLoaded){
		UnloadMesh(batch->mesh);
	}
	if (batch->materialLoaded){
		UnloadMaterial(batch->material);
	}
	*batch = (polygonTriangleBatch){0};
}

static bool loadPolygonTriangleBatch(polygonTriangleBatch* batch,
	int capacity, int shapeVertexCount){
	if (capacity <= 0 || shapeVertexCount < 3 ||
		shapeVertexCount > SIM_POLYGON_MAX_VERTICES){
		return false;
	}

	const int verticesPerParticle = 3*(shapeVertexCount - 2);
	const size_t maximumVertexCount =
		(size_t)INT_MAX/(3*sizeof(float));
	if ((size_t)capacity >
		maximumVertexCount/(size_t)verticesPerParticle){
		TraceLog(LOG_ERROR, "Polygon render batch is too large");
		return false;
	}
	const size_t vertexCount = (size_t)capacity*(size_t)verticesPerParticle;

	Mesh mesh = {0};
	mesh.vertexCount = (int)vertexCount;
	mesh.triangleCount = mesh.vertexCount/3;
	mesh.vertices = calloc(vertexCount*3, sizeof(float));
	mesh.colors = malloc(vertexCount*4*sizeof(unsigned char));
	if (mesh.vertices == NULL || mesh.colors == NULL){
		free(mesh.vertices);
		free(mesh.colors);
		TraceLog(LOG_ERROR, "Unable to allocate the polygon particle render batch");
		return false;
	}
	memset(mesh.colors, 255, vertexCount*4*sizeof(unsigned char));
	UploadMesh(&mesh, true);

	Material material = LoadMaterialDefault();
	if (material.maps == NULL){
		UnloadMesh(mesh);
		TraceLog(LOG_ERROR, "Unable to allocate the polygon particle material");
		return false;
	}

	batch->mesh = mesh;
	batch->material = material;
	batch->capacity = capacity;
	batch->shapeVertexCount = shapeVertexCount;
	batch->verticesPerParticle = verticesPerParticle;
	batch->count = 0;
	batch->meshLoaded = true;
	batch->materialLoaded = true;
	return true;
}

static bool ensurePolygonTriangleBatch(polygonTriangleBatch* batch,
	int capacity, int shapeVertexCount){
	if (batch->meshLoaded && batch->materialLoaded &&
		batch->capacity == capacity &&
		batch->shapeVertexCount == shapeVertexCount){
		return true;
	}
	unloadPolygonTriangleBatch(batch);
	return loadPolygonTriangleBatch(batch, capacity, shapeVertexCount);
}

static void setPolygonTriangleBatchVertex(polygonTriangleBatch* batch,
	int vertex, Vector2 point, Color color){
	float* position = batch->mesh.vertices + (size_t)vertex*3;
	position[0] = point.x;
	position[1] = point.y;
	position[2] = 0.0f;

	unsigned char* vertexColor = batch->mesh.colors + (size_t)vertex*4;
	vertexColor[0] = color.r;
	vertexColor[1] = color.g;
	vertexColor[2] = color.b;
	vertexColor[3] = color.a;
}

static Vector2 polygonParticleScreenVertex(Vector2 canonical,
	const particlePrediction* predicted, double cosine, double sine,
	double factor, double scale){
	const double localX = scale*predicted->rad*canonical.x;
	const double localY = scale*predicted->rad*canonical.y;
	const double worldX = predicted->x + cosine*localX - sine*localY;
	const double worldY = predicted->y + sine*localX + cosine*localY;
	return (Vector2){
		(float)(worldX*factor),
		(float)((Ly - worldY)*factor),
	};
}

static bool appendPolygonParticle(polygonTriangleBatch* batch,
	const Vector2* canonical, int shapeVertexCount,
	const particlePrediction* predicted, double angle, double factor,
	Color color){
	if (!batch->meshLoaded || !batch->materialLoaded || canonical == NULL ||
		predicted == NULL || batch->count >= batch->capacity ||
		shapeVertexCount != batch->shapeVertexCount){
		return false;
	}

	const double cosine = cos(angle);
	const double sine = sin(angle);
	const Vector2 first = polygonParticleScreenVertex(canonical[0], predicted,
		cosine, sine, factor, 1.0);
	int vertex = batch->count*batch->verticesPerParticle;
	for (int triangle = 1; triangle < shapeVertexCount - 1; triangle++){
		const Vector2 second = polygonParticleScreenVertex(canonical[triangle],
			predicted, cosine, sine, factor, 1.0);
		const Vector2 third = polygonParticleScreenVertex(canonical[triangle + 1],
			predicted, cosine, sine, factor, 1.0);
		setPolygonTriangleBatchVertex(batch, vertex++, first, color);
		setPolygonTriangleBatchVertex(batch, vertex++, second, color);
		setPolygonTriangleBatchVertex(batch, vertex++, third, color);
	}
	batch->count++;
	return true;
}

static void drawPolygonTriangleBatch(polygonTriangleBatch* batch){
	if (!batch->meshLoaded || !batch->materialLoaded || batch->count == 0){
		return;
	}
	const int vertexCount = batch->count*batch->verticesPerParticle;
	UpdateMeshBuffer(batch->mesh, 0, batch->mesh.vertices,
		vertexCount*3*(int)sizeof(float), 0);
	UpdateMeshBuffer(batch->mesh, 3, batch->mesh.colors,
		vertexCount*4*(int)sizeof(unsigned char), 0);
	Mesh visibleMesh = batch->mesh;
	visibleMesh.vertexCount = vertexCount;
	visibleMesh.triangleCount = vertexCount/3;
	DrawMesh(visibleMesh, batch->material, MatrixIdentity());
}
#endif

typedef struct velocityArrowShape2D velocityArrowShape2D;
struct velocityArrowShape2D{
	Vector2 start;
	Vector2 shaftEnd;
	Vector2 tip;
	Vector2 normal;
	float outerShaftHalfWidth;
	float outerHeadHalfWidth;
	float innerShaftHalfWidth;
	float innerHeadHalfWidth;
};

static void unloadVelocityArrowBatch(velocityArrowBatch* batch){
	if (batch->loaded){
		UnloadMesh(batch->mesh);
	}
	*batch = (velocityArrowBatch){0};
}

static bool loadVelocityArrowBatch(velocityArrowBatch* batch, int capacity){
	if (capacity <= 0){
		return false;
	}
	const size_t vertexCount = (size_t)capacity*VELOCITY_ARROW_VERTICES;
	Mesh mesh = {0};
	mesh.vertexCount = (int)vertexCount;
	mesh.triangleCount = (int)(vertexCount/3);
	mesh.vertices = calloc(vertexCount*3, sizeof(float));
	mesh.colors = malloc(vertexCount*4*sizeof(unsigned char));
	if (mesh.vertices == NULL || mesh.colors == NULL){
		free(mesh.vertices);
		free(mesh.colors);
		TraceLog(LOG_ERROR, "Unable to allocate the velocity-arrow render batch");
		return false;
	}
	memset(mesh.colors, 255, vertexCount*4*sizeof(unsigned char));
	UploadMesh(&mesh, true);
	batch->mesh = mesh;
	batch->capacity = capacity;
	batch->count = 0;
	batch->loaded = true;
	return true;
}

static bool ensureVelocityArrowBatch(velocityArrowBatch* batch, int capacity){
	if (batch->loaded && batch->capacity == capacity){
		return true;
	}
	unloadVelocityArrowBatch(batch);
	return loadVelocityArrowBatch(batch, capacity);
}

static void setVelocityArrowVertex(velocityArrowBatch* batch, int vertex,
	Vector2 point, Color color){
	float* position = batch->mesh.vertices + (size_t)vertex*3;
	position[0] = point.x;
	position[1] = point.y;
	position[2] = 0;
	unsigned char* vertexColor = batch->mesh.colors + (size_t)vertex*4;
	vertexColor[0] = color.r;
	vertexColor[1] = color.g;
	vertexColor[2] = color.b;
	vertexColor[3] = color.a;
}

static void appendVelocityArrowTriangle(velocityArrowBatch* batch, int* vertex,
	Vector2 a, Vector2 b, Vector2 c, Color color){
	setVelocityArrowVertex(batch, (*vertex)++, a, color);
	setVelocityArrowVertex(batch, (*vertex)++, b, color);
	setVelocityArrowVertex(batch, (*vertex)++, c, color);
}

static bool velocityArrowShape(double vx, double vy, Vector2 center,
	float radius, float zoom, velocityArrowShape2D* shape){
	const double speedSquared = vx*vx + vy*vy;
	if (shape == NULL || !isfinite(speedSquared) || speedSquared <= 1e-20){
		return false;
	}
	const double speed = sqrt(speedSquared);
	const Vector2 direction = {
		(float)(vx/speed),
		(float)(-vy/speed),
	};
	const float safeZoom = fmaxf(zoom, 0.01f);
	const float lengthScale = Clamp(1.55f + 0.70f*(float)log1p(speed),
		1.55f, 4.75f);
	const float length = fmaxf(radius*lengthScale, 8.0f/safeZoom);
	const float headLength = fminf(0.36f*length,
		fmaxf(1.05f*radius, 4.0f/safeZoom));
	const Vector2 normal = {-direction.y, direction.x};
	shape->start = Vector2Subtract(center,
		Vector2Scale(direction, 0.12f*radius));
	shape->tip = Vector2Add(center, Vector2Scale(direction, length));
	shape->shaftEnd = Vector2Subtract(shape->tip,
		Vector2Scale(direction, headLength));
	shape->normal = normal;
	shape->outerShaftHalfWidth = fmaxf(0.10f*radius, 0.80f/safeZoom);
	shape->outerHeadHalfWidth = fmaxf(0.34f*radius, 2.40f/safeZoom);
	shape->innerShaftHalfWidth = 0.56f*shape->outerShaftHalfWidth;
	shape->innerHeadHalfWidth = 0.68f*shape->outerHeadHalfWidth;
	return true;
}

static void appendVelocityArrowLayer(velocityArrowBatch* batch, int* vertex,
	const velocityArrowShape2D* shape, float shaftHalfWidth,
	float headHalfWidth, Color color){
	const Vector2 shaftOffset = Vector2Scale(shape->normal, shaftHalfWidth);
	const Vector2 headOffset = Vector2Scale(shape->normal, headHalfWidth);
	const Vector2 startMinus = Vector2Subtract(shape->start, shaftOffset);
	const Vector2 startPlus = Vector2Add(shape->start, shaftOffset);
	const Vector2 endMinus = Vector2Subtract(shape->shaftEnd, shaftOffset);
	const Vector2 endPlus = Vector2Add(shape->shaftEnd, shaftOffset);
	appendVelocityArrowTriangle(batch, vertex,
		startMinus, startPlus, endPlus, color);
	appendVelocityArrowTriangle(batch, vertex,
		startMinus, endPlus, endMinus, color);
	appendVelocityArrowTriangle(batch, vertex,
		Vector2Subtract(shape->shaftEnd, headOffset),
		Vector2Add(shape->shaftEnd, headOffset), shape->tip, color);
}

static bool appendVelocityArrow2D(velocityArrowBatch* batch, double vx,
	double vy, Vector2 center, float radius, float zoom){
	if (!batch->loaded || batch->count >= batch->capacity){
		return false;
	}
	velocityArrowShape2D shape;
	if (!velocityArrowShape(vx, vy, center, radius, zoom, &shape)){
		return false;
	}
	int vertex = batch->count*VELOCITY_ARROW_VERTICES;
	const Color outline = {236, 239, 244, 235};
	const Color fill = {255, 183, 45, 255};
	appendVelocityArrowLayer(batch, &vertex, &shape,
		shape.outerShaftHalfWidth, shape.outerHeadHalfWidth, outline);
	appendVelocityArrowLayer(batch, &vertex, &shape,
		shape.innerShaftHalfWidth, shape.innerHeadHalfWidth, fill);
	batch->count++;
	return true;
}

static void drawVelocityArrow2DImmediate(double vx, double vy, Vector2 center,
	float radius, float zoom){
	velocityArrowShape2D shape;
	if (!velocityArrowShape(vx, vy, center, radius, zoom, &shape)) return;
	const Color outline = {236, 239, 244, 235};
	const Color fill = {255, 183, 45, 255};
	const Vector2 outerHead = Vector2Scale(shape.normal,
		shape.outerHeadHalfWidth);
	const Vector2 innerHead = Vector2Scale(shape.normal,
		shape.innerHeadHalfWidth);
	DrawLineEx(shape.start, shape.shaftEnd,
		2.0f*shape.outerShaftHalfWidth, outline);
	DrawTriangle(Vector2Subtract(shape.shaftEnd, outerHead),
		Vector2Add(shape.shaftEnd, outerHead), shape.tip, outline);
	DrawLineEx(shape.start, shape.shaftEnd,
		2.0f*shape.innerShaftHalfWidth, fill);
	DrawTriangle(Vector2Subtract(shape.shaftEnd, innerHead),
		Vector2Add(shape.shaftEnd, innerHead), shape.tip, fill);
}

static void drawVelocityArrowBatch(velocityArrowBatch* batch){
	if (!batch->loaded || batch->count == 0) return;
	const int vertexCount = batch->count*VELOCITY_ARROW_VERTICES;
	UpdateMeshBuffer(batch->mesh, 0, batch->mesh.vertices,
		vertexCount*3*(int)sizeof(float), 0);
	UpdateMeshBuffer(batch->mesh, 3, batch->mesh.colors,
		vertexCount*4*(int)sizeof(unsigned char), 0);
	Mesh visibleMesh = batch->mesh;
	visibleMesh.vertexCount = vertexCount;
	visibleMesh.triangleCount = vertexCount/3;
	DrawMesh(visibleMesh, arrowMaterial, MatrixIdentity());
}
#endif

#if THREE_D
static const char* sphereVertexShader3D =
	#if defined(PLATFORM_WEB)
	"#version 300 es\n"
	"precision highp float;\n"
	#else
	"#version 330\n"
	#endif
	"layout(location = 0) in vec3 vertexPosition;\n"
	"layout(location = 2) in vec3 vertexNormal;\n"
	"layout(location = 9) in mat4 instanceTransform;\n"
	"uniform mat4 mvp;\n"
	"out vec3 sphereNormal;\n"
	"out vec3 sphereWorldPosition;\n"
	"out vec4 sphereColor;\n"
	"void main(void)\n"
	"{\n"
	"    mat3 scale = mat3(instanceTransform);\n"
	"    vec3 center = instanceTransform[3].xyz;\n"
	"    vec3 worldPosition = center + scale*vertexPosition;\n"
	"    sphereNormal = normalize(scale*vertexNormal);\n"
	"    sphereWorldPosition = worldPosition;\n"
	"    sphereColor = vec4(instanceTransform[0].w,\n"
	"                       instanceTransform[1].w,\n"
	"                       instanceTransform[2].w,\n"
	"                       instanceTransform[3].w);\n"
	"    gl_Position = mvp*vec4(worldPosition, 1.0);\n"
	"}\n";

static const char* sphereFragmentShader3D =
	#if defined(PLATFORM_WEB)
	"#version 300 es\n"
	"precision highp float;\n"
	#else
	"#version 330\n"
	#endif
	"in vec3 sphereNormal;\n"
	"in vec3 sphereWorldPosition;\n"
	"in vec4 sphereColor;\n"
	"uniform int lightingMode;\n"
	"uniform vec3 cameraPosition;\n"
	"uniform vec3 backgroundColor;\n"
	"uniform float sceneRadius;\n"
	"out vec4 finalColor;\n"
	"void main(void)\n"
	"{\n"
	"    vec3 normal = normalize(sphereNormal);\n"
	"    vec3 keyLight = normalize(vec3(-0.35, 0.80, 0.48));\n"
	"    vec3 fillLight = normalize(vec3(0.55, -0.20, -0.80));\n"
	"    float key = max(dot(normal, keyLight), 0.0);\n"
	"    float fill = max(dot(normal, fillLight), 0.0);\n"
	"    if (lightingMode == 0)\n"
	"    {\n"
	"        float brightness = min(1.0, 0.42 + 0.50*key + 0.10*fill);\n"
	"        finalColor = vec4(sphereColor.rgb*brightness, sphereColor.a);\n"
	"        return;\n"
	"    }\n"
	"    vec3 viewVector = cameraPosition - sphereWorldPosition;\n"
	"    vec3 viewDirection = viewVector*\n"
	"        inversesqrt(max(dot(viewVector, viewVector), 1.0e-8));\n"
	"    float hemisphere = 0.5 + 0.5*normal.y;\n"
	"    vec3 ambient = mix(vec3(0.16, 0.18, 0.23),\n"
	"                       vec3(0.32, 0.35, 0.39), hemisphere);\n"
	"    float wrappedKey = max((dot(normal, keyLight) + 0.16)/1.16, 0.0);\n"
	"    vec3 halfVector = keyLight + viewDirection;\n"
	"    vec3 halfDirection = halfVector*\n"
	"        inversesqrt(max(dot(halfVector, halfVector), 1.0e-8));\n"
	"    float specular = pow(max(dot(normal, halfDirection), 0.0), 32.0);\n"
	"    specular *= smoothstep(0.0, 0.18, key);\n"
	"    vec3 litColor = sphereColor.rgb*\n"
	"        (ambient + vec3(0.56*wrappedKey + 0.11*fill));\n"
	"    litColor += vec3(1.00, 0.93, 0.82)*(0.20*specular);\n"
	"    if (lightingMode >= 2)\n"
	"    {\n"
	"        float viewFacing = clamp(dot(normal, viewDirection), 0.0, 1.0);\n"
	"        float edgeSeparation = smoothstep(0.02, 0.62, viewFacing);\n"
	"        litColor *= mix(0.72, 1.0, edgeSeparation);\n"
	"        float rim = pow(1.0 - viewFacing, 3.0);\n"
	"        litColor += vec3(0.10, 0.15, 0.22)*rim*(0.35 + 0.65*key);\n"
	"        float relativeDepth = sqrt(max(dot(viewVector, viewVector), 0.0))/\n"
	"            max(sceneRadius, 0.001);\n"
	"        float haze = 0.12*smoothstep(2.25, 6.0, relativeDepth);\n"
	"        litColor = mix(litColor, backgroundColor, haze);\n"
	"    }\n"
	"    finalColor = vec4(clamp(litColor, 0.0, 1.0), sphereColor.a);\n"
	"}\n";

static void sphereResolution3D(int particleCount, int* rings, int* slices){
	switch (guiLive.sphere_quality){
		case SIM_SPHERE_QUALITY_LOW:
			*rings = 4; *slices = 6; return;
		case SIM_SPHERE_QUALITY_MEDIUM:
			*rings = 8; *slices = 12; return;
		case SIM_SPHERE_QUALITY_HIGH:
			*rings = 16; *slices = 24; return;
		case SIM_SPHERE_QUALITY_ULTRA:
			*rings = 24; *slices = 36; return;
		case SIM_SPHERE_QUALITY_AUTO:
		default:
			break;
	}
	/* Keep the total triangle budget roughly constant as systems get larger. */
	if (particleCount <= 2500){
		*rings = 16;
		*slices = 24;
	}
	else if (particleCount <= 5000){
		*rings = 12;
		*slices = 16;
	}
	else if (particleCount <= 10000){
		*rings = 8;
		*slices = 12;
	}
	else if (particleCount <= 20000){
		*rings = 6;
		*slices = 8;
	}
	else{
		*rings = 4;
		*slices = 6;
	}
}

static void unloadSphereBatch(sphereBatch* batch){
	if (batch->materialLoaded){
		/* raylib's UnloadMaterial also releases its non-default shader. */
		UnloadMaterial(batch->material);
	}
	else if (batch->shaderLoaded){
		UnloadShader(batch->shader);
	}
	if (batch->meshLoaded){
		UnloadMesh(batch->mesh);
	}
	free(batch->transforms);
	*batch = (sphereBatch){0};
}

static bool loadSphereBatch(sphereBatch* batch, int capacity,
	int rings, int slices){
	if (capacity <= 0){
		return false;
	}

	sphereBatch loaded = {
		.capacity = capacity,
		.rings = rings,
		.slices = slices,
	};
	loaded.transforms = malloc((size_t)capacity*sizeof(Matrix));
	if (loaded.transforms == NULL){
		TraceLog(LOG_ERROR, "Unable to allocate the 3D particle batch");
		return false;
	}

	loaded.mesh = GenMeshSphere(1.0f, rings, slices);
	loaded.meshLoaded = (loaded.mesh.vertices != NULL);
	if (!loaded.meshLoaded){
		TraceLog(LOG_ERROR, "Unable to create the 3D sphere mesh");
		unloadSphereBatch(&loaded);
		return false;
	}
	if (loaded.mesh.vaoId == 0){
		UploadMesh(&loaded.mesh, false);
	}
	if (loaded.mesh.vaoId == 0){
		TraceLog(LOG_ERROR, "Unable to upload the 3D sphere mesh");
		unloadSphereBatch(&loaded);
		return false;
	}

	loaded.shader = LoadShaderFromMemory(sphereVertexShader3D,
		sphereFragmentShader3D);
	loaded.shaderLoaded = (loaded.shader.id != 0);
	if (!loaded.shaderLoaded || (loaded.shader.locs == NULL)){
		TraceLog(LOG_ERROR, "Unable to load the instanced sphere shader");
		unloadSphereBatch(&loaded);
		return false;
	}
	loaded.shader.locs[SHADER_LOC_MATRIX_MVP] = GetShaderLocation(
		loaded.shader, "mvp");
	loaded.shader.locs[SHADER_LOC_MATRIX_MODEL] = GetShaderLocationAttrib(
		loaded.shader, "instanceTransform");
	loaded.lightingModeLocation = GetShaderLocation(loaded.shader,
		"lightingMode");
	loaded.cameraPositionLocation = GetShaderLocation(loaded.shader,
		"cameraPosition");
	loaded.backgroundColorLocation = GetShaderLocation(loaded.shader,
		"backgroundColor");
	loaded.sceneRadiusLocation = GetShaderLocation(loaded.shader,
		"sceneRadius");
	if ((loaded.shader.locs[SHADER_LOC_MATRIX_MVP] < 0) ||
		(loaded.shader.locs[SHADER_LOC_MATRIX_MODEL] < 0) ||
		(loaded.lightingModeLocation < 0) ||
		(loaded.cameraPositionLocation < 0) ||
		(loaded.backgroundColorLocation < 0) ||
		(loaded.sceneRadiusLocation < 0)){
		TraceLog(LOG_ERROR, "The instanced sphere shader is missing required inputs");
		unloadSphereBatch(&loaded);
		return false;
	}

	loaded.material = LoadMaterialDefault();
	loaded.materialLoaded = (loaded.material.maps != NULL);
	if (!loaded.materialLoaded){
		TraceLog(LOG_ERROR, "Unable to create the 3D sphere material");
		unloadSphereBatch(&loaded);
		return false;
	}
	loaded.material.shader = loaded.shader;
	loaded.count = 0;
	loaded.loaded = true;
	*batch = loaded;
	return true;
}

static bool ensureSphereBatch(sphereBatch* batch, int capacity){
	int rings = 0;
	int slices = 0;
	sphereResolution3D(capacity, &rings, &slices);
	if (batch->loaded && (batch->capacity == capacity) &&
		(batch->rings == rings) && (batch->slices == slices)){
		return true;
	}
	unloadSphereBatch(batch);
	return loadSphereBatch(batch, capacity, rings, slices);
}

static bool appendSphereBatch(sphereBatch* batch, Vector3 center,
	float radius, Color color){
	if (!batch->loaded || (batch->count >= batch->capacity)){
		return false;
	}
	const float inverse255 = 1.0f/255.0f;
	batch->transforms[batch->count] = (Matrix){
		.m0 = radius,
		.m5 = radius,
		.m10 = radius,
		.m12 = center.x,
		.m13 = center.y,
		.m14 = center.z,
		/* DrawMeshInstanced uploads matrices verbatim. The affine bottom row is
		 * unused by the shader, so it carries each particle's RGBA color. */
		.m3 = color.r*inverse255,
		.m7 = color.g*inverse255,
		.m11 = color.b*inverse255,
		.m15 = color.a*inverse255,
	};
	batch->count++;
	return true;
}

static void drawSphereBatch(sphereBatch* batch, const Camera3D* camera){
	if (!batch->loaded || (batch->count == 0) || (camera == NULL)){
		return;
	}
	int lightingMode = (int)guiLive.lighting_mode;
	if (lightingMode < SIM_LIGHTING_SIMPLE ||
		lightingMode > SIM_LIGHTING_DEPTH_ENHANCED){
		lightingMode = SIM_LIGHTING_STUDIO;
	}
	const Color background = colorFromRgba(guiLive.background_rgba);
	const Vector3 backgroundColor = {
		background.r/255.0f, background.g/255.0f, background.b/255.0f,
	};
	const float sceneRadius = boxRadius3D();
	SetShaderValue(batch->shader, batch->lightingModeLocation,
		&lightingMode, SHADER_UNIFORM_INT);
	SetShaderValue(batch->shader, batch->cameraPositionLocation,
		&camera->position, SHADER_UNIFORM_VEC3);
	SetShaderValue(batch->shader, batch->backgroundColorLocation,
		&backgroundColor, SHADER_UNIFORM_VEC3);
	SetShaderValue(batch->shader, batch->sceneRadiusLocation,
		&sceneRadius, SHADER_UNIFORM_FLOAT);
	DrawMeshInstanced(batch->mesh, batch->material,
		batch->transforms, batch->count);
}
#endif

#if !THREE_D
static bool particleIntersectsViewport(float left, float top, float right, float bottom, Rectangle viewport){
	return !((right < viewport.x) || (left > viewport.x + viewport.width) ||
		(bottom < viewport.y) || (top > viewport.y + viewport.height));
}
#endif

static int structureFactorPopulation(void){
	int population = N - (int)(N*fractionSmallN);
	if (population < 0) population = 0;
	if (population > N) population = N;
	return population;
}

static int validatedStructureFactorSampleCap(void){
	int cap = guiLive.structure_factor_sample_cap;
	if (cap <= 0){
		cap = STRUCTURE_FACTOR_DEFAULT_SAMPLE_CAP;
	}
	if (cap > STRUCTURE_FACTOR_HARD_SAMPLE_CAP){
		cap = STRUCTURE_FACTOR_HARD_SAMPLE_CAP;
	}
	return cap;
}

static int desiredStructureFactorSampleCount(void){
	const int population = structureFactorPopulation();
	const int sampleCap = validatedStructureFactorSampleCap();
	return (population < sampleCap) ? population : sampleCap;
}

static bool allocateStructureFactorPositions(void){
	const int desiredCount = desiredStructureFactorSampleCount();
	if (desiredCount == 0){
		structSampleCount = 0;
		return false;
	}
	if ((positions != NULL) && (desiredCount <= structSampleCapacity)){
		structSampleCount = desiredCount;
		structRejectedSampleCount = -1;
		return true;
	}
	if (desiredCount == structRejectedSampleCount){
		return false;
	}
	position* resized = realloc(positions, (size_t)desiredCount*sizeof(position));
	if (resized == NULL){
		structRejectedSampleCount = desiredCount;
		TraceLog(LOG_ERROR,
			"Unable to grow the structure-factor sample buffer to %d particles",
			desiredCount);
		return false;
	}
	positions = resized;
	structSampleCapacity = desiredCount;
	structSampleCount = desiredCount;
	structRejectedSampleCount = -1;
	return true;
}

static void captureStructureFactorPositions(void){
	const int population = structureFactorPopulation();
	if ((positions == NULL) || (structSampleCount <= 0) || (population <= 0)){
		return;
	}
	for (int sample = 0; sample < structSampleCount; sample++){
		const int particleIndex = (int)(((long long)sample*population)/structSampleCount);
		const particlePrediction predicted = graphicsParticlePrediction(particles + particleIndex);
		positions[sample].x = predicted.x;
		positions[sample].y = predicted.y;
	}
}

static bool structureFactorJobIsRunning(void){
	#if defined(PLATFORM_WEB)
	return structJobRunning;
	#else
	pthread_mutex_lock(&structMutex);
	const bool running = structJobRunning;
	pthread_mutex_unlock(&structMutex);
	return running;
	#endif
}

static void releaseStructureFactorPositions(void){
	#if defined(PLATFORM_WEB)
	structJobCancelRequested = true;
	#else
	pthread_mutex_lock(&structMutex);
	if (structJobRunning){
		structJobCancelRequested = true;
	}
	pthread_mutex_unlock(&structMutex);
	#endif
	awaitStructFactor();
	free(positions);
	positions = NULL;
	structSampleCount = 0;
	structSampleCapacity = 0;
	structRejectedSampleCount = -1;
}

#if !THREE_D
#if POLYGON_PARTICLES
static bool loadPolygonCanonicalVertices(Vector2* vertices, int vertexCount){
	if (vertices == NULL || vertexCount < 3 ||
		vertexCount > SIM_POLYGON_MAX_VERTICES){
		return false;
	}
	for (int vertex = 0; vertex < vertexCount; vertex++){
		double x = 0;
		double y = 0;
		if (!polygon_particles_vertex(vertex, &x, &y) ||
			!isfinite(x) || !isfinite(y)){
			return false;
		}
		vertices[vertex] = (Vector2){(float)x, (float)y};
	}
	return true;
}

static void drawPolygonSelectionOutline(const particle* selectedParticle,
	const particlePrediction* selectedPrediction, const Vector2* canonical,
	int shapeVertexCount, double factor, float zoom){
	if (selectedParticle == NULL || selectedPrediction == NULL ||
		canonical == NULL || shapeVertexCount < 3){
		return;
	}
	double angle = polygon_particles_angle_at(selectedParticle->num, t);
	if (!isfinite(angle)) angle = 0;
	const double cosine = cos(angle);
	const double sine = sin(angle);
	Vector2 outline[SIM_POLYGON_MAX_VERTICES];
	for (int vertex = 0; vertex < shapeVertexCount; vertex++){
		outline[vertex] = polygonParticleScreenVertex(canonical[vertex],
			selectedPrediction, cosine, sine, factor, 1.09);
	}
	const Color selectionColor = {245, 247, 250, 255};
	const float safeZoom = fmaxf(zoom, 0.01f);
	const float thickness = fmaxf(1.6f/safeZoom,
		(float)(0.055*selectedPrediction->rad*factor));
	for (int vertex = 0; vertex < shapeVertexCount; vertex++){
		DrawLineEx(outline[vertex], outline[(vertex + 1)%shapeVertexCount],
			thickness, selectionColor);
	}
}

static void drawPolygonParticlesAndBox(double factor, state* screenState,
	float zoom, Rectangle viewport){
	particle* particleUnderClick = screenState->particleUnderClick;
	const bool predictionsReady = refreshRenderPredictions();
	const bool colorsReady = refreshParticleColorCache(screenState);
	const int shapeVertexCount = polygon_particles_vertex_count();
	Vector2 canonical[SIM_POLYGON_MAX_VERTICES] = {0};
	const bool shapeReady = loadPolygonCanonicalVertices(canonical,
		shapeVertexCount);
	const bool particleBatchReady = (N == 0) || (shapeReady &&
		ensurePolygonTriangleBatch(&particlePolygonBatch, N,
			shapeVertexCount));
	const bool arrowBatchReady = !screenState->arrow ||
		(arrowMaterialLoaded && ensureVelocityArrowBatch(&arrowBatch, N));
	particlePolygonBatch.count = 0;
	arrowBatch.count = 0;

	for (int i = N - 1; i >= 0; i--){
		const particlePrediction predicted = predictionsReady ?
			renderPredictions[i] : graphicsParticlePrediction(particles + i);
		if (colorsReady){
			screenState->particleColor = colorSelect(
				(screenState->colorFunctionArray[i] - colorCache.min)/
					(colorCache.max - colorCache.min), screenState->colorArray);
		}
		else{
			screenState->particleColor =
				screenState->particleArray[particles[i].type];
		}

		const float left = (float)((predicted.x - predicted.rad)*factor);
		const float top = (float)((Ly - predicted.y - predicted.rad)*factor);
		const float diameter = (float)(2*predicted.rad*factor);
		if (!particleIntersectsViewport(left, top, left + diameter,
			top + diameter, viewport)){
			continue;
		}

		if (particleBatchReady && shapeReady){
			double angle = polygon_particles_angle_at(i, t);
			if (!isfinite(angle)) angle = 0;
			appendPolygonParticle(&particlePolygonBatch, canonical,
				shapeVertexCount, &predicted, angle, factor,
				screenState->particleColor);
		}
		else if (shapeReady){
			/* Allocation/GPU failure is exceptional, but it must not make every
			 * particle disappear. Keep a correct immediate-mode fallback. */
			double angle = polygon_particles_angle_at(i, t);
			if (!isfinite(angle)) angle = 0;
			const double cosine = cos(angle);
			const double sine = sin(angle);
			Vector2 outline[SIM_POLYGON_MAX_VERTICES];
			for (int vertex = 0; vertex < shapeVertexCount; vertex++){
				outline[vertex] = polygonParticleScreenVertex(canonical[vertex],
					&predicted, cosine, sine, factor, 1.0);
			}
			DrawTriangleFan(outline, shapeVertexCount,
				screenState->particleColor);
		}

		if (screenState->arrow && arrowBatchReady){
			appendVelocityArrow2D(&arrowBatch, predicted.vx, predicted.vy,
				(Vector2){(float)(predicted.x*factor),
					(float)((Ly - predicted.y)*factor)},
				(float)(predicted.rad*factor), zoom);
		}
	}

	if (particleBatchReady && shapeReady){
		drawPolygonTriangleBatch(&particlePolygonBatch);
	}
	if (screenState->arrow && arrowBatchReady){
		drawVelocityArrowBatch(&arrowBatch);
	}
	else if (screenState->arrow){
		for (int i = N - 1; i >= 0; i--){
			const particlePrediction predicted = predictionsReady ?
				renderPredictions[i] : graphicsParticlePrediction(particles + i);
			const float radius = (float)(predicted.rad*factor);
			const float centerX = (float)(predicted.x*factor);
			const float centerY = (float)((Ly - predicted.y)*factor);
			if (!particleIntersectsViewport(centerX - radius, centerY - radius,
				centerX + radius, centerY + radius, viewport)) continue;
			drawVelocityArrow2DImmediate(predicted.vx, predicted.vy,
				(Vector2){centerX, centerY}, radius, zoom);
		}
	}
	if (shapeReady && guiLive.show_selection_outline &&
		screenState->selected && particleUnderClick != NULL){
		const particlePrediction selected =
			cachedGraphicsParticlePrediction(particleUnderClick);
		drawPolygonSelectionOutline(particleUnderClick, &selected, canonical,
			shapeVertexCount, factor, zoom);
	}
	if (guiLive.show_boundary && addWallx){
		const Color boundaryColor = colorFromRgba(guiLive.boundary_rgba);
		DrawLineEx((Vector2){0, 0}, (Vector2){0, Ly*factor},
			0.3*factor/zoom, boundaryColor);
		DrawLineEx((Vector2){Lx*factor, 0}, (Vector2){Lx*factor, Ly*factor},
			0.3*factor/zoom, boundaryColor);
	}
	if (guiLive.show_boundary && addWally){
		const Color boundaryColor = colorFromRgba(guiLive.boundary_rgba);
		DrawLineEx((Vector2){0, 0}, (Vector2){Lx*factor, 0},
			0.3*factor/zoom, boundaryColor);
		DrawLineEx((Vector2){0, Ly*factor}, (Vector2){Lx*factor, Ly*factor},
			0.3*factor/zoom, boundaryColor);
	}
	else if (guiLive.show_boundary && addCircularWall){
		DrawRing((Vector2){Lx/2*factor, Lx/2*factor}, Lx/2*factor,
			(Lx/2 + 1)*factor, 0, 360, 360,
			colorFromRgba(guiLive.boundary_rgba));
	}
}
#endif

void drawParticlesAndBox(double factor, state* screenState, float zoom, Rectangle viewport){
	#if POLYGON_PARTICLES
	if (polygon_particles_requested()){
		drawPolygonParticlesAndBox(factor, screenState, zoom, viewport);
		return;
	}
	#endif
	particle* particleUnderClick = screenState->particleUnderClick;
	const bool predictionsReady = refreshRenderPredictions();
	const bool colorsReady = refreshParticleColorCache(screenState);
	const bool particleBatchReady = proceduralCircleMaterialLoaded &&
		ensureCircleBatch(&particleCircleBatch, N);
	const bool arrowBatchReady = !screenState->arrow ||
		(arrowMaterialLoaded && ensureVelocityArrowBatch(&arrowBatch, N));
	particleCircleBatch.count = 0;
	arrowBatch.count = 0;

	for (int i = N - 1; i >= 0; i--){
		const particlePrediction predicted = predictionsReady ? renderPredictions[i] :
			graphicsParticlePrediction(particles + i);
		if (colorsReady){
			screenState->particleColor = colorSelect(
				(screenState->colorFunctionArray[i] - colorCache.min)/
				(colorCache.max - colorCache.min), screenState->colorArray);
		}
		else{
			screenState->particleColor = screenState->particleArray[particles[i].type];
		}

		const float left = (float)((predicted.x - predicted.rad)*factor);
		const float top = (float)((Ly - predicted.y - predicted.rad)*factor);
		const float diameter = (float)(2*predicted.rad*factor);
		if (!particleIntersectsViewport(left, top, left + diameter, top + diameter, viewport)){
			continue;
		}

		const Vector2 particleCorners[4] = {
			{left, top}, {left, top + diameter},
			{left + diameter, top + diameter}, {left + diameter, top},
		};
		if (!particleBatchReady || !appendCircleQuad(&particleCircleBatch,
			particleCorners, screenState->particleColor)){
			DrawCircleV((Vector2){left + 0.5f*diameter, top + 0.5f*diameter},
				0.5f*diameter, screenState->particleColor);
		}

		if (screenState->arrow && arrowBatchReady){
			appendVelocityArrow2D(&arrowBatch, predicted.vx, predicted.vy,
				(Vector2){(float)(predicted.x*factor),
					(float)((Ly - predicted.y)*factor)},
				(float)(predicted.rad*factor), zoom);
		}
	}

	if (particleBatchReady){
		drawCircleBatch(&particleCircleBatch);
	}
	if (screenState->arrow && arrowBatchReady){
		drawVelocityArrowBatch(&arrowBatch);
	}
	else if (screenState->arrow){
		for (int i = N - 1; i >= 0; i--){
			const particlePrediction predicted = predictionsReady ? renderPredictions[i] :
				graphicsParticlePrediction(particles + i);
			const float radius = (float)(predicted.rad*factor);
			const float centerX = (float)(predicted.x*factor);
			const float centerY = (float)((Ly - predicted.y)*factor);
			if (!particleIntersectsViewport(centerX - radius, centerY - radius,
				centerX + radius, centerY + radius, viewport)) continue;
			drawVelocityArrow2DImmediate(predicted.vx, predicted.vy,
				(Vector2){centerX, centerY}, radius, zoom);
		}
	}
	if (guiLive.show_selection_outline && screenState->selected &&
		particleUnderClick != NULL){
		const particlePrediction selected = cachedGraphicsParticlePrediction(particleUnderClick);
		const Color selectionColor = {245, 247, 250, 255};
		DrawRing((Vector2){selected.x*factor, (Ly - selected.y)*factor},
			1.03*selected.rad*factor, 1.13*selected.rad*factor,
			0, 360, 48, selectionColor);
	}
	if (guiLive.show_boundary && addWallx){
		const Color boundaryColor = colorFromRgba(guiLive.boundary_rgba);
		DrawLineEx((Vector2){0, 0}, (Vector2){0, Ly*factor} , 0.3*factor/zoom, boundaryColor);
		DrawLineEx((Vector2){Lx*factor, 0}, (Vector2){Lx*factor, Ly*factor} , 0.3*factor/zoom, boundaryColor);
	}
	if (guiLive.show_boundary && addWally){
		const Color boundaryColor = colorFromRgba(guiLive.boundary_rgba);
		DrawLineEx((Vector2){0, 0}, (Vector2){Lx*factor, 0} , 0.3*factor/zoom, boundaryColor);
		DrawLineEx((Vector2){0, Ly*factor}, (Vector2){Lx*factor, Ly*factor} , 0.3*factor/zoom, boundaryColor);
	}
	#if !THREE_D
	else if (guiLive.show_boundary && addCircularWall){
		DrawRing((Vector2){Lx/2*factor, Lx/2*factor}, Lx/2*factor, (Lx/2 + 1)*factor, 0, 360, 360, colorFromRgba(guiLive.boundary_rgba));
	}
	#endif
}
#endif

#if THREE_D
static void drawBoxBoundary3D(Color color){
	const Vector3 corners[8] = {
		{0, 0, 0}, {(float)Lx, 0, 0},
		{0, (float)Ly, 0}, {(float)Lx, (float)Ly, 0},
		{0, 0, (float)Lz}, {(float)Lx, 0, (float)Lz},
		{0, (float)Ly, (float)Lz}, {(float)Lx, (float)Ly, (float)Lz},
	};
	static const int edges[12][2] = {
		{0, 1}, {0, 2}, {1, 3}, {2, 3},
		{4, 5}, {4, 6}, {5, 7}, {6, 7},
		{0, 4}, {1, 5}, {2, 6}, {3, 7},
	};
	const float radius = fmaxf(0.002f*boxRadius3D(), 0.0005f);
	for (int edge = 0; edge < 12; edge++){
		DrawCylinderEx(corners[edges[edge][0]], corners[edges[edge][1]],
			radius, radius, 8, color);
	}
}

static void drawSphericalBoundary3D(Color color){
	const Vector3 center = {
		(float)(0.5*Lx), (float)(0.5*Ly), (float)(0.5*Lz),
	};
	DrawSphereWires(center, (float)(0.5*Lx), 24, 48, color);
}

static void drawWallGrid3D(int axis, float coordinate, Color color){
	const int divisions = 5;
	for (int i = 0; i <= divisions; i++){
		const float u = (float)i/(float)divisions;
		if (axis == 0){
			const float y = (float)(u*Ly);
			const float z = (float)(u*Lz);
			DrawLine3D((Vector3){coordinate, y, 0},
				(Vector3){coordinate, y, (float)Lz}, color);
			DrawLine3D((Vector3){coordinate, 0, z},
				(Vector3){coordinate, (float)Ly, z}, color);
		}
		else if (axis == 1){
			const float x = (float)(u*Lx);
			const float z = (float)(u*Lz);
			DrawLine3D((Vector3){x, coordinate, 0},
				(Vector3){x, coordinate, (float)Lz}, color);
			DrawLine3D((Vector3){0, coordinate, z},
				(Vector3){(float)Lx, coordinate, z}, color);
		}
		else{
			const float x = (float)(u*Lx);
			const float y = (float)(u*Ly);
			DrawLine3D((Vector3){x, 0, coordinate},
				(Vector3){x, (float)Ly, coordinate}, color);
			DrawLine3D((Vector3){0, y, coordinate},
				(Vector3){(float)Lx, y, coordinate}, color);
		}
	}
}

static bool drawVelocityArrow3D(const CameraFrustum3D* frustum,
	const particlePrediction* predicted){
	if (predicted == NULL || !isfinite(predicted->rad) || predicted->rad <= 0)
		return false;
	const double speedSquared = predicted->vx*predicted->vx +
		predicted->vy*predicted->vy + predicted->vz*predicted->vz;
	if (!isfinite(speedSquared) || speedSquared <= 1e-20) return false;
	const double speed = sqrt(speedSquared);
	const Vector3 direction = {
		(float)(predicted->vx/speed),
		(float)(predicted->vy/speed),
		(float)(predicted->vz/speed),
	};
	const float radius = (float)predicted->rad;
	const float lengthScale = Clamp(2.20f + 0.80f*(float)log1p(speed),
		2.20f, 5.50f);
	const float length = radius*lengthScale;
	const float headLength = fminf(0.34f*length, 1.15f*radius);
	const Vector3 center = {
		(float)predicted->x, (float)predicted->y, (float)predicted->z,
	};
	const Vector3 tip = Vector3Add(center, Vector3Scale(direction, length));
	if (!camera_frustum_3d_intersects_arrow(frustum, center, tip,
		0.30f*radius)){
		return false;
	}
	const Vector3 shaftEnd = Vector3Subtract(tip,
		Vector3Scale(direction, headLength));
	const Color shaftColor = {236, 239, 244, 245};
	const Color headColor = {255, 181, 30, 255};
	/* A real shaft and cone remain readable while orbiting, unlike the old
	   one-pixel line.  Low side counts keep the opt-in overlay inexpensive. */
	DrawCylinderEx(center, shaftEnd, 0.095f*radius, 0.095f*radius,
		6, shaftColor);
	DrawCylinderEx(shaftEnd, tip, 0.30f*radius, 0.0f,
		8, headColor);
	return true;
}

static void drawParticlesAndBox3D(const Camera3D* camera,
	const CameraFrustum3D* frustum, state* screenState,
	bool predictionsReady){
	const bool colorsReady = refreshParticleColorCache(screenState);
	const bool sphereBatchReady = ensureSphereBatch(&particleSphereBatch, N);
	particleSphereBatch.count = 0;
	guiMetrics.visible_particle_count = 0;
	guiMetrics.frustum_culled_particle_count = 0;
	const Color boundaryColor = colorFromRgba(guiLive.boundary_rgba);
	if (guiLive.show_boundary){
		if (addCircularWall){
			drawSphericalBoundary3D(boundaryColor);
		}
		else{
			drawBoxBoundary3D(boundaryColor);
		}
	}
	if (guiLive.show_boundary && addWallx){
		const Color wallColor = Fade(boundaryColor, 0.55f);
		drawWallGrid3D(0, 0.0f, wallColor);
		drawWallGrid3D(0, (float)Lx, wallColor);
	}
	if (guiLive.show_boundary && addWally){
		const Color wallColor = Fade(boundaryColor, 0.55f);
		drawWallGrid3D(1, 0.0f, wallColor);
		drawWallGrid3D(1, (float)Ly, wallColor);
	}
	if (guiLive.show_boundary && addWallz){
		const Color wallColor = Fade(boundaryColor, 0.55f);
		drawWallGrid3D(2, 0.0f, wallColor);
		drawWallGrid3D(2, (float)Lz, wallColor);
	}

	if (guiLive.show_axes){
		const float axisLength = (float)(0.12*fmin(Lx, fmin(Ly, Lz)));
		const Vector3 axisOrigin = addCircularWall ?
			(Vector3){(float)(0.5*Lx), (float)(0.5*Ly), (float)(0.5*Lz)} :
			(Vector3){0};
		DrawLine3D(axisOrigin,
			(Vector3){axisOrigin.x + axisLength, axisOrigin.y, axisOrigin.z}, RED);
		DrawLine3D(axisOrigin,
			(Vector3){axisOrigin.x, axisOrigin.y + axisLength, axisOrigin.z}, GREEN);
		DrawLine3D(axisOrigin,
			(Vector3){axisOrigin.x, axisOrigin.y, axisOrigin.z + axisLength}, BLUE);
	}

	for (int i = N - 1; i >= 0; i--){
		const particlePrediction predicted = predictionsReady ? renderPredictions[i] :
			graphicsParticlePrediction(particles + i);
		const Vector3 center = {
			(float)predicted.x, (float)predicted.y, (float)predicted.z,
		};
		if (!camera_frustum_3d_intersects_sphere(frustum, center,
			predicted.rad)){
			guiMetrics.frustum_culled_particle_count++;
			continue;
		}
		guiMetrics.visible_particle_count++;
		Color particleColor;
		if (colorsReady){
			particleColor = colorSelect(
				(screenState->colorFunctionArray[i] - colorCache.min)/
					(colorCache.max - colorCache.min), screenState->colorArray);
		}
		else{
			particleColor = screenState->particleArray[particles[i].type];
		}
		if (sphereBatchReady){
			appendSphereBatch(&particleSphereBatch, center,
				(float)predicted.rad, particleColor);
		}
		else{
			DrawSphereEx(center, (float)predicted.rad, 4, 8, particleColor);
		}
	}
	if (sphereBatchReady){
		drawSphereBatch(&particleSphereBatch, camera);
	}
	if (screenState->arrow){
		for (int i = N - 1; i >= 0; i--){
			const particlePrediction predicted = predictionsReady ? renderPredictions[i] :
				graphicsParticlePrediction(particles + i);
			drawVelocityArrow3D(frustum, &predicted);
		}
	}

	if (guiLive.show_selection_outline && screenState->selected &&
		(screenState->particleUnderClick != NULL)){
		const particlePrediction selected = cachedGraphicsParticlePrediction(
			screenState->particleUnderClick);
		const Vector3 center = {
			(float)selected.x, (float)selected.y, (float)selected.z,
		};
		const Color selectionColor = {245, 247, 250, 255};
		const float outlineRadius = (float)(1.08*selected.rad);
		if (camera_frustum_3d_intersects_sphere(frustum, center,
			outlineRadius)){
			DrawSphereWires(center, outlineRadius, 8, 12, selectionColor);
		}
	}
}
#endif

static void setGraphicalStructureFactorEnabled(bool enabled){
	if (enabled == structFactorActivated){
		return;
	}
	if (enabled){
		if (!allocateStructureFactorPositions()){
			releaseStructureFactorPositions();
			guiLive.graphical_structure_factor = false;
			return;
		}
		structFactorActivated = true;
		captureStructureFactorPositions();
		asyncStructFactor();
		lastStructLaunchTime = GetTime();
	}
	else{
		releaseStructureFactorPositions();
		structFactorActivated = false;
	}
}

static void updateGraphicalStructureFactor(void){
	setGraphicalStructureFactorEnabled(guiLive.graphical_structure_factor);
	if (!structFactorActivated){
		return;
	}
	if (finishedStructComputation()){
		awaitStructFactor();
		normalizeStruct();
		UpdateTexture(texture, structFactorGraphics);
		lastStructLaunchTime = GetTime();
	}
	bool sampleCountChanged = false;
	if (!structureFactorJobIsRunning() &&
		(structSampleCount != desiredStructureFactorSampleCount())){
		const int previousSampleCount = structSampleCount;
		if (allocateStructureFactorPositions()){
			sampleCountChanged = (structSampleCount != previousSampleCount);
		}
	}
	const double refreshInterval = validatedRenderRefreshSeconds(
		guiLive.structure_factor_refresh_seconds,
		STRUCTURE_FACTOR_DEFAULT_REFRESH_SECONDS);
	if (!structureFactorJobIsRunning() &&
		(sampleCountChanged || (guiLive.running &&
		GetTime() - lastStructLaunchTime >= refreshInterval))){
		captureStructureFactorPositions();
		asyncStructFactor();
		lastStructLaunchTime = GetTime();
	}
}

static void updateAppearance(state* screenState){
	const int previousMode = screenState->colorParam2;
	const int requestedMode = (int)guiLive.color_mode;
	screenState->particleArray[0] = colorFromRgba(guiLive.type0_rgba);
	screenState->particleArray[1] = colorFromRgba(guiLive.type1_rgba);
	screenState->arrow = guiLive.velocity_arrows;

	switch (guiLive.palette){
		case SIM_PALETTE_VIRIDIS: screenState->colorArray = Viridis; break;
		case SIM_PALETTE_COPPER: screenState->colorArray = Copper; break;
		case SIM_PALETTE_PLASMA:
		default: screenState->colorArray = Plasma; break;
	}

	screenState->colorParam2 = requestedMode;
	switch (guiLive.color_mode){
		case SIM_COLOR_BY_COLLISIONS:
			screenState->colorFunction = &colorCollision;
			break;
		case SIM_COLOR_BY_SPEED:
			screenState->colorFunction = &colorVelocity;
			break;
		#if !THREE_D
		case SIM_COLOR_BY_HEX_ORDER:
			nSym = 6;
			screenState->colorFunction = &colorBOOP;
			break;
		case SIM_COLOR_BY_SQUARE_ORDER:
			nSym = 4;
			screenState->colorFunction = &colorBOOP;
			break;
		#endif
		case SIM_COLOR_BY_RADIUS:
			screenState->colorFunction = &colorRadius;
			break;
		case SIM_COLOR_BY_CHARGE:
			screenState->colorFunction = &colorCharge;
			break;
		case SIM_COLOR_BY_TYPE:
		default:
			screenState->colorParam2 = 0;
			break;
	}

	if (screenState->colorParam2 != 0){
		ensureColorFunctionArray(screenState);
	}
	else if (screenState->colorFunctionArray != NULL){
		free(screenState->colorFunctionArray);
		screenState->colorFunctionArray = NULL;
		screenState->colorFunctionCapacity = 0;
	}
	if (previousMode != screenState->colorParam2){
		invalidateParticleColorCache();
	}
}

static GuiParticleSnapshot selectedParticleSnapshot(const state* screenState){
	GuiParticleSnapshot result = {0};
	if (!screenState->selected || screenState->particleUnderClick == NULL){
		return result;
	}
	const particle* p = screenState->particleUnderClick;
	const particlePrediction predicted = graphicsParticlePrediction(p);
	result.valid = true;
	result.index = p->num;
	result.type = p->type;
	result.charge = p->charge;
	result.collisions = (uint64_t)p->coll;
	result.x = predicted.x;
	result.y = predicted.y;
	result.vx = predicted.vx;
	result.vy = predicted.vy;
	#if THREE_D
	result.z = predicted.z;
	result.vz = predicted.vz;
	#endif
	result.radius = predicted.rad;
	result.mass = p->m;
	result.last_collision_time = p->lastColl;
	return result;
}

double graphicsSimulationTimeLimit(double currentSimulationTime){
	synchronizeSimulationPacing(true, currentSimulationTime);
	if (!pacingState.enabled)
		return HUGE_VAL;
	/* A manual single-step can put the event clock ahead of its pacing target.
	 * Never move time backwards; wall time will naturally catch up. */
	return fmax(currentSimulationTime, pacingState.targetSimulationTime);
}

void graphicsRebaseSimulationPacing(double currentSimulationTime){
	/* Growth time parameterizes an initializer; it is not elapsed physical
	 * dynamics.  Begin pacing production from the completed configuration
	 * instead of making it wait for wall time to catch up with preparation. */
	pacingState.initialized = false;
	synchronizeSimulationPacing(guiLive.running, currentSimulationTime);
}

static const GuiMetrics* refreshGuiMetrics(void){
	const double now = GetTime();
	synchronizeSimulationPacing(guiLive.running, t);
	guiMetrics.simulation_time_pacing = pacingState.enabled;
	guiMetrics.target_simulation_time_rate = pacingState.rate;
	guiMetrics.simulation_time_lag = pacingState.lag;
	if (lastMetricsWallTime >= 0 && now - lastMetricsWallTime < 0.2){
		guiMetrics.frames_per_second = (float)GetFPS();
		return &guiMetrics;
	}
	const double elapsed = (lastMetricsWallTime >= 0) ?
		fmax(now - lastMetricsWallTime, 1e-9) : 0;
	const unsigned long events = ncol + ncross;
	guiMetrics.dimensions = THREE_D ? 3 : 2;
	guiMetrics.particle_count = N;
	guiMetrics.cells_x = Nxcells;
	guiMetrics.cells_y = Nycells;
	#if THREE_D
	guiMetrics.cells_z = Nzcells;
	guiMetrics.box_lz = Lz;
	#else
	guiMetrics.cells_z = 1;
	guiMetrics.box_lz = 0;
	#endif
	guiMetrics.queued_events = 0;
	guiMetrics.event_capacity = 2*N + suppSizeEvent;
	if (eventList != NULL){
		for (int index = 0; index < guiMetrics.event_capacity; index++){
			if (eventList[index] != NULL && eventList[index]->q >= 0){
				guiMetrics.queued_events++;
			}
		}
	}
	guiMetrics.simulation_time = t;
	guiMetrics.box_lx = Lx;
	guiMetrics.box_ly = Ly;
	guiMetrics.packing_fraction = simulationCurrentPackingFraction(t);
	if (!isfinite(guiMetrics.packing_fraction))
		guiMetrics.packing_fraction = phi;
	guiMetrics.target_packing_fraction = phi;
	guiMetrics.pressure = pressure;
	/* Pressure is accumulated by the headless thermodynamic sampling path,
	 * which graphical runs intentionally do not schedule. */
	guiMetrics.pressure_valid = false;
	guiMetrics.bath_temperature = T;
	guiMetrics.collisions = (uint64_t)ncol;
	guiMetrics.crossings = (uint64_t)ncross;
	guiMetrics.frames_per_second = (float)GetFPS();
	guiMetrics.growing = freeFly == &freeFlyGrow
		#if POLYGON_PARTICLES && !THREE_D
		|| polygon_particles_growing()
		#endif
		;

	double translationalEnergy = 0;
	for (int i = 0; i < N; i++){
		const particlePrediction predicted = graphicsParticlePrediction(particles + i);
		double speedSquared = predicted.vx*predicted.vx + predicted.vy*predicted.vy;
		#if THREE_D
		speedSquared += predicted.vz*predicted.vz;
		#endif
		translationalEnergy += 0.5*particles[i].m*speedSquared;
	}
	double energy = translationalEnergy;
	#if POLYGON_PARTICLES && !THREE_D
	energy += polygon_particles_rotational_energy();
	#endif
	guiMetrics.total_energy = energy;
	guiMetrics.energy_per_particle = (N > 0) ? energy/N : 0;
	guiMetrics.temperature = (N > 0) ?
		2.0*translationalEnergy/((THREE_D ? 3.0 : 2.0)*N) : 0;
	if (elapsed > 0){
		guiMetrics.event_rate = (events - lastMetricsEvents)/elapsed;
		guiMetrics.simulation_time_rate =
			(t - lastMetricsSimulationTime)/elapsed;
	}
	lastMetricsWallTime = now;
	lastMetricsSimulationTime = t;
	lastMetricsEvents = events;
	return &guiMetrics;
}

static void requestedSceneTargetSize(int* width, int* height){
	const double supersampling = Clamp((float)guiLive.scene_supersampling,
		0.5f, 3.0f);
	#if defined(PLATFORM_WEB)
	/* raylib 5.5 reports a fixed 1x DPI scale on PLATFORM_WEB (and logs a
	 * warning whenever GetWindowScaleDPI() is queried).  Use the browser's
	 * actual backing-store ratio so the particle scene is still rendered with
	 * enough samples on Retina/HiDPI displays. */
	const double browserDpi = emscripten_get_device_pixel_ratio();
	const double dpiX = (isfinite(browserDpi) && browserDpi > 0) ?
		Clamp(browserDpi, 1.0, 4.0) : 1.0;
	const double dpiY = dpiX;
	#else
	const Vector2 dpi = GetWindowScaleDPI();
	const double dpiX = (isfinite(dpi.x) && dpi.x > 0) ?
		Clamp(dpi.x, 1.0f, 4.0f) : 1.0;
	const double dpiY = (isfinite(dpi.y) && dpi.y > 0) ?
		Clamp(dpi.y, 1.0f, 4.0f) : 1.0;
	#endif
	const float logicalWidth = (guiViewport.width > 16) ?
		guiViewport.width : (float)GetScreenWidth();
	const float logicalHeight = (guiViewport.height > 16) ?
		guiViewport.height : (float)GetScreenHeight();
	/* ImGui reports logical points, while this texture is sampled into the
	 * native framebuffer. Include the display scale so 1.0 means one scene
	 * pixel per physical display pixel. Quantization avoids reallocating a
	 * large GPU target for every pixel of a dock-splitter drag. */
	*width = ((int)ceilf((float)(logicalWidth*dpiX*supersampling)/16.0f))*16;
	*height = ((int)ceilf((float)(logicalHeight*dpiY*supersampling)/16.0f))*16;
	if (*width > 6144) *width = 6144;
	if (*height > 6144) *height = 6144;
	if (*width < 64) *width = 64;
	if (*height < 64) *height = 64;
}

static bool renderSimulationScene(window* screenWindow, state* screenState){
	int width = 0;
	int height = 0;
	requestedSceneTargetSize(&width, &height);
	if (!ensureSceneTarget(width, height)){
		return false;
	}
	BeginTextureMode(sceneTarget);
	ClearBackground(colorFromRgba(guiLive.background_rgba));
	#if THREE_D
	ensureCamera3D(screenWindow);
	screenWindow->cam3D.fovy = (float)guiLive.camera_fov_degrees;
	const bool predictionsReady = refreshRenderPredictions();
	const double maximumParticleRadius = predictionsReady ?
		maximumPredictedParticleRadius3D : NAN;
	const cameraClipPlanes3D clipPlanes = computeCameraClipPlanes3D(
		&screenWindow->cam3D, screenState->arrow, maximumParticleRadius);
	rlSetClipPlanes(clipPlanes.nearDistance, clipPlanes.farDistance);
	const double exactNearDistance = rlGetCullDistanceNear();
	const double exactFarDistance = rlGetCullDistanceFar();
	const double sceneAspect = (double)sceneTargetWidth/sceneTargetHeight;
	CameraFrustum3D frustum = {0};
	camera_frustum_3d_build_perspective(&frustum, &screenWindow->cam3D,
		sceneAspect, exactNearDistance, exactFarDistance);
	BeginMode3D(screenWindow->cam3D);
	drawParticlesAndBox3D(&screenWindow->cam3D, &frustum, screenState,
		predictionsReady);
	EndMode3D();
	#else
	if (camera2DNeedsReset){
		resetCamera2D(screenWindow);
	}
	const Vector2 cornerA = GetScreenToWorld2D((Vector2){0, 0}, screenWindow->cam);
	const Vector2 cornerB = GetScreenToWorld2D(
		(Vector2){(float)sceneTargetWidth, (float)sceneTargetHeight},
		screenWindow->cam);
	const Rectangle particleViewport = {
		.x = fminf(cornerA.x, cornerB.x),
		.y = fminf(cornerA.y, cornerB.y),
		.width = fabsf(cornerB.x - cornerA.x),
		.height = fabsf(cornerB.y - cornerA.y),
	};
	BeginMode2D(screenWindow->cam);
	drawParticlesAndBox(screenWindow->factor, screenState,
		screenWindow->cam.zoom, particleViewport);
	EndMode2D();
	#endif
	EndTextureMode();
	return true;
}

static void clearParticleSelection(state* screenState){
	clearPendingParticlePress();
	screenState->selected = 0;
	screenState->leftClicked = 0;
	screenState->particleUnderClick = NULL;
}

static void updateSelectedParticle(const GuiParticleUpdate* update,
	state* screenState){
	if (gui_imgui_simulation_faulted()){
		return;
	}
	if (update == NULL || update->index < 0 || update->index >= N){
		return;
	}
	particle* p = particles + update->index;
	#if (STRESS || HEAT)
	beginTransportVelocityBatch();
	#endif
	freeFly(p);
	const double oldVx = p->vx;
	const double oldVy = p->vy;
	#if THREE_D
	const double oldVz = p->vz;
	#else
	const double oldVz = 0;
	#endif
	p->vx = update->vx;
	p->vy = update->vy;
	#if THREE_D
	p->vz = update->vz;
	#endif
	#if (STRESS || HEAT)
	recordTransportVelocityChange(p, oldVx, oldVy, oldVz);
	#else
	(void)oldVx;
	(void)oldVy;
	(void)oldVz;
	#endif
	#if POLYGON_PARTICLES && !THREE_D
	if (polygon_particles_active()){
		if (!polygon_particles_set_mass(update->index, update->mass)) return;
	} else {
		p->m = update->mass;
	}
	#else
	p->m = update->mass;
	#endif
	p->coll++;
	removeEventFromQueue(eventList[N + p->num]);
	removeEventFromQueue(eventList[p->num]);
	crossingEvent(p->num);
	collisionEvent(p->num);
	physicalQ();
	live_analysis_reset_generation();
	screenState->particleUnderClick = p;
	screenState->selected = 1;
	invalidateParticleColorCache();
}

static void handleGuiCommand(const GuiCommand* command, window* screenWindow,
	state* screenState){
	if (command == NULL){
		return;
	}
	switch (command->type){
		case GUI_COMMAND_APPLY_CONFIG: {
			SimulationConfig current;
			char message[512] = {0};
			if (hostAppliedConfigValid){
				current = hostAppliedConfig;
			}
			else{
				simulation_config_capture_globals(&current);
			}
			const SimulationReconfigureAction action =
				simulation_reconfigure_core(&current, &command->config,
					message, sizeof(message));

			if (action == SIM_RECONFIGURE_FAILED){
				char retainedMessage[1024];
				snprintf(retainedMessage, sizeof(retainedMessage),
					"Live update rolled back. The previous valid configuration "
					"and particle state were retained; no partial changes were "
					"applied.%s%s",
					message[0] != '\0' ? "\n\n" : "",
					message[0] != '\0' ? message : "");
				TraceLog(LOG_ERROR, "%s", retainedMessage);
				gui_imgui_acknowledge_applied_config(
					NULL, GUI_CONFIG_APPLY_REJECTED, retainedMessage);
				break;
			}

			if (action != SIM_RECONFIGURE_RESTART_REQUIRED){
				/* In-place updates retain particle addresses, selection, camera,
				 * drag state, time, and counters.  Only derived displays whose
				 * meaning changed with the physics are refreshed. */
				gui_imgui_acknowledge_applied_config(
					&command->config, GUI_CONFIG_APPLY_SUCCEEDED, NULL);
				hostAppliedConfig = command->config;
				hostAppliedConfigValid = true;
				if (action != SIM_RECONFIGURE_UNCHANGED){
					live_analysis_reset_generation();
					invalidateParticleColorCache();
					lastMetricsWallTime = -1;
				}
				if (simulationHasFault()){
					graphicsPauseSimulation(screenState);
					graphicsReportSimulationFault(simulationFaultMessage());
				}
				break;
			}

			/* State-defining edits and runtime constraints which cannot contain
			 * the current particles intentionally start a new trajectory.  Reuse
			 * the explicit restart machinery, and tell the user why preservation
			 * was not possible. */
			char restartNotice[1024];
			snprintf(restartNotice, sizeof(restartNotice),
				"Configuration applied by rebuilding the simulation at t = 0. "
				"The previous particle state could not be preserved.%s%s",
				message[0] != '\0' ? "\n\n" : "",
				message[0] != '\0' ? message : "");
			char error[512] = {0};
			clearPendingParticlePress();
			if (activeParticleDrag.active){
				finishParticleDrag(screenState);
			}
			if (simulation_restart_core(&command->config, error, sizeof(error))){
				if (structFactorActivated){
					releaseStructureFactorPositions();
					structFactorActivated = false;
				}
				clearParticleSelection(screenState);
				free(screenState->colorFunctionArray);
				screenState->colorFunctionArray = NULL;
				screenState->colorFunctionCapacity = 0;
				SimulationConfig effective;
				simulation_config_capture_globals(&effective);
				effective.rendering = command->config.rendering;
				hostAppliedConfig = effective;
				hostAppliedConfigValid = true;
				gui_imgui_acknowledge_applied_config(
					&effective, GUI_CONFIG_APPLY_RESTARTED, restartNotice);
				screenState->running = command->config.rendering.start_running;
				guiLive.running = screenState->running;
				if (error[0] != '\0'){
					TraceLog(LOG_WARNING, "%s", error);
					graphicsReportRuntimeError(error);
				}
				if (simulationHasFault()){
					graphicsPauseSimulation(screenState);
					graphicsReportSimulationFault(simulationFaultMessage());
				}
				screenWindow->factor = sceneTargetHeight/Ly;
				#if THREE_D
				resetCamera3D(screenWindow);
				#else
				camera2DNeedsReset = true;
				#endif
				invalidateParticleColorCache();
				lastMetricsWallTime = -1;
			}
			else{
				char retainedMessage[1024];
				snprintf(retainedMessage, sizeof(retainedMessage),
					"Live update required a fresh initial state, but rebuilding "
					"it failed. The previous valid configuration and simulation "
					"state were retained; no partial changes were applied.%s%s",
					error[0] != '\0' ? "\n\n" : "",
					error[0] != '\0' ? error : "");
				TraceLog(LOG_ERROR, "%s", retainedMessage);
				gui_imgui_acknowledge_applied_config(
					NULL, GUI_CONFIG_APPLY_REJECTED, retainedMessage);
			}
			break;
		}
		case GUI_COMMAND_RESTART_APPLIED: {
			char error[512] = {0};
			clearPendingParticlePress();
			if (activeParticleDrag.active){
				finishParticleDrag(screenState);
			}
			if (simulation_restart_core(&command->config, error, sizeof(error))){
				if (structFactorActivated){
					releaseStructureFactorPositions();
					structFactorActivated = false;
				}
				clearParticleSelection(screenState);
				free(screenState->colorFunctionArray);
				screenState->colorFunctionArray = NULL;
				screenState->colorFunctionCapacity = 0;
				SimulationConfig effective;
				simulation_config_capture_globals(&effective);
				effective.rendering = command->config.rendering;
				hostAppliedConfig = effective;
				hostAppliedConfigValid = true;
				gui_imgui_acknowledge_applied_config(
					&effective, GUI_CONFIG_APPLY_RESTARTED, NULL);
				screenState->running = command->config.rendering.start_running;
				guiLive.running = screenState->running;
				if (error[0] != '\0'){
					TraceLog(LOG_WARNING, "%s", error);
					graphicsReportRuntimeError(error);
				}
				if (simulationHasFault()){
					graphicsPauseSimulation(screenState);
					graphicsReportSimulationFault(simulationFaultMessage());
				}
				screenWindow->factor = sceneTargetHeight/Ly;
				#if THREE_D
				resetCamera3D(screenWindow);
				#else
				camera2DNeedsReset = true;
				#endif
				invalidateParticleColorCache();
				lastMetricsWallTime = -1;
			}
			else{
				TraceLog(LOG_ERROR, "Simulation restart failed: %s", error);
				char retainedMessage[1024];
				snprintf(retainedMessage, sizeof(retainedMessage),
					"Restart failed. The previous valid configuration and simulation "
					"state were retained; no partial changes were applied.%s%s",
					error[0] != '\0' ? "\n\n" : "",
					error[0] != '\0' ? error : "");
				gui_imgui_acknowledge_applied_config(
					NULL, GUI_CONFIG_APPLY_REJECTED, retainedMessage);
			}
			break;
		}
		case GUI_COMMAND_SAVE_SNAPSHOT:
			saveSnapshotGraphics();
			break;
		case GUI_COMMAND_STEP_ONE_EVENT:
			screenState->stepRequested = true;
			break;
		case GUI_COMMAND_UPDATE_PARTICLE:
			updateSelectedParticle(&command->particle, screenState);
			break;
		case GUI_COMMAND_CLEAR_SELECTION:
			clearParticleSelection(screenState);
			break;
		case GUI_COMMAND_RESET_VIEW:
			#if THREE_D
			resetCamera3D(screenWindow);
			#else
			resetCamera2D(screenWindow);
			#endif
			break;
		case GUI_COMMAND_TOGGLE_FULLSCREEN:
			toggleFullscreenWindow(screenWindow);
			break;
		case GUI_COMMAND_NONE:
		default:
			break;
	}
}


void draw(int argc, char *argv[], window* screenWindow, state* screenState){
	(void)argc;
	(void)argv;
	if (!imguiInitialized){
		SimulationConfig applied;
		simulation_config_capture_globals(&applied);
		hostAppliedConfig = applied;
		hostAppliedConfigValid = true;
		imguiInitialized = gui_imgui_init(&applied);
		guiLive.running = screenState->running;
		if (imguiInitialized && pendingRuntimeError[0] != '\0'){
			if (pendingRuntimeFault){
				gui_imgui_report_fault(pendingRuntimeError);
			}
			else{
				gui_imgui_report_error(pendingRuntimeError);
			}
			pendingRuntimeError[0] = '\0';
			pendingRuntimeFault = false;
		}
		if (imguiInitialized && pendingTransientNotice[0] != '\0'){
			gui_imgui_report_notice(pendingTransientNotice);
			pendingTransientNotice[0] = '\0';
		}
	}
	/* Match the headless stop condition without closing the interactive
	 * workspace: reaching tmax pauses the event loop so the final state can be
	 * inspected or restarted.  Apply this after first-frame GUI initialization,
	 * which imports the pre-existing run state. */
	if (tmax < never && t > tmax){
		graphicsPauseSimulation(screenState);
	}

	updateAppearance(screenState);
	updateGraphicalStructureFactor();
	bool sceneReady = false;
	if (guiLive.show_scene){
		sceneReady = renderSimulationScene(screenWindow, screenState);
	}
	const GuiParticleSnapshot selected = selectedParticleSnapshot(screenState);
	const GuiMetrics* metrics = refreshGuiMetrics();

	BeginDrawing();
	ClearBackground(colorFromRgba(guiLive.background_rgba));
	GuiFrameResult frame = {0};
	if (imguiInitialized){
		frame = gui_imgui_frame(sceneReady ? &sceneTarget : NULL,
			structFactorActivated ? &texture : NULL,
			metrics, &selected, &guiLive);
		if (frame.viewport.width > 0 && frame.viewport.height > 0){
			guiViewport = frame.viewport;
		}
	}
	else if (sceneReady){
		DrawTexturePro(sceneTarget.texture,
			(Rectangle){0, 0, (float)sceneTargetWidth, (float)-sceneTargetHeight},
			(Rectangle){0, 0, (float)GetScreenWidth(), (float)GetScreenHeight()},
			(Vector2){0, 0}, 0, WHITE);
		guiViewport = (GuiViewport){0, 0, (float)GetScreenWidth(),
			(float)GetScreenHeight(), true, true};
	}
	/* On the web, DOM callbacks update raylib's current input state between
	 * frames and EndDrawing() advances it to "previous".  Consume scene input
	 * first so clicks, releases and wheel motion are not erased.  This ordering
	 * also gives native builds the viewport/capture state from this GUI frame. */
	getInput(screenWindow, screenState);
	EndDrawing();

	screenState->running = guiLive.running;
	synchronizeSimulationPacing(screenState->running, t);
	#if !defined(PLATFORM_WEB)
	if (guiLive.target_fps != appliedTargetFps && guiLive.target_fps > 0){
		SetTargetFPS(guiLive.target_fps);
		appliedTargetFps = guiLive.target_fps;
	}
	#else
	/* The browser callback performs non-blocking frame capping.  Calling
	 * SetTargetFPS() here would make raylib sleep on the JavaScript thread. */
	#endif
	handleGuiCommand(&frame.command, screenWindow, screenState);
}

int graphicsTargetFrameRate(void){
	return guiLive.target_fps > 0 ? guiLive.target_fps : 0;
}

void graphicsPauseSimulation(state* screenState){
	if (screenState != NULL){
		screenState->running = false;
		screenState->stepRequested = false;
	}
	guiLive.running = false;
}

void graphicsReportRuntimeError(const char* message){
	if (message == NULL || message[0] == '\0'){
		return;
	}
	if (imguiInitialized){
		gui_imgui_report_error(message);
		return;
	}
	snprintf(pendingRuntimeError, sizeof(pendingRuntimeError), "%s", message);
	pendingRuntimeFault = false;
}

void graphicsReportTransientNotice(const char* message){
	if (message == NULL || message[0] == '\0'){
		return;
	}
	if (imguiInitialized){
		gui_imgui_report_notice(message);
		return;
	}
	snprintf(pendingTransientNotice, sizeof(pendingTransientNotice), "%s",
		message);
}

void graphicsReportSimulationFault(const char* message){
	if (message == NULL || message[0] == '\0'){
		return;
	}
	if (imguiInitialized){
		gui_imgui_report_fault(message);
		return;
	}
	snprintf(pendingRuntimeError, sizeof(pendingRuntimeError), "%s", message);
	pendingRuntimeFault = true;
}

#if defined(PLATFORM_WEB)
static void downloadBrowserFile(const char* path){
	if (path == NULL || path[0] == '\0'){
		return;
	}
	EM_ASM({
		const path = UTF8ToString($0);
		try {
			const bytes = FS.readFile(path);
			const blob = new Blob([bytes], {type: 'text/plain;charset=utf-8'});
			const url = URL.createObjectURL(blob);
			const anchor = document.createElement('a');
			anchor.href = url;
			anchor.download = path.split('/').pop() || 'EDMD-snapshot.dump';
			anchor.style.display = 'none';
			document.body.appendChild(anchor);
			anchor.click();
			anchor.remove();
			setTimeout(() => URL.revokeObjectURL(url), 1000);
		} catch (error) {
			console.error('Unable to download EDMD snapshot', error);
		}
	}, path);
}
#endif

void saveSnapshotGraphics(){
	for (int i = 0; i < N; i++){
		freeFly(particles + i);
	}
	char snapshotPath[350] = {0};
	if (!customSnapshotName(snapshotPath, sizeof(snapshotPath))){
		graphicsReportRuntimeError("Unable to construct a snapshot filename");
		return;
	}
	FILE* periodicStream = fichier;
	FILE* snapshotStream = fopen(snapshotPath, "w");
	if (snapshotStream == NULL){
		char message[512];
		snprintf(message, sizeof(message),
			"Unable to save snapshot '%s': %s", snapshotPath, strerror(errno));
		TraceLog(LOG_ERROR, "%s", message);
		graphicsReportRuntimeError(message);
		return;
	}
	fichier = snapshotStream;
	saveTXT();
	const bool snapshotClosed = fclose(snapshotStream) == 0;
	if (!snapshotClosed){
		char message[512];
		snprintf(message, sizeof(message),
			"Unable to finish writing snapshot '%s': %s",
			snapshotPath, strerror(errno));
		TraceLog(LOG_ERROR, "%s", message);
		graphicsReportRuntimeError(message);
	}
	fichier = periodicStream;
	#if defined(PLATFORM_WEB)
	if (snapshotClosed){
		downloadBrowserFile(snapshotPath);
	}
	#endif
}

double colorVelocity(particle* p){
	const particlePrediction predicted = cachedGraphicsParticlePrediction(p);
	double speedSquared = predicted.vx*predicted.vx + predicted.vy*predicted.vy;
	#if THREE_D
	speedSquared += predicted.vz*predicted.vz;
	#endif
	return sqrt(speedSquared);
}

double colorCollision(particle* p){
	return (double)p->coll;
}

double colorRadius(particle* p){
	return cachedGraphicsParticlePrediction(p).rad;
}


double colorCharge(particle* p){
	return p->charge;
}

#if !THREE_D
double colorBOOP(particle* p1){
	double treshold = 2;
	if ((addWell) && (U > 0)){
		treshold *= sig;
	}
	if (p1->type != 1){
		return 0;
	}
	int count = 0;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	double re = 0;
	double im = 0;
	const particlePrediction p1Now = cachedGraphicsParticlePrediction(p1);
	/* The predicate below has a center-distance cutoff of
	 * sqrt(treshold)*(r1 + r2), not one linked-cell width. */
	const double maximumCutoff = sqrt(treshold)*
		(p1Now.rad + graphicsMaximumParticleRadius());
	const int rangeX = graphicsCellNeighborRange(maximumCutoff, Lx, Nxcells);
	const int rangeY = graphicsCellNeighborRange(maximumCutoff, Ly, Nycells);
	const int numberX = graphicsNeighborCellCount(Nxcells, rangeX);
	const int numberY = graphicsNeighborCellCount(Nycells, rangeY);
	for (int yi = 0; yi < numberY; yi++){
		const int neighborY = graphicsNeighborCellCoordinate(Y, Nycells, rangeY, yi);
		for (int xi = 0; xi < numberX; xi++){
			const int neighborX = graphicsNeighborCellCoordinate(X, Nxcells, rangeX, xi);
			particle* p2 = cellList[neighborY*Nxcells + neighborX];
			while (p2 != NULL){
				if ((p1->num != p2->num) && (p2->type == 1)){
					const particlePrediction p2Now = cachedGraphicsParticlePrediction(p2);
					double dy = p1Now.y - p2Now.y;
					double dx = p1Now.x - p2Now.x;
					graphicsPBCXY(&dx, &dy);
					const double distanceSquared = dx*dx + dy*dy;
					const double radiusSum = p1Now.rad + p2Now.rad;
					if ((distanceSquared > 0) &&
						(distanceSquared < treshold*radiusSum*radiusSum)){
						count++;
						const double inverseDistance = 1.0/sqrt(distanceSquared);
						const double directionReal = dx*inverseDistance;
						const double directionImaginary = dy*inverseDistance;
						double harmonicReal = 1;
						double harmonicImaginary = 0;
						for (int harmonic = 0; harmonic < nSym; harmonic++){
							const double nextReal = harmonicReal*directionReal -
								harmonicImaginary*directionImaginary;
							harmonicImaginary = harmonicReal*directionImaginary +
								harmonicImaginary*directionReal;
							harmonicReal = nextReal;
						}
						re += harmonicReal;
						im += harmonicImaginary;
					}
				}
				p2 = p2->nxt;
			}
		}
	}
	if (count >= 1)
		return (re*re + im*im)/count;
	return 0;
}
#endif

Color colorSelect(double value, Color* colorArray){
	if (!isfinite(value)) value = 0;
	if (value < 0) value = 0;
	if (value > 1) value = 1;
	return colorArray[(int)(value*color_size)];
}

#if !defined(PLATFORM_WEB)
static bool structureFactorJobWasCancelled(unsigned long generation){
	pthread_mutex_lock(&structMutex);
	const bool cancelled = structThreadPoolStopping ||
		structJobCancelRequested || generation != structJobGeneration;
	pthread_mutex_unlock(&structMutex);
	return cancelled;
}
#endif

static bool computeStructureFactorRow(int row, unsigned long generation){
	if (row < 0 || row >= qN || structSampleCount <= 0 || positions == NULL){
		return false;
	}
	const double qStep = qxG[1] - qxG[0];
	double rowReal[qN] = {0};
	double rowImaginary[qN] = {0};
	for (int n = 0; n < structSampleCount; n++){
		#if !defined(PLATFORM_WEB)
		if ((n & 63) == 0 && structureFactorJobWasCancelled(generation)){
			return false;
		}
		#else
		(void)generation;
		#endif
		const position* p = positions + n;
		const double firstPhase = qxG[row]*p->x + qxG[0]*p->y;
		const double phaseStep = qStep*p->y;
		double phaseReal = cos(firstPhase);
		double phaseImaginary = sin(firstPhase);
		const double stepReal = cos(phaseStep);
		const double stepImaginary = sin(phaseStep);

		for (int column = 0; column < qN; column++){
			rowReal[column] += phaseReal;
			rowImaginary[column] += phaseImaginary;
			const double nextReal = phaseReal*stepReal -
				phaseImaginary*stepImaginary;
			phaseImaginary = phaseReal*stepImaginary +
				phaseImaginary*stepReal;
			phaseReal = nextReal;
		}
	}

	for (int column = 0; column < qN; column++){
		const double magnitude = (rowReal[column]*rowReal[column] +
			rowImaginary[column]*rowImaginary[column])/structSampleCount;
		structFactorGraphics[row][column] =
			(float)log(magnitude*magnitude + 1.0);
		if ((row <= qN/2 + 3) && (row >= qN/2 - 3) &&
			(column >= qN/2 - 3) && (column <= qN/2 + 3)){
			structFactorGraphics[row][column] = 0;
		}
	}
	return true;
}

#if !defined(PLATFORM_WEB)
void* computeStructureFactorGraphic(void* arg){
	threadArg* threadArgument = (threadArg*)arg;
	unsigned long completedGeneration = 0;

	for (;;){
		pthread_mutex_lock(&structMutex);
		while (!structThreadPoolStopping && (completedGeneration == structJobGeneration)){
			pthread_cond_wait(&structJobReady, &structMutex);
		}
		if (structThreadPoolStopping){
			pthread_mutex_unlock(&structMutex);
			break;
		}
		completedGeneration = structJobGeneration;
		pthread_mutex_unlock(&structMutex);

		bool cancelled = false;
		for (int i = threadArgument->start;
			i < threadArgument->end && !cancelled; i++){
			cancelled = !computeStructureFactorRow(i, completedGeneration);
		}

		pthread_mutex_lock(&structMutex);
		structWorkersRemaining--;
		if (structWorkersRemaining == 0){
			pthread_cond_broadcast(&structJobFinished);
		}
		pthread_mutex_unlock(&structMutex);
	}
	return NULL;
}
#endif

void normalizeStruct(){
	float max = 0;
	for (int i = 0; i < qN; i++){
		for (int j = 0; j < qN; j++){
			float value = structFactorGraphics[i][j];
			if (!isfinite(value) || (value < 0)){
				value = 0;
				structFactorGraphics[i][j] = 0;
			}
			if (value > max){
				max = value;
			}
		}
	}
	if (!isfinite(max) || (max <= 0)){
		return;
	}

	for (int i = 0; i < qN; i++){
		for (int j = 0; j < qN; j++){
			structFactorGraphics[i][j] = structFactorGraphics[i][j]/max;
		}
	}
}

#if !defined(PLATFORM_WEB)
void threadPoolInit(){
	const int iterationsPerThread = qN/NUM_THREADS;
	const int remainingIterations = qN%NUM_THREADS;
	for (int i = 0; i < NUM_THREADS; i++) {
		threadArgs[i].start = i*iterationsPerThread;
		threadArgs[i].end = (i + 1)*iterationsPerThread;
		if (i == NUM_THREADS - 1){
			threadArgs[i].end += remainingIterations;
		}
		if (pthread_create(&threads[i], NULL, computeStructureFactorGraphic,
			(void*)&threadArgs[i]) != 0){
			TraceLog(LOG_FATAL, "Unable to start structure-factor worker %d", i);
		}
	}
	structThreadPoolStarted = true;
}

bool finishedStructComputation(){
	pthread_mutex_lock(&structMutex);
	const bool finished = structJobRunning && (structWorkersRemaining == 0);
	pthread_mutex_unlock(&structMutex);
	return finished;
}

void asyncStructFactor(){
	if (!structThreadPoolStarted || (structSampleCount <= 0) || (positions == NULL)){
		return;
	}
	pthread_mutex_lock(&structMutex);
	if (!structJobRunning){
		structWorkersRemaining = NUM_THREADS;
		structJobRunning = true;
		structJobCancelRequested = false;
		structJobGeneration++;
		pthread_cond_broadcast(&structJobReady);
	}
	pthread_mutex_unlock(&structMutex);
}

void awaitStructFactor(){
	if (!structThreadPoolStarted){
		return;
	}
	pthread_mutex_lock(&structMutex);
	while (structJobRunning && (structWorkersRemaining > 0)){
		pthread_cond_wait(&structJobFinished, &structMutex);
	}
	structJobRunning = false;
	pthread_mutex_unlock(&structMutex);
}

static void shutdownStructureFactorThreadPool(void){
	if (!structThreadPoolStarted){
		return;
	}
	pthread_mutex_lock(&structMutex);
	structThreadPoolStopping = true;
	structJobCancelRequested = true;
	pthread_cond_broadcast(&structJobReady);
	pthread_mutex_unlock(&structMutex);
	awaitStructFactor();
	for (int i = 0; i < NUM_THREADS; i++){
		pthread_join(threads[i], NULL);
	}
	structThreadPoolStarted = false;
}
#else
/* Static hosting should not require SharedArrayBuffer or COOP/COEP.  Spread
 * structure-factor rows over animation frames instead of spawning pthreads;
 * the numerical kernel and result are identical to the desktop workers. */
void threadPoolInit(){
	structWebNextRow = 0;
	structJobRunning = false;
	structJobCancelRequested = false;
	structThreadPoolStopping = false;
	structThreadPoolStarted = true;
}

bool finishedStructComputation(){
	if (!structThreadPoolStarted || !structJobRunning){
		return false;
	}
	if (structJobCancelRequested || structThreadPoolStopping){
		structJobRunning = false;
		structWebNextRow = 0;
		return false;
	}

	const double frameStart = GetTime();
	do{
		if (!computeStructureFactorRow(structWebNextRow, structJobGeneration)){
			structJobRunning = false;
			structWebNextRow = 0;
			return false;
		}
		structWebNextRow++;
	} while (structWebNextRow < qN && GetTime() - frameStart < 0.0025);

	return structWebNextRow >= qN;
}

void asyncStructFactor(){
	if (!structThreadPoolStarted || structJobRunning ||
		structSampleCount <= 0 || positions == NULL){
		return;
	}
	structJobCancelRequested = false;
	structJobGeneration++;
	structWebNextRow = 0;
	structJobRunning = true;
}

void awaitStructFactor(){
	if (!structThreadPoolStarted){
		return;
	}
	structJobRunning = false;
	structJobCancelRequested = false;
	structWebNextRow = 0;
}

static void shutdownStructureFactorThreadPool(void){
	if (!structThreadPoolStarted){
		return;
	}
	structThreadPoolStopping = true;
	structJobCancelRequested = true;
	awaitStructFactor();
	structThreadPoolStarted = false;
}
#endif



window graphicalInit(){
	/* Keep the executable's initial engine state and the GUI's default
	 * model sourced from one schema.  Command-line arguments are applied later
	 * by constantInit(), before the GUI captures the effective configuration. */
	SimulationConfig defaults;
	const BuildCapabilities capabilities = simulation_build_capabilities();
	simulation_config_set_defaults(&defaults, &capabilities);
	simulation_config_apply_globals(&defaults);
	Camera2D cam = {0};
    cam.zoom = 1;
	#if THREE_D
	Camera3D cam3D = {
		.position = {0},
		.target = {0},
		.up = {0, 1, 0},
		.fovy = 45.0f,
		.projection = CAMERA_PERSPECTIVE,
	};
	#endif

	for (int i = 0; i < qN; i++){
		qxG[i] = -10 + 20.0*i/(double)qN;
	}
	
	threadPoolInit();
	/* Request a native-resolution framebuffer on Retina/HiDPI displays.
	 * rlImGui uses GetWindowScaleDPI() for its framebuffer scale and the scene
	 * has its own bounded supersampled target, so both UI text and 3D edges stay
	 * sharp without changing their logical size. */
	SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_WINDOW_HIGHDPI |
		FLAG_MSAA_4X_HINT);
	
	#if THREE_D
	InitWindow(1800, 900, "EDMD - 3D system");
	#else
	InitWindow(1800, 900, "EDMD");
	#endif
	{
		#if defined(PLATFORM_WEB)
		const float browserDpi = (float)emscripten_get_device_pixel_ratio();
		const Vector2 startupDpi = {
			(isfinite(browserDpi) && browserDpi > 0) ? browserDpi : 1.0f,
			(isfinite(browserDpi) && browserDpi > 0) ? browserDpi : 1.0f,
		};
		#else
		const Vector2 startupDpi = GetWindowScaleDPI();
		#endif
		TraceLog(LOG_INFO,
			"EDMD: logical %dx%d, framebuffer %dx%d, DPI scale %.2fx%.2f",
			GetScreenWidth(), GetScreenHeight(), GetRenderWidth(), GetRenderHeight(),
			startupDpi.x, startupDpi.y);
	}
	#if !defined(PLATFORM_WEB)
	SetWindowMinSize(800, 450);
	const int display = GetCurrentMonitor();
	const int monitorWidth = GetMonitorWidth(display);
	const int monitorHeight = GetMonitorHeight(display);
	const int availableWidth = (monitorWidth > 80) ? monitorWidth - 80 : monitorWidth;
	const int availableHeight = (monitorHeight > 100) ? monitorHeight - 100 : monitorHeight;
	const int initialWidth = (1800 < availableWidth) ? 1800 : availableWidth;
	const int initialHeight = (900 < availableHeight) ? 900 : availableHeight;
	if ((initialWidth != GetScreenWidth()) || (initialHeight != GetScreenHeight())){
		SetWindowSize(initialWidth, initialHeight);
	}
	const Vector2 monitorPosition = GetMonitorPosition(display);
	SetWindowPosition((int)monitorPosition.x + (monitorWidth - initialWidth)/2,
		(int)monitorPosition.y + (monitorHeight - initialHeight)/2);
	SetTargetFPS(144);
	#endif
	image = (Image){
		.data = structFactorGraphics,
		.width = qN,
		.height = qN,
		.mipmaps = 1,
		.format = PIXELFORMAT_UNCOMPRESSED_R32,
	};
	texture = LoadTextureFromImage(image);
	
	#if !THREE_D
	arrowMaterial = LoadMaterialDefault();
	arrowMaterialLoaded = (arrowMaterial.maps != NULL);
	loadProceduralCircleMaterial();
	#endif

	const int initialScreenWidth = GetScreenWidth();
	const int initialScreenHeight = GetScreenHeight();
    return (window){
		.screenWidth = initialScreenWidth,
		.screenHeight = initialScreenHeight,
		.factor = 1.0,
		.start = initialScreenHeight,
		.xGUI = initialScreenWidth/1800.0,
		.yGUI = initialScreenHeight/900.0,
        .cam = cam,
		#if THREE_D
		.cam3D = cam3D,
		.cameraBoxLx = -1,
		.cameraBoxLy = -1,
		.cameraBoxLz = -1,
		#endif
    };
}

state GUIinit(){
	return (state){
		.leftClicked = 0,
		.selected = 0,
		.particleUnderClick = NULL,
		.particleArray = {{66, 150, 250, 255}, {242, 157, 73, 255}},
		.particleColor = BLACK,
		.colorEditing = false,
		.colorParam = 0,
		.colorEditing2 = false,
		.colorParam2 = 0,
		.arrow = false,
		.colorArray = Plasma,
		.colorFunction = colorCollision,
		.colorFunctionArray = NULL,
		.colorFunctionCapacity = 0,
		.running = true,
		.stepRequested = false,
		.wallMoving = false,
		.spacePressed = false,
		.editVx = false,
		.editVy = false,
		#if THREE_D
		.editVz = false,
		#endif
		.editM = false,
		.wallEditing = false,
		.wallParam = 0,
	};
}

void graphicsFree(state* screenState){
	if (structFactorActivated)
		releaseStructureFactorPositions();
	shutdownStructureFactorThreadPool();
	#if THREE_D
	unloadSphereBatch(&particleSphereBatch);
	#endif
	if (sceneTarget.texture.id != 0){
		UnloadRenderTexture(sceneTarget);
		sceneTarget = (RenderTexture2D){0};
		sceneTargetWidth = 0;
		sceneTargetHeight = 0;
	}
	if (imguiInitialized){
		gui_imgui_shutdown();
		imguiInitialized = false;
	}
	#if !THREE_D
	#if POLYGON_PARTICLES
	unloadPolygonTriangleBatch(&particlePolygonBatch);
	#endif
	unloadCircleBatch(&particleCircleBatch);
	unloadVelocityArrowBatch(&arrowBatch);
	if (proceduralCircleMaterialLoaded){
		UnloadMaterial(proceduralCircleMaterial);
		proceduralCircleMaterial = (Material){0};
		proceduralCircleMaterialLoaded = false;
	}
	if (arrowMaterialLoaded){
		UnloadMaterial(arrowMaterial);
		arrowMaterial = (Material){0};
		arrowMaterialLoaded = false;
	}
	#endif
	free(renderPredictions);
	renderPredictions = NULL;
	renderPredictionCapacity = 0;
	UnloadTexture(texture);
	free(screenState->colorFunctionArray);
	screenState->colorFunctionArray = NULL;
	screenState->colorFunctionCapacity = 0;
	invalidateParticleColorCache();
}
