#define RAYGUI_IMPLEMENTATION
#include "raygui.h" 
#include "color.h"
#include "graphics.h"
#include "EDMD.h"
#include "tex_data.h"  
#include "raymath.h"
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#define qN 180
#define SIZE 400
#define NUM_THREADS 4
#define STRUCTURE_FACTOR_MAX_SAMPLES 4096
#define STRUCTURE_FACTOR_REFRESH_SECONDS 0.5
#define COLOR_REFRESH_SECONDS (1.0/30.0)
#if !THREE_D
#define BOOP_COLOR_REFRESH_SECONDS 0.2
int nSym = 4;
#endif

double qxG[qN] = {0};
float structFactorGraphics[qN][qN] = {0};
bool structFactorActivated = 0;
position* positions;
pthread_t threads[NUM_THREADS];
threadArg threadArgs[NUM_THREADS];

typedef struct spriteBatch spriteBatch;
struct spriteBatch{
	Mesh mesh;
	int capacity;
	int count;
	bool loaded;
};

static spriteBatch particleBatch = {0};
static spriteBatch arrowBatch = {0};
#if THREE_D
typedef struct sphereBatch sphereBatch;
struct sphereBatch{
	Mesh mesh;
	Material material;
	Shader shader;
	Matrix* transforms;
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
static Material spriteMaterial = {0};
static Texture2D defaultMaterialTexture = {0};
static bool spriteMaterialLoaded = false;
static particlePrediction* renderPredictions = NULL;
static int renderPredictionCapacity = 0;

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

static pthread_mutex_t structMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t structJobReady = PTHREAD_COND_INITIALIZER;
static pthread_cond_t structJobFinished = PTHREAD_COND_INITIALIZER;
static unsigned long structJobGeneration = 0;
static int structWorkersRemaining = 0;
static int structSampleCount = 0;
static bool structJobRunning = false;
static bool structThreadPoolStopping = false;
static bool structThreadPoolStarted = false;
static double lastStructLaunchTime = -STRUCTURE_FACTOR_REFRESH_SECONDS;

Image image;
Texture2D texture;

Texture2D circleTexture;
Texture2D arrowTexture;

#if THREE_D
static RenderTexture2D sceneTarget = {0};
static int sceneTargetWidth = 0;
static int sceneTargetHeight = 0;
#endif

static particlePrediction graphicsParticlePrediction(const particle* p);
static void invalidateParticleColorCache(void);

#if THREE_D
static void resetCamera3D(window* screenWindow);
static void ensureCamera3D(window* screenWindow);
static bool ensureSceneTarget3D(int width, int height);
static void sceneTargetSize3D(int logicalWidth, int logicalHeight,
	int* pixelWidth, int* pixelHeight);
static void drawParticlesAndBox3D(state* screenState);
#endif

#if !THREE_D
static void graphicsPBCXY(double* dx, double* dy){
	PBC(dx, dy);
}
#endif

int getWhatsUnderClick(window* screenWindow, state* screenState){
	const Vector2 mouse = GetMousePosition();
	if ((mouse.x < 0) || (mouse.x >= screenWindow->start) ||
		(mouse.y < 0) || (mouse.y >= GetScreenHeight())){
		return -1;
	}

	#if THREE_D
	ensureCamera3D(screenWindow);
	int scenePixelWidth = 0;
	int scenePixelHeight = 0;
	sceneTargetSize3D(screenWindow->start, GetScreenHeight(),
		&scenePixelWidth, &scenePixelHeight);
	const Vector2 sceneMouse = {
		mouse.x*scenePixelWidth/(float)screenWindow->start,
		mouse.y*scenePixelHeight/(float)GetScreenHeight(),
	};
	const Ray ray = GetScreenToWorldRayEx(sceneMouse, screenWindow->cam3D,
		scenePixelWidth, scenePixelHeight);
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
	Vector2 pos = GetScreenToWorld2D(mouse, screenWindow->cam);
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

		for (int j = -1; j <= 1; j++){
			for (int k = -1; k <= 1; k++){
				screenState->particleUnderClick =  cellList[PBCcellY(Y + j) * Nxcells + PBCcellX(X + k)];
				while (screenState->particleUnderClick != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					const particlePrediction predicted = graphicsParticlePrediction(screenState->particleUnderClick);
					double dx = predicted.x - x;
					double dy = predicted.y - y;
					graphicsPBCXY(&dx, &dy);
					if (dx*dx + dy*dy < predicted.rad*predicted.rad){
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

#if THREE_D
static float boxRadius3D(void){
	const double diagonalSquared = Lx*Lx + Ly*Ly + Lz*Lz;
	return (float)(0.5*sqrt(fmax(diagonalSquared, 1.0)));
}

static void resetCamera3D(window* screenWindow){
	const Vector3 center = {
		(float)(0.5*Lx), (float)(0.5*Ly), (float)(0.5*Lz),
	};
	const Vector3 viewDirection = Vector3Normalize((Vector3){1.25f, 0.9f, 1.4f});
	const float distance = fmaxf(3.0f*boxRadius3D(), 5.0f);
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
}

static void ensureCamera3D(window* screenWindow){
	if ((screenWindow->cameraBoxLx <= 0) ||
		(fabs(screenWindow->cameraBoxLx - Lx) > 1e-9) ||
		(fabs(screenWindow->cameraBoxLy - Ly) > 1e-9) ||
		(fabs(screenWindow->cameraBoxLz - Lz) > 1e-9)){
		resetCamera3D(screenWindow);
	}
}

static bool ensureSceneTarget3D(int width, int height){
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
		TraceLog(LOG_ERROR, "Unable to allocate the 3D scene render target");
		return false;
	}
	sceneTargetWidth = width;
	sceneTargetHeight = height;
	SetTextureFilter(sceneTarget.texture, TEXTURE_FILTER_BILINEAR);
	return true;
}

static void sceneTargetSize3D(int logicalWidth, int logicalHeight,
	int* pixelWidth, int* pixelHeight){
	const int screenWidth = GetScreenWidth();
	const int screenHeight = GetScreenHeight();
	float scaleX = (screenWidth > 0) ?
		(float)GetRenderWidth()/(float)screenWidth : 1.0f;
	float scaleY = (screenHeight > 0) ?
		(float)GetRenderHeight()/(float)screenHeight : 1.0f;

	/* Match Retina/HiDPI framebuffer pixels. On a 1x display, modest
	 * supersampling keeps thin 3D edges from becoming stair-stepped. */
	scaleX = fmaxf(scaleX, 1.5f);
	scaleY = fmaxf(scaleY, 1.5f);
	*pixelWidth = (int)ceilf(logicalWidth*scaleX);
	*pixelHeight = (int)ceilf(logicalHeight*scaleY);
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

void getInput(window* screenWindow, state* screenState){
	#if !THREE_D
	particle* particleUnderClick = screenState->particleUnderClick;
	#endif

	if (IsWindowResized()){
		updateWindowLayout(screenWindow, !IsWindowFullscreen());
	}

	if (IsKeyPressed(KEY_F)){
		const int display = GetCurrentMonitor();
		if (IsWindowFullscreen()){
			ToggleFullscreen();
			SetWindowSize(screenWindow->screenWidth,
				screenWindow->screenHeight);
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
	ensureCamera3D(screenWindow);
	const Vector2 mouse = GetMousePosition();
	const bool mouseInScene = (mouse.x >= 0) && (mouse.x < screenWindow->start) &&
		(mouse.y >= 0) && (mouse.y < GetScreenHeight());
	const float wheel = GetMouseWheelMove();
	if (mouseInScene && (wheel != 0)){
		zoomCamera3D(screenWindow, wheel);
	}
	if (mouseInScene && IsMouseButtonDown(MOUSE_BUTTON_RIGHT)){
		orbitCamera3D(screenWindow, GetMouseDelta());
	}
	if (IsKeyPressed(KEY_R)){
		resetCamera3D(screenWindow);
	}
	if (mouseInScene && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)){
		const int result = getWhatsUnderClick(screenWindow, screenState);
		screenState->selected = (result == 1);
		screenState->leftClicked = (result == 1);
		if (result != 1){
			screenState->particleUnderClick = NULL;
		}
	}
	if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)){
		screenState->leftClicked = 0;
	}
	#else
	float wheel = GetMouseWheelMove();
	if (wheel != 0){
		Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);
		screenWindow->cam.offset = GetMousePosition();
		screenWindow->cam.target = mouseWorldPos;

		// evil trick
		double speed = 0.0625f*100/Ly;
		screenWindow->cam.zoom += wheel * speed*(screenWindow->cam.zoom*screenWindow->cam.zoom);
		if (screenWindow->cam.zoom < 0.2)
			screenWindow->cam.zoom = 0.2;
		if (screenWindow->cam.zoom > 10)
			screenWindow->cam.zoom = 10;
	}

	if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)){
		Vector2 delta = GetMouseDelta();
		float scale = -1./screenWindow->cam.zoom;
		delta.x = delta.x*scale;
		delta.y = delta.y*scale;
		screenWindow->cam.target.x += delta.x;
		screenWindow->cam.target.y += delta.y;
	}

	if (IsKeyDown(KEY_R)){
		screenWindow->cam.zoom = 1;
		screenWindow->cam.offset = (Vector2){0, 0};
		screenWindow->cam.target = (Vector2){0, 0};
	}
	if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)){
		
		if (screenState->leftClicked > 0){
			freeFly(particleUnderClick);
			invalidateParticleColorCache();
			Vector2 pos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);

			particleUnderClick->coll++;
		
			double x = pos.x/screenWindow->factor;
			double y = Ly - pos.y/screenWindow->factor; 
			double dx = (particleUnderClick->x - x);
			double dy = (particleUnderClick->y - y);
			graphicsPBCXY(&dx, &dy);

			particleUnderClick->vx -= dx*0.3 + 0.3*particleUnderClick->vx;

			particleUnderClick->vy -= dy*0.3 + 0.3*particleUnderClick->vy;

			
			removeEventFromQueue(eventList[particleUnderClick->num]);
			crossingEvent(particleUnderClick->num);

			removeEventFromQueue(eventList[N + particleUnderClick->num]);
			collisionEvent(particleUnderClick->num);
		}
		else{
			int resultat = getWhatsUnderClick(screenWindow, screenState);
			if (resultat == 2){
				screenState->wallMoving = true;
				Vector2 pos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);
				if (Lx < pos.x/screenWindow->factor)
					Lx = pos.x/screenWindow->factor;
			}
			if (resultat == 3){
				screenState->wallMoving = true;
				Vector2 pos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);
				if (Ly < (Ly - pos.y/screenWindow->factor)){
					Ly = (Ly - pos.y/screenWindow->factor);
					screenWindow->cam.target.y += fabs(pos.y);
				}
			}
			if (resultat == 1){
				GuiLock();
				screenState->leftClicked = 1;
				screenState->selected = 1;
			}
			else if (resultat == 0){
				screenState->selected = 0;
			}
		}
	}
	else{
		GuiUnlock();
		screenState->leftClicked = 0;
		if (screenState->wallMoving){
			screenState->wallMoving = false;
			invalidateParticleColorCache();
			for(int i = 0; i < N; i++){
				particle* p = particles + i;
				freeFly(p);
				if (addWell){
					free(p->particlesInWell);
				}
				
			}
			free(cellList);
			boxConstantHelper();
			cellListInit();
			if (addWell)
				pairsInit();
			
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				removeEventFromQueue(eventList[N + p->num]);
				removeEventFromQueue(eventList[p->num]);
				crossingEvent(p->num);
				collisionEvent(p->num);
			}
		}
	}
	#endif

	if (IsKeyPressed(KEY_SPACE) && !screenState->spacePressed){
		screenState->spacePressed = true;
		screenState->running = !screenState->running;
	}
	
	if (IsKeyReleased(KEY_SPACE)){
		screenState->spacePressed = false;
	}
	if (IsKeyReleased(KEY_S)){
		saveSnapshotGraphics();
	}
}

void GuiSliderBarDouble(Rectangle bounds, const char *textLeft, const char *textRight, double *value, double minValue, double maxValue){
	const float oldFloatValue = (float)*value;
	float floatValue = oldFloatValue;
	GuiSliderBar(bounds, textLeft, textRight, &floatValue, minValue, maxValue);
	if (floatValue != oldFloatValue){
		*value = (double)floatValue;
	}
}

static void invalidateParticleColorCache(void){
	colorCache.valid = false;
}

static particlePrediction graphicsParticlePrediction(const particle* p){
	particlePrediction prediction;
	if (freeFly == freeFlyGrow){
		prediction = predictFreeFlyGrow(p, t);
	}
	else{
		prediction = predictFreeFlyNormal(p, t);
	}
	PBCpostX(&prediction.x);
	PBCpostY(&prediction.y);
	#if THREE_D
	PBCpostZ(&prediction.z);
	#endif
	return prediction;
}

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
	for (int i = 0; i < N; i++){
		renderPredictions[i] = graphicsParticlePrediction(particles + i);
	}
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

static bool refreshParticleColorCache(state* screenState){
	if ((screenState->colorParam2 == 0) || !ensureColorFunctionArray(screenState)){
		return false;
	}

	const double now = GetTime();
	double refreshInterval = COLOR_REFRESH_SECONDS;
	#if !THREE_D
	if ((screenState->colorParam2 == 3) || (screenState->colorParam2 == 4)){
		refreshInterval = BOOP_COLOR_REFRESH_SECONDS;
	}
	#endif
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

static void unloadSpriteBatch(spriteBatch* batch){
	if (batch->loaded){
		UnloadMesh(batch->mesh);
	}
	*batch = (spriteBatch){0};
}

static bool loadSpriteBatch(spriteBatch* batch, int capacity){
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
	for (int sprite = 0; sprite < capacity; sprite++){
		memcpy(mesh.texcoords + (size_t)sprite*12, quadUV, sizeof(quadUV));
	}
	memset(mesh.colors, 255, vertexCount*4*sizeof(unsigned char));

	UploadMesh(&mesh, true);
	batch->mesh = mesh;
	batch->capacity = capacity;
	batch->count = 0;
	batch->loaded = true;
	return true;
}

static bool ensureSpriteBatch(spriteBatch* batch, int capacity){
	if (batch->loaded && (batch->capacity == capacity)){
		return true;
	}
	unloadSpriteBatch(batch);
	return loadSpriteBatch(batch, capacity);
}

static void setBatchVertex(spriteBatch* batch, int vertex, float x, float y, Color color){
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

static bool appendBatchQuad(spriteBatch* batch, const Vector2 corners[4], Color color){
	if (!batch->loaded || (batch->count >= batch->capacity)){
		return false;
	}

	const int firstVertex = batch->count*6;
	static const int cornerOrder[6] = {0, 1, 2, 0, 2, 3};
	for (int i = 0; i < 6; i++){
		const Vector2 corner = corners[cornerOrder[i]];
		setBatchVertex(batch, firstVertex + i, corner.x, corner.y, color);
	}
	batch->count++;
	return true;
}

static void drawSpriteBatch(spriteBatch* batch, Texture2D texture){
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
	spriteMaterial.maps[MATERIAL_MAP_DIFFUSE].texture = texture;
	DrawMesh(visibleMesh, spriteMaterial, MatrixIdentity());
}

#if THREE_D
static const char* sphereVertexShader3D =
	"#version 330\n"
	"layout(location = 0) in vec3 vertexPosition;\n"
	"layout(location = 2) in vec3 vertexNormal;\n"
	"layout(location = 9) in mat4 instanceTransform;\n"
	"uniform mat4 mvp;\n"
	"out vec3 sphereNormal;\n"
	"out vec4 sphereColor;\n"
	"void main(void)\n"
	"{\n"
	"    mat3 scale = mat3(instanceTransform);\n"
	"    vec3 center = instanceTransform[3].xyz;\n"
	"    vec3 worldPosition = center + scale*vertexPosition;\n"
	"    sphereNormal = normalize(scale*vertexNormal);\n"
	"    sphereColor = vec4(instanceTransform[0].w,\n"
	"                       instanceTransform[1].w,\n"
	"                       instanceTransform[2].w,\n"
	"                       instanceTransform[3].w);\n"
	"    gl_Position = mvp*vec4(worldPosition, 1.0);\n"
	"}\n";

static const char* sphereFragmentShader3D =
	"#version 330\n"
	"in vec3 sphereNormal;\n"
	"in vec4 sphereColor;\n"
	"out vec4 finalColor;\n"
	"void main(void)\n"
	"{\n"
	"    vec3 normal = normalize(sphereNormal);\n"
	"    vec3 keyLight = normalize(vec3(-0.35, 0.80, 0.48));\n"
	"    vec3 fillLight = normalize(vec3(0.55, -0.20, -0.80));\n"
	"    float key = max(dot(normal, keyLight), 0.0);\n"
	"    float fill = max(dot(normal, fillLight), 0.0);\n"
	"    float brightness = min(1.0, 0.30 + 0.62*key + 0.12*fill);\n"
	"    finalColor = vec4(sphereColor.rgb*brightness, sphereColor.a);\n"
	"}\n";

static void sphereResolution3D(int particleCount, int* rings, int* slices){
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
	if ((loaded.shader.locs[SHADER_LOC_MATRIX_MVP] < 0) ||
		(loaded.shader.locs[SHADER_LOC_MATRIX_MODEL] < 0)){
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

static void drawSphereBatch(sphereBatch* batch){
	if (!batch->loaded || (batch->count == 0)){
		return;
	}
	DrawMeshInstanced(batch->mesh, batch->material,
		batch->transforms, batch->count);
}
#endif

static bool particleIntersectsViewport(float left, float top, float right, float bottom, Rectangle viewport){
	return !((right < viewport.x) || (left > viewport.x + viewport.width) ||
		(bottom < viewport.y) || (top > viewport.y + viewport.height));
}

static int structureFactorPopulation(void){
	int population = N - (int)(N*fractionSmallN);
	if (population < 0) population = 0;
	if (population > N) population = N;
	return population;
}

static int desiredStructureFactorSampleCount(void){
	const int population = structureFactorPopulation();
	return (population < STRUCTURE_FACTOR_MAX_SAMPLES) ? population : STRUCTURE_FACTOR_MAX_SAMPLES;
}

static bool allocateStructureFactorPositions(void){
	structSampleCount = desiredStructureFactorSampleCount();
	if (structSampleCount == 0){
		positions = NULL;
		return false;
	}
	positions = calloc((size_t)structSampleCount, sizeof(position));
	if (positions == NULL){
		structSampleCount = 0;
		TraceLog(LOG_ERROR, "Unable to allocate the structure-factor sample");
		return false;
	}
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
	pthread_mutex_lock(&structMutex);
	const bool running = structJobRunning;
	pthread_mutex_unlock(&structMutex);
	return running;
}

static void releaseStructureFactorPositions(void){
	awaitStructFactor();
	free(positions);
	positions = NULL;
	structSampleCount = 0;
}

void drawParticlesAndBox(double factor, state* screenState, float zoom, Rectangle viewport){
	particle* particleUnderClick = screenState->particleUnderClick;
	const bool predictionsReady = refreshRenderPredictions();
	const bool colorsReady = refreshParticleColorCache(screenState);
	const bool particleBatchReady = ensureSpriteBatch(&particleBatch, N);
	const bool arrowBatchReady = !screenState->arrow ||
		(particleBatchReady && ensureSpriteBatch(&arrowBatch, N));
	particleBatch.count = 0;
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
		if (particleBatchReady){
			appendBatchQuad(&particleBatch, particleCorners, screenState->particleColor);
		}
		else{
			DrawTexturePro(circleTexture,
				(Rectangle){0.0f, 0.0f, (float)circleTexture.width, (float)circleTexture.height},
				(Rectangle){left, top, diameter, diameter}, (Vector2){0, 0}, 0.0f,
				screenState->particleColor);
		}

		// Draw velocity arrow if enabled
		if (screenState->arrow) {
			double vx = predicted.vx;
			double vy = predicted.vy;
			double speed_squared = vx*vx + vy*vy;
			
			if (speed_squared > 1e-20) {
				const double inverseSpeed = 1.0/sqrt(speed_squared);
				const float sine = (float)(vx*inverseSpeed);
				const float cosine = (float)(vy*inverseSpeed);
				const float halfArrow = (float)(predicted.rad*factor*0.5);
				const float centerX = (float)(predicted.x*factor);
				const float centerY = (float)((Ly - predicted.y)*factor);
				const Vector2 localCorners[4] = {
					{-halfArrow, -halfArrow}, {-halfArrow, halfArrow},
					{halfArrow, halfArrow}, {halfArrow, -halfArrow},
				};
				Vector2 arrowCorners[4];
				for (int corner = 0; corner < 4; corner++){
					arrowCorners[corner].x = centerX + localCorners[corner].x*cosine - localCorners[corner].y*sine;
					arrowCorners[corner].y = centerY + localCorners[corner].x*sine + localCorners[corner].y*cosine;
				}
				if (arrowBatchReady){
					appendBatchQuad(&arrowBatch, arrowCorners, BLACK);
				}
			}
		}
	}

	if (particleBatchReady){
		drawSpriteBatch(&particleBatch, circleTexture);
	}
	if (screenState->arrow && arrowBatchReady){
		drawSpriteBatch(&arrowBatch, arrowTexture);
	}
	else if (screenState->arrow){
		for (int i = N - 1; i >= 0; i--){
			const particlePrediction predicted = predictionsReady ? renderPredictions[i] :
				graphicsParticlePrediction(particles + i);
			const double speedSquared = predicted.vx*predicted.vx + predicted.vy*predicted.vy;
			if (speedSquared <= 1e-20) continue;
			const float radius = (float)(predicted.rad*factor);
			const float centerX = (float)(predicted.x*factor);
			const float centerY = (float)((Ly - predicted.y)*factor);
			if (!particleIntersectsViewport(centerX - radius, centerY - radius,
				centerX + radius, centerY + radius, viewport)) continue;
			DrawTexturePro(arrowTexture,
				(Rectangle){0.0f, 0.0f, (float)arrowTexture.width, (float)arrowTexture.height},
				(Rectangle){centerX, centerY, radius, radius},
				(Vector2){radius*0.5f, radius*0.5f},
				(float)(atan2(predicted.vx, predicted.vy)*180.0/M_PI), BLACK);
		}
	}
	if (screenState->selected){
		const particlePrediction selected = cachedGraphicsParticlePrediction(particleUnderClick);
		DrawRing((Vector2){selected.x*factor, (Ly - selected.y)*factor},
			0.7*selected.rad*factor, selected.rad*factor, 0, 360, 20, BLACK);
	}
	if (addWallx){
		DrawLineEx((Vector2){0, 0}, (Vector2){0, Ly*factor} , 0.3*factor/zoom, BLACK);
		DrawLineEx((Vector2){Lx*factor, 0}, (Vector2){Lx*factor, Ly*factor} , 0.3*factor/zoom, BLACK);
	}
	if (addWally){
		DrawLineEx((Vector2){0, 0}, (Vector2){Lx*factor, 0} , 0.3*factor/zoom, BLACK);
		DrawLineEx((Vector2){0, Ly*factor}, (Vector2){Lx*factor, Ly*factor} , 0.3*factor/zoom, BLACK);
	}
	#if !THREE_D
	else if (addCircularWall){
		DrawRing((Vector2){Lx/2*factor, Lx/2*factor}, Lx/2*factor, (Lx/2 + 1)*factor, 0, 360, 360, BLACK);
	}
	#endif
}

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
		else{
			const float x = (float)(u*Lx);
			const float z = (float)(u*Lz);
			DrawLine3D((Vector3){x, coordinate, 0},
				(Vector3){x, coordinate, (float)Lz}, color);
			DrawLine3D((Vector3){0, coordinate, z},
				(Vector3){(float)Lx, coordinate, z}, color);
		}
	}
}

static void drawParticlesAndBox3D(state* screenState){
	const bool predictionsReady = refreshRenderPredictions();
	const bool colorsReady = refreshParticleColorCache(screenState);
	const bool sphereBatchReady = ensureSphereBatch(&particleSphereBatch, N);
	particleSphereBatch.count = 0;
	drawBoxBoundary3D((Color){105, 115, 125, 255});
	if (addWallx){
		const Color wallColor = Fade(BLACK, 0.45f);
		drawWallGrid3D(0, 0.0f, wallColor);
		drawWallGrid3D(0, (float)Lx, wallColor);
	}
	if (addWally){
		const Color wallColor = Fade(BLACK, 0.45f);
		drawWallGrid3D(1, 0.0f, wallColor);
		drawWallGrid3D(1, (float)Ly, wallColor);
	}

	const float axisLength = (float)(0.12*fmin(Lx, fmin(Ly, Lz)));
	DrawLine3D((Vector3){0}, (Vector3){axisLength, 0, 0}, RED);
	DrawLine3D((Vector3){0}, (Vector3){0, axisLength, 0}, GREEN);
	DrawLine3D((Vector3){0}, (Vector3){0, 0, axisLength}, BLUE);

	for (int i = N - 1; i >= 0; i--){
		const particlePrediction predicted = predictionsReady ? renderPredictions[i] :
			graphicsParticlePrediction(particles + i);
		Color particleColor;
		if (colorsReady){
			particleColor = colorSelect(
				(screenState->colorFunctionArray[i] - colorCache.min)/
					(colorCache.max - colorCache.min), screenState->colorArray);
		}
		else{
			particleColor = screenState->particleArray[particles[i].type];
		}
		const Vector3 center = {
			(float)predicted.x, (float)predicted.y, (float)predicted.z,
		};
		if (sphereBatchReady){
			appendSphereBatch(&particleSphereBatch, center,
				(float)predicted.rad, particleColor);
		}
		else{
			DrawSphereEx(center, (float)predicted.rad, 4, 8, particleColor);
		}

		if (screenState->arrow){
			const double speedSquared = predicted.vx*predicted.vx +
				predicted.vy*predicted.vy + predicted.vz*predicted.vz;
			if (speedSquared > 1e-20){
				const float arrowLength = (float)(2.5*predicted.rad);
				const float inverseSpeed = (float)(1.0/sqrt(speedSquared));
				const Vector3 end = {
					center.x + arrowLength*(float)predicted.vx*inverseSpeed,
					center.y + arrowLength*(float)predicted.vy*inverseSpeed,
					center.z + arrowLength*(float)predicted.vz*inverseSpeed,
				};
				DrawLine3D(center, end, BLACK);
			}
		}
	}
	if (sphereBatchReady){
		drawSphereBatch(&particleSphereBatch);
	}

	if (screenState->selected && (screenState->particleUnderClick != NULL)){
		const particlePrediction selected = cachedGraphicsParticlePrediction(
			screenState->particleUnderClick);
		const Vector3 center = {
			(float)selected.x, (float)selected.y, (float)selected.z,
		};
		DrawSphereWires(center, (float)(1.08*selected.rad), 8, 12, BLACK);
	}
}
#endif

void draw(int argc, char *argv[], window* screenWindow, state* screenState){
    int start = screenWindow->start;
	double xGUI = screenWindow->xGUI;
    double yGUI = screenWindow->yGUI;
	#if !THREE_D
    double factor = screenWindow->factor;
	#endif
	particle* particleUnderClick = screenState->particleUnderClick;
	char name[200];

	#if THREE_D
	ensureCamera3D(screenWindow);
	int scenePixelWidth = 0;
	int scenePixelHeight = 0;
	sceneTargetSize3D(start, GetScreenHeight(),
		&scenePixelWidth, &scenePixelHeight);
	const bool sceneReady = ensureSceneTarget3D(scenePixelWidth,
		scenePixelHeight);
	if (sceneReady){
		BeginTextureMode(sceneTarget);
		ClearBackground((Color){245, 247, 250, 255});
		BeginMode3D(screenWindow->cam3D);
		drawParticlesAndBox3D(screenState);
		EndMode3D();
		EndTextureMode();
	}
	#endif

	BeginDrawing();

	ClearBackground(RAYWHITE);
	#if THREE_D
	if (sceneReady){
		DrawTexturePro(sceneTarget.texture,
			(Rectangle){0, 0, (float)sceneTargetWidth, (float)-sceneTargetHeight},
			(Rectangle){0, 0, (float)start, (float)GetScreenHeight()},
			(Vector2){0, 0}, 0.0f, WHITE);
	}
	else{
		DrawText("Unable to create the 3D viewport", 20, 20, 20, MAROON);
	}
	DrawText("Right drag: orbit   Wheel: zoom   R: reset view",
		10, 10, 18, Fade(BLACK, 0.7f));
	#else
	const Vector2 viewportCornerA = GetScreenToWorld2D((Vector2){0, 0}, screenWindow->cam);
	const Vector2 viewportCornerB = GetScreenToWorld2D(
		(Vector2){(float)screenWindow->start, (float)GetScreenHeight()}, screenWindow->cam);
	const Rectangle particleViewport = {
		.x = fminf(viewportCornerA.x, viewportCornerB.x),
		.y = fminf(viewportCornerA.y, viewportCornerB.y),
		.width = fabsf(viewportCornerB.x - viewportCornerA.x),
		.height = fabsf(viewportCornerB.y - viewportCornerA.y),
	};

	BeginMode2D(screenWindow->cam);
	drawParticlesAndBox(factor, screenState, screenWindow->cam.zoom, particleViewport);
	
	EndMode2D();
	#endif

	DrawRectangle(start, 0, GetScreenWidth() - start, GetScreenHeight(), Fade((Color){220, 220, 220, 255}, 1.f));

	if (screenState->selected && (particleUnderClick != NULL)){
		freeFly(particleUnderClick);
		int changes = 0;
		changes += doubleBox(&(particleUnderClick->vx), "vx", &screenState->editVx, (Rectangle){ start + 40*xGUI, 800*yGUI, 120*xGUI, 24*yGUI});
		changes += doubleBox(&(particleUnderClick->vy), "vy", &screenState->editVy, (Rectangle){ start + 40*xGUI, 824*yGUI, 120*xGUI, 24*yGUI});
		#if THREE_D
		changes += doubleBox(&(particleUnderClick->vz), "vz", &screenState->editVz, (Rectangle){ start + 40*xGUI, 848*yGUI, 120*xGUI, 24*yGUI});
		changes += doubleBox(&(particleUnderClick->m), "m", &screenState->editM, (Rectangle){ start + 40*xGUI, 872*yGUI, 120*xGUI, 24*yGUI});
		#else
		changes += doubleBox(&(particleUnderClick->m), "m", &screenState->editM, (Rectangle){ start + 40*xGUI, 848*yGUI, 120*xGUI, 24*yGUI});
		#endif
		if (changes){
			invalidateParticleColorCache();
			particleUnderClick->coll++;

			removeEventFromQueue(eventList[N + particleUnderClick->num]);
			removeEventFromQueue(eventList[particleUnderClick->num]);
			crossingEvent(particleUnderClick->num);
			collisionEvent(particleUnderClick->num);
		}
	}

	if (screenState->colorEditing2) 
		GuiLock();

	int dirtyColorParam2 = screenState->colorParam2;
	#if THREE_D
	const char* colorModeNames = "Uniform; Coll. Based; Vel. Based; Radius Based; Charge Based";
	#else
	const char* colorModeNames = "Uniform; Coll. Based; Vel. Based; Hex. Based; Square. Based; Radius Based; Charge Based";
	#endif
	if (GuiDropdownBox((Rectangle){start  + 100*xGUI, 755*yGUI, 120*xGUI, 15*yGUI }, colorModeNames, &screenState->colorParam2, screenState->colorEditing2)){
		
		screenState->colorEditing2 = !screenState->colorEditing2;
		if (screenState->colorParam2 == 0){
			screenState->particleColor = MAROON;
		}
		else if (screenState->colorParam2 == 1){
			screenState->colorFunction = &colorCollision;
		}
		else if (screenState->colorParam2 == 2){
			screenState->colorFunction = &colorVelocity;
		}
		#if THREE_D
		else if (screenState->colorParam2 == 3){
			screenState->colorFunction = &colorRadius;
		}
		else if (screenState->colorParam2 == 4){
			screenState->colorFunction = &colorCharge;
		}
		#else
		else if (screenState->colorParam2 == 3){
			nSym = 6;
			screenState->colorFunction = &colorBOOP;
		}
		else if (screenState->colorParam2 == 4){
			nSym = 4;
			screenState->colorFunction = &colorBOOP;
		}
		else if (screenState->colorParam2 == 5){
			nSym = 4;
			screenState->colorFunction = &colorRadius;
		}
		else if (screenState->colorParam2 == 6){
			screenState->colorFunction = &colorCharge;
		}
		#endif
		if (dirtyColorParam2 != screenState->colorParam2){
			invalidateParticleColorCache();
		}
		if (screenState->colorParam2 != 0){
			ensureColorFunctionArray(screenState);
		}
		else if (dirtyColorParam2 != 0){
			free(screenState->colorFunctionArray);
			screenState->colorFunctionArray = NULL;
			screenState->colorFunctionCapacity = 0;
		}
	}

	bool tempPoly = polydispersity;
	GuiCheckBox((Rectangle){start + 100*xGUI, 715*yGUI, 40*xGUI, 40*yGUI}, "Polydisp.", &polydispersity);
	if (polydispersity != tempPoly){
		if (polydispersity){
			//small hack to get the structure factor right!
			fractionSmallN = 0;
		}
		if (structFactorActivated)
			releaseStructureFactorPositions();
		freeArrays();

		reset(argc, argv, &screenWindow->factor, screenState);
		
	}

	GuiCheckBox((Rectangle){start + 200*xGUI, 715*yGUI, 40*xGUI, 40*yGUI}, "Arrow", &screenState->arrow);

	if (screenState->leftClicked == 0)
		GuiUnlock();

	if (screenState->colorParam2){
		if (screenState->colorEditing) 
			GuiLock();

		if (GuiDropdownBox((Rectangle){start  + 220*xGUI, 755*yGUI, 120*xGUI, 15*yGUI }, "Plasma;Viridis;Copper", &screenState->colorParam, screenState->colorEditing)){
			screenState->colorEditing = !screenState->colorEditing;
			if (screenState->colorParam == 0){
				screenState->colorArray = Plasma;
			}
			else if (screenState->colorParam == 1){
				screenState->colorArray = Viridis;
			}
			else{
				screenState->colorArray = Copper;
			}
		}	

		if (screenState->leftClicked == 0)
			GuiUnlock();
	}		
	
	const bool oldDamping = damping;
	bool requestedDamping = damping;
	GuiCheckBox((Rectangle){start + 406*xGUI, 40*yGUI, 40*xGUI, 40*yGUI},
		"Damping", &requestedDamping);
	if (oldDamping != requestedDamping){
		/* Commit every particle with the old law before changing event semantics. */
		for(int i = 0; i < N; i++){
			freeFly(particles + i);
		}
		/* Remove the complete old-law schedule before rebuilding any new-law event. */
		for(int i = 0; i < N; i++){
			particle* p = particles + i;
			removeEventFromQueue(eventList[N + p->num]);
			removeEventFromQueue(eventList[p->num]);
		}
		#if INTRUDER
		removeEventFromQueue(eventList[2*N + 8]);
		#endif
		damping = requestedDamping;
		for(int i = 0; i < N; i++){
			particle* p = particles + i;
			crossingEvent(p->num);
			collisionEvent(p->num);
		}
		#if INTRUDER
		if (freeFly == &freeFlyNormal){
			collisionEventIntruder();
		}
		#endif
		invalidateParticleColorCache();
	}
	else{
		damping = requestedDamping;
	}

	


	int dirtyNoise = noise;
	GuiToggleGroup((Rectangle){start  + 100*xGUI, 40*yGUI, 100*xGUI, 40*yGUI}, "No Thermostat;Langevin;Vel. Rescale", &noise); 
	
	if (dirtyNoise != noise){

		if (noise == 0){
			removeEventFromQueue(eventList[2*N + 2]);
		}
		if (dirtyNoise == 0){
			addEventNoise(t + dtnoise);
		}
	}	

	if (noise != 0 || damping != 0){
		const double oldGamma = gamm;
		sprintf(name, "%.3f", gamm);
		GuiSliderBarDouble(
			(Rectangle){ start + 100*xGUI, 130*yGUI, 505*xGUI, 40*yGUI },
			"Gamma", name, &gamm, 0.001f, 0.1f);
		if (oldGamma != gamm){
			if (damping){
				for(int i = 0; i < N; i++){
					particle* p = particles + i;
					freeFly(p);
					removeEventFromQueue(eventList[N + p->num]);
					removeEventFromQueue(eventList[p->num]);
					crossingEvent(p->num);
					collisionEvent(p->num);
				}
			}
			invalidateParticleColorCache();
		}
		if (noise == 1){
			sprintf(name, "%.3f", T);
			GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 80*yGUI, 505*xGUI, 40*yGUI}, "Temperature", name, &T, 0.0000000001f, 1.f);
		}
	}
	

	bool dirtyWell = addWell;
	GuiCheckBox((Rectangle){ start  + 100*xGUI, 180*yGUI, 40*xGUI, 40*yGUI}, "Potential", &addWell);
	if (dirtyWell != addWell){
		invalidateParticleColorCache();
		for(int i = 0; i < N; i++){
			particle* p = particles + i;
			freeFly(p);
			if (addWell == 0){
				free(p->particlesInWell);
			}
			
		}
		
		free(cellList);
		boxConstantHelper();
		cellListInit();
		if (addWell)
			pairsInit();
		
		for (int i = 0; i < N; i++){
			particle* p = particles + i;
			removeEventFromQueue(eventList[N + p->num]);
			removeEventFromQueue(eventList[p->num]);
			crossingEvent(p->num);
			collisionEvent(p->num);
		}
	}
	
	if (addWell){
		const double oldPotential = U;
		sprintf(name, "%.6f", U);
		GuiSliderBarDouble((Rectangle){ start + 100*xGUI, 220*yGUI, 505*xGUI, 40*yGUI }, "U", name, &U, -1.f, 1.f);
		if (U != oldPotential){
			invalidateParticleColorCache();
		}
		const double oldPotentialRadius = sig;
		sprintf(name, "%.3f", sig);
		GuiSliderBarDouble(
			(Rectangle){ start  + 100*xGUI, 270*yGUI, 505*xGUI, 40*yGUI},
			"Pot. rad.", name, &sig, 1.01f, 2.5f);
		if (sig != oldPotentialRadius){
			free(cellList);
			boxConstantHelper();
			cellListInit();
			for(int i = 0; i < N; i++){
				particle *p = particles + i;
				freeFly(p);
				free(p->particlesInWell);
			} 
			pairsInit();
			
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				removeEventFromQueue(eventList[N + p->num]);
				removeEventFromQueue(eventList[p->num]);
				crossingEvent(p->num);
				collisionEvent(p->num);
				

			}
			invalidateParticleColorCache();
		}
		
		GuiCheckBox((Rectangle){ start  + 200*xGUI, 180*yGUI, 40*xGUI, 40*yGUI}, "Charged", &charged);
		if (charged){
			const double oldPositiveChargeFraction = proportionPositivelyCharged;
			sprintf(name, "%.2f", proportionPositivelyCharged);
			GuiSliderBarDouble(
				(Rectangle){ start + 400*xGUI, 180*yGUI, 50*xGUI, 40*yGUI },
				"%+", name, &proportionPositivelyCharged, 0.f, 0.5f);
			const double oldNeutralChargeFraction = proportionNeutralyCharged;
			sprintf(name, "%.2f", proportionNeutralyCharged);
			GuiSliderBarDouble(
				(Rectangle){ start + 500*xGUI, 180*yGUI, 50*xGUI, 40*yGUI },
				"%n", name, &proportionNeutralyCharged, 0.f, 0.5f);
			if ((oldPositiveChargeFraction != proportionPositivelyCharged) ||
				(oldNeutralChargeFraction != proportionNeutralyCharged)){
				for (int i = 0; i < N; i++){
					double randomTemp = drand(0, 1);
					if (randomTemp < proportionPositivelyCharged){
						particles[i].charge = 1;
					}
					else if (randomTemp < proportionNeutralyCharged + proportionPositivelyCharged){
						particles[i].charge = 0;
					}
					else{
						particles[i].charge = -1;
					}
				}
				invalidateParticleColorCache();
			}
		}
	}

	bool dirtyField = addField;
	GuiCheckBox((Rectangle){ start  + 100*xGUI, 320*yGUI, 40*xGUI, 40*yGUI}, "Const. Field", &addField);
	if (dirtyField != addField){
		for(int i = 0; i < N; i++){			
			particle* p = particles + i;
			freeFly(p);
			removeEventFromQueue(eventList[N + p->num]);
			removeEventFromQueue(eventList[p->num]);
			crossingEvent(p->num);
			collisionEvent(p->num);
		}
		invalidateParticleColorCache();
	}

	GuiCheckBox((Rectangle){start + 250*xGUI, 320*yGUI, 40*xGUI, 40*yGUI}, "Delta", &addDelta);
	#if !THREE_D
	GuiCheckBox((Rectangle){start + 400*xGUI, 320*yGUI, 40*xGUI, 40*yGUI}, "Tangent Delta", &addDeltaTangent);
	#endif
	if (addDelta || addDeltaTangent){
		sprintf(name, "%.3f", delta);
		GuiSliderBarDouble((Rectangle){ start + 100*xGUI, 360*yGUI, 505*xGUI, 20*yGUI}, "Delta", name, &delta, 0, 1);
	}

	if (addField){
		const double oldField = field;
		sprintf(name, "%.3f", field);
		GuiSliderBarDouble(
			(Rectangle){ start + 100*xGUI, 380*yGUI, 505*xGUI, 20*yGUI},
			"g", name, &field, -0.03, -0.0001);
		if (oldField != field){
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				freeFly(p);
				removeEventFromQueue(eventList[N + p->num]);
				removeEventFromQueue(eventList[p->num]);
				crossingEvent(p->num);
				collisionEvent(p->num);

			}
			invalidateParticleColorCache();
		}
	}
	
	sprintf(name, "%.3f", res);
	GuiSliderBarDouble((Rectangle){start + 100*xGUI, 580*yGUI, 200*xGUI, 40*yGUI}, "Coeff of res.", name, &res, 0.f, 1.f);

	float requestedParticleCount = (float)N;
	sprintf(name, "%d", N);
	GuiSliderBar((Rectangle){start + 100*xGUI, 615*yGUI, 200*xGUI, 25*yGUI},
		"N. of particles", name, &requestedParticleCount, 50.f, 50000.f);
	const int newParticleCount = (int)requestedParticleCount;
	if (newParticleCount != N){
		if (structFactorActivated)
			releaseStructureFactorPositions();
		freeArrays();
		free(screenState->colorFunctionArray);
		screenState->colorFunctionArray = NULL;
		screenState->colorFunctionCapacity = 0;
		N = newParticleCount;
		reset(argc, argv, &screenWindow->factor, screenState);
		if (screenState->colorParam2)
			ensureColorFunctionArray(screenState);
	}

	const double oldSizeRatio = sizeratio;
	const double oldSmallFraction = fractionSmallN;
	if (!polydispersity){
		
		sprintf(name, "%.3lf", sizeratio);
		GuiSliderBarDouble(
			(Rectangle){start + 100*xGUI, 640*yGUI, 200*xGUI, 25*yGUI},
			"Sizeratio", name, &sizeratio, 0.25, 1);


		sprintf(name, "%.3lf", fractionSmallN);
		GuiSliderBarDouble(
			(Rectangle){start + 100*xGUI, 665*yGUI, 200*xGUI, 25*yGUI},
			"Frac. Of Small.", name, &fractionSmallN, 0, 1);
	}
	

	const double oldPackingFraction = phi;
	sprintf(name, "%.3lf", phi);
	GuiSliderBarDouble(
		(Rectangle){start + 100*xGUI, 690*yGUI, 200*xGUI, 25*yGUI},
		"Packing fraction", name, &phi, 0.1, 0.88);
	if ((oldPackingFraction != phi) || (oldSmallFraction != fractionSmallN) ||
		(oldSizeRatio != sizeratio)){
	/*
		if ((phiTemp > phi) && (t > 1/vr)){
			double mult = sqrt(phiTemp/phi);
			Lx = Lx*mult;
			Ly = Lx;
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				freeFly(p);
				p->x = p->x*mult;
				p->y = p->y*mult;
			}
			printf("Lx = %lf\n", Lx);
			
			free(cellList);
			boxConstantHelper();
			cellListInit();
			factor = screenHeight/Ly;
			if (addWell){
				for (int i = 0; i < N; i++){
					particle* p = particles + i;
					free(particles[i].particlesInWell);
				}
				pairsInit();
			}
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				p->coll++;
				removeEventFromQueue(eventList[N + p->num]);
				removeEventFromQueue(eventList[p->num]);
				crossingEvent(p->num);
				collisionEvent(p->num);
			}
		}
		else{
		*/
		if (structFactorActivated)
			releaseStructureFactorPositions();
		freeArrays();

		reset(argc, argv, &screenWindow->factor, screenState);
        
		//}
	}
	



	int dirtyWallParam = screenState->wallParam;
	#if THREE_D
	const char* wallModeNames = "No Wall; Horiz. Wall; Vert. Wall; Square Walls";
	#else
	const char* wallModeNames = "No Wall; Horiz. Wall; Vert. Wall; Square Walls; Circl. Wall";
	#endif
	if (GuiDropdownBox((Rectangle){start + 100*xGUI, 400*yGUI, 120*xGUI, 24*yGUI }, wallModeNames, &screenState->wallParam, screenState->wallEditing)){
		screenState->wallEditing = !screenState->wallEditing;
		if (dirtyWallParam != screenState->wallParam){
			invalidateParticleColorCache();
			if (screenState->wallParam == 0){
				addWallx = 0;
				addWally = 0;
				addCircularWall = 0;
			}
			else if (screenState->wallParam == 1){
				addWallx = 0;
				addWally = 1;
				addCircularWall = 0;
			}
			else if (screenState->wallParam == 2){
				addWallx = 1;
				addWally = 0;
				addCircularWall = 0;
			}
			else if (screenState->wallParam == 3){
				addWallx = 1;
				addWally = 1;
				addCircularWall = 0;
			}
			#if !THREE_D
			else if (screenState->wallParam == 4){
				addWallx = 0;
				addWally = 0;
				addCircularWall = 1;
			}
			#endif
			if (screenState->wallParam != 0){
				
				if (structFactorActivated)
					releaseStructureFactorPositions();
				freeArrays();
				reset(argc, argv, &screenWindow->factor, screenState);
			}
			else{
				if (addWell){
					//To take care of the boundary conditions
					for (int i = 0; i < N; i++){
						particle* p = particles + i;
						freeFly(p);
						free(particles[i].particlesInWell);
					}
					pairsInit();
				}
				for (int i = 0; i < N; i++){
					particle* p = particles + i;
					freeFly(p);
					p->coll++;
					removeEventFromQueue(eventList[N + p->num]);
					removeEventFromQueue(eventList[p->num]);
					crossingEvent(p->num);
					collisionEvent(p->num);
				}
			}
		}
	}	

	
	bool tempStruct = structFactorActivated;
	#if THREE_D
	const char* structureFactorLabel = "XY Struct. Factor";
	#else
	const char* structureFactorLabel = "Struct. Factor";
	#endif
	GuiCheckBox((Rectangle){start + 500*xGUI, 40*yGUI, 40*xGUI, 40*yGUI}, structureFactorLabel, &structFactorActivated);
	if (tempStruct != structFactorActivated){
		if (structFactorActivated){
			if (allocateStructureFactorPositions()){
				captureStructureFactorPositions();
				asyncStructFactor();
				lastStructLaunchTime = GetTime();
			}
			else{
				structFactorActivated = false;
			}
		}
		else{
			releaseStructureFactorPositions();
		}
	}
	if (structFactorActivated){
		if (finishedStructComputation()){
			awaitStructFactor();
			normalizeStruct();
			UpdateTexture(texture, structFactorGraphics);
			lastStructLaunchTime = GetTime();
		}
		if (screenState->running && !structureFactorJobIsRunning() &&
			(GetTime() - lastStructLaunchTime >= STRUCTURE_FACTOR_REFRESH_SECONDS)){
			captureStructureFactorPositions();
			asyncStructFactor();
			lastStructLaunchTime = GetTime();
		}
		DrawTexturePro(texture,
            (Rectangle){ 0, 0, qN, qN },
            (Rectangle){ 400 + start*xGUI, 450*yGUI, SIZE*xGUI, SIZE*yGUI},
            (Vector2){ 0, 0 }, 0.0f, (Color){255, 255, 255, 255});
	}


	DrawFPS(GetScreenWidth() - 100, 10);
	
	EndDrawing();
}

void saveSnapshotGraphics(){
	for (int i = 0; i < N; i++){
		freeFly(particles + i);
	}
	customName();
	fichier = fopen(fileName, "w");
	saveTXT();
	fclose(fichier);
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
	int count = 0;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	double re = 0;
	double im = 0;
	const particlePrediction p1Now = cachedGraphicsParticlePrediction(p1);
	for (int j = -1; j <= 1; j++){
		for (int k = -1; k <= 1; k++){
			particle* p2 =  cellList[PBCcellY(Y + j) * Nxcells + PBCcellX(X + k)];
			while (p2 != NULL){
				if ((p1->num != p2->num) && (p1->type == 1) && (p2->type == 1)){
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

void* computeStructureFactorGraphic(void* arg){
	threadArg* threadArgument = (threadArg*)arg;
	unsigned long completedGeneration = 0;
	const double qStep = qxG[1] - qxG[0];

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

		for (int i = threadArgument->start; i < threadArgument->end; i++){
			double rowReal[qN] = {0};
			double rowImaginary[qN] = {0};
			for (int n = 0; n < structSampleCount; n++){
				const position* p = positions + n;
				const double firstPhase = qxG[i]*p->x + qxG[0]*p->y;
				const double phaseStep = qStep*p->y;
				double phaseReal = cos(firstPhase);
				double phaseImaginary = sin(firstPhase);
				const double stepReal = cos(phaseStep);
				const double stepImaginary = sin(phaseStep);

				for (int j = 0; j < qN; j++){
					rowReal[j] += phaseReal;
					rowImaginary[j] += phaseImaginary;
					const double nextReal = phaseReal*stepReal - phaseImaginary*stepImaginary;
					phaseImaginary = phaseReal*stepImaginary + phaseImaginary*stepReal;
					phaseReal = nextReal;
				}
			}

			for (int j = 0; j < qN; j++){
				const double magnitude = (rowReal[j]*rowReal[j] +
					rowImaginary[j]*rowImaginary[j])/structSampleCount;
				structFactorGraphics[i][j] = (float)log(magnitude*magnitude + 1.0);
				if ((i <= qN/2 + 3) && (i >= qN/2 - 3) &&
					(j >= qN/2 - 3) && (j <= qN/2 + 3)){
					structFactorGraphics[i][j] = 0;
				}
			}
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
	awaitStructFactor();
	pthread_mutex_lock(&structMutex);
	structThreadPoolStopping = true;
	pthread_cond_broadcast(&structJobReady);
	pthread_mutex_unlock(&structMutex);
	for (int i = 0; i < NUM_THREADS; i++){
		pthread_join(threads[i], NULL);
	}
	structThreadPoolStarted = false;
}

void reset(int argc, char *argv[], double* factor, state* screenState){
	invalidateParticleColorCache();
	t = 0;
	ncol = 0;
	ncross = 0;
	paulTime = 0;
	actualPaulList = 0;
	screenState->selected = 0;
	screenState->leftClicked = 0;
	screenState->particleUnderClick = NULL;
	version = 1;

	
	constantInit(argc, argv);
	particlesInit();
	cellListInit();
	//stupid allocation to be able to change sig and sigtemp in the graphical mode even though we don't need it yet!
	if (addWell)
		pairsInit();
	
	eventListInit();
	physicalQ();
	format();
	if (structFactorActivated){
		if (!allocateStructureFactorPositions()){
			structFactorActivated = false;
		}
	}
	
	*factor = GetScreenHeight()/Ly;
}

int doubleBox(double* ptr, char* text, bool* activate, Rectangle position){
	char textBuffer[128] = {0};
	char textBuffer2[128] = {0};
	sprintf(textBuffer, "%lf", *ptr);
	sprintf(textBuffer2, "%lf", *ptr);
	if (GuiTextBox(position, textBuffer, sizeof(textBuffer), *activate)){
		*activate = !(*activate);
	}
	Rectangle textBounds = {0};
	textBounds.width = (float)GetTextWidth(text);
	textBounds.height = (float)GuiGetStyle(DEFAULT, TEXT_SIZE);
	textBounds.x = position.x - textBounds.width - GuiGetStyle(SLIDER, TEXT_PADDING);
	textBounds.y = position.y + position.height/2 - GuiGetStyle(DEFAULT, TEXT_SIZE)/2;

	GuiDrawText(text, textBounds, TEXT_ALIGN_RIGHT, Fade(GetColor(GuiGetStyle(SLIDER, TEXT + (STATE_NORMAL*3))), guiAlpha));

	if (strcmp(textBuffer, textBuffer2) == 0){
		return 0;
	}
	int conversion = sscanf(textBuffer, "%lf", ptr);
	//awful dumb trick because sscanf "says" that the conversion is successful even when the number is -0

	if (*ptr == -0.0){
			//for some reason *ptr = 0 doesn't work :) the compiler optimize it out with -Ofast
            *ptr = 0.000000000000000001;
    }

	if (conversion == 1){
		return 1;
	}
	else{
		sscanf(textBuffer2, "%lf", ptr);
		return 0;
	}
}

window graphicalInit(){

	firstScreen = 0;
	load = 0;
	dtnoise = 1;
	dtime = 1;
	N = 10000;
	phi = 0.5;
	field = -0.001;
	res = 0.6;
	delta = 1;
	gamm = 0.025;
	sig = 1.5;
	U = 0.15;
	T = 0.01;
	Einit = 0.01;

	load = 0;
	S1 = 0;
	Hex = 0;
	coexistenceOld = 0;
	coexistence = 0;

	interfaceThermo = 0;
	interfaceThermoTemp = 0;

	clusterThermo = 0;
	dumpCluster = 0;
	killCluster = 0;
	strucThermo = 0;
	areaThermo = 0;
	boopThermo = 0;
	pcfThermo = 0;
	pcfBondOrderThermo = 0;
	pcfg6Thermo = 0;
	critical = 0;
	snapshotCritical = 0;

	aspectRatio = 1;
	#if THREE_D
	/* These controls have only planar implementations. */
	addCircularWall = 0;
	addDeltaTangent = false;
	res = 1;
	damping = false;
	#endif
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
	SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_WINDOW_HIGHDPI);
	
	#if THREE_D
	InitWindow(1800, 900, "EDMD - 3D system");
	#else
	InitWindow(1800, 900, "EDMD");
	#endif
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
		const Vector2 monitorPosition = GetMonitorPosition(display);
		SetWindowPosition((int)monitorPosition.x + (monitorWidth - initialWidth)/2,
			(int)monitorPosition.y + (monitorHeight - initialHeight)/2);
	}
	SetTargetFPS(144);
	image = (Image){
		.data = structFactorGraphics,
		.width = qN,
		.height = qN,
		.mipmaps = 1,
		.format = PIXELFORMAT_UNCOMPRESSED_R32,
	};
	texture = LoadTextureFromImage(image);
	
	// Load circle texture from embedded data
	Image circleImage = LoadImageFromMemory(".png", src_tex_png, src_tex_png_len);
	circleTexture = LoadTextureFromImage(circleImage);
	UnloadImage(circleImage);  

	Image arrowImage = LoadImage("src/arrow.png");
	arrowTexture = LoadTextureFromImage(arrowImage);
	UnloadImage(arrowImage);  

	spriteMaterial = LoadMaterialDefault();
	defaultMaterialTexture = spriteMaterial.maps[MATERIAL_MAP_DIFFUSE].texture;
	spriteMaterialLoaded = true;

	const int initialScreenWidth = GetScreenWidth();
	const int initialScreenHeight = GetScreenHeight();
    return (window){
		.screenWidth = initialScreenWidth,
		.screenHeight = initialScreenHeight,
		.factor = initialScreenHeight/Ly,
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
		.particleArray = {GRAY, MAROON},
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
	if (sceneTarget.texture.id != 0){
		UnloadRenderTexture(sceneTarget);
		sceneTarget = (RenderTexture2D){0};
		sceneTargetWidth = 0;
		sceneTargetHeight = 0;
	}
	#endif
	unloadSpriteBatch(&particleBatch);
	unloadSpriteBatch(&arrowBatch);
	free(renderPredictions);
	renderPredictions = NULL;
	renderPredictionCapacity = 0;
	if (spriteMaterialLoaded){
		spriteMaterial.maps[MATERIAL_MAP_DIFFUSE].texture = defaultMaterialTexture;
		UnloadMaterial(spriteMaterial);
		spriteMaterialLoaded = false;
	}
	UnloadTexture(circleTexture);
	UnloadTexture(arrowTexture);
	UnloadTexture(texture);
	free(screenState->colorFunctionArray);
	screenState->colorFunctionArray = NULL;
	screenState->colorFunctionCapacity = 0;
	invalidateParticleColorCache();
}
