#include "pcf.h"
#include "EDMD.h"
#include "mersenne.c" 
#include "quartic.c" 
#include "interface.c"
#include "critical.c"
#include "struc.c"
#include "cluster.c"
#include "boop.h"
#include "jc_voronoi.h"
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<stdio.h>
#include<time.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<stdbool.h>
#include<float.h>
#include<limits.h>
#include<stdint.h>
#include<errno.h>
#include<stdarg.h>
#include"parser.h"
#include<getopt.h>
#include<time.h>
#include "voronoi_edmd.h"
#include "osmosis.h"
#include "field.h"

#if POLYGON_PARTICLES && !THREE_D
#include "polygon_particles.h"
#endif

#include "intruder.h"



#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#ifndef NONADDITIVE
#define NONADDITIVE 1
#endif

#if NONADDITIVE == 0
    #define DIST(r1, r2) (sqrt(4 * (r1) * (r2)))
	#define DIST2(r1, r2) ((4 * (r1) * (r2)))
#else
    #define DIST(r1, r2) ((r1) + (r2))
	#define DIST2(r1, r2) (((r1) + (r2))*((r1) + (r2)))
#endif

#ifndef G
#define G 0
#endif
#if G
#include "graphics.h"
#include "gui_imgui.h"
#include "live_analysis.h"
#include "raylib.h"
#if !defined(PLATFORM_WEB)
#include <pthread.h>
#else
#include <emscripten/emscripten.h>
#endif
#endif

#ifndef STRESS
#define STRESS 0
#endif
#ifndef HEAT
#define HEAT 0
#endif

#define CGREEN printf("\033[1;32m")
#define CWHITE printf("\033[0;37m")
#define CCYAN printf("\033[1;96m")

#define ERROR 0


#ifndef PAUL_DT_SCALE
#define PAUL_DT_SCALE 5.0
#endif

#ifndef EDMD_VERIFY_QUEUE
#define EDMD_VERIFY_QUEUE 0
#endif

#ifndef EDMD_VERIFY_OVERLAPS
#define EDMD_VERIFY_OVERLAPS 0
#endif

#ifndef EDMD_VERIFY_COLLISION_ACCOUNTING
#define EDMD_VERIFY_COLLISION_ACCOUNTING 0
#endif

/* Focused scheduler tests construct an already prepared particle state in
 * memory. Keep that test state distinct from production restart modes. */
#ifndef EDMD_TEST_PREPARED_STATE
#define EDMD_TEST_PREPARED_STATE 0
#endif

//Values if no load file
int N = 500;
double phi = 0.38;
double sizeratio = 0.5;
double fractionSmallN = 0;
double aspectRatio = 1.0;
int controlLenght = 0;


int Nbig, Nsmall, Nx, Ny;
#if THREE_D
int Nz;
#endif
double px, py, E, pressure, dxr, dyr, dist, active, Ep, MSD;
double COM_x = 0;
double COM_y = 0;
#if THREE_D
double COM_z = 0;
#endif
double Lx, Ly, Lz;
double lastScreen = 0;
double lastCollNum = 0;

double phig = 0.6;
double phil = 0.2;
int Nlx, Nly;

int Nxcells, Nycells;
int collisionCellNeighborRange = 1;
static double collisionGridMaxRadius = 1.0;
#if POLYGON_PARTICLES && !THREE_D
int polygonCollisionCellNeighborRange = 1;
#endif
double halfLx, halfLy; 
double cellxSize, cellySize; 
double halfxC, halfyC; 
double cellxFac, cellyFac; 
#if THREE_D
	int Nzcells;
	double halfLz, cellzSize, halfzC, cellzFac, dzr, pz;
	double collTermZ = 0;
	double collTermXZ = 0;
	double collTermZX = 0;
	double collTermYZ = 0;
	double collTermZY = 0;
#endif

const int suppSizeEvent = 10;

double t = 0;

double collTermX = 0;
double collTermY = 0;
double collTermXY = 0;
double collTermYX = 0;
#if STRESS
double sumKineticX = 0;
double sumKineticY = 0;
double sumKineticXY = 0;
double sumKineticYX = 0;

double runningKineticX = 0;
double runningKineticY = 0;
double runningKineticXY = 0;
double runningKineticYX = 0;
#if THREE_D
double sumKineticZ = 0;
double sumKineticXZ = 0;
double sumKineticZX = 0;
double sumKineticYZ = 0;
double sumKineticZY = 0;
double runningKineticZ = 0;
double runningKineticXZ = 0;
double runningKineticZX = 0;
double runningKineticYZ = 0;
double runningKineticZY = 0;
#endif
#endif

#if HEAT
double sumHeatX = 0;
double sumHeatY = 0;
double runningHeatX = 0;
double runningHeatY = 0;
double collHeatX = 0;
double collHeatY = 0;
#if THREE_D
double sumHeatZ = 0;
double runningHeatZ = 0;
double collHeatZ = 0;
double collisionEnergySource = 0;
#endif
/* Extra third-order peculiar-velocity moments let the Lees-Edwards heat
 * current be integrated exactly between events.  Keeping these aggregates
 * avoids an O(N) particle sweep at every collision. */
static double shearHeatCx2Vy = 0;
static double shearHeatVy3 = 0;
static double shearHeatVyVz2 = 0;
static double shearHeatCxVy2 = 0;
static double shearHeatVy2Vz = 0;
static double shearHeatCxVyVz = 0;
#endif

#if (STRESS || HEAT)
double lastCollTime = 0;
#endif
double dp = 0;
double cumulativeWallImpulse = 0;
double dpLeft = 0;
double dpRight = 0;
double dpMid = 0;

// Temperature/energy/heat-current profile along x.
int tempProfileThermo = 0;
#define TEMP_PROFILE_N_SLABS 40
int tempProfileNSlabs = TEMP_PROFILE_N_SLABS;
double tempSlabWidth = 0;
double JxProfile[TEMP_PROFILE_N_SLABS];
double JyProfile[TEMP_PROFILE_N_SLABS];
#if THREE_D
double JzProfile[TEMP_PROFILE_N_SLABS];
#endif
double energyProfile[TEMP_PROFILE_N_SLABS];
double countTempProfile[TEMP_PROFILE_N_SLABS];

unsigned long int ncol = 0;
unsigned long int ncross = 0;
char fileName[350];
char thermoName[350];
char pressureFieldName[350];
char interfaceName[350];
char osmosisName[350];
char strucName[350];
char transverseCorrelationName[350];
char clusterName[350];
char pcfName[350];
char pcfBondOrderName[350];
char pcfg6Name[350];
char velocityDistributionName[350];
char tempProfileName[350];
FILE* file;

/* The restart source has to outlive command-line parsing so the graphical
 * frontend can stage and reapply it without relying on getopt state. */
char initialFilePath[350] = "restart.dump";



int load = 0;
int S1 = 0;
int Hex = 0;
int MembraneInit = 0;
int hexRadius = 0;

int osmosisThermo = 0;
double sigma = 2;

int interfaceThermo = 0;
int interfaceThermoTemp = 0;

int velocityDistributionThermo = 0;

int profileThermo = 0;
char profileName[350];
#define PROFILE_BIN_CAPACITY 4096
int nProfileBins = 100;

int clusterThermo = 0;
int dumpCluster = 0 ;
int killCluster = 0;
int clusterMode = 0;
int clusterNeighborThreshold = 3;
double clusterCutoff = 3.55;


int strucThermo = 0;
int doFQT = 0;
double qmax = 0.3;

int bubbleStop = 0;

int areaThermo = 0;

int boopThermo = 0;

int pcfThermo = 0;
int pcfBondOrderThermo = 0;
int pcfg6Thermo = 0;

int critical = 0;
int snapshotCritical = 0;
double fracCritical = -1./6.;



double tmax = G ? never : 30000.0;
double firstScreen = 100;
double dtime = 100000;
int periodicDumps = G ? 0 : 1;

double dtimeThermo = 1000;
double firstThermo = 1000;
int Nthermo = 1;


int    logSpacing = 0;
double base = 1.01;


double dr_pressure = 1.0;

//if -1, screenshot will be taken at constant interval of dtimeThermo
double nextScreen = -1;



#if G || defined(EDMD_TEST_MUTABLE_CONFIG)
bool addWell = false;
bool addField = false;
int noise = 0;
int addWally = 0;
int addWallx = 0;
int addWallz = 0;
int addMidWall = 0;
int addCircularWall = 0;
bool damping = false;
bool addDelta = false;
bool addDeltaTangent = false;
bool addDoubleDelta = false;
bool addDoubleDeltaTangent = false;
bool polydispersity = false;
bool thermoWall = 0;
bool charged = 0;
bool addShear = 0;
bool addShearWall = 0;
int addEvolvingDelta = 0;
int addExpo = 0;
int addEnergy = 0;
int addActiveDelta = 0;
int liquidliquid = 0;
#else
/* Headless physics is configured here, not by make.  Keeping these values
 * const lets the compiler erase branches belonging to disabled models. */
const bool addWell = false;
const bool addField = false;
const int noise = 0;
const int addWally = 0;
const int addWallx = 0;
const int addWallz = 0;
const int addMidWall = 0;
const int addCircularWall = 0;
const bool damping = false;
const bool addDelta = false;
const bool addDeltaTangent = false;
const bool addDoubleDelta = false;
const bool addDoubleDeltaTangent = false;
const bool polydispersity = false;
const bool thermoWall = false;
const bool charged = false;
const bool addShear = false;
const bool addShearWall = false;
const int addEvolvingDelta = 0;
const int addExpo = 0;
const int addEnergy = 0;
const int addActiveDelta = 0;
const int liquidliquid = 0;
#endif

void (*collisionEvent)(int);
void (*doTheCollision)(void);
void (*freeFly)(particle*);
void (*doTheWall)(void);
void (*crossingEvent)(int);

//if noise use euler solver
#if G
int euler = 0;
#else
const int euler = 0;
#endif

//save a file with thermodynamics quantities each dtimeThermo
int ther = 1;
int msd = 0;

// Graphical restarts may change output columns without recompiling.
#if G
int reduce = 0;
#else
//reduce the number of properties dumped into the dump files
const int reduce = 0;
#endif
double* posxInitial = NULL;
double* posyInitial = NULL;
#if THREE_D
double* poszInitial = NULL;
#endif

//value for delta model
double delta = 0.5;

//values for double delta model
double deltaM = 2;
double deltam = 0.1;
double ts = 10;

//value for evolving delta model
double tau = 5;
double deltaMax = 0.1;

//values for exponential model
double vo = 3.5;
double ao = 1.5;

//values for square potential model
double sig = 1.6;
double U = -0.8;

//Value of the energy input at collision

#if G
double randomInjection = 0;
int simplifyInjection = 0;
#else
const double randomInjection = 0;
const int simplifyInjection = 0;
#endif
double deltaE = 50;
double beta = 10;
double taur = 3;
double additionalEnergy = 0.01;


//value for the field, must be < 0
double field = -0.02; 

//Initial temperature
double Einit = 1;

//coeff of restitution of the wall
double resW = 1.;
double resT = 1.;

//coeff of restitution of particles
double res = 1;
double resBeta = 0.;

double deltaOmega = 100;

//parameter if noise or damping
double gamm = 0.02;
double T = 0.05;
//time between kicks
double dtnoise = 0.5;

double strainRate = -0.01;
/* Lees-Edwards displacement at time t is
 *   (shearStrainPhase + strainRate*(t - dynamicsStart))*Ly.
 * Keeping the phase separate lets a runtime strain-rate edit preserve the
 * current image displacement instead of making every seam jump. */
static double shearStrainPhase = 0.0;
#if G
double shearWallSpeed = 0.1;
#else
const double shearWallSpeed = 0.1;
#endif

double a = 2;
double b = 5;



double proportionPositivelyCharged = 0.2;
double proportionNeutralyCharged = 0.1;

#if INTRUDER

intruder intru;
double res_intruder = 1;
double delta_intruder = 0.0;

double bathRad = 1;


double Im = 1000000000000000;
double M = 10000000000000000;
double angle = 0;
double omega = 0.0;
double VX = 0.0;
double VY = 0.0;

#if 0
enum intruderShape shape = TRIANGLE;
double paramIntruder[4] = {30, 30};
#elif 0
enum intruderShape shape = SQUARE;
double paramIntruder[4] = {40};
#elif 1
#define sqrt2 1.41421356237
enum intruderShape shape = POLYGON;
double paramIntruder[4] = {10*sqrt2, 4};
#elif 0
enum intruderShape shape = GEN_WHEEL;
double paramIntruder[4] = {20, 40, 5, 0};
//enum intruderShape shape = COMPLEX_WHEEL;
//double paramIntruder[4] = {12.0, 6.5, 2.0, 0.35};
#else
enum intruderShape shape = CROSS;
double paramIntruder[2] = {20, 4};
#endif

const int rethermalize_after_intruder_collision = 0;


double displayRad = 0.;

static IntruderKM km = {0};
static char intruderKMName[350];
int km_enabled = 1;
static double intruderImpulseX = 0.0;
static double intruderImpulseY = 0.0;
static double intruderTorqueZ = 0.0;

double getIntruderKineticModelBurnTime(void){
	return km.burn_time;
}

void setIntruderKineticModelBurnTime(double burnTime){
	km.burn_time = burnTime;
}

#endif


#if G
int addEnergyThermo = 1;
int addPressureThermo = 1;
int thermoPressureField = 0;
int addCollisionNumberThermo = 1;
int addEnergyCenterOfMassThermo = 0;
int thermoKurt = 0;
int thermoSonine = 1;
int unwrap = 0;
#else
const int addEnergyThermo = 1;
const int addPressureThermo = 1;
const int thermoPressureField = 0;
const int addCollisionNumberThermo = 1;
const int addEnergyCenterOfMassThermo = 0;
const int thermoKurt = 0;
const int thermoSonine = 1;
const int unwrap = 0;
#endif


int automaticGrowthRate = 1;
double vr;
int version = 1;
/* Output numbering is independent from the random seed in graphical builds.
 * Keeping it separate also means probing for an unused filename never changes
 * the state reproduced by Restart. */
int outputVersion = 1;

double paulTime = 0;
int actualPaulList = 0;
double dtPaul;
int paulListN;

double t1 = 0;
FILE *fichier;
FILE *thermo;
FILE *interface;
FILE *strucFile;
FILE *clusterFile;
FILE *osmosisFile;
FILE *velocityFile;
FILE *profileFile;
FILE *tempProfileFile;
node *root;
static node *treeMin;
static node *eventStorage;
/* Number of live nodes in the tree and Paul buckets.  Keeping this explicitly
 * avoids an O(N) event-storage scan every time the current tree is empty. */
static size_t scheduledEventCount;
#define EVENT_NODE(index) (eventStorage + (index))
Dump* dump;
//array containing the particles
particle* particles;
//2D array of pointer
particle** cellList;
node** eventList;
node** eventPaul;
//variable containing the nextEvent
node* nextEvent;
typedef struct outputSamplingState outputSamplingState;
struct outputSamplingState{
	int logarithmicIndex;
	int temperatureProfileCount;
	int thermoCount;
	double energyAverage;
	double pressureAverage;
	double pressureXAverage;
	double pressureYAverage;
	double pressureXYAverage;
	double pressureYXAverage;
	double wallPressureAverage;
	int wallPressureCount;
	double midWallPressureAverage;
	double membranePressureDifferenceAverage;
	int midWallPressureCount;
	double velocityProfile[PROFILE_BIN_CAPACITY];
	int velocityProfileCount[PROFILE_BIN_CAPACITY];
	#if THREE_D
	double pressureZAverage;
	double pressureXZAverage;
	double pressureZXAverage;
	double pressureYZAverage;
	double pressureZYAverage;
	#endif
	#if HEAT
	double heatXAverage;
	double heatYAverage;
	#if THREE_D
	double heatZAverage;
	double energySourceRateAverage;
	#endif
	#endif
	double sonineAverage;
	int bubbleCount;
};

static outputSamplingState outputSamples = {0};
static bool periodicDumpStreamOpen = false;

#if G && POLYGON_PARTICLES && !THREE_D
enum {
	GRAPHICAL_RECOVERY_SNAPSHOT_COUNT = 4,
	GRAPHICAL_RECOVERY_AUXILIARY_COUNT = 9
};

typedef struct graphicalRecoverySnapshot graphicalRecoverySnapshot;
struct graphicalRecoverySnapshot{
	bool valid;
	unsigned long generation;
	particle* particleIdentity;
	int particleCount;
	double boxLx, boxLy;
	int wallX, wallY, middleWall, circularWall;
	int thermostatMode;
	double thermostatInterval;
	double time;
	double lastOutputTime;
	double lastCollisionCount;
	unsigned long collisions;
	unsigned long crossings;
	double virialX, virialY, virialXY, virialYX;
	double wallImpulse;
	double sampledWallImpulse;
	double sampledLeftImpulse;
	double sampledRightImpulse;
	double sampledMiddleImpulse;
	#if STRESS
	double kineticX, kineticY, kineticXY, kineticYX;
	double runningKineticX, runningKineticY;
	double runningKineticXY, runningKineticYX;
	#endif
	#if HEAT
	double heatX, heatY, runningHeatX, runningHeatY;
	double collisionHeatX, collisionHeatY;
	double shearCx2Vy, shearVy3, shearVyVz2;
	double shearCxVy2, shearVy2Vz, shearCxVyVz;
	#endif
	#if STRESS || HEAT
	double lastTransportTime;
	#endif
	outputSamplingState samples;
	unsigned long randomState[NMT];
	int randomIndex;
	particle* particleState;
	void* polygonState;
	size_t polygonStateBytes;
	node auxiliary[GRAPHICAL_RECOVERY_AUXILIARY_COUNT];
	bool auxiliaryQueued[GRAPHICAL_RECOVERY_AUXILIARY_COUNT];
};

static graphicalRecoverySnapshot graphicalRecovery[
	GRAPHICAL_RECOVERY_SNAPSHOT_COUNT] = {{0}};
static int graphicalRecoveryNext = 0;
static int graphicalRecoveryCount = 0;
static unsigned long graphicalRecoveryGeneration = 0;
static double graphicalRecoveryLastWallTime = -DBL_MAX;
static double graphicalRecoveryLastSimulationTime = -DBL_MAX;
static unsigned long graphicalRecoveryLastCollisionCount = 0;

static void captureGraphicalRecoveryCheckpoint(const state* screenState);
static bool tryGraphicalPolygonRecovery(state* screenState);
static void releaseGraphicalRecoveryCheckpoints(void);
#endif

#if (STRESS || HEAT)
void initializeTransportAccumulators(void);
#endif

#if EDMD_VERIFY_OVERLAPS
static void verifyNoOverlaps(void);
#endif

#if G
static void simulationSetFault(const char* format, ...){
	if (simulationHasFault())
		return;
	char message[512];
	va_list arguments;
	va_start(arguments, format);
	(void)vsnprintf(message, sizeof(message), format, arguments);
	va_end(arguments);
	message[sizeof(message) - 1] = '\0';
	simulationReportFault(message);
}

/* Execute exactly one interactive frame.  Native builds call this from their
 * ordinary loop; the browser registers the same function with Emscripten so
 * control returns to the browser between frames. */
static bool graphicalFrame(int argc, char** argv, window* screenWindow,
	state* screenState){
	#if !defined(PLATFORM_WEB)
	if (WindowShouldClose()){
		return false;
	}
	#endif
	#if POLYGON_PARTICLES && !THREE_D
	captureGraphicalRecoveryCheckpoint(screenState);
	#endif

	if ((screenState->running || screenState->stepRequested) &&
		!simulationHasFault() && !screenState->wallMoving){
		/* Keep rendering/input responsive even when many potential events are
		 * required per unit of simulated time. */
		const bool singleStep = screenState->stepRequested;
		screenState->stepRequested = false;
		const double simulationFrameStart = GetTime();
		/* A fixed 1/180 s budget is appropriate for a native 144 Hz window,
		 * but leaves two thirds of a browser's ordinary 60 Hz frame idle.  Use
		 * the measured presentation interval when it is slower, while retaining
		 * a hard 1/90 s ceiling so input and rendering cannot be starved. */
		const bool configurationGrowth = freeFly == &freeFlyGrow
			#if POLYGON_PARTICLES && !THREE_D
			|| polygon_particles_growing()
			#endif
			;
		double simulationFrameBudget = 1.0/180.0;
		const double previousFrameDuration = GetFrameTime();
		if (isfinite(previousFrameDuration) &&
			previousFrameDuration > 1.0/120.0){
			simulationFrameBudget = fmin(1.0/90.0,
				fmax(1.0/180.0, 0.70*previousFrameDuration));
		}
		if (configurationGrowth){
			/* Preparation is a transient computation, not the simulation being
			 * observed.  Spend a larger slice on it while retaining enough render
			 * cadence for progress, input, and browser event handling. */
			#if defined(PLATFORM_WEB)
			simulationFrameBudget = fmax(simulationFrameBudget, 1.0/75.0);
			#else
			simulationFrameBudget = fmax(simulationFrameBudget, 1.0/45.0);
			#endif
		}
		do{
			nextEvent = findNextEvent();
			if (nextEvent == NULL){
				simulationSetFault(
					"No finite event remains at t=%.9g. The simulation is paused; "
					"configuration and inspection remain available. Change the setup or restart to continue.",
					t);
				TraceLog(LOG_WARNING, "%s", simulationFaultMessage());
				graphicsPauseSimulation(screenState);
				graphicsReportSimulationFault(simulationFaultMessage());
				break;
			}
			const bool preparingConfiguration = freeFly == &freeFlyGrow
				#if POLYGON_PARTICLES && !THREE_D
				|| polygon_particles_growing()
				#endif
				;
			if (!singleStep && !preparingConfiguration){
				const double pacedTimeLimit = graphicsSimulationTimeLimit(t);
				if (nextEvent->t > pacedTimeLimit){
					/* Particle state is evaluated lazily, so advancing only the
					 * clock to the real-time pacing target is exact. */
					t = pacedTimeLimit;
					#if EDMD_VERIFY_OVERLAPS && POLYGON_PARTICLES && !THREE_D
					if (polygon_particles_active()){
						verifyNoOverlaps();
						if (simulationHasFault()){
							if (tryGraphicalPolygonRecovery(screenState)){
								graphicsReportTransientNotice(
									"A small polygon overlap was rolled back; "
									"collision search is temporarily using finer precision.");
							}
							else{
								graphicsPauseSimulation(screenState);
								graphicsReportSimulationFault(
									simulationFaultMessage());
							}
						}
					}
					#endif
					break;
				}
			}

			t = nextEvent->t;
			removeEventFromQueue(nextEvent);
			switch(nextEvent->type){
				case COLLISION: doTheCollision(); break;
				case WALL: doTheWall(); break;
				case CELLCROSS:
					ncross++;
					doTheCrossing();
					break;
				#if POLYGON_PARTICLES && !THREE_D
				case POLYGON_CELLCROSS:
					ncross++;
					doTheCrossingPolygon();
					break;
				#endif
				case ADDINGNOISE: addNoise(); break;
				case SCREENSHOT: takeAScreenshot(); break;
				case THERMO:
					takeAThermo();
					fflush(stdout);
					break;
				case IN: doIn(); break;
				case OUT: doOut(); break;
				case GROWSTOP:
					stopGrow();
					/* Growth is deliberately unpaced. Production starts a new
					 * physical-time epoch at this completed configuration. */
					graphicsRebaseSimulationPacing(t);
					#if POLYGON_PARTICLES && !THREE_D
					/* Random-growth runs enter polygon dynamics at this event.
					 * Capture the freshly valid production state before another
					 * event in the same render frame can require recovery. */
					captureGraphicalRecoveryCheckpoint(screenState);
					#endif
					break;
				#if INTRUDER
				case INTRUDER_COLLISION: doTheCollisionIntruder(); break;
				#endif
				#if POLYGON_PARTICLES && !THREE_D
					case POLYGON_RECHECK: polygon_particles_do_recheck(); break;
				#endif
			}
			#if EDMD_VERIFY_OVERLAPS && POLYGON_PARTICLES && !THREE_D
				if (!simulationHasFault() && polygon_particles_active())
					verifyNoOverlaps();
			#endif
			if (simulationHasFault()){
				#if POLYGON_PARTICLES && !THREE_D
				if (tryGraphicalPolygonRecovery(screenState)){
					graphicsReportTransientNotice(
						"A small polygon overlap was rolled back; collision "
						"search is temporarily using finer precision.");
					if (singleStep)
						break;
					continue;
				}
				#endif
				graphicsPauseSimulation(screenState);
				graphicsReportSimulationFault(simulationFaultMessage());
				break;
			}
		} while (!singleStep && !simulationHasFault() && t <= tmax &&
			GetTime() - simulationFrameStart < simulationFrameBudget);
	}

	draw(argc, argv, screenWindow, screenState);
	return true;
}

static void graphicalShutdown(state* screenState){
	#if INTRUDER
	km_finalize_and_dump(km, intruderKMName, intru);
	#endif
	printClose();
	outputLifecycleClose();
	#if POLYGON_PARTICLES && !THREE_D
	releaseGraphicalRecoveryCheckpoints();
	#endif
	graphicsFree(screenState);
	CloseWindow();
	freeArrays();
}

#if defined(PLATFORM_WEB)
typedef struct browserGraphicalContext browserGraphicalContext;
struct browserGraphicalContext{
	int argc;
	char** argv;
	window screenWindow;
	state screenState;
	bool shutdown;
	int frameRate;
	double nextRenderTime;
};

static browserGraphicalContext browserContext = {0};

static void browserGraphicalFrame(void* opaque){
	browserGraphicalContext* context = (browserGraphicalContext*)opaque;
	if (context == NULL || context->shutdown){
		emscripten_cancel_main_loop();
		return;
	}
	const int requestedFrameRate = graphicsTargetFrameRate();
	const double now = GetTime();
	if (requestedFrameRate != context->frameRate){
		context->frameRate = requestedFrameRate;
		context->nextRenderTime = now;
	}
	if (requestedFrameRate > 0){
		const double interval = 1.0/requestedFrameRate;
		const double tolerance = fmin(0.0015, 0.1*interval);
		if (now + tolerance < context->nextRenderTime){
			return;
		}
		context->nextRenderTime += interval;
		if (context->nextRenderTime < now){
			context->nextRenderTime = now + interval;
		}
	}
	if (!graphicalFrame(context->argc, context->argv,
		&context->screenWindow, &context->screenState)){
		context->shutdown = true;
		graphicalShutdown(&context->screenState);
		emscripten_cancel_main_loop();
	}
}
#endif
#endif

static const double headlessKineticEnergyCutoff = 1e-10;

/* A low kinetic-energy state is terminal only when no enabled mechanism can
 * replenish it.  In particular, a driven Delta system must not be stopped
 * merely because it passes through a cold transient. */
static bool unforcedDissipationActive(void){
	const bool dissipativeBoundary =
		(addWallx || addWally || addWallz || addMidWall || addCircularWall) &&
		(resW < 1 || resT < 1);
	const bool dissipationActive = damping || res < 1 || dissipativeBoundary;
	const bool kineticEnergyCanRecover = noise || addDelta || addDeltaTangent ||
		addDoubleDelta || addDoubleDeltaTangent || addEvolvingDelta ||
		addExpo || addEnergy || addActiveDelta || thermoWall || addShear ||
		addShearWall || addField || charged || addWell;
	return dissipationActive && !kineticEnergyCanRecover;
}

#if !G
static bool headlessKineticEnergyDepleted(void){
	return unforcedDissipationActive() && isfinite(E) &&
		E < headlessKineticEnergyCutoff;
}

static void reportHeadlessKineticEnergyStop(void){
	printf("\nTotal kinetic energy %.9g fell below %.1e at t=%.9g; "
		"stopping the freely cooling headless simulation.\n",
		E, headlessKineticEnergyCutoff, t);
}
#endif


int main(int argc, char *argv[]){
	#if G
	#if defined(PLATFORM_WEB)
	browserContext.argc = argc;
	browserContext.argv = argv;
	browserContext.screenWindow = graphicalInit();
	browserContext.screenState = GUIinit();
	window* screenWindow = &browserContext.screenWindow;
	state* screenState = &browserContext.screenState;
	#else
	window nativeScreenWindow = graphicalInit();
	state nativeScreenState = GUIinit();
	window* screenWindow = &nativeScreenWindow;
	state* screenState = &nativeScreenState;
	#endif
	#endif


	constantInit(argc, argv);
	/* The version/seed shown by the GUI and command line is the actual seed of
	 * the initial state as well as every explicit restart.  Previously the
	 * first graphical run was seeded from wall-clock time, while the GUI
	 * claimed seed 1, so Restart did not reproduce what was shown. */
	init_genrand((unsigned long)version);
	srand((unsigned int)version);
	configureTransverseVelocityCorrelationSampling(&dtimeThermo, &logSpacing);
	runningCheck();
	particlesInit();
	adaptCellGridToParticles();
	cellListInit();
	if (addWell && U != 0.0){
		pairsInit();
	}
	physicalQ();
	#if (STRESS || HEAT)
	initializeTransportAccumulators();
	#endif
	format();
	char outputError[1024] = {0};
	if (!outputLifecycleOpen(outputError, sizeof(outputError))){
		fprintf(stderr, "WARNING: %s\n", outputError);
		#if G
		graphicsReportRuntimeError(outputError);
		#endif
	}
	eventListInit();
	#if G
	if (simulationHasFault()){
		graphicsPauseSimulation(screenState);
		graphicsReportSimulationFault(simulationFaultMessage());
	}
	#endif

	#if G == 1
	screenWindow->factor = GetScreenHeight()/Ly;
	#if defined(PLATFORM_WEB)
	/* A zero frame rate delegates scheduling to requestAnimationFrame.  The
	 * context is static because the callback outlives main()'s stack frame. */
	emscripten_set_main_loop_arg(browserGraphicalFrame, &browserContext, 0, true);
	#else
	while (graphicalFrame(argc, argv, screenWindow, screenState)){}
	graphicalShutdown(screenState);
	#endif
	return 0;
	#else
	while (t <= tmax){
		nextEvent = findNextEvent();
		if (nextEvent == NULL){
				fprintf(stderr,
					"No finite event remains at t=%.17g; stopping.\n", t);
			break;
		}

		t = nextEvent->t;
		removeEventFromQueue(nextEvent);
		switch(nextEvent->type){
			case COLLISION: doTheCollision(); break;
			case WALL: doTheWall(); break;
			case CELLCROSS:
				ncross++;
				doTheCrossing();
				break;
			#if POLYGON_PARTICLES && !THREE_D
			case POLYGON_CELLCROSS:
				ncross++;
				doTheCrossingPolygon();
				break;
			#endif
			case ADDINGNOISE: addNoise(); break;
			case SCREENSHOT:
				takeAScreenshot();
				if (headlessKineticEnergyDepleted()){
					reportHeadlessKineticEnergyStop();
					goto exiting;
				}
				break;
			case THERMO:
				takeAThermo();
				fflush(stdout);
				if (headlessKineticEnergyDepleted()){
					saveTXT();
					reportHeadlessKineticEnergyStop();
					goto exiting;
				}
				break;
			case IN: doIn(); break;
			case OUT: doOut(); break;
			case GROWSTOP: stopGrow(); break;
			#if INTRUDER
			case INTRUDER_COLLISION: doTheCollisionIntruder(); break;
			#endif
			#if POLYGON_PARTICLES && !THREE_D
			case POLYGON_RECHECK: polygon_particles_do_recheck(); break;
			#endif
		}
		#if EDMD_VERIFY_OVERLAPS && POLYGON_PARTICLES && !THREE_D
			if (polygon_particles_active())
				verifyNoOverlaps();
		#endif
	}

exiting:
	#if INTRUDER
	km_finalize_and_dump(km, intruderKMName, intru);
	#endif
	printInfo();
	printClose();
	outputLifecycleClose();
	freeArrays();
	return 0;
	#endif
}



/* ------------------------------------------/
/											 /
/	INITIALISATION OF ARRAYS AND CONSTANTS 	 /
/											 /
/-------------------------------------------*/

#if !THREE_D
/* `particle.rad` remains the polygon circumradius so the linked-cell stencil
 * stays conservative. Box sizing and native polygon preparation use the area
 * and boundary of the body that is actually simulated and rendered. */
static double unitParticleArea2D(void){
	#if POLYGON_PARTICLES
	if (polygon_particles_requested()){
		const double area = polygon_particles_unit_area();
		if (isfinite(area) && area > 0.0)
			return area;
	}
	#endif
	return M_PI;
}
#endif

/* --------------------------------------------
	Initializes the number of particles N
	(for the densely packed config) and
	constants required for the cellList.
-------------------------------------------- */


#if THREE_D
static int nearestPositiveEven(double value){
	if (!isfinite(value) || value > INT_MAX - 2)
		return 0;
	int result = (int)llround(value/2.0)*2;
	return result < 2 ? 2 : result;
}
#endif

bool hexagonalLatticeDimensions(int requestedCount, double requestedAspect,
	int* nx, int* ny, int* nz){
	if (nx == NULL || ny == NULL || requestedCount <= 0 ||
		!isfinite(requestedAspect) || requestedAspect <= 0)
		return false;

	#if THREE_D
	if (nz == NULL || requestedCount < 4)
		return false;
	/* An HCP layer has Ly/Lx=(sqrt(3)/2) Ny/Nx and the AB stack has
	 * Lz/Lx=sqrt(2/3) Nz/Nx.  Search a small neighbourhood of the
	 * continuous solution because both periodic offsets require even counts. */
	const double layerY = sqrt(3.0)/2.0;
	const double layerZ = sqrt(2.0/3.0);
	const double estimate = cbrt(requestedCount/
		(sqrt(2.0)*requestedAspect*requestedAspect));
	if (!isfinite(estimate) || estimate > INT_MAX - 65)
		return false;
	const int centerX = (int)fmax(1.0, llround(estimate));
	const int firstX = centerX > 64 ? centerX - 64 : 1;
	const int lastX = centerX + 64;
	double bestScore = DBL_MAX;
	long long bestCountError = LLONG_MAX;
	int bestX = 0, bestY = 0, bestZ = 0;
	for (int candidateX = firstX; candidateX <= lastX; candidateX++){
		const int centerY = nearestPositiveEven(
			requestedAspect*candidateX/layerY);
		if (centerY == 0)
			continue;
		for (int yOffset = -4; yOffset <= 4; yOffset += 2){
			const int candidateY = centerY + yOffset;
			if (candidateY < 2)
				continue;
			const int shapeZ = nearestPositiveEven(
				requestedAspect*candidateX/layerZ);
			const int countZ = nearestPositiveEven((double)requestedCount/
				((double)candidateX*candidateY));
			const int centersZ[2] = {shapeZ, countZ};
			for (int zCenterIndex = 0; zCenterIndex < 2; zCenterIndex++){
				if (centersZ[zCenterIndex] == 0)
					continue;
				for (int zOffset = -4; zOffset <= 4; zOffset += 2){
					const int candidateZ = centersZ[zCenterIndex] + zOffset;
					if (candidateZ < 2)
						continue;
					/* The realized lattice must fit both the particle/event index and
					 * this intermediate multiplication.  Checking against INT_MAX a
					 * factor at a time also avoids signed long-long overflow for
					 * pathological aspect ratios. */
					if (candidateX > INT_MAX/candidateY)
						continue;
					const int layerCount = candidateX*candidateY;
					if (candidateZ > INT_MAX/layerCount)
						continue;
					const long long count = (long long)layerCount*candidateZ;
					if (count <= 0)
						continue;
					const long long countError = llabs(count - requestedCount);
					const double aspectY = layerY*candidateY/candidateX;
					const double aspectZ = layerZ*candidateZ/candidateX;
					const double score = 4.0*(double)countError/requestedCount +
						fabs(aspectY/requestedAspect - 1.0) +
						fabs(aspectZ/requestedAspect - 1.0);
					if (score < bestScore ||
						(score == bestScore && countError < bestCountError)){
						bestScore = score;
						bestCountError = countError;
						bestX = candidateX;
						bestY = candidateY;
						bestZ = candidateZ;
					}
				}
			}
		}
	}
	if (bestX == 0)
		return false;
	*nx = bestX;
	*ny = bestY;
	*nz = bestZ;
	#else
	/* Preserve the established triangular-lattice sizing in 2D, while
	 * rejecting the former zero-row division for very small configurations. */
	double yValue = sqrt((double)requestedCount)*requestedAspect;
	if (!isfinite(yValue) || yValue > INT_MAX - 1)
		return false;
	int resultY = (int)yValue;
	if (resultY & 1)
		resultY++;
	if (resultY < 2)
		resultY = 2;
	int resultX = requestedCount/resultY;
	if (resultX < 1)
		return false;
	*nx = resultX;
	*ny = resultY;
	if (nz != NULL)
		*nz = 1;
	#endif
	return true;
}

static void hexagonalLatticeSitePosition(int ix, int iy, int iz,
	double* x, double* y, double* z){
	const double dx = Lx/Nx;
	const double dy = Ly/Ny;
	#if THREE_D
	const int layer = iz & 1;
	*x = fmod((ix + 0.5*(iy & 1) + 0.5*layer)*dx, Lx);
	*y = fmod(iy*dy + layer*dy/3.0, Ly);
	*z = iz*(Lz/Nz);
	#else
	(void)iz;
	*x = ix*dx + (iy & 1)*(dx/2.0);
	*y = iy*dy;
	if (z != NULL)
		*z = 0;
	#endif
}

static bool hexagonalLatticeSiteInsideHole(double x, double y, double z){
	double radiusSquared = (x - halfLx)*(x - halfLx) +
		(y - halfLy)*(y - halfLy);
	#if THREE_D
	radiusSquared += (z - halfLz)*(z - halfLz);
	#else
	(void)z;
	#endif
	return radiusSquared < (double)hexRadius*hexRadius;
}

#if POLYGON_PARTICLES && !THREE_D
enum { POLYGON_CELL_HEADS_PER_PARTICLE_LIMIT = 128 };
#ifndef EDMD_POLYGON_ABSOLUTE_CELL_HEAD_LIMIT
#define EDMD_POLYGON_ABSOLUTE_CELL_HEAD_LIMIT 4194304UL
#endif

static bool polygonCellHeadCountForRange(int range, size_t* count){
	if (count == NULL || range < 1 ||
		!(collisionGridMaxRadius > 0.0) ||
		!isfinite(collisionGridMaxRadius) ||
		!(Lx > 0.0) || !(Ly > 0.0) ||
		!isfinite(Lx) || !isfinite(Ly)){
		return false;
	}
	const double targetCellSize = 2*collisionGridMaxRadius/range;
	double nxValue = floor(Lx/targetCellSize);
	double nyValue = floor(Ly/targetCellSize);
	if (nxValue < 1.0) nxValue = 1.0;
	if (nyValue < 1.0) nyValue = 1.0;
	if (nxValue > INT_MAX || nyValue > INT_MAX)
		return false;
	const size_t nx = (size_t)nxValue;
	const size_t ny = (size_t)nyValue;
	if (ny > SIZE_MAX/nx)
		return false;
	*count = nx*ny;
	return true;
}

/* Empirical dense-CCD tuning: range 2 is best near phi=.1, range 6 near
   phi=.375, and range 12 near phi=.8.  Bound cell-head memory independently
   of particle count, but never below the range-1 grid required for exactness
   in a very dilute/large box. */
static int polygonCellNeighborRangeForBox(void){
	if (!polygon_particles_requested())
		return 1;
	int range = 2;
	if (isfinite(phi) && phi > 0.0){
		const double rounded = floor(16.0*phi + 0.5);
		if (rounded >= 2.0 && rounded <=
			POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE){
			range = (int)rounded;
		}
		else if (rounded > POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE){
			range = POLYGON_PARTICLES_MAX_CELL_NEIGHBOR_RANGE;
		}
	}

	size_t rangeOneHeads = 0;
	if (!polygonCellHeadCountForRange(1, &rangeOneHeads))
		return 1;
	size_t proportionalLimit = SIZE_MAX;
	if (N > 0 && (size_t)N <=
		SIZE_MAX/POLYGON_CELL_HEADS_PER_PARTICLE_LIMIT){
		proportionalLimit =
			(size_t)N*POLYGON_CELL_HEADS_PER_PARTICLE_LIMIT;
	}
	size_t headLimit = proportionalLimit;
	if (headLimit > (size_t)EDMD_POLYGON_ABSOLUTE_CELL_HEAD_LIMIT)
		headLimit = (size_t)EDMD_POLYGON_ABSOLUTE_CELL_HEAD_LIMIT;
	if (headLimit < rangeOneHeads)
		headLimit = rangeOneHeads;

	while (range > 1){
		size_t heads = 0;
		if (polygonCellHeadCountForRange(range, &heads) &&
			heads <= headLimit){
			return range;
		}
		range--;
	}
	return 1;
}
#endif

// Convenient function for graphical part of the program
void boxConstantHelper(){
	halfLx = Lx/2;
	halfLy = Ly/2;
	#if THREE_D
	halfLz = Lz/2;
	#endif
	#if POLYGON_PARTICLES && !THREE_D
	polygonCollisionCellNeighborRange = 1;
	#endif
	if (addWell && U != 0.0){
		/* A finer potential grid reduces the large 3^d candidate-cell volume.
		   The full +/-2 stencil remains conservative for a 2*sig cutoff. */
		#if THREE_D
		collisionCellNeighborRange = 2;
		#else
		collisionCellNeighborRange = (sig >= 2) ? 2 : 1;
		#endif
		double targetCellSize = 2*sig*collisionGridMaxRadius/
			collisionCellNeighborRange;
		Nycells = (int)(Ly/targetCellSize);
		Nxcells = (int)(Lx/targetCellSize);
		#if THREE_D
		Nzcells = (int)(Lz/targetCellSize);
		#endif
	}
	else{
		/* Keep the disk stencil and its hot stack frames unchanged. Polygon CCD
		   uses a separately bounded, density-adaptive fine grid. */
		collisionCellNeighborRange = 1;
		int gridNeighborRange = 1;
		#if POLYGON_PARTICLES && !THREE_D
		polygonCollisionCellNeighborRange =
			polygonCellNeighborRangeForBox();
		gridNeighborRange = polygonCollisionCellNeighborRange;
		#endif
		double targetCellSize = 2*collisionGridMaxRadius/
			gridNeighborRange;
		Nycells = (int)(Ly/targetCellSize);
		Nxcells = (int)(Lx/targetCellSize);
		#if THREE_D
		Nzcells = (int)(Lz/targetCellSize);
		#endif
	}
	if (Nxcells < 1) Nxcells = 1;
	if (Nycells < 1) Nycells = 1;
	#if THREE_D
	if (Nzcells < 1) Nzcells = 1;
	#endif

	cellxSize = Lx/Nxcells;
	cellySize = Ly/Nycells;
	halfxC = cellxSize/2;
	halfyC = cellySize/2;
	cellxFac = 1/cellxSize;
	cellyFac = 1/cellySize;
	#if THREE_D
	cellzSize = Lz/Nzcells;
	halfzC = cellzSize/2;
	cellzFac = 1/cellzSize;
	#endif
	dtPaul = PAUL_DT_SCALE/(double)N;

    paulListN = N;
}

/* File-based initializers may contain particles larger than the built-in
 * unit-radius species. Size the grid from the state that was actually loaded;
 * otherwise a collision cutoff may span more cells than the stencil. Keeping
 * the lower bound at one also covers zero-radius particles that later grow to
 * the built-in target radius. */
void adaptCellGridToParticles(){
	double maxRadius = 1.0;
	for (int i = 0; i < N; i++){
		if (isfinite(particles[i].rad) && particles[i].rad > maxRadius){
			maxRadius = particles[i].rad;
		}
		if (isfinite(vr) && vr != 0 && isfinite(particles[i].vr)){
			double targetRadius = particles[i].rad + fabs(particles[i].vr/vr);
			if (targetRadius > maxRadius)
				maxRadius = targetRadius;
		}
	}
	collisionGridMaxRadius = maxRadius;
	boxConstantHelper();
}


void constantInit(int argc, char *argv[]){

	struct option longopt[] = {
		{"load", required_argument, NULL, 'l'},
		{"number", required_argument, NULL, 'N'},
		{"phi", required_argument, NULL, 'p'},
		{"res", required_argument, NULL, 'r'},
		{"delta", required_argument, NULL, 'd'},
		{"gamma", required_argument, NULL, 'g'},
		{"time", required_argument, NULL, 't'},
		{"dt", required_argument, NULL, 'D'},
		{"initial-energy", required_argument, NULL, 'E'},
		{"potential", required_argument, NULL, 'U'},
		{"well-radius", required_argument, NULL, 'R'},
		{"field", required_argument, NULL, 'f'},
		{"ts", required_argument, NULL, 's'},
		{"xs", required_argument, NULL, 'x'},
		{"sizeratio", required_argument, NULL, 'q'},
		{"temperature", required_argument, NULL, 'T'},
		{"aspect", required_argument, NULL, 'a'},
		{"Ly", required_argument, NULL, 'L'},
		{"Lx", required_argument, NULL, 'X'},
		{"deltaE", required_argument, NULL, 'e'},
		{"additionalEnergy", required_argument, NULL, 'A'},
		{"dtimeThermo", required_argument, NULL, 'o'},
		{"dr_pressure", required_argument, NULL, 1001},
		{"phil", required_argument, NULL, 'P'},
		{"phig", required_argument, NULL, 'G'},
		{"version", required_argument, NULL, 'v'},
		{"strainRate", required_argument, NULL, 'S'},
		{"cluster-mode", required_argument, NULL, 'C'},
		{"cluster-neighbor-threshold", required_argument, NULL, 'H'},
		{"temp-profile", no_argument, NULL, 1007},
		{"temp-slabs", required_argument, NULL, 1008},
		#if !G
		{"msd", no_argument, NULL, 1010},
		{"pair-correlation", no_argument, NULL, 1011},
		{"bond-order-pair-correlation", no_argument, NULL, 1012},
		{"g6-pair-correlation", no_argument, NULL, 1013},
		{"structure-factor", required_argument, NULL, 1014},
		{"velocity-distribution", no_argument, NULL, 1015},
		#endif
		{"intruderMass", required_argument, NULL, 'M'},
		{"intruderDelta", required_argument, NULL, 'z'},
		{"intruderRes", required_argument, NULL, 'Y'},
		{"intruderKMBurn", required_argument, NULL, 'B'},
		{NULL, 0, NULL, 0}
	};
	
	
	int c;

	/* Command-line arguments are a one-time input source. GUI restarts call
	 * this function with argc == 0 and only rebuild derived geometry. */
	if (argc > 0 && argv != NULL){
		optind = 1;
	}
	while (argc > 0 && argv != NULL &&
		(c = getopt_long(argc, argv, "l:N:p:r:d:g:t:D:E:U:R:f:s:x:q:T:a:L:X:e:A:o:P:G:v:S:C:H:M:B:z:Y:n:", longopt, NULL)) != -1){
		switch(c){
			case 'l':
				strncpy(initialFilePath, optarg, sizeof(initialFilePath) - 1);
				initialFilePath[sizeof(initialFilePath) - 1] = '\0';
				/* The only supported restart format is the full dump format,
				   which preserves velocities and complete 3D state. */
				load = 2;
				break;
			case 'N':
				load = 0;
				sscanf(optarg, "%d", &N);
				break;
			case 'p':
				load = 0;
				sscanf(optarg, "%lf", &phi);
				break;
			case 'x':
				load = 0;
				sscanf(optarg, "%lf", &fractionSmallN);
				break;
			case 'a':
				load = 0;
				sscanf(optarg, "%lf", &aspectRatio);
				break;
			case 'q':
				load = 0;
				sscanf(optarg, "%lf", &sizeratio);
				break;
			case 'r':
				sscanf(optarg, "%lf", &res);
				break;
			case 'd':
				if (!(addDelta || addDeltaTangent || addDoubleDelta ||
					addDoubleDeltaTangent || addEvolvingDelta || addActiveDelta))
					printf("WARNING:\033[0;31m Changing Delta while Delta model is not enabled!\033[0m\n");
				sscanf(optarg, "%lf", &delta);
				break;
			case 'g':
				if (!damping && !noise)
					printf("WARNING:\033[0;31m Changing gamma while noise or damping is not enabled!\033[0m\n");
				sscanf(optarg, "%lf", &gamm);
				break;
			case 't':
				sscanf(optarg, "%lf", &tmax);
				break;
			case 'D':
				sscanf(optarg, "%lf", &dtime);
				break;
			case 'E':
				sscanf(optarg, "%lf", &Einit);
				break;
			case 'U':
				if (addWell != 1)
					printf("WARNING:\033[0;31m Changing potential while it is not enabled!\033[0m\n");
				sscanf(optarg, "%lf", &U);
				break;
			case 'R':
				if (addWell != 1)
					printf("WARNING:\033[0;31m Changing potential radius while it is not enabled!\033[0m\n");
				sscanf(optarg, "%lf", &sig);
				break;
			case 'f':
				if (addField != 1)
					printf("WARNING:\033[0;31m Changing field while it is not enabled!\033[0m\n");
				sscanf(optarg, "%lf", &field);
				break;
			case 's':
				sscanf(optarg, "%lf", &ts);
				break;
			case 'T':
				sscanf(optarg, "%lf", &T);
				break;
			case 'L':
				controlLenght = 1;
				load = 0;
				sscanf(optarg, "%lf", &Ly);
				break;
			case 'X':
				controlLenght = 2;
				load = 0;
				sscanf(optarg, "%lf", &Lx);
				break;
			case 'e':
				sscanf(optarg, "%lf", &deltaE);
				break;
			case 'A':
				sscanf(optarg, "%lf", &additionalEnergy);
				break;
			case 'o':
				sscanf(optarg, "%lf", &dtimeThermo);
				break;
			case 1001:
				sscanf(optarg, "%lf", &dr_pressure);
				if (dr_pressure <= 0){
					printf("WARNING:\033[0;31m dr_pressure must be > 0. Forcing dr_pressure = 1.0\033[0m\n");
					dr_pressure = 1.0;
				}
				break;
			case 'P':
				sscanf(optarg, "%lf", &phil);
				break;
			case 'G':
				sscanf(optarg, "%lf", &phig);
				break;
			case 'v':
				sscanf(optarg, "%d", &version);
				outputVersion = version;
				init_genrand(version);
				break;
			case 'S':
				sscanf(optarg, "%lf", &strainRate);
				break;
			case 'C':
				sscanf(optarg, "%d", &clusterMode);
				if (clusterMode != 0 && clusterMode != 1){
					printf("WARNING:\033[0;31m clusterMode should be 0 (dense) or 1 (dilute). Forcing dense.\033[0m\n");
					clusterMode = 0;
				}
				break;
			case 'H':
				sscanf(optarg, "%d", &clusterNeighborThreshold);
				if (clusterNeighborThreshold < 0){
					printf("WARNING:\033[0;31m clusterNeighborThreshold must be >= 0. Forcing 0.\033[0m\n");
					clusterNeighborThreshold = 0;
				}
				break;
			case 1007:
				tempProfileThermo = 1;
				break;
			case 1008:
				sscanf(optarg, "%d", &tempProfileNSlabs);
				if (tempProfileNSlabs < 1){
					tempProfileNSlabs = 1;
				}
				if (tempProfileNSlabs > TEMP_PROFILE_N_SLABS){
					printf("WARNING:\033[0;31m tempProfileNSlabs is capped at %d.\033[0m\n", TEMP_PROFILE_N_SLABS);
					tempProfileNSlabs = TEMP_PROFILE_N_SLABS;
				}
				break;
			#if !G
			case 1010:
				msd = 1;
				break;
			case 1011:
				pcfThermo = 1;
				break;
			case 1012:
				pcfBondOrderThermo = 1;
				break;
			case 1013:
				pcfg6Thermo = 1;
				break;
			case 1014:
				if (strcmp(optarg, "density") == 0 || strcmp(optarg, "1") == 0)
					strucThermo = 1;
				else if (strcmp(optarg, "velocity") == 0 || strcmp(optarg, "2") == 0)
					strucThermo = 2;
				else{
					fprintf(stderr,
						"Unknown structure-factor mode '%s'; use density or velocity.\n",
						optarg);
					exit(2);
				}
				break;
			case 1015:
				velocityDistributionThermo = 1;
				break;
			#endif
			#if INTRUDER
			case 'M':
				sscanf(optarg, "%lf", &M);
				break;
			case 'B':
				km.burn_time = atof(optarg);
				break;
			case 'z':
				sscanf(optarg, "%lf", &delta_intruder);
				break;
			case 'Y':
				sscanf(optarg, "%lf", &res_intruder);
				break;
			case 'n':
				sscanf(optarg, "%lf", &paramIntruder[1]);
			#endif
		}
	}

	#if !THREE_D
	const double unitParticleArea = unitParticleArea2D();
	#endif

	if (load == 2){
		file = fopen(initialFilePath, "r");
		if (file == NULL){
			fprintf(stderr, "Unable to open initial dump '%s'\n", initialFilePath);
			exit(EXIT_FAILURE);
		}
		dump = dump_open(initialFilePath,'r');
		N = get_natoms(dump);
		char hboxx;
		Lx = get_boxx(&hboxx, 1, dump);
		char hboxy;
		Ly = get_boxy(&hboxy, 1, dump);
		#if THREE_D
		char hboxz;
		Lz = get_boxz(&hboxz, 1, dump);
		#endif
	}
	else if (S1){
		Nsmall = N/2;
		Nlx = (int)sqrt(Nsmall/aspectRatio);
		Nly = (int)sqrt(Nsmall*aspectRatio);
		Nsmall = Nlx*Nly;
		Nbig = Nsmall;
		N = Nbig + Nsmall;
		#if THREE_D
		Lx = sqrt(Nsmall*M_PI*(sizeratio*sizeratio + 1)/
			(phi*aspectRatio));
		#else
		Lx = sqrt(Nsmall*unitParticleArea*(sizeratio*sizeratio + 1)/
			(phi*aspectRatio));
		#endif
		Ly = Lx*aspectRatio;
	}
	else if (Hex){
		int latticeZ = 1;
		if (!hexagonalLatticeDimensions(N, aspectRatio, &Nx, &Ny,
			&latticeZ)){
			fprintf(stderr,
				"Unable to derive a complete %dD hexagonal lattice for N=%d and aspect=%.17g\n",
				THREE_D ? 3 : 2, N, aspectRatio);
			exit(EXIT_FAILURE);
		}
		#if THREE_D
		Nz = latticeZ;
		N = Nx*Ny*Nz;
		/* HCP nearest-neighbour spacing.  At phi=pi/(3 sqrt(2)) this
		 * evaluates to the hard-sphere diameter, as it should. */
		const double spacing = cbrt((4.0/3.0)*M_PI/
			(phi*(sqrt(2.0)/2.0)));
		Lx = Nx*spacing;
		Ly = Ny*(sqrt(3.0)/2.0)*spacing;
		Lz = Nz*sqrt(2.0/3.0)*spacing;
		#else
		N = Nx*Ny;
		Ly = sqrt((sqrt(3.0)/2.0)*((double)Ny/Nx)*
			unitParticleArea*N/phi);
		Lx = (2.0/sqrt(3.0))*((double)Nx/Ny)*Ly;
		#endif
		Nbig = N;
		Nsmall = 0;
		if (hexRadius != 0){
			for (int i = 0; i < Nx; i++){
				for (int j = 0; j < Ny; j++){
					#if THREE_D
					for (int k = 0; k < Nz; k++){
					#else
					for (int k = 0; k < 1; k++){
					#endif
						double x, y, z;
						hexagonalLatticeSitePosition(i, j, k, &x, &y, &z);
						if (hexagonalLatticeSiteInsideHole(x, y, z))
							N--;
					}
				}
			}
			Nbig = N;
		}
		if (N <= 0){
			fprintf(stderr,
				"The central hexagonal-lattice hole removes every lattice site\n");
			exit(EXIT_FAILURE);
		}
	}
	else if (addMidWall){
		double val = (phig/phil)*2/(sizeratio*sizeratio);
		fractionSmallN = val/(1 + val);
		#if THREE_D
		Lx = sqrt(M_PI*N*(1-fractionSmallN)/(0.5*aspectRatio*phil));
		#else
		Lx = sqrt(unitParticleArea*N*(1-fractionSmallN)/
			(0.5*aspectRatio*phil));
		#endif
		Ly = Lx*aspectRatio;
		Nsmall = N*fractionSmallN;
		Nbig = N - Nsmall;
		#if THREE_D
		phi = (Nbig + Nsmall*(sizeratio*sizeratio))/(Lx*Ly)*M_PI;
		#else
		phi = (Nbig + Nsmall*(sizeratio*sizeratio))/(Lx*Ly)*
			unitParticleArea;
		#endif
	}
	else{
		#if THREE_D
		double r2 = (1 + fractionSmallN*(sizeratio*sizeratio*sizeratio - 1));
		#else
		double r2 = (1 + fractionSmallN*(sizeratio*sizeratio - 1));
		#endif
		if (polydispersity){
			#if THREE_D
			r2 = b*tgamma(1 + 3/a)*tgamma(b)/tgamma(1 + b + 3/a);
			#else
			r2 = b*tgamma(1 + 2/a)*tgamma(b)/tgamma(1 + b + 2/a);
			#endif
			Ly = Lx;
		}
		if (addCircularWall){
			#if THREE_D
			/* The packing fraction in a spherical container is
			   sum(radius^3)/(L/2)^3. */
			Lx = cbrt(8*N*r2/phi);
			Ly = Lx;
			Lz = Lx;
			#else
			Lx = sqrt((4.0/M_PI)*unitParticleArea*N/phi*r2);
			Ly = Lx;
			#endif
		}
		else{
			if (controlLenght == 0){
				#if THREE_D
				Lx = pow((4./3.)*M_PI*N/(phi*aspectRatio*aspectRatio)*r2, 1./3.);
				Lz = Lx*aspectRatio;
				#else
				Lx = sqrt(unitParticleArea*N/(phi*aspectRatio)*r2);
				#endif
				Ly = Lx*aspectRatio;
				
			}
			// NO 3d
			else if (controlLenght == 1){
				Lx = Ly/aspectRatio;
				Lz = Ly;
			}
			else{
				Ly = Lx*aspectRatio;
				#if THREE_D
				N = (int)(phi*Lx*Ly/(M_PI*r2));
				phi = M_PI*N*r2/(Lx*Ly);
				#else
				N = (int)(phi*Lx*Ly/(unitParticleArea*r2));
				phi = unitParticleArea*N*r2/(Lx*Ly);
				#endif
				Lz = Ly;
			}

		}
		
	}

	#if INTRUDER
	intruderInit(&intru, Lx/2, Ly/2, VX, VY, angle, omega, M, Im, shape, paramIntruder, km_enabled, &km);

	

	if (Hex){
		for (int i = 0; i < Nx; i++){
			for (int j = 0; j < Ny; j++){
				double X = i*Lx/Nx + (j%2)*(Lx/(2*Nx));
				double Y = j*Ly/Ny;
				particle p;
				p.x = X;
				p.y = Y;
				p.rad = 1;
				if (particleOutsideIntruder(&p, &intru)){
					N--;
				}
			}
		}
	}
	#endif

	if ((critical > 0 || snapshotCritical > 0) && (fracCritical < 0)){
		fracCritical = 2*aspectRatio;
	}

	boxConstantHelper();
}

static void appendOutputError(char* error, size_t errorCapacity,
	const char* format, ...){
	if (error == NULL || errorCapacity == 0){
		return;
	}
	size_t used = strlen(error);
	if (used >= errorCapacity - 1){
		return;
	}
	if (used != 0){
		int written = snprintf(error + used, errorCapacity - used, "; ");
		if (written < 0 || (size_t)written >= errorCapacity - used){
			return;
		}
		used += (size_t)written;
	}
	va_list arguments;
	va_start(arguments, format);
	(void)vsnprintf(error + used, errorCapacity - used, format, arguments);
	va_end(arguments);
	error[errorCapacity - 1] = '\0';
}

static FILE* openOutputStream(const char* path, const char* description,
	char* error, size_t errorCapacity){
	errno = 0;
	FILE* stream = fopen(path, "w");
	if (stream == NULL){
		appendOutputError(error, errorCapacity, "Unable to open %s '%s': %s",
			description, path, strerror(errno));
	}
	return stream;
}

bool outputThermoSamplingRequested(void){
	return ther || thermoPressureField || tempProfileThermo || profileThermo ||
		interfaceThermo || osmosisThermo || strucThermo ||
		transverseVelocityCorrelationIsEnabled() ||
		(clusterThermo && (ther || dumpCluster || killCluster > 0)) ||
		pcfThermo || pcfBondOrderThermo || pcfg6Thermo ||
		velocityDistributionThermo || bubbleStop ||
		(ther && (critical || boopThermo || msd)) ||
		/* Preserve the headless low-energy termination even when ordinary
		 * thermodynamic output and periodic snapshots are disabled. */
		(!G && unforcedDissipationActive());
}

static bool outputSamplingNeedsSynchronizedParticles(void){
	return damping || msd || tempProfileThermo || profileThermo ||
		interfaceThermo || osmosisThermo || strucThermo ||
		transverseVelocityCorrelationIsEnabled() || clusterThermo ||
		pcfThermo || pcfBondOrderThermo || pcfg6Thermo ||
		velocityDistributionThermo || critical || bubbleStop;
}

static void reportOutputRuntimeFailure(const char* format, ...){
	char message[768];
	va_list arguments;
	va_start(arguments, format);
	(void)vsnprintf(message, sizeof(message), format, arguments);
	va_end(arguments);
	message[sizeof(message) - 1] = '\0';
	fprintf(stderr, "WARNING: %s\n", message);
	#if G
	graphicsReportRuntimeError(message);
	#endif
}

void resetOutputSamplingState(void){
	memset(&outputSamples, 0, sizeof(outputSamples));
	memset(JxProfile, 0, sizeof(JxProfile));
	memset(JyProfile, 0, sizeof(JyProfile));
	#if THREE_D
	memset(JzProfile, 0, sizeof(JzProfile));
	#endif
	memset(energyProfile, 0, sizeof(energyProfile));
	memset(countTempProfile, 0, sizeof(countTempProfile));
	free(posxInitial);
	free(posyInitial);
	#if THREE_D
	free(poszInitial);
	poszInitial = NULL;
	COM_z = 0.0;
	#endif
	posxInitial = NULL;
	posyInitial = NULL;
	COM_x = 0.0;
	COM_y = 0.0;
	MSD = 0.0;
	/* These accumulators describe the interval beginning at lastScreen.  A
	 * phase change (notably growth -> production) must not leak its collision
	 * impulses into the first production sample.  The monotonic live-analysis
	 * counter deliberately remains untouched. */
	collTermX = 0.0;
	collTermY = 0.0;
	collTermXY = 0.0;
	collTermYX = 0.0;
	#if THREE_D
	collTermZ = 0.0;
	collTermXZ = 0.0;
	collTermZX = 0.0;
	collTermYZ = 0.0;
	collTermZY = 0.0;
	#endif
	dp = 0.0;
	dpLeft = 0.0;
	dpRight = 0.0;
	dpMid = 0.0;
	lastScreen = t;
}

void outputLifecycleClose(void){
	periodicDumpStreamOpen = false;
	if (fichier != NULL){
		fclose(fichier);
		fichier = NULL;
	}
	if (thermo != NULL){
		fclose(thermo);
		thermo = NULL;
	}
	if (interface != NULL){
		fclose(interface);
		interface = NULL;
	}
	if (strucFile != NULL){
		fclose(strucFile);
		strucFile = NULL;
	}
	if (clusterFile != NULL){
		fclose(clusterFile);
		clusterFile = NULL;
	}
	if (osmosisFile != NULL){
		fclose(osmosisFile);
		osmosisFile = NULL;
	}
	if (velocityFile != NULL){
		fclose(velocityFile);
		velocityFile = NULL;
	}
	if (profileFile != NULL){
		fclose(profileFile);
		profileFile = NULL;
	}
	if (tempProfileFile != NULL){
		fclose(tempProfileFile);
		tempProfileFile = NULL;
	}
	pressure_field_close_file();
	pressure_field_free();
	freeInterface();
	free_osmosis();
	freeStructureFactor();
	finalizeTransverseVelocityCorrelation();
	resetOutputSamplingState();
}

static void writeThermoHeader(void){
	if (thermo == NULL){
		return;
	}
		fprintf(thermo, "t ");
		if (addCollisionNumberThermo){
			fprintf(thermo, "Ncol ");
		}
		if (addEnergyThermo){
			fprintf(thermo, "E ");
		}
		#if THREE_D
		if (addPressureThermo){
			fprintf(thermo, "p px py pz pxy pyx pxz pzx pyz pzy ");
		}
		#if HEAT
		fprintf(thermo, "Jx Jy Jz Edot ");
		#endif
		#else
		if (addPressureThermo){
			fprintf(thermo, "p px py pxy pyx ");
		}
		#if HEAT
		fprintf(thermo, "Jx Jy ");
		#endif
		#endif
		if (addCircularWall || addWallx || addWally || addWallz){
			fprintf(thermo, "pwall ");
		}
		if (addMidWall){
			fprintf(thermo, "pmid dp packLeft packRight tempSolventLeft tempSolventRight tempSolute ");
		}
		if (addWell){
			fprintf(thermo, "U ");
		}
		if (msd){
			fprintf(thermo, "MSD ");
		}
		if (addEnergyCenterOfMassThermo){
			fprintf(thermo, "ECOM ");
		}
		if (critical){
			#if THREE_D
			fprintf(thermo, "dense1 dense2 dense3 dense4 dilute1 dilute2 dilute3 dilute4 Edense1 Edense2 Edense3 Edense4 Edilute1 Edilute2 Edilute3 Edilute4 ");
			#else
			fprintf(thermo, "dense1 dense2 dilute1 dilute2 Edense1 Edense2 Edilute1 Edilute2 ");
			#endif
		}
		if (clusterThermo){
			fprintf(thermo, "largestCluster ");
		}
		if (boopThermo){
			fprintf(thermo, "q6 ");
		}
		#if INTRUDER
			fprintf(thermo, "intruderX intruderY intruderVX intruderVY intruderAngle intruderOmega intruderFx intruderFy intruderTorque ");
		#endif
		if (thermoKurt){
			fprintf(thermo, "K1 K3 ratioKurt ");
		}
		if (thermoSonine){
			fprintf(thermo, "a2 ");
		}
		fprintf(thermo, "\n");
	}

bool outputLifecycleOpen(char* error, size_t errorCapacity){
	bool success = true;
	if (error != NULL && errorCapacity > 0){
		error[0] = '\0';
	}
	outputLifecycleClose();
	resetOutputSamplingState();
	if (!periodicDumps && !outputThermoSamplingRequested()){
		/* The graphical live-analysis path is deliberately memory-only.  Avoid
		 * even creating/probing an output directory until the user explicitly
		 * asks for a manual snapshot. */
		return true;
	}
	customName();
	if (fileName[0] == '\0'){
		appendOutputError(error, errorCapacity,
			"Unable to construct output filenames");
		return false;
	}

	if (periodicDumps){
		fichier = openOutputStream(fileName, "periodic dump", error,
			errorCapacity);
		periodicDumpStreamOpen = fichier != NULL;
		success = periodicDumpStreamOpen && success;
		if (!periodicDumpStreamOpen){
			periodicDumps = 0;
		}
	}
	if (ther){
		thermo = openOutputStream(thermoName, "thermodynamic stream", error,
			errorCapacity);
		if (thermo != NULL){
			writeThermoHeader();
		}
		else{
			success = false;
			#if G
			ther = 0;
			#endif
		}
	}
	if (thermoPressureField){
		if (!pressure_field_init(Lx, Ly, dr_pressure, pressureFieldName)){
			appendOutputError(error, errorCapacity,
				"Unable to initialize pressure-field output '%s'",
				pressureFieldName);
			success = false;
			#if G
			thermoPressureField = 0;
			#endif
		}
	}
	if (interfaceThermo){
		interface = openOutputStream(interfaceName, "interface stream", error,
			errorCapacity);
		if (interface != NULL && interfaceThermoTemp == 0){
			if (liquidliquid == 0){
				#if THREE_D
				initLocalDensity(Lx, Ly, Lz, 8, 8, 8, N, interface, liquidliquid, 0);
				#else
				initLocalDensity(Lx, Ly, 8, 8, N, interface, liquidliquid, 0);
				#endif
			}
			else{
				#if THREE_D
				initLocalDensity(Lx, Ly, Lz, 4, 4, 4, N, interface, liquidliquid, 0);
				#else
				initLocalDensity(Lx, Ly, 4, 4, N, interface, liquidliquid, 0);
				#endif
			}
		}
		else if (interface != NULL){
			#if THREE_D
			initLocalDensity(Lx, Ly, Lz, 4, 4, 4, N, interface, liquidliquid, 1);
			#else
			initLocalDensity(Lx, Ly, 4, 4, N, interface, liquidliquid, 1);
			#endif
		}
		if (interface == NULL){
			success = false;
			interfaceThermo = 0;
			interfaceThermoTemp = 0;
		}
	}
	if (osmosisThermo){
		osmosisFile = openOutputStream(osmosisName, "osmosis stream", error,
			errorCapacity);
		if (osmosisFile == NULL || !init_osmosis(Lx, Ly, sigma, osmosisFile)){
			appendOutputError(error, errorCapacity,
				"Unable to initialize osmosis output '%s'", osmosisName);
			if (osmosisFile != NULL){
				fclose(osmosisFile);
				osmosisFile = NULL;
			}
			free_osmosis();
			success = false;
			osmosisThermo = 0;
		}
	}
	if (strucThermo){
		strucFile = openOutputStream(strucName, "structure-factor stream",
			error, errorCapacity);
		if (strucFile == NULL ||
			!initStructureFactor(qmax, Lx, Ly, N, strucFile, doFQT)){
			appendOutputError(error, errorCapacity,
				"Unable to initialize structure-factor output '%s'", strucName);
			if (strucFile != NULL){
				fclose(strucFile);
				strucFile = NULL;
			}
			freeStructureFactor();
			success = false;
			strucThermo = 0;
		}
	}
	if (transverseVelocityCorrelationIsEnabled()){
		if (!initTransverseVelocityCorrelation(Lx, Ly, N,
			transverseCorrelationName)){
			appendOutputError(error, errorCapacity,
				"Unable to initialize TVACF output '%s'",
				transverseCorrelationName);
			success = false;
			int enabled, modeM, modeN, samples, checkpoint;
			double span;
			getTransverseVelocityCorrelationConfig(&enabled, &modeM, &modeN,
				&samples, &span, &checkpoint);
			setTransverseVelocityCorrelationConfig(0, modeM, modeN, samples,
				span, checkpoint);
		}
	}
	if (dumpCluster && clusterThermo){
		clusterFile = openOutputStream(clusterName, "cluster histogram", error,
			errorCapacity);
		success = clusterFile != NULL && success;
		if (clusterFile == NULL){
			dumpCluster = 0;
		}
	}
	if (velocityDistributionThermo){
		velocityFile = openOutputStream(velocityDistributionName,
			"velocity distribution", error, errorCapacity);
		if (velocityFile != NULL){
			#if THREE_D
			fprintf(velocityFile, "vx vy vz\n");
			#else
			fprintf(velocityFile, "vx vy\n");
			#endif
		}
		else{
			success = false;
			velocityDistributionThermo = 0;
		}
	}
	if (profileThermo){
		profileFile = openOutputStream(profileName, "velocity profile", error,
			errorCapacity);
		if (profileFile != NULL){
			for (int i = 0; i < nProfileBins; i++){
				fprintf(profileFile, "%lf ", (i + 0.5)*Ly/nProfileBins);
			}
			fprintf(profileFile, "\n");
		}
		else{
			success = false;
			profileThermo = 0;
		}
	}
	if (tempProfileThermo){
		tempProfileFile = openOutputStream(tempProfileName,
			"temperature/heat profile", error, errorCapacity);
		tempSlabWidth = Lx/((double)tempProfileNSlabs);
		if (tempProfileFile){
			#if THREE_D
			fprintf(tempProfileFile,
				"# x-slab-centres; blocks: kinetic-energy Jx Jy Jz\n");
			#else
			fprintf(tempProfileFile,
				"# x-slab-centres; blocks: kinetic-energy Jx Jy\n");
			#endif
			for (int i = 0; i < tempProfileNSlabs; i++){
				fprintf(tempProfileFile, "%lf ", (i + 0.5)*tempSlabWidth);
			}
			fprintf(tempProfileFile, "\n");
		}
		else{
			success = false;
			tempProfileThermo = 0;
		}
	}
	return success;
}

void initThermo(){
	char error[1024] = {0};
	if (!outputLifecycleOpen(error, sizeof(error))){
		fprintf(stderr, "WARNING: %s\n", error);
	}
}


/* ------------------------------------------
	Initializes the array containing all
	the particles and their properties
------------------------------------------ */
static double randomCoordinateWithMargin(double lower, double upper,
		double preferredMargin){
	/* Several wall initializers historically sampled [2.5, L-2.5]
	 * directly.  Very small boxes (or the half-box used by the membrane)
	 * can make that interval empty, which feeds reversed bounds to drand().
	 * Preserve the usual clearance whenever it fits and continuously reduce
	 * it for a narrow interval; the overlap/fault machinery can then report
	 * a physically impossible dense state instead of corrupting positions. */
	if (!isfinite(lower) || !isfinite(upper) || upper <= lower)
		return 0.5*(lower + upper);
	const double width = upper - lower;
	double margin = isfinite(preferredMargin) && preferredMargin > 0.0
		? preferredMargin : 0.0;
	const double largestSafeMargin = width*(0.5 - 8.0*DBL_EPSILON);
	if (margin > largestSafeMargin)
		margin = largestSafeMargin;
	const double minimum = lower + margin;
	const double maximum = upper - margin;
	return maximum > minimum ? drand(minimum, maximum) : 0.5*(lower + upper);
}

void particlesInit(){

	particles = calloc(N, sizeof(particle));
	if (particles == NULL){
		fprintf(stderr, "Particle allocation failed for N=%d\n", N);
		exit(EXIT_FAILURE);
	}

	
	

	int i = 0;
	if (load == 2){
		int nframes = dump->nframes;
		jump_to_frame(nframes - 1, dump);

		double* x = calloc(N, sizeof(double));
		double* y = calloc(N, sizeof(double));
		double* vx = calloc(N, sizeof(double));
		double* vy = calloc(N, sizeof(double));
		#if THREE_D
		double* z = calloc(N, sizeof(double));
		double* vz = calloc(N, sizeof(double));
		#endif
		double* rad = calloc(N, sizeof(double));
		double* m = calloc(N, sizeof(double));
		int* type = calloc(N, sizeof(int));
		if ((x == NULL) || (y == NULL) || (vx == NULL) || (vy == NULL) ||
			(rad == NULL) || (m == NULL) || (type == NULL)
			#if THREE_D
			|| (z == NULL) || (vz == NULL)
			#endif
		){
			fprintf(stderr, "Restart-buffer allocation failed for N=%d\n", N);
			free(x);
			free(y);
			free(vx);
			free(vy);
			#if THREE_D
			free(z);
			free(vz);
			#endif
			free(rad);
			free(m);
			free(type);
			exit(EXIT_FAILURE);
		}

		get_doubleatomprop("x", x, N, dump);
		get_doubleatomprop("y", y, N, dump);
		#if THREE_D
		get_doubleatomprop("z", z, N, dump);
		#endif
		/* Input requirements are independent from the format selected for future
		 * output. GUI users may load a full restart and request reduced dumps for
		 * the new run; consulting `reduce` here previously left mass/radius/type
		 * arrays at zero and corrupted that otherwise valid restart. */
		get_doubleatomprop("vx", vx, N, dump);
		get_doubleatomprop("vy", vy, N, dump);
		#if THREE_D
		get_doubleatomprop("vz", vz, N, dump);
		#endif
		get_doubleatomprop("radius", rad, N, dump);
		get_doubleatomprop("m", m, N, dump);
		get_intatomprop("type", type, N, dump);
		
		
		particle* p;
		for (i = 0; i < N; i++){
		    p = particles + i;
		    p->x = x[i];
			p->y = y[i];
			PBCpostX(&(p->x));
			PBCpostY(&(p->y));
			#if THREE_D
			p->z = z[i];
			PBCpostZ(&(p->z));
			#endif
			p->type = type[i];
			p->vx = vx[i];
			p->vy = vy[i];
			#if THREE_D
			p->vz = vz[i];
			#endif
			p->rad = rad[i];
			p->m = m[i];
		    Nbig += p->type;
			p->num = i;
			p->t = 0;
			p->coll = 0;
            p->lastColl = 0;
		}
		fclose(file);
		Nsmall = N - Nbig;
		sizeratio = particles[N - 1].rad;
		free(x);
		free(y);
		free(vx);
		free(vy);
		#if THREE_D
		free(z);
		free(vz);
		#endif
		free(type);
		free(m);
		free(rad);
		dump_close(dump);
	}
	else if (S1){
		particle* p;
		int k = 0;
	    float dx = Lx/Nlx;
		float dy = Ly/Nly;
		for (i = 0; i < Nlx; i++){
			for (int j = 0; j < Nly; j++){
				p = particles + k;
				p->x = i*dx;
				p->y = j*dy;
				p->type = 1;
				p->rad = 1;
				p->m = 1;
				double U1 = drand(0, 1);
				double U2 = drand(0, 1);
				p->vx = sqrt(-2*log(U1)/p->m*Einit)*cos(2*M_PI*U2);
				p->vy = sqrt(-2*log(U1)/p->m*Einit)*sin(2*M_PI*U2);
				p->num = k;
				p->t = 0;
				p->coll = 0;
				p->lastColl = 0;
				k++;
			}
		}
		
		for (i = 0; i < Nlx; i++){
			for (int j = 0; j < Nly; j++){
				p = particles + k;
				p->x = dx/2 + i*dx;
				p->y = dy/2 + j*dy;
				p->type = 0;
				p->rad = sizeratio;
				p->m = sizeratio*sizeratio*sizeratio;
				double U1 = drand(0, 1);
				double U2 = drand(0, 1);
				p->vx = sqrt(-2*log(U1)/p->m*Einit)*cos(2*M_PI*U2);
				p->vy = sqrt(-2*log(U1)/p->m*Einit)*sin(2*M_PI*U2);
				p->num = k;
				p->t = 0;
				p->coll = 0;
				p->lastColl = 0;
				k++;
			}
		}
		normalizePhysicalQ();
	}
	else if (Hex){
		particle* p;
		int particleIndex = 0;
		for (i = 0; i < Nx; i++){
			for (int j = 0; j < Ny; j++){
				#if THREE_D
				for (int layer = 0; layer < Nz; layer++){
				#else
				for (int layer = 0; layer < 1; layer++){
				#endif
				p = particles + particleIndex;
				double siteZ = 0;
				hexagonalLatticeSitePosition(i, j, layer,
					&p->x, &p->y, &siteZ);
				#if THREE_D
				p->z = siteZ;
				#endif
				if (hexRadius != 0){
					if (hexagonalLatticeSiteInsideHole(p->x, p->y, siteZ)){
						continue;
					}
				}
				#if INTRUDER
				if (particleOutsideIntruder(p, &intru)){
					continue;
				}
				#endif
				if (liquidliquid){
					if (particleIndex < N/2)
						p->type = 1;
					else
						p->type = 0;
				}
				else{
					p->type = 1;
				}
				p->rad = 1;
				p->m = 1;
				double U1 = drand(0, 1);
				double U2 = drand(0, 1);
				p->vx = sqrt(-2*log(U1)/p->m*Einit)*cos(2*M_PI*U2);
				p->vy = sqrt(-2*log(U1)/p->m*Einit)*sin(2*M_PI*U2);
				#if THREE_D
				double U3 = drand(0, 1);
				double U4 = drand(0, 1);
				p->vz = sqrt(-2*log(U3)/p->m*Einit)*cos(2*M_PI*U4);
				#endif
				if (hexRadius != 0){
					double X = p->x;
					double Y = p->y;
					double dist = sqrt((X - Lx/2)*(X - Lx/2) + (Y - Ly/2)*(Y - Ly/2));
					/* A 3D spherical hole can retain a site on its z axis.  Its
					 * projected XY radius is then zero, so the legacy planar swirl
					 * has no defined tangent and must not divide by zero. */
					if (dist > DBL_MIN){
						double velocity_magnitude = 2*hexRadius/dist;
						double angle = atan2(Y - Ly/2, X - Lx/2) + M_PI/2;
						p->vx += velocity_magnitude*cos(angle);
						p->vy += velocity_magnitude*sin(angle);
					}
				}
				p->num = particleIndex;
				p->t = 0;
				p->coll = 0;
				p->lastColl = 0;
				particleIndex++;
				}
			}
		}
		if (particleIndex != N){
			fprintf(stderr,
				"Hexagonal lattice count mismatch: built %d particles, expected %d\n",
				particleIndex, N);
			exit(EXIT_FAILURE);
		}
		normalizePhysicalQ();
	}
	else{
		
		double EinitGrow = 0.05;
		particle* p;
		if (automaticGrowthRate){
			optimizeGrowConstant();
		}
		double *targetRad = calloc(N, sizeof(double));

		if (polydispersity){
			for (int i = 0; i < N; i++)	
				targetRad[i] = pow(1 - pow(1 - drand(0, 1), 1/b), 1/a);
			qsort(targetRad, N, sizeof(double), cmpDouble);
		}
		else{
			Nsmall= (int)(N*fractionSmallN);
			Nbig = N - Nsmall;
		}
		for (int i = 0; i < N; i++){
			p = particles + i;
			if (addCircularWall){
				do{
					p->x = drand(0, Lx);
					p->y = drand(0, Ly);
					#if THREE_D
					p->z = drand(0, Lz);
					#endif
				}
				#if THREE_D
				while ((p->x - halfLx)*(p->x - halfLx) +
					(p->y - halfLy)*(p->y - halfLy) +
					(p->z - halfLz)*(p->z - halfLz) > halfLx*halfLx);
				#else
				while ((p->x - halfLx)*(p->x - halfLx) +
					(p->y - halfLy)*(p->y - halfLy) > halfLx*halfLx);
				#endif
			}
			else{
				if (addMidWall && i < Nbig){
					p->x = randomCoordinateWithMargin(0.0, Lx/2.0, 2.5);
				}
				else{
					p->x = (addWallx || addMidWall)
						? randomCoordinateWithMargin(0.0, Lx, 2.5)
						: drand(0.0, Lx);
				}
				p->y = addWally
					? randomCoordinateWithMargin(0.0, Ly, 2.5)
					: drand(0.0, Ly);
				#if THREE_D
				p->z = addWallz
					? randomCoordinateWithMargin(0.0, Lz, 2.5)
					: drand(0.0, Lz);
				#endif
			}
			p->num = i;
			p->rad = 0;
			p->crossX = 0;
			p->crossY = 0;
			#if THREE_D
			p->crossZ = 0;
			#endif

			if (polydispersity){
				p->vr = vr*targetRad[i];
				p->type = 1;
				p->m = pow(targetRad[i], 3);
			}
			else{
				if (i < Nbig){
					p->type = 1;
					p->m = 1;
					p->vr = vr;
					
				}
				else{
					p->type = 0;
					p->m = pow(sizeratio, 3);
					p->vr = vr*sizeratio;
				}
			}

			#if TANGENTIAL
			p->J = 1;
			#endif
			p->t = 0;
			p->coll = 0;

			double U1 = drand(0, 1);
			double U2 = drand(0, 1);
			p->vx = sqrt(-2*log(U1)/p->m*EinitGrow)*cos(2*M_PI*U2);
			p->vy = sqrt(-2*log(U1)/p->m*EinitGrow)*sin(2*M_PI*U2);
			#if THREE_D
			/* A separate Gaussian keeps the initial 3D state isotropic. */
			double U3 = drand(0, 1);
			double U4 = drand(0, 1);
			p->vz = sqrt(-2*log(U3)/p->m*EinitGrow)*cos(2*M_PI*U4);
			#endif

			p->lastColl = 0;
			
			// Needed only if charged == 1/true..
			double randomTemp = drand(0, 1);
			if (randomTemp < proportionPositivelyCharged){
				p->charge = 1;
			}
			else if (randomTemp < proportionNeutralyCharged + proportionPositivelyCharged){
				p->charge = 0;
			}
			else{
				p->charge = -1;
			}
			#if INTRUDER
				p->vr = vr*bathRad;
			#endif
			#if POLYGON_PARTICLES && !THREE_D
			if (polygon_particles_requested()){
				/* SAT requires a strictly positive scale. Start the actual body
				 * almost point-like and reduce its rate so t=1/vr still lands
				 * exactly on the configured species radius. */
				const double targetRadius = p->vr/vr;
				const double initialFraction = 1.0e-6;
				p->rad = targetRadius*initialFraction;
				p->vr = vr*(targetRadius - p->rad);
			}
			#endif
			
		}
		if (polydispersity){
			free(targetRad);
		}
	}


	if (load){
		if (addCircularWall){
			#if THREE_D
			phi = (sizeratio*sizeratio*sizeratio*Nsmall + Nbig)/
				(halfLx*halfLx*halfLx);
			#else
			phi = unitParticleArea2D()*
				(sizeratio*sizeratio*Nsmall + Nbig)/
				(M_PI*halfLx*halfLx);
			#endif
		}
		else{
			#if THREE_D
			phi = (4./3.)*M_PI*(sizeratio*sizeratio*sizeratio*Nsmall + Nbig)/
				(Lx*Ly*Lz);
			#else
			phi = unitParticleArea2D()*
				(sizeratio*sizeratio*Nsmall + Nbig)/(Lx*Ly);
			#endif
		}
	}

}


/* ---------------------------------------------
	Initializes the list of particles close!
--------------------------------------------- */

static double potentialPairStepEnergy(const particle* first,
	const particle* second);

/* Returns the distinct periodic cell coordinates in the active collision stencil.
   Keeping this list unique matters when a short box has only one or two cells:
   the old loops visited the same linked list as many as 27 times in 3D. */
static int neighboringCellCoordinates(int center, int count, int coordinates[5]){
	/* The ordinary hard-disk grid always uses the +/-1 stencil.  Keep its
	 * overwhelmingly common case branch-only: building the same three entries
	 * through wrapping and duplicate-search loops is paid for every collision
	 * reschedule.  The one- and two-cell cases retain the distinct-coordinate
	 * guarantee needed by very small periodic boxes. */
	if (collisionCellNeighborRange == 1){
		if (count <= 1){
			coordinates[0] = 0;
			return 1;
		}
		coordinates[0] = center == 0 ? count - 1 : center - 1;
		coordinates[1] = center;
		if (count == 2)
			return 2;
		coordinates[2] = center + 1 == count ? 0 : center + 1;
		return 3;
	}

	int number = 0;
	for (int offset = -collisionCellNeighborRange;
		offset <= collisionCellNeighborRange; offset++){
		int coordinate = center + offset;
		while (coordinate < 0)
			coordinate += count;
		while (coordinate >= count)
			coordinate -= count;
		int duplicate = 0;
		for (int i = 0; i < number; i++){
			if (coordinates[i] == coordinate){
				duplicate = 1;
				break;
			}
		}
		if (!duplicate)
			coordinates[number++] = coordinate;
	}
	return number;
}

void pairsInit(){
	for(int i = 0; i < N; i++){
		particle* p1 = particles + i;
		p1->particlesInWell = NULL;
		p1->numberOfParticlesInWell = 0;
		p1->particlesInWellCapacity = 0;

		int X = p1->cell[0];
		int Y = p1->cell[1];
		int neighboringX[5];
		int neighboringY[5];
		int numberX = neighboringCellCoordinates(X, Nxcells, neighboringX);
		int numberY = neighboringCellCoordinates(Y, Nycells, neighboringY);
		#if THREE_D
		int Z = p1->cell[2];
		int neighboringZ[5];
		int numberZ = neighboringCellCoordinates(Z, Nzcells, neighboringZ);
		#else
		int numberZ = 1;
		#endif

		for (int zi = 0; zi < numberZ; zi++){
			for (int yi = 0; yi < numberY; yi++){
				for (int xi = 0; xi < numberX; xi++){
				#if THREE_D
					particle* p2 = cellList[neighboringZ[zi]*Nycells*Nxcells +
						neighboringY[yi]*Nxcells + neighboringX[xi]];
				#else
					particle* p2 = cellList[neighboringY[yi]*Nxcells + neighboringX[xi]];
				#endif
				while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
						if (p1->num != p2->num &&
							potentialPairStepEnergy(p1, p2) != 0.0){
						double dx = p2->x - p1->x;
						double dy = p2->y - p1->y;
						double cutoff = sig*(p1->rad + p2->rad);
						#if THREE_D
						double dz = p2->z - p1->z;
						PBC(&dx, &dy, &dz);
						if (dx*dx + dy*dy + dz*dz < cutoff*cutoff)
							addParticleInWellList(p1, p2->num);
						#else
						PBC(&dx, &dy);
						if (dx*dx + dy*dy < cutoff*cutoff)
							addParticleInWellList(p1, p2->num);
						#endif
					}
					p2 = p2->nxt;
				}
				}
			}
		}
	}

	fflush(stdout);
}

/* ---------------------------------------------
	Initializes the cellList. cellList[X*Nycells*Y]
	points toward a single particle contained
	in the cell {X, Y}, all the others particles
	in the same cell are reachable by a
	double linked list containing cellList[X*Nycells*Y]
--------------------------------------------- */
void cellListInit(){
	if (Nxcells <= 0 || Nycells <= 0
		#if THREE_D
		|| Nzcells <= 0
		#endif
	){
		fprintf(stderr, "Invalid cell-list dimensions: %d x %d",
			Nxcells, Nycells);
		#if THREE_D
		fprintf(stderr, " x %d", Nzcells);
		#endif
		fprintf(stderr, "\n");
		exit(EXIT_FAILURE);
	}
	size_t cellCount = (size_t)Nxcells;
	if ((size_t)Nycells > SIZE_MAX/cellCount){
		fprintf(stderr, "Cell-list dimensions overflow: %d x %d\n",
			Nxcells, Nycells);
		exit(EXIT_FAILURE);
	}
	cellCount *= (size_t)Nycells;
	#if THREE_D
	if ((size_t)Nzcells > SIZE_MAX/cellCount){
		fprintf(stderr, "Cell-list dimensions overflow: %d x %d x %d\n",
			Nxcells, Nycells, Nzcells);
		exit(EXIT_FAILURE);
	}
	cellCount *= (size_t)Nzcells;
	#endif
	/* Cell indices are stored and combined as int throughout the event code. */
	if (cellCount > (size_t)INT_MAX ||
		cellCount > SIZE_MAX/sizeof(particle*)){
		fprintf(stderr, "Cell list is too large (%zu cells)\n", cellCount);
		exit(EXIT_FAILURE);
	}
	cellList = calloc(cellCount, sizeof(particle*));
	if (cellList == NULL){
		fprintf(stderr, "Cell-list allocation failed for %zu cells\n", cellCount);
		exit(EXIT_FAILURE);
	}


	for (int i = 0; i < N; i++){
		particles[i].nxt = NULL;
		particles[i].prv = NULL;
		addToCell(i);
	}
}


/* ---------------------------------------------
	Initializes the eventList and the binary
	search tree. cellCrossing for the particle
	i is contained in eventList[i], the first
	collision for particle i is contained
	in eventList[N + i], some space is left
	in eventList for other events or utilitaries.

	We then calculate the next (possible)
	crossings and collisions for each particles
	at the end of the function
--------------------------------------------- */
static bool initializationNeedsGrowthPhase(void){
	#if EDMD_TEST_PREPARED_STATE
	return false;
	#else
	return load == 0 && Hex == 0 && S1 == 0;
	#endif
}

static double postGrowthDynamicsStartTime(void){
	return initializationNeedsGrowthPhase() ? 1/vr : 0;
}

static void selectNormalDynamicsHandlers(void){
	#if POLYGON_PARTICLES && !THREE_D
	if (polygon_particles_requested()){
		polygon_particles_activate();
		collisionEvent = &polygon_particles_collision_event;
		doTheCollision = &polygon_particles_do_collision;
		freeFly = &polygon_particles_free_fly;
		doTheWall = &polygon_particles_do_wall;
		crossingEvent = &polygon_particles_crossing_event;
		return;
	}
	#endif
	collisionEvent = &collisionEventNormal;
	doTheCollision = &doTheCollisionNormal;
	freeFly = &freeFlyNormal;
	doTheWall = &doTheWallNormal;
	crossingEvent = &crossingEventNormal;
}

void eventListInit(){
	/* A newly constructed run always begins from the configured zero-strain
	 * reference.  Incremental updates do not call this initializer. */
	shearStrainPhase = 0.0;
	#if G && POLYGON_PARTICLES && !THREE_D
	graphicalRecoveryGeneration++;
	if (graphicalRecoveryGeneration == 0)
		graphicalRecoveryGeneration = 1;
	#endif

	if ((N <= 0) || (suppSizeEvent < 0) ||
		(N > (INT_MAX - suppSizeEvent)/2)){
		fprintf(stderr, "Invalid particle count for event storage: N=%d\n", N);
		exit(EXIT_FAILURE);
	}
	int totalSize = 2*N + suppSizeEvent;
	scheduledEventCount = 0;
	/* Random preparations contain zero-radius particles which must grow.
	 * Lattice initializers are already at their final radii. */
	const bool initializationNeedsGrowth = initializationNeedsGrowthPhase();
	const double growthDuration = postGrowthDynamicsStartTime();
	#if POLYGON_PARTICLES && !THREE_D
	if (!polygon_particles_prepare()){
		fprintf(stderr, "Polygon-particle state allocation failed\n");
		exit(EXIT_FAILURE);
	}
	#endif
	eventList = calloc(totalSize, sizeof(node*));
	eventPaul = calloc(paulListN + 1, sizeof(node*));
	eventStorage = calloc(totalSize, sizeof(node));

	if ((eventList == NULL) || (eventPaul == NULL) || (eventStorage == NULL)){
		fprintf(stderr, "Event-list allocation failed\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < totalSize; i++){
		eventList[i] = eventStorage + i;
		EVENT_NODE(i)->q = -1;
	}

	root = EVENT_NODE(2*N);
	root->t = never + 1; //arbitrary value, root->rgt will always be NULL.
	root->rgt = NULL;
	root->lft = NULL;
	root->top = NULL;
	treeMin = NULL;

	
	double firstThermoEvent = firstThermo;
	double firstDumpEvent = nextScreen >= 0.0 ? nextScreen : firstScreen;
	if (firstDumpEvent < 0.0){
		firstDumpEvent = tmax - 1.0;
	}
	if (initializationNeedsGrowth){
		if (firstDumpEvent < growthDuration){
			firstDumpEvent = growthDuration + firstDumpEvent;
		}
		if (firstThermoEvent < growthDuration){
			firstThermoEvent = growthDuration + firstThermoEvent;
		}
	}
	if (outputThermoSamplingRequested()){
		addEventThermo(firstThermoEvent + 0.000000001);
	}
	if (initializationNeedsGrowth){
		#if G != 1
		tmax += growthDuration;
		#endif
		addEventGrow(growthDuration);
		#if POLYGON_PARTICLES && !THREE_D
		if (polygon_particles_requested()){
			polygon_particles_begin_growth();
			collisionEvent = &polygon_particles_collision_event;
			doTheCollision = &polygon_particles_do_collision;
			freeFly = &polygon_particles_free_fly;
			doTheWall = &polygon_particles_do_wall;
			crossingEvent = &polygon_particles_crossing_event;
		}
		else
		#endif
		{
			collisionEvent = &collisionEventGrow;
			doTheCollision = &doTheCollisionGrow;
			freeFly = &freeFlyGrow;
			doTheWall = &doTheWallGrow;
			crossingEvent = &crossingEventGrow;
		}

	}
	else{
		selectNormalDynamicsHandlers();
	}
	/* The interactive renderer has its own wall-clock frame loop. Periodic
	 * simulation-time snapshots are opt-in because each one synchronizes all
	 * particles and can dominate dense/potential runs. */
	if (periodicDumps && firstDumpEvent < never){
		addEventScreenshot(firstDumpEvent);
	}
	if (noise){
		/* Preparation velocities are deliberately discarded at GROWSTOP.
		 * Deferring the thermostat avoids injecting unrelated noise into the
		 * homothetic polygon packing stage. */
		double firstNoiseTime = dtnoise;
		#if POLYGON_PARTICLES && !THREE_D
		if (polygon_particles_growing())
			firstNoiseTime = growthDuration + dtnoise;
		#endif
		addEventNoise(firstNoiseTime);
	}
	for (int i = 0; i < N; i++){
		EVENT_NODE(i)->i = i;
		EVENT_NODE(N + i)->i = i;
		crossingEvent(i);
		collisionEvent(i);
	}
	#if INTRUDER
	if (!initializationNeedsGrowth)
		collisionEventIntruder();
	#endif
	

}

/* Bring the value-only, lazily advanced particle records to the current
 * global simulation time before changing a trajectory law or topology. */
void synchronizeParticlesAtCurrentTime(void){
	if (particles == NULL || freeFly == NULL)
		return;
	for (int i = 0; i < N; i++)
		freeFly(particles + i);
}

/* Replace only the two predictions owned by each particle.  Auxiliary
 * SCREENSHOT/THERMO/NOISE/GROWSTOP nodes remain in their original calendar
 * slots with their original deadlines. */
void rescheduleParticleEventsAtCurrentTime(void){
	if (particles == NULL || eventList == NULL || eventStorage == NULL ||
		crossingEvent == NULL || collisionEvent == NULL)
		return;
	for (int i = 0; i < N; i++){
		removeEventFromQueue(EVENT_NODE(i));
		removeEventFromQueue(EVENT_NODE(N + i));
	}
	for (int i = 0; i < N; i++){
		crossingEvent(i);
		collisionEvent(i);
	}
	#if INTRUDER
	if (freeFly == &freeFlyNormal){
		removeEventFromQueue(EVENT_NODE(2*N + 8));
		collisionEventIntruder();
	}
	#endif
	nextEvent = NULL;
}

/* A thermostat edit affects no pair or crossing prediction.  Replacing this
 * one auxiliary node avoids disturbing output cadence or the growth stop. */
void rescheduleThermostatEventAtCurrentTime(void){
	if (eventList == NULL || eventStorage == NULL)
		return;
	removeEventFromQueue(EVENT_NODE(2*N + 2));
	if (noise){
		double nextNoiseTime = t + dtnoise;
		#if POLYGON_PARTICLES && !THREE_D
		if (polygon_particles_growing()){
			const double growthStopTime = EVENT_NODE(2*N + 5)->t;
			if (isfinite(growthStopTime) && growthStopTime >= t)
				nextNoiseTime = growthStopTime + dtnoise;
		}
		#endif
		addEventNoise(nextNoiseTime);
	}
	nextEvent = NULL;
}

static void clearPotentialMemberships(void){
	if (particles == NULL)
		return;
	for (int i = 0; i < N; i++){
		free(particles[i].particlesInWell);
		particles[i].particlesInWell = NULL;
		particles[i].numberOfParticlesInWell = 0;
		particles[i].particlesInWellCapacity = 0;
	}
}

/* Species/charge rules alter only which nearby pairs own a square-well
 * membership.  Keep the linked-cell allocation intact for these edits. */
void rebuildPotentialMembershipsAtCurrentTime(void){
	clearPotentialMemberships();
	if (particles != NULL && addWell && U != 0.0)
		pairsInit();
}

/* Rebuild the structures whose geometry depends on the square-well cutoff.
 * Particle allocation and every physical/runtime counter remain owned by the
 * current run.  The caller reschedules particle events after this returns. */
void rebuildSpatialIndexAtCurrentTime(void){
	if (particles == NULL)
		return;
	const int oldNxcells = Nxcells;
	const int oldNycells = Nycells;
	#if THREE_D
	const int oldNzcells = Nzcells;
	#endif
	clearPotentialMemberships();
	adaptCellGridToParticles();
	const bool gridDimensionsChanged = oldNxcells != Nxcells ||
		oldNycells != Nycells
		#if THREE_D
		|| oldNzcells != Nzcells
		#endif
	;
	if (gridDimensionsChanged){
		free(cellList);
		cellList = NULL;
		cellListInit();
	}
	else{
		bool ownershipChanged = false;
		for (int i = 0; i < N && !ownershipChanged; i++){
			ownershipChanged = particles[i].cell[0] !=
				coordToCell(particles[i].x, 1) ||
				particles[i].cell[1] != coordToCell(particles[i].y, 0)
				#if THREE_D
				|| particles[i].cell[2] != coordToCell(particles[i].z, 2)
				#endif
			;
		}
		if (ownershipChanged){
			size_t cellCount = (size_t)Nxcells*(size_t)Nycells
				#if THREE_D
				*(size_t)Nzcells
				#endif
			;
			memset(cellList, 0, cellCount*sizeof(*cellList));
			for (int i = 0; i < N; i++){
				particles[i].nxt = NULL;
				particles[i].prv = NULL;
				addToCell(i);
			}
		}
	}
	if (addWell && U != 0.0)
		pairsInit();
}


void freeArrays(){
	#if POLYGON_PARTICLES && !THREE_D
	polygon_particles_release();
	#endif
	/* Ownership, rather than the currently selected model, decides what is
	 * released. This keeps potential -> hard-particle restarts safe. */
	if (particles != NULL){
		for (int i = 0; i < N; i++){
			free(particles[i].particlesInWell);
			particles[i].particlesInWell = NULL;
		}	
	}
	free(particles);
	particles = NULL;
	free(cellList);
	cellList = NULL;
	free(eventStorage);
	eventStorage = NULL;
	free(eventList);
	eventList = NULL;
	free(eventPaul);
	eventPaul = NULL;
	scheduledEventCount = 0;
	free(posxInitial);
	posxInitial = NULL;
	free(posyInitial);
	posyInitial = NULL;
	#if THREE_D
	free(poszInitial);
	poszInitial = NULL;
	#endif

	if (interfaceThermo){
		freeInterface();
	}
	pressure_field_free();

}
/* --------------------------/
/							 /
/	CellList utilitaries 	 /
/							 /
/---------------------------*/


/* --------------------------------------
   	Adds Particle i to the CellList
-------------------------------------- */
void addToCell(int i){

	particle* p = particles + i;
	int X = coordToCell(p->x, 1);
	int Y = coordToCell(p->y, 0);
	


	p->cell[0] = X;
	p->cell[1] = Y;
	p->prv = NULL;

	#if THREE_D
	int Z = coordToCell(p->z, 2);
	p->cell[2] = Z;
	p->nxt = cellList[Z*Nycells*Nxcells + Y*Nxcells + X];
	cellList[Z*Nycells*Nxcells + Y*Nxcells + X] = p;
	#else
	p->nxt = cellList[Y*Nxcells + X];
	cellList[Y*Nxcells + X] = p;
	#endif
	
	if (p->nxt)
		p->nxt->prv = p;

}

/* ----------------------------------------
   	Removes Particle i from the CellList
-----------------------------------------*/
void removeFromCell(int i){
	particle* p = particles + i;
	if (p->prv == NULL)
		#if THREE_D
		cellList[p->cell[2]*Nycells*Nxcells + p->cell[1]*Nxcells + p->cell[0]] = p->nxt;
		#else
		cellList[p->cell[1]*Nxcells + p->cell[0]] = p->nxt;
		#endif
	else
		p->prv->nxt = p->nxt;
	if (p->nxt != NULL)
		p->nxt->prv = p->prv;
}


int coordToCell(double a, int x){
	if (x == 1){
		int coordinate = (int)(a*cellxFac);
		if (coordinate < 0) return 0;
		if (coordinate >= Nxcells) return Nxcells - 1;
		return coordinate;
	}
	#if THREE_D
	else if (x == 2){
		int coordinate = (int)(a*cellzFac);
		if (coordinate < 0) return 0;
		if (coordinate >= Nzcells) return Nzcells - 1;
		return coordinate;
	}
	#endif
	int coordinate = (int)(a*cellyFac);
	if (coordinate < 0) return 0;
	if (coordinate >= Nycells) return Nycells - 1;
	return coordinate;
}


int PBCcellX(int a){
	if (a < 0)
		return a + Nxcells;
	else if (a >= Nxcells)
		return a - Nxcells;
	return a;
}

int PBCcellY(int a){
	if (a < 0)
		return a + Nycells;
	else if (a >= Nycells)
		return a - Nycells;
	return a;
}
#if THREE_D
int PBCcellZ(int a){
	if (a < 0)
		return a + Nzcells;
	else if (a >= Nzcells)
		return a - Nzcells;
	return a;
}
#endif

/* ---------------------------------------/
/										  /
/	Calendar/Event queue utilitaries 	  /
/										  /
/----------------------------------------*/

/* -ffast-math may fold isfinite() to true.  Inspecting the IEEE-754 exponent
 * keeps NaN/Inf sentinels out of the calendar in optimized builds too. */
static bool eventTimeCanBeQueued(double eventTime){
	uint64_t bits = 0;
	memcpy(&bits, &eventTime, sizeof(bits));
	const bool finite =
		(bits & UINT64_C(0x7ff0000000000000)) !=
		UINT64_C(0x7ff0000000000000);
	return finite && eventTime < never;
}

static void clearEventQueueLinks(node* eventNode){
	eventNode->lft = NULL;
	eventNode->rgt = NULL;
	eventNode->top = NULL;
}

/* Classify one detached, finite event without changing the live-node count.
 * Overflow redistribution uses this path because those nodes never leave the
 * queue while their bucket assignment is refreshed. */
static void queueDetachedEvent(node* toAdd){
	double dt = toAdd->t - paulTime;

	if (dt < dtPaul){
		toAdd->q = actualPaulList;
		addEventToTree(toAdd);
		return;
	}

	int paulListIndex;
	const double horizon = dtPaul*paulListN;
	if (dt >= horizon){
		paulListIndex = paulListN;
	}
	else{
		/* Floating division immediately below the horizon can round to N.
		 * Without this clamp the modulo result is the current bucket even
		 * though the node is stored in a linked list, so removal treats it as
		 * a tree node and corrupts the queue. */
		int bucketOffset = (int)(dt/dtPaul);
		if (bucketOffset < 1)
			bucketOffset = 1;
		if (bucketOffset >= paulListN)
			bucketOffset = paulListN - 1;
		paulListIndex = actualPaulList + bucketOffset;
		if (paulListIndex >= paulListN)
			paulListIndex -= paulListN;
	}

	toAdd->q = paulListIndex;
	toAdd->lft = NULL;
	toAdd->top = NULL;
	toAdd->rgt = eventPaul[paulListIndex];
	if (toAdd->rgt != NULL)
		toAdd->rgt->lft = toAdd;
	eventPaul[paulListIndex] = toAdd;
}

/* --------------------------------
	Adds an event to the queue
-------------------------------- */
void addEventToQueue(node* toAdd){
	if (toAdd == NULL)
		return;

	/* Replacing an already scheduled fixed-storage node must never leave a
	 * second link to that same node in the calendar. */
	if (toAdd->q >= 0)
		removeEventFromQueue(toAdd);
	clearEventQueueLinks(toAdd);
	toAdd->q = -1;

	if (!eventTimeCanBeQueued(toAdd->t))
		return;

	queueDetachedEvent(toAdd);
	scheduledEventCount++;
}

/* --------------------------------
	Adds an event to the BST
-------------------------------- */
void addEventToTree(node* toAdd){
	node* walker = root;
	while (1){
		if (toAdd->t < walker->t){
			if (walker->lft == NULL){
				walker->lft = toAdd;
				break;
			}
			walker = walker->lft;
		}
		else{
			if (walker->rgt == NULL){
				walker->rgt = toAdd;
				break;
			}
			walker = walker->rgt;
		}
	}
	toAdd->lft = NULL;
	toAdd->rgt = NULL;
	toAdd->top = walker;
	if ((treeMin == NULL) || (toAdd->t < treeMin->t))
		treeMin = toAdd;
}


/* -----------------------------------------
	Removes an event from the event queue
----------------------------------------- */
void removeEventFromQueue(node* toRemove){
	if (toRemove == NULL)
		return;
	if (toRemove->q < 0)
		return;

	if (toRemove->q != actualPaulList){
		if (toRemove->rgt != NULL)
			toRemove->rgt->lft = toRemove->lft;
		if (toRemove->lft != NULL)
			toRemove->lft->rgt = toRemove->rgt;
		else
			eventPaul[toRemove->q] = toRemove->rgt;
	}
	else
		removeEventFromTree(toRemove);

	toRemove->q = -1;
	clearEventQueueLinks(toRemove);
	if (scheduledEventCount > 0)
		scheduledEventCount--;
	#if EDMD_VERIFY_QUEUE
	else{
		fprintf(stderr, "Event-queue count underflow at t=%.17g\n", t);
		exit(3);
	}
	#endif
}



/* -----------------------------------
	Removes an event from the BST
------------------------------------ */
void removeEventFromTree(node* toRemove){
	node* top = toRemove->top;
	node* replacement;

	if (toRemove == treeMin){
		if (toRemove->rgt != NULL){
			treeMin = toRemove->rgt;
			while (treeMin->lft != NULL)
				treeMin = treeMin->lft;
		}
		else
			treeMin = (top == root) ? NULL : top;
	}

	if ((toRemove->lft == NULL) && (toRemove->rgt == NULL))
		replacement = NULL;
	else if (toRemove->lft == NULL){
		replacement = toRemove->rgt;
		replacement->top = top;
	}
	else if (toRemove->rgt == NULL){
		replacement = toRemove->lft;
		replacement->top = top;
	}
	else{
		replacement = toRemove->rgt;
		while (replacement->lft != NULL)
			replacement = replacement->lft;

		if (replacement->top != toRemove){
			replacement->top->lft = replacement->rgt;
			if (replacement->rgt != NULL)
				replacement->rgt->top = replacement->top;
			toRemove->rgt->top = replacement;
			replacement->rgt = toRemove->rgt;
		}
		toRemove->lft->top = replacement;
		replacement->lft = toRemove->lft;
		replacement->top = top;
	}

	if (top->lft == toRemove)
		top->lft = replacement;
	else
		top->rgt = replacement;
}


static void moveCurrentPaulBucketToTree(void){
	node* temp = eventPaul[actualPaulList];
	eventPaul[actualPaulList] = NULL;
	while (temp != NULL){
		node* next = temp->rgt;
		clearEventQueueLinks(temp);
		temp->q = actualPaulList;
		addEventToTree(temp);
		temp = next;
	}
}

static void redistributePaulOverflow(void){
	node* temp = eventPaul[paulListN];
	eventPaul[paulListN] = NULL;
	while (temp != NULL){
		node* next = temp->rgt;
		temp->q = -1;
		clearEventQueueLinks(temp);
		if (eventTimeCanBeQueued(temp->t)){
			queueDetachedEvent(temp);
		}
		else if (scheduledEventCount > 0){
			/* Defensive only: new invalid events are rejected on admission. */
			scheduledEventCount--;
		}
		temp = next;
	}
}

static bool regularPaulBucketsAreEmpty(void){
	for (int index = 0; index < paulListN; index++){
		if (eventPaul[index] != NULL)
			return false;
	}
	return true;
}

static double earliestPaulOverflowTime(void){
	double earliest = never;
	for (node* eventNode = eventPaul[paulListN]; eventNode != NULL;
		eventNode = eventNode->rgt){
		if (eventNode->t < earliest)
			earliest = eventNode->t;
	}
	return earliest;
}

/* ---------------------------------------------
	Advance to the next occupied Paul bucket
--------------------------------------------- */

void addNextPaulEvent(){
	/* Before the next wrap no overflow event can enter the active horizon, so
	 * empty buckets can be skipped in one update without changing ordering. */
	const int stepsToWrap = paulListN - actualPaulList;
	for (int steps = 1; steps < stepsToWrap; steps++){
		int candidate = actualPaulList + steps;
		if (eventPaul[candidate] != NULL){
			actualPaulList = candidate;
			paulTime += steps*dtPaul;
			moveCurrentPaulBucketToTree();
			return;
		}
	}

	/* If overflow is the whole buffered queue, its minimum is the global next
	 * event.  Rebase the empty calendar there instead of walking potentially
	 * billions of empty five-unit epochs. */
	if (regularPaulBucketsAreEmpty() && eventPaul[paulListN] != NULL){
		double earliest = earliestPaulOverflowTime();
		if (earliest > paulTime)
			paulTime = earliest;
		redistributePaulOverflow();
		return;
	}

	/* A wrap may make an overflow event earlier than regular buckets on the
	 * other side of the ring.  Stop at that boundary, redistribute first, and
	 * let the next call choose the now-visible minimum bucket. */
	actualPaulList = 0;
	paulTime += stepsToWrap*dtPaul;
	redistributePaulOverflow();
	moveCurrentPaulBucketToTree();
}


/* ---------------------------------------------
	Finds the next event according to the BST
--------------------------------------------- */
node* findNextEvent(){
	while (treeMin == NULL){
		if (scheduledEventCount == 0)
			return NULL;
		addNextPaulEvent();
	}
	#if EDMD_VERIFY_QUEUE
		node* verifiedMin = root->lft;
		while ((verifiedMin != NULL) && (verifiedMin->lft != NULL))
			verifiedMin = verifiedMin->lft;
		if (verifiedMin != treeMin){
			fprintf(stderr, "Event-tree minimum invariant failed at t=%.17g\n", t);
			exit(3);
		}
	#endif
	return treeMin;
}

#if G && POLYGON_PARTICLES && !THREE_D
static bool graphicalRecoveryConfigurationMatches(
	const graphicalRecoverySnapshot* snapshot){
	return snapshot != NULL && snapshot->valid &&
		snapshot->generation == graphicalRecoveryGeneration &&
		snapshot->particleIdentity == particles &&
		snapshot->particleCount == N &&
		snapshot->boxLx == Lx && snapshot->boxLy == Ly &&
		snapshot->wallX == addWallx && snapshot->wallY == addWally &&
		snapshot->middleWall == addMidWall &&
		snapshot->circularWall == addCircularWall &&
		snapshot->thermostatMode == noise &&
		snapshot->thermostatInterval == dtnoise;
}

static bool reserveGraphicalRecoveryStorage(
	graphicalRecoverySnapshot* snapshot, size_t polygonBytes){
	if (snapshot == NULL || N <= 0 || polygonBytes == 0)
		return false;
	if (snapshot->particleCount != N || snapshot->particleState == NULL){
		particle* resized = realloc(snapshot->particleState,
			(size_t)N*sizeof(*resized));
		if (resized == NULL)
			return false;
		snapshot->particleState = resized;
	}
	if (snapshot->polygonStateBytes != polygonBytes ||
		snapshot->polygonState == NULL){
		void* resized = realloc(snapshot->polygonState, polygonBytes);
		if (resized == NULL)
			return false;
		snapshot->polygonState = resized;
		snapshot->polygonStateBytes = polygonBytes;
	}
	return true;
}

static void captureGraphicalRecoveryCheckpoint(const state* screenState){
	/* Growth checkpoints require phase flags and preparation progress baselines
	 * in addition to the ordinary rigid state. Keep recovery production-only;
	 * GROWSTOP clears this ring and immediately seeds a fresh valid snapshot. */
	if (!polygon_particles_active() || polygon_particles_growing() ||
		particles == NULL || eventStorage == NULL ||
		eventList == NULL || eventPaul == NULL || simulationHasFault() ||
		(screenState != NULL && screenState->wallMoving))
		return;
	const double wallTime = GetTime();
	bool configurationChanged = true;
	if (graphicalRecoveryCount > 0){
		const int latest = (graphicalRecoveryNext +
			GRAPHICAL_RECOVERY_SNAPSHOT_COUNT - 1) %
			GRAPHICAL_RECOVERY_SNAPSHOT_COUNT;
		configurationChanged =
			!graphicalRecoveryConfigurationMatches(&graphicalRecovery[latest]);
	}
	if (!configurationChanged &&
		wallTime - graphicalRecoveryLastWallTime < 0.10)
		return;
	if (!configurationChanged &&
		wallTime - graphicalRecoveryLastWallTime < 0.25 &&
		ncol >= graphicalRecoveryLastCollisionCount &&
		ncol - graphicalRecoveryLastCollisionCount < 256)
		return;
	if (!configurationChanged && t == graphicalRecoveryLastSimulationTime &&
		ncol == graphicalRecoveryLastCollisionCount)
		return;

	graphicalRecoverySnapshot* snapshot =
		&graphicalRecovery[graphicalRecoveryNext];
	snapshot->valid = false;
	const size_t polygonBytes = polygon_particles_checkpoint_bytes();
	if (!reserveGraphicalRecoveryStorage(snapshot, polygonBytes))
		return;

	/* Keep the engine's lazy trajectories lazy in the common path.  The
	 * particle record includes its trajectory origin (`lastColl`), while the
	 * polygon sidecar includes the corresponding angular origin.  Together
	 * they are an exact snapshot at snapshot->time without advancing every
	 * particle on each render-frame checkpoint.  The rare restore path
	 * synchronizes the copied values once before rebuilding cell ownership. */
	memcpy(snapshot->particleState, particles,
		(size_t)N*sizeof(*snapshot->particleState));
	if (!polygon_particles_capture_checkpoint(snapshot->polygonState,
		polygonBytes))
		return;

	snapshot->generation = graphicalRecoveryGeneration;
	snapshot->particleIdentity = particles;
	snapshot->particleCount = N;
	snapshot->boxLx = Lx;
	snapshot->boxLy = Ly;
	snapshot->wallX = addWallx;
	snapshot->wallY = addWally;
	snapshot->middleWall = addMidWall;
	snapshot->circularWall = addCircularWall;
	snapshot->thermostatMode = noise;
	snapshot->thermostatInterval = dtnoise;
	snapshot->time = t;
	snapshot->lastOutputTime = lastScreen;
	snapshot->lastCollisionCount = lastCollNum;
	snapshot->collisions = ncol;
	snapshot->crossings = ncross;
	snapshot->virialX = collTermX;
	snapshot->virialY = collTermY;
	snapshot->virialXY = collTermXY;
	snapshot->virialYX = collTermYX;
	snapshot->wallImpulse = cumulativeWallImpulse;
	snapshot->sampledWallImpulse = dp;
	snapshot->sampledLeftImpulse = dpLeft;
	snapshot->sampledRightImpulse = dpRight;
	snapshot->sampledMiddleImpulse = dpMid;
	#if STRESS
	snapshot->kineticX = sumKineticX;
	snapshot->kineticY = sumKineticY;
	snapshot->kineticXY = sumKineticXY;
	snapshot->kineticYX = sumKineticYX;
	snapshot->runningKineticX = runningKineticX;
	snapshot->runningKineticY = runningKineticY;
	snapshot->runningKineticXY = runningKineticXY;
	snapshot->runningKineticYX = runningKineticYX;
	#endif
	#if HEAT
	snapshot->heatX = sumHeatX;
	snapshot->heatY = sumHeatY;
	snapshot->runningHeatX = runningHeatX;
	snapshot->runningHeatY = runningHeatY;
	snapshot->collisionHeatX = collHeatX;
	snapshot->collisionHeatY = collHeatY;
	snapshot->shearCx2Vy = shearHeatCx2Vy;
	snapshot->shearVy3 = shearHeatVy3;
	snapshot->shearVyVz2 = shearHeatVyVz2;
	snapshot->shearCxVy2 = shearHeatCxVy2;
	snapshot->shearVy2Vz = shearHeatVy2Vz;
	snapshot->shearCxVyVz = shearHeatCxVyVz;
	#endif
	#if STRESS || HEAT
	snapshot->lastTransportTime = lastCollTime;
	#endif
	snapshot->samples = outputSamples;
	memcpy(snapshot->randomState, mt, sizeof(snapshot->randomState));
	snapshot->randomIndex = mti;
	for (int offset = 1; offset <= GRAPHICAL_RECOVERY_AUXILIARY_COUNT;
		offset++){
		const node* source = EVENT_NODE(2*N + offset);
		snapshot->auxiliary[offset - 1] = *source;
		snapshot->auxiliaryQueued[offset - 1] = source->q >= 0;
	}
	snapshot->valid = true;
	graphicalRecoveryNext = (graphicalRecoveryNext + 1) %
		GRAPHICAL_RECOVERY_SNAPSHOT_COUNT;
	if (graphicalRecoveryCount < GRAPHICAL_RECOVERY_SNAPSHOT_COUNT)
		graphicalRecoveryCount++;
	graphicalRecoveryLastWallTime = wallTime;
	graphicalRecoveryLastSimulationTime = t;
	graphicalRecoveryLastCollisionCount = ncol;
}

static void restoreGraphicalRecoveryCalendar(
	const graphicalRecoverySnapshot* snapshot){
	const int totalSize = 2*N + suppSizeEvent;
	memset(eventPaul, 0, (size_t)(paulListN + 1)*sizeof(*eventPaul));
	for (int index = 0; index < totalSize; ++index){
		eventList[index] = EVENT_NODE(index);
		memset(EVENT_NODE(index), 0, sizeof(*EVENT_NODE(index)));
		EVENT_NODE(index)->q = -1;
	}
	root = EVENT_NODE(2*N);
	root->t = never + 1.0;
	root->q = -1;
	treeMin = NULL;
	scheduledEventCount = 0;
	actualPaulList = 0;
	paulTime = t;
	nextEvent = NULL;
	for (int offset = 1; offset <= GRAPHICAL_RECOVERY_AUXILIARY_COUNT;
		offset++){
		if (!snapshot->auxiliaryQueued[offset - 1])
			continue;
		const node* saved = &snapshot->auxiliary[offset - 1];
		node* restored = EVENT_NODE(2*N + offset);
		restored->t = saved->t;
		restored->collActual = saved->collActual;
		restored->i = saved->i;
		restored->j = saved->j;
		restored->type = saved->type;
		addEventToQueue(restored);
	}
	for (int index = 0; index < N; ++index){
		EVENT_NODE(index)->i = index;
		EVENT_NODE(N + index)->i = index;
		crossingEvent(index);
		collisionEvent(index);
		if (simulationHasFault())
			break;
	}
	nextEvent = NULL;
}

static bool restoreGraphicalRecoverySnapshot(
	const graphicalRecoverySnapshot* snapshot){
	if (!graphicalRecoveryConfigurationMatches(snapshot) ||
		snapshot->particleState == NULL || snapshot->polygonState == NULL ||
		snapshot->polygonStateBytes != polygon_particles_checkpoint_bytes() ||
		snapshot->randomIndex < 0 || snapshot->randomIndex > NMT + 1)
		return false;
	if (!polygon_particles_restore_checkpoint(snapshot->polygonState,
		snapshot->polygonStateBytes))
		return false;
	t = snapshot->time;
	for (int index = 0; index < N; ++index){
		free(particles[index].particlesInWell);
		particles[index] = snapshot->particleState[index];
		particles[index].nxt = NULL;
		particles[index].prv = NULL;
		particles[index].particlesInWell = NULL;
		particles[index].numberOfParticlesInWell = 0;
		particles[index].particlesInWellCapacity = 0;
	}
	/* Convert the restored lazy trajectories to snapshot->time only on this
	 * exceptional path.  Their saved cell links describe the later failed
	 * state, so discard those links and force an unconditional ownership
	 * rebuild after synchronization. */
	synchronizeParticlesAtCurrentTime();
	for (int index = 0; index < N; ++index){
		particles[index].cell[0] = -1;
		particles[index].cell[1] = -1;
	}
	simulationClearFault();
	polygon_particles_enable_recovery_precision(1.0);
	rebuildSpatialIndexAtCurrentTime();
	restoreGraphicalRecoveryCalendar(snapshot);
	if (simulationHasFault())
		return false;

	char verification[512];
	const double tolerance = 1.0e-10*(1.0 + fabs(Lx) + fabs(Ly));
	if (!polygon_particles_verify_local_no_overlaps(t, tolerance, verification,
		sizeof(verification)))
		return false;

	ncol = snapshot->collisions;
	ncross = snapshot->crossings;
	collTermX = snapshot->virialX;
	collTermY = snapshot->virialY;
	collTermXY = snapshot->virialXY;
	collTermYX = snapshot->virialYX;
	cumulativeWallImpulse = snapshot->wallImpulse;
	dp = snapshot->sampledWallImpulse;
	dpLeft = snapshot->sampledLeftImpulse;
	dpRight = snapshot->sampledRightImpulse;
	dpMid = snapshot->sampledMiddleImpulse;
	lastScreen = snapshot->lastOutputTime;
	lastCollNum = snapshot->lastCollisionCount;
	#if STRESS
	sumKineticX = snapshot->kineticX;
	sumKineticY = snapshot->kineticY;
	sumKineticXY = snapshot->kineticXY;
	sumKineticYX = snapshot->kineticYX;
	runningKineticX = snapshot->runningKineticX;
	runningKineticY = snapshot->runningKineticY;
	runningKineticXY = snapshot->runningKineticXY;
	runningKineticYX = snapshot->runningKineticYX;
	#endif
	#if HEAT
	sumHeatX = snapshot->heatX;
	sumHeatY = snapshot->heatY;
	runningHeatX = snapshot->runningHeatX;
	runningHeatY = snapshot->runningHeatY;
	collHeatX = snapshot->collisionHeatX;
	collHeatY = snapshot->collisionHeatY;
	shearHeatCx2Vy = snapshot->shearCx2Vy;
	shearHeatVy3 = snapshot->shearVy3;
	shearHeatVyVz2 = snapshot->shearVyVz2;
	shearHeatCxVy2 = snapshot->shearCxVy2;
	shearHeatVy2Vz = snapshot->shearVy2Vz;
	shearHeatCxVyVz = snapshot->shearCxVyVz;
	#endif
	#if STRESS || HEAT
	lastCollTime = snapshot->lastTransportTime;
	#endif
	outputSamples = snapshot->samples;
	memcpy(mt, snapshot->randomState, sizeof(snapshot->randomState));
	mti = snapshot->randomIndex;
	physicalQ();
	nextEvent = NULL;
	return true;
}

static bool tryGraphicalPolygonRecovery(state* screenState){
	(void)screenState;
	if (!simulationFaultIsRecoverableOverlap() ||
		!polygon_particles_active() || polygon_particles_growing() ||
		graphicalRecoveryCount <= 0)
		return false;
	char originalFault[512];
	(void)snprintf(originalFault, sizeof(originalFault), "%s",
		simulationFaultMessage());
	for (int age = 0; age < graphicalRecoveryCount; ++age){
		const int index = (graphicalRecoveryNext +
			GRAPHICAL_RECOVERY_SNAPSHOT_COUNT - 1 - age) %
			GRAPHICAL_RECOVERY_SNAPSHOT_COUNT;
		graphicalRecoverySnapshot* snapshot = &graphicalRecovery[index];
		if (!graphicalRecoveryConfigurationMatches(snapshot))
			continue;
		/* Never replay the same failed checkpoint indefinitely. Older ring
		 * entries remain available if the first replay also proves unsafe. */
		const bool restored = restoreGraphicalRecoverySnapshot(snapshot);
		snapshot->valid = false;
		if (restored){
			simulationClearFault();
			live_analysis_reset_generation();
			graphicalRecoveryLastWallTime = -DBL_MAX;
			graphicalRecoveryLastSimulationTime = -DBL_MAX;
			graphicalRecoveryLastCollisionCount = ncol;
			return true;
		}
	}
	simulationClearFault();
	simulationReportRecoverableOverlapFault(originalFault);
	return false;
}

static void releaseGraphicalRecoveryCheckpoints(void){
	for (int index = 0; index < GRAPHICAL_RECOVERY_SNAPSHOT_COUNT; ++index){
		free(graphicalRecovery[index].particleState);
		free(graphicalRecovery[index].polygonState);
		memset(&graphicalRecovery[index], 0,
			sizeof(graphicalRecovery[index]));
	}
	graphicalRecoveryNext = 0;
	graphicalRecoveryCount = 0;
	graphicalRecoveryLastWallTime = -DBL_MAX;
	graphicalRecoveryLastSimulationTime = -DBL_MAX;
	graphicalRecoveryLastCollisionCount = 0;
}
#endif

/* -------------------------------/
/								  /
/	Event handling utilitaries 	  /
/								  /
/--------------------------------*/

/* ---------------------------------------------
	Finds the next crossing of particle i by
	calculating which cell border will be
	reached first
--------------------------------------------- */



void crossingEventGrow(int i){

	const particle* p = particles + i;
	double travelTimeX = never;
	double travelTimeY = never;
	int xx = 2;
	int yy = 1;
	

	// xx and yy indicates the direction to the new cell for example xx = 1 means that we will go from cell Cell[X][Y] to Cell[X - 1][Y]
	//printf("p: %d %lf |||", p.num, p.vx);
	if (p->vx < 0){
		xx = 1;
		if (!((addWallx || addCircularWall) && p->cell[0] == 0))
			travelTimeX = PBCinsideCellX(
				p->cell[0]*cellxSize - p->x)/p->vx;
	}
	else if (p->vx > 0){
		xx = 2;
		if (!((addWallx || addCircularWall) &&
			p->cell[0] == Nxcells - 1))
			travelTimeX = PBCinsideCellX(
				(1 + p->cell[0])*cellxSize - p->x)/p->vx;
	}

		

	if (p->vy < 0){
		yy = 3;
		if (!((addWally || addCircularWall) && p->cell[1] == 0))
			travelTimeY = PBCinsideCellY(
				p->cell[1]*cellySize - p->y)/p->vy;
	}
	else if (p->vy > 0){
		yy = 4;
		if (!((addWally || addCircularWall) &&
			p->cell[1] == Nycells - 1))
			travelTimeY = PBCinsideCellY(
				(1 + p->cell[1])*cellySize - p->y)/p->vy;
	}
	
	#if THREE_D
	double travelTimeZ = never;
	int zz = 6;
	if (p->vz < 0){
		zz = 5;
		if (!((addWallz || addCircularWall) && p->cell[2] == 0))
			travelTimeZ = PBCinsideCellZ(
				p->cell[2]*cellzSize - p->z)/p->vz;
	}
	else if (p->vz > 0){
		zz = 6;
		if (!((addWallz || addCircularWall) &&
			p->cell[2] == Nzcells - 1))
			travelTimeZ = PBCinsideCellZ(
				(1 + p->cell[2])*cellzSize - p->z)/p->vz;
	}
	if (travelTimeX <= travelTimeY && travelTimeX <= travelTimeZ) {
   		addCrossingEvent(i, xx, t + travelTimeX);
	} 
	else if (travelTimeY <= travelTimeX && travelTimeY <= travelTimeZ) {
		addCrossingEvent(i, yy, t + travelTimeY);
	} 
	else {
		addCrossingEvent(i, zz, t + travelTimeZ);
	}
	#else
	if (travelTimeX < travelTimeY)
		addCrossingEvent(i, xx, t + travelTimeX);
	else
		addCrossingEvent(i, yy, t + travelTimeY);
	#endif

}

void crossingEventNormal(int i){

	const particle* p = particles + i;
	double travelTimeX = never;
	double travelTimeY = never;

	int xx = 2;
	int yy = 1;

	// xx and yy indicates the direction to the new cell for example xx = 1 means that we will go from cell Cell[X][Y] to Cell[X - 1][Y]
	//printf("p: %d %lf |||", p.num, p.vx);
	if (p->vx < 0){
		xx = 1;
		if (!((addWallx || addCircularWall) && p->cell[0] == 0))
			travelTimeX = logTime(PBCinsideCellX(
				p->cell[0]*cellxSize - p->x)/p->vx);
	}
	else if (p->vx > 0){
		xx = 2;
		if (!((addWallx || addCircularWall) &&
			p->cell[0] == Nxcells - 1))
			travelTimeX = logTime(PBCinsideCellX(
				(1 + p->cell[0])*cellxSize - p->x)/p->vx);
	}


	if (addField){
		double vy = p->vy;
		double dy = PBCinsideCellY(p->y - (1 + p->cell[1])*cellySize);

		if ((dy < 0) && (vy > sqrt(2*field*dy))){
			travelTimeY = (-vy + sqrt(vy*vy - 2*field*dy))/field;
			yy = 4;
			
		}
		else{
			dy = PBCinsideCellY(p->y - p->cell[1]*cellySize);
			travelTimeY = (-vy - sqrt(vy*vy - 2*field*dy))/field;
			yy = 3;
		}	
	}	
	else{
		if (p->vy < 0){
			yy = 3;
			if (!((addWally || addCircularWall) && p->cell[1] == 0))
				travelTimeY = logTime(PBCinsideCellY(
					p->cell[1]*cellySize - p->y)/p->vy);
		}
		else if (p->vy > 0){
			yy = 4;
			if (!((addWally || addCircularWall) &&
				p->cell[1] == Nycells - 1))
				travelTimeY = logTime(PBCinsideCellY(
					(1 + p->cell[1])*cellySize - p->y)/p->vy);
		}
	}

	#if THREE_D
	double travelTimeZ = never;
	int zz = 6;
	if (p->vz < 0){
		zz = 5;
		if (!((addWallz || addCircularWall) && p->cell[2] == 0))
			travelTimeZ = logTime(PBCinsideCellZ(
				p->cell[2]*cellzSize - p->z)/p->vz);
	}
	else if (p->vz > 0){
		zz = 6;
		if (!((addWallz || addCircularWall) &&
			p->cell[2] == Nzcells - 1))
			travelTimeZ = logTime(PBCinsideCellZ(
				(1 + p->cell[2])*cellzSize - p->z)/p->vz);
	}
	if (travelTimeX <= travelTimeY && travelTimeX <= travelTimeZ) {
		addCrossingEvent(i, xx, t + travelTimeX);
	} 
	else if (travelTimeY <= travelTimeX && travelTimeY <= travelTimeZ) {
		addCrossingEvent(i, yy, t + travelTimeY);
	} 
	else {
		addCrossingEvent(i, zz, t + travelTimeZ);
	}
	#else

	//printf("\n Now checking cellList: %lf %lf\n", travelTimeX, travelTimeY);
	if (travelTimeX < travelTimeY)
		addCrossingEvent(i, xx, t + travelTimeX);
	else
		addCrossingEvent(i, yy, t + travelTimeY);
	#endif
}
/* ---------------------------------------------
	Creates the crossing event and adds it
	to the BST.
--------------------------------------------- */
void addCrossingEvent(int i, int info, double tCross){
	node* toAdd;

	toAdd = EVENT_NODE(i);
	toAdd->j = info;
	toAdd->type = CELLCROSS;
	toAdd->t = tCross;

	addEventToQueue(toAdd);
}

/* -----------------------------------
	Performs the crossing event
----------------------------------- */
static inline
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
void relocateParticleAfterCrossing(int i, particle* p){
	//updates the cell in which the particle is in after the cellcrossing
	removeFromCell(i);
	int newCell;

	switch (nextEvent->j){
		case 1:
			newCell = p->cell[0] - 1;
				if (newCell == -1){
				p->cell[0] = Nxcells - 1;
				p->crossX -= 1;
			}
			else
				p->cell[0] = newCell;
			break;
		case 2:
			newCell = p->cell[0] + 1;
			if (newCell == Nxcells){
				p->cell[0] = 0;
				p->crossX += 1;
			}
			else
				p->cell[0] = newCell;
			break;
		case 3:
			newCell = p->cell[1] - 1;
			if (newCell == -1){
				p->cell[1] = Nycells - 1;
				p->crossY -= 1;
				if (addShear && t >= postGrowthDynamicsStartTime()){
					/* The lower image is represented at y=Ly while it traverses
					 * the last cell.  Coupled with the -strainRate*Ly velocity
					 * remap, this keeps peculiar vx continuous across the seam. */
					p->y = Ly;
					shearCorrection(&(p->x), &(p->vx), -1, p);
						/* The image-frame velocity jump invalidates collision events
						 * owned by other particles that reference this particle. */
						p->coll++;
						p->cell[0] = coordToCell(p->x, 1);
				}
			}
			else
				p->cell[1] = newCell;
			break;
		case 4:
			newCell = p->cell[1] + 1;
			if (newCell == Nycells){
				p->cell[1] = 0;
				p->crossY += 1;
					if (addShear && t >= postGrowthDynamicsStartTime()){
						shearCorrection(&(p->x), &(p->vx), 1, p);
						p->coll++;
						p->cell[0] = coordToCell(p->x, 1);
				}
			}
			else
				p->cell[1] = newCell;
			break;
		#if THREE_D
		case 5:
			newCell = p->cell[2] - 1;
			if (newCell == -1){
				p->cell[2] = Nzcells - 1;
				p->crossZ -= 1;
			}
			else
				p->cell[2] = newCell;
			break;
		case 6:
			newCell = p->cell[2] + 1;
			if (newCell == Nzcells){
				p->cell[2] = 0;
				p->crossZ += 1;
			}
			else
				p->cell[2] = newCell;
			break;
		#endif
	}

	p->prv = NULL;
	#if THREE_D
	p->nxt = cellList[p->cell[2]*Nycells*Nxcells + p->cell[1]*Nxcells + p->cell[0]];
	cellList[p->cell[2]*Nycells*Nxcells + p->cell[1]*Nxcells + p->cell[0]] = p;
	#else
	p->nxt = cellList[p->cell[1]*Nxcells + p->cell[0]];
	cellList[p->cell[1]*Nxcells + p->cell[0]] = p;
	#endif
	if (p->nxt != NULL)
		p->nxt->prv = p;
}

void doTheCrossing(){
	int i = nextEvent->i;
	particle* p = particles + i;
	const int normalPhase = (freeFly == &freeFlyNormal);
	if (normalPhase)
		freeFlyNormal(p);
	else
		freeFlyGrow(p);

	relocateParticleAfterCrossing(i, p);

	if (normalPhase){
		crossingEventNormal(i);
		removeEventFromQueue(EVENT_NODE(N + i));
		collisionEventNormal(i);
	}
	else{
		crossingEventGrow(i);
		removeEventFromQueue(EVENT_NODE(N + i));
		collisionEventGrow(i);
	}
}

#if POLYGON_PARTICLES && !THREE_D
void doTheCrossingPolygon(void){
	const int i = nextEvent->i;
	particle* p = particles + i;
	polygon_particles_free_fly(p);
	relocateParticleAfterCrossing(i, p);
	polygon_particles_crossing_event(i);
	removeEventFromQueue(EVENT_NODE(N + i));
	polygon_particles_collision_event(i);
}
#endif



/* --------------------------------------
	Calculates the time of collision
	of p1 and p2
--------------------------------------- */
particlePrediction predictParticle(const particle* p, double targetTime,
	enum particleFlightPhysics physics, double parameter){
	const double dt = targetTime - p->t;
	particlePrediction prediction = {
		.x = p->x,
		.y = p->y,
		.vx = p->vx,
		.vy = p->vy,
		.rad = p->rad,
		#if THREE_D
		.z = p->z,
		.vz = p->vz,
		#endif
	};

	if (physics == CONSTANT_FIELD_FLIGHT){
		prediction.x += dt*prediction.vx;
		prediction.y += dt*prediction.vy + 0.5*parameter*dt*dt;
		prediction.vy += dt*parameter;
		#if THREE_D
		prediction.z += dt*prediction.vz;
		#endif
	}
	else if (physics == LINEAR_DAMPING_FLIGHT){
		const double decay = exp(-parameter*dt);
		prediction.x += prediction.vx/parameter*(1 - decay);
		prediction.y += prediction.vy/parameter*(1 - decay);
		prediction.vx *= decay;
		prediction.vy *= decay;
		#if THREE_D
		prediction.z += prediction.vz/parameter*(1 - decay);
		prediction.vz *= decay;
		#endif
	}
	else if (physics == GROWTH_FLIGHT){
		prediction.x += dt*prediction.vx;
		prediction.y += dt*prediction.vy;
		prediction.rad += dt*p->vr;
		#if THREE_D
		prediction.z += dt*prediction.vz;
		#endif
	}
	else{
		prediction.x += dt*prediction.vx;
		prediction.y += dt*prediction.vy;
		#if THREE_D
		prediction.z += dt*prediction.vz;
		#endif
	}

	return prediction;
}

particlePrediction predictFreeFlyNormal(const particle* p, double targetTime){
	if (addField == 1){
		return predictParticle(p, targetTime, CONSTANT_FIELD_FLIGHT, field);
	}
	if (damping == 1){
		return predictParticle(p, targetTime, LINEAR_DAMPING_FLIGHT, gamm);
	}
	return predictParticle(p, targetTime, BALLISTIC_FLIGHT, 0);
}

particlePrediction predictFreeFlyGrow(const particle* p, double targetTime){
	return predictParticle(p, targetTime, GROWTH_FLIGHT, 0);
}

double simulationCurrentPackingFraction(double queryTime){
	if (particles == NULL || N <= 0 || !isfinite(queryTime))
		return NAN;
	const bool growing = freeFly == &freeFlyGrow
		#if POLYGON_PARTICLES && !THREE_D
		|| polygon_particles_growing()
		#endif
		;
	double radiusPowerSum = 0.0;
	double compensation = 0.0;
	for (int index = 0; index < N; index++){
		double radius = particles[index].rad;
		if (growing)
			radius = predictFreeFlyGrow(particles + index, queryTime).rad;
		if (!isfinite(radius) || radius < 0.0)
			return NAN;
		#if THREE_D
		const double contribution = radius*radius*radius;
		#else
		const double contribution = radius*radius;
		#endif
		/* Compensated accumulation keeps the displayed value stable for broad
		 * polydispersity and large graphical populations. */
		const double adjusted = contribution - compensation;
		const double updated = radiusPowerSum + adjusted;
		compensation = (updated - radiusPowerSum) - adjusted;
		radiusPowerSum = updated;
	}

	#if THREE_D
	if (addCircularWall){
		if (!(halfLx > 0.0) || !isfinite(halfLx)) return NAN;
		return radiusPowerSum/(halfLx*halfLx*halfLx);
	}
	const double boxVolume = Lx*Ly*Lz;
	if (!(boxVolume > 0.0) || !isfinite(boxVolume)) return NAN;
	return ((4.0/3.0)*M_PI)*radiusPowerSum/boxVolume;
	#else
	const double containerArea = addCircularWall
		? M_PI*halfLx*halfLx : Lx*Ly;
	#if POLYGON_PARTICLES
	SimParticleShapeConfig particleShape;
	polygon_particles_get_config(&particleShape);
	return simulation_particle_shape_packing_fraction_2d(
		&particleShape, radiusPowerSum, containerArea);
	#else
	if (!(containerArea > 0.0) || !isfinite(containerArea)) return NAN;
	return M_PI*radiusPowerSum/containerArea;
	#endif
	#endif
}

double collisionTimeGrow(const particle* p1, const particle* p2){

	double lat2 = t - p2->t;
	
	double dvx = p2->vx - p1->vx;
	double dvy = p2->vy - p1->vy;
	double dvr = p1->vr + p2->vr;
	
	double dx = (p2->x + lat2*p2->vx) - p1->x;
	double dy = (p2->y + lat2*p2->vy) - p1->y;


	double dr = DIST(p1->rad, p2->rad + lat2*p2->vr);

	#if THREE_D
	double dz = (p2->z + lat2*p2->vz) - p1->z;
	double dvz = p2->vz - p1->vz;
	PBC(&dx, &dy, &dz);
	double b = dx*dvx + dy*dvy + dz*dvz - dvr*dr;
	double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
	double distOfSquare = dx*dx + dy*dy + dz*dz;
	#else
	PBC(&dx, &dy);
	double b = dx*dvx + dy*dvy - dvr*dr;
	double v2 = dvx*dvx + dvy*dvy;
	double distOfSquare = dx*dx + dy*dy;
	#endif
	



	double det = b*b - (v2 - dvr*dvr)*(distOfSquare - dr*dr);




	if (det < 0)
		return never;
 
	double plus = (-b + sqrt(det))/(v2 - dvr*dvr);
	double minus = (-b - sqrt(det))/(v2 - dvr*dvr);
	
	//printf("%.10lf %.10lf\n", plus, minus);
		
	if (((minus > 0) && (plus > 0) && (minus < plus)) || ((minus > 0) && (plus < 0))){
		return minus;
	}
	//hacky
	else if (((minus > 0) && (plus > 0.000000001) && (plus < minus)) || ((minus < 0) && (plus > 0.000000001))){
		return plus;
	}
	if (distOfSquare - dr*dr < -0.01)
		return -never;
	
	return never;

	

}

typedef struct normalPairState normalPairState;
struct normalPairState{
	double dx, dy, dvx, dvy;
	#if THREE_D
	double dz, dvz;
	#endif
	double distanceSquared;
	double relativeSpeedSquared;
	double radialVelocity;
};

static normalPairState normalPairStateAtCurrentTime(const particle* p1,
	const particle* p2){
	const particlePrediction p2Now = predictFreeFlyNormal(p2, t);
	normalPairState state = {
		.dvx = p2Now.vx - p1->vx,
		.dvy = p2Now.vy - p1->vy,
		.dx = p2Now.x - p1->x,
		.dy = p2Now.y - p1->y,
		#if THREE_D
		.dvz = p2Now.vz - p1->vz,
		.dz = p2Now.z - p1->z,
		#endif
	};

	if (addShear){
		correctDistances(p1, p2, p2Now.x, p2Now.vx, &state.dx, &state.dvx);
	}

	#if THREE_D
	PBC(&state.dx, &state.dy, &state.dz);
	state.distanceSquared = state.dx*state.dx + state.dy*state.dy + state.dz*state.dz;
	state.relativeSpeedSquared = state.dvx*state.dvx + state.dvy*state.dvy + state.dvz*state.dvz;
	state.radialVelocity = state.dx*state.dvx + state.dy*state.dvy + state.dz*state.dvz;
	#else
	PBC(&state.dx, &state.dy);
	state.distanceSquared = state.dx*state.dx + state.dy*state.dy;
	state.relativeSpeedSquared = state.dvx*state.dvx + state.dvy*state.dvy;
	state.radialVelocity = state.dx*state.dvx + state.dy*state.dvy;
	#endif
	return state;
}

static double collisionTimeFromPairState(const normalPairState* state,
	double contactDistanceSquared){
	double gapSquared = state->distanceSquared - contactDistanceSquared;
	if (gapSquared < -0.01)
		return -never;
	//no collision.
	if (state->radialVelocity > 0 || state->relativeSpeedSquared <= 1e-30)
		return never;
	/* Roundoff at a just-committed contact can leave a minute negative gap.
	 * The quadratic's entry root is then negative even though the pair is
	 * approaching, so returning `never` lets that harmless numerical overlap
	 * grow into a real missed collision. Commit the contact immediately. */
	if (gapSquared <= 0)
		return 0;
	double det = state->radialVelocity*state->radialVelocity -
		state->relativeSpeedSquared*gapSquared;
	if (det < 0)
		return never;
	/* q/(sqrt(b^2-v^2 q)-b) is the same smaller root as
	 * (-b-sqrt(...))/v^2, but avoids catastrophic cancellation as q -> 0. */
	double displacementTime = gapSquared/
		(sqrt(det) - state->radialVelocity);
	if (displacementTime < 0 && displacementTime > -1e-12)
		displacementTime = 0;
	return displacementTime >= 0 ? logTime(displacementTime) : never;
}

static double wellEntryTimeFromPairState(const normalPairState* state,
	double wellDistanceSquared){
	if (state->radialVelocity > 0 || state->relativeSpeedSquared <= 1e-30)
		return never;
	double gapSquared = state->distanceSquared - wellDistanceSquared;
	double det = state->radialVelocity*state->radialVelocity -
		state->relativeSpeedSquared*gapSquared;
	if (det < 0)
		return never;
	/* Rationalized entry root. The direct (-b-sqrt(det))/v^2 expression
	 * loses most of its precision for a pair just outside the shell. */
	double denominator = sqrt(det) - state->radialVelocity;
	if (denominator <= 0)
		return never;
	double displacementTime = gapSquared/denominator;
	if (displacementTime < 0 && displacementTime > -1e-12)
		displacementTime = 0;
	return displacementTime >= 0 ? logTime(displacementTime) : never;
}

static double wellExitTimeFromPairState(const normalPairState* state,
	double wellDistanceSquared){
	if (state->relativeSpeedSquared <= 1e-30)
		return never;
	double gapSquared = state->distanceSquared - wellDistanceSquared;
	double det = state->radialVelocity*state->radialVelocity -
		state->relativeSpeedSquared*gapSquared;
	if (det < 0)
		return never;
	double displacementTime = (-state->radialVelocity + sqrt(det))/
		state->relativeSpeedSquared;
	return displacementTime >= 0 ? logTime(displacementTime) : never;
}

double collisionTimeNormal(const particle* p1, const particle* p2){
	/* Hard-particle scheduling dominates dense disk runs.  Build only the
	 * quantities needed by each rejection stage: more than half of ordinary
	 * linked-cell candidates are already moving apart, so calculating their
	 * relative-speed norm and discriminant eagerly is wasted work.  Keep the
	 * overlap check before the radial rejection so a separating invalid state
	 * is still reported exactly as it was by collisionTimeFromPairState(). */
	const particlePrediction p2Now = predictFreeFlyNormal(p2, t);
	double dx = p2Now.x - p1->x;
	double dy = p2Now.y - p1->y;
	double dvx = p2Now.vx - p1->vx;
	double dvy = p2Now.vy - p1->vy;
	#if THREE_D
	double dz = p2Now.z - p1->z;
	double dvz = p2Now.vz - p1->vz;
	#endif

	if (addShear){
		correctDistances(p1, p2, p2Now.x, p2Now.vx, &dx, &dvx);
	}

	#if THREE_D
	PBC(&dx, &dy, &dz);
	const double distanceSquared = dx*dx + dy*dy + dz*dz;
	const double radialVelocity = dx*dvx + dy*dvy + dz*dvz;
	#else
	PBC(&dx, &dy);
	const double distanceSquared = dx*dx + dy*dy;
	const double radialVelocity = dx*dvx + dy*dvy;
	#endif
	const double gapSquared = distanceSquared - DIST2(p1->rad, p2->rad);
	if (gapSquared < -0.01)
		return -never;
	if (radialVelocity > 0)
		return never;

	#if THREE_D
	const double relativeSpeedSquared =
		dvx*dvx + dvy*dvy + dvz*dvz;
	#else
	const double relativeSpeedSquared = dvx*dvx + dvy*dvy;
	#endif
	if (relativeSpeedSquared <= 1e-30)
		return never;
	if (gapSquared <= 0)
		return 0;
	const double det = radialVelocity*radialVelocity -
		relativeSpeedSquared*gapSquared;
	if (det < 0)
		return never;
	double displacementTime = gapSquared/
		(sqrt(det) - radialVelocity);
	if (displacementTime < 0 && displacementTime > -1e-12)
		displacementTime = 0;
	return displacementTime >= 0 ? logTime(displacementTime) : never;

}

static double scheduledCollisionTimeNormal(const particle* p1, const particle* p2){
	const double prediction = collisionTimeNormal(p1, p2);
	if (prediction != -never)
		return prediction;
	if (bubbleStop){
		#if G
		simulationSetFault(
			"The bubble/overlap stopping criterion was reached at t=%.9g. "
			"Event processing is locked; adjust the configuration and restart the simulation.",
			t);
		return never;
		#else
		printf("\nBubble detected! Stopping simulation...\n");
		if (thermo != NULL){
			fprintf(thermo, "-1\n");
		}
		exit(4);
		#endif
	}
	#if G
	simulationSetFault(
		"Particles %d and %d overlap at t=%.9g. The simulation was paused "
		"before processing more events. Correct the setup and restart to build a fresh event queue.",
		p1->num, p2->num, t);
	return never;
	#else
	printf("\nERROR:\033[0;31m Overlap detected during SIMULATION between particle %d and particle %d.\033[0m\n",
		p1->num, p2->num);
	saveTXT();
	exit(3);
	#endif
}

static double scheduledCollisionTimeGrow(const particle* p1, const particle* p2){
	const double prediction = collisionTimeGrow(p1, p2);
	if (prediction != -never)
		return prediction;
	#if G
	simulationSetFault(
		"Particles %d and %d overlap during radius growth at t=%.9g. "
		"The simulation was paused and locked; reduce density/growth or choose a new initial state, then restart.",
		p1->num, p2->num, t);
	return never;
	#else
	printf("\nERROR:\033[0;31m Overlap detected during GROWTH between particle %d and particle %d.\033[0m\n",
		p1->num, p2->num);
	exit(3);
	#endif
}

#if !THREE_D
/* Scan one linked-cell head for the narrow, no-potential 2D disk path.  Keeping
 * this helper inline lets collisionEventNormal spell out its fixed 3x3
 * stencil without rebuilding coordinate arrays or carrying variable trip
 * counts through the hottest scheduler loop. */
static inline int scanHardDiskCellNormal(const particle* p1, particle* p2,
	double* bestTime, int* bestPartner){
	int improved = 0;
	while (p2 != NULL){
		if (p1 != p2){
			const double candidateTime =
				scheduledCollisionTimeNormal(p1, p2);
			if (*bestTime > candidateTime){
				*bestTime = candidateTime;
				*bestPartner = p2->num;
				improved = 1;
			}
		}
		p2 = p2->nxt;
	}
	return improved;
}
#endif

//time to enter in the square well: https://faculty.biu.ac.il/~rapaport/papers/09b-ptp.pdf
double sphereTime(const particle* p1, const particle* p2){
	const normalPairState state = normalPairStateAtCurrentTime(p1, p2);
	double wellDistance = sig*(p1->rad + p2->rad);
	return wellEntryTimeFromPairState(&state, wellDistance*wellDistance);
}


//time to leave the square well: https://faculty.biu.ac.il/~rapaport/papers/09b-ptp.pdf
double separateTime(const particle* p1, const particle* p2){
	const normalPairState state = normalPairStateAtCurrentTime(p1, p2);
	double wellDistance = sig*(p1->rad + p2->rad);
	return wellExitTimeFromPairState(&state, wellDistance*wellDistance);
}

static double scheduledPotentialEventTime(const particle* p1, const particle* p2,
	int insideWell, enum event* eventType){
	/* A zero-height step is exactly a hard-particle system.  Avoid generating
	 * redundant IN/OUT events (and the former same-time IN/OUT churn) while
	 * retaining the user's potential setting in the configuration/output name. */
	if (U == 0.0 || potentialPairStepEnergy(p1, p2) == 0.0){
		*eventType = COLLISION;
		return scheduledCollisionTimeNormal(p1, p2);
	}
	const normalPairState state = normalPairStateAtCurrentTime(p1, p2);
	double wellDistance = sig*(p1->rad + p2->rad);
	if (!insideWell){
		*eventType = IN;
		return wellEntryTimeFromPairState(&state, wellDistance*wellDistance);
	}

	double collisionTime = collisionTimeFromPairState(&state, DIST2(p1->rad, p2->rad));
	if (collisionTime == -never){
		if (bubbleStop){
			#if G
			simulationSetFault(
				"The bubble/overlap stopping criterion was reached at t=%.9g. "
				"Adjust the configuration and restart the simulation.", t);
			return never;
			#else
			printf("\nBubble detected! Stopping simulation...\n");
			if (thermo != NULL){
				fprintf(thermo, "-1\n");
			}
			exit(4);
			#endif
		}
		#if G
		simulationSetFault(
			"Particles %d and %d overlap inside the active potential at t=%.9g. "
			"The simulation was paused and locked; correct the setup and restart.",
			p1->num, p2->num, t);
		return never;
		#else
		printf("\nERROR:\033[0;31m Overlap detected during SIMULATION between particle %d and particle %d.\033[0m\n",
			p1->num, p2->num);
		saveTXT();
		exit(3);
		#endif
	}
	double exitTime = wellExitTimeFromPairState(&state, wellDistance*wellDistance);
	if (collisionTime < exitTime){
		*eventType = COLLISION;
		return collisionTime;
	}
	*eventType = OUT;
	return exitTime;
}
/* ---------------------------------------------
	Finds the next collision of particle i
	by calculating collision time of i with
	its surrouding (using the cellList)
--------------------------------------------- */
#if INTRUDER
void collisionEventIntruder(){
	double dt = never;
	int index = 0;

	for (int i = 0; i < N; i++){
		const particle* p1 = particles + i;
		double time = collisionTimeIntruderAt(&intru, p1, t);
		if ((time < dt) && (time > INTRUDER_TOL)){
			dt = time;
			index = i;
		}
	}
	addCollisionIntruderEvent(-1, index, t + dt);
}
#endif

static double smallestPositiveQuadraticExitRoot(double A, double B, double C,
	double surfaceTolerance){
	/* Solve A*t^2 + B*t + C = 0 for the first transition from C <= 0
	 * (inside) to C > 0 (outside).  Surface classification is deliberately
	 * supplied by the geometry.  Folding velocity-sized A/B coefficients into
	 * that tolerance used to classify deeply interior high-speed particles as
	 * already touching the wall. */
	const bool onSurface = fabs(C) <= surfaceTolerance;
	if (C > surfaceTolerance)
		return 0;

	const double beta = 0.5*B;
	if (onSurface){
		if (beta > 0)
			return 0;
		if (beta == 0)
			return A > 0 ? 0 : never;
		/* The zero root is a re-entry point.  A positive quadratic term is
		 * required for a later exit on the far side. */
		return A > 0 ? -2.0*beta/A : never;
	}

	if (A == 0)
		return beta > 0 ? -C/(2.0*beta) : never;

	double discriminant = beta*beta - A*C;
	if (discriminant < 0){
		const double discriminantTolerance = 64.0*DBL_EPSILON*
			(beta*beta + fabs(A*C) + DBL_MIN);
		if (discriminant >= -discriminantTolerance)
			discriminant = 0;
		else
			return never;
	}
	const double squareRoot = sqrt(discriminant);
	if (beta > 0){
		/* Rationalized outward root: accurate even an ulp below the wall. */
		const double denominator = squareRoot + beta;
		return denominator > 0 ? -C/denominator : never;
	}
	if (A > 0)
		return (squareRoot - beta)/A;
	return never;
}

double circularWallTimeNormal(const particle* p){
	double rx = p->x - halfLx;
	double ry = p->y - halfLy;
	double radius = halfLx - p->rad;
	double positionSquared = rx*rx + ry*ry;
	double speedSquared = p->vx*p->vx + p->vy*p->vy;
	double radialVelocity = rx*p->vx + ry*p->vy;
	#if THREE_D
	double rz = p->z - halfLz;
	positionSquared += rz*rz;
	speedSquared += p->vz*p->vz;
	radialVelocity += rz*p->vz;
	#endif
	double C = positionSquared - radius*radius;
	const double surfaceTolerance = 8.0*DBL_EPSILON*
		(positionSquared + radius*radius + 1.0);
	if (addField){
		if (C > surfaceTolerance)
			return 0;
		double quartic = field*field/4;
		double cubic = field*p->vy;
		double quadratic = speedSquared + field*ry;
		double linear = 2*radialVelocity;
		if (fabs(C) <= surfaceTolerance){
			/* At a surface point, inspect the first non-zero derivative of the
			 * trajectory. A tangent particle at the top of a sphere may accelerate
			 * inward under gravity, while the same state at the bottom accelerates
			 * outward and must collide immediately. */
			double initialDirection = linear != 0 ? linear :
				(quadratic != 0 ? quadratic :
				(cubic != 0 ? cubic : quartic));
			if (initialDirection > 0)
				return 0;
			if (initialDirection == 0)
				return never;
		}
		const double minimumRoot = fabs(C) <= surfaceTolerance ?
			64*DBL_EPSILON : 0;
		double result = smallestRootAfter(quartic, cubic, quadratic,
			linear, C, minimumRoot);
		return isfinite(result) ? result : never;
	}
	double displacementTime = smallestPositiveQuadraticExitRoot(
		speedSquared, 2*radialVelocity, C, surfaceTolerance);
	return displacementTime == never ? never : logTime(displacementTime);
}

double circularWallTimeGrow(const particle* p){
	double rx = p->x - halfLx;
	double ry = p->y - halfLy;
	double radius = halfLx - p->rad;
	double positionSquared = rx*rx + ry*ry;
	double speedSquared = p->vx*p->vx + p->vy*p->vy;
	double radialVelocity = rx*p->vx + ry*p->vy;
	#if THREE_D
	double rz = p->z - halfLz;
	positionSquared += rz*rz;
	speedSquared += p->vz*p->vz;
	radialVelocity += rz*p->vz;
	#endif
	double A = speedSquared - p->vr*p->vr;
	double B = 2*(radialVelocity + radius*p->vr);
	double C = positionSquared - radius*radius;
	const double surfaceTolerance = 8.0*DBL_EPSILON*
		(positionSquared + radius*radius + 1.0);
	return smallestPositiveQuadraticExitRoot(A, B, C, surfaceTolerance);
}

void collisionEventNormal(int i){
	int finalPartner = 0;
	double dt = never;
	enum event type = COLLISION;
	particle* p1 = particles + i;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	#if THREE_D
	int Z = p1->cell[2];
	#endif
	int xy = 2;
	double dtTemp;
	enum event typeTemp;



	if (addWally){
		if (addField){
			if (Y == 0){
				dt = (-p1->vy - sqrt(p1->vy*p1->vy - 2*field*(p1->y - p1->rad)))/field;
				type = WALL;
				xy = 1;
			}
			if ((p1->vy > 0)){
				double dt2 = (-p1->vy + sqrt(p1->vy*p1->vy - 2*field*(p1->y - (Ly - p1->rad))))/field;
				if (!isnan(dt2) && (dt2 < dt)){
					dt = dt2;
					type = WALL;
					xy = 1;
				}
			}
		}
		else{
			if ((Y == 0) && (p1->vy < 0)){
				dt = logTime(-(p1->y - p1->rad)/p1->vy);
				type = WALL;
				xy = 1;
			}
			else if ((Y == Nycells - 1) && (p1->vy > 0)){

				dt = logTime((Ly - p1->y - p1->rad)/p1->vy);
				type = WALL;
				xy = 1;
			}
		}
	}
	if (addWallx){
		if ((X == 0) && (p1->vx < 0)){

			dtTemp = logTime(-(p1->x - p1->rad)/p1->vx);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 0;
			}
		}

		else if ((X == Nxcells - 1) && (p1->vx > 0)){

			dtTemp = logTime((Lx - p1->x - p1->rad)/p1->vx);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 0;
			}
		}
	}
	#if THREE_D
	if (addWallz){
		if ((Z == 0) && (p1->vz < 0)){
			dtTemp = logTime(-(p1->z - p1->rad)/p1->vz);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 3;
			}
		}
		else if ((Z == Nzcells - 1) && (p1->vz > 0)){
			dtTemp = logTime((Lz - p1->z - p1->rad)/p1->vz);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 3;
			}
		}
	}
	#endif

	if ((addMidWall) && (p1->type == 1)){
		if ((p1->x > Lx/2 - 4) && (p1->x < Lx/2) && (p1->vx > 0)){
			dtTemp = logTime((Lx/2 - p1->x - p1->rad)/p1->vx);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 0;
			}
		}
		else if ((p1->x < Lx/2 + 4) && (p1->x > Lx/2) && (p1->vx < 0)){
		dtTemp = logTime(((p1->x - p1->rad) - Lx/2)/(-p1->vx));
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 0;
			}
		}
	}


	if (addCircularWall){
		dt = circularWallTimeNormal(p1);
		type = WALL;
		xy = 2;
	}
	#if THREE_D
	int neighboringZ[5];
	int numberZ = neighboringCellCoordinates(Z, Nzcells, neighboringZ);
	#else
	int numberZ = 1;
	#endif
	int usedFixedDiskStencil = 0;
	#if !THREE_D
	if (!addWell && !addShear && !addField && !addWallx && !addWally &&
		!addWallz &&
		!addMidWall && !addCircularWall && collisionCellNeighborRange == 1 &&
		Nxcells > 2 && Nycells > 2){
		const int left = X == 0 ? Nxcells - 1 : X - 1;
		const int right = X + 1 == Nxcells ? 0 : X + 1;
		const int lowerRow = (Y == 0 ? Nycells - 1 : Y - 1)*Nxcells;
		const int middleRow = Y*Nxcells;
		const int upperRow = (Y + 1 == Nycells ? 0 : Y + 1)*Nxcells;
		int foundCollision = 0;
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[lowerRow + left], &dt, &finalPartner);
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[lowerRow + X], &dt, &finalPartner);
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[lowerRow + right], &dt, &finalPartner);
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[middleRow + left], &dt, &finalPartner);
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[middleRow + X], &dt, &finalPartner);
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[middleRow + right], &dt, &finalPartner);
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[upperRow + left], &dt, &finalPartner);
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[upperRow + X], &dt, &finalPartner);
		foundCollision |= scanHardDiskCellNormal(p1,
			cellList[upperRow + right], &dt, &finalPartner);
		if (foundCollision){
			type = COLLISION;
		}
		usedFixedDiskStencil = 1;
	}
	#endif
	if (!usedFixedDiskStencil){
		int neighboringX[5];
		int neighboringY[5];
		int numberX = neighboringCellCoordinates(X, Nxcells, neighboringX);
		int numberY = neighboringCellCoordinates(Y, Nycells, neighboringY);
		for (int zi = 0; zi < numberZ; zi++){
			for (int yi = 0; yi < numberY; yi++){
				for (int xi = 0; xi < numberX; xi++){
				#if THREE_D
					particle* p2 = cellList[neighboringZ[zi]*Nycells*Nxcells +
						neighboringY[yi]*Nxcells + neighboringX[xi]];
				#else
					particle* p2 = cellList[neighboringY[yi]*Nxcells + neighboringX[xi]];
				#endif
					while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
						if (p1 != p2){
							if (addWell){
								dtTemp = scheduledPotentialEventTime(p1, p2,
									isParticleInWellList(p1, p2->num), &typeTemp);
							}
							else{
								dtTemp = scheduledCollisionTimeNormal(p1, p2);
								typeTemp = COLLISION;
							}
							if (dt > dtTemp){ //get min
								finalPartner = p2->num;
								dt = dtTemp;
								type = typeTemp;
							}
						}
						p2 = p2->nxt;
					}
				}
			}
		}
	}
	if (addShear){
		if (Y == 0){
			for (int zi = 0; zi < numberZ; zi++){
			for (int k = 0; k < Nxcells; k++){
				#if THREE_D
				particle* p2 = cellList[neighboringZ[zi]*Nycells*Nxcells +
					(Nycells - 1)*Nxcells + k];
				#else
				particle* p2 = cellList[(Nycells - 1)*Nxcells + k];
				#endif
				while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (p1 != p2){
						if (addWell){
							dtTemp = scheduledPotentialEventTime(p1, p2,
								isParticleInWellList(p1, p2->num), &typeTemp);
						}
						else{
							dtTemp = scheduledCollisionTimeNormal(p1, p2);
							typeTemp = COLLISION;
						}

						if (dt > dtTemp){ //get min
							finalPartner = p2->num;
							dt = dtTemp;
							type = typeTemp;
						}
					}
					p2 = p2->nxt;
				}
			}
			}
		}
		else if (Y == Nycells - 1){
			for (int zi = 0; zi < numberZ; zi++){
			for (int k = 0; k < Nxcells; k++){
				#if THREE_D
				particle* p2 = cellList[neighboringZ[zi]*Nycells*Nxcells + k];
				#else
				particle* p2 = cellList[k];
				#endif
				while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (p1 != p2){
						if (addWell){
							dtTemp = scheduledPotentialEventTime(p1, p2,
								isParticleInWellList(p1, p2->num), &typeTemp);
						}
						else{
							dtTemp = scheduledCollisionTimeNormal(p1, p2);
							typeTemp = COLLISION;
						}

						if (dt > dtTemp){ //get min
							finalPartner = p2->num;
							dt = dtTemp;
							type = typeTemp;
						}
					}
					p2 = p2->nxt;
				}
			}
			}
		}
	}
	#if INTRUDER
		double time = collisionTimeIntruderAt(&intru, p1, t);
		//printf("%lf\n", time);
		if ((time < dt) && (time > INTRUDER_TOL)){
			addCollisionIntruderEvent(i, -1, t + time);
			return;
		}
		//printf("Done checking intruder\n\n");
	#endif
	
	if (type == COLLISION){
		addCollisionEvent(i, finalPartner, t + dt);
	}
	else if (type == OUT)
		addOutEvent(i, finalPartner, t + dt);
	else if (type == IN)
		addInEvent(i, finalPartner, t + dt);
	else{
		addWallEvent(i, xy, t + dt);
	}

}

void collisionEventGrow(int i){
	int finalPartner = 0;
	double dt = 10000000;
	enum event type = COLLISION;
	particle* p1 = particles + i;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	#if THREE_D
	int Z = p1->cell[2];
	#endif
	int xy = 2;
	double dtTemp;
	enum event typeTemp;


	double vrParticle = p1->vr;
	if (addWally){
		double vyPlus = p1->vy + vrParticle;
		double vyMinus = p1->vy - vrParticle;
		if ((Y == 0) && (vyMinus < 0)){
			dt = -(p1->y - p1->rad)/vyMinus;
			type = WALL;
			xy = 1;
			}
		else if ((Y == Nycells - 1) && (vyPlus > 0)){

			dt = (Ly - p1->y - p1->rad)/vyPlus;
			type = WALL;
			xy = 1;
		}
	}
	if (addWallx){
		double vxPlus = p1->vx + vrParticle;
		double vxMinus = p1->vx - vrParticle;
		if ((X == 0) && (vxMinus < 0)){
			dtTemp = -(p1->x - p1->rad)/vxMinus;
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 0;
			}
		}

		else if ((X == Nxcells - 1) && (vxPlus > 0)){

			dtTemp = (Lx - p1->x - p1->rad)/vxPlus;
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 0;
			}
		}
	}
	#if THREE_D
	if (addWallz){
		double vzPlus = p1->vz + vrParticle;
		double vzMinus = p1->vz - vrParticle;
		if ((Z == 0) && (vzMinus < 0)){
			dtTemp = -(p1->z - p1->rad)/vzMinus;
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 3;
			}
		}
		else if ((Z == Nzcells - 1) && (vzPlus > 0)){
			dtTemp = (Lz - p1->z - p1->rad)/vzPlus;
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 3;
			}
		}
	}
	#endif

	if (addCircularWall){
		dt = circularWallTimeGrow(p1);
		type = WALL;
		xy = 2;
	}


	if (addMidWall && p1->type == 1){
		double vxPlus = p1->vx + vrParticle;
		double vxMinus = p1->vx - vrParticle;
		if ((p1->x > Lx/2 - 4) && (p1->x < Lx/2) && (vxPlus > 0)){
			dtTemp = (Lx/2 - p1->x - p1->rad)/(vxPlus);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 0;
			}
		}
		else if ((p1->x < Lx/2 + 4) && (p1->x > Lx/2) && (vxMinus < 0)){
			dtTemp = ((p1->x - p1->rad) - Lx/2)/(-vxMinus);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
				xy = 0;
			}
		}
	}
	int neighboringX[5];
	int neighboringY[5];
	int numberX = neighboringCellCoordinates(X, Nxcells, neighboringX);
	int numberY = neighboringCellCoordinates(Y, Nycells, neighboringY);
	#if THREE_D
	int neighboringZ[5];
	int numberZ = neighboringCellCoordinates(Z, Nzcells, neighboringZ);
	#else
	int numberZ = 1;
	#endif
	for (int zi = 0; zi < numberZ; zi++){
		for (int yi = 0; yi < numberY; yi++){
			for (int xi = 0; xi < numberX; xi++){
			#if THREE_D
				particle* p2 = cellList[neighboringZ[zi]*Nycells*Nxcells +
					neighboringY[yi]*Nxcells + neighboringX[xi]];
			#else
				particle* p2 = cellList[neighboringY[yi]*Nxcells + neighboringX[xi]];
			#endif
				while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (p1 != p2){
						dtTemp = scheduledCollisionTimeGrow(p1, p2);
						typeTemp = COLLISION;
						

						if (dt > dtTemp){ //get min
							finalPartner = p2->num;
							dt = dtTemp;
							type = typeTemp;
						}
					}
					p2 = p2->nxt;
				}
			}
		}
	}

	if (type == COLLISION){
		addCollisionEvent(i, finalPartner, t + dt);
	}
	else
		addWallEvent(i, xy, t + dt);

}


/* ---------------------------------------------
	Updates the Well List
--------------------------------------------- */
static int particleWellLowerBound(const particle* p, int num){
	int lower = 0;
	int upper = p->numberOfParticlesInWell;
	while (lower < upper){
		const int middle = lower + (upper - lower)/2;
		if (p->particlesInWell[middle] < num)
			lower = middle + 1;
		else
			upper = middle;
	}
	return lower;
}

void addParticleInWellList(particle* p1, int num){
	/* Keep the compact neighbour array ordered.  Collision prediction asks
	 * whether each candidate pair is inside the square well far more often
	 * than IN/OUT events mutate the set, so binary lookup is substantially
	 * cheaper than rescanning a dense 3D neighbour shell every time. */
	const int lower = particleWellLowerBound(p1, num);
	if ((lower < p1->numberOfParticlesInWell) &&
		(p1->particlesInWell[lower] == num))
		return;
	if (p1->numberOfParticlesInWell == p1->particlesInWellCapacity){
		int newCapacity = p1->particlesInWellCapacity ? 2*p1->particlesInWellCapacity : 8;
		int* resized = realloc(p1->particlesInWell, (size_t)newCapacity*sizeof(*resized));
		if (resized == NULL){
			fprintf(stderr, "Unable to grow the potential-neighbor list for particle %d\n", p1->num);
			exit(EXIT_FAILURE);
		}
		p1->particlesInWell = resized;
		p1->particlesInWellCapacity = newCapacity;
	}
	memmove(p1->particlesInWell + lower + 1,
		p1->particlesInWell + lower,
		(size_t)(p1->numberOfParticlesInWell - lower)*
			sizeof(*p1->particlesInWell));
	p1->particlesInWell[lower] = num;
	p1->numberOfParticlesInWell++;
}

void removeParticleInWellList(particle* p1, int num){
	const int lower = particleWellLowerBound(p1, num);
	if ((lower >= p1->numberOfParticlesInWell) ||
		(p1->particlesInWell[lower] != num))
		return;
	p1->numberOfParticlesInWell--;
	memmove(p1->particlesInWell + lower,
		p1->particlesInWell + lower + 1,
		(size_t)(p1->numberOfParticlesInWell - lower)*
			sizeof(*p1->particlesInWell));
}

int isParticleInWellList(const particle* p1, int num){
	/* A short sequential scan is faster and more cache-friendly for the small
	 * shells common in 2D and dilute 3D systems.  Use the ordered representation
	 * only once the shell is large enough for logarithmic lookup to win. */
	if (p1->numberOfParticlesInWell <= 16){
		for (int i = 0; i < p1->numberOfParticlesInWell; i++){
			if (p1->particlesInWell[i] == num)
				return 1;
		}
		return 0;
	}
	const int lower = particleWellLowerBound(p1, num);
	return (lower < p1->numberOfParticlesInWell) &&
		(p1->particlesInWell[lower] == num);
}


/* ---------------------------------------------
	Creates the collision event and adds it
	to the BST.
--------------------------------------------- */
#if INTRUDER
void addCollisionIntruderEvent(int i, int j, double tColl) {
	node* toAdd;

	if ((i >= 0) && (tColl > EVENT_NODE(i)->t))
		return;
	
	if (i == -1){
		toAdd = EVENT_NODE(2*N + 8);
		toAdd->collActual = particles[j].coll;
	}
	else{
		toAdd = EVENT_NODE(N + i);
		toAdd->collActual = intru.coll;
	}
	toAdd->j = j;
	toAdd->type = INTRUDER_COLLISION;
	toAdd->t = tColl;

	addEventToQueue(toAdd);
}
#endif


void addCollisionEvent(int i, int j, double tColl){
	node* toAdd;

	if (tColl > EVENT_NODE(i)->t)
		return;

	toAdd = EVENT_NODE(N + i);
	toAdd->j = j;
	toAdd->type = COLLISION;
	toAdd->t = tColl;
	toAdd->collActual = particles[j].coll;

	addEventToQueue(toAdd);
}

void addInEvent(int i, int j, double tColl){
	node* toAdd;

	if (tColl > EVENT_NODE(i)->t)
		return;

	toAdd = EVENT_NODE(N + i);
	toAdd->j = j;
	toAdd->type = IN;
	toAdd->t = tColl;
	toAdd->collActual = particles[j].coll; //the collision trick of the article
	addEventToQueue(toAdd);
}

void addOutEvent(int i, int j, double tColl){
	node* toAdd;

	if (tColl > EVENT_NODE(i)->t)
		return;

	toAdd = EVENT_NODE(N + i);
	toAdd->j = j;
	toAdd->type = OUT;
	toAdd->t = tColl;
	toAdd->collActual = particles[j].coll; //the collision trick of the article
	addEventToQueue(toAdd);
}

void addWallEvent(int i, int xy, double tColl){
	node* toAdd;

	if (tColl > EVENT_NODE(i)->t)
		return;
	if (!eventTimeCanBeQueued(tColl))
		return;
	/* A zero-radius growing particle can reach the outer cell face and the
	 * physical wall at exactly the same time.  The crossing node is inserted
	 * first, so it used to wrap the particle to the opposite side and delete
	 * the wall collision.  Once the wall is the selected interaction it owns
	 * the trajectory up to that time and will schedule a fresh crossing after
	 * reflection. */
	removeEventFromQueue(EVENT_NODE(i));

	toAdd = EVENT_NODE(N + i);
	toAdd->type = WALL;
	toAdd->t = tColl;
	toAdd->j = xy;
	addEventToQueue(toAdd);
}

static void accumulateCollisionVirial(double impulseX, double impulseY,
	double impulseZ, double separationX, double separationY, double separationZ){
	collTermX += impulseX*separationX;
	collTermY += impulseY*separationY;
	collTermXY += impulseX*separationY;
	collTermYX += impulseY*separationX;
	#if THREE_D
	collTermZ += impulseZ*separationZ;
	collTermXZ += impulseX*separationZ;
	collTermZX += impulseZ*separationX;
	collTermYZ += impulseY*separationZ;
	collTermZY += impulseZ*separationY;
	#else
	(void)impulseZ;
	(void)separationZ;
	#endif
}

static bool shearDynamicsActive(void){
	return addShear && t >= postGrowthDynamicsStartTime();
}

static double runtimeShearStrain(void){
	return shearStrainPhase +
		strainRate*(t - postGrowthDynamicsStartTime());
}

void resetRuntimeShearPhase(void){
	shearStrainPhase = 0.0;
}

void rebaseRuntimeShear(bool wasActive, double oldRate,
	bool isActive, double newRate){
	const double start = postGrowthDynamicsStartTime();
	if (t < start){
		/* Growth is unsheared.  Keep the configured shear origin at the exact
		 * transition to normal dynamics. */
		shearStrainPhase = 0.0;
		return;
	}
	const double strainNow = wasActive
		? shearStrainPhase + oldRate*(t - start)
		: 0.0;
	shearStrainPhase = isActive
		? strainNow - newRate*(t - start)
		: 0.0;
}

static double shearPeculiarVelocityX(double storedVx, double y){
	return shearDynamicsActive() ?
		storedVx + strainRate*(y - halfLy) : storedVx;
}

#if (STRESS || HEAT)
static double transportVelocityX(const particle* p, double storedVx){
	return shearPeculiarVelocityX(storedVx, p->y);
}

static void advanceTransportAccumulators(){
	double interval = t - lastCollTime;
	if (interval < 0)
		interval = 0;
	#if STRESS
	if (shearDynamicsActive()){
		/* Between events the peculiar x velocity is
		 * c_x(s)=c_x(0)+strainRate*v_y*s.  Integrate its polynomial
		 * moments exactly instead of treating them as piecewise constant. */
		const double dt2 = interval*interval;
		const double dt3 = dt2*interval;
		const double rate = strainRate;
		const double oldX = sumKineticX;
		const double oldXY = sumKineticXY;
		const double oldYX = sumKineticYX;
		runningKineticX += oldX*interval + rate*oldXY*dt2 +
			rate*rate*sumKineticY*dt3/3.0;
		runningKineticY += sumKineticY*interval;
		runningKineticXY += oldXY*interval +
			0.5*rate*sumKineticY*dt2;
		runningKineticYX += oldYX*interval +
			0.5*rate*sumKineticY*dt2;
		sumKineticX = oldX + 2.0*rate*oldXY*interval +
			rate*rate*sumKineticY*dt2;
		sumKineticXY = oldXY + rate*sumKineticY*interval;
		sumKineticYX = oldYX + rate*sumKineticY*interval;
		#if THREE_D
		const double oldXZ = sumKineticXZ;
		const double oldZX = sumKineticZX;
		runningKineticZ += sumKineticZ*interval;
		runningKineticXZ += oldXZ*interval +
			0.5*rate*sumKineticYZ*dt2;
		runningKineticZX += oldZX*interval +
			0.5*rate*sumKineticZY*dt2;
		runningKineticYZ += sumKineticYZ*interval;
		runningKineticZY += sumKineticZY*interval;
		sumKineticXZ = oldXZ + rate*sumKineticYZ*interval;
		sumKineticZX = oldZX + rate*sumKineticZY*interval;
		#endif
	}
	else{
		runningKineticX += sumKineticX*interval;
		runningKineticY += sumKineticY*interval;
		runningKineticXY += sumKineticXY*interval;
		runningKineticYX += sumKineticYX*interval;
	#if THREE_D
		runningKineticZ += sumKineticZ*interval;
		runningKineticXZ += sumKineticXZ*interval;
		runningKineticZX += sumKineticZX*interval;
		runningKineticYZ += sumKineticYZ*interval;
		runningKineticZY += sumKineticZY*interval;
	#endif
	}
	#endif
	#if HEAT
	if (shearDynamicsActive()){
		const double dt2 = interval*interval;
		const double dt3 = dt2*interval;
		const double dt4 = dt3*interval;
		const double rate = strainRate;
		const double rate2 = rate*rate;
		const double coefficientX1 = 0.5*rate*
			(3.0*shearHeatCx2Vy + shearHeatVy3 + shearHeatVyVz2);
		const double coefficientX2 = 1.5*rate2*shearHeatCxVy2;
		const double coefficientX3 = 0.5*rate2*rate*shearHeatVy3;
		const double coefficientY1 = rate*shearHeatCxVy2;
		const double coefficientY2 = 0.5*rate2*shearHeatVy3;
		runningHeatX += sumHeatX*interval + 0.5*coefficientX1*dt2 +
			coefficientX2*dt3/3.0 + 0.25*coefficientX3*dt4;
		runningHeatY += sumHeatY*interval + 0.5*coefficientY1*dt2 +
			coefficientY2*dt3/3.0;
		sumHeatX += coefficientX1*interval + coefficientX2*dt2 +
			coefficientX3*dt3;
		sumHeatY += coefficientY1*interval + coefficientY2*dt2;
	#if THREE_D
		const double coefficientZ1 = rate*shearHeatCxVyVz;
		const double coefficientZ2 = 0.5*rate2*shearHeatVy2Vz;
		runningHeatZ += sumHeatZ*interval + 0.5*coefficientZ1*dt2 +
			coefficientZ2*dt3/3.0;
		sumHeatZ += coefficientZ1*interval + coefficientZ2*dt2;
	#endif
		const double oldCx2Vy = shearHeatCx2Vy;
		const double oldCxVy2 = shearHeatCxVy2;
		const double oldCxVyVz = shearHeatCxVyVz;
		shearHeatCx2Vy = oldCx2Vy + 2.0*rate*oldCxVy2*interval +
			rate2*shearHeatVy3*dt2;
		shearHeatCxVy2 = oldCxVy2 + rate*shearHeatVy3*interval;
		shearHeatCxVyVz = oldCxVyVz +
			rate*shearHeatVy2Vz*interval;
	}
	else{
		runningHeatX += sumHeatX*interval;
		runningHeatY += sumHeatY*interval;
	#if THREE_D
		runningHeatZ += sumHeatZ*interval;
	#endif
	}
	#endif
	lastCollTime = t;
}

static void updateParticleTransport(const particle* p, double oldVx, double oldVy,
	double oldVz){
	const double oldTransportVx = transportVelocityX(p, oldVx);
	const double newTransportVx = transportVelocityX(p, p->vx);
	#if STRESS
	sumKineticX += p->m*(newTransportVx*newTransportVx -
		oldTransportVx*oldTransportVx);
	sumKineticY += p->m*(p->vy*p->vy - oldVy*oldVy);
	sumKineticXY += p->m*(newTransportVx*p->vy -
		oldTransportVx*oldVy);
	sumKineticYX += p->m*(p->vy*newTransportVx -
		oldVy*oldTransportVx);
	#if THREE_D
	sumKineticZ += p->m*(p->vz*p->vz - oldVz*oldVz);
	sumKineticXZ += p->m*(newTransportVx*p->vz -
		oldTransportVx*oldVz);
	sumKineticZX += p->m*(p->vz*newTransportVx -
		oldVz*oldTransportVx);
	sumKineticYZ += p->m*(p->vy*p->vz - oldVy*oldVz);
	sumKineticZY += p->m*(p->vz*p->vy - oldVz*oldVy);
	#endif
	#endif
	#if HEAT
	double oldEnergy = 0.5*p->m*(oldTransportVx*oldTransportVx +
		oldVy*oldVy
		#if THREE_D
		+ oldVz*oldVz
		#endif
	);
	double newEnergy = 0.5*p->m*(newTransportVx*newTransportVx +
		p->vy*p->vy
		#if THREE_D
		+ p->vz*p->vz
		#endif
	);
	sumHeatX += newEnergy*newTransportVx - oldEnergy*oldTransportVx;
	sumHeatY += newEnergy*p->vy - oldEnergy*oldVy;
	#if THREE_D
	const double newVz = p->vz;
	sumHeatZ += newEnergy*newVz - oldEnergy*oldVz;
	#else
	const double newVz = 0.0;
	#endif
	if (shearDynamicsActive()){
		const double mass = p->m;
		shearHeatCx2Vy += mass*(newTransportVx*newTransportVx*p->vy -
			oldTransportVx*oldTransportVx*oldVy);
		shearHeatVy3 += mass*(p->vy*p->vy*p->vy -
			oldVy*oldVy*oldVy);
		shearHeatVyVz2 += mass*(p->vy*newVz*newVz -
			oldVy*oldVz*oldVz);
		shearHeatCxVy2 += mass*(newTransportVx*p->vy*p->vy -
			oldTransportVx*oldVy*oldVy);
		shearHeatVy2Vz += mass*(p->vy*p->vy*newVz -
			oldVy*oldVy*oldVz);
		shearHeatCxVyVz += mass*(newTransportVx*p->vy*newVz -
			oldTransportVx*oldVy*oldVz);
	}
	#else
	(void)oldVz;
	#endif
}

#if HEAT
static double recordPairHeat(const particle* pi, const particle* pj,
	double oldVix, double oldViy, double oldViz,
	double oldVjx, double oldVjy, double oldVjz,
	double separationX, double separationY, double separationZ,
	int recordNonconservativeSource){
	const double oldTransportVix = transportVelocityX(pi, oldVix);
	const double oldTransportVjx = transportVelocityX(pj, oldVjx);
	const double newTransportVix = transportVelocityX(pi, pi->vx);
	const double newTransportVjx = transportVelocityX(pj, pj->vx);
	double oldEi = 0.5*pi->m*(oldTransportVix*oldTransportVix + oldViy*oldViy
		#if THREE_D
		+ oldViz*oldViz
		#endif
	);
	double oldEj = 0.5*pj->m*(oldTransportVjx*oldTransportVjx + oldVjy*oldVjy
		#if THREE_D
		+ oldVjz*oldVjz
		#endif
	);
	double newEi = 0.5*pi->m*(newTransportVix*newTransportVix + pi->vy*pi->vy
		#if THREE_D
		+ pi->vz*pi->vz
		#endif
	);
	double newEj = 0.5*pj->m*(newTransportVjx*newTransportVjx + pj->vy*pj->vy
		#if THREE_D
		+ pj->vz*pj->vz
		#endif
	);
	double transfer = 0.5*((newEi - oldEi) - (newEj - oldEj));
	collHeatX += transfer*separationX;
	collHeatY += transfer*separationY;
	#if THREE_D
	collHeatZ += transfer*separationZ;
	if (recordNonconservativeSource){
		/* A Lees image remapping changes stored vx while leaving peculiar
		 * velocity continuous.  Conversely, shear work can change peculiar
		 * kinetic energy in an otherwise elastic collision.  The explicit
		 * nonconservative source diagnostic should therefore retain the raw
		 * collision-energy balance, while heat transfer uses peculiar energy. */
		const double rawOldEi = 0.5*pi->m*(oldVix*oldVix + oldViy*oldViy +
			oldViz*oldViz);
		const double rawOldEj = 0.5*pj->m*(oldVjx*oldVjx + oldVjy*oldVjy +
			oldVjz*oldVjz);
		const double rawNewEi = 0.5*pi->m*(pi->vx*pi->vx + pi->vy*pi->vy +
			pi->vz*pi->vz);
		const double rawNewEj = 0.5*pj->m*(pj->vx*pj->vx + pj->vy*pj->vy +
			pj->vz*pj->vz);
		collisionEnergySource += (rawNewEi - rawOldEi) +
			(rawNewEj - rawOldEj);
	}
	#else
	(void)oldViz;
	(void)oldVjz;
	(void)separationZ;
	(void)recordNonconservativeSource;
	#endif
	return transfer;
}
#endif

static void rebuildInstantaneousTransportSums(){
	#if STRESS
	sumKineticX = sumKineticY = sumKineticXY = sumKineticYX = 0;
	#if THREE_D
	sumKineticZ = sumKineticXZ = sumKineticZX = sumKineticYZ = sumKineticZY = 0;
	#endif
	#endif
	#if HEAT
	sumHeatX = sumHeatY = 0;
	#if THREE_D
	sumHeatZ = 0;
	#endif
	shearHeatCx2Vy = shearHeatVy3 = shearHeatVyVz2 = 0;
	shearHeatCxVy2 = shearHeatVy2Vz = shearHeatCxVyVz = 0;
	#endif
	for (int i = 0; i < N; i++){
		particle* p = particles + i;
		#if (STRESS || HEAT)
		const double transportVx = transportVelocityX(p, p->vx);
		#endif
		#if STRESS
		sumKineticX += p->m*transportVx*transportVx;
		sumKineticY += p->m*p->vy*p->vy;
		sumKineticXY += p->m*transportVx*p->vy;
		sumKineticYX += p->m*p->vy*transportVx;
		#if THREE_D
		sumKineticZ += p->m*p->vz*p->vz;
		sumKineticXZ += p->m*transportVx*p->vz;
		sumKineticZX += p->m*p->vz*transportVx;
		sumKineticYZ += p->m*p->vy*p->vz;
		sumKineticZY += p->m*p->vz*p->vy;
		#endif
		#endif
		#if HEAT
		double energy = 0.5*p->m*(transportVx*transportVx + p->vy*p->vy
			#if THREE_D
			+ p->vz*p->vz
			#endif
		);
		sumHeatX += energy*transportVx;
		sumHeatY += energy*p->vy;
		#if THREE_D
		const double transportVz = p->vz;
		sumHeatZ += energy*transportVz;
		#else
		const double transportVz = 0.0;
		#endif
		if (shearDynamicsActive()){
			const double mass = p->m;
			shearHeatCx2Vy += mass*transportVx*transportVx*p->vy;
			shearHeatVy3 += mass*p->vy*p->vy*p->vy;
			shearHeatVyVz2 += mass*p->vy*transportVz*transportVz;
			shearHeatCxVy2 += mass*transportVx*p->vy*p->vy;
			shearHeatVy2Vz += mass*p->vy*p->vy*transportVz;
			shearHeatCxVyVz += mass*transportVx*p->vy*transportVz;
		}
		#endif
	}
	lastCollTime = t;
}

void initializeTransportAccumulators(){
	#if STRESS
	runningKineticX = runningKineticY = runningKineticXY = runningKineticYX = 0;
	#if THREE_D
	runningKineticZ = runningKineticXZ = runningKineticZX = runningKineticYZ = runningKineticZY = 0;
	#endif
	#endif
	#if HEAT
	runningHeatX = runningHeatY = collHeatX = collHeatY = 0;
	#if THREE_D
	runningHeatZ = collHeatZ = collisionEnergySource = 0;
	#endif
	#endif
	rebuildInstantaneousTransportSums();
}

/* Close the old model's transport interval and establish the instantaneous
 * kinetic/heat moments for a newly selected flight/shear model.  Historical
 * integrals and collision transfers intentionally remain untouched. */
void rebaseTransportAccumulatorsAfterModelChange(void){
	advanceTransportAccumulators();
	rebuildInstantaneousTransportSums();
}
#endif

void beginTransportVelocityBatch(void){
	#if (STRESS || HEAT)
	advanceTransportAccumulators();
	#endif
}

void recordTransportVelocityChange(const particle* p,
	double oldVx, double oldVy, double oldVz){
	#if (STRESS || HEAT)
	updateParticleTransport(p, oldVx, oldVy, oldVz);
	#else
	(void)p;
	(void)oldVx;
	(void)oldVy;
	(void)oldVz;
	#endif
}

typedef struct velocitySnapshot velocitySnapshot;
struct velocitySnapshot{
	double x, y, z;
};

static velocitySnapshot captureVelocity(const particle* p){
	return (velocitySnapshot){
		.x = p->vx,
		.y = p->vy,
		#if THREE_D
		.z = p->vz,
		#else
		.z = 0,
		#endif
	};
}

static void recordCapturedVelocity(const particle* p, velocitySnapshot old){
	recordTransportVelocityChange(p, old.x, old.y, old.z);
}

typedef struct planarWallContact planarWallContact;
struct planarWallContact{
	double coordinate;
	double outwardNormal;
};

static void considerPlanarWallContact(double particleCoordinate,
	double coordinate, double outwardNormal, planarWallContact* closest,
	double* closestDistance){
	const double distance = fabs(particleCoordinate - coordinate);
	if (distance < *closestDistance){
		*closestDistance = distance;
		closest->coordinate = coordinate;
		closest->outwardNormal = outwardNormal;
	}
}

/* Recover which face produced an axis-only wall event.  Legacy events store
 * only the axis, not the side.  Choosing the closest exact contact surface is
 * robust for outer walls and for either face of the optional center membrane. */
static planarWallContact closestPlanarWallContact(const particle* p, int axis){
	planarWallContact closest = {0};
	double closestDistance = DBL_MAX;
	if (axis == 0){
		considerPlanarWallContact(p->x, p->rad, -1, &closest,
			&closestDistance);
		considerPlanarWallContact(p->x, Lx - p->rad, 1, &closest,
			&closestDistance);
		if (addMidWall && p->type == 1){
			considerPlanarWallContact(p->x, halfLx - p->rad, 1,
				&closest, &closestDistance);
			considerPlanarWallContact(p->x, halfLx + p->rad, -1,
				&closest, &closestDistance);
		}
	}
	else if (axis == 1){
		considerPlanarWallContact(p->y, p->rad, -1, &closest,
			&closestDistance);
		considerPlanarWallContact(p->y, Ly - p->rad, 1, &closest,
			&closestDistance);
	}
	#if THREE_D
	else if (axis == 3){
		considerPlanarWallContact(p->z, p->rad, -1, &closest,
			&closestDistance);
		considerPlanarWallContact(p->z, Lz - p->rad, 1, &closest,
			&closestDistance);
	}
	#endif
	return closest;
}

double confiningWallNormalVelocityAfter(double normalVelocity,
	double inwardSpeed, double restitution){
	/* The boundary velocity along its outward normal is -inwardSpeed. */
	return -restitution*(normalVelocity + inwardSpeed) - inwardSpeed;
}

static bool projectParticleToCircularWall(particle* p,
	double* nx, double* ny, bool* outsideRepair
	#if THREE_D
	, double* nz
	#endif
){
	const double accessibleRadius = halfLx - p->rad;
	if (!(accessibleRadius >= 0) || !isfinite(accessibleRadius))
		return false;
	double rx = p->x - halfLx;
	double ry = p->y - halfLy;
	#if THREE_D
	double rz = p->z - halfLz;
	const double distance = hypot(hypot(rx, ry), rz);
	#else
	const double distance = hypot(rx, ry);
	#endif
	if (!isfinite(distance))
		return false;
	const double distanceTolerance = 64.0*DBL_EPSILON*
		(distance + accessibleRadius + 1.0);
	if (outsideRepair != NULL)
		*outsideRepair = distance > accessibleRadius + distanceTolerance;
	if (!(distance > DBL_MIN)){
		/* A genuine wall event cannot occur at the center unless the accessible
		 * radius vanished.  Do not teleport a malformed positive-radius state
		 * from the center to an arbitrary point on the wall. */
		if (accessibleRadius > distanceTolerance)
			return false;
		rx = 1;
		ry = 0;
		#if THREE_D
		rz = 0;
		#endif
	}
	else{
		rx /= distance;
		ry /= distance;
		#if THREE_D
		rz /= distance;
		#endif
	}
	*nx = rx;
	*ny = ry;
	#if THREE_D
	*nz = rz;
	#endif
	p->x = halfLx + accessibleRadius*rx;
	p->y = halfLy + accessibleRadius*ry;
	#if THREE_D
	p->z = halfLz + accessibleRadius*rz;
	#endif
	return true;
}

/* Corrective projection can move a malformed outside state across one or
 * more cell faces.  The linked-list membership must follow the projected
 * coordinates before new crossing and collision events are constructed. */
static void repairParticleCellAfterCircularProjection(int i){
	particle* p = particles + i;
	const int projectedX = coordToCell(p->x, 1);
	const int projectedY = coordToCell(p->y, 0);
	bool changed = projectedX != p->cell[0] || projectedY != p->cell[1];
	#if THREE_D
	const int projectedZ = coordToCell(p->z, 2);
	changed = changed || projectedZ != p->cell[2];
	#endif
	if (changed){
		removeFromCell(i);
		addToCell(i);
	}
}


void doTheWallGrow(){
	int i = nextEvent->i;
	int xy = nextEvent->j;
	ncol++;
	particle* pi = particles + i;
	freeFlyGrow(pi);

	double vrParticle = pi->vr;


	if (xy == 0 || xy == 1
		#if THREE_D
		|| xy == 3
		#endif
	){
		const planarWallContact contact = closestPlanarWallContact(pi, xy);
		double* velocity = xy == 0 ? &pi->vx : &pi->vy;
		#if THREE_D
		if (xy == 3)
			velocity = &pi->vz;
		#endif
		if (xy == 0)
			pi->x = contact.coordinate;
		else if (xy == 1)
			pi->y = contact.coordinate;
		#if THREE_D
		else
			pi->z = contact.coordinate;
		#endif
		const double normalVelocity = contact.outwardNormal*(*velocity);
		const double relativeNormalVelocity = normalVelocity + vrParticle;
		const double velocityTolerance = 64*DBL_EPSILON*
			(fabs(normalVelocity) + fabs(vrParticle) + 1);
		if (relativeNormalVelocity > velocityTolerance){
			const double normalAfter = confiningWallNormalVelocityAfter(
				normalVelocity, vrParticle, resW);
			*velocity += (normalAfter - normalVelocity)*contact.outwardNormal;
		}
	}
	else{
		double nx, ny;
		#if THREE_D
		double nz;
		if (projectParticleToCircularWall(pi, &nx, &ny, NULL, &nz)){
		#else
		if (projectParticleToCircularWall(pi, &nx, &ny, NULL)){
		#endif
			repairParticleCellAfterCircularProjection(i);
			double normalVelocity = pi->vx*nx + pi->vy*ny;
			#if THREE_D
			normalVelocity += pi->vz*nz;
			#endif
			const double relativeNormalVelocity = normalVelocity + vrParticle;
			const double velocityTolerance = 64*DBL_EPSILON*
				(fabs(normalVelocity) + fabs(vrParticle) + 1);
			if (relativeNormalVelocity > velocityTolerance){
				const double normalAfter = confiningWallNormalVelocityAfter(
					normalVelocity, vrParticle, resW);
				const double velocityChange = normalAfter - normalVelocity;
				pi->vx += velocityChange*nx;
				pi->vy += velocityChange*ny;
				#if THREE_D
				pi->vz += velocityChange*nz;
				#endif
			}
		}
	}

	pi->coll++;
	removeEventFromQueue(EVENT_NODE(i));
	crossingEventGrow(i);
	collisionEventGrow(i);
}

static void thermalizePlanarWall(particle* p, int normalAxis,
	double outwardNormal){
	const double normalSpeed = sqrt(-2*T*log(genrand_real3())/p->m);
	const double u1 = genrand_real3();
	const double u2 = genrand_real3();
	const double gaussianScale = sqrt(T/p->m)*sqrt(-2*log(u1));
	const double tangent1 = gaussianScale*cos(2*M_PI*u2);
	#if THREE_D
	const double tangent2 = gaussianScale*sin(2*M_PI*u2);
	#endif
	if (normalAxis == 0){
		p->vx = -outwardNormal*normalSpeed;
		p->vy = tangent1;
		#if THREE_D
		p->vz = tangent2;
		#endif
	}
	else if (normalAxis == 1){
		p->vy = -outwardNormal*normalSpeed;
		p->vx = tangent1;
		#if THREE_D
		p->vz = tangent2;
		#endif
	}
	#if THREE_D
	else{
		p->vz = -outwardNormal*normalSpeed;
		p->vx = tangent1;
		p->vy = tangent2;
	}
	#endif
}

void doTheWallNormal(){
	int i = nextEvent->i;
	int xy = nextEvent->j;
	ncol++;
	particle* pi = particles + i;
	freeFlyNormal(pi);
	double wallOldVx = pi->vx;
	double wallOldVy = pi->vy;
	#if THREE_D
	double wallOldVz = pi->vz;
	#endif
	#if (STRESS || HEAT)
	advanceTransportAccumulators();
	#endif

	if (xy == 0 || xy == 1
		#if THREE_D
		|| xy == 3
		#endif
	){
		const planarWallContact contact = closestPlanarWallContact(pi, xy);
		if (xy == 0)
			pi->x = contact.coordinate;
		else if (xy == 1)
			pi->y = contact.coordinate;
		#if THREE_D
		else
			pi->z = contact.coordinate;
		#endif
		double oldNormalComponent = xy == 0 ? wallOldVx : wallOldVy;
		#if THREE_D
		if (xy == 3)
			oldNormalComponent = wallOldVz;
		#endif
		const double normalVelocity = contact.outwardNormal*oldNormalComponent;
		const double velocityTolerance = 64*DBL_EPSILON*
			(fabs(normalVelocity) + 1);
		const bool fieldPressesIntoWall = xy == 1 && addField &&
			field*contact.outwardNormal > 0;
		const bool wallImpact = normalVelocity > velocityTolerance ||
			(fabs(normalVelocity) <= velocityTolerance && fieldPressesIntoWall);
		if (wallImpact){
			if (thermoWall){
				thermalizePlanarWall(pi, xy == 3 ? 2 : xy,
					contact.outwardNormal);
			}
			else{
				double normalAfter = -resW*fmax(normalVelocity, 0);
				/* Prevent the finite-time wall-collapse sequence of an inelastic
				 * particle settling under a constant field.  Below this tiny cutoff
				 * impacts are treated elastically at the cutoff velocity. */
				if (fieldPressesIntoWall){
					const double collapseCutoff = 1e-3*
						sqrt(fabs(field)*fmax(pi->rad, 1.0));
					if (-normalAfter < collapseCutoff)
						normalAfter = -collapseCutoff;
				}
				if (xy == 0){
					pi->vx = contact.outwardNormal*normalAfter;
					pi->vy = resT*wallOldVy;
					#if THREE_D
					pi->vz = resT*wallOldVz;
					#endif
				}
				else if (xy == 1){
					pi->vy = contact.outwardNormal*normalAfter;
					pi->vx = resT*wallOldVx;
					#if THREE_D
					pi->vz = resT*wallOldVz;
					#endif
				}
				#if THREE_D
				else{
					pi->vz = contact.outwardNormal*normalAfter;
					pi->vx = resT*wallOldVx;
					pi->vy = resT*wallOldVy;
				}
				#endif
			}
			if (xy == 1 && addShearWall){
				/* Couette wall: prescribe x; z keeps the restitution/thermal rule. */
				pi->vx = contact.outwardNormal > 0 ? shearWallSpeed :
					-shearWallSpeed;
				if (addDelta)
					pi->vy += delta*sign(pi->vy);
			}
		}
		double newNormalComponent = xy == 0 ? pi->vx : pi->vy;
		#if THREE_D
		if (xy == 3)
			newNormalComponent = pi->vz;
		#endif
		const double normalImpulse = fabs(pi->m*
			(newNormalComponent - oldNormalComponent));
		dp += normalImpulse;
		cumulativeWallImpulse += normalImpulse;
		if (xy == 0 && addMidWall){
			const double surfaceTolerance = 64*DBL_EPSILON*(Lx + 1);
			if (fabs(contact.coordinate - pi->rad) <= surfaceTolerance)
				dpLeft += normalImpulse;
			else if (fabs(contact.coordinate - (Lx - pi->rad)) <=
				surfaceTolerance)
				dpRight += normalImpulse;
			else
				dpMid += normalImpulse;
		}
	}
	else{
		double nx, ny;
		bool outsideRepair = false;
		#if THREE_D
		double nz;
		const bool validSurface = projectParticleToCircularWall(
			pi, &nx, &ny, &outsideRepair, &nz);
		#else
		const bool validSurface = projectParticleToCircularWall(
			pi, &nx, &ny, &outsideRepair);
		#endif
		if (validSurface){
			repairParticleCellAfterCircularProjection(i);
			#if THREE_D
			double normalVelocity = wallOldVx*nx + wallOldVy*ny + wallOldVz*nz;
			double tangentSpeedSquared =
				wallOldVx*wallOldVx + wallOldVy*wallOldVy + wallOldVz*wallOldVz -
				normalVelocity*normalVelocity;
			#else
			double normalVelocity = wallOldVx*nx + wallOldVy*ny;
			double tangentSpeedSquared =
				wallOldVx*wallOldVx + wallOldVy*wallOldVy -
				normalVelocity*normalVelocity;
			#endif
			if (tangentSpeedSquared < 0)
				tangentSpeedSquared = 0;
			const double velocityTolerance = 64*DBL_EPSILON*
				(fabs(normalVelocity) + sqrt(tangentSpeedSquared) + 1);
			/* An inward velocity means this was a roundoff re-entry event from a
			 * center a few ulps outside. Projection repairs it without reflecting
			 * the particle back out of the sphere. */
			const bool wallImpact = normalVelocity >= -velocityTolerance &&
				!(outsideRepair && normalVelocity < 0);
			if (wallImpact && thermoWall == 0){
				const double tangentScale = resT;
				double inwardSpeed = resW*fmax(normalVelocity, 0);
				/* Curved walls with zero/very small normal restitution otherwise
				 * undergo an infinite collision sequence. This standard collapse
				 * cutoff preserves the requested law above 0.1% of tangential speed. */
				double collapseCutoff = 1e-3*fabs(tangentScale)*
					sqrt(tangentSpeedSquared);
				if (inwardSpeed < collapseCutoff)
					inwardSpeed = collapseCutoff;
				pi->vx = tangentScale*(wallOldVx - normalVelocity*nx) -
					inwardSpeed*nx;
				pi->vy = tangentScale*(wallOldVy - normalVelocity*ny) -
					inwardSpeed*ny;
				#if THREE_D
				pi->vz = tangentScale*(wallOldVz - normalVelocity*nz) -
					inwardSpeed*nz;
				#endif
				/* A constant field tangent to a curved surface can produce an
				 * immediate t^4 re-contact even when its normal component is zero.
				 * Test the complete post-impact trajectory, not only a.n, before
				 * applying the small collapse regularization. */
				if (addField && circularWallTimeNormal(pi) == 0){
					const double fieldCutoff = 1e-3*sqrt(
						fabs(field)*fmax(pi->rad, 1.0));
					if (fieldCutoff > inwardSpeed){
						const double additionalInwardSpeed =
							fieldCutoff - inwardSpeed;
						pi->vx -= additionalInwardSpeed*nx;
						pi->vy -= additionalInwardSpeed*ny;
						#if THREE_D
						pi->vz -= additionalInwardSpeed*nz;
						#endif
					}
				}
			}
			else if (wallImpact){
			/* Flux-weighted Maxwellian normal speed and Maxwellian tangential speed. */
			double u = genrand_real3();
			double normalSpeed = sqrt(-2*T*log(u)/pi->m);
			double u1 = genrand_real3();
			double u2 = genrand_real3();
			double gaussian1 = sqrt(-2*log(u1))*cos(2*M_PI*u2);
			#if THREE_D
			double gaussian2 = sqrt(-2*log(u1))*sin(2*M_PI*u2);
			double ax = fabs(nz) < 0.9 ? 0 : 1;
			double ay = 0;
			double az = fabs(nz) < 0.9 ? 1 : 0;
			double tx = ay*nz - az*ny;
			double ty = az*nx - ax*nz;
			double tz = ax*ny - ay*nx;
			double tangentNorm = sqrt(tx*tx + ty*ty + tz*tz);
			tx /= tangentNorm;
			ty /= tangentNorm;
			tz /= tangentNorm;
			double sx = ny*tz - nz*ty;
			double sy = nz*tx - nx*tz;
			double sz = nx*ty - ny*tx;
			double tangentScale = sqrt(T/pi->m);
			pi->vx = -normalSpeed*nx + tangentScale*(gaussian1*tx + gaussian2*sx);
			pi->vy = -normalSpeed*ny + tangentScale*(gaussian1*ty + gaussian2*sy);
			pi->vz = -normalSpeed*nz + tangentScale*(gaussian1*tz + gaussian2*sz);
			#else
			double tangentSpeed = sqrt(T/pi->m)*gaussian1;
			pi->vx = -normalSpeed*nx - tangentSpeed*ny;
			pi->vy = -normalSpeed*ny + tangentSpeed*nx;
			#endif
			}
			double normalImpulse = pi->m*((pi->vx - wallOldVx)*nx + (pi->vy - wallOldVy)*ny
			#if THREE_D
			+ (pi->vz - wallOldVz)*nz
			#endif
			);
			const double impulseMagnitude = fabs(normalImpulse);
			dp += impulseMagnitude;
			cumulativeWallImpulse += impulseMagnitude;
		}
	}
	#if (STRESS || HEAT)
	updateParticleTransport(pi, wallOldVx, wallOldVy,
		#if THREE_D
		wallOldVz
		#else
		0
		#endif
	);
	#endif

	pi->coll++;
	removeEventFromQueue(EVENT_NODE(i));
	crossingEventNormal(i);
	collisionEventNormal(i);
}
/* -----------------------------------
	Performs the collision event
----------------------------------- */
void doTheCollisionGrow(){

	double delta = 0.1*0;
	double res = 0.75*0 + 1;

	int i = nextEvent->i;
	int j = nextEvent->j;

	particle* pi = particles + i;
	particle* pj = particles + j;


	freeFlyGrow(pi);
	if (nextEvent->collActual != pj->coll){ //the collision trick of the article
		collisionEventGrow(i);
        return;
    }

    ncol++;
	freeFlyGrow(pj);

	pi->coll++;
	pj->coll++;

	double dx = pj->x - pi->x;
	double dy = pj->y - pi->y;
	

	double dvx = pj->vx - pi->vx;
	double dvy = pj->vy - pi->vy;
	double dvr = pi->vr + pj->vr;

	#if THREE_D
	double dz = pj->z - pi->z;
	double dvz = pj->vz - pi->vz;
	PBC(&dx, &dy, &dz);
	dist = sqrt(dx*dx + dy*dy + dz*dz);
	dzr = dz/dist;
	#else
	PBC(&dx, &dy);
	dist = sqrt(dx*dx + dy*dy);
	#endif
	dxr = dx/dist;
	dyr = dy/dist;
	

	double invMass = 1/(pi->m + pj->m);
	#if THREE_D
	double collKernel = invMass*(1 + res)*(dxr*dvx + dyr*dvy + dzr*dvz - dvr);
	pi->vz += collKernel*pj->m*dzr - 2*pj->m*invMass*delta*dzr;
	pj->vz -= collKernel*pi->m*dzr - 2*pi->m*invMass*delta*dzr;
	#else
	double collKernel = invMass*(1 + res)*(dxr*dvx + dyr*dvy - dvr);
	#endif

	pi->vx += collKernel*pj->m*dxr - 2*pj->m*invMass*delta*dxr;
	pi->vy += collKernel*pj->m*dyr - 2*pj->m*invMass*delta*dyr;
	pj->vx -= collKernel*pi->m*dxr - 2*pi->m*invMass*delta*dxr;
	pj->vy -= collKernel*pi->m*dyr - 2*pi->m*invMass*delta*dyr;
	

	

    pi->lastColl = t;
    pj->lastColl = t;

	//recomputes crossing and collision of pi and pj
	removeEventFromQueue(EVENT_NODE(i));
	removeEventFromQueue(EVENT_NODE(j));
	crossingEventGrow(i);
	crossingEventGrow(j);

	removeEventFromQueue(EVENT_NODE(N + j));
	collisionEventGrow(i);
	collisionEventGrow(j);
}

void doTheCollisionNormal(){
	int i = nextEvent->i;
	int j = nextEvent->j;

	particle* pi = particles + i;
	particle* pj = particles + j;


	freeFlyNormal(pi);


	if (nextEvent->collActual != pj->coll){ //the collision trick of the article
		collisionEventNormal(i);
        return;
    }

    ncol++;
	freeFlyNormal(pj);

	pi->coll++;
	pj->coll++;

	#if (STRESS || HEAT)
	advanceTransportAccumulators();
	double oldvix = pi->vx;
	double oldviy = pi->vy;
	double oldvjx = pj->vx;
	double oldvjy = pj->vy;
	#if THREE_D
	double oldviz = pi->vz;
	double oldvjz = pj->vz;
	#else
	double oldviz = 0;
	double oldvjz = 0;
	#endif
	#endif
	#if EDMD_VERIFY_COLLISION_ACCOUNTING && !TANGENTIAL
	const double verifyOldVix = pi->vx;
	const double verifyOldViy = pi->vy;
	const double verifyOldVjx = pj->vx;
	const double verifyOldVjy = pj->vy;
	#if THREE_D
	const double verifyOldViz = pi->vz;
	const double verifyOldVjz = pj->vz;
	const double verifyOldVirialZ = collTermZ;
	const double verifyOldVirialXZ = collTermXZ;
	const double verifyOldVirialZX = collTermZX;
	const double verifyOldVirialYZ = collTermYZ;
	const double verifyOldVirialZY = collTermZY;
	#else
	const double verifyOldViz = 0;
	const double verifyOldVjz = 0;
	#endif
	const double verifyOldVirialX = collTermX;
	const double verifyOldVirialY = collTermY;
	const double verifyOldVirialXY = collTermXY;
	const double verifyOldVirialYX = collTermYX;
	#endif

	double dx = pj->x - pi->x;
	double dy = pj->y - pi->y;
	
	double dvx = pj->vx - pi->vx;
	double dvy = pj->vy - pi->vy;

	if (addShear){
		correctDistances(pi, pj, pj->x, pj->vx, &dx, &dvx);
	}
	double invMass = 1/(pi->m + pj->m);
	double effectiveRestitution = res;
	if (addExpo){
		double relativeSpeedSquared = dvx*dvx + dvy*dvy;
		#if THREE_D
		relativeSpeedSquared += (pj->vz - pi->vz)*(pj->vz - pi->vz);
		#endif
		effectiveRestitution = resCoeff(sqrt(relativeSpeedSquared));
	}

	#if THREE_D
	double dz = pj->z - pi->z;
	double dvz = pj->vz - pi->vz;
	PBC(&dx, &dy, &dz);
	double collKernel = invMass*(1 + effectiveRestitution)*(dx*dvx + dy*dvy + dz*dvz);
	#else
	PBC(&dx, &dy);
	double collKernel = invMass*(1 + effectiveRestitution)*(dx*dvx + dy*dvy);
	#endif
	
	double distSquared = DIST2(pi->rad, pj->rad);
	double funkyFactor = collKernel/distSquared;

	#if !TANGENTIAL
	if (addDelta || addDoubleDelta || addEvolvingDelta || addDeltaTangent || addDoubleDeltaTangent){
		/* Keep the configured injection amplitude declarative.  The double and
		 * evolving laws choose a per-collision amplitude; writing that value back
		 * to the global `delta` made an unrelated live GUI edit capture internal
		 * collision state as though it were user configuration. */
		double collisionDelta = delta;
		//double tau = -ts*log(1 - drand(0, 1));
		if (addDoubleDelta){
			if ((t - pi->lastColl > ts) && (t - pj->lastColl > ts)){
				collisionDelta = deltam;
			}
			else{
				collisionDelta = deltaM;
			}
		}
		if (addDoubleDeltaTangent){
			if ((t - pi->lastColl < ts) && (t - pj->lastColl < ts)){
				collisionDelta = deltam;
			}
			else{
				collisionDelta = deltaM;
			}
		}
		else if (addEvolvingDelta){
			double dti = t - pi->lastColl;
			double dtj = t - pj->lastColl;
			collisionDelta = deltaMax*pow((2 - (exp(-dtj/tau) + exp(-dti/tau)))/2, 3);
		}

		/* This event was predicted on the hard contact surface, whose distance
		 * is already known from the two radii.  Recomputing hypot(dx,dy[,dz])
		 * here adds a square root to every Delta collision and only re-measures
		 * accumulated flight roundoff. */
		const double inverseContactDistance = 1.0/DIST(pi->rad, pj->rad);
		double nx = dx*inverseContactDistance;
		double ny = dy*inverseContactDistance;
		#if THREE_D
		double nz = dz*inverseContactDistance;
		#endif
		double deltaDirectionX = nx;
		double deltaDirectionY = ny;
		#if THREE_D
		double deltaDirectionZ = nz;
		#endif
		double deltaFactor = -2*collisionDelta*invMass;
		if (addDeltaTangent || addDoubleDeltaTangent){
			/* Fixed-axis chirality: t = u_z x n.  It is deliberately not
			   normalized, so polar contacts remain smooth and nonsingular. */
			deltaDirectionX = -ny;
			deltaDirectionY = nx;
			#if THREE_D
			deltaDirectionZ = 0;
			#endif
			deltaFactor = 2*collisionDelta*invMass;
		}

		double dviX = funkyFactor*pj->m*dx + pj->m*deltaFactor*deltaDirectionX;
		double dviY = funkyFactor*pj->m*dy + pj->m*deltaFactor*deltaDirectionY;
		double dvjX = -funkyFactor*pi->m*dx - pi->m*deltaFactor*deltaDirectionX;
		double dvjY = -funkyFactor*pi->m*dy - pi->m*deltaFactor*deltaDirectionY;
		#if THREE_D
		double dviZ = funkyFactor*pj->m*dz + pj->m*deltaFactor*deltaDirectionZ;
		double dvjZ = -funkyFactor*pi->m*dz - pi->m*deltaFactor*deltaDirectionZ;
		#else
		double dviZ = 0;
		#endif
		pi->vx += dviX;
		pi->vy += dviY;
		pj->vx += dvjX;
		pj->vy += dvjY;
		#if THREE_D
		pi->vz += dviZ;
		pj->vz += dvjZ;
		#endif

		double impulseX = pi->m*dviX;
		double impulseY = pi->m*dviY;
		double impulseZ = pi->m*dviZ;
		accumulateCollisionVirial(impulseX, impulseY, impulseZ, dx, dy,
			#if THREE_D
			dz
			#else
			0
			#endif
		);
		#if EDMD_VERIFY_COLLISION_ACCOUNTING
		if (addDeltaTangent || addDoubleDeltaTangent){
			const double tx = -ny;
			const double ty = nx;
			const double tz = 0;
			const double oldGx = verifyOldVjx - verifyOldVix;
			const double oldGy = verifyOldVjy - verifyOldViy;
			const double oldGz = verifyOldVjz - verifyOldViz;
			const double normalSpeed = oldGx*nx + oldGy*ny
				#if THREE_D
				+ oldGz*nz
				#endif
			;
			const double expectedGx = oldGx - (1 + effectiveRestitution)*normalSpeed*nx - 2*collisionDelta*tx;
			const double expectedGy = oldGy - (1 + effectiveRestitution)*normalSpeed*ny - 2*collisionDelta*ty;
			const double expectedGz = oldGz
				#if THREE_D
				- (1 + effectiveRestitution)*normalSpeed*nz
				#endif
				- 2*collisionDelta*tz;
			const double newGx = pj->vx - pi->vx;
			const double newGy = pj->vy - pi->vy;
			const double newGz =
				#if THREE_D
				pj->vz - pi->vz
				#else
				0
				#endif
			;
			const double oldEnergy = 0.5*pi->m*(verifyOldVix*verifyOldVix + verifyOldViy*verifyOldViy
				#if THREE_D
				+ verifyOldViz*verifyOldViz
				#endif
			) + 0.5*pj->m*(verifyOldVjx*verifyOldVjx + verifyOldVjy*verifyOldVjy
				#if THREE_D
				+ verifyOldVjz*verifyOldVjz
				#endif
			);
			const double newEnergy = 0.5*pi->m*(pi->vx*pi->vx + pi->vy*pi->vy
				#if THREE_D
				+ pi->vz*pi->vz
				#endif
			) + 0.5*pj->m*(pj->vx*pj->vx + pj->vy*pj->vy
				#if THREE_D
				+ pj->vz*pj->vz
				#endif
			);
			const double reducedMass = pi->m*pj->m/(pi->m + pj->m);
			const double tangentSpeed = oldGx*tx + oldGy*ty + oldGz*tz;
			const double tangentNormSquared = tx*tx + ty*ty + tz*tz;
			const double expectedEnergyChange = 0.5*reducedMass*
				(effectiveRestitution*effectiveRestitution - 1)*normalSpeed*normalSpeed -
				2*reducedMass*collisionDelta*tangentSpeed +
				2*reducedMass*collisionDelta*collisionDelta*tangentNormSquared;
			const double oldMomentumX = pi->m*verifyOldVix + pj->m*verifyOldVjx;
			const double oldMomentumY = pi->m*verifyOldViy + pj->m*verifyOldVjy;
			const double oldMomentumZ = pi->m*verifyOldViz + pj->m*verifyOldVjz;
			const double newMomentumX = pi->m*pi->vx + pj->m*pj->vx;
			const double newMomentumY = pi->m*pi->vy + pj->m*pj->vy;
			const double newMomentumZ =
				#if THREE_D
				pi->m*pi->vz + pj->m*pj->vz
				#else
				0
				#endif
			;
			const double scale = 1 + fabs(oldEnergy) + fabs(newEnergy) +
				fabs(oldMomentumX) + fabs(oldMomentumY) + fabs(oldMomentumZ);
			int invalid = !isfinite(newEnergy) ||
				fabs(newGx - expectedGx) > 1e-7*scale ||
				fabs(newGy - expectedGy) > 1e-7*scale ||
				fabs(newGz - expectedGz) > 1e-7*scale ||
				fabs(newEnergy - oldEnergy - expectedEnergyChange) > 1e-7*scale ||
				fabs(newMomentumX - oldMomentumX) > 1e-9*scale ||
				fabs(newMomentumY - oldMomentumY) > 1e-9*scale ||
				fabs(newMomentumZ - oldMomentumZ) > 1e-9*scale ||
				fabs((collTermX - verifyOldVirialX) - impulseX*dx) > 1e-9*scale ||
				fabs((collTermY - verifyOldVirialY) - impulseY*dy) > 1e-9*scale ||
				fabs((collTermXY - verifyOldVirialXY) - impulseX*dy) > 1e-9*scale ||
				fabs((collTermYX - verifyOldVirialYX) - impulseY*dx) > 1e-9*scale;
			#if THREE_D
			invalid = invalid ||
				fabs((collTermZ - verifyOldVirialZ) - impulseZ*dz) > 1e-9*scale ||
				fabs((collTermXZ - verifyOldVirialXZ) - impulseX*dz) > 1e-9*scale ||
				fabs((collTermZX - verifyOldVirialZX) - impulseZ*dx) > 1e-9*scale ||
				fabs((collTermYZ - verifyOldVirialYZ) - impulseY*dz) > 1e-9*scale ||
				fabs((collTermZY - verifyOldVirialZY) - impulseZ*dy) > 1e-9*scale;
			#endif
			if (invalid){
				fprintf(stderr, "Tangential collision accounting invariant failed at t=%.17g for %d/%d\n",
					t, i, j);
				exit(3);
			}
		}
		#endif
		double collTermXLocal = impulseX*dx;
		double collTermYLocal = impulseY*dy;
		double collTermXYLocal = impulseX*dy;
		double collTermYXLocal = impulseY*dx;
		if (thermoPressureField){
			pressure_field_add_collision(pi, dx, dy, collTermXLocal, collTermYLocal, collTermXYLocal, collTermYXLocal);
		}

	}
	else if (addEnergy){
		#if THREE_D
		double b = dx*dvx + dy*dvy + dz*dvz;
		#else
		double b = dx*dvx + dy*dvy;
		#endif
		double distSquared = DIST2(pi->rad, pj->rad);
		double dti = t - pi->lastColl;
		double dtj = t - pj->lastColl;
		
		double deltaETotal;
		if (simplifyInjection)
			deltaETotal = 2*additionalEnergy + deltaE*(2 + sign(dti - taur) + sign(dtj - taur));
		else
			deltaETotal = 2*additionalEnergy + deltaE*(pow(1 - exp(-dtj/taur), beta) + pow(1 - exp(-dti/taur), beta));
		if (randomInjection){
			deltaETotal *= drand(0.9, 1./0.9);
		}

		double funkyFactor2 = (b - sqrt(effectiveRestitution*effectiveRestitution*b*b +
			2*distSquared*deltaETotal*(pi->m + pj->m)/(pi->m*pj->m)))/(2*distSquared);
		double dviX = 2*pj->m*invMass*funkyFactor2*dx;
		double dviY = 2*pj->m*invMass*funkyFactor2*dy;
		double impulseX = pi->m*dviX;
		double impulseY = pi->m*dviY;
		#if THREE_D
		double dviZ = 2*pj->m*invMass*funkyFactor2*dz;
		double impulseZ = pi->m*dviZ;
		#else
		double impulseZ = 0;
		#endif
		double collTermXLocal = impulseX*dx;
		double collTermYLocal = impulseY*dy;
		double collTermXYLocal = impulseX*dy;
		double collTermYXLocal = impulseY*dx;
		accumulateCollisionVirial(impulseX, impulseY,
			#if THREE_D
			impulseZ, dx, dy, dz
			#else
			impulseZ, dx, dy, 0
			#endif
		);
		if (thermoPressureField){
			pressure_field_add_collision(pi, dx, dy, collTermXLocal, collTermYLocal, collTermXYLocal, collTermYXLocal);
		}
		//double E = 0.5*(pi->m*(pi->vx*pi->vx + pi->vy*pi->vy) + pj->m*(pj->vx*pj->vx + pj->vy*pj->vy ));
		pi->vx += dviX;
		pi->vy += dviY;
		pj->vx -= 2*pi->m*invMass*funkyFactor2*dx;
		pj->vy -= 2*pi->m*invMass*funkyFactor2*dy;
		#if THREE_D
		pi->vz += dviZ;
		pj->vz -= 2*pi->m*invMass*funkyFactor2*dz;
		#endif
		//printf("%lf %lf\n", deltaETotal, - E + 0.5*(pi->m*(pi->vx*pi->vx + pi->vy*pi->vy) + pj->m*(pj->vx*pj->vx + pj->vy*pj->vy)));
	}
	else if (addActiveDelta){
		double rij2 = dx*dx + dy*dy;
		#if THREE_D
		rij2 += dz*dz;
		#endif
		double rij = sqrt(rij2);
		double nx = dx/rij;
		double ny = dy/rij;
		double gdotn = dvx*nx + dvy*ny;
		#if THREE_D
		double nz = dz/rij;
		gdotn += dvz*nz;
		double dvt_z = dvz - gdotn*nz;
		#endif
		double dvt_x = dvx - gdotn*nx;
		double dvt_y = dvy - gdotn*ny;

		double dvi_x = funkyFactor*pj->m*dx + 2*delta*pj->m*invMass*dvt_x;
		double dvi_y = funkyFactor*pj->m*dy + 2*delta*pj->m*invMass*dvt_y;
		double dvj_x = -funkyFactor*pi->m*dx - 2*delta*pi->m*invMass*dvt_x;
		double dvj_y = -funkyFactor*pi->m*dy - 2*delta*pi->m*invMass*dvt_y;
		#if THREE_D
		double dvi_z = funkyFactor*pj->m*dz + 2*delta*pj->m*invMass*dvt_z;
		double dvj_z = -funkyFactor*pi->m*dz - 2*delta*pi->m*invMass*dvt_z;
		#endif

		pi->vx += dvi_x;
		pi->vy += dvi_y;
		pj->vx += dvj_x;
		pj->vy += dvj_y;
		#if THREE_D
		pi->vz += dvi_z;
		pj->vz += dvj_z;
		#endif

		accumulateCollisionVirial(pi->m*dvi_x, pi->m*dvi_y,
			#if THREE_D
			pi->m*dvi_z, dx, dy, dz
			#else
			0, dx, dy, 0
			#endif
		);
	}
	else {
	#endif
		double collTerm = pi->m*pj->m*funkyFactor;
		double collTermXLocal = collTerm*dx*dx;
		double collTermYLocal = collTerm*dy*dy;
		double collTermXYLocal = collTerm*dx*dy;
		double collTermYXLocal = collTerm*dy*dx;

		accumulateCollisionVirial(collTerm*dx, collTerm*dy,
			#if THREE_D
			collTerm*dz, dx, dy, dz
			#else
			0, dx, dy, 0
			#endif
		);
		if (thermoPressureField){
			pressure_field_add_collision(pi, dx, dy, collTermXLocal, collTermYLocal, collTermXYLocal, collTermYXLocal);
		}
		

		pi->vx += funkyFactor*pj->m*dx;
		pi->vy += funkyFactor*pj->m*dy;
		pj->vx -= funkyFactor*pi->m*dx;
		pj->vy -= funkyFactor*pi->m*dy; 
		#if THREE_D
		pi->vz += funkyFactor*pj->m*dz;
		pj->vz -= funkyFactor*pi->m*dz;
		#endif
		
		
		#if TANGENTIAL
		dist = sqrt(dx*dx + dy*dy);
		double dxr = dx/dist;
		double dyr = dy/dist;
		double gnx = (dvx*dxr + dvy*dyr)*dxr;
		double gny = (dvx*dxr + dvy*dyr)*dyr;
		double gtx = dvx - gnx + dist/2*(dyr*(pi->omega + pj->omega));
		double gty = dvy - gny - dist/2*(dxr*(pi->omega + pj->omega));
		double q = pi->J/(pi->m*dist*dist/4.);
		

		
		
		pi->vx += (1 - resBeta)/(2*(1 + 1/q))*gtx;
		pi->vy += (1 - resBeta)/(2*(1 + 1/q))*gty;
		pj->vx -= (1 - resBeta)/(2*(1 + 1/q))*gtx;
		pj->vy -= (1 - resBeta)/(2*(1 + 1/q))*gty;
		if (addDeltaTangent){
			double dxrTemp = dxr;
			dxr = dyr;  
			dyr = -dxrTemp;
			pi->vx += - 2*pj->m*invMass*delta*dxr;
			pi->vy += - 2*pj->m*invMass*delta*dyr;
			pj->vx -= - 2*pi->m*invMass*delta*dxr;
			pj->vy -= - 2*pi->m*invMass*delta*dyr;
		}

		pi->omega += (1 - resBeta)/(1 + q)*1/(dist)*(gty*dxr - gtx*dyr) + deltaOmega;
		pj->omega += (1 - resBeta)/(1 + q)*1/(dist)*(gty*dxr - gtx*dyr) + deltaOmega;
		#endif
	#if !TANGENTIAL
	}
	#endif

	if (thermoPressureField){
		pressure_field_add_kinetic(pi);
		pressure_field_add_kinetic(pj);
	}
	#if (STRESS || HEAT)
	updateParticleTransport(pi, oldvix, oldviy, oldviz);
	updateParticleTransport(pj, oldvjx, oldvjy, oldvjz);
	#endif
	#if HEAT
	double collisionalHeatTransfer = recordPairHeat(pi, pj,
		oldvix, oldviy, oldviz, oldvjx, oldvjy, oldvjz, dx, dy,
		#if THREE_D
		dz
		#else
		0
		#endif
		, 1
	);
	if (tempProfileThermo && tempSlabWidth > 0){
		double meanX = 0.5*(pi->x + pj->x);
		PBCpostX(&meanX);
		int indexSlab = (int)(meanX/tempSlabWidth);
		if (indexSlab < 0)
			indexSlab = 0;
		else if (indexSlab >= tempProfileNSlabs)
			indexSlab = tempProfileNSlabs - 1;
		JxProfile[indexSlab] += dx*collisionalHeatTransfer/dtimeThermo;
		JyProfile[indexSlab] += dy*collisionalHeatTransfer/dtimeThermo;
		#if THREE_D
		JzProfile[indexSlab] += dz*collisionalHeatTransfer/dtimeThermo;
		#endif
	}
	#endif

    pi->lastColl = t;
    pj->lastColl = t;
	//recomputes crossing and collision of pi and pj
	removeEventFromQueue(EVENT_NODE(i));
	removeEventFromQueue(EVENT_NODE(j));
	crossingEventNormal(i);
	crossingEventNormal(j);

	removeEventFromQueue(EVENT_NODE(N + j));
	collisionEventNormal(i);
	collisionEventNormal(j);

}


#if INTRUDER

void doTheCollisionIntruder(){


	int i = nextEvent->i;
	int j = nextEvent->j;

	intruder* P = &intru;
	particle* p;

	if (j == -1){
		p = particles + i;
		freeFlyNormal(p);
		freeFlyIntruder(P);
		if (nextEvent->collActual != P->coll){
			collisionEventNormal(i);
			return;
		}
	}
	else{
		p = particles + j;
		freeFlyNormal(p);
		freeFlyIntruder(P);
		if (nextEvent->collActual != p->coll){
			collisionEventIntruder();
			return;
		}
	}

	P->coll++;
	p->coll++;
	ncol++;

	double px = p->x;
	double py = p->y;
	polygon* poly = &P->real_shape;

	double bestDist2 = never;
	double nx = 0.0;
	double ny = 0.0;
	int bestEdge = -1;

	for (int e = 0; e < poly->n_vertices; e++){
		int e2 = (e + 1) % poly->n_vertices;

		double ax = poly->vertices[e][0];
		double ay = poly->vertices[e][1];
		double abx = poly->vertices[e2][0] - ax;
		double aby = poly->vertices[e2][1] - ay;
		PBC(&abx, &aby);

		double apx = px - ax;
		double apy = py - ay;
		PBC(&apx, &apy);

		double ab2 = abx*abx + aby*aby;
		if (ab2 <= 0.0)
			continue;

		double s = (apx*abx + apy*aby)/ab2;
		if (s < 0.0) s = 0.0;
		if (s > 1.0) s = 1.0;

		double qx = ax + s*abx;
		double qy = ay + s*aby;

		double dx = px - qx;
		double dy = py - qy;
		PBC(&dx, &dy);

		double dist2 = dx*dx + dy*dy;
		if (dist2 < bestDist2){
			bestDist2 = dist2;
			bestEdge = e;
			if (dist2 > 1e-24){
				double dist = sqrt(dist2);
				nx = dx/dist;
				ny = dy/dist;
			}
		}
	}

	if (bestEdge >= 0 && bestDist2 <= 1e-24){
		nx = poly->normals[bestEdge][0];
		ny = poly->normals[bestEdge][1];
		double n2 = nx*nx + ny*ny;
		if (n2 > 0.0){
			double invn = 1.0/sqrt(n2);
			nx *= invn;
			ny *= invn;
		}
	}

	if (bestEdge < 0 || (nx == 0.0 && ny == 0.0)){
		printf("Error: no normal found for collision!\n");
		exit(3);
	}

	p->normalAtCollision[0] = nx;
	p->normalAtCollision[1] = ny;

	double contactX = px - p->rad*nx;
	double contactY = py - p->rad*ny;

	double rx = contactX - P->x;
	double ry = contactY - P->y;
	PBC(&rx, &ry);

	double vcx = P->vx - P->omega*ry;
	double vcy = P->vy + P->omega*rx;

	double rvx = p->vx - vcx;
	double rvy = p->vy - vcy;
	double vn = rvx*nx + rvy*ny; // g_n

	// Tangent at contact: t = z x n
	double tx = -ny;
	double ty =  nx;
	double vt = rvx*tx + rvy*ty; // g_t (currently unused by the chosen tangential rule)
	(void)vt;

	if (vn < vn_tol){

		double invPMass = 1.0/p->m;
		double invMassIntruder = 1.0/P->M;
		double Im = P->Im;

		double invInertiaIntruder = 1.0/Im;

		double aeff = invPMass + invMassIntruder;
		double kappa_n = rx*ny - ry*nx;      // (r x n)_z
		double kappa_t = rx*ty - ry*tx;      // (r x t)_z = r . n in 2D

		double lambda_n = aeff + (kappa_n*kappa_n)*invInertiaIntruder;
		double lambda_t = aeff + (kappa_t*kappa_t)*invInertiaIntruder;

		if (lambda_n > 0.0 && lambda_t > 0.0){
			// Tangential impulse (chosen model): J_t = 2*Delta / lambda_t
			double Jt = (2.0*delta_intruder)/lambda_t;

			// Normal impulse with coupling:
			// J_n = (1+alpha)/lambda_n * g_n - (kappa_n*kappa_t)/(I*lambda_n) * J_t
			double Jn = ((1.0 + res_intruder)/lambda_n)*vn
					- (kappa_n*kappa_t)*invInertiaIntruder*Jt/lambda_n;

			// Full impulse vector J = J_n n + J_t t
			double Jx = Jn*nx + Jt*tx;
			double Jy = Jn*ny + Jt*ty;
			intruderImpulseX += Jx;
			intruderImpulseY += Jy;
			intruderTorqueZ += (rx*Jy - ry*Jx);

			double cphi = cos(P->angle);
			double sphi = sin(P->angle);
			double Jbx = Jx*cphi + Jy*sphi;
			double Jby = -Jx*sphi + Jy*cphi;

			
			// --- Online KM estimator update (collision-to-collision, lag=1) ---
			if (km.enabled && km.has_prev){
			    double dtFlight = t - km.last_t;
			    // Only use intervals that start after the burn time (stationary sampling)
			    if (dtFlight > 1e-12 && km.last_t >= km.burn_time){
			        double int_vbx = 0.0, int_vby = 0.0;
			        km_integrate_vbody(km.vx0, km.vy0, km.phi0, km.om0, dtFlight, &int_vbx, &int_vby);
			
			        // Design vector uses integrals directly: x = [dt, ∫Vb_x dt, ∫Vb_y dt, ∫Ω dt]
			        double xvec[4];
			        xvec[0] = dtFlight;
			        xvec[1] = int_vbx;
			        xvec[2] = int_vby;
			        xvec[3] = km.om0 * dtFlight;
			
			        // Response is the inertial body increment at the collision: y = [ΔVb_x, ΔVb_y, ΔΩ]
			        double yvec[3];
			        yvec[0] = Jbx * invMassIntruder;                 // J_body/M
			        yvec[1] = Jby * invMassIntruder;
			        yvec[2] = (rx*Jy - ry*Jx) * invInertiaIntruder;  // (r×J)/I
			
			        km.nsamples++;
			        km.sum_dt += dtFlight;
				    double Etrans = 0.5 * P->M * (km.vx0*km.vx0 + km.vy0*km.vy0);
				    double Erot = 0.5 * P->Im * (km.om0*km.om0);
				    km.sum_Etrans_dt += Etrans * dtFlight;
				    km.sum_Erot_dt += Erot * dtFlight;
			        km.sum_dU[0] += yvec[0];
			        km.sum_dU[1] += yvec[1];
			        km.sum_dU[2] += yvec[2];
			
			        km.sum_intVb[0] += int_vbx;
			        km.sum_intVb[1] += int_vby;
			        km.sum_intOm_dt += km.om0 * dtFlight;
			        km.sum_intOmVb[0] += km.om0 * int_vbx;
			        km.sum_intOmVb[1] += km.om0 * int_vby;
			
			        double w = 1.0/dtFlight; // weights = 1/dt
			        for (int ii = 0; ii < 4; ++ii){
			            for (int jj = 0; jj < 4; ++jj){
			                km.XtWX[ii][jj] += w * xvec[ii] * xvec[jj];
			                km.XTX[ii][jj]  +=     xvec[ii] * xvec[jj];
			            }
			            for (int kk = 0; kk < 3; ++kk){
			                km.XtWY[ii][kk] += w * xvec[ii] * yvec[kk];
			                km.XTY[ii][kk]  +=     xvec[ii] * yvec[kk];
			            }
			        }
			        for (int kk = 0; kk < 3; ++kk){
			            for (int ll = 0; ll < 3; ++ll){
			                km.YTY[kk][ll] += yvec[kk] * yvec[ll];
			            }
			        }
			    }
			}
			
			// Update velocities using:
			// v' = v - J/m, V' = V + J/M, Omega' = Omega + (r x J)_z / I
			p->vx -= Jx*invPMass;
			p->vy -= Jy*invPMass;

			P->vx += Jx*invMassIntruder;
			P->vy += Jy*invMassIntruder;
			P->omega += (rx*Jy - ry*Jx)*invInertiaIntruder;
			if (rethermalize_after_intruder_collision){
				rethermalize_bath_particle_at_intruder_contact(p, P, contactX, contactY, rx, ry, nx, ny, tx, ty, T);
			}
			// Update KM baseline state at the collision time (start of next free flight)
			// MUST be inside the collision block: non-collision events (vn >= vn_tol)
			// must NOT reset the baseline, otherwise flight intervals get fragmented
			// and sum_intVb / sum_dt (= Ubar) becomes biased.
			if (km.enabled){
			    km.last_t = t;
			    km.phi0 = P->angle;
			    km.vx0 = P->vx;
			    km.vy0 = P->vy;
			    km.om0 = P->omega;
			    km.has_prev = 1;
			}
		}
	}
	
	p->lastColl = t;


	removeEventFromQueue(EVENT_NODE(p->num));
	crossingEvent(p->num);

	if (j == -1){
		removeEventFromQueue(EVENT_NODE(2*N + 8));
	}
	else{
		removeEventFromQueue(EVENT_NODE(N + j));
	}
	collisionEventNormal(p->num);
	collisionEventIntruder();
}



#endif

static double potentialPairStepEnergy(const particle* first,
	const particle* second){
	if (charged){
		int chargeProduct = first->charge*second->charge;
		if (chargeProduct < 0)
			return -fabs(U);
		if (chargeProduct > 0)
			return fabs(U);
		return 0;
	}
	if (liquidliquid && first->type != second->type)
		return 0;
	return U;
}

bool potentialStepNormalVelocities(double firstMass, double secondMass,
	double firstVelocity, double secondVelocity, double deltaPotential,
	double* firstVelocityAfter, double* secondVelocityAfter){
	if (firstVelocityAfter == NULL || secondVelocityAfter == NULL ||
		!isfinite(firstMass) || !isfinite(secondMass) ||
		firstMass <= 0 || secondMass <= 0){
		return false;
	}
	double totalMass = firstMass + secondMass;
	double reducedMass = firstMass*secondMass/totalMass;
	double momentum = firstMass*firstVelocity + secondMass*secondVelocity;
	double relativeVelocity = firstVelocity - secondVelocity;
	double relativeSquaredAfter = relativeVelocity*relativeVelocity -
		2*deltaPotential/reducedMass;
	/* At the exact barrier height the transmitted radial speed would be zero
	 * on the discontinuity. Preserve the former strict-threshold behavior and
	 * reflect, avoiding a pair stranded on the shell. */
	bool transmitted = isfinite(relativeSquaredAfter) &&
		relativeSquaredAfter > 0;
	double relativeVelocityAfter = transmitted ?
		copysign(sqrt(relativeSquaredAfter), relativeVelocity) :
		-relativeVelocity;
	*firstVelocityAfter = momentum/totalMass +
		secondMass*relativeVelocityAfter/totalMass;
	*secondVelocityAfter = momentum/totalMass -
		firstMass*relativeVelocityAfter/totalMass;
	return transmitted;
}

void doIn(){

	int i = nextEvent->i;
	int j = nextEvent->j;

	particle* pi = particles + i;
	particle* pj = particles + j;


	freeFly(pi); //might put it on top?
	if (nextEvent->collActual != pj->coll){ //the collision trick of the article
        collisionEvent(i);
        return;
    }
	ncol++;



	freeFly(pj);


	pi->coll++;
	pj->coll++;
	#if (STRESS || HEAT)
	advanceTransportAccumulators();
	double oldvix = pi->vx;
	double oldviy = pi->vy;
	double oldvjx = pj->vx;
	double oldvjy = pj->vy;
	#if THREE_D
	double oldviz = pi->vz;
	double oldvjz = pj->vz;
	#else
	double oldviz = 0;
	double oldvjz = 0;
	#endif
	#endif

	double stepEnergy = potentialPairStepEnergy(pi, pj);
	if (stepEnergy != 0){

		double dx = pj->x - pi->x;
		double dy = pj->y - pi->y;
		#if THREE_D
		double dz =  pj->z - pi->z;
		PBC(&dx, &dy, &dz);
		double dr = sqrt(dx*dx + dy*dy + dz*dz);
		double dzr = dz/dr;
		#else
		PBC(&dx, &dy);
		double dr = sqrt(dx*dx + dy*dy);
		#endif

		double dxr = dx/dr;
		double dyr = dy/dr;


		//https://scholarworks.rit.edu/cgi/viewcontent.cgi?article=5982&context=theses
		#if THREE_D
		double vi = pi->vx*dxr + pi->vy*dyr + pi->vz*dzr;
		double vj = pj->vx*dxr + pj->vy*dyr + pj->vz*dzr;
		#else
		double vi = pi->vx*dxr + pi->vy*dyr;
		double vj = pj->vx*dxr + pj->vy*dyr;
		#endif

		double mi = pi->m;
		double mj = pj->m;
		double vif, vjf;
		bool transmitted = potentialStepNormalVelocities(mi, mj, vi, vj,
			stepEnergy, &vif, &vjf);

		if (transmitted){
			addParticleInWellList(pi, j);
			addParticleInWellList(pj, i);
		}

		
		/* The virial helper takes the impulse itself and multiplies by r.
		   Keep the factor of dr in the separation, not in the impulse. */
		double impulseMagnitude = -pj->m*(vjf - vj);
		accumulateCollisionVirial(impulseMagnitude*dxr, impulseMagnitude*dyr,
			#if THREE_D
			impulseMagnitude*dzr, dx, dy, dz
			#else
			0, dx, dy, 0
			#endif
		);

		double dvjx = (vjf - vj)*dxr;
		double dvjy = (vjf - vj)*dyr;
		#if THREE_D
		double dvjz = (vjf - vj)*dzr;
		pi->vz += (vif - vi)*dzr;
		pj->vz += dvjz;
		#endif
		pi->vx += (vif - vi)*dxr;
		pi->vy += (vif - vi)*dyr;
		pj->vx += dvjx;
		pj->vy += dvjy;
		#if (STRESS || HEAT)
		updateParticleTransport(pi, oldvix, oldviy, oldviz);
		updateParticleTransport(pj, oldvjx, oldvjy, oldvjz);
		#endif
		#if HEAT
		recordPairHeat(pi, pj, oldvix, oldviy, oldviz, oldvjx, oldvjy, oldvjz,
			dx, dy,
			#if THREE_D
			dz
			#else
			0
			#endif
			, 0
		);
		#endif
	}
	else{
		addParticleInWellList(pi, j);
		addParticleInWellList(pj, i);
	}


	//recomputes crossing and collision of pi and pj
	removeEventFromQueue(EVENT_NODE(i));
	removeEventFromQueue(EVENT_NODE(j));
	crossingEvent(i);
	crossingEvent(j);

	removeEventFromQueue(EVENT_NODE(N + j));
	collisionEvent(i);
	collisionEvent(j);



}

//Leave the square well or bounce on it
void doOut(){

	int i = nextEvent->i;
	int j = nextEvent->j;

	particle* pi = particles + i;
	particle* pj = particles + j;



	freeFly(pi); //might put it on top?
	if (nextEvent->collActual != pj->coll){ //the collision trick of the article
        collisionEvent(i);
        return;
    }
	ncol++;

	freeFly(pj);

	pi->coll++;
	pj->coll++;
	#if (STRESS || HEAT)
	advanceTransportAccumulators();
	double oldvix = pi->vx;
	double oldviy = pi->vy;
	double oldvjx = pj->vx;
	double oldvjy = pj->vy;
	#if THREE_D
	double oldviz = pi->vz;
	double oldvjz = pj->vz;
	#else
	double oldviz = 0;
	double oldvjz = 0;
	#endif
	#endif

	double stepEnergy = potentialPairStepEnergy(pi, pj);
	if (stepEnergy != 0){
		double dx = pj->x - pi->x;
		double dy = pj->y - pi->y;
		#if THREE_D
		double dz =  pj->z - pi->z;
		PBC(&dx, &dy, &dz);
		double dr = sqrt(dx*dx + dy*dy + dz*dz);
		double dzr = dz/dr;
		#else
		PBC(&dx, &dy);
		double dr = sqrt(dx*dx + dy*dy);
		#endif

		double dxr = dx/dr;
		double dyr = dy/dr;


		// https://scholarworks.rit.edu/cgi/viewcontent.cgi?article=5982&context=theses
		#if THREE_D
		double vi = pi->vx*dxr + pi->vy*dyr + pi->vz*dzr;
		double vj = pj->vx*dxr + pj->vy*dyr + pj->vz*dzr;
		#else
		double vi = pi->vx*dxr + pi->vy*dyr;
		double vj = pj->vx*dxr + pj->vy*dyr;
		#endif

		double mi = pi->m;
		double mj = pj->m;
		double vif, vjf;
		bool transmitted = potentialStepNormalVelocities(mi, mj, vi, vj,
			-stepEnergy, &vif, &vjf);
		if (transmitted){
			removeParticleInWellList(pi, j);
			removeParticleInWellList(pj, i);
		}

		double impulseMagnitude = -pj->m*(vjf - vj);
		accumulateCollisionVirial(impulseMagnitude*dxr, impulseMagnitude*dyr,
			#if THREE_D
			impulseMagnitude*dzr, dx, dy, dz
			#else
			0, dx, dy, 0
			#endif
		);

		double dvjx = (vjf - vj)*dxr;
		double dvjy = (vjf - vj)*dyr;
		#if THREE_D
		double dvjz = (vjf - vj)*dzr;
		pi->vz += (vif - vi)*dzr;
		pj->vz += dvjz;
		#endif
		pi->vx += (vif - vi)*dxr;
		pi->vy += (vif - vi)*dyr;
		pj->vx += dvjx;
		pj->vy += dvjy;
		#if (STRESS || HEAT)
		updateParticleTransport(pi, oldvix, oldviy, oldviz);
		updateParticleTransport(pj, oldvjx, oldvjy, oldvjz);
		#endif
		#if HEAT
		recordPairHeat(pi, pj, oldvix, oldviy, oldviz, oldvjx, oldvjy, oldvjz,
			dx, dy,
			#if THREE_D
			dz
			#else
			0
			#endif
			, 0
		);
		#endif
	}
	else{
		removeParticleInWellList(pi, j);
		removeParticleInWellList(pj, i);
	}
	//recomputes crossing and collision of pi and pj
	removeEventFromQueue(EVENT_NODE(i));
	removeEventFromQueue(EVENT_NODE(j));
	crossingEvent(i);
	crossingEvent(j);

	removeEventFromQueue(EVENT_NODE(N + j));
	collisionEvent(i);
	collisionEvent(j);
}

/* ---------------------------------------------
	Creates the screenshot event and adds it
	to the BST.
--------------------------------------------- */
void addEventScreenshot(double tscreen){
	node* toAdd;

	toAdd = EVENT_NODE(2*N + 1);
	toAdd->type = SCREENSHOT;
	toAdd->t = tscreen;
	toAdd->j = 0; //whatever.

	addEventToQueue(toAdd);
}

void addEventThermo(double tscreen){
	node* toAdd;

	toAdd = EVENT_NODE(2*N + 3);
	toAdd->type = THERMO;
	toAdd->t = tscreen;
	toAdd->j = 0; //whatever.

	addEventToQueue(toAdd);
}

/* ---------------------------------------------
	Updates the position of the particles
	at actual time and take a screenshot of
	the system. Also schedules a new
	screenshot event.
--------------------------------------------- */
#if EDMD_VERIFY_OVERLAPS
static void verifyNoOverlaps(void){
	#if POLYGON_PARTICLES && !THREE_D
	if (polygon_particles_active()){
		char message[512];
		const double tolerance = 1.0e-10*
			(1.0 + fabs(Lx) + fabs(Ly));
		if (!polygon_particles_verify_no_overlaps(
				t, tolerance, message, sizeof(message))){
			#if G
				simulationReportRecoverableOverlapFault(message);
				return;
			#else
				fprintf(stderr, "%s\n", message);
				exit(3);
			#endif
		}
		return;
	}
	#endif
	for (int i = 0; i < N; i++){
		if (addCircularWall){
			double dxWall = particles[i].x - halfLx;
			double dyWall = particles[i].y - halfLy;
			double distance2 = dxWall*dxWall + dyWall*dyWall;
			#if THREE_D
			double dzWall = particles[i].z - halfLz;
			distance2 += dzWall*dzWall;
			#endif
			double accessibleRadius = halfLx - particles[i].rad;
			if (distance2 > accessibleRadius*accessibleRadius + 1e-8){
				#if G
				simulationSetFault(
					"Confinement verification found particle %d outside the wall at t=%.9g. "
					"The simulation is locked until restart.", i, t);
				return;
				#else
				fprintf(stderr,
					"Confinement invariant failed at t=%.17g for particle %d: %.17g > %.17g\n",
					t, i, distance2, accessibleRadius*accessibleRadius);
				exit(3);
				#endif
			}
		}
		for (int j = i + 1; j < N; j++){
			double dx = particles[j].x - particles[i].x;
			double dy = particles[j].y - particles[i].y;
			#if THREE_D
				double dz = particles[j].z - particles[i].z;
				PBC(&dx, &dy, &dz);
				double separation2 = dx*dx + dy*dy + dz*dz;
			#else
				PBC(&dx, &dy);
				double separation2 = dx*dx + dy*dy;
			#endif
			double contact2 = DIST2(particles[i].rad, particles[j].rad);
			if (separation2 < contact2 - 1e-8){
				#if G
				simulationSetFault(
					"Overlap verification found particles %d and %d intersecting at t=%.9g. "
					"The simulation is locked until restart.", i, j, t);
				return;
				#else
				fprintf(stderr,
					"Overlap invariant failed at t=%.17g for particles %d and %d: %.17g < %.17g\n",
					t, i, j, separation2, contact2);
				exit(3);
				#endif
			}
		}
	}
}
#endif

void takeAScreenshot(){

	for (int i = 0; i < N; i++){
		freeFly(particles + i);
	}
	#if EDMD_VERIFY_OVERLAPS
		verifyNoOverlaps();
		#if G
		if (simulationHasFault())
			return;
		#endif
	#endif
	#if INTRUDER
	freeFlyIntruder(&intru);
	#endif
	physicalQ();
	if (periodicDumpStreamOpen && dtime > 0 && dtime < never){
		addEventScreenshot(t + dtime);
	}
	if (periodicDumpStreamOpen && fichier != NULL){
		saveTXT();
	}
	
	printInfo();
	fflush(stdout);
}

void takeAThermo(){
	if (outputSamplingNeedsSynchronizedParticles()){
		for (int i = 0; i < N; i++){
			freeFly(particles + i);
		}
	}
	if (msd && posxInitial == NULL){
			posxInitial = calloc((unsigned int)N, sizeof(double));
			posyInitial = calloc((unsigned int)N, sizeof(double));
			#if THREE_D
			poszInitial = calloc((unsigned int)N, sizeof(double));
			#endif
			if (posxInitial == NULL || posyInitial == NULL
				#if THREE_D
				|| poszInitial == NULL
				#endif
			){
				free(posxInitial);
				free(posyInitial);
				#if THREE_D
				free(poszInitial);
				poszInitial = NULL;
				#endif
				posxInitial = NULL;
				posyInitial = NULL;
				reportOutputRuntimeFailure(
					"Unable to allocate mean-squared-displacement reference positions");
				#if G
				msd = 0;
				#endif
			}
			else{
			for (int i = 0; i < N; i++){
				COM_x += particles[i].x + particles[i].crossX*Lx;
				COM_y += particles[i].y + particles[i].crossY*Ly;
				posxInitial[i] = particles[i].x + particles[i].crossX*Lx;
				posyInitial[i] = particles[i].y + particles[i].crossY*Ly;
				#if THREE_D
				COM_z += particles[i].z + particles[i].crossZ*Lz;
				poszInitial[i] = particles[i].z + particles[i].crossZ*Lz;
				#endif
			}
			COM_x = COM_x/N;
			COM_y = COM_y/N;
			#if THREE_D
			COM_z = COM_z/N;
			#endif
			}
	}
	physicalQ();
	saveThermo();
	saveTempProfile();

	double nextTime = t + dtimeThermo;
	if (logSpacing){
		nextTime = t + dtimeThermo*pow(base, outputSamples.logarithmicIndex);
		outputSamples.logarithmicIndex++;
	}
	if (outputThermoSamplingRequested()){
		addEventThermo(nextTime);
	}
	printInfo();

}

void printInfo(){
	int size = 29;
	const bool indefinite = tmax >= never;
	int place = indefinite ? 0 : (int)(size*t/tmax);
	if (indefinite)
		printf("\r│\033[1;32m%.2e\033[0;37m (unbounded) ", t);
	else
		printf("\r│\033[1;32m%.2e\033[0;37m/%.2e ", t, tmax);
	CCYAN;
	for (int i = 0; i < size; i++){
		if (indefinite)
			printf("·");
		else if (i < place - 1)
			printf("█");
		else if (i == (place - 1)){
			if (i == size - 1)
				printf("█");
			else
				printf("🭬");
			}
		else
			printf(" ");
	}
	CWHITE;
	


	printf(" │ n° coll = \033[1;32m%.2e\033[0;37m  Energy = \033[1;32m%.8lf\033[0;37m  Pressure = \033[1;32m%.3lf\033[0;37m │", (float)ncol, E/N, pressure);

}

void stopGrow(){
	particle* p;
	for (int i = 0; i < N; i++){
		p = particles + i;
		freeFly(p);
		p->coll = 0;
		ncol = 0;
		ncross = 0;
		p->crossX = 0;
		p->crossY = 0;
		#if THREE_D
		p->crossZ = 0;
		#endif
		//#p->rad = 0.001;
	}
	#if POLYGON_PARTICLES && !THREE_D
	if (polygon_particles_growing()){
		/* Growth CCD intentionally stops at a small conservative pre-contact
		 * tolerance.  Pay for one exact local SAT pass at handoff, rather than
		 * machine-precision bisection at every preparation collision. */
		char verification[512];
		if (!polygon_particles_verify_local_no_overlaps(t, 1.0e-10,
			verification, sizeof(verification))){
			simulationReportFault(verification[0] != '\0' ? verification :
				"Polygon preparation failed its final overlap check because the final linked-cell state is invalid.");
			return;
		}
		polygon_particles_finish_growth();
		#if G
		releaseGraphicalRecoveryCheckpoints();
		#endif
	}
	#endif
	selectNormalDynamicsHandlers();
	normalizePhysicalQ();
	
	if (addWell && U != 0.0){
		for (int i  = 0; i < N; i++){
			p = particles + i;
			free(p->particlesInWell);
		} 
		pairsInit();
	}
	

	
	for (int i = 0; i < N; i++){
		p = particles + i;
		p->coll++;
		removeEventFromQueue(EVENT_NODE(N + p->num));
		removeEventFromQueue(EVENT_NODE(p->num));
		crossingEvent(p->num);
		collisionEvent(p->num);
	}
	#if INTRUDER
	collisionEventIntruder();
	#endif
	#if (STRESS || HEAT)
	initializeTransportAccumulators();
	#endif
	resetOutputSamplingState();
	#if G && POLYGON_PARTICLES && !THREE_D
	/* Discard preparation-era pressure/energy histories and wall-impulse
	 * baselines. The first production sample establishes fresh observables. */
	live_analysis_reset_generation();
	#endif
}

void addEventGrow(double time){
	node* toAdd;
	toAdd = EVENT_NODE(2*N + 5);
	toAdd->type = GROWSTOP;
	toAdd->t = time;
	toAdd->j = 0; 
	addEventToQueue(toAdd);
}
/* ---------------------------------------------
	Creates the event responsible for
	adding energy into the system
--------------------------------------------- */
void addEventNoise(double tnoise){
	node* toAdd;
	toAdd = EVENT_NODE(2*N + 2);
	toAdd->type = ADDINGNOISE;
	toAdd->t = tnoise;
	toAdd->j = 0; //whatever.

	addEventToQueue(toAdd);
}

/* ---------------------------------------------
	Kicks every particles of the system
--------------------------------------------- */

bool velocityRescaleFactor(double kineticEnergy, int particleCount,
	double targetTemperature, double* factor){
	if (factor == NULL)
		return false;
	*factor = 0;
	if (particleCount <= 0 || !isfinite(kineticEnergy) ||
		kineticEnergy <= 0 || !isfinite(targetTemperature) ||
		targetTemperature <= 0)
		return false;
	double ratio = targetTemperature/(
		kineticEnergy/(double)particleCount);
	if (!isfinite(ratio) || ratio <= 0)
		return false;
	double value = sqrt(ratio);
	if (!isfinite(value) || value <= 0)
		return false;
	*factor = value;
	return true;
}

void addNoise(){
	//double k = 2*M_PI*3/Lx;
	if (noise == 2 || noise == 4){
		physicalQ();
	}
	double velocityRescale = 1;
	if (noise == 2 && !velocityRescaleFactor(E, N, T, &velocityRescale)){
		#if G
		simulationSetFault(
			"Global velocity rescaling cannot continue at t=%.9g: kinetic "
			"energy and target temperature must both be finite and strictly "
			"positive (E=%.9g, T=%.9g). Change the setup and restart.",
			t, E, T);
		return;
		#else
		fprintf(stderr,
			"Global velocity rescaling failed at t=%.17g: kinetic energy and "
			"target temperature must both be finite and strictly positive "
			"(E=%.17g, T=%.17g).\n", t, E, T);
		exit(EXIT_FAILURE);
		#endif
	}
	#if INTRUDER
	freeFlyIntruder(&intru);
	#endif
	beginTransportVelocityBatch();
	if (noise == 4){
		for (int i = 0; i < N; i++){
			particle* p = particles + i;
			freeFly(p);
			const velocitySnapshot old = captureVelocity(p);
			p->vx -= px/(N*p->m);
			p->vy -= py/(N*p->m);
			#if THREE_D
			p->vz -= pz/(N*p->m);
			#endif
			recordCapturedVelocity(p, old);
			p->coll++;
		}
	}
	else if (noise == 8){
		for (int i = 0; i < N; i++){
			particle* p = particles + i;
			p->coll++;

			freeFly(p);
			const velocitySnapshot old = captureVelocity(p);
			double u1 = genrand_real3();
			double u2 = genrand_real3();
			double a = sqrt(-2*log(u1));
			double b = 2*M_PI*u2;
			double std = sqrt(T/p->m);
			p->vx = std*a*cos(b);
			p->vy = std*a*sin(b);
			#if THREE_D
			double u3 = genrand_real3();
			double u4 = genrand_real3();
			p->vz = std*sqrt(-2*log(u3))*cos(2*M_PI*u4);
			#endif
			recordCapturedVelocity(p, old);
		}
	}
	else{
		for (int i = 0; i < N; i++){

			particle* p = particles + i; //(int)(genrand_int32()%N)
			p->coll++;
			//updates position
			freeFly(p);
			const velocitySnapshot old = captureVelocity(p);
			//adds the kick
			if (noise == 1){
				randomGaussian(p, T);
				#if POLYGON_PARTICLES && !THREE_D
				if (polygon_particles_active() &&
					!polygon_particles_apply_angular_langevin(
						i, T, gamm, dtnoise, euler)) {
					simulationReportFault(
						"The polygon Langevin thermostat could not produce a finite angular update. Use the analytical OU integrator and finite positive gamma/interval values, then restart.");
					return;
				}
				#endif
			}
			else{
				p->vx *= velocityRescale;
				p->vy *= velocityRescale;
				#if THREE_D
				p->vz *= velocityRescale;
				#endif
			}
			recordCapturedVelocity(p, old);
			//p->vx += 0.03*sin(k*p->x - 0.1*t)*dtnoise;
			//p->vy += 0.3*sin(k*p->y - 0.1*t)*dtnoise;

			//Calculate the new collisions and the new cell crossings.
		}
	}
	for (int j = 0; j < N; j++){
		removeEventFromQueue(EVENT_NODE(j));
		crossingEvent(j);

		removeEventFromQueue(EVENT_NODE(N + j));
		collisionEvent(j);
	}
	#if INTRUDER
	if (freeFly == &freeFlyNormal){
		removeEventFromQueue(EVENT_NODE(2*N + 8));
		collisionEventIntruder();
	}
	#endif
	addEventNoise(t + dtnoise);
}


/* -------------------------------/
/								  /
/     Miscellaneous functions  	  /
/								  /
/--------------------------------*/


void freeFlyNormal(particle* p){
	const particlePrediction prediction = predictFreeFlyNormal(p, t);
	p->t = t;
	p->x = prediction.x;
	p->y = prediction.y;
	p->vx = prediction.vx;
	p->vy = prediction.vy;
	#if THREE_D
	p->z = prediction.z;
	p->vz = prediction.vz;
	#endif

	if (!addWallx && !addCircularWall)
		PBCpostX(&(p->x));
	if (!addWally && !addCircularWall){
		/* A lower Lees crossing deliberately represents the new image at
		 * y=Ly with vy<0.  Preserve that endpoint if another same-time event
		 * synchronizes the particle before it has moved into the last cell. */
		if (!(shearDynamicsActive() && p->y == Ly && p->vy < 0))
			PBCpostY(&(p->y));
	}
	#if THREE_D
	if (!addWallz && !addCircularWall)
		PBCpostZ(&(p->z));
	#endif
}

void freeFlyGrow(particle* p){
	const particlePrediction prediction = predictFreeFlyGrow(p, t);
	p->t = t;
	p->x = prediction.x;
	p->y = prediction.y;
	p->rad = prediction.rad;
	#if THREE_D
	p->z = prediction.z;
	#endif

	if (!addWallx && !addCircularWall)
		PBCpostX(&(p->x));
	if (!addWally && !addCircularWall)
		PBCpostY(&(p->y));
	#if THREE_D
	if (!addWallz && !addCircularWall)
		PBCpostZ(&(p->z));
	#endif
}




void saveTXT(){
	if (fichier == NULL){
		return;
	}
	int num_particles = N;
	#if INTRUDER
	num_particles = N + countFakeParticlesForIntruder(&intru, displayRad);
	#endif
	if (reduce){
		#if THREE_D
		fprintf(fichier, "ITEM: TIMESTEP\n%lf\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n0 %lf\n0 %lf\n0 %lf\nITEM: ATOMS id x y z vx vy vz\n", t, num_particles, Lx, Ly, Lz);
		for(int i = 0; i < N; i++){
			double X = particles[i].x;
			double Y = particles[i].y;
			double Z = particles[i].z;
			if (unwrap){
				X += particles[i].crossX*Lx;
				Y += particles[i].crossY*Ly;
				Z += particles[i].crossZ*Lz;
			}
			double vx = particles[i].vx;
			if (addShear)
				vx += strainRate*(particles[i].y - halfLy);
			fprintf(fichier, "%d %.2lf %.2lf %.2lf %.6lf %.6lf %.6lf\n",
				i, X, Y, Z, vx, particles[i].vy, particles[i].vz);
		}
		#else
		fprintf(fichier, "ITEM: TIMESTEP\n%lf\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n0 %lf\n0 %lf\n0 0\nITEM: ATOMS id x y vx vy\n", t, num_particles, Lx, Ly);
		for(int i = 0; i < N; i++){
			double X = particles[i].x;
			double Y = particles[i].y;
			if (unwrap){
				X += particles[i].crossX*Lx;
				Y += particles[i].crossY*Ly;
			}
			double vx = particles[i].vx;
			if (addShear){
				vx += strainRate*(particles[i].y - halfLy);
			}
			
			fprintf(fichier, "%d %.2lf %.2lf %.2lf %.2lf\n", i, X, Y, vx, particles[i].vy);
		}
		#endif
	}
	else{
		boop_data* boop = NULL;
		double* area = NULL;
		int* particleCluster = NULL;
		if (clusterThermo){
			findClusters(particles, N, clusterCutoff, cellList, Nxcells, Nycells);
			particleCluster = malloc(N * sizeof(int));
			for (int i = 0; i < numClusters; i++) {
				for (int j = 0; j < clusters[i].size; j++) {
					particleCluster[clusters[i].particles[j]->num] = i;
				}
			}
		}
		#if THREE_D
		fprintf(fichier, "ITEM: TIMESTEP\n%.17g\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n0 %.17g\n0 %.17g\n0 %.17g\nITEM: ATOMS id type x y z vx vy vz radius m coll\n", t, num_particles, Lx, Ly, Lz);

		#else
		fprintf(fichier, "ITEM: TIMESTEP\n%.17g\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n0 %.17g\n0 %.17g\n0 0\nITEM: ATOMS id type x y vx vy radius m coll", t, num_particles, Lx, Ly);
		if (boopThermo == 1){
			fprintf(fichier, " q5 q6 q7 argq6 neighbors");
			boop = computeBOOPVoronoi(particles, N, Lx, Ly);
		}
		else if (boopThermo == 2){
			fprintf(fichier, " q5 q6 q7 argq6 neighbors");
			boop = computeBOOPCutoff(particles, N, 2.5, cellList, Nxcells);
		}
		if (areaThermo){
			fprintf(fichier, " packingFraction");
			area = get_particle_voronoi_area(particles, N, Lx, Ly);
		}
		
		#if TANGENTIAL
		fprintf(fichier, " w");
		#endif
		fprintf(fichier, "\n");
		#endif
		double shift = 0;
		if (snapshotCritical){
			shift = getCOMx(particles, N, Lx);
		}
		for(int i = 0; i < N; i++){
			double X = particles[i].x;
			int type = particles[i].type;
			if (snapshotCritical){
				X = particles[i].x + (Lx / 4 - shift);
				if (X < 0) {
           			X += Lx;
				}
				else if (X > Lx){
					X -= Lx;
				}
				double dfrac = Lx*fracCritical/8.;
				double left = Lx/4;
				double right = 3*Lx/4;
				
				int firstx = ( (X > (left - dfrac)) && (X < (left + dfrac)) );
				int seconx = ( (X > (right - dfrac)) && (X < (right + dfrac)) );
				int firsty = (particles[i].y < (Ly/2));
				int secony = (particles[i].y > (Ly/2));
				#if THREE_D
				int firstz = (particles[i].z < (Lz/2));
				int seconz = (particles[i].z > (Lz/2));
				type = (firstx && firsty && firstz) ? 2 :
					   (firstx && secony && firstz) ? 3 :
					   (firstx && firsty && seconz) ? 4 :
					   (firstx && secony && seconz) ? 5 :
					   (seconx && firsty && firstz) ? 6 :
					   (seconx && secony && firstz) ? 7 :
					   (seconx && firsty && seconz) ? 8 :
					   (seconx && secony && seconz) ? 9 : type;
				#else
				type = (firstx && firsty) ? 2 :
					   (firstx && secony) ? 3 :
					   (seconx && firsty) ? 4 :
					   (seconx && secony) ? 5 : type;
				#endif
			}
			else if (clusterThermo){
				type = particleCluster[i];
			}

			#if THREE_D
			double Y = particles[i].y;
			double Z = particles[i].z;
			if (unwrap){
				X += particles[i].crossX*Lx;
				Y += particles[i].crossY*Ly;
				Z += particles[i].crossZ*Lz;
			}
			fprintf(fichier, "%d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %d\n", i, particles[i].type, X, Y, Z, particles[i].vx, particles[i].vy, particles[i].vz, particles[i].rad, particles[i].m, type);
			#else
			double Y = particles[i].y;
			if (unwrap){
				Y += particles[i].crossY*Ly;
				X += particles[i].crossX*Lx;
			}
			fprintf(fichier, "%d %d %.17g %.17g %.17g %.17g %.17g %.17g %d", i, particles[i].type, X , Y, particles[i].vx, particles[i].vy, particles[i].rad, particles[i].m, type);
			#if TANGENTIAL
			fprintf(fichier, " %lf", particles[i].omega);
			#endif
			if (boopThermo){
				fprintf(fichier, " %.2lf %.2lf %.2lf %.2lf %d", boop[i].q5, boop[i].q6, boop[i].q7, boop[i].q6_arg, boop[i].neighbors);
			}
			if (areaThermo){
				fprintf(fichier, " %lf", M_PI*particles[i].rad*particles[i].rad/area[i]);
			}
			fprintf(fichier, "\n");
			#endif
		}
		if (clusterThermo){
			freeClusters();
			free(particleCluster);
		}
		if (boopThermo){
			free(boop);
		}
		if (areaThermo){
			free(area);
		}
	}
	#if INTRUDER
		printIntruderShapeFromSmallParticles(fichier, &intru, displayRad, N);

	#endif
	fflush(fichier);
}


// Writes the x-resolved energy and heat-current profile into tempProfileFile.
// Each block contains: mean kinetic energy per particle, Jx, Jy[, Jz].
void saveTempProfile(){
	if (!tempProfileThermo){
		return;
	}
	if (tempProfileFile == NULL){
		return;
	}
	if (tempSlabWidth <= 0){
		tempSlabWidth = Lx/((double)tempProfileNSlabs);
	}

	outputSamples.temperatureProfileCount++;

	for (int i = 0; i < N; i++){
		freeFly(particles + i);
	}

	for (int i = 0; i < N; i++){
		particle* p = particles + i;
		double x = p->x;
		PBCpostX(&x);
		int slab = (int)(x/tempSlabWidth);
		if (slab < 0){
			slab = 0;
		}
		else if (slab >= tempProfileNSlabs){
			slab = tempProfileNSlabs - 1;
		}

		double vx = p->vx;
		if (addShear){
			vx += strainRate*(p->y - halfLy);
		}
		double ke = 0.5*p->m*(vx*vx + p->vy*p->vy);
		#if THREE_D
		ke += 0.5*p->m*(p->vz*p->vz);
		#endif

		energyProfile[slab] += ke;
		countTempProfile[slab] += 1;
		JxProfile[slab] += ke*vx;
		JyProfile[slab] += ke*p->vy;
		#if THREE_D
		JzProfile[slab] += ke*p->vz;
		#endif
	}

	if (outputSamples.temperatureProfileCount >= Nthermo){
		for (int i = 0; i < tempProfileNSlabs; i++){
			double e = 0;
			if (countTempProfile[i] > 0){
				e = energyProfile[i]/countTempProfile[i];
			}
			fprintf(tempProfileFile, "%.4e ", e);
			energyProfile[i] = 0;
			countTempProfile[i] = 0;
		}
		fprintf(tempProfileFile, "\n");

		double slabVolume = Ly*tempSlabWidth
			#if THREE_D
			*Lz
			#endif
		;
		for (int i = 0; i < tempProfileNSlabs; i++){
			fprintf(tempProfileFile, "%.4e ",
				JxProfile[i]/(slabVolume*outputSamples.temperatureProfileCount));
			JxProfile[i] = 0;
		}
		fprintf(tempProfileFile, "\n");

		for (int i = 0; i < tempProfileNSlabs; i++){
			fprintf(tempProfileFile, "%.4e ",
				JyProfile[i]/(slabVolume*outputSamples.temperatureProfileCount));
			JyProfile[i] = 0;
		}
		fprintf(tempProfileFile, "\n");

		#if THREE_D
		for (int i = 0; i < tempProfileNSlabs; i++){
			fprintf(tempProfileFile, "%.4e ",
				JzProfile[i]/(slabVolume*outputSamples.temperatureProfileCount));
			JzProfile[i] = 0;
		}
		fprintf(tempProfileFile, "\n");
		#endif
		fflush(tempProfileFile);
		outputSamples.temperatureProfileCount = 0;
	}
}
static double outerWallMeasure(void){
	double measure = 0.0;
	#if THREE_D
	if (addCircularWall)
		measure += 4*M_PI*halfLx*halfLx;
	if (addWallx)
		measure += 2*Ly*Lz;
	if (addWally)
		measure += 2*Lx*Lz;
	if (addWallz)
		measure += 2*Lx*Ly;
	#else
	if (addCircularWall)
		measure += M_PI*Lx;
	if (addWallx)
		measure += 2*Ly;
	if (addWally)
		measure += 2*Lx;
	#endif
	return measure;
}

void saveThermo(){
	physicalQ();
	if (transverseVelocityCorrelationIsEnabled()){
		updateTransverseVelocityCorrelation(particles, t);
	}
	double pressureX = 0;
	double pressureY = 0;
	double pressureXY = 0;
	double pressureYX = 0;

	#if INTRUDER
	freeFlyIntruder(&intru);
	#endif
	#if HEAT
	double Jx = 0;
	double Jy = 0;
	#if THREE_D
	double Jz = 0;
	double energySourceRate = 0;
	#endif
	#endif
	double energyX = 0, energyY = 0, energyXY = 0, energyYX = 0;
	#if THREE_D
	double pressureZ = 0;
	double pressureXZ = 0, pressureZX = 0, pressureYZ = 0, pressureZY = 0;
	double energyZ = 0, energyXZ = 0, energyZX = 0, energyYZ = 0, energyZY = 0;
	double vol = addCircularWall ? (4./3.)*M_PI*halfLx*halfLx*halfLx : Lx*Ly*Lz;
	#else
	double area = Lx*Ly;
	if (addCircularWall){
		area = M_PI*halfLx*halfLx; 
	}
	#endif
	double deltaTime = t - lastScreen;
	lastScreen = t;
	#if (STRESS || HEAT)
	advanceTransportAccumulators();
	#endif

	#if STRESS
	energyX = runningKineticX/deltaTime;
	energyY = runningKineticY/deltaTime;
	energyXY = runningKineticXY/deltaTime;
	energyYX = runningKineticYX/deltaTime;
	runningKineticX = 0;
	runningKineticY = 0;
	runningKineticXY = 0;
	runningKineticYX = 0;
	#if THREE_D
	energyZ = runningKineticZ/deltaTime;
	energyXZ = runningKineticXZ/deltaTime;
	energyZX = runningKineticZX/deltaTime;
	energyYZ = runningKineticYZ/deltaTime;
	energyZY = runningKineticZY/deltaTime;
	runningKineticZ = 0;
	runningKineticXZ = 0;
	runningKineticZX = 0;
	runningKineticYZ = 0;
	runningKineticZY = 0;
	#endif
	#else
	for (int i = 0; i < N; i++){
		double vx = particles[i].vx;
		double vy = particles[i].vy;
		if (shearDynamicsActive()){
			double currentY = particles[i].y +
				(t - particles[i].t)*particles[i].vy;
			/* A particle cannot pass its queued crossing event, but tolerate a
			 * one-box roundoff excursion without choosing the wrong shear image. */
			if (currentY < 0)
				currentY += Ly;
			else if (currentY > Ly)
				currentY -= Ly;
			vx = shearPeculiarVelocityX(vx, currentY);
		}
		energyX += particles[i].m*vx*vx;
		energyY += particles[i].m*vy*vy;
		energyXY += particles[i].m*vx*vy;
		energyYX += particles[i].m*vy*vx;
		#if THREE_D
		double vz = particles[i].vz;
		energyZ += particles[i].m*vz*vz;
		energyXZ += particles[i].m*vx*vz;
		energyZX += particles[i].m*vz*vx;
		energyYZ += particles[i].m*vy*vz;
		energyZY += particles[i].m*vz*vy;
		#endif
	}
	#endif
	#if HEAT
	double heatKineticX = runningHeatX/deltaTime;
	double heatKineticY = runningHeatY/deltaTime;
	double heatCollX = collHeatX/deltaTime;
	double heatCollY = collHeatY/deltaTime;
	#if THREE_D
	double heatKineticZ = runningHeatZ/deltaTime;
	double heatCollZ = collHeatZ/deltaTime;
	Jx = (heatKineticX + heatCollX)/vol;
	Jy = (heatKineticY + heatCollY)/vol;
	Jz = (heatKineticZ + heatCollZ)/vol;
	energySourceRate = collisionEnergySource/(vol*deltaTime);
	#else
	Jx = (heatKineticX + heatCollX)/area;
	Jy = (heatKineticY + heatCollY)/area;
	#endif
	runningHeatX = 0;
	runningHeatY = 0;
	collHeatX = 0;
	collHeatY = 0;
	#if THREE_D
	runningHeatZ = 0;
	collHeatZ = 0;
	collisionEnergySource = 0;
	#endif
	#endif
	#if THREE_D
	pressureX = -collTermX/(vol*deltaTime) + energyX/vol;
	pressureY = -collTermY/(vol*deltaTime) + energyY/vol;
	pressureZ = -collTermZ/(vol*deltaTime) + energyZ/vol;
	pressureXY = -collTermXY/(vol*deltaTime) + energyXY/vol;
	pressureYX = -collTermYX/(vol*deltaTime) + energyYX/vol;
	pressureXZ = -collTermXZ/(vol*deltaTime) + energyXZ/vol;
	pressureZX = -collTermZX/(vol*deltaTime) + energyZX/vol;
	pressureYZ = -collTermYZ/(vol*deltaTime) + energyYZ/vol;
	pressureZY = -collTermZY/(vol*deltaTime) + energyZY/vol;
	pressure = (pressureX + pressureY + pressureZ)/3;
	collTermZ = 0;
	collTermXZ = 0;
	collTermZX = 0;
	collTermYZ = 0;
	collTermZY = 0;
	#else
	pressureX = -1*collTermX/(area*deltaTime) + energyX/area;
	pressureY = -1*collTermY/(area*deltaTime) + energyY/area;
	pressureXY = -1*collTermXY/(area*deltaTime) + energyXY/area;
	pressureYX = -1*collTermYX/(area*deltaTime) + energyYX/area;
	pressure = (-1/(deltaTime))*(collTermX + collTermY)/(2*area) + E/area;
	#endif
	collTermX = 0;
	collTermY = 0;
	collTermXY = 0;
	collTermYX = 0;

	/* Wall impulses use the same sampling cadence as the virial tensor.  The
	 * old code retained `dp` until the Nthermo output row, then divided the
	 * multi-sample impulse by only the final sample interval.  Form one
	 * pressure value per interval and average those values just like the bulk
	 * pressure components above. */
	const double wallMeasure = outerWallMeasure();
	if (wallMeasure > 0.0 && deltaTime > 0.0){
		double outerImpulse = dp;
		if (addMidWall)
			outerImpulse -= dpMid;
		if (outerImpulse < 0.0 && outerImpulse > -1e-12*(fabs(dp) + 1.0))
			outerImpulse = 0.0;
		if (outerImpulse >= 0.0){
			outputSamples.wallPressureAverage +=
				outerImpulse/(deltaTime*wallMeasure);
			outputSamples.wallPressureCount++;
		}
	}
	if (addMidWall && Ly > 0.0 && deltaTime > 0.0){
		outputSamples.midWallPressureAverage += dpMid/(deltaTime*Ly);
		outputSamples.membranePressureDifferenceAverage +=
			-(dpRight - dpLeft)/(deltaTime*Ly);
		outputSamples.midWallPressureCount++;
	}
	/* `cumulativeWallImpulse` is intentionally not reset here: it is the
	 * independent monotonic signal used by the live-analysis window. */
	dp = 0.0;
	dpMid = 0.0;
	dpRight = 0.0;
	dpLeft = 0.0;

	if (thermoPressureField){
		pressure_field_step(t, deltaTime);
	}
	outputSamples.energyAverage += E/N;
	outputSamples.pressureAverage += pressure;
	outputSamples.pressureXAverage += pressureX;
	outputSamples.pressureYAverage += pressureY;
	outputSamples.pressureXYAverage += pressureXY;
	outputSamples.pressureYXAverage += pressureYX;
	if (profileThermo){
		for (int i = 0; i < N; i++){
			int bin = (int)(particles[i].y/(Ly/nProfileBins));
			if (bin == nProfileBins){
				bin = nProfileBins - 1;
			}
			double vx = particles[i].vx;
			outputSamples.velocityProfile[bin] += vx;
			outputSamples.velocityProfileCount[bin]++;
		}

	}
	#if THREE_D
	outputSamples.pressureZAverage += pressureZ;
	outputSamples.pressureXZAverage += pressureXZ;
	outputSamples.pressureZXAverage += pressureZX;
	outputSamples.pressureYZAverage += pressureYZ;
	outputSamples.pressureZYAverage += pressureZY;
	#endif
	#if HEAT
	outputSamples.heatXAverage += Jx;
	outputSamples.heatYAverage += Jy;
	#if THREE_D
	outputSamples.heatZAverage += Jz;
	outputSamples.energySourceRateAverage += energySourceRate;
	#endif
	#endif

	if (thermoSonine){
		double vx_av = 0;
		for (int i = 0; i < N; i++){
			double vx = particles[i].vx;
			if (addShear){
				vx += strainRate*(particles[i].y - halfLy);
			}
			vx_av += vx;
		}
		vx_av /= N;

		double v_2_av = 0;
		double v_4_av = 0;
		for (int i = 0; i < N; i++){
			double vx = particles[i].vx;
			if (addShear){
				vx += strainRate*(particles[i].y - halfLy);
			}
			double dvx = vx - vx_av;
			double dvx2 = dvx*dvx;
			v_2_av += dvx2;
			v_4_av += dvx2*dvx2;
		}
		v_2_av /= N;
		v_4_av /= N;

		outputSamples.sonineAverage += (v_2_av > 0.0) ?
			v_4_av/(3.0*v_2_av*v_2_av) - 1.0 : 0.0;
	}
	bool clusterResultsReady = false;
	int largestCluster = 0;
	if (t != 0 && clusterThermo &&
		(ther || dumpCluster || killCluster > 0)){
		findClusters(particles, N, clusterCutoff, cellList, Nxcells, Nycells);
		clusterResultsReady = true;
		largestCluster = numClusters > 0 ? clusters[0].size : 0;
		if (killCluster > 0 && largestCluster > killCluster){
			printf("\nNucleation detected!\nExiting...");
			exit(3);
		}
	}

	outputSamples.thermoCount++;
	if (outputSamples.thermoCount % Nthermo == 0){
		if (t != 0 && ther && thermo != NULL){
			fprintf(thermo, "%lf ", t);
			if (addCollisionNumberThermo){
				fprintf(thermo, "%ld ", ncol);
			}
			if (addEnergyThermo){
				fprintf(thermo, "%lf ",
					outputSamples.energyAverage/outputSamples.thermoCount);
			}
			#if THREE_D
			if (addPressureThermo){
				fprintf(thermo, "%lf %lf %lf %lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf ",
					outputSamples.pressureAverage/outputSamples.thermoCount,
					outputSamples.pressureXAverage/outputSamples.thermoCount,
					outputSamples.pressureYAverage/outputSamples.thermoCount,
					outputSamples.pressureZAverage/outputSamples.thermoCount,
					outputSamples.pressureXYAverage/outputSamples.thermoCount,
					outputSamples.pressureYXAverage/outputSamples.thermoCount,
					outputSamples.pressureXZAverage/outputSamples.thermoCount,
					outputSamples.pressureZXAverage/outputSamples.thermoCount,
					outputSamples.pressureYZAverage/outputSamples.thermoCount,
					outputSamples.pressureZYAverage/outputSamples.thermoCount);
			}
			#if HEAT
			fprintf(thermo, "%.10lf %.10lf %.10lf %.10lf ",
				outputSamples.heatXAverage/outputSamples.thermoCount,
				outputSamples.heatYAverage/outputSamples.thermoCount,
				outputSamples.heatZAverage/outputSamples.thermoCount,
				outputSamples.energySourceRateAverage/outputSamples.thermoCount);
			#endif
			#else
			if (addPressureThermo){
				fprintf(thermo, "%lf %lf %lf %.10lf %.10lf ",
					outputSamples.pressureAverage/outputSamples.thermoCount,
					outputSamples.pressureXAverage/outputSamples.thermoCount,
					outputSamples.pressureYAverage/outputSamples.thermoCount,
					outputSamples.pressureXYAverage/outputSamples.thermoCount,
					outputSamples.pressureYXAverage/outputSamples.thermoCount);
			}
			#if HEAT
			fprintf(thermo, "%.10lf %.10lf ",
				outputSamples.heatXAverage/outputSamples.thermoCount,
				outputSamples.heatYAverage/outputSamples.thermoCount);
			#endif
			#endif
			if (addCircularWall || addWallx || addWally || addWallz){
				fprintf(thermo, "%lf ", outputSamples.wallPressureCount > 0 ?
					outputSamples.wallPressureAverage/
					outputSamples.wallPressureCount : 0.0);
			}
			if (addMidWall){
				int Nleft = 0;
				int Nright = 0;
				int Nsolute = 0;
				double Tsolute = 0;
				double TsolventL = 0;
				double TsolventR = 0;
				for (int i = 0; i < N; i++){
					if (particles[i].type == 0){
						Nleft += (particles[i].x < halfLx) ? 1 : 0;
						Nright += (particles[i].x > halfLx) ? 1 : 0;
						double temp = 0.5*particles[i].m*(particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy);
						TsolventL += (particles[i].x < halfLx) ? temp : 0;
						TsolventR += (particles[i].x > halfLx) ? temp : 0;
					}
					else{
						Nsolute++;
						Tsolute += 0.5*particles[i].m*(particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy);
					}
				}
				double ratio = M_PI/(halfLx*Ly);
				fprintf(thermo, "%lf %lf %lf %lf %lf %lf %lf ",
					outputSamples.midWallPressureCount > 0 ?
						outputSamples.midWallPressureAverage/
						outputSamples.midWallPressureCount : 0.0,
					outputSamples.midWallPressureCount > 0 ?
						outputSamples.membranePressureDifferenceAverage/
						outputSamples.midWallPressureCount : 0.0,
					ratio*Nleft*sizeratio*sizeratio,
					ratio*Nright*sizeratio*sizeratio,
					TsolventL/Nleft, TsolventR/Nright, Tsolute/(Nsolute));
			}
			if (addWell){
				Ep = 0;
				for (int i = 0; i < N; i ++){
					particle* pi = particles + i;
					for (int k = 0; k < pi->numberOfParticlesInWell; k++){
						int j = pi->particlesInWell[k];
						if (j <= i || j < 0 || j >= N)
							continue;
						Ep += potentialPairStepEnergy(pi, particles + j);
					}
				}
				fprintf(thermo, "%lf ", Ep/N);
			}
			if (msd){
				fprintf(thermo, "%lf ", MSD/N);
			}
			if (addEnergyCenterOfMassThermo){
				double ex = 0;
				double ey = 0;
				#if THREE_D
				double ez = 0;
				#endif
				for (int i = 0; i < N; i++){
					ex += particles[i].m*particles[i].vx;
					ey += particles[i].m*particles[i].vy;
					#if THREE_D
					ez += particles[i].m*particles[i].vz;
					#endif
				}
				fprintf(thermo, "%lf ", 0.5*(ex*ex + ey*ey
					#if THREE_D
					+ ez*ez
					#endif
					)/N);
			}
			if (critical){
				int count[8];
				double energyCount[8];
				double dfrac = Lx*fracCritical/8.;
				countConditions(particles, N, Lx, Ly, count, energyCount, dfrac);
				#if THREE_D
				fprintf(thermo, "%d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf ", count[0], count[1], count[2], count[3], count[4], count[5], count[6], count[7], energyCount[0], energyCount[1], energyCount[2], energyCount[3], energyCount[4], energyCount[5], energyCount[6], energyCount[7]);
				#else
				fprintf(thermo, "%d %d %d %d %lf %lf %lf %lf ", count[0], count[1], count[2], count[3], energyCount[0], energyCount[1], energyCount[2], energyCount[3]);
				#endif
			}
			if (clusterThermo){
				fprintf(thermo, "%d ", largestCluster);
			}
			if (boopThermo){
				boop_data* boop = NULL;
				if (boopThermo == 1){
					boop = computeBOOPVoronoi(particles, N, Lx, Ly);
				}
				else if (boopThermo == 2){
					boop = computeBOOPCutoff(particles, N, 2.5, cellList, Nxcells);
				}
				
				double q6 = 0;
				for (int i = 0; i < N; i++){
					q6 += boop[i].q6;
				}
				fprintf(thermo, "%lf ", q6/N);
				free(boop);
			}
			#if INTRUDER
			double vx_body = intru.vx*cos(intru.angle) + intru.vy*sin(intru.angle);
			double vy_body = -intru.vx*sin(intru.angle) + intru.vy*cos(intru.angle);
			double invDeltaTime = (deltaTime > 0.0) ? (1.0 / deltaTime) : 0.0;
			double meanFx = intruderImpulseX * invDeltaTime;
			double meanFy = intruderImpulseY * invDeltaTime;
			double meanTorque = intruderTorqueZ * invDeltaTime;
			fprintf(thermo, "%lf %lf %lf %lf %lf %lf %lf %lf %lf ", intru.x, intru.y, vx_body, vy_body, intru.angle, intru.omega, meanFx, meanFy, meanTorque);
			intruderImpulseX = 0.0;
			intruderImpulseY = 0.0;
			intruderTorqueZ = 0.0;
			#endif
			if (thermoKurt){
				double ax1 = 0.0;
				double ax3 = 0.0;
				double vx2 = 0.0;

				for (int i = 0; i < N; i++){
					double vx = particles[i].vx;   // better: subtract mean flow first if nonzero
					ax1 += fabs(vx);
					ax3 += fabs(vx)*fabs(vx)*fabs(vx);
					vx2 += vx*vx;
				}

				ax1 /= N;
				ax3 /= N;
				vx2 /= N;

				double v_th = (vx2 > 0.0) ? sqrt(vx2) : 0.0;   // one-component thermal speed

				double K1 = (v_th > 0.0) ? ax1 / (2.0 * v_th) : 0.0;
				double K3 = (v_th > 0.0) ? ax3 / (2.0 * v_th * v_th * v_th) : 0.0;

				double ratio = (K1 > 0.0) ? K3 / (2.0 * K1) : 0.0;

				fprintf(thermo, "%lf %lf %lf ", K1, K3, ratio);
			}
			if (thermoSonine){
				fprintf(thermo, "%lf ", outputSamples.sonineAverage/
					outputSamples.thermoCount);
			}
			fprintf(thermo, "\n");
			fflush(thermo);
		}
		if (profileThermo && profileFile != NULL){
			for (int i = 0; i < nProfileBins; i++){
				int samples = outputSamples.velocityProfileCount[i];
				fprintf(profileFile, "%lf ", samples > 0 ?
					outputSamples.velocityProfile[i]/samples : 0.0);
			}
			fprintf(profileFile, "\n");
			fflush(profileFile);
		}
		outputSamples.thermoCount = 0;
		outputSamples.energyAverage = 0;
		outputSamples.pressureAverage = 0;
		outputSamples.pressureXAverage = 0;
		outputSamples.pressureYAverage = 0;
		outputSamples.pressureXYAverage = 0;
		outputSamples.pressureYXAverage = 0;
		outputSamples.wallPressureAverage = 0;
		outputSamples.wallPressureCount = 0;
		outputSamples.midWallPressureAverage = 0;
		outputSamples.membranePressureDifferenceAverage = 0;
		outputSamples.midWallPressureCount = 0;
		#if THREE_D
		outputSamples.pressureZAverage = 0;
		outputSamples.pressureXZAverage = 0;
		outputSamples.pressureZXAverage = 0;
		outputSamples.pressureYZAverage = 0;
		outputSamples.pressureZYAverage = 0;
		#endif
		#if HEAT
		outputSamples.heatXAverage = 0;
		outputSamples.heatYAverage = 0;
		#if THREE_D
		outputSamples.heatZAverage = 0;
		outputSamples.energySourceRateAverage = 0;
		#endif
		#endif
		outputSamples.sonineAverage = 0;
		if (profileThermo){
			memset(outputSamples.velocityProfile, 0,
				sizeof(outputSamples.velocityProfile));
			memset(outputSamples.velocityProfileCount, 0,
				sizeof(outputSamples.velocityProfileCount));
		}
		
	}

	if (t != 0){	
		if (liquidliquid && interfaceThermo && interface != NULL){
			computeInterfacesPos(particles, 0.1);
			fflush(interface);
		}
		else if (interfaceThermo && interface != NULL){
			saveDensityCoarse(particles);
			fflush(interface);
		}
		if (strucThermo && strucFile != NULL){
			//double start_time = omp_get_wtime();
			saveStructureFactor(particles, strucThermo);
			//double end_time = omp_get_wtime();
			//double time_spent = end_time - start_time;
			//printf("\n Time struc: %f seconds\n", time_spent);
			fflush(strucFile);
		}
		if (dumpCluster && clusterResultsReady && clusterFile != NULL){
			saveHistrogramCluster(clusterFile);
			fflush(clusterFile);
		}
		if (pcfBondOrderThermo){
			if (!save_pcf_boop(pcfName, pcfBondOrderName, particles, N, 2,
				fmin(Lx, Ly)/2, Lx, Ly, phi)){
				reportOutputRuntimeFailure(
					"Unable to write pair/bond-order correlation outputs");
				pcfBondOrderThermo = 0;
			}
		}
		else if (pcfThermo){
			if (!save_pcf(pcfName, particles, N, 0.1, fmin(Lx, Ly)/2,
				Lx, Ly)){
				reportOutputRuntimeFailure(
					"Unable to write pair-correlation output '%s'", pcfName);
				pcfThermo = 0;
			}
		}
		if (pcfg6Thermo){
			if (!save_pcf_g6(pcfg6Name, particles, N, 2,
				fmin(Lx, Ly)/2, Lx, Ly)){
				reportOutputRuntimeFailure(
					"Unable to write g6 correlation output '%s'", pcfg6Name);
				pcfg6Thermo = 0;
			}
		}
		if (osmosisThermo && osmosisFile != NULL){
			save_osmosis(particles, N);
		}
		if (velocityDistributionThermo && velocityFile != NULL){
			for (int i = 0; i < N; i++){
				#if THREE_D
				fprintf(velocityFile, "%.4lf %.4lf %.4lf\n",
					particles[i].vx, particles[i].vy, particles[i].vz);
				#else
				fprintf(velocityFile, "%.4lf %.4lf\n",
					particles[i].vx, particles[i].vy);
				#endif
			}
			fflush(velocityFile);
		}
	}
	if (clusterResultsReady){
		freeClusters();
	}
	if (bubbleStop){
		if (0){
			double Z = (1 + phi*phi/8)/pow(1 - phi, 2);
			double dZdphi = 2*Z/(1 - phi) + (phi/4)/pow(1 - phi, 2);
			double compressibility = 1/(Z + phi*dZdphi);
			double re = 0, im = 0;
			double qx = 4*M_PI/Lx;
			double qy = 4*M_PI/Ly;
			for (int n = 0; n < N; n++){
					particle* p = particles + n;
					double qr = qx*p->x + qy*p->y;

					re += cos(qr);
					im += sin(qr);
				}
			
			double structFactor = (re*re + im*im)/N;
			if (structFactor > 10*compressibility){
				outputSamples.bubbleCount++;
			}
			else{
				outputSamples.bubbleCount = 0;
			}
		}
		double energyTreshold = 2*delta*delta/(1 - res*res);
		if (damping == 1){
			double g_contact = (((1 + phi*phi/8 - phi*phi*phi*phi/10)/((1 - phi)*(1 - phi))) - 1)/(2*phi);
			double freq = 4*phi*sqrt(1/M_PI)*g_contact;
			double epsi = 4*gamm/(freq);
			energyTreshold = pow((epsi + sqrt(epsi*epsi + 4*(1-res*res))/(2*(1 - res*res))), 2);
		}
		if (E/N > 2*energyTreshold){
			outputSamples.bubbleCount++;
		}
		else{
			outputSamples.bubbleCount = 0;
		}
		if (outputSamples.bubbleCount >= 2 || (E/N > 4*energyTreshold)){
			printf("\nBubble detected! Stopping simulation...\n");
			if (thermo != NULL){
				fprintf(thermo, "-1\n");
			}
			exit(4);
		}
	}

}

double drand(double min, double max){
    return (genrand()*(max - min)) + min;
}

double sign(double x){
	return (x < 0) ? -1 : (x > 0);
}



inline double logTime(double time){
	if (damping != 0){
		double var = (1 - time*gamm);
		if (var < 0){
			return never;
		}
		else{
			return -log(var)/gamm;
		}
	}
	else
		return time;
}

void normalizePhysicalQ(){
	if (addCircularWall != 1){
		physicalQ();
		for (int i = 0; i < N; i++){
			particles[i].vx -= px/(N*particles[i].m);
			particles[i].vy -= py/(N*particles[i].m);
			#if THREE_D
			particles[i].vz -= pz/(N*particles[i].m);
			#endif
		}
	}
	else{
		/* Remove the rigid-body rotation of a confined sample.  In 3D this
		   requires the full inertia tensor, rather than only Lz. */
		#if THREE_D
		double inertia[3][3] = {{0}};
		double angularMomentum[3] = {0};
		for (int i = 0; i < N; i++){
			particle* p = particles + i;
			double rx = p->x - halfLx;
			double ry = p->y - halfLy;
			double rz = p->z - halfLz;
			double r2 = rx*rx + ry*ry + rz*rz;
			angularMomentum[0] += p->m*(ry*p->vz - rz*p->vy);
			angularMomentum[1] += p->m*(rz*p->vx - rx*p->vz);
			angularMomentum[2] += p->m*(rx*p->vy - ry*p->vx);
			inertia[0][0] += p->m*(r2 - rx*rx);
			inertia[1][1] += p->m*(r2 - ry*ry);
			inertia[2][2] += p->m*(r2 - rz*rz);
			inertia[0][1] -= p->m*rx*ry;
			inertia[0][2] -= p->m*rx*rz;
			inertia[1][2] -= p->m*ry*rz;
		}
		inertia[1][0] = inertia[0][1];
		inertia[2][0] = inertia[0][2];
		inertia[2][1] = inertia[1][2];
		double determinant =
			inertia[0][0]*(inertia[1][1]*inertia[2][2] - inertia[1][2]*inertia[2][1]) -
			inertia[0][1]*(inertia[1][0]*inertia[2][2] - inertia[1][2]*inertia[2][0]) +
			inertia[0][2]*(inertia[1][0]*inertia[2][1] - inertia[1][1]*inertia[2][0]);
		if (fabs(determinant) > 1e-20){
			double inverse[3][3];
			inverse[0][0] = (inertia[1][1]*inertia[2][2] - inertia[1][2]*inertia[2][1])/determinant;
			inverse[0][1] = (inertia[0][2]*inertia[2][1] - inertia[0][1]*inertia[2][2])/determinant;
			inverse[0][2] = (inertia[0][1]*inertia[1][2] - inertia[0][2]*inertia[1][1])/determinant;
			inverse[1][0] = (inertia[1][2]*inertia[2][0] - inertia[1][0]*inertia[2][2])/determinant;
			inverse[1][1] = (inertia[0][0]*inertia[2][2] - inertia[0][2]*inertia[2][0])/determinant;
			inverse[1][2] = (inertia[0][2]*inertia[1][0] - inertia[0][0]*inertia[1][2])/determinant;
			inverse[2][0] = (inertia[1][0]*inertia[2][1] - inertia[1][1]*inertia[2][0])/determinant;
			inverse[2][1] = (inertia[0][1]*inertia[2][0] - inertia[0][0]*inertia[2][1])/determinant;
			inverse[2][2] = (inertia[0][0]*inertia[1][1] - inertia[0][1]*inertia[1][0])/determinant;
			double omega[3] = {0};
			for (int row = 0; row < 3; row++)
				for (int column = 0; column < 3; column++)
					omega[row] += inverse[row][column]*angularMomentum[column];
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				double rx = p->x - halfLx;
				double ry = p->y - halfLy;
				double rz = p->z - halfLz;
				p->vx -= omega[1]*rz - omega[2]*ry;
				p->vy -= omega[2]*rx - omega[0]*rz;
				p->vz -= omega[0]*ry - omega[1]*rx;
			}
		}
		#else
		double angularMomentum = 0;
		double inertia = 0;
		for (int i = 0; i < N; i++){
			particle* p = particles + i;
			double rx = p->x - halfLx;
			double ry = p->y - halfLy;
			angularMomentum += p->m*(rx*p->vy - ry*p->vx);
			inertia += p->m*(rx*rx + ry*ry);
		}
		if (inertia > 0){
			double omega = angularMomentum/inertia;
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				double rx = p->x - halfLx;
				double ry = p->y - halfLy;
				p->vx += omega*ry;
				p->vy -= omega*rx;
			}
		}
		#endif
	}


	physicalQ();
	#if POLYGON_PARTICLES && !THREE_D
	if (polygon_particles_active()){
		/* Polygon spin is initialized from its own angular-temperature
		   setting.  Preserve the historical meaning of Einit as the target
		   translational energy instead of counting rotation in a scale factor
		   that is applied only to vx/vy. */
		const double translationalEnergy =
			E - polygon_particles_rotational_energy();
		if (translationalEnergy > 0 && Einit > 0){
			const double velocityScale =
				sqrt(translationalEnergy/N/Einit);
			for (int i = 0; i < N; i++){
				particles[i].vx /= velocityScale;
				particles[i].vy /= velocityScale;
			}
		}
		return;
	}
	#endif
	for (int i = 0; i < N; i++){
		particles[i].vx /= sqrt(E/N/Einit);
		particles[i].vy /= sqrt(E/N/Einit);
		#if THREE_D
		particles[i].vz /= sqrt(E/N/Einit);
		#endif
	}
}

void optimizeGrowConstant(){
	double baseGrow;
	#if THREE_D
	double criticalPhi = 0.6;
	#else
	double criticalPhi = 0.8;
	#endif
	if (addCircularWall){
		baseGrow = 0.03;
	}
	else{
		baseGrow = 0.1;
	}
	#if POLYGON_PARTICLES && !THREE_D
	if (polygon_particles_requested()){
		/* The legacy density slowdown was tuned for disks close to their
		 * envelope jamming point. Actual polygons need enough expansion-driven
		 * contact motion to rotate and unjam; slowing by 8x instead creates
		 * resolution-level chatter. */
		vr = baseGrow;
		return;
	}
	#endif
	if (polydispersity != 0){
		if ((sizeratio > 0.5) && (fractionSmallN < 0.9) && (phi > 0.82)){
			baseGrow *= 0.1;
		}
	}
	if (phi < criticalPhi){
		vr = baseGrow;
	}
	else{
		vr = baseGrow*pow(criticalPhi/phi, 30);
	}
	// Try to homogenize the system!
	if (clusterThermo || phi > 0.7){
		vr /= 8;
	}
	#if INTRUDER
	if (bathRad < 0.8){
		vr = 0.1;
	}
	#endif
}


void randomGaussian(particle* p, double T){
	double u1 = genrand_real3();
	double u2 = genrand_real3();
	double a = sqrt(-2*log(u1));
	double b = 2*M_PI*u2;
	#if THREE_D
	double u3 = genrand_real3();
	double u4 = genrand_real3();
	double gaussianZ = sqrt(-2*log(u3))*cos(2*M_PI*u4);
	#endif

	if (damping == 0){
		if (euler){
			double std = sqrt(2*gamm*T*dtnoise)/p->m;
			p->vx += std*a*cos(b) - sign(p->vx)*gamm*dtnoise/p->m;
			p->vy += std*a*sin(b) - sign(p->vy)*gamm*dtnoise/p->m;
			#if THREE_D
			p->vz += std*gaussianZ - sign(p->vz)*gamm*dtnoise/p->m;
			#endif
		}
		else{
			double c = exp(-gamm*dtnoise);
			double std = sqrt(T*(1 - c*c)/p->m);
			p->vx = std*a*cos(b) + (p->vx)*c;
			p->vy = std*a*sin(b) + (p->vy)*c;
			#if THREE_D
			p->vz = std*gaussianZ + p->vz*c;
			#endif
			//p->vx = p->vx - gamm*p->vx*dtnoise/p->m - gamm*pow(p->vx, 3)*dtnoise/p->m;
			//p->vy = p->vy - gamm*p->vy*dtnoise/p->m - gamm*pow(p->vy, 3)*dtnoise/p->m;
		}
	}
	else{
		double std = sqrt(2*gamm*T*dtnoise/p->m);
		p->vx += std*a*cos(b);
		p->vy += std*a*sin(b);
		#if THREE_D
		p->vz += std*gaussianZ;
		#endif
	}
}

void shearCorrection(double* x, double* vx, int sign, particle* p){
	#if (STRESS || HEAT)
	advanceTransportAccumulators();
	#endif
	(void)p;
	*x += sign*fmod(runtimeShearStrain()*Ly, Lx);
	*vx += sign*strainRate*Ly;
	/* vx and y change image together at the seam, so the peculiar velocity
	 * vx + strainRate*(y-halfLy) is continuous.  No kinetic/heat impulse is
	 * recorded for this coordinate remapping. */
	PBCpostX(x);
}

void correctDistances(const particle* p1, const particle* p2,
	double p2x, double p2vx, double* dx, double* dvx){
	if (p1->cell[1] == 0 && p2->cell[1] == Nycells - 1){
		double x = p2x + fmod(runtimeShearStrain()*Ly, Lx);
		double vx = p2vx + strainRate*Ly;
		PBCpostX(&x);
		*dx = x - p1->x;
		*dvx = vx - p1->vx;
		}
	if (p1->cell[1] == Nycells - 1 && p2->cell[1] == 0){
		double x = p2x - fmod(runtimeShearStrain()*Ly, Lx);
		double vx = p2vx - strainRate*Ly;
		PBCpostX(&x);
		*dx = x - p1->x;
		*dvx = vx - p1->vx;
	}
}

#if THREE_D
void PBC(double* dx, double* dy, double* dz){
	if ((!addWallx) && (!addCircularWall)){
		if (*dx >= halfLx){
			*dx -= Lx;
		}
		else if (*dx < -halfLx){
			*dx += Lx;
		}			
	}
	if ((!addWally) && (!addCircularWall)){
		if (*dy >= halfLy){
			*dy -= Ly;
		}
		else if (*dy < -halfLy){
			*dy += Ly;
		}
	}
	
	if (!addWallz && !addCircularWall){
		if (*dz >= halfLz){
			*dz -= Lz;
		}
		else if (*dz < -halfLz){
			*dz += Lz;
		}
	}
}
#else
void PBC(double* dx, double* dy){
	if ((!addWallx) && (!addCircularWall)){
		if (*dx >= halfLx){
			*dx -= Lx;
		}
		else if (*dx < -halfLx){
			*dx += Lx;
		}			
	}
	if ((!addWally) && (!addCircularWall)){
		if (*dy >= halfLy){
			*dy -= Ly;
		}
		else if (*dy < -halfLy){
			*dy += Ly;
		}
	}
}
#endif

int cmpDouble(const void * a, const void * b){
    double A = *(double*) a;
   	double B = *(double*) b;
    return (int)sign(B - A);
}

double PBCinsideCellX(double dx){
	if (dx >= halfLx)
		return dx - Lx;
	else if (dx < -halfLx)
		return dx + Lx;
	return dx;
}

double PBCinsideCellY(double dy){
	if (dy >= halfLy)
		return dy - Ly;
	else if (dy < -halfLy)
		return dy + Ly;
	return dy;
}

#if THREE_D
double PBCinsideCellZ(double dz){
	if (dz >= halfLz)
		return dz - Lz;
	else if (dz < -halfLz)
		return dz + Lz;
	return dz;
}
#endif

void PBCpostX(double* val){
	if (*val < 0)
		*val += Lx;
	else if (*val >= Lx)
		*val -= Lx;
}
void PBCpostY(double* val){
	if (*val < 0)
		*val += Ly;
	else if (*val >= Ly)
		*val -= Ly;
}
void PBCpostZ(double* val){
	if (*val < 0)
		*val += Lz;
	else if (*val >= Lz)
		*val -= Lz;
}


void physicalQ(){
	px = 0;
	py = 0;
	#if THREE_D
	pz = 0;
	#endif
	E = 0;
	active = 0;
	for (int i = 0; i < N; i++){
		px += particles[i].vx*particles[i].m;
		py += particles[i].vy*particles[i].m;
		#if THREE_D
		pz += particles[i].vz*particles[i].m;
		double vx = particles[i].vx;
		if (shearDynamicsActive()){
			double currentY = particles[i].y +
				(t - particles[i].t)*particles[i].vy;
			if (currentY < 0)
				currentY += Ly;
			else if (currentY > Ly)
				currentY -= Ly;
			vx = shearPeculiarVelocityX(vx, currentY);
		}
		E += 0.5*particles[i].m*(vx*vx + particles[i].vy*particles[i].vy +
			particles[i].vz*particles[i].vz);
		if (vx*vx + particles[i].vy*particles[i].vy +
			particles[i].vz*particles[i].vz > 1e-6){
		#else
		double vx = particles[i].vx;
		if (shearDynamicsActive()) {
			double currentY = particles[i].y +
				(t - particles[i].t)*particles[i].vy;
			if (currentY < 0)
				currentY += Ly;
			else if (currentY > Ly)
				currentY -= Ly;
			vx = shearPeculiarVelocityX(vx, currentY);
		}
		E += 0.5 * particles[i].m * (vx * vx + particles[i].vy * particles[i].vy);
		#if TANGENTIAL
		E += 0.5*particles[i].J*pow(particles[i].omega, 2);
		#endif
		if (vx*vx + particles[i].vy*particles[i].vy > 1e-6){
		#endif
			active += 1.;
		}
	}
	#if POLYGON_PARTICLES && !THREE_D
	E += polygon_particles_rotational_energy();
	#endif
	active /= N;
	if (msd != 0 && posxInitial != NULL && posyInitial != NULL &&
		#if THREE_D
		poszInitial != NULL &&
		#endif
		t >= firstThermo){
		MSD = 0;
		double COMt_x = 0;
		double COMt_y = 0;
		#if THREE_D
		double COMt_z = 0;
		#endif
		for (int i = 0; i < N; i++){
			COMt_x += particles[i].x + particles[i].crossX*Lx;
			COMt_y += particles[i].y + particles[i].crossY*Ly;
			#if THREE_D
			COMt_z += particles[i].z + particles[i].crossZ*Lz;
			#endif
		}
		COMt_x = COMt_x/N;
		COMt_y = COMt_y/N;
		#if THREE_D
		COMt_z = COMt_z/N;
		#endif
		for (int i = 0; i < N; i++){
			MSD += pow(particles[i].x + particles[i].crossX*Lx -
				posxInitial[i] - (COMt_x - COM_x), 2) +
				pow(particles[i].y + particles[i].crossY*Ly -
				posyInitial[i] - (COMt_y - COM_y), 2)
				#if THREE_D
				+ pow(particles[i].z + particles[i].crossZ*Lz -
					poszInitial[i] - (COMt_z - COM_z), 2)
				#endif
				;
		}
	}
	
}

static bool appendOutputName(char* destination, size_t capacity,
	const char* format, ...){
	size_t used = strlen(destination);
	if (used >= capacity - 1){
		return false;
	}
	va_list arguments;
	va_start(arguments, format);
	int written = vsnprintf(destination + used, capacity - used,
		format, arguments);
	va_end(arguments);
	return written >= 0 && (size_t)written < capacity - used;
}

static void makeVersionedOutputName(char* destination, size_t capacity,
	const char* base, int outputVersion, const char* suffix){
	/* customName() limits base to 300 bytes, leaving ample room for a signed
	 * integer version and the longest output suffix. */
	const int baseLimit = capacity > 48 ? (int)(capacity - 48) : 0;
	snprintf(destination, capacity, "%.*s%d%.*s", baseLimit, base,
		outputVersion, 24, suffix);
}

static bool buildOutputBaseName(char* destination, size_t capacity){
	#if defined(_WIN32)
	_mkdir("dump/");
	#else 
    mkdir("dump/", 0777);
	#endif

	int initialNameLength = snprintf(destination, capacity,
		"dump/N_%dres_%.3lfphi_%.6lfq_%.3lfrat_%.3lfLx_%.3lfLy_%.3lf",
		N, res, phi, (double)Nsmall/N, sizeratio, Lx, Ly);
	bool completeName = initialNameLength >= 0 &&
		initialNameLength < (int)capacity;
	#if THREE_D
		completeName = appendOutputName(destination, capacity,
			"Lz_%.3lf", Lz) && completeName;
	#endif	
	if (damping || noise){
		completeName = appendOutputName(destination, capacity,
			"gamma_%.8lfT_%.3lfdt_%.3lf",
			gamm, T, dtnoise) && completeName;
	};
	if (addDelta || addDeltaTangent || addActiveDelta){
		completeName = appendOutputName(destination, capacity,
			"delta_%.8lf", delta) && completeName;
	}
	if (addDoubleDelta || addDoubleDeltaTangent){
		completeName = appendOutputName(destination, capacity,
			"deltam_%.4lfdeltaM_%.4lfts_%.4lf",
			deltam, deltaM, ts) && completeName;
	}

	if (addEnergy){
		completeName = appendOutputName(destination, capacity,
			"deltaE_%.3lfDeltaE_%.3lfts_%.3lfbeta_%.3lf",
			additionalEnergy, deltaE, taur, beta) && completeName;
	}
	if (addMidWall){
		completeName = appendOutputName(destination, capacity,
			"phiS_%.6lfphis_%.6lf", phil, phig) &&
			completeName;
	}
	if (addShear){
		completeName = appendOutputName(destination, capacity,
			"shear_%.6lf", strainRate) && completeName;
	}
	#if INTRUDER
		completeName = appendOutputName(destination, capacity,
			"Im_%lfM_%.5lfresI_%.3lfdeltaI_%.3lf",
			intru.Im, intru.M, res_intruder, delta_intruder) && completeName;
		if (shape == GEN_WHEEL){
			completeName = appendOutputName(destination, capacity,
				"Rinside_%lfRoutside_%lfNspikes_%dChi_%lf",
				paramIntruder[0], paramIntruder[1], (int)paramIntruder[2],
				paramIntruder[3]) && completeName;
		}
		if (shape == POLYGON){
			completeName = appendOutputName(destination, capacity,
				"Nside_%dR_%lf",
				(int)paramIntruder[1], paramIntruder[0]) && completeName;
		}
	#endif

	completeName = appendOutputName(destination, capacity, "v_") && completeName;
	/* Reserve suffix space even when a very descriptive parameter set was
	 * selected. This also prevents different version numbers from truncating to
	 * the same path in the existence loop below. */
	if (strlen(destination) > 300){
		destination[300] = '\0';
		completeName = false;
	}
	if (!completeName){
		fprintf(stderr, "Output filename was shortened to fit safely\n");
	}
	return destination[0] != '\0';
}

static bool versionedOutputExists(const char* baseFileName,
	int candidateVersion){
	static const char* suffixes[] = {
		".dump", ".thermo", ".pressure_field", ".interface", ".osmosis",
		".struc", ".tvacf", ".cluster", ".pcfBO", ".pcf", ".pcfg6",
		".vdist", ".profile", ".temp", ".intruderKM"
	};
	char candidate[350];
	for (size_t index = 0; index < sizeof(suffixes)/sizeof(suffixes[0]);
		index++){
		makeVersionedOutputName(candidate, sizeof(candidate), baseFileName,
			candidateVersion, suffixes[index]);
		if (access(candidate, F_OK) == 0){
			return true;
		}
	}
	return false;
}

static int availableOutputVersion(const char* baseFileName){
	int candidateVersion = outputVersion;
	while (versionedOutputExists(baseFileName, candidateVersion)){
		if (candidateVersion == INT_MAX){
			break;
		}
		candidateVersion++;
	}
	return candidateVersion;
}

bool customSnapshotName(char* destination, size_t capacity){
	char baseFileName[350] = {0};
	if (destination == NULL || capacity == 0 ||
		!buildOutputBaseName(baseFileName, sizeof(baseFileName))){
		return false;
	}
	makeVersionedOutputName(destination, capacity, baseFileName,
		availableOutputVersion(baseFileName), ".dump");
	return destination[0] != '\0';
}

void customName(){
	char baseFileName[350] = {0};
	if (!buildOutputBaseName(baseFileName, sizeof(baseFileName))){
		fileName[0] = '\0';
		return;
	}
	const int resolvedVersion = availableOutputVersion(baseFileName);

	makeVersionedOutputName(fileName, sizeof(fileName), baseFileName,
		resolvedVersion, ".dump");

	makeVersionedOutputName(thermoName, sizeof(thermoName), baseFileName,
		resolvedVersion, ".thermo");
	if (thermoPressureField){
		makeVersionedOutputName(pressureFieldName, sizeof(pressureFieldName),
			baseFileName, resolvedVersion, ".pressure_field");
	}


	#if INTRUDER
		makeVersionedOutputName(intruderKMName, sizeof(intruderKMName),
			baseFileName, resolvedVersion, ".intruderKM");
	#endif


	if (interfaceThermo) {
		makeVersionedOutputName(interfaceName, sizeof(interfaceName), baseFileName,
			resolvedVersion, ".interface");
	}
	if (osmosisThermo) {
		makeVersionedOutputName(osmosisName, sizeof(osmosisName), baseFileName,
			resolvedVersion, ".osmosis");
	}

	if (strucThermo) {
		makeVersionedOutputName(strucName, sizeof(strucName), baseFileName,
			resolvedVersion, ".struc");
	}
	if (transverseVelocityCorrelationIsEnabled()) {
		makeVersionedOutputName(transverseCorrelationName,
			sizeof(transverseCorrelationName), baseFileName, resolvedVersion,
			".tvacf");
	}

	if (dumpCluster && clusterThermo) {
		makeVersionedOutputName(clusterName, sizeof(clusterName), baseFileName,
			resolvedVersion, ".cluster");
	}

	if (pcfBondOrderThermo) {
		makeVersionedOutputName(pcfBondOrderName, sizeof(pcfBondOrderName),
			baseFileName, resolvedVersion, ".pcfBO");
		makeVersionedOutputName(pcfName, sizeof(pcfName), baseFileName,
			resolvedVersion,
			".pcf");
	}
	else if (pcfThermo){
		makeVersionedOutputName(pcfName, sizeof(pcfName), baseFileName,
			resolvedVersion, ".pcf");
	}
	if (pcfg6Thermo){
		makeVersionedOutputName(pcfg6Name, sizeof(pcfg6Name), baseFileName,
			resolvedVersion, ".pcfg6");
	}
	if (velocityDistributionThermo){
		makeVersionedOutputName(velocityDistributionName,
			sizeof(velocityDistributionName), baseFileName, resolvedVersion,
			".vdist");
	}
	if (profileThermo){
		makeVersionedOutputName(profileName, sizeof(profileName), baseFileName,
			resolvedVersion, ".profile");
	}
	if (tempProfileThermo){
		makeVersionedOutputName(tempProfileName, sizeof(tempProfileName),
			baseFileName, resolvedVersion, ".temp");
	}
}

double resCoeff(double v){
	double temp = ao*exp(-v/vo);
	if (temp < 0.01)
		return 1;
	return temp;

}

void format(){
	printf("BOUNDARY CONDITIONS: ");
	CGREEN;
	if (addCircularWall){
		#if THREE_D
		printf("SPHERICAL WALL");
		#else
		printf("CIRCULAR WALL");
		#endif
	}
	else if ((addWallx) || (addWally) || (addWallz)){
		if (resW == 1)
			printf("HARD WALLS");
		else
			printf("DISSIPATIVE WALLS");
		}
	else
		printf("PERIODICS");
	CWHITE;
	printf("\nTHERMOSTAT: ");
	CGREEN;
	if (noise)
		printf("ON");
	else
		printf("OFF");
	CWHITE;
	printf("\nDISSIPATION: ");
	CGREEN;
	if ((res < 1) || addExpo || damping)
		printf("ON");
	else
		printf("OFF");
	CWHITE;
	printf("\nCOLLISION ENERGY INPUT: ");
	CGREEN;
	if ((addDelta || addDeltaTangent) && (addExpo))
		printf("Δ AND e^");
	else if (addExpo)
		printf("e^");
	else if (addDeltaTangent)
		printf("TANGENTIAL Δ (u_z x n)");
	else if (addDelta)
		printf("NORMAL Δ");
	else
		printf("NONE");
	CWHITE;
	printf("\nBOX SHAPE: ");
	CGREEN;
	if (addCircularWall){
		#if THREE_D
		printf("SPHERE");
		#else
		printf("CIRCLE");
		#endif
	}
	#if THREE_D
	else if ((Lx == Ly) && (Ly == Lz))
		printf("CUBE");
	else
		printf("CUBOID");
	#else
	else if (Lx == Ly)
		printf("SQUARE");
	else
		printf("RECTANGULAR");
	#endif
	CWHITE;
	printf("\n┌───────────────────────────────────────────────────────────────────────────────────────────────────────┐\n│\t N = \033[1;32m%d\033[0;37m Δt = \033[1;32m%.3lf\033[0;37m α = \033[1;32m%.3lf\033[0;37m γ = \033[1;32m%.3lf\033[0;37m T = \033[1;32m%.3lf\033[0;37m Φ = \033[1;32m%.3lf\033[0;37m rat = \033[1;32m%.3lf\033[0;37m L = \033[1;32m%.3lf\033[0;37m\t\t│ \n├────────────────────────────────────────────────┬──────────────────────────────────────────────────────┤\n", N, dtnoise, res, gamm, T, phi, sizeratio, Lx);
}

void printClose(){
	printf("\n");
	printf("└────────────────────────────────────────────────┴──────────────────────────────────────────────────────┘\n");
	if (bubbleStop){
		double energyTreshold = 2*delta*delta/(1 - res*res);
		if (damping == 1){
			double g_contact = (((1 + phi*phi/8 - phi*phi*phi*phi/10)/((1 - phi)*(1 - phi))) - 1)/(2*phi);
			double freq = 4*phi*sqrt(1/M_PI)*g_contact;
			double epsi = 4*gamm/(freq);
			energyTreshold = pow((epsi + sqrt(epsi*epsi + 4*(1-res*res))/(2*(1 - res*res))), 2);
		}
		if (E/N > 2*energyTreshold){
			printf("\nBubble detected! Stopping simulation...\n");
			if (thermo != NULL){
				fprintf(thermo, "-1\n");
			}
			exit(4);
		}
	}
}

void runningCheck(){
	double stop = 0;
	if (res > 1){
		printf("ERROR:\033[1;31m Coefficient of restitution greater than 1!\033[0m\n");
		stop = 1;
	}
	if (unforcedDissipationActive()){
		if (G){
			printf("WARNING:\033[0;31m Dissipative system without energy input; it will cool freely.\033[0m\n");
		}
		else{
			printf("WARNING:\033[0;31m Dissipative system without energy input; "
				"the headless run will stop below total kinetic energy %.1e.\033[0m\n",
				headlessKineticEnergyCutoff);
		}
	}

	if (field > 0){
		printf("ERROR:\033[1;31m Field must be negative!\033[0m\n");
		stop = 1;
	}
	if ((addDelta) && (res == 1)){
		printf("ERROR:\033[1;31m Delta model used without dissipation!\033[0m\n");
		stop = 1;
	}
	if ((addExpo) && (ao <= 1)){
		printf("ERROR:\033[1;31m Exponential model used without energy input!\033[0m\n");
		stop = 1;
	}
	if ((noise) && ((res == 1) && (addDelta == 0) && (addExpo == 0)))
		printf("WARNING:\033[0;31m Thermostat used even though the system is not dissipative!\033[0m\n");
	if ((noise) && ((addDelta) || (addExpo)))
		printf("WARNING:\033[0;31m Thermostat used with energy input based on collision!\033[0m\n");
	if ((damping) && (euler)){
		printf("ERROR:\033[1;31m Can't use damping == 1 with Euler == 1\033[0m\n");
		stop = 1;
	}
	if ((addShear) && (charged || addField || addWell || addWallx || addWally || addWallz || addCircularWall || damping)){
		printf("ERROR:\033[1;31m Can't use addShear == 1 with charged/Field/Well/Walls/damping == 1\033[0m\n");
		stop = 1;
	}

	printf("\n");

	if ((ERROR) && (stop))
		exit(3);
}
