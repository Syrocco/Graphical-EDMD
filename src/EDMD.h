#ifndef EDMD
#define EDMD

#include <stdbool.h>
#include <stdio.h>

#if defined(__GNUC__) || defined(__clang__)
#define EDMD_NOINLINE __attribute__((noinline))
#define EDMD_PURE __attribute__((pure))
#else
#define EDMD_NOINLINE
#define EDMD_PURE
#endif


#define never 100000000000000000000000000.0

#ifndef THREE_D
#    define THREE_D 0
#endif

#ifndef TANGENTIAL
#    define TANGENTIAL 0
#endif

#ifndef INTRUDER
#    define INTRUDER 0
#endif

enum event{
	CELLCROSS,
	COLLISION, 
	SCREENSHOT,
	THERMO,
	ADDINGNOISE,
	WALL, 
	OUT, IN,
	GROWSTOP,
	#if INTRUDER
	INTRUDER_COLLISION,
	#endif
};

typedef struct arguments arguments;
struct arguments{
	int argc;
	char **argv;
};

typedef struct particle particle;
struct particle{
	double x, y, vx, vy, rad, t, vr;
	particle *nxt, *prv;
	double m, lastColl;
	unsigned long int coll; //coll = counter of collision at collision
	int* particlesInWell;
	#if THREE_D
	int cell[3], crossZ;
	double vz, z;
	#else
	int cell[2];
	#endif
	int num, type, numberOfParticlesInWell, particlesInWellCapacity;
	int crossX, crossY, synchro;
	int charge;
	#if TANGENTIAL
    double omega;     // Angular velocity
    double J;         // Moment of inertia
    #endif
	#if INTRUDER
	double normalAtCollision[2];
	#endif

};

/* A value-only particle state evaluated at an explicit simulation time. */
typedef struct particlePrediction particlePrediction;
struct particlePrediction{
	double x, y, vx, vy, rad;
	#if THREE_D
	double z, vz;
	#endif
};

enum particleFlightPhysics{
	BALLISTIC_FLIGHT,
	CONSTANT_FIELD_FLIGHT,
	LINEAR_DAMPING_FLIGHT,
	GROWTH_FLIGHT,
};

typedef struct node node;
struct node{
	node *lft, *rgt, *top;
	double t;
	unsigned long int collActual; //collActual = counter of collision at event prediciton
	int i, j, q;
	enum event type;
};

/* A graphical run keeps a damaged state available for inspection instead of
 * terminating the process.  Event processing remains locked until a fresh
 * simulation is constructed.  Headless builds retain their fail-fast exit
 * behavior. */
void simulationClearFault(void);
bool simulationHasFault(void);
const char* simulationFaultMessage(void);



void doOut();
void doIn();

void addWallEvent(int i, int xy, double tColl);
void doTheWallGrow();
void doTheWallNormal();

void constantInit(int argc, char *argv[]);
void particlesInit();
void adaptCellGridToParticles();
void cellListInit();
void boxConstantHelper();
/* Derive complete periodic triangular/HCP lattice counts.  In 2D nz is set
 * to one; in 3D ny and nz are even so the row and AB-layer offsets close
 * exactly through the periodic boundaries. */
bool hexagonalLatticeDimensions(int requestedCount, double requestedAspect,
	int* nx, int* ny, int* nz);
void eventListInit();
/* Incremental graphical reconfiguration helpers.  They preserve particles,
 * simulation time, physical counters, output streams, and unrelated
 * auxiliary events. */
void synchronizeParticlesAtCurrentTime(void);
void rescheduleParticleEventsAtCurrentTime(void);
void rescheduleThermostatEventAtCurrentTime(void);
void rebuildPotentialMembershipsAtCurrentTime(void);
void rebuildSpatialIndexAtCurrentTime(void);
void rebaseRuntimeShear(bool wasActive, double oldRate,
	bool isActive, double newRate);
void resetRuntimeShearPhase(void);
void freeArrays();
void initThermo();
bool outputLifecycleOpen(char* error, size_t errorCapacity);
void outputLifecycleClose(void);
void resetOutputSamplingState(void);
bool outputThermoSamplingRequested(void);

/* Total normal impulse delivered by physical wall events in the current run.
 * Unlike the legacy `dp` sampling accumulator, this value is never cleared by
 * thermodynamic output.  Consumers can therefore form independent intervals
 * without racing the file-output cadence. */
extern double cumulativeWallImpulse;

void addToCell(int i);
void removeFromCell(int i);
int coordToCell(double a, int x);
int PBCcellX(int a);
int PBCcellY(int a);
int PBCcellZ(int a);
void shearCorrection(double* x, double* vx, int sign, particle* p);
void correctDistances(const particle* p1, const particle* p2,
	double p2x, double p2vx, double* dx, double* dvx);

void addEventToQueue(node* toAdd);
void addEventToTree(node* toAdd);
void removeEventFromQueue(node* toRemove);
void removeEventFromTree(node* toRemove);
void addNextPaulEvent();
node* findNextEvent();



EDMD_NOINLINE void crossingEventGrow(int i);
EDMD_NOINLINE void crossingEventNormal(int i);
void addCrossingEvent(int i, int info, double tCross);
void doTheCrossing();

EDMD_PURE double collisionTimeNormal(const particle* p1, const particle* p2);
EDMD_PURE double collisionTimeGrow(const particle* p1, const particle* p2);
double circularWallTimeNormal(const particle* p);
double circularWallTimeGrow(const particle* p);
/* Normal velocity after impact with a confining boundary moving inward at
 * `inwardSpeed`.  Positive input points out of the accessible region. */
EDMD_PURE double confiningWallNormalVelocityAfter(double normalVelocity,
	double inwardSpeed, double restitution);
EDMD_NOINLINE void collisionEventGrow(int i);
EDMD_NOINLINE void collisionEventNormal(int i);
void addCollisionEvent(int i, int j, double tColl);
void doTheCollisionNormal();
void doTheCollisionGrow();
void doTheCollisionIntruder();

void addCollisionIntruderEvent(int i, int j, double tColl);

void addEventScreenshot(double tscreen);
void takeAScreenshot();

void addEventThermo(double tscreen);
void takeAThermo();

void addEventGrow(double time);
void stopGrow();

void addEventNoise(double tnoise);
void addNoise();
/* Computes the multiplicative factor used by the global velocity-rescale
 * thermostat.  A false result means that rescaling would divide by zero or
 * introduce a non-finite velocity. */
bool velocityRescaleFactor(double kineticEnergy, int particleCount,
	double targetTemperature, double* factor);

void setTransverseVelocityCorrelationConfig(int enabled, int modeM, int modeN,
	int samples, double timeSpan, int checkpointInterval);
void getTransverseVelocityCorrelationConfig(int* enabled, int* modeM,
	int* modeN, int* samples, double* timeSpan, int* checkpointInterval);
int transverseVelocityCorrelationIsEnabled(void);
double transverseVelocityCorrelationSampleInterval(void);
void configureTransverseVelocityCorrelationSampling(double* thermoInterval,
	int* useLogSpacing);
int initTransverseVelocityCorrelation(double Lx, double Ly, int numParticles,
	const char* outputName);
void updateTransverseVelocityCorrelation(const particle* particles,
	double sampleTime);
void finalizeTransverseVelocityCorrelation(void);
int initStructureFactor(double qmax, double Lx, double Ly, int num, FILE* file,
	double doFQT);
void freeStructureFactor(void);

EDMD_NOINLINE void freeFlyGrow(particle* p);
EDMD_NOINLINE void freeFlyNormal(particle* p);
/* parameter is acceleration for field flight and gamma for damping flight. */
EDMD_PURE particlePrediction predictParticle(const particle* p,
	double targetTime, enum particleFlightPhysics physics, double parameter);
EDMD_PURE particlePrediction predictFreeFlyGrow(const particle* p, double targetTime);
EDMD_PURE particlePrediction predictFreeFlyNormal(const particle* p, double targetTime);
void saveTXT();
void saveThermo();
void saveTempProfile();
double sign(double x);
double logTime(double time);
double drand(double min, double max);
void randomGaussian(particle* p, double T);
void optimizeGrowConstant();
#if THREE_D
void PBC(double* dx, double* dy, double* dz);
#else
void PBC(double* dx, double* dy);
#endif
void PBCpostX(double* val);
void PBCpostY(double* val);
void PBCpostZ(double* val);
double PBCinsideCellX(double dx);
double PBCinsideCellY(double dy);
double PBCinsideCellZ(double dz);
void physicalQ();
/* Rebuilds time-integrated stress/heat state after particle initialization. */
#if (STRESS || HEAT)
void initializeTransportAccumulators(void);
void rebaseTransportAccumulatorsAfterModelChange(void);
/* Bracket a discrete velocity update which is external to a pair/wall
 * collision (thermostats, swaps, or a paused GUI edit).  They are no-ops in
 * builds without STRESS/HEAT, so callers do not need preprocessor branches. */
void beginTransportVelocityBatch(void);
void recordTransportVelocityChange(const particle* p,
	double oldVx, double oldVy, double oldVz);
#endif
int cmpDouble(const void * a, const void * b);
void normalizePhysicalQ();
void customName();
bool customSnapshotName(char* destination, size_t capacity);
double resCoeff(double v);
void format();
void printInfo();
void printClose();
void runningCheck();

EDMD_PURE double separateTime(const particle* p1, const particle* p2);
EDMD_PURE double sphereTime(const particle* p1, const particle* p2);
/* Solve one radial square-step crossing. deltaPotential is the potential
 * energy after the crossing minus the value before it. The return value is
 * true for transmission and false for an elastic reflection. */
bool potentialStepNormalVelocities(double firstMass, double secondMass,
	double firstVelocity, double secondVelocity, double deltaPotential,
	double* firstVelocityAfter, double* secondVelocityAfter);

void addOutEvent(int i, int j, double tColl);
void addInEvent(int i, int j, double tColl);

void pairsInit();

void addParticleInWellList(particle* p1, int num);
void removeParticleInWellList(particle* p1, int num);
int isParticleInWellList(const particle* p1, int num);

void collisionEventIntruder();
#if INTRUDER
double getIntruderKineticModelBurnTime(void);
void setIntruderKineticModelBurnTime(double burnTime);
#endif
#endif
