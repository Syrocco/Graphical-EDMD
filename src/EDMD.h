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

#define THREE_D 0
#ifndef THREE_D
#    define THREE_D 0
#endif

#define TANGENTIAL 0

#define INTRUDER 0

enum event{
	CELLCROSS,
	COLLISION, 
	SCREENSHOT,
	THERMO,
	ADDINGNOISE,
	WALL, 
	UPDATE,
	OUT, IN,
	GROWSTOP,
	UMBRELLA,
	MULLER_PLATHE,
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
	int num, type, numberOfParticlesInWell, crossX, crossY, synchro;
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



void doOut();
void doIn();

void addWallEvent(int i, int xy, double tColl);
void doTheWallGrow();
void doTheWallNormal();

void constantInit(int argc, char *argv[]);
void particlesInit();
void cellListInit();
void boxConstantHelper();
void eventListInit();
void freeArrays();
void initThermo();

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
void addEventMullerPlathe(double tswap);
void takeAThermo();

void addEventGrow(double time);
void stopGrow();

void addEventNoise(double tnoise);
void addNoise();

void addEventUpdate(double tupdate);
void updateT();

void addEventUmbrella(double tscreen);
void doUmbrella();

void initMullerPlathe();
void doMullerPlatheSwap();

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
int cmpDouble(const void * a, const void * b);
void normalizePhysicalQ();
void customName();
int mygetline(char* str, FILE* f);
double resCoeff(double v);
void format();
void printInfo();
void printClose();
void runningCheck();

EDMD_PURE double separateTime(const particle* p1, const particle* p2);
EDMD_PURE double sphereTime(const particle* p1, const particle* p2);

void addOutEvent(int i, int j, double tColl);
void addInEvent(int i, int j, double tColl);

void pairsInit();

void addParticleInWellList(particle* p1, int num);
void removeParticleInWellList(particle* p1, int num);
int isParticleInWellList(particle* p1, int num);

void collisionEventIntruder();
#endif
