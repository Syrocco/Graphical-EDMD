#ifndef EDMD
#define EDMD

#include <stdbool.h>
#include <stdio.h>

#define THREE_D 0
#ifndef THREE_D
#    define THREE_D 0
#endif

#define TANGENTIAL 0

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
};

typedef struct arguments arguments;
struct arguments{
	int argc;
	char **argv;
};

typedef struct particle particle;
struct particle{
	double rad, x, y, vx, vy, m, lastColl, t, vr;
	particle *prv, *nxt;
	int num, type, numberOfParticlesInWell, crossX, crossY, synchro;
	int* particlesInWell;
	unsigned long int coll; //coll = counter of collision at collision
	int charge;
	#if THREE_D
	double vz, z;
	int cell[3], crossZ;
	#else
	int cell[2];
	#endif
	#if TANGENTIAL
    double omega;     // Angular velocity
    double J;         // Moment of inertia
    #endif

};

typedef struct node node;
struct node{
	node *lft, *rgt, *top;
	int i, j, q; 
	double t;
	enum event type;

	unsigned long int collActual; //collActual = counter of collision at event prediciton
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
int PBCcellX(double a);
int PBCcellY(double a);
int PBCcellZ(double a);
void shearCorrection(double* x, double* vx, int sign, particle* p);
void correctDistances(particle* p1, particle* p2, double lat2, double* dx, double* dvx);

void addEventToQueue(node* toAdd);
void addEventToTree(node* toAdd);
void removeEventFromQueue(node* toRemove);
void removeEventFromTree(node* toRemove);
void addNextPaulEvent();
node* findNextEvent();



void crossingEventGrow(int i);
void crossingEventNormal(int i);
void addCrossingEvent(int i, int info, double tCross);
void doTheCrossing();

double collisionTimeNormal(particle* p1, particle* p2);
double collisionTimeGrow(particle* p1, particle* p2);
void collisionEventGrow(int i);
void collisionEventNormal(int i);
void addCollisionEvent(int i, int j, double tColl);
void doTheCollisionNormal();
void doTheCollisionGrow();

void addEventScreenshot(double tscreen);
void takeAScreenshot();

void addEventThermo(double tscreen);
void takeAThermo();

void addEventGrow(double time);
void stopGrow();

void addEventNoise(double tnoise);
void addNoise();

void addEventUpdate(double tupdate);
void updateT();

void addEventUmbrella(double tscreen);
void doUmbrella();

void freeFlyGrow(particle* p);
void freeFlyNormal(particle* p);
void saveTXT();
void saveThermo();
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

double separateTime(particle* p1, particle* p2);
double sphereTime(particle* p1, particle* p2);

void addOutEvent(int i, int j, double tColl);
void addInEvent(int i, int j, double tColl);

void pairsInit();

void addParticleInWellList(particle* p1, int num);
void removeParticleInWellList(particle* p1, int num);
int isParticleInWellList(particle* p1, int num);

#endif






