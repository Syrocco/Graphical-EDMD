#ifndef EDMD
#define EDMD

#include <stdbool.h>
#include <stdio.h>

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
	int num, type, numberOfParticlesInWell, cell[2], crossX, crossY, synchro;
	int* particlesInWell;

	unsigned long int coll; //coll = counter of collision at collision

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


void addToCell(int i);
void removeFromCell(int i);
int coordToCell(double a, int x);
int PBCcellX(double a);
int PBCcellY(double a);

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


void freeFlyGrow(particle* p);
void freeFlyNormal(particle* p);
void saveTXT();
void saveThermo();
void initializeThermo();
double sign(double x);
double logTime(double time);
double drand(double min, double max);
void randomGaussian(particle* p);
void optimizeGrowConstant();
void PBC(double* dx, double* dy);
void PBCpostX(double* val);
void PBCpostY(double* val);
double PBCinsideCellX(double dx);
double PBCinsideCellY(double dy);
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






