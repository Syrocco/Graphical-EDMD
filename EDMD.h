#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<stdio.h>
#include<time.h>
#include<sys/stat.h>
#include<sys/types.h>

typedef struct particle particle;
struct particle{
	double rad, x, y, vx, vy, m;
	particle *prv, *nxt;
	int cell[2]; //very stupid
	int num, type;
	unsigned long int coll; //coll = counter of collision at collision
	double t;
	int crossX, crossY;
	double lastColl;
	unsigned int* particlesInWell;
	int numberOfParticlesInWell;
};

typedef struct node node;
struct node{
	node *lft, *rgt, *top;
	int i, j; 
	int type;
	unsigned long int collActual; //collActual = counter of collision at event prediciton
	double t;
	int q;
};

void doOut();
void doIn();

void addWallEvent(int i, int xy, double tColl);
void doTheWall();

void constantInit(int argc, char *argv[]);
void particlesInit();
void cellListInit();
void eventListInit();
void freeArrays();


void addToCell(int i);
void removeFromCell(int i);
int coordToCell(double a, int x);
int PBCcell(double a, int x);


void addEventToQueue(node* toAdd);
void addEventToTree(node* toAdd);
void removeEventFromQueue(node* toRemove);
void removeEventFromTree(node* toRemove);
void addNextPaulEvent();
node* findNextEvent();



void crossingEvent(int i);
void addCrossingEvent(int i, int info, double tCross);
void doTheCrossing();

double collisionTime(particle* p1, particle* p2);
void collisionEvent(int i);
void addCollisionEvent(int i, int j, double tColl);
void doTheCollision();

void addEventScreenshot(double tscreen);
void takeAScreenshot();

void addEventThermo(double tscreen);
void takeAThermo();

void addEventNoise(double tnoise);
void addNoise();

void addEventUpdate(double tupdate);
void updateT();


void freeFly(particle* p);
void saveTXT();
void saveThermo();
double logTime(double time);
double drand(double min, double max);
void randomGaussian(particle* p);
void PBC(double* dx, double* dy);
void PBCpost(double* val, int x);
double PBCinsideCell(double dx, int x);
void physicalQ();
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


