#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<stdio.h>
#include<time.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<pthread.h>

#include <getopt.h>
#include<time.h>

#define G 1

#if G
#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h" 
#include<stdbool.h>
#include "color.h"

#endif

#if G
typedef struct position position;
struct position{
	double x, y;
};

typedef struct threadArg threadArg;
struct threadArg{
    int start;
    int end;
};
#endif

typedef struct particle particle;
struct particle{
	double rad, x, y, vx, vy, m, lastColl, t;
	particle *prv, *nxt;
	int num, type, numberOfParticlesInWell, cell[2], crossX, crossY, synchro;
	unsigned int* particlesInWell;

	unsigned long int coll; //coll = counter of collision at collision

};

typedef struct node node;
struct node{
	node *lft, *rgt, *top;
	int i, j, type, q; 
	double t;

	unsigned long int collActual; //collActual = counter of collision at event prediciton
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
void takeAScreenshot(int argc, char *argv[]);

void addEventThermo(double tscreen);
void takeAThermo();

void addEventNoise(double tnoise);
void addNoise();

void addEventUpdate(double tupdate);
void updateT();


void freeFly(particle* p);
void saveTXT();
void saveThermo();
void initializeThermo();
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

#if G
void getInput();
void drawParticles();
void draw(int argc, char *argv[]);
void GuiSliderBarDouble(Rectangle bounds, const char *textLeft, const char *textRight, double *value, double minValue, double maxValue);
int getParticleUnderClick();
void addEventInput(double tin);
Color colorSelect(double value);
double colorVelocity(particle* p);
double colorCollision(particle* p);
double colorBOOP(particle* p1);
void* computeStructureFactor(void* arg);
void normalizeStruct();
void awaitStructFactor();
void asyncStructFactor();
void threadPoolInit();
void reset(int argc, char *argv[]);
int doubleBox(double* ptr, char* text, bool* activate, Rectangle position);
#endif






