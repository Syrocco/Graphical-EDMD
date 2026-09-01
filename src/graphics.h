#include "raylib.h"
#include "EDMD.h"

extern double Lx;
extern double Ly;
#if THREE_D
extern double Lz;
#endif
extern int N;
extern particle* particles;
extern node** eventList;
extern int addWallx;
extern int addWally;
extern int addWallz;
extern int addMidWall;
extern int addCircularWall;
extern int noise;
extern double delta;
extern bool addDelta;
extern bool addDoubleDelta;
extern bool addDeltaTangent;
extern bool addDoubleDeltaTangent;
extern int addEvolvingDelta;
extern int addEnergy;
extern int addActiveDelta;
extern double ts;
extern double deltaM;
extern bool damping;
extern bool addWell;
extern bool addField;
extern bool polydispersity;
extern bool thermoWall;
extern double vr;
extern int automaticGrowthRate;
extern double dtnoise;
extern double t;
extern double gamm;
extern double T;
extern double field;
extern double res;
extern double aspectRatio;
extern double U;
extern bool charged;
extern double proportionPositivelyCharged;
extern double proportionNeutralyCharged;
extern double sig;
extern double phi;
extern double sizeratio;
extern double fractionSmallN;
extern particle** cellList;
extern int Nxcells;
extern int Nycells;
extern long unsigned int ncol;
extern long unsigned int ncross;
extern double paulTime;
extern int actualPaulList;
extern double dtime;
extern double Einit;
extern int load;
extern double firstScreen;

extern void (*collisionEvent)(int);
extern void (*doTheCollision)(void);
extern void (*freeFly)(particle*);
extern void (*doTheWall)(void);
extern void (*crossingEvent)(int);


extern int S1;
extern int Hex;
extern int interfaceThermo;
extern int interfaceThermoTemp;
extern int clusterThermo;
extern int dumpCluster;
extern int killCluster;

extern int strucThermo;
extern int areaThermo;
extern int boopThermo;
extern int pcfThermo;
extern int pcfBondOrderThermo;
extern int pcfg6Thermo;
extern int critical;
extern int snapshotCritical;
extern int version;
extern int outputVersion;
extern int periodicDumps;
extern int ther;

extern char fileName[350];
extern FILE* fichier;

typedef struct window window;
struct window{
    int screenWidth;
    int screenHeight;
    double factor;
    int start;
    double xGUI;
    double yGUI;
    Camera2D cam;
    #if THREE_D
    Camera3D cam3D;
    double cameraBoxLx;
    double cameraBoxLy;
    double cameraBoxLz;
    #endif
};

typedef struct state state;
struct state{
    int leftClicked;
    int selected;
    particle* particleUnderClick;
    Color particleArray[2];
    Color particleColor;
    bool colorEditing;
    int colorParam;
    bool colorEditing2;
    int colorParam2;
    bool arrow;
    Color* colorArray;
    double (*colorFunction)(particle*);
    double* colorFunctionArray;
    int colorFunctionCapacity;
    bool running;
    bool stepRequested;
    bool wallMoving;
    bool spacePressed;
    bool editVx;
    bool editVy;
    #if THREE_D
    bool editVz;
    #endif
    bool editM;
    bool wallEditing;
    int wallParam;
};

typedef struct position position;
struct position{
	double x, y;
};

typedef struct threadArg threadArg;
struct threadArg{
    int start;
    int end;
};

void getInput(window* screenWindow, state* screenState);
void draw(int argc, char *argv[], window* screenWindow, state* screenState);

#if !THREE_D
void drawParticlesAndBox(double factor, state* screenState, float zoom, Rectangle viewport);
#endif
int getWhatsUnderClick(window* screenWindow, state* screenState);
void addEventInput(double tin);
Color colorSelect(double value, Color* colorArray);
double colorVelocity(particle* p);
double colorCollision(particle* p);
double colorRadius(particle* p);
#if !THREE_D
double colorBOOP(particle* p1);
#endif
double colorCharge(particle* p);
void* computeStructureFactorGraphic(void* arg);
void normalizeStruct();
void awaitStructFactor();
void asyncStructFactor();
void threadPoolInit();
bool finishedStructComputation();
window graphicalInit();
state GUIinit();
void graphicsFree(state* screenState);
/* Keeps the event-loop run flag and ImGui toolbar state in sync. */
void graphicsPauseSimulation(state* screenState);
/* Upper simulation time reachable at the current wall-clock instant.  An
 * unpaced run returns +infinity.  The event loop may advance the global clock
 * up to this limit when the next discrete event lies later. */
double graphicsSimulationTimeLimit(double currentSimulationTime);
/* Requested render cap.  Browser hosts use it to skip animation callbacks
 * without blocking the JavaScript thread. */
int graphicsTargetFrameRate(void);
void graphicsReportRuntimeError(const char* message);
void graphicsReportSimulationFault(const char* message);
void saveSnapshotGraphics();
