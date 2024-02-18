#include "raylib.h"
#include "EDMD.h"

extern double Lx;
extern double Ly;
extern int N;
extern particle* particles;
extern node** eventList;
extern int addWallx;
extern int addWally;
extern int addCircularWall;
extern int noise;
extern double delta;
extern bool addDelta;
extern bool addDoubleDelta;
extern double ts;
extern double deltaM;
extern bool damping;
extern bool addWell;
extern bool addField;
extern bool polydispersity;
extern bool thermoWall;
extern double vr;
extern double dtnoise;
extern double t;
extern double gamm;
extern double T;
extern double field;
extern double res;
extern double U;
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

typedef struct window window;
struct window{
    int screenWidth;
    int screenHeight;
    double factor;
    int start;
    double xGUI;
    double yGUI;
    Camera2D cam;
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
    Color* colorArray;
    double (*colorFunction)(particle*);
    double* colorFunctionArray;
    bool running;
    bool wallMoving;
    bool spacePressed;
    bool editVx;
    bool editVy;
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
    int resting;
    int start;
    int end;
};

void getInput(window* screenWindow, state* screenState);
void draw(int argc, char *argv[], window* screenWindow, state* screenState);

void drawParticlesAndBox(double factor, state* screenState, float zoom);
void GuiSliderBarDouble(Rectangle bounds, const char *textLeft, const char *textRight, double *value, double minValue, double maxValue);
int getWhatsUnderClick(window* screenWindow, state* screenState);
void addEventInput(double tin);
Color colorSelect(double value, Color* colorArray);
double colorVelocity(particle* p);
double colorCollision(particle* p);
double colorRadius(particle* p);
double colorBOOP(particle* p1);
void* computeStructureFactor(void* arg);
void normalizeStruct();
void awaitStructFactor();
void asyncStructFactor();
void threadPoolInit();
bool finishedStructComputation();
void reset(int argc, char* argv[], double* factor, state* screenState);
int doubleBox(double* ptr, char* text, bool* activate, Rectangle position);
window graphicalInit();
state GUIinit();
void graphicsFree(state* screenState);