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
extern bool addWell;
extern bool addField;
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
extern particle*** cellList;
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
extern bool running;

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

typedef struct position position;
struct position{
	double x, y;
};

typedef struct threadArg threadArg;
struct threadArg{
    int start;
    int end;
};

void getInput(window* screenWindow);
void draw(int argc, char *argv[], window* screenWindow);

void drawParticlesAndBox(double factor);
void GuiSliderBarDouble(Rectangle bounds, const char *textLeft, const char *textRight, double *value, double minValue, double maxValue);
int getParticleUnderClick(window* screenWindow);
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
window graphicalInit();