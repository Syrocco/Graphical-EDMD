#include "EDMD.h"
#include "mersenne.c" 
#include "quartic.c" 
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<stdio.h>
#include<time.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<stdbool.h>
#include"parser.h"

#include<getopt.h>
#include<time.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


#ifndef G
#define G 0
#endif
#define SIGN 1
#if G
#include "graphics.h"
#include "raylib.h"
#include<pthread.h>
#endif




#define CGREEN printf("\033[1;32m")
#define CWHITE printf("\033[0;37m")
#define CCYAN printf("\033[1;96m")

#define ERROR 0
#define never 100000000000000000000000000.0

//Values if no load file
int N = 100;
double phi = 0.05;
double sizeratio = 0.5;
double fractionSmallN = 0;
double aspectRatio = 1;
int controlLenght = 0;

int Nbig, Nsmall;
double px, py, E, pressure, dxr, dyr, dist, active, Ep, MSD;
double COM_x = 0;
double COM_y = 0;
double Lx, Ly;
double lastScreen = 0;
double lastCollNum = 0;



int Nycells;
int Nxcells;
double halfLx; // = Lx/2
double halfLy; // = Ly/2
double cellxSize; // L/Ncells
double cellySize; // L/Ncells
double halfxC; // = cellSize/2;
double cellxFac; // = 1/cellSize
double halfyC; // = cellSize/2;
double cellyFac; // = 1/cellSize

//for now we only need 3 more events
const int suppSizeEvent = 10;

double t = 0;

double collTerm = 0;
double dp = 0;
double dpLeft = 0;
double dpRight = 0;
double dpMid = 0;

unsigned long int ncol = 0;
unsigned long int ncross = 0;
char fileName[350];
char thermoName[350];
char buffer[255];
FILE* file;


//load a configuration if == 1 (if nothing specified it loads from data.txt)
int load = 2;
int S1 = 0;

//duration of simulation
double tmax = 1100;  
//time between each screenshots
double dtime = 100;
double firstScreen = 0;

double dtimeThermo = 100;
double firstThermo = 0;

//if -1, screenshot will be taken at constant interval of dtimeThermo
double nextScreen = -1;



#if G
bool addWell = true;
bool addField = false;
int noise = 1;
int addWally = 0;
int addWallx = 0;
int addCircularWall = 0;
bool damping = false;
bool addDelta = false;
bool addDoubleDelta = false;
const int addEvolvingDelta = 0;
const int addExpo = 0;
const int addEnergy = 0;
bool polydispersity = false;
bool thermoWall = 0;
bool charged = 0;
bool addShear = 0;
bool addShearWall = 0;
#else
const int addWell = 0;
const int addField = 0;
const int noise = 0;
const int addWally = 0;
const int addWallx = 0;
const int addCircularWall = 0;
const int damping = 0;
const int addDelta = 1;
const int addEnergy = 0;
const int addDoubleDelta = 0;
const int addEvolvingDelta = 0;
const int addExpo = 0;
const int polydispersity = 0;
const int thermoWall = 0;
const int charged = 0;
const int addShear = 0;
const int addShearWall = 0;
#endif

//update the thermostat temperature to reach a wanted temperature for the system
const int updating = 0;


//add a wall at x = Lx/2
const int addMidWall = 0;



void (*collisionEvent)(int);
void (*doTheCollision)(void);
void (*freeFly)(particle*);
void (*doTheWall)(void);
void (*crossingEvent)(int);

//if noise use euler solver
const int euler = 0;

//save a file with thermodynamics quantities each dtimeThermo
const int ther = 1;

//reduce the number of properties dumped into the dump files
const int reduce = 0;
const int msd = 0;
double* posxInitial = NULL;
double* posyInitial = NULL;

//value for delta model
double delta = 0.022;

//values for double delta model
double deltaM = 0.05;
double deltam = 0;
double ts = 4;

//value for evolving delta model
double tau = 5;
double deltaMax = 0.1;

//values for exponential model
double vo = 3.5;
double ao = 1.5;

//values for square potential model
double sig = 1.125;
double U = -2;

//Value of the energy input at collision
double deltaE = 50;
double beta = 6;
double taur = 3;
double additionalEnergy = 0.05;

//value for the field, must be < 0
double field = -0.1; 

//Initial temperature
double Einit = 0.03;

//coeff of restitution of the wall
double resW = 1;

//coeff of restitution of particles
double res = 0.95;

//parameter if noise or damping
double gamm = 0.01;
double T = 0;
double expE = 1;
//time between kicks
double dtnoise = 1;


double strainRate = 0.01;

double a = 2;
double b = 5;

double proportionPositivelyCharged = 0.2;
double proportionNeutralyCharged = 0.1;

const int addEnergyThermo = 1;
const int addPressureThermo = 1;
const int addCollisionNumberThermo = 1;
const int addEnergyCenterOfMassThermo = 0;

//parameter to reach a given temperature out of equilibrium with langevin dynamics
double updateTime = 100;
double updateTimeMax = 2000;
double vr;


double sizeRat;


double paulTime = 0;
int actualPaulList = 0;
double dtPaul;
int paulListN;

double t1 = 0;
FILE *fichier;
FILE *thermo;
node *root;
Dump* dump;
//array containing the particles
particle* particles;
//2D array of pointer
particle** cellList;
node** eventList;
node** eventPaul;
//variable containing the nextEvent
node* nextEvent;

void* computeEvolution(void *arg){
	
	#if G
	window screenWindow = graphicalInit();
	state screenState = GUIinit();
	#endif

	arguments* argument = (arguments*)arg;
	int argc = argument->argc;
	char** argv = argument->argv;

	constantInit(argc, argv);
	runningCheck();
	particlesInit();
	cellListInit();
	if (addWell){
		pairsInit();
	}
	eventListInit();
	physicalQ();
	format();
	#if G != 1
	customName();
	fichier = fopen(fileName, "w");
	initThermo();
	while (t <= tmax){
	#else
	screenWindow.factor = GetScreenHeight()/Ly;
	while (!WindowShouldClose()){
		if (!((!screenState.running) || (screenState.wallMoving))){
	#endif

			nextEvent = findNextEvent();
			//double dtEventM = nextEvent->t - t;
			//if (dtEventM < 0.){
			//	printf("%lf %d %d %d | ", dtEventM, nextEvent->type, nextEvent->i, nextEvent->j);
			//}
			
			t = nextEvent->t;
			//printf("%lf %lf %lf %lf\n", particles[0].x, particles[0].y, particles[1].x, particles[1].y);
			removeEventFromQueue(nextEvent);
			switch(nextEvent->type){
				case COLLISION:
					doTheCollision();
					break;
				case WALL:
					doTheWall();
					break;
				case CELLCROSS:
					ncross++;
					doTheCrossing();
					break;
				case ADDINGNOISE:
					addNoise();
					break;
				case SCREENSHOT:
					takeAScreenshot();
					#if G
					draw(argc, argv, &screenWindow, &screenState);
					getInput(&screenWindow, &screenState);
					fflush(stdout);
					#endif

					#if G != 1
					if (noise == 0 && E < 1e-10){// && addWell == 0){
						goto exiting;
					}
					#endif
					break;
				case THERMO:
					takeAThermo();
					fflush(stdout);
					#if G != 1
					if (noise == 0 && E < 1e-10){// && addWell == 0){
						saveTXT();
						goto exiting;
					}
					#endif
					break;
				case UPDATE:
					updateT();
					break;

				case IN:
					doIn();
					break;
				case OUT:
					doOut();
					break;
				case GROWSTOP:
					stopGrow();			
					break;
			}
		#if G
		}
		else{
			draw(argc, argv, &screenWindow, &screenState);
			getInput(&screenWindow, &screenState);
		}
		#endif
	}

	exiting:
		printInfo();
		printClose();
		#if G != 1
		fclose(fichier);
		if (ther)
			fclose(thermo);
		#else
		graphicsFree(&screenState);
		#endif
		freeArrays();
	#if G
	pthread_exit(NULL);
	#else
	return NULL;
	#endif
}


int main(int argc, char *argv[]){

	
	
	//init_genrand(time(NULL));
	init_genrand(666);
	
	#if G
	pthread_t mainThread;
	pthread_create(&mainThread, NULL, computeEvolution,&(arguments){argc, argv});
	pthread_join(mainThread, NULL);
	#else
	computeEvolution(&(arguments){argc, argv});
	#endif
	return 0;
}


/* ------------------------------------------/
/											 /
/	INITIALISATION OF ARRAYS AND CONSTANTS 	 /
/											 /
/-------------------------------------------*/

/* --------------------------------------------
	Initializes the number of particles N
	(for the densely packed config) and
	constants required for the cellList.
-------------------------------------------- */


// Convenient function for graphical part of the program
void boxConstantHelper(){
	halfLx = Lx/2;
	halfLy = Ly/2;
	if (addWell){
		Nycells = (int)(halfLy/sig);
		Nxcells = (int)(halfLx/sig);
	}
	else{
		Nycells = (int)(halfLy);
		Nxcells = (int)(halfLx);
	}

	cellxSize = Lx/Nxcells;
	cellySize = Ly/Nycells;
	halfxC = cellxSize/2;
	halfyC = cellxSize/2;
	cellxFac = 1/cellxSize;
	cellyFac = 1/cellySize;

	dtPaul = 100/(double)N;

    paulListN = N;
}


void constantInit(int argc, char *argv[]){

	char filename[200];
	sprintf(filename, "data.txt");

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
		{NULL, 0, NULL, 0}
	};
	
	
	int c;

	while ((c = getopt_long(argc, argv, "l:N:p:r:d:g:t:D:E:U:R:f:s:x:q:T:a:L:X:e:", longopt, NULL)) != -1){
		switch(c){
			case 'l':
				strcpy(filename, optarg);
				int len = strlen(filename);
				if (len > 0 && filename[len - 1] == 'p') {
					load = 2;
				} else {
					load = 1;
				}
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
				if (addDelta != 1)
					printf("WARNING:\033[0;31m Changing Delta while Delta model is not enabled!\033[0m\n");
				sscanf(optarg, "%lf", &delta);
				break;
			case 'g':
				if ((damping != 1) || (noise != 1))
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
		}
	}

	
	if (load == 1){
		file = fopen(filename, "r");
		mygetline(buffer, file);
		int nValues = sscanf(buffer, "%d %lf %lf\n", &N, &Lx, &Ly);
		if (nValues == 2)
			Ly = Lx;
	}
	else if (load == 2){
		file = fopen(filename, "r");
		dump = dump_open(filename,'r');		
		N = get_natoms(dump);
		char hboxx;
		Lx = get_boxx(&hboxx, 1, dump);
		char hboxy;
		Ly = get_boxy(&hboxy, 1, dump);
	}
	else if (S1){
		Nsmall = N/2;
		Nsmall = (int)sqrt(Nsmall);
		Nsmall = Nsmall*Nsmall;
		Nbig = Nsmall;
		N = Nbig + Nsmall;
		Lx = sqrt(Nsmall*M_PI*(sizeratio*sizeratio + 1)/phi);
		Ly = Lx;
	}
	else{
		double r2 = (1 + fractionSmallN*(sizeratio*sizeratio - 1));
		if (polydispersity){
			r2 = b*tgamma(1 + 2/a)*tgamma(b)/tgamma(1 + b + 2/a);
			Ly = Lx;
		}
		if (addCircularWall){
			Lx = sqrt(4*N/phi*r2);
			Ly = Lx;
		}
		else{
			if (controlLenght == 0){
				Lx = sqrt(M_PI*N/(phi*aspectRatio)*r2);
				Ly = Lx*aspectRatio;
			}
			else if (controlLenght == 1){
				Lx = Ly/aspectRatio;
				N = (int)(phi*Lx*Ly/(M_PI*r2));
				phi = M_PI*N*r2/(Lx*Ly);
			}
			else{
				Ly = Lx*aspectRatio;
				N = (int)(phi*Lx*Ly/(M_PI*r2));
				phi = M_PI*N*r2/(Lx*Ly);
			}

		}
		
	}

	

	if (nextScreen == -1)
		nextScreen = dtime;
	if (firstScreen == -1)
		firstScreen = tmax - 1;
	
	boxConstantHelper();
}

void initThermo(){
	if ((ther)){
		thermo = fopen(thermoName, "w");
		fprintf(thermo, "t ");
		if (addCollisionNumberThermo){
			fprintf(thermo, "Ncol ");
		}
		if (addEnergyThermo){
			fprintf(thermo, "E ");
		}
		if (addPressureThermo){
			fprintf(thermo, "p ");
		}
		if (addCircularWall || addWallx || addWally){
			fprintf(thermo, "pwall ");
		}
		if (addMidWall){
			fprintf(thermo, "pmid dp ");
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
		fprintf(thermo, "\n");
	}	
}

/* ------------------------------------------
	Initializes the array containing all
	the particles and their properties
------------------------------------------ */
void particlesInit(){

	particles = calloc(N, sizeof(particle));
	if (particles == NULL)
		printf("Memory alloc failed");

	
	

	int i = 0;
	if (load == 1){
		particle* p;
		for (i = 0; i < N; i++){
		    p = particles + i;
		    mygetline(buffer, file);
		    sscanf(buffer, "%d %lf %lf %lf %lf\n", &(p->type), &(p->x), &(p->y), &(p->rad), &(p->m));
		    Nbig += p->type;
			double U1 = drand(0, 1);
			double U2 = drand(0, 1);
            p->vx = sqrt(-2*log(U1)/p->m*Einit)*cos(2*M_PI*U2);
            p->vy = sqrt(-2*log(U1)/p->m*Einit)*sin(2*M_PI*U2);
			p->num = i;
			p->t = 0;
			p->coll = 0;
            p->lastColl = 0;
		}
		fclose(file);
		Nsmall = N - Nbig;
		sizeratio = particles[N-1].rad;
		normalizePhysicalQ();
	}
	else if (load == 2){
		int nframes = dump->nframes;
		jump_to_frame(nframes - 1, dump);

		double* x = calloc(N, sizeof(double));
		double* y = calloc(N, sizeof(double));
		double* vx = calloc(N, sizeof(double));
		double* vy = calloc(N, sizeof(double));
		double* rad = calloc(N, sizeof(double));
		double* m = calloc(N, sizeof(double));
		int* type = calloc(N, sizeof(int));

		get_doubleatomprop("x", x, N, dump);
		get_doubleatomprop("y", y, N, dump);
		get_doubleatomprop("vx", vx, N, dump);
		get_doubleatomprop("vy", vy, N, dump);
		get_doubleatomprop("radius", rad, N, dump);
		get_doubleatomprop("m", m, N, dump);
		get_intatomprop("type", type, N, dump);

		
		

		particle* p;
		for (i = 0; i < N; i++){
		    p = particles + i;
			p->type = type[i];
			p->x = x[i];
			p->y = y[i];
			p->vx = vx[i];
			p->vy = vy[i];
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
		sizeratio = particles[N-1].rad;
		free(x);
		free(y);
		free(vx);
		free(vy);
		free(type);
		free(m);
		free(rad);
		dump_close(dump);
	}
	else if (S1){
		particle* p;
		int Nsqrt = (int)sqrt(Nsmall);
		int k = 0;
		float dx = Lx/Nsqrt;
		for (i = 0; i < Nsqrt; i++){
			for (int j = 0; j < Nsqrt; j++){
				p = particles + k;
				p->x = i*dx;
				p->y = j*dx;
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
		
		for (i = 0; i < Nsqrt; i++){
			for (int j = 0; j < Nsqrt; j++){
				p = particles + k;
				p->x = dx/2 + i*dx;
				p->y = dx/2 + j*dx;
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
	else{
		
		double EinitGrow = 0.05;
		particle* p;
		optimizeGrowConstant();
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
					p->x = drand(0.1, Lx - 0.1);
					p->y = drand(0.1, Ly - 0.1);
				}
				while ((p->x - halfLx)*(p->x - halfLx) + (p->y - halfLy)*(p->y - halfLy) > 0.7*halfLx*halfLx); //0.2 because so cool effect !!
			}
			else{
				p->x = drand(2.5, Lx - 2.5); 
				p->y = drand(2.5, Ly - 2.5);
			}
			p->num = i;
			p->rad = 0;
			p->crossX = 0;
			p->crossY = 0;

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


			p->t = 0;
			p->coll = 0;

			double U1 = drand(0, 1);
			double U2 = drand(0, 1);
			p->vx = sqrt(-2*log(U1)/p->m*EinitGrow)*cos(2*M_PI*U2);
			p->vy = sqrt(-2*log(U1)/p->m*EinitGrow)*sin(2*M_PI*U2);

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
			
		}
		if (polydispersity){
			free(targetRad);
		}
	}
	

	

	if (load){
		if (addCircularWall)
			phi = (sizeratio*sizeratio*Nsmall + Nbig)/(halfLx*halfLx);
		else
			phi = M_PI*(sizeratio*sizeratio*Nsmall + Nbig)/(Lx*Ly);
	}
}


/* ---------------------------------------------
	Initializes the list of particles close!
--------------------------------------------- */

void pairsInit(){

	

	particle* p1;

	for(int i = 0; i < N; i++){
		p1 = particles + i;
		p1->particlesInWell = calloc(30, sizeof(unsigned int));
		p1->numberOfParticlesInWell = 0;

		int X = p1->cell[0];
		int Y = p1->cell[1];

		for (int j = -1; j <= 1; j++){
			for (int k = -1; k <= 1; k++){
				particle* p2 = cellList[PBCcellX(X + j)*Nycells + PBCcellY(Y + k)];
				//printf("id: %d vs: %d | PBC %d %d | Nxcells: %d Nycells: %d\n", PBCcellX(X + j)*Nycells + PBCcellY(Y + k), Nxcells*Nycells, PBCcellX(X + j), PBCcellX(Y + k), Nxcells, Nycells);
				while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (p1->num != p2->num){
						double dx = p2->x - p1->x;
						double dy = p2->y - p1->y;
						PBC(&dx, &dy);
						if (sqrt(dx*dx + dy*dy) < sig*(p1->rad + p2->rad) && (isParticleInWellList(p1, p2->num) == 0)) //last condition to prevent double counting for small system
							addParticleInWellList(p1, p2->num);
					}
					p2 = p2->nxt;
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

	cellList = calloc(Nxcells*Nycells, sizeof(particle*));


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
void eventListInit(){

	int totalSize = 2*N + suppSizeEvent;
	eventList = calloc(totalSize, sizeof(node*));
	eventPaul = calloc(paulListN + 1, sizeof(node*));

	for (int i = 0; i < totalSize; i++)
		eventList[i] = calloc(1, sizeof(node));

	if (eventList == NULL)
		printf("malloc failed");


	root = eventList[2*N];
	root->t = never + 1; //arbitrary value, root->rgt will always be NULL.
	root->rgt = NULL;
	root->lft = NULL;
	root->top = NULL;

	
	
	#if G != 1
	
	if (load == 0 && S1 == 0){
		if (firstScreen < 1/vr){
			firstScreen = 1/vr + firstScreen;
		}
		if (firstThermo < 1/vr){
			firstThermo = 1/vr + firstThermo;
		}
	}
	
	addEventThermo(firstThermo + 0.000000001);
	
	
	#endif
	if (load == 0 && S1 == 0){
		tmax += 1/vr;
		addEventGrow(1/vr);
		collisionEvent = &collisionEventGrow;
		doTheCollision = &doTheCollisionGrow;
		freeFly = &freeFlyGrow;
		doTheWall = &doTheWallGrow;
		crossingEvent = &crossingEventGrow;

	}
	else{
		collisionEvent = &collisionEventNormal;
		doTheCollision = &doTheCollisionNormal;
		freeFly = &freeFlyNormal;
		doTheWall = &doTheWallNormal;
		crossingEvent = &crossingEventNormal;
	}
	addEventScreenshot(firstScreen);
	if (noise)
		addEventNoise(dtnoise);
	if (updating)
		addEventUpdate(updateTime);

	for (int i = 0; i < N; i++){
		eventList[i]->i = i;
		eventList[N + i]->i = i;
		crossingEvent(i);
		collisionEvent(i);
	}
	

}


void freeArrays(){
	int totalSize = 2*N + suppSizeEvent;
	if (addWell){
		for (int i = 0; i < N; i++){
			free(particles[i].particlesInWell);
		}	
	}
	free(particles);
	free(cellList);
	for (int i = 0; i < totalSize; i++){
		free(eventList[i]);
	}
	free(eventList);
	free(eventPaul);

	

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
	p->nxt = cellList[X*Nycells + Y];
	cellList[X*Nycells + Y] = p;
	if (p->nxt)
		p->nxt->prv = p;

}

/* ----------------------------------------
   	Removes Particle i from the CellList
-----------------------------------------*/
void removeFromCell(int i){
	particle* p = particles + i;
	if (p->prv == NULL)
		cellList[p->cell[0]*Nycells + p->cell[1]] = p->nxt;
	else
		p->prv->nxt = p->nxt;
	if (p->nxt != NULL)
		p->nxt->prv = p->prv;
}


int coordToCell(double a, int x){
	if (x)
		return (int)(a*cellxFac);
	return (int)(a*cellyFac);
}


int PBCcellX(double a){
	if (a < 0)
		return a + Nxcells;
	else if (a >= Nxcells)
		return a - Nxcells;
	return a;
}

int PBCcellY(double a){
	if (a < 0)
		return a + Nycells;
	else if (a >= Nycells)
		return a - Nycells;
	return a;
}


/* ---------------------------------------/
/										  /
/	Calendar/Event queue utilitaries 	  /
/										  /
/----------------------------------------*/

/* --------------------------------
	Adds an event to the queue
-------------------------------- */
void addEventToQueue(node* toAdd){

	double dt = toAdd->t - paulTime;

	if (dt < dtPaul){
		toAdd->q = actualPaulList;
		addEventToTree(toAdd);
	}
	else{
		int paulListIndex = actualPaulList + dt/dtPaul; //could put long int
		if (paulListIndex >= paulListN){
			paulListIndex -= paulListN;
			if (paulListIndex > actualPaulList - 1)
				paulListIndex = paulListN;
		}
		else if (paulListIndex < 0)
			paulListIndex = paulListN;

		toAdd->q = paulListIndex;
		toAdd->lft = NULL;
		toAdd->rgt = eventPaul[paulListIndex];
		if (toAdd->rgt != NULL)
			toAdd->rgt->lft = toAdd;
		eventPaul[paulListIndex] = toAdd;
	}
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
			else{
				walker = walker->lft;
			}
		}
		else{
			if (walker->rgt == NULL){
				walker->rgt = toAdd;
				break;
			}
			else{
				walker = walker->rgt;
			}
		}
	}
	toAdd->lft = NULL;
	toAdd->rgt = NULL;
	toAdd->top = walker;
}


/* -----------------------------------------
	Removes an event from the event queue
----------------------------------------- */
void removeEventFromQueue(node* toRemove){

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
}



/* -----------------------------------
	Removes an event from the BST
------------------------------------ */
void removeEventFromTree(node* toRemove){

    node* top = toRemove->top;
    node* n;

	if ((toRemove->lft == NULL) && (toRemove->rgt == NULL))
		n = NULL;

	else if (toRemove->lft == NULL){
		n = toRemove->rgt;
		n->top = top;
	}
	else if (toRemove->rgt == NULL){
		n = toRemove->lft;
		n->top = top;
	}
	//http://www.mathcs.emory.edu/~cheung/Courses/171/Syllabus/9-BinTree/BST-delete2.html
	else{
		n = toRemove->rgt;

		//look for successor
		while (n->lft != NULL)
			n = n->lft;

		if (n->top != toRemove) //if (...) then we move successor to toRemove and link what need to be linked
		{
			n->top->lft = n->rgt;
			if (n->rgt)
				n->rgt->top = n->top;
			toRemove->rgt->top = n;
			n->rgt = toRemove->rgt;
		}
		toRemove->lft->top = n;
		n->lft = toRemove->lft;
		n->top = top;
	}
	//link new node to tree
	if (top->lft == toRemove)
		top->lft = n;
	else
		top->rgt = n;
}


/* ---------------------------------------------
	Add the buffer Paul's list into the tree
--------------------------------------------- */

void addNextPaulEvent(){
	actualPaulList++;
	paulTime += dtPaul;
	if (actualPaulList == paulListN){
		actualPaulList = 0;
		node* temp = eventPaul[paulListN];
		eventPaul[paulListN] = NULL;
		while (temp != NULL){
			node* temp2 = temp->rgt;
			addEventToQueue(temp);
			temp = temp2;
		}
	}
	node* temp = eventPaul[actualPaulList];
	while (temp != NULL){
			node* temp2 = temp->rgt;
			addEventToTree(temp);
			temp = temp2;
	}
	eventPaul[actualPaulList] = NULL;
}


/* ---------------------------------------------
	Finds the next event according to the BST
--------------------------------------------- */
node* findNextEvent(){
	node* chosenONE = root->lft;
	
	while (chosenONE == NULL){
		addNextPaulEvent();
		chosenONE = root->lft;
	}
	while (chosenONE->lft != NULL)
		chosenONE = chosenONE->lft;
	return chosenONE;
}

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

	particle p = particles[i];
	double travelTimeX = 10000000;
	double travelTimeY = 10000000;

	int xx;
	int yy = 1;

	// xx and yy indicates the direction to the new cell for example xx = 1 means that we will go from cell Cell[X][Y] to Cell[X - 1][Y]
	//printf("p: %d %lf |||", p.num, p.vx);
	if (p.vx < 0){
		travelTimeX = PBCinsideCellX(p.cell[0]*cellxSize - p.x)/p.vx;
		xx = 1;
		
	}
	else{
		travelTimeX = PBCinsideCellX((1 + p.cell[0])*cellxSize - p.x)/p.vx;
		xx = 2;
		
	}

		

	if (p.vy < 0){
		travelTimeY = PBCinsideCellY(p.cell[1]*cellySize - p.y)/p.vy;
		yy = 3;
	}
	else{
		travelTimeY = PBCinsideCellY((1 + p.cell[1])*cellySize - p.y)/p.vy;
		yy = 4;
	}
	

	if (travelTimeX < travelTimeY)
		addCrossingEvent(i, xx, t + travelTimeX);
	else
		addCrossingEvent(i, yy, t + travelTimeY);
}

void crossingEventNormal(int i){

	particle p = particles[i];
	double travelTimeX = 10000000;
	double travelTimeY = 10000000;

	int xx;
	int yy = 1;

	// xx and yy indicates the direction to the new cell for example xx = 1 means that we will go from cell Cell[X][Y] to Cell[X - 1][Y]
	//printf("p: %d %lf |||", p.num, p.vx);
	if (p.vx < 0){
		travelTimeX = logTime(PBCinsideCellX(p.cell[0]*cellxSize - p.x)/p.vx);
		xx = 1;
		
	}
	else{
		travelTimeX = logTime(PBCinsideCellX((1 + p.cell[0])*cellxSize - p.x)/p.vx);
		xx = 2;
		
	}


	if (addField){
		double vy = p.vy;
		double dy = PBCinsideCellY(p.y - (1 + p.cell[1])*cellySize);
		if (vy > sqrt(2*field*dy)){
			travelTimeY = (-vy + sqrt(vy*vy - 2*field*dy))/field;
			yy = 4;
			
		}
		else{
			dy = PBCinsideCellY(p.y - p.cell[1]*cellySize);
			travelTimeY = (-vy - sqrt(vy*vy - 2*field*dy))/field;
			yy = 3;
		}	
	}

		
	else{
		if (p.vy < 0){
			travelTimeY = logTime(PBCinsideCellY(p.cell[1]*cellySize - p.y)/p.vy);
			yy = 3;
		}
		else{
			travelTimeY = logTime(PBCinsideCellY((1 + p.cell[1])*cellySize - p.y)/p.vy);
			yy = 4;
		}
	}

	if (travelTimeX < travelTimeY)
		addCrossingEvent(i, xx, t + travelTimeX);
	else
		addCrossingEvent(i, yy, t + travelTimeY);
}
/* ---------------------------------------------
	Creates the crossing event and adds it
	to the BST.
--------------------------------------------- */
void addCrossingEvent(int i, int info, double tCross){
	node* toAdd;

	toAdd = eventList[i];
	toAdd->j = info;
	toAdd->type = CELLCROSS;
	toAdd->t = tCross;

	addEventToQueue(toAdd);
}

/* -----------------------------------
	Performs the crossing event
----------------------------------- */
void doTheCrossing(){
	int i = nextEvent->i;
	particle* p = particles + i;
	freeFly(p);
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
				if ((addShear && t >= 1/vr) || (addShear && load)){
					shearCorrection(&(p->x), &(p->vx), -1);
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
				if ((addShear && t >= 1/vr) || (addShear && load)){
					shearCorrection(&(p->x), &(p->vx), 1);
					p->cell[0] = coordToCell(p->x, 1); 
				}
			}
			else
				p->cell[1] = newCell;
			break;
	}

	p->prv = NULL;
	p->nxt = cellList[p->cell[0]*Nycells + p->cell[1]];
	cellList[p->cell[0]*Nycells + p->cell[1]] = p;
	if (p->nxt != NULL)
		p->nxt->prv = p;

	//p->coll++;
	crossingEvent(i);
	removeEventFromQueue(eventList[N + i]);
	collisionEvent(i);
}



/* --------------------------------------
	Calculates the time of collision
	of p1 and p2
--------------------------------------- */
double collisionTimeGrow(particle* p1, particle* p2){

	double lat2 = t - p2->t;
	
	double dvx = p2->vx - p1->vx;
	double dvy = p2->vy - p1->vy;
	double dvr = p1->vr + p2->vr;
	

	double dx = (p2->x + lat2*p2->vx) - p1->x;
	double dy = (p2->y + lat2*p2->vy) - p1->y;
	double dr = 2*sqrt(p1->rad*(p2->rad + lat2*p2->vr));
	PBC(&dx, &dy);
	double b = dx*dvx + dy*dvy - dvr*dr;

	
	


	double v2 = dvx*dvx + dvy*dvy;
	double distOfSquare = dx*dx + dy*dy;
	double det = pow(b, 2) - (v2 - dvr*dvr)*(distOfSquare - dr*dr);




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
		printf("%lf %lf\n", minus, plus);
		return plus;
	}
	if (distOfSquare - dr*dr< -0.01){
		printf("\nERROR:\033[0;31m Overlaps detected during GROWTH between particle %d and particle %d ! Position: %lf %lf %lf %lf, Speed: %lf %lf %lf %lf, Radius: %lf %lf, distance: \033[0m\n", p1->num, p2->num, p1->x, p1->y, p2->x, p2->y, p1->vx, p1->vy, p2->vx, p2->vy, p1->rad, p2->rad);
		exit(3);
	}
	
	return never;

	

}

double collisionTimeNormal(particle* p1, particle* p2){
	

	if ((damping == 1) || (addField == 1)){
		freeFly(p2); //utterly retarded evil dumb trick. With damping == 0, no damping, v = cst, we can calculate dx by (x + lat2*vx).
	}                //with damping, it would be harder, so I just update the position of p2. 

	double lat2 = t - p2->t;

	double dvx = p2->vx - p1->vx;
	double dvy = p2->vy - p1->vy;

	double dx = (p2->x + lat2*p2->vx) - p1->x;
	double dy = (p2->y + lat2*p2->vy) - p1->y;

	if (addShear){
		correctDistances(p1, p2, lat2, &dx, &dvx);
	}

	PBC(&dx, &dy);
	
	double b = dx*dvx + dy*dvy;

	//no collision.
	if (b > 0)
		return never;


	double v2 = dvx*dvx + dvy*dvy;
	double distOfSquare = 4*p1->rad*p2->rad;
	distOfSquare = dx*dx + dy*dy - 4*p1->rad*p2->rad;
	double det = pow(b, 2) - v2*distOfSquare;


	//to delete when confident with life decisions...
	if (distOfSquare < -0.01){
		printf("\nERROR:\033[0;31m Overlaps detected during SIMULATION between particle %d and particle %d ! Position: %lf %lf %lf %lf, Speed: %lf %lf %lf %lf, Radius: %lf %lf, distance: %lf \033[0m\n", p1->num, p2->num, p1->x, p1->y, p2->x, p2->y, p1->vx, p1->vy, p2->vx, p2->vy, p1->rad, p2->rad, dx*dx + dy*dy);
		saveTXT();
		exit(3);
	}

	if (det < 0)
		return never;
	return logTime((-b - sqrt(det))/v2);

}

//time to enter in the square well: https://faculty.biu.ac.il/~rapaport/papers/09b-ptp.pdf
double sphereTime(particle* p1, particle* p2){

	if ((damping == 1) || (addField == 1)){
		freeFly(p2); //utterly retarded evil dumb trick. With damping == 0, no damping, v = cst, we can calculate dx by (x + lat2*vx).
	}                //with damping, it would be harder, so I just update the position of p2. 

	double lat2 = t - p2->t;

	double dvx = p2->vx - p1->vx;
	double dvy = p2->vy - p1->vy;
	double dx = (p2->x + lat2*p2->vx) - p1->x;
	double dy = (p2->y + lat2*p2->vy) - p1->y;
	PBC(&dx, &dy);
	double b = dx*dvx + dy*dvy;

	//no collision.
	if (b > 0)
		return never;


	double v2 = dvx*dvx + dvy*dvy;
	double distOfSquare = dx*dx + dy*dy - sig*(p1->rad + p2->rad)*sig*(p1->rad + p2->rad);
	double det = b*b - v2*distOfSquare;


	if (det < 0)
		return never;

	return logTime((-b - sqrt(det))/v2);
}


//time to leave the square well: https://faculty.biu.ac.il/~rapaport/papers/09b-ptp.pdf
double separateTime(particle* p1, particle* p2){

	if ((damping == 1)|| (addField == 1)){
		freeFly(p2); //utterly retarded evil dumb trick. With damping == 0, no damping, v = cst, we can calculate dx by (x + lat2*vx).
	}                //with damping, it would be harder, so I just update the position of p2. 

	double lat2 = t - p2->t;

	double dvx = p2->vx - p1->vx;
	double dvy = p2->vy - p1->vy;
	double dx = (p2->x + lat2*p2->vx) - p1->x;
	double dy = (p2->y + lat2*p2->vy) - p1->y;
	PBC(&dx, &dy);
	double b = dx*dvx + dy*dvy;



	double v2 = dvx*dvx + dvy*dvy;
	double distOfSquare = dx*dx + dy*dy - sig*(p1->rad + p2->rad)*sig*(p1->rad + p2->rad);
	double det = b*b - v2*distOfSquare;


	if (det < 0)
		return never;

	double temp = (- b + sqrt(det))/v2;
	return logTime(temp);
}
/* ---------------------------------------------
	Finds the next collision of particle i
	by calculating collision time of i with
	its surrouding (using the cellList)
--------------------------------------------- */
void collisionEventNormal(int i){
	int finalPartner = 0;
	double dt = never;
	enum event type = COLLISION;
	particle* p1 = particles + i;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	int xy = 2;
	double dtTemp, dtTemp2;
	int typeTemp;



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
			}
			xy = 0;
		}

		else if ((X == Nxcells - 1) && (p1->vx > 0)){

			dtTemp = logTime((Lx - p1->x - p1->rad)/p1->vx);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
			}
			xy = 0;
		}
	}

	if ((addMidWall) && (p1->type == 1)){
		if ((p1->x > Lx/2 - 4) && (p1->x < Lx/2) && (p1->vx > 0)){
			dtTemp = logTime((Lx/2 - p1->x - p1->rad)/p1->vx);
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
			}
			xy = 0;
		}
		else if ((p1->x < Lx/2 + 4) && (p1->x > Lx/2) && (p1->vx < 0)){
		dtTemp = logTime(((p1->x - p1->rad) - Lx/2)/(-p1->vx));
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
			}
			xy = 0;
		}
	}


	if (addCircularWall){
		if (addField){
			double y0 = p1->y - halfLy;
			double x0 = p1->x - halfLx;
			double rad = p1->rad;
			double e = x0*x0 + y0*y0 - (halfLx - rad)*(halfLx - rad);
			//if (e > -100){ //arbitrary, have to think about it
				double v0x = p1->vx;
				double v0y = p1->vy;
				

				double a = field*field/4;
				double b = field*v0y;
				double c = v0x*v0x + v0y*v0y + field*y0;
				double d = 2*v0x*x0 + 2*v0y*y0;
				
				dt = smallestRoot(a, b, c, d, e);
				xy = 2;	
				type = WALL;
			//}	
		}	
		else{
			double y0 = p1->y - halfLx;
			double x0 = p1->x - halfLx;
			double rad = p1->rad;
			double C = (x0*x0 + y0*y0) - (halfLx - rad)*(halfLx - rad);
			//if (C > -100){
				double v0x = p1->vx;
				double v0y = p1->vy;
				

				double A = (v0x*v0x + v0y*v0y);
				double B = 2*(v0x*x0 + v0y*y0);
				
				dt  = (-B + sqrt(B*B - 4*A*C))/(2*A);
				type = WALL;
				xy = 2;
			//}
		}
	}

	
	for (int j = -1; j <= 1; j++){
		for (int k = -1; k <= 1; k++){
			particle* p2 = cellList[PBCcellX(X + j)*Nycells + PBCcellY(Y + k)];
			while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
				if (p1->num != p2->num){
					if (addWell){
						if (isParticleInWellList(p1, p2->num)){
							dtTemp = collisionTimeNormal(p1, p2);
							dtTemp2 = separateTime(p1, p2);
							if (dtTemp < dtTemp2){
								typeTemp = COLLISION;
							}
							else{
								dtTemp = dtTemp2;
								typeTemp = OUT;
							}
						}
						else{
							dtTemp = sphereTime(p1, p2);
							typeTemp = IN;
						}
					}
					else{
						dtTemp = collisionTimeNormal(p1, p2);
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
	if (addShear){
		if (Y == 0){
			for (int k = 0; k < Nxcells; k++){
				particle* p2 = cellList[k*Nycells + Nycells - 1];
				while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (p1->num != p2->num){
						if (addWell){
							if (isParticleInWellList(p1, p2->num)){
								dtTemp = collisionTimeNormal(p1, p2);
								dtTemp2 = separateTime(p1, p2);
								if (dtTemp < dtTemp2){
									typeTemp = COLLISION;
								}
								else{
									dtTemp = dtTemp2;
									typeTemp = OUT;
								}
							}
							else{
								dtTemp = sphereTime(p1, p2);
								typeTemp = IN;
							}
						}
						else{
							dtTemp = collisionTimeNormal(p1, p2);
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
		else if (Y == Nycells - 1){
			for (int k = 0; k < Nxcells; k++){
				particle* p2 = cellList[k*Nycells];
				while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (p1->num != p2->num){
						if (addWell){
							if (isParticleInWellList(p1, p2->num)){
								dtTemp = collisionTimeNormal(p1, p2);
								dtTemp2 = separateTime(p1, p2);
								if (dtTemp < dtTemp2){
									typeTemp = COLLISION;
								}
								else{
									dtTemp = dtTemp2;
									typeTemp = OUT;
								}
							}
							else{
								dtTemp = sphereTime(p1, p2);
								typeTemp = IN;
							}
						}
						else{
							dtTemp = collisionTimeNormal(p1, p2);
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
	int xy = 2;
	double dtTemp;
	int typeTemp;


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
			}
			xy = 0;
		}

		else if ((X == Nxcells - 1) && (vxPlus > 0)){

			dtTemp = (Lx - p1->x - p1->rad)/vxPlus;
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
			}
			xy = 0;
		}
	}

	if (addCircularWall){
		
			double y0 = p1->y - halfLx;
			double x0 = p1->x - halfLx;
			double v0x = p1->vx;
			double v0y = p1->vy;
			double rad = p1->rad;

			double A = v0x*v0x + v0y*v0y - vrParticle*vrParticle;
			double B = 2*(v0x*x0 + v0y*y0 + (halfLx - rad)*vrParticle);
			double C = (x0*x0 + y0*y0) - (halfLx - rad)*(halfLx - rad);
			dt  = (-B + sqrt(B*B - 4*A*C))/(2*A);
			type = WALL;
			xy = 2;

	}

	for (int j = -1; j <= 1; j++){
		for (int k = -1; k <= 1; k++){
			particle* p2 = cellList[PBCcellX(X + j)*Nycells + PBCcellY(Y + k)];
			while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
				if (p1->num != p2->num){
					dtTemp = collisionTimeGrow(p1, p2);
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
	//printf("%lf %d\n", dt, type);
	if (type == COLLISION){
		addCollisionEvent(i, finalPartner, t + dt);
	}
	else
		addWallEvent(i, xy, t + dt);

}


/* ---------------------------------------------
	Updates the Well List
--------------------------------------------- */
void addParticleInWellList(particle* p1, int num){
	p1->particlesInWell[p1->numberOfParticlesInWell] = num;
	p1->numberOfParticlesInWell += 1;
}

void removeParticleInWellList(particle* p1, int num){
	for (int i = 0; i < p1->numberOfParticlesInWell; i++){
		if (p1->particlesInWell[i] == num){
			p1->numberOfParticlesInWell -= 1;
			if (p1->numberOfParticlesInWell != 0)
				p1->particlesInWell[i] = p1->particlesInWell[p1->numberOfParticlesInWell];
		}
	}
}

int isParticleInWellList(particle* p1, int num){
	for (int i = 0; i < p1->numberOfParticlesInWell; i++){
		if (p1->particlesInWell[i] == num)
			return 1;
	}
	return 0;
}


/* ---------------------------------------------
	Creates the collision event and adds it
	to the BST.
--------------------------------------------- */
void addCollisionEvent(int i, int j, double tColl){
	node* toAdd;

	toAdd = eventList[N + i];
	toAdd->j = j;
	toAdd->type = COLLISION;
	toAdd->t = tColl;
	toAdd->collActual = particles[j].coll; //the collision trick of the article


	addEventToQueue(toAdd);
}

void addInEvent(int i, int j, double tColl){
	node* toAdd;

	toAdd = eventList[N + i];
	toAdd->j = j;
	toAdd->type = IN;
	toAdd->t = tColl;
	toAdd->collActual = particles[j].coll; //the collision trick of the article
	addEventToQueue(toAdd);
}

void addOutEvent(int i, int j, double tColl){
	node* toAdd;

	toAdd = eventList[N + i];
	toAdd->j = j;
	toAdd->type = OUT;
	toAdd->t = tColl;
	toAdd->collActual = particles[j].coll; //the collision trick of the article
	addEventToQueue(toAdd);
}

void addWallEvent(int i, int xy, double tColl){
	node* toAdd;

	toAdd = eventList[N + i];
	toAdd->type = WALL;
	toAdd->t = tColl;
	toAdd->j = xy;

	addEventToQueue(toAdd);
}


void doTheWallGrow(){
	int i = nextEvent->i;
	int xy = nextEvent->j;
	ncol++;
	particle* pi = particles + i;
	freeFly(pi);

	double vrParticle = pi->vr;

	//hacky
	if (xy == 0){
		if (pi->x > Lx - 2){
			pi->vx = -pi->vx - 2*vrParticle;
		}
		else{
			pi->vx = -pi->vx + 2*vrParticle;
		}
	}
	else if (xy == 1){

		if (pi->y > Ly - 2){
			pi->vy = -pi->vy - 2*vrParticle;
		}
		else{
			pi->vy = -pi->vy + 2*vrParticle;
		}


	}
	//Very hacky
	else{
		double theta = atan2(pi->y - halfLx, pi->x - halfLx);
		double nx = cos(theta);
		double ny = sin(theta);
		double vx = pi->vx;
		double vy = pi->vy;
		
		if (0){
			double mix = -(1 + resW)*(sign(vx)*4*vrParticle*nx + sign(vy)*4*vrParticle*ny);
			pi->vx = mix*nx;
			pi->vy = mix*ny;
		}
		else{
			double mix = -(1 + resW)*(pi->vx*nx + pi->vy*ny);
			pi->vx += (mix - 2*vrParticle)*nx;
			pi->vy += (mix - 2*vrParticle)*ny;
		}
	}

	pi->coll++;
	removeEventFromQueue(eventList[i]);
	crossingEvent(i);
	collisionEvent(i);
}

void doTheWallNormal(){
	int i = nextEvent->i;
	int xy = nextEvent->j;
	ncol++;
	particle* pi = particles + i;
	freeFly(pi);

	if (xy == 0){
		dp += fabs((resW + 1)*pi->m*pi->vx);
		if (addMidWall){
			if (pi->vx < 0){
				if (pi->x - 2 < 0){
					dpLeft += fabs((resW + 1)*pi->m*pi->vx);
				}
				else{
					dpMid += fabs((resW + 1)*pi->m*pi->vx);
				}
			}
			else{
				if (pi->x + 2 > Lx){
					dpRight += fabs((resW + 1)*pi->m*pi->vx);
				}
				else{
					dpMid += fabs((resW + 1)*pi->m*pi->vx);
				}
			}
		}
		if (thermoWall == 1){
			pi->vx = -sign(pi->vx)*drand(0, T/pi->m);
		}
		else{
			pi->vx = -resW*pi->vx;
		}
	}
	else if (xy == 1){
		dp += fabs((resW + 1)*pi->m*pi->vy);
		/*
		if (pi->y + 2*pi->rad > Ly){
			pi->vy = drand(-3, 0);
		}
		else{
			pi->vy = drand(0, 5);
		}
		*/
		if (thermoWall){
			pi->vy = -sign(pi->vy)*drand(0, T/pi->m);
		}
		
		else{
			pi->vy = -resW*pi->vy;
		}
		if (addShearWall){
			if (pi->y + 2*pi->rad > Ly){
				pi->vx = 0.1;
			}
			else{
				pi->vx = -0.1;
			}
			if (addDelta){
				pi->vy += delta*sign(pi->vy); 
			}
		}
		//pi->vy = drand(-50, 50);
	}
	else{
		double theta = atan2(pi->y - halfLx, pi->x - halfLx);
		double nx = cos(theta);
		double ny = sin(theta);
		double mix = -(1 + resW)*(pi->vx*nx + pi->vy*ny);
		if (thermoWall == 0){
			pi->vx += mix*nx;
			pi->vy += mix*ny;
		}
		else{
			pi->vx = -nx*drand(0, T);
			pi->vy = -ny*drand(0, T);
		}
		dp += fabs(pi->m*mix);
	}

	pi->coll++;
	removeEventFromQueue(eventList[i]);
	crossingEvent(i);
	collisionEvent(i);
}
/* -----------------------------------
	Performs the collision event
----------------------------------- */
void doTheCollisionGrow(){

	double delta = 0.1;
	double res = 0.75;

	int i = nextEvent->i;
	int j = nextEvent->j;

	particle* pi = particles + i;
	particle* pj = particles + j;


	freeFly(pi);
	if (nextEvent->collActual != pj->coll){ //the collision trick of the article
        collisionEvent(i);
        return;
    }

    ncol++;
	freeFly(pj);

	pi->coll++;
	pj->coll++;

	double dx = pj->x - pi->x;
	double dy = pj->y - pi->y;

	double dvx = pj->vx - pi->vx;
	double dvy = pj->vy - pi->vy;
	double dvr = pi->vr + pj->vr;;


	PBC(&dx, &dy);
	dist = sqrt(dx*dx + dy*dy);
	dxr = dx/dist;
	dyr = dy/dist;

	double invMass = 1/(pi->m + pj->m);
	double collKernel = invMass*(1 + res)*(dxr*dvx + dyr*dvy - dvr);



	

	pi->vx += collKernel*pj->m*dxr - 2*pj->m*invMass*delta*dxr;
	pi->vy += collKernel*pj->m*dyr - 2*pj->m*invMass*delta*dyr;
	pj->vx -= collKernel*pi->m*dxr - 2*pi->m*invMass*delta*dxr;
	pj->vy -= collKernel*pi->m*dyr - 2*pi->m*invMass*delta*dyr;

	/*
	pi->vx = - 2*vr*dxr;
	pi->vy = - 2*vr*dyr;
	pj->vx = + 2*vr*dxr;
	pj->vy = + 2*vr*dyr;
	
	if ((sqrt(pi->vx*pi->vx + pi->vy*pi->vy) < vr) || (sqrt(pj->vx*pj->vx + pj->vy*pj->vy) < vr)){
		printf("%lf et %lf | ", sqrt(pi->vx*pi->vx + pi->vy*pi->vy), sqrt(pj->vx*pj->vx + pj->vy*pj->vy));
	}

	if (pow(pi->vx - pj->vx, 2) + pow(pi->vy - pj->vy, 2) < vr*vr){
		printf("WHO");
	}*/

	

    pi->lastColl = t;
    pj->lastColl = t;

	//recomputes crossing and collision of pi and pj
	removeEventFromQueue(eventList[i]);
	removeEventFromQueue(eventList[j]);
	crossingEvent(i);
	crossingEvent(j);

	removeEventFromQueue(eventList[N + j]);
	collisionEvent(i);
	collisionEvent(j);
}

void doTheCollisionNormal(){

	int i = nextEvent->i;
	int j = nextEvent->j;

	particle* pi = particles + i;
	particle* pj = particles + j;


	freeFly(pi); 
	if (nextEvent->collActual != pj->coll){ //the collision trick of the article
        collisionEvent(i);
        return;
    }

    ncol++;
	freeFly(pj);

	pi->coll++;
	pj->coll++;

	double dx = pj->x - pi->x;
	double dy = pj->y - pi->y;

	

	double dvx = pj->vx - pi->vx;
	double dvy = pj->vy - pi->vy;


	if (addShear){
		correctDistances(pi, pj, 0, &dx, &dvx);
	}

	PBC(&dx, &dy);

	if (addExpo)
		res = resCoeff(sqrt(dvx*dvx + dvy*dvy));
	
	double invMass = 1/(pi->m + pj->m);
	double collKernel = invMass*(1 + res)*(dx*dvx + dy*dvy);
	collTerm += pi->m*pj->m*collKernel;
	double funkyFactor = collKernel/(4*pi->rad*pj->rad);

	if ((addDelta) ||(addDoubleDelta || (addEvolvingDelta))){
		//double tau = -ts*log(1 - drand(0, 1));
		if (addDoubleDelta){

			if ((t - pi->lastColl > ts) && (t - pj->lastColl > ts)){
				delta = deltam;
			}
			else{
				delta = deltaM;
			}
		}
		else if (addEvolvingDelta){
			double dti = t - pi->lastColl;
			double dtj = t - pj->lastColl;
			//delta = deltaMax*(2 - (exp(-dtj/tau) + exp(-dti/tau)))/2;
			//delta = deltaMax*(pow(dtj/tau, 3) + pow(dti/tau, 3));
			delta = deltaMax*pow((2 - (exp(-dtj/tau) + exp(-dti/tau)))/2, 3);
		}

		dist = sqrt(dx*dx + dy*dy);
		collTerm -= pi->m*pj->m*invMass*2*delta*dist;
		dxr = dx/dist;
		dyr = dy/dist;

		pi->vx += funkyFactor*pj->m*dx - 2*pj->m*invMass*delta*dxr;
		pi->vy += funkyFactor*pj->m*dy - 2*pj->m*invMass*delta*dyr;
		pj->vx -= funkyFactor*pi->m*dx - 2*pi->m*invMass*delta*dxr;
		pj->vy -= funkyFactor*pi->m*dy - 2*pi->m*invMass*delta*dyr;
	}
	else if (addEnergy){
		double b = dx*dvx + dy*dvy;
		double distSquared = 4*pi->rad*pj->rad;
		//double dti = t - pi->lastColl;
		//double dtj = t - pj->lastColl;
		double deltaETotal = additionalEnergy;// + deltaE*(pow(1 - exp(-dtj/taur), beta) + pow(1 - exp(-dti/taur), beta));

		double funkyFactor2 = (b - sqrt(res*res*b*b + 4*distSquared*deltaETotal))/(2*distSquared); //lacking mass xDeltaE

		pi->vx += funkyFactor2*pj->m*dx;
		pi->vy += funkyFactor2*pj->m*dy;
		pj->vx -= funkyFactor2*pi->m*dx;
		pj->vy -= funkyFactor2*pi->m*dy;

	}
	else{

		pi->vx += funkyFactor*pj->m*dx;
		pi->vy += funkyFactor*pj->m*dy;
		pj->vx -= funkyFactor*pi->m*dx;
		pj->vy -= funkyFactor*pi->m*dy;
	}


    pi->lastColl = t;
    pj->lastColl = t;
	//recomputes crossing and collision of pi and pj
	removeEventFromQueue(eventList[i]);
	removeEventFromQueue(eventList[j]);
	crossingEvent(i);
	crossingEvent(j);

	removeEventFromQueue(eventList[N + j]);
	collisionEvent(i);
	collisionEvent(j);
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

	// Insanely hacky!
	if ((pi->charge*pj->charge != 0) || (charged == 0)) {

		double dx = pj->x - pi->x;
		double dy = pj->y - pi->y;

		PBC(&dx, &dy);

		double dxr = dx/sqrt(dx*dx + dy*dy);
		double dyr = dy/sqrt(dx*dx + dy*dy);


		//https://scholarworks.rit.edu/cgi/viewcontent.cgi?article=5982&context=theses
		double vi = pi->vx*dxr + pi->vy*dyr;
		double vj = pj->vx*dxr + pj->vy*dyr;
		double mi = pi->m;
		double mj = pj->m;
		double vif, vjf;

		double temporary = mi*mj*(vi - vj)*(vi - vj);

		if (((charged == 1) && (pi->charge*pj->charge < 0)) || ((charged == 0) && (U < 0))){
			addParticleInWellList(pi, j);
			addParticleInWellList(pj, i);
			double atroce = (mi*vi + mj*vj);
			double horrible = sqrt(mi*mi*mj*mj*(vi - vj)*(vi - vj) + 2*mi*mj*(mi + mj)*fabs(U));
			vif = (mi*atroce + horrible)/(mi*(mi + mj));
			vjf = (mj*atroce - horrible)/(mj*(mi + mj));
		}
		else{
			if (fabs(U) < -temporary/(2*(mi + mj))){
				double atroce = (mi*vi + mj*vj);
				double horrible = sqrt(mi*mj*temporary - 2*mi*mj*(mi + mj)*fabs(U));
				vif = (mi*atroce - horrible)/(mi*(mi + mj));
				vjf = (mj*atroce + horrible)/(mj*(mi + mj));
				addParticleInWellList(pi, j);
				addParticleInWellList(pj, i);

			}
			else{ //bounces on the square well
				vif = ((mi - mj)*vi + 2*mj*vj)/(mi + mj);
				vjf = ((mj - mi)*vj + 2*mi*vi)/(mi + mj);
			}
		}


		double dvjx = (vjf - vj)*dxr;
		double dvjy = (vjf - vj)*dyr;
		collTerm += dvjy*dy + dvjx*dx;

		pi->vx += (vif - vi)*dxr;
		pi->vy += (vif - vi)*dyr;
		pj->vx += dvjx;
		pj->vy += dvjy;
	}
	else{
		addParticleInWellList(pi, j);
		addParticleInWellList(pj, i);
	}


	//recomputes crossing and collision of pi and pj
	removeEventFromQueue(eventList[i]);
	removeEventFromQueue(eventList[j]);
	crossingEvent(i);
	crossingEvent(j);

	removeEventFromQueue(eventList[N + j]);
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

	// Insanely hacky!
	if ((pi->charge*pj->charge != 0) || (charged == 0)){
		double dx = pj->x - pi->x;
		double dy = pj->y - pi->y;

		PBC(&dx, &dy);

		double dxr = dx/sqrt(dx*dx + dy*dy);
		double dyr = dy/sqrt(dx*dx + dy*dy);


		// https://scholarworks.rit.edu/cgi/viewcontent.cgi?article=5982&context=theses
		double vi = pi->vx*dxr + pi->vy*dyr;
		double vj = pj->vx*dxr + pj->vy*dyr;
		double mi = pi->m;
		double mj = pj->m;
		double temporary = mi*mj*(vi - vj)*(vi - vj);
		double vif, vjf;
		if (((charged == 1) && (pi->charge*pj->charge < 0)) || ((charged == 0) && (U < 0))){
			if (fabs(U) < temporary/(2*(mi + mj))){ //leaves the square well
				double atroce = (mi*vi + mj*vj);
				double horrible = sqrt(mi*mj*temporary - 2*mi*mj*(mi + mj)*fabs(U));
				vif = (mi*atroce - horrible)/(mi*(mi + mj));
				vjf = (mj*atroce + horrible)/(mj*(mi + mj));
				removeParticleInWellList(pi, j);
				removeParticleInWellList(pj, i);

			}
			else{ //bounces on the square well
				vif = ((mi - mj)*vi + 2*mj*vj)/(mi + mj);
				vjf = ((mj - mi)*vj + 2*mi*vi)/(mi + mj);
			}
		}
		else{
			double atroce = (mi*vi + mj*vj);
			double horrible = sqrt(mi*mj*temporary + 2*mi*mj*(mi + mj)*fabs(U));
			vif = (mi*atroce - horrible)/(mi*(mi + mj));
			vjf = (mj*atroce + horrible)/(mj*(mi + mj));
			removeParticleInWellList(pi, j);
			removeParticleInWellList(pj, i);
		}

		double dvjx = (vjf - vj)*dxr;
		double dvjy = (vjf - vj)*dyr;
		collTerm += dvjy*dy + dvjx*dx;

		pi->vx += (vif - vi)*dxr;
		pi->vy += (vif - vi)*dyr;
		pj->vx += dvjx;
		pj->vy += dvjy;
	}
	else{
		removeParticleInWellList(pi, j);
		removeParticleInWellList(pj, i);
	}
	//recomputes crossing and collision of pi and pj
	removeEventFromQueue(eventList[i]);
	removeEventFromQueue(eventList[j]);
	crossingEvent(i);
	crossingEvent(j);

	removeEventFromQueue(eventList[N + j]);
	collisionEvent(i);
	collisionEvent(j);
}

/* ---------------------------------------------
	Creates the screenshot event and adds it
	to the BST.
--------------------------------------------- */
void addEventScreenshot(double tscreen){
	node* toAdd;

	toAdd = eventList[2*N + 1];
	toAdd->type = SCREENSHOT;
	toAdd->t = tscreen;
	toAdd->j = 0; //whatever.

	addEventToQueue(toAdd);
}

void addEventThermo(double tscreen){
	node* toAdd;

	toAdd = eventList[2*N + 3];
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
void takeAScreenshot(){

	for (int i = 0; i < N; i++){
		freeFly(particles + i);
	}
	
	addEventScreenshot(t + dtime);
	#if G == 0
	saveTXT();
	#else
	physicalQ(); //don't remember why i did this.
	#endif
	
	printInfo();
}

void takeAThermo(){

	if ((damping == 1) || (msd)){
		for (int i = 0; i < N; i++){
			freeFly(particles + i);
		}
	}
	if ((fabs(t - firstThermo) < 0.0001) && (firstThermo != 0)){
		if (msd){
			posxInitial = calloc((unsigned int)N, sizeof(double));
			posyInitial = calloc((unsigned int)N, sizeof(double));
			for (int i = 0; i < N; i++){
				COM_x += particles[i].x + particles[i].crossX*Lx;
				COM_y += particles[i].y + particles[i].crossY*Ly;
				posxInitial[i] = particles[i].x + particles[i].crossX*Lx;
				posyInitial[i] = particles[i].y + particles[i].crossY*Ly;
			}
			COM_x = COM_x/N;
			COM_y = COM_y/N;
		}
		
	}
	physicalQ();
	saveThermo();

	addEventThermo(t + dtimeThermo);
	printInfo();

}

void printInfo(){
	int size = 29;
	int place = (int)(size*t/tmax);
	printf("\r│\033[1;32m%.2e\033[0;37m/%.2e ", t, tmax);
	CCYAN;
	for (int i = 0; i < size; i++){
		if (i < place - 1)
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
	}
	collisionEvent = &collisionEventNormal;
	doTheCollision = &doTheCollisionNormal;
	freeFly = &freeFlyNormal;
	doTheWall = &doTheWallNormal;
	crossingEvent = &crossingEventNormal;
	normalizePhysicalQ();
	
	if (addWell){
		for (int i  = 0; i < N; i++){
			p = particles + i;
			free(p->particlesInWell);
		} 
		pairsInit();
	}
	

	
	for (int i = 0; i < N; i++){
		p = particles + i;
		p->coll++;
		removeEventFromQueue(eventList[N + p->num]);
		collisionEvent(p->num);
		removeEventFromQueue(eventList[p->num]);
		crossingEvent(p->num);
	}
}

void addEventGrow(double time){
	node* toAdd;
	toAdd = eventList[2*N + 5];
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
	toAdd = eventList[2*N + 2];
	toAdd->type = ADDINGNOISE;
	toAdd->t = tnoise;
	toAdd->j = 0; //whatever.

	addEventToQueue(toAdd);
}

/* ---------------------------------------------
	Kicks every particles of the system
--------------------------------------------- */

void addNoise(){
	//double k = 2*M_PI*3/Lx;
	if (noise == 2 || noise == 4){
		physicalQ();
	}
	if (noise == 4){
		for (int i = 0; i < N; i++){
			freeFly(particles + i);
			particles[i].vx -= px/(N*particles[i].m);
			particles[i].vy -= py/(N*particles[i].m);
			particles[i].coll++;
		}
		for (int j = 0; j < N; j++){
			removeEventFromQueue(eventList[j]);
			crossingEvent(j);

			removeEventFromQueue(eventList[N + j]);
			collisionEvent(j);
		}
	}
	else if (noise == 3){
		for (int i = 0; i < N; i++){
			freeFly(particles + i);
			particles[i].vx += (particles[i].y - halfLy)/(particles[i].m)*dtnoise*strainRate*0.1;
			particles[i].coll++;
		}
		for (int j = 0; j < N; j++){
			removeEventFromQueue(eventList[j]);
			crossingEvent(j);

			removeEventFromQueue(eventList[N + j]);
			collisionEvent(j);
		}
	}
	else{
		for (int i = 0; i < N; i++){

			particle* p = particles + i; //(int)(genrand_int32()%N)
			int j = p->num;
			p->coll++;
			//updates position
			freeFly(p);
			//adds the kick
			if (noise == 1){
				randomGaussian(p);
			}
			else{
				p->vx /= sqrt(E/N/T);
				p->vy /= sqrt(E/N/T);
			}
			//p->vx += 0.03*sin(k*p->x - 0.1*t)*dtnoise;
			//p->vy += 0.3*sin(k*p->y - 0.1*t)*dtnoise;

			//Calculate the new collisions and the new cell crossings.
			removeEventFromQueue(eventList[j]);
			crossingEvent(j);

			removeEventFromQueue(eventList[N + j]);
			collisionEvent(j);
		}
	}
	addEventNoise(t + dtnoise);
}


void addEventUpdate(double tupdate){
	node* toAdd;
	toAdd = eventList[2*N + 5];
	toAdd->type = UPDATE;
	toAdd->t = tupdate;
	toAdd->j = 0; //whatever.

	addEventToQueue(toAdd);
}



void updateT(){
	physicalQ();
	T = T*exp(expE - E/N);
	if (t + updateTime < updateTimeMax)
		addEventUpdate(t + updateTime);
}



/* -------------------------------/
/								  /
/     Miscellaneous functions  	  /
/								  /
/--------------------------------*/


void freeFlyNormal(particle* p){
    double dt = t - p->t;
    p->t = t;
	if (addField == 1){
		p->x += dt*p->vx;
		p->y += dt*p->vy + 0.5*field*dt*dt;
		p->vy += dt*field;
	}
	else if (damping == 1){
		double ex = exp(-gamm*dt);
		p->x += p->vx/gamm*(1 - ex);
		p->y += p->vy/gamm*(1 - ex);
		p->vx *= ex;
		p->vy *= ex;
	}
	else{
		p->x += dt*p->vx;
		p->y += dt*p->vy;
		
	}

	PBCpostX(&(p->x));
	PBCpostY(&(p->y));
}

void freeFlyGrow(particle* p){
    double dt = t - p->t;
    p->t = t;

	p->x += dt*p->vx;
	p->y += dt*p->vy;
	p->rad += dt*p->vr;

	PBCpostX(&(p->x));
	PBCpostY(&(p->y));
}

void saveTXT(){
	if (reduce){
		fprintf(fichier, "ITEM: TIMESTEP\n%lf\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n0 %lf\n0 %lf\n0 0\nITEM: ATOMS id x y vx vy\n", t, N, Lx, Ly);
		for(int i = 0; i < N; i++){
			fprintf(fichier, "%d %.2lf %.2lf %.6lf %.6lf\n", i, particles[i].x, particles[i].y, particles[i].vx, particles[i].vy);
		}
	}
	else{
	fprintf(fichier, "ITEM: TIMESTEP\n%lf\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n0 %lf\n0 %lf\n0 0\nITEM: ATOMS id type x y vx vy radius m coll\n", t, N, Lx, Ly);
		for(int i = 0; i < N; i++){
			
            int synchro = 0;
            if (t - particles[i].lastColl > ts)
                synchro = 1;
			fprintf(fichier, "%d %d %.3lf %lf %lf %lf %lf %lf %d\n", i, particles[i].type, particles[i].x, particles[i].y, particles[i].vx, particles[i].vy, particles[i].rad, particles[i].m, synchro);
		}
	}
	fflush(fichier);
}

void initializeThermo(){
	if ((addWallx) || (addWally) || (addCircularWall)){

		if (addMidWall){
			fprintf(thermo, "t E pm p\n");
		}
		else{
			fprintf(thermo, "t E p pW\n");
		}
	}
	else{
		if (addWell){
			fprintf(thermo, "t E Ep p\n");
		}
		else{
			fprintf(thermo, "t E p\n");
		}
	}
}


void saveThermo(){

	double area = Lx*Ly;
	double deltaTime = t - lastScreen;
	lastScreen = t;
	if (addCircularWall){
		area = M_PI*halfLx*halfLx; 
	}
	
    pressure = (-1/(deltaTime))*collTerm/(2*area) + E/area;
	collTerm = 0;


	if (t != 0){
		fprintf(thermo, "%lf ", t);
		if (addCollisionNumberThermo){
			fprintf(thermo, "%ld ", ncol);
		}
		if (addEnergyThermo){
			fprintf(thermo, "%lf ", E/N);
		}
		if (addPressureThermo){
			fprintf(thermo, "%lf ", pressure);
		}
		if (addCircularWall || addWallx || addWally){
			double perimeter = 0;
			if (addCircularWall){
				perimeter = M_PI*Lx;
			}
			if (addWallx){
				perimeter += 2*Ly;
			}
			if (addWally){
				perimeter += 2*Lx;
			}
			fprintf(thermo, "%lf ", dp/(deltaTime*perimeter));
			dp = 0;
		}
		if (addMidWall){
			fprintf(thermo, "%lf %lf", dpMid/(deltaTime*Ly),(dpRight - dpLeft)/(deltaTime*Ly)) ;
			dpMid = 0;
			dpRight = 0;
			dpLeft = 0;
		}
		if (addWell){
			Ep = 0;
			for (int i = 0; i < N; i ++){
				Ep -= U*particles[i].numberOfParticlesInWell/2;
			}
			fprintf(thermo, "%lf ", Ep/N);
		}
		if (msd){
			fprintf(thermo, "%lf ", MSD/N);
		}
		if (addEnergyCenterOfMassThermo){
			double ex = 0;
			double ey = 0;
			for (int i = 0; i < N; i++){
				ex += particles[i].m*particles[i].vx;
				ey += particles[i].m*particles[i].vy;
			}
			fprintf(thermo, "%lf ", 0.5*(ex*ex + ey*ey)/N);
		}
		fprintf(thermo, "\n");
	}
	fflush(thermo);

}

double drand(double min, double max){
    return (genrand()*(max - min)) + min;
}

double sign(double x){
	return (x < 0) ? -1 : (x > 0);
}



inline double logTime(double time){
	if (damping == 1){
		double var = (1 - time*gamm);
		if (var < 0)
			return never;
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
		}
	}
	else{

		double dvT = 0;
		double dvM = 0;
		for (int i = 0; i < N; i ++){
			particle* p = particles + i;
			double r = sqrt((p->x - halfLx)*(p->x - halfLx) + (p->y - halfLx)*(p->y - halfLx));
			dvT -= p->m*(-p->vx*(p->y - halfLx) + p->vy*(p->x - halfLx));
			dvM += p->m*r;
		}

		double dv = dvT/dvM;

		for(int i = 0; i < N; i++){
			particle* p = particles + i;
			double r = sqrt((p->x - halfLx)*(p->x - halfLx) + (p->y - halfLx)*(p->y - halfLx));
			p->vx += -dv*(p->y - halfLx)/r;
			p->vy += dv*(p->x - halfLx)/r;
		}
	}


	physicalQ();
	for (int i = 0; i < N; i++){
		particles[i].vx /= sqrt(E/N/Einit);
		particles[i].vy /= sqrt(E/N/Einit);
	}
}

void optimizeGrowConstant(){
	double baseGrow;
	double criticalPhi = 0.8;
	if (addCircularWall){
		baseGrow = 0.03;
	}
	else{
		baseGrow = 0.1;
	}
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
}


void randomGaussian(particle* p){
	double u1 = genrand_real3();
	double u2 = genrand_real3();
	double a = sqrt(-2*log(u1));
	double b = 2*M_PI*u2;

	if (damping == 0){
		if (euler){
			double std = sqrt(2*gamm*T*dtnoise)/p->m;
			p->vx += std*a*cos(b) - sign(p->vx)*gamm*dtnoise/p->m;
			p->vy += std*a*sin(b) - sign(p->vy)*gamm*dtnoise/p->m;
		}
		else{
			double c = exp(-gamm*dtnoise);
			double std = sqrt(T*(1 - c*c)/p->m);
			p->vx = std*a*cos(b) + (p->vx)*c;
			p->vy = std*a*sin(b) + (p->vy)*c;
		}
	}
	else{
		double std = sqrt(2*gamm*T*dtnoise/p->m);
		p->vx += std*a*cos(b);
		p->vy += std*a*sin(b);
	}
}

void shearCorrection(double* x, double* vx, int sign){
	*x += sign*(fmod(strainRate*(t - 1/vr)*Ly, Lx));
	*vx += sign*strainRate*Ly;
	PBCpostX(x);
}

void correctDistances(particle* p1, particle* p2, double lat2, double* dx, double* dvx){
	if (p1->cell[1] == 0 && p2->cell[1] == Nycells - 1){
			double x = p2->x + lat2*p2->vx;
			double vx = p2->vx;
			shearCorrection(&x, &vx, 1);
			*dx = x - p1->x;
			*dvx = vx - p1->vx;
		}
	if (p1->cell[1] == Nycells - 1 && p2->cell[1] == 0){
		double x = p2->x + lat2*p2->vx;
		double vx = p2->vx;
		shearCorrection(&x, &vx, -1);
		*dx = x - p1->x;
		*dvx = vx - p1->vx;
	}
}

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


void physicalQ(){
	px = 0;
	py = 0;
	E = 0;
	active = 0;
	for (int i = 0; i < N; i++){
		px += particles[i].vx*particles[i].m;
		py += particles[i].vy*particles[i].m;
		E += 0.5*(pow(particles[i].vx, 2)*particles[i].m + pow(particles[i].vy, 2)*particles[i].m);
		if (particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy > 1e-6){
			active += 1.;
		}
	}
	active /= N;
	if (msd != 0 && t >= firstThermo){
		MSD = 0;
		double COMt_x = 0;
		double COMt_y = 0;
		for (int i = 0; i < N; i++){
			COMt_x += particles[i].x + particles[i].crossX*Lx;
			COMt_y += particles[i].y + particles[i].crossY*Ly;
		}
		COMt_x = COMt_x/N;
		COMt_y = COMt_y/N;
		for (int i = 0; i < N; i++){
			MSD += pow(particles[i].x + particles[i].crossX*Lx - posxInitial[i] - (COMt_x - COM_x), 2) + pow(particles[i].y + particles[i].crossY*Ly - posyInitial[i] - (COMt_y - COM_y), 2);
		}
	}
	
}

void customName(){
	
	#if defined(_WIN32)
	_mkdir("dump/");
	#else 
	mkdir("dump/", 0777);
	#endif
	int v = 1;
	sprintf(fileName, "dump/N_%ddtnoise_%.3lfres_%.3lfgamma_%.3lfT_%.3lfphi_%.6lfrat_%.3lfvo_%.3lfao_%.3lfdelta_%.3lfLx_%.3lfLy_%.3lfq_%.3lfv_%d.dump", N, dtnoise, res, gamm, T, phi, sizeratio, ts, deltaE, delta, Lx, Ly, (double)Nsmall/N, v);
	while (access(fileName, F_OK) == 0){
		v += 1;
			sprintf(fileName, "dump/N_%ddtnoise_%.3lfres_%.3lfgamma_%.3lfT_%.3lfphi_%.6lfrat_%.3lfvo_%.3lfao_%.3lfdelta_%.3lfLx_%.3lfLy_%.3lfq_%.3lfv_%d.dump", N, dtnoise, res, gamm, T, phi, sizeratio, ts, deltaE, delta, Lx, Ly, (double)Nsmall/N, v);
	}
    snprintf(thermoName, sizeof(fileName), "dump/N_%ddtnoise_%.3lfres_%.3lfgamma_%.3lfT_%.3lfphi_%.6lfrat_%.3lfvo_%.3lfao_%.3lfdelta_%.3lfLx_%.3lfLy_%.3lfq_%.3lfv_%d.thermo", N, dtnoise, res, gamm, T, phi, sizeratio, ts, deltaE, delta, Lx, Ly, (double)Nsmall/N, v);
}

int mygetline(char* str, FILE* f){
  int comment = 1;
  while (comment)
  {
    if (!fgets(str, 255, f)) return -1;
    if (str[0] != '#') comment = 0;
  }
  return 0;
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
	if ((addWallx) || (addWally)){
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
	if ((res < 1) || (addExpo))
		printf("ON");
	else
		printf("OFF");
	CWHITE;
	printf("\nCOLLISION ENERGY INPUT: ");
	CGREEN;
	if ((addDelta) && (addExpo))
		printf("Δ AND e^");
	else if (addExpo)
		printf("e^");
	else if (addDelta)
		printf("Δ");
	else
		printf("NONE");
	CWHITE;
	printf("\nBOX SHAPE: ");
	CGREEN;
	if (Lx == Ly)
		printf("SQUARE");
	else
		printf("RECTANGULAR");
	CWHITE;
	printf("\n┌───────────────────────────────────────────────────────────────────────────────────────────────────────┐\n│\t N = \033[1;32m%d\033[0;37m Δt = \033[1;32m%.3lf\033[0;37m α = \033[1;32m%.3lf\033[0;37m γ = \033[1;32m%.3lf\033[0;37m T = \033[1;32m%.3lf\033[0;37m Φ = \033[1;32m%.3lf\033[0;37m rat = \033[1;32m%.3lf\033[0;37m L = \033[1;32m%.3lf\033[0;37m\t\t│ \n├────────────────────────────────────────────────┬──────────────────────────────────────────────────────┤\n", N, dtnoise, res, gamm, T, phi, sizeratio, Lx);
}

void printClose(){
	printf("\n");
	printf("└────────────────────────────────────────────────┴──────────────────────────────────────────────────────┘\n");
}

void runningCheck(){
	double stop = 0;
	if (res > 1){
		printf("ERROR:\033[1;31m Coefficient of restitution greater than 1!\033[0m\n");
		stop = 1;
	}
	if (( ( ((addWally) || (addWallx) || (addCircularWall) || (addMidWall)) && (resW < 1) ) || (res < 1) || (damping == 1)) && (((addDelta == 0) && (addDoubleDelta == 0) && (addEvolvingDelta == 0)) && (addExpo == 0) && (noise == 0))){
		printf("ERROR:\033[1;31m Dissipative system without energy input!\033[0m\n");
		stop = 1;
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
	if ((addShear) && (charged || addField || addWell || addWallx || addWally || addCircularWall || damping)){
		printf("ERROR:\033[1;31m Can't use addShear == 1 with charged/Field/Well/Walls/damping == 1\033[0m\n");
		stop = 1;
	}

	printf("\n");

	if ((ERROR) && (stop))
		exit(3);
}


