#include"EDMD.h"
#include"mersenne.c"
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#define CELLCROSS 0
#define COLLISION 1
#define SCREENSHOT 2
#define THERMO 6
#define ADDINGNOISE 3
#define WALL 4
#define UPDATE 5
#define OUT 8
#define IN 7



#define GREEN printf("\033[1;32m")
#define WHITE printf("\033[0;37m")
#define CYAN printf("\033[1;96m")

#define ERROR 1


int N, Nbig, Nsmall;
double px, py, E, pressure, dxr, dyr, dist, active;
double lastScreen = 0;
double lastCollNum = 0;



int Nycells = -1;
int Nxcells = -1;
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

double on = 1;
double t = 0;

double collTerm = 0;
double dp = 0;
double dpLeft = 0;
double dpRight = 0;
double dpMid = 0;

unsigned long int ncol = 0;
unsigned long int ncolss = 0;
unsigned long int ncolsb = 0;
unsigned long int ncolbb = 0;
unsigned long int ncross = 0;
char fileName[350];
char thermoName[350];
char buffer[255];
FILE* file;


//load a configuration if == 1 (if nothing specified it loads from data.txt)
const int load = 1;

const int tinyOnes = 0;

//will be overwritten if:    load = 1     below.
double Lx = 12000;
double Ly = 12000;

//spacing bewteen big particles in the init state if no loading data given
double dL = 2;

//values of radius and mass if load == 0 (so they will be overwritten if we load from a file)
double rad1 = 1.;
double rad2 = 0.457;
double m1 = 1;
double m2 = 0.06;

//duration of simulation
double tmax = 1000000;
//time between each screenshots
double dtime = 10;
double dtimeThermo = 10;
double firstScreen = 0;

//if -1, screenshot will be taken at constant interval of dtimeThermo
double nextScreen = -1;





//add a noise
const int noise = 0;

//update the thermostat temperature to reach a wanted temperature for the system
const int updating = 0;

//add damping
const int damping = 0;

//activate the delta model
const int addDelta = 0;

//activate double delta model
const int addDoubleDelta = 0;

//activate evoling delta model
const int addEvolvingDelta = 0;



//activate the exponential model
const int addExpo = 0;

//add a wall at y = Ly and y = 0
const int addWally = 0;
//add a wall at x = Lx and x = 0
const int addWallx = 0;
//add a wall at x = Lx/2
const int addMidWall = 0;
//add a circular wall
const int addCircularWall = 0;

//add square potential
const int addWell = 0;

//if noise use euler solver
const int euler = 0;

//save a file with thermodynamics quantities each dtimeThermo
const int ther = 1;

//reduce the number of properties dumped into the dump files
const int reduce = 0;


//value for delta model
double delta = 0.03;

//values for double delta model
double deltaM = 0.05;
double deltam = 0;
double ts = 6;

//value for evolving delta model
double tau = 5;
double deltaMax = 0.02;

//values for exponential model
double vo = 1;
double ao = 1.3;

//values for square potential model
double sig = 1.6;
double U = 10.56;

//Initial temperature
double Einit = 0.1;

//coeff of restitution of the wall
double resW = 1;

//coeff of restitution of particles
double res = 1;

//parameter if noise or damping
double gamm = 0.1;
double T = 0.056;
double expE = 1;
//time between kicks
double dtnoise = 0.3;


//parameter to reach a given temperature out of equilibrium with langevin dynamics
double updateTime = 100;
double updateTimeMax = 2000;



double sizeRat;
double phi;



double paulTime = 0;
int actualPaulList = 0;
double dtPaul;
int paulListN;

double t1 = 0;
FILE *fichier;
FILE *thermo;
node *root;
//array containing the particles
particle* particles;
//2D array of pointer
particle*** cellList;
node** eventList;
node** eventPaul;
//variable containing the nextEvent
node* nextEvent;


int main(int argc, char *argv[]){

	init_genrand(time(NULL));
	//init_genrand(666);

	constantInit(argc, argv);
	runningCheck();
	particlesInit();
	fichier = fopen(fileName, "w");
	cellListInit();
	if (addWell)
		pairsInit();
	eventListInit();
	physicalQ();
	format();

	if ((ther) || (addWallx) || (addWally))
		thermo = fopen(thermoName, "w");
	
	while (t <= tmax){
		nextEvent = findNextEvent();
		t = nextEvent->t;

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
   				fflush(stdout);
				if (noise == 0 && E < 1e-10 && addWell == 0){
					goto exiting;
				}
				break;
			case THERMO:
				takeAThermo();
				fflush(stdout);
				if (noise == 0 && E < 1e-10 && addWell == 0){
					saveTXT();
					goto exiting;
				}
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
		}
	}

	exiting:
		printInfo();
		printClose();
		fclose(fichier);
		if (ther)
			fclose(thermo);
		freeArrays();
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
void constantInit(int argc, char *argv[]){
	if (load){
		char filename[200];
		if (argc == 5){
		    double temp1, temp2, temp3;
            sscanf(argv[1], "%lf", &temp1);
            sscanf(argv[2], "%lf", &temp2);
            sscanf(argv[3], "%lf", &temp3);
            char filename[200];
            sprintf(filename, "%.7lf%.7lf%.7lf%.7lf.txt",temp1, temp2, temp3, 0.);
		}
		else if (argc == 4){
		    double temp1, temp2, temp3;
            sscanf(argv[1], "%lf", &temp1);
            sscanf(argv[2], "%lf", &temp2);
            sscanf(argv[3], "%lf", &temp3);
            char filename[200];
            sprintf(filename, "%.7lf%.7lf%.7lf.txt",temp1, temp2, temp3);
		}
		else{
			sprintf(filename, "data.txt");
			if (argc == 2)
				sscanf(argv[1], "%lf", &res);

			if (argc == 3){
				sscanf(argv[1], "%lf", &res);
				sscanf(argv[2], "%lf", &deltaM);
			}
		}
		printf("%s\n", fileName);
		file = fopen(filename, "r");

		mygetline(buffer, file);
		sscanf(buffer, "%d %lf %lf\n", &N, &Lx, &Ly);
		if (Ly == 12000)
			Ly = Lx;
	}
	else{
	//brute force (ie stupid) computation of N

		if (argc == 2){
			sscanf(argv[1], "%lf", &dL);
		}

		int i = 0;
		double y = 0.1 + rad1;
		while (y < Ly - rad1){
			for (double x = rad1; x < Lx - rad1; x += dL + 2*rad1){
				i++;
			}
			y += dL + 2*rad1;
		}

		Nbig = i;

		if (tinyOnes){
			y = 0.5*dL + 2*rad1;

			while (y < Ly - rad2){
				for (double x = 0.5*dL + 2*rad1; x < Lx - rad2; x += dL + 2*rad1){
					i++;
				}
				y += dL + 2*rad1;
			}
		}
		N = i;
		Nsmall = N - Nbig;
	}

	if (nextScreen == -1)
		nextScreen = dtime;
	if (firstScreen == -1)
		firstScreen = tmax - 1;
	halfLx = Lx/2;
	halfLy = Ly/2;
	if ((Nxcells == -1) || (Nycells == -1)){
		if (addWell){
			Nycells = (int)(halfLy/sig);
			Nxcells = (int)(halfLx/sig);
		}
		else{
			Nycells = (int)(halfLy);
			Nxcells = (int)(halfLx);
		}
	}
	cellxSize = Lx/Nxcells;
	cellySize = Ly/Nycells;
	halfxC = cellxSize/2;
	halfyC = cellxSize/2;
	cellxFac = 1/cellxSize;
	cellyFac = 1/cellySize;

	dtPaul = 1000/(double)N;

    paulListN = N;



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
	if (load){
		particle* p;
		for (i = 0; i < N; i++){
		    p = particles + i;
		    mygetline(buffer, file);
		    sscanf(buffer, "%d %lf %lf %lf %lf\n", &(p->type), &(p->x), &(p->y), &(p->rad), &(p->m));
		    Nbig += p->type;
			double U1 = drand(0, 1);
			double U2 = drand(0, 1);
            p->vx = sqrt(-2*log(U1)/p->m*Einit)*cos(2*M_PI*U2) ;
            p->vy = sqrt(-2*log(U1)/p->m*Einit)*sin(2*M_PI*U2) ;
			p->num = i;
			p->t = 0;
			p->coll = 0;
            p->lastColl = 0;
		}
		fclose(file);
		Nsmall = N - Nbig;
	}
	else{
		double vxTot, vyTot;
		//big particles
		double y = 0.1 + rad1;
		while (y < Ly - rad1){
			for (double x = rad1; x < Lx - rad1; x += dL + 2*rad1){
				particles[i].x = x;
				particles[i].y = y;
				particles[i].vx = drand(-2, 2);
				particles[i].vy = drand(-2, 2);
				vxTot += particles[i].vx;
				vyTot += particles[i].vy;
				particles[i].num = i;
				particles[i].rad = rad1;
				particles[i].m = m1;
				particles[i].type = 1;
				particles[i].t = 0;
				particles[i].coll = 0;
				i++;
			}
			y+= dL + 2*rad1;
		}

		if (tinyOnes){
		//tiny ones surrounded by big ones
			y = 0.5*dL + 2*rad1;
			while (y < Ly - rad2){
				for (double x = 0.5*dL + 2*rad1; x < Lx - rad2; x += dL + 2*rad1){
					particles[i].x = x;
					particles[i].y = y;
					particles[i].vx = drand(-2,2);
					particles[i].vy = drand(-2,2);
					vxTot += particles[i].vx;
					vyTot += particles[i].vy;
					particles[i].num = i;
					particles[i].rad = rad2;
					particles[i].m = m2;
					particles[i].type = 0;
					particles[i].t = 0;
					particles[i].coll = 0;
					i++;
				}
				y += dL + 2*rad1;
			}
		}
	}

	sizeRat = particles[N-1].rad;
	if (addCircularWall)
		phi = (sizeRat*sizeRat*Nsmall + Nbig)/(halfLx*halfLx);
	else
		phi = M_PI*(sizeRat*sizeRat*Nsmall + Nbig)/(Lx*Ly);
	customName();



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


/* ---------------------------------------------
	Initializes the list of particles close!
--------------------------------------------- */

void pairsInit(){

	
	printf("Generating pairs for well potential...");

	fflush(stdout);

	

	particle* p1;

	for(int i = 0; i < N; i++){
		p1 = particles + i;
		p1->particlesInWell = calloc(20, sizeof(unsigned int));
		p1->numberOfParticlesInWell = 0;
		int X = p1->cell[0];
		int Y = p1->cell[1];

		for (int j = -1; j <= 1; j++){
			for (int k = -1; k <= 1; k++){
				particle* p2 = cellList[PBCcell(X + j, 1)][PBCcell(Y + k, 0)];
				while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (p1->num != p2->num){
						double dx = p2->x - p1->x;
						double dy = p2->y - p1->y;
						PBC(&dx, &dy);
						if (sqrt(dx*dx + dy*dy) < sig*(p1->rad + p2->rad))
							addParticleInWellList(p1, p2->num);
					}
					p2 = p2->nxt;
				}
			}
		}
	}

	printf("Done!\n");

	fflush(stdout);
}

/* ---------------------------------------------
	Initializes the cellList. cellList[X][Y]
	points toward a single particle contained
	in the cell {X, Y}, all the others particles
	in the same cell are reachable by a
	double linked list containing cellList[X][Y]
--------------------------------------------- */
void cellListInit(){

	cellList = calloc(Nxcells, sizeof(particle**));
	for (int i = 0; i < Nxcells; i++){
		cellList[i] = calloc(Nycells, sizeof(particle*));
	}


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
	root->t = 1000000001; //arbitrary value, root->rgt will always be NULL.
	root->rgt = NULL;
	root->lft = NULL;
	root->top = NULL;


	for (int i = 0; i < N; i++){
		eventList[i]->i = i;
		eventList[N + i]->i = i;
		crossingEvent(i);
		collisionEvent(i);
	}

	addEventScreenshot(firstScreen);
	addEventThermo(dtimeThermo);

	if (noise)
		addEventNoise(dtnoise);
	if (updating)
		addEventUpdate(updateTime);
}


void freeArrays(){
	int totalSize = 2*N + suppSizeEvent;
	if (addWell){
		for (int i = 0; i < N; i++){
			free(particles[i].particlesInWell);
		}	
	}
	free(particles);
	for (int i = 0; i < Nxcells; i++){
			free(cellList[i]);
	}
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
	p->nxt = cellList[X][Y];
	cellList[X][Y] = p;
	if (p->nxt)
		p->nxt->prv = p;

}

/* ----------------------------------------
   	Removes Particle i from the CellList
-----------------------------------------*/
void removeFromCell(int i){
	particle* p = particles + i;
	if (p->prv == NULL)
		cellList[p->cell[0]][p->cell[1]] = p->nxt;
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


int PBCcell(double a, int x){
	if (x){
		if (a < 0)
			return a + Nxcells;
		else if (a >= Nxcells)
			return a - Nxcells;
		return a;
	}
	else{
		if (a < 0)
			return a + Nycells;
		else if (a >= Nycells)
			return a - Nycells;
		return a;
	}
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



void crossingEvent(int i){

	particle p = particles[i];
	double travelTimeX = 10000000;
	double travelTimeY = 10000000;

	int xx;
	int yy = 1;

	// xx and yy indicates the direction to the new cell for example xx = 1 means that we will go from cell Cell[X][Y] to Cell[X - 1][Y]
	if (p.vx < 0){
		travelTimeX = logTime(PBCinsideCell(p.cell[0]*cellxSize - p.x, 1)/p.vx);
		xx = 1;
	}
	else{
		travelTimeX = logTime(PBCinsideCell((1 + p.cell[0])*cellxSize - p.x, 1)/p.vx);
		xx = 2;
	}
	if (p.vy < 0){
		travelTimeY = logTime(PBCinsideCell(p.cell[1]*cellySize - p.y, 0)/p.vy);
		yy = 3;
	}
	else{
		travelTimeY = logTime(PBCinsideCell((1 + p.cell[1])*cellySize - p.y, 0)/p.vy);
		yy = 4;
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
			}
			else
				p->cell[1] = newCell;
			break;
		case 4:
			newCell = p->cell[1] + 1;
			if (newCell == Nycells){
				p->cell[1] = 0;
				p->crossY += 1;
			}
			else
				p->cell[1] = newCell;
			break;
	}

	p->prv = NULL;
	p->nxt = cellList[p->cell[0]][p->cell[1]];
	cellList[p->cell[0]][p->cell[1]] = p;
	if (p->nxt != NULL)
		p->nxt->prv = p;


	crossingEvent(i);
	removeEventFromQueue(eventList[N + i]);
	collisionEvent(i);
}



/* --------------------------------------
	Calculates the time of collision
	of p1 and p2
--------------------------------------- */
double collisionTime(particle* p1, particle* p2){

	//return 10000000; //to deactivate collisions

	if (damping == 1){
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
		return 100000000;


	double v2 = dvx*dvx + dvy*dvy;
	double distOfSquare = 4*p1->rad*p2->rad;
	distOfSquare = dx*dx + dy*dy - 4*p1->rad*p2->rad;
	double det = pow(b, 2) - v2*distOfSquare;


	//to delete when confident with life decisions...
	if (distOfSquare < -0.01){
		printf("\nERROR:\033[0;31m Overlaps detected!\033[0m\n");
		exit(3);
	}

	if (det < 0)
		return 100000000;

	double timeOfCollision = (-b - sqrt(det))/v2;

	return logTime(timeOfCollision);
}

//time to enter in the square well: https://faculty.biu.ac.il/~rapaport/papers/09b-ptp.pdf
double sphereTime(particle* p1, particle* p2){

	if (damping == 1){
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
		return 100000000;


	double v2 = dvx*dvx + dvy*dvy;
	double distOfSquare = dx*dx + dy*dy -  sig*(p1->rad + p2->rad)*sig*(p1->rad + p2->rad);
	double det = b*b - v2*distOfSquare;


	if (det < 0)
		return 100000000;

	return logTime(- b - sqrt(det))/v2;
}


//time to leave the square well: https://faculty.biu.ac.il/~rapaport/papers/09b-ptp.pdf
double separateTime(particle* p1, particle* p2){

	if (damping == 1){
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
		return 10000000000;

	double temp = (- b + sqrt(det))/v2;
	if (temp < 0)
		return 10000000000;
	else
		return logTime(temp);
}
/* ---------------------------------------------
	Finds the next collision of particle i
	by calculating collision time of i with
	its surrouding (using the cellList)
--------------------------------------------- */
void collisionEvent(int i){
	int finalPartner = 0;
	double dt = 10000000;
	int type = COLLISION;
	particle* p1 = particles + i;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	int xy = 2;
	double dtTemp, dtTemp2;
	int typeTemp;


	if (addWally){

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
		
			double y0 = p1->y;
			double x0 = p1->x;
			double v0x = p1->vx;
			double v0y = p1->vy;
			double rad = p1->rad;
			
			x0 = x0 - halfLx;
			y0 = y0 - halfLx;
			double A = (v0x*v0x + v0y*v0y);
			double B = 2*(v0x*x0 + v0y*y0);
			double C = (x0*x0 + y0*y0) - (halfLx - rad)*(halfLx - rad);
			dt  = (-B + sqrt(B*B - 4*A*C))/(2*A);
			
			type = WALL;
			xy = 2;
			

	}

	for (int j = -1; j <= 1; j++){
		for (int k = -1; k <= 1; k++){
			particle* p2 = cellList[PBCcell(X + j, 1)][PBCcell(Y + k, 0)];
			while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
				if (p1->num != p2->num){
					if (addWell){
						if (isParticleInWellList(p1, p2->num)){
							dtTemp = collisionTime(p1, p2);
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
						dtTemp = collisionTime(p1, p2);
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

	if (type == COLLISION){
		addCollisionEvent(i, finalPartner, t + dt);
	}
	else if (type == OUT)
		addOutEvent(i, finalPartner, t + dt);
	else if (type == IN)
		addInEvent(i, finalPartner, t + dt);
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

void doTheWall(){
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
		pi->vx = -resW*pi->vx;
	}
	else if (xy == 1){
		dp += fabs((resW + 1)*pi->m*pi->vy);
		pi->vy = -resW*pi->vy;
	}
	else{
		double theta = atan2(pi->y - halfLx, pi->x - halfLx);
		double nx = cos(theta);
		double ny = sin(theta);
		double mix = -(1 + resW)*(pi->vx*nx + pi->vy*ny);

		pi->vx += mix*nx;
		pi->vy += mix*ny;
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
void doTheCollision(){

	int i = nextEvent->i;
	int j = nextEvent->j;

	particle* pi = particles + i;
	particle* pj = particles + j;


	freeFly(pi); //might put it on top?
	if (nextEvent->collActual != pj->coll){ //the collision trick of the article
        collisionEvent(i);
        return;
    }

    if ((pi->type == 1) && (pj->type == 1)){
    	ncolbb += 1;
    }
    else if ((pi->type == 0) && (pj->type == 0)){
    	ncolss += 1;
    }
    else{
    	ncolsb += 1;
    }
    ncol++;
	freeFly(pj);

	pi->coll++;
	pj->coll++;

	double dx = pj->x - pi->x;
	double dy = pj->y - pi->y;

	double dvx = pj->vx - pi->vx;
	double dvy = pj->vy - pi->vy;



	PBC(&dx, &dy);



	if (addExpo)
		res = resCoeff(sqrt(dvx*dvx + dvy*dvy));

	double invMass = 1/(pi->m + pj->m);
	double collKernel = invMass*(1 + res)*(dx*dvx + dy*dvy);



	collTerm += pi->m*pj->m*collKernel;
	if ((addDelta) ||(addDoubleDelta || (addEvolvingDelta))){
		//double tau = -ts*log(1 - drand(0, 1));
		if (addDoubleDelta){
			if ((t - pi->lastColl > ts) && (t - pj->lastColl > ts))
				delta = deltam;
			else
				delta = deltaM;
		}
		else if (addEvolvingDelta){
			double dti = t - pi->lastColl;
			double dtj = t - pj->lastColl;
			delta = deltaMax*(2 - (exp(-dtj/tau) + exp(-dti/tau)))/2;
		}

		dist = sqrt(dx*dx + dy*dy);
		collTerm -= pi->m*pj->m*invMass*2*delta*dist;
		dxr = dx/dist;
		dyr = dy/dist;
	}


	double funkyFactor = collKernel/(4*pi->rad*pj->rad);




	if ((addDoubleDelta) || (addDelta) || (addEvolvingDelta)){

		pi->vx += funkyFactor*pj->m*dx - 2*pj->m*invMass*delta*dxr;
		pi->vy += funkyFactor*pj->m*dy - 2*pj->m*invMass*delta*dyr;
		pj->vx -= funkyFactor*pi->m*dx - 2*pi->m*invMass*delta*dxr;
		pj->vy -= funkyFactor*pi->m*dy - 2*pi->m*invMass*delta*dyr;
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

	if (U > 0){
		addParticleInWellList(pi, j);
		addParticleInWellList(pj, i);
		double atroce = (mi*vi + mj*vj);
		double horrible = sqrt(mi*mi*mj*mj*(vi - vj)*(vi - vj) + 2*mi*mj*(mi + mj)*U);
		vif = (mi*atroce + horrible)/(mi*(mi + mj));
		vjf = (mj*atroce - horrible)/(mj*(mi + mj));
	}
	else{
		if (U > -temporary/(2*(mi + mj))){
			double atroce = (mi*vi + mj*vj);
			double horrible = sqrt(mi*mj*temporary + 2*mi*mj*(mi + mj)*U);
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

	if (U > 0){
		if (U < temporary/(2*(mi + mj))){ //leaves the square well
			double atroce = (mi*vi + mj*vj);
			double horrible = sqrt(mi*mj*temporary - 2*mi*mj*(mi + mj)*U);
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
		double horrible = sqrt(mi*mj*temporary - 2*mi*mj*(mi + mj)*U);
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
	physicalQ();
	saveTXT();
	if (on){
		on = 0;
		addEventScreenshot(t + dtime);
		return;
	}
	on = 1;
	addEventScreenshot(t + nextScreen);
	printInfo();
}

void takeAThermo(){

	if (damping == 1){
		for (int i = 0; i < N; i++){
			freeFly(particles + i);
		}
	}
	physicalQ();
	saveThermo();

	addEventThermo(t + dtimeThermo);
	printInfo();
}

void printInfo(){
	int size = 29;
	int place = (int)(size*t)/(int)tmax;
	printf("\r│\033[1;32m%.2e\033[0;37m/%.2e ", t, tmax);
	CYAN;
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
	WHITE;
	
	/*
	if (addWell){
		for (int i = 0; i < N; i ++){
			for (int j = 0; j < i; j++){
				E -= U*pairs[i*N + j];
			} 
		}
	}
	*/

	printf(" │ n° coll = \033[1;32m%.2e\033[0;37m  Energy = \033[1;32m%.3lf\033[0;37m  Pressure = \033[1;32m%.3lf\033[0;37m │", (float)ncol, E/N, pressure);

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
	for (int i = 0; i < N; i++){
		particle* p = particles + i; //(int)(genrand_int32()%N)
		int j = p->num;
		p->coll++;
		//updates position
		freeFly(p);
		//adds the kick
		randomGaussian(p);


		//Calculate the new collisions and the new cell crossings.
		removeEventFromQueue(eventList[j]);
		crossingEvent(j);

		removeEventFromQueue(eventList[N + j]);
		collisionEvent(j);
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


void freeFly(particle* p){
    double dt = t - p->t;
    p->t = t;
	if (damping == 0){
		p->x += dt*p->vx;
		p->y += dt*p->vy;
	}
	else{
		double ex = exp(-gamm*dt);
		p->x += p->vx/gamm*(1 - ex);
		p->y += p->vy/gamm*(1 - ex);
		p->vx *= ex;
		p->vy *= ex;
	}

	PBCpost(&(p->x), 1);
	PBCpost(&(p->y), 0);
}


void saveTXT(){
	if (reduce){
		fprintf(fichier, "ITEM: TIMESTEP\n%lf\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS xy xz yz\n0 %lf 0\n0 %lf 0\n0 2 0\nITEM: ATOMS id type x y radius\n", t, N, Lx, Ly);
		for(int i = 0; i < N; i++){
			fprintf(fichier, "%d %d %lf %lf %lf\n", i, particles[i].type, particles[i].x, particles[i].y, particles[i].rad);
		}
	}
	else{
	fprintf(fichier, "ITEM: TIMESTEP\n%lf\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS xy xz yz\n0 %lf 0\n0 %lf 0\n0 2 0\nITEM: ATOMS id type x y vx vy radius m\n", t, N, Lx, Ly);
		for(int i = 0; i < N; i++){
			/*
            int synchro = 0;
            if (t - particles[i].lastColl > ts)
                synchro = 1;*/
			fprintf(fichier, "%d %d %lf %lf %lf %lf %lf %lf\n", i, particles[i].type, particles[i].x, particles[i].y, particles[i].vx, particles[i].vy, particles[i].rad, particles[i].m);
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
		if (ther){
			double deltaColl = ncol - lastCollNum;
			lastCollNum = ncol;
			if ((addWallx) || (addWally) || (addCircularWall)){

				if (addMidWall){

					fprintf(thermo, "%lf %lf %lf %lf\n", t, E/N, dpMid/((t - t1)*Ly), (dpRight - dpLeft)/((t - t1)*Ly));
					dpMid = 0;
					dpRight = 0;
					dpLeft = 0;
				}
				else{
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
					fprintf(thermo, "%lf %lf %lf %lf %lf %lf %lf %lf\n", t, ncolss/deltaTime, ncolbb/deltaTime, ncolsb/deltaTime, deltaColl/deltaTime, E/N, pressure, dp/((t - t1)*perimeter));

					dp = 0;
				}
			t1 = t;
			}
			else{
				fprintf(thermo, "%lf %lf %lf %lf %lf %.10lf %lf %lf\n", t, ncolss/deltaTime, ncolbb/deltaTime, ncolsb/deltaTime, deltaColl/deltaTime, E/N, pressure, active);
			}
			ncolss = 0;
			ncolsb = 0;
			ncolbb = 0;
		}
	}
}

double drand(double min, double max){
    return (genrand()*(max - min)) + min;
}

double sign(double x){
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}



inline double logTime(double time){
	if (damping == 1){
		double var = (1 - time*gamm);
		if (var < 0)
			return 1000000000000;
		else
			return  -log(var)/gamm;
	}
	else
		return time;
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


void PBC(double* dx, double* dy){
	if (*dx >= halfLx)
		*dx -= Lx;
	else if (*dx < -halfLx)
		*dx += Lx;
	if (*dy >= halfLy)
		*dy -= Ly;
	else if (*dy < -halfLy)
		*dy += Ly;
}

double PBCinsideCell(double dx, int x){
	if (x){
		if (dx >= halfLx)
			return dx - Lx;
		else if (dx < -halfLx){
			return dx + Lx;

			}
		return dx;
	}
	else{
		if (dx >= halfLy)
			return dx - Ly;
		else if (dx < -halfLy)
			return dx + Ly;
	return dx;
	}
}

void PBCpost(double* val, int x){
	if (x){
		if (*val < 0)
			*val += Lx;
		else if (*val >= Lx)
			*val -= Lx;
	}
	else{
		if (*val < 0)
			*val += Ly;
		else if (*val >= Ly)
			*val -= Ly;
	}
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
}

void customName(){
	mkdir("dump/", 0777);
	int v = 1;
	sprintf(fileName, "dump/N_%ddtnoise_%.3lfres_%.3lfgamma_%.3lfT_%.3lfphi_%.6lfrat_%.3lfvo_%.3lfao_%.3lfdelta_%.3lfLx_%.3lfLy_%.3lfq_%.3lfdL_%.3lfv_%d.dump", N, dtnoise, res, gamm, T, phi, sizeRat, vo, ao, deltaM, Lx, Ly, (double)Nsmall/N, dL, v);
	while (access(fileName, F_OK) == 0){
		v += 1;
			sprintf(fileName, "dump/N_%ddtnoise_%.3lfres_%.3lfgamma_%.3lfT_%.3lfphi_%.6lfrat_%.3lfvo_%.3lfao_%.3lfdelta_%.3lfLx_%.3lfLy_%.3lfq_%.3lfdL_%.3lfv_%d.dump", N, dtnoise, res, gamm, T, phi, sizeRat, vo, ao, deltaM, Lx, Ly, (double)Nsmall/N, dL, v);
	}
    snprintf(thermoName, sizeof(fileName), "dump/N_%ddtnoise_%.3lfres_%.3lfgamma_%.3lfT_%.3lfphi_%.6lfrat_%.3lfvo_%.3lfao_%.3lfdelta_%.3lfLx_%.3lfLy_%.3lfq_%.3lfdL_%.3lfv_%d.thermo", N, dtnoise, res, gamm, T, phi, sizeRat, vo, ao, deltaM, Lx, Ly, (double)Nsmall/N, dL, v);
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
	GREEN;
	if ((addWallx) || (addWally)){
		if (resW == 1)
			printf("HARD WALLS");
		else
			printf("DISSIPATIVE WALLS");
		}
	else
		printf("PERIODICS");
	WHITE;
	printf("\nTHERMOSTAT: ");
	GREEN;
	if (noise)
		printf("ON");
	else
		printf("OFF");
	WHITE;
	printf("\nDISSIPATION: ");
	GREEN;
	if ((res < 1) || (addExpo))
		printf("ON");
	else
		printf("OFF");
	WHITE;
	printf("\nCOLLISION ENERGY INPUT: ");
	GREEN;
	if ((addDelta) && (addExpo))
		printf("Δ AND e^");
	else if (addExpo)
		printf("e^");
	else if (addDelta)
		printf("Δ");
	else
		printf("NONE");
	WHITE;
	printf("\nBOX SHAPE: ");
	GREEN;
	if (Lx == Ly)
		printf("SQUARE");
	else
		printf("RECTANGULAR");
	WHITE;
	printf("\n┌───────────────────────────────────────────────────────────────────────────────────────────────────────┐\n│\t N = \033[1;32m%d\033[0;37m Δt = \033[1;32m%.3lf\033[0;37m α = \033[1;32m%.3lf\033[0;37m γ = \033[1;32m%.3lf\033[0;37m T = \033[1;32m%.3lf\033[0;37m Φ = \033[1;32m%.3lf\033[0;37m rat = \033[1;32m%.3lf\033[0;37m L = \033[1;32m%.3lf\033[0;37m\t\t│ \n├────────────────────────────────────────────────┬──────────────────────────────────────────────────────┤\n", N, dtnoise, res, gamm, T, phi, sizeRat, Lx);
}

void printClose(){
	printf("\n");
	printf("└────────────────────────────────────────────────┴──────────────────────────────────────────────────────┘\n");
}

void runningCheck(){
	double stop = 0;
	if ((res > 1) && (addDelta == 0)){
		printf("ERROR:\033[1;31m Coefficient of restitution greater than 1!\033[0m\n");
		stop = 1;
	}
	if (( ( ((addWally) || (addWallx) || (addCircularWall)) && (resW < 1) ) || (res < 1)) && (((addDelta == 0) && (addDoubleDelta == 0) && (addEvolvingDelta == 0)) && (addExpo == 0) && (noise == 0))){
		printf("ERROR:\033[1;31m Dissipative system without energy input!\033[0m\n");
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

	printf("\n");

	if ((ERROR) && (stop))
		exit(3);
}


