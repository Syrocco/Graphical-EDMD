#include"EDMD.h"
#include"mersenne.c"

#define CELLCROSS 0
#define COLLISION 1
#define SCREENSHOT 2
#define ADDINGNOISE 3
#define WALL 4
#define UPDATE 5

#define GREEN printf("\033[1;32m")
#define WHITE printf("\033[0;37m")
#define CYAN printf("\033[1;96m")

#define ERROR 1


int N, Nbig, Nsmall;
double px, py, E, pressure, dxr, dyr, dist;
double lastScreen = 0;
double lastCollNum = 0;

//will be overwritten if:    load = 1     below.
double Lx = 12000;
double Ly = 12000;
//spacing bewteen big particles in the init state
double dL = 2;


double t = 0;
double tmax = 30000;
//time between each screenshots
double dtime = 1;

double firstScreen = 29999;
double nextScreen = -1;
double on = 1;

double updateTime = 100;
double updateTimeMax = 2000;
unsigned long int ncol = 0;
unsigned long int ncolss = 0;
unsigned long int ncolsb = 0;
unsigned long int ncolbb = 0;
unsigned long int ncross = 0;
char fileName[350];
char thermoName[350];
char buffer[255];
FILE* file;

double rad1 = 1.;
double rad2 = 0.457;
double m1 = 1;
double m2 = 0.06;

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

const int tinyOnes = 0;
const int noise = 1;
const int updating = 0;

const int addDelta = 1;
const int addExpo = 0;
const int addWall = 0;
const int euler = 0;
const int ther = 1;
const int reduce = 0;
double dp = 0;

double delta = 0;
double deltaM = 0.05;
double deltam = 0;
double ts = 4;
double vo = 1;
double ao = 1.3;

double Einit = 0.3;
double resW = 1;
double res = 0.4;
double gamm = 0.003;
double T = 0;
double expE = 1;
//time between kicks
double dtnoise = 0.1;
double collTerm = 0;

double sizeRat;
double phi;

const int load = 1;

double paulTime = 0;
int actualPaulList = 0;
double dtPaul;
int paulListN;

double t1 = 0;
FILE *fichier;
FILE *thermo;
node* root;
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
	eventListInit();
	physicalQ();
	format();

	if ((ther) || (addWall))
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
				break;
			case UPDATE:
				updateT();
				break;

		}
	}
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
            sscanf(argv[4], "%lf", &deltaM);
            char filename[200];
            sprintf(filename, "%.7lf%.7lf%.7lf.txt",temp1, temp2, temp3);
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
				sscanf(argv[2], "%lf", &tmax);
			}
		}
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
		Nycells = (int)(halfLy);
		Nxcells = (int)(halfLx);
	}
	cellxSize = Lx/Nxcells;
	cellySize = Ly/Nycells;
	halfxC = cellxSize/2;
	halfyC = cellxSize/2;
	cellxFac = 1/cellxSize;
	cellyFac = 1/cellySize;

	dtPaul = 1/(double)N;

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
			p->vx = drand(-1, 1);
			p->vy = drand(-1, 1);
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
	phi = M_PI*(sizeRat*sizeRat*Nsmall + Nbig)/(Lx*Ly);
	customName();
	physicalQ();
	for (int i = 0; i < N; i++){
		particles[i].vx -= px/(N*particles[i].m);
		particles[i].vy -= py/(N*particles[i].m);
	}
	physicalQ();
	for (int i = 0; i < N; i++){
		particles[i].vx /= sqrt(E/N/Einit);
		particles[i].vy /= sqrt(E/N/Einit);
	}
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

	if (noise)
		addEventNoise(dtnoise);
	if (updating)
		addEventUpdate(updateTime);
}


void freeArrays(){
	int totalSize = 2*N + suppSizeEvent;

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
		else if (a >= Nycells )
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

	if (p.vx < 0){
		travelTimeX = PBCinsideCell(p.cell[0]*cellxSize - p.x, 1)/p.vx;
		xx = 1;
	}
	else{
		travelTimeX = PBCinsideCell((1 + p.cell[0])*cellxSize - p.x, 1)/p.vx;
		xx = 2;
	}
	if (p.vy < 0){
		travelTimeY = PBCinsideCell(p.cell[1]*cellySize - p.y, 0)/p.vy;
		yy = 3;
	}
	else{
		travelTimeY = PBCinsideCell((1 + p.cell[1])*cellySize - p.y, 0)/p.vy;
		yy = 4;
	}
	//break ___ if ((travelTimeX < 0) || (travelTimeY < 0))
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

	double lat2 = t - p2->t;

	double dvx = p2->vx - p1->vx;
	double dvy = p2->vy - p1->vy;
	double dx = (p2->x + lat2*p2->vx) - p1->x;
	double dy = (p2->y + lat2*p2->vy) - p1->y;
	PBC(&dx, &dy);
	double b = dx*dvx + dy*dvy;

	//no collision.
	if (b > 0)
		return 100000;


	double v2 = dvx*dvx + dvy*dvy;
	double distOfSquare = 4*p1->rad*p2->rad;
	distOfSquare = dx*dx + dy*dy - 4*p1->rad*p2->rad;
	double det = pow(b, 2) - v2*distOfSquare;



	//to delete when confident with life decisions...
	if (distOfSquare < -0.0001){
		printf("ERROR:\033[0;31m Overlaps detected!\033[0m\n");
		exit(3);

	}

	if (det < 0)
		return 100000;

	return (- b - sqrt(det))/v2;
}


/* ---------------------------------------------
	Finds the next collision of particle i
	by calculating collision time of i with
	its surrouding (using the cellList)
--------------------------------------------- */
void collisionEvent(int i){
	int finalPartner = 0;
	double dt = 1000000;
	int type = COLLISION;
	particle* p1 = particles + i;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	int xy = 2;
	double dtTemp;

	if (addWall){
		if ((Y == 0) && (p1->vy < 0)){
			dt = -(p1->y - p1->rad)/p1->vy;
			type = WALL;
			xy = 1;
		}
		else if ((Y == Nycells - 1) && (p1->vy > 0)){

			dt = (Ly - p1->y - p1->rad)/p1->vy;
			type = WALL;
			xy = 1;
		}


		if ((X == 0) && (p1->vx < 0)){

			dtTemp = -(p1->x - p1->rad)/p1->vx;
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
			}
			xy = 0;
		}
		else if ((X == Nxcells - 1) && (p1->vx > 0)){

			dtTemp = (Lx - p1->x - p1->rad)/p1->vx;
			if (dt > dtTemp){
				type = WALL;
				dt = dtTemp;
			}
			xy = 0;
		}
	}

	for (int j = -1; j <= 1; j++){
		for (int k = -1; k <= 1; k++){
			particle* p2 = cellList[PBCcell(X + j, 1)][PBCcell(Y + k, 0)];
			while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
				if (p1->num != p2->num){
			 		double dtTemp = collisionTime(p1, p2);
			 		if (dt > dtTemp){ //get min
						finalPartner = p2->num;
						dt = dtTemp;
						type = COLLISION;
					}
				}
				p2 = p2->nxt;
	 		}
	 	}
	}

	if (type == COLLISION)
		addCollisionEvent(i, finalPartner, t + dt);
	else
		addWallEvent(i, xy, t + dt);
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

/* ---------------------------------------------
	Creates the collision event and adds it
	to the BST.
--------------------------------------------- */
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
		pi->vx = -resW*pi->vx;
	}
	else{
		dp += fabs((resW + 1)*pi->m*pi->vy);
		pi->vy = -resW*pi->vy;
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
    	//res = 0.9;
    }
    else if ((pi->type == 0) && (pj->type == 0)){
    	ncolss += 1;
    	//res = 0.1;
    }
    else{
    	ncolsb += 1;
 	 	//res = 0.95;
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
	if (addDelta){
		double tau = -ts*log(1 - drand(0, 1));
        if ((t - pi->lastColl > tau) && (t - pj->lastColl > tau))
            delta = deltam;
        else
            delta = deltaM;
		dist = sqrt(dx*dx + dy*dy);
		collTerm -= pi->m*pj->m*invMass*2*delta*dist;
		dxr = dx/dist;
		dyr = dy/dist;

	}



	double funkyFactor = collKernel/(4*pi->rad*pj->rad);




	if (addDelta){

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

void printInfo(){
	int size = 29;
	physicalQ();
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
    p->x += dt*p->vx;
    p->y += dt*p->vy;
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
	fprintf(fichier, "ITEM: TIMESTEP\n%lf\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS xy xz yz\n0 %lf 0\n0 %lf 0\n0 2 0\nITEM: ATOMS id type x y xu yu vx vy radius m z synchro\n", t, N, Lx, Ly);
		for(int i = 0; i < N; i++){
            int synchro = 0;
            if (t - particles[i].lastColl > -ts*log(1 - drand(0, 1)))
                synchro = 1;
			fprintf(fichier, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", i, particles[i].type, particles[i].x, particles[i].y, particles[i].x + particles[i].crossX*Lx, particles[i].y + particles[i].crossY*Ly, particles[i].vx, particles[i].vy, particles[i].rad, particles[i].m, particles[i].rad, synchro);
		}
	}
    physicalQ();

	double deltaTime = t - lastScreen;
	lastScreen = t;
    pressure = (-1/(deltaTime))*collTerm/(2*Lx*Ly) + E/(Lx*Ly);
    collTerm = 0;

    if (t != 0){
		if (ther){
			double deltaColl = ncol - lastCollNum;
			lastCollNum = ncol;
			if (addWall){
				fprintf(thermo, "%lf %lf %lf %lf %lf %lf %lf %lf\n", t, ncolss/deltaTime, ncolbb/deltaTime, ncolsb/deltaTime, deltaColl/deltaTime, E/N, pressure, dp/((t - t1)*2*(Lx + Ly)));
				printf("virial %lf et wall %lf \n", pressure, dp/((t - t1)*2*(Lx + Ly)));
				t1 = t;
				dp = 0;
			}
			else{
				fprintf(thermo, "%lf %lf %lf %lf %lf %lf %lf\n", t, ncolss/deltaTime, ncolbb/deltaTime, ncolsb/deltaTime, deltaColl/deltaTime, E/N, pressure);
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

void randomGaussian(particle* p){
	double u1 = genrand_real3();
	double u2 = genrand_real3();
	double a = sqrt(-2*log(u1));
	double b = 2*M_PI*u2;
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
		else if (dx < -halfLx)
			return dx + Lx;
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
	for (int i = 0; i < N; i++){
		px += particles[i].vx*particles[i].m;
		py += particles[i].vy*particles[i].m;
		E += 0.5*(pow(particles[i].vx, 2)*particles[i].m + pow(particles[i].vy, 2)*particles[i].m);
	}
}

void customName(){
	mkdir("dump/", 0777);
    snprintf(fileName, sizeof(fileName), "dump/N_%ddtnoise_%.3lfres_%.3lfgamma_%.3lfT_%.3lfphi_%.6lfrat_%.3lfvo_%.3lfao_%.3lfdelta_%.3lfLx_%.3lfLy_%.3lfq_%.3lfdL_%.3lf.dump", N, dtnoise, res, gamm, T, phi, sizeRat, vo, ao, deltaM, Lx, Ly, (double)Nsmall/N, dL);
    snprintf(thermoName, sizeof(fileName), "dump/N_%ddtnoise_%.3lfres_%.3lfgamma_%.3lfT_%.3lfphi_%.6lfrat_%.3lfvo_%.3lfao_%.3lfdelta_%.3lfLx_%.3lfLy_%.3lfq_%.3lfdL_%.3lf.thermo", N, dtnoise, res, gamm, T, phi, sizeRat, vo, ao, deltaM, Lx, Ly, (double)Nsmall/N, dL);
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
	if (addWall){
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
	if (( ( (addWall) && (resW < 1) ) || (res < 1)) && ((addDelta == 0) && (addExpo == 0) && (noise == 0))){
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

	printf("\n");

	if ((ERROR) && (stop))
		exit(3);
}


