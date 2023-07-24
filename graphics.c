#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h" 
#include "color.h"
#include "graphics.h"
#include "EDMD.h"
#include <pthread.h>

#define qN 200
#define SIZE 400
#define NUM_THREADS 6

int screenWidth = 1800;
int screenHeight = 900;
double factor;
int start;
double xGUI = 1;
double yGUI = 1;

int leftClicked = 0;
int selected = 0;
particle* particleUnderClick = NULL;
Camera2D cam = {0};

Color particleArray[2] = {GRAY, MAROON};
Color particleColor;
bool colorEditing = false;
int colorParam = 0;
bool colorEditing2 = false;
int colorParam2 = 0;
Color* colorArray = Plasma;
double (*colorFunction)(particle*);
double* colorFunctionArray;

double qx[qN] = {0};
float structFactor[qN][qN] = {0};
bool structFactorActivated = 0;
position* positions;
int counter = 0;
pthread_t threads[NUM_THREADS];
threadArg threadArgs[NUM_THREADS];
bool running = true;
bool spacePressed = false;

bool editVx = false;
bool editVy = false;
bool editM = false;

bool wallEditing = false;
int wallParam = 0;

int nSym = 6;

Image image;
Texture2D texture;



int getParticleUnderClick(){
	Vector2 pos = GetScreenToWorld2D(GetMousePosition(), cam);
	if (GetMousePosition().x > start){
		return 2;
	}
	double x = pos.x/factor;
	double y = Ly - pos.y/factor; 
	if ((x < Lx) && (y < Ly)){
		int X = coordToCell(x, 1);
		int Y = coordToCell(y, 0);

		for (int j = -1; j <= 1; j++){
			for (int k = -1; k <= 1; k++){
				particleUnderClick = cellList[PBCcell(X + j, 1)][PBCcell(Y + k, 0)];
				while (particleUnderClick != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (pow((particleUnderClick->x - x), 2) + pow(particleUnderClick->y - y, 2) < particleUnderClick->rad*particleUnderClick->rad){
						return 1;
					}
					particleUnderClick = particleUnderClick->nxt;
				}
			}
		}
	}
	return 0;
}



void getInput(){

	if (IsKeyPressed(KEY_SPACE) && !spacePressed){
		spacePressed = true;
		running = !running;
		
	}
	
	if (IsKeyReleased(KEY_SPACE)){
		spacePressed = false;
	}

	 if (IsWindowResized() && !IsWindowFullscreen()){
            screenWidth = GetScreenWidth();
            screenHeight = GetScreenHeight();
			factor = GetScreenHeight()/Lx;
            start = GetScreenHeight();
			xGUI = GetScreenWidth()/1800.;
			yGUI = GetScreenHeight()/900.;
        }

        // check for alt + enter
 		if (IsKeyPressed(KEY_F))
 		{
            // see what display we are on right now
 			int display = GetCurrentMonitor();
 
            
            if (IsWindowFullscreen())
            {
                // if we are full screen, then go back to the windowed size
                SetWindowSize(screenWidth, screenHeight);
            }
            else
            {
                // if we are not full screen, set the window size to match the monitor we are on
                SetWindowSize(GetMonitorWidth(display), GetMonitorHeight(display));
            }
 
            // toggle the state
 			ToggleFullscreen();
 		}



	if (IsKeyPressed(KEY_F)){
		if (!IsWindowFullscreen()){
			int display = GetCurrentMonitor();
			factor = GetMonitorHeight(display)/Lx;
			start = GetMonitorHeight(display);
			xGUI = GetMonitorWidth(display)/1800.;
			yGUI = GetMonitorHeight(display)/900.;
			
			SetWindowSize(GetMonitorWidth(display), GetMonitorHeight(display));
			
			//ToggleFullscreen();
		}    
    }

	float wheel = GetMouseWheelMove();
	if (wheel != 0){
		
		Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), cam);
		

		cam.offset = GetMousePosition();

		cam.target = mouseWorldPos;

		// evil trick
		double speed = 0.0625f*100/Lx;
		
		cam.zoom += wheel * speed*(cam.zoom*cam.zoom);
		if (cam.zoom < 0.2)
			cam.zoom = 0.2;
		if (cam.zoom > 10)
			cam.zoom = 10;
	}

	if (IsKeyDown(KEY_R)){
		cam.zoom = 1;
		cam.offset = (Vector2){0, 0};
		cam.target = (Vector2){0, 0};
	}

	if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)){
		
		if (leftClicked > 0){
			freeFly(particleUnderClick);
			Vector2 pos = GetScreenToWorld2D(GetMousePosition(), cam);

			particleUnderClick->coll++;
		
			double x = pos.x/factor;
			double y = Ly - pos.y/factor; 
			double dx = (particleUnderClick->x - x);
			double dy = (particleUnderClick->y - y);
			PBC(&dx, &dy);

			particleUnderClick->vx -= dx*0.3 + 0.3*particleUnderClick->vx;

			particleUnderClick->vy -= dy*0.3 + 0.3*particleUnderClick->vy;
			
			
			removeEventFromQueue(eventList[particleUnderClick->num]);
			crossingEvent(particleUnderClick->num);

			removeEventFromQueue(eventList[N + particleUnderClick->num]);
			collisionEvent(particleUnderClick->num);
		}
		else{
			int resultat = getParticleUnderClick();

			if (resultat == 1){
				GuiLock();
				leftClicked = 1;
				selected = 1;
			}
			else if (resultat == 0){
				selected = 0;
			}
		}
	}
	else{
		GuiUnlock();
		leftClicked = 0;
	}

}

void GuiSliderBarDouble(Rectangle bounds, const char *textLeft, const char *textRight, double *value, double minValue, double maxValue){
	float floatValue = *value;
	GuiSliderBar(bounds, textLeft, textRight, &floatValue, minValue, maxValue);
	*value = (double)floatValue;
}

void drawParticlesAndBox(){
	double min = 1000000000000;
	double max = -1;

	if (colorParam2 != 0){
		
		for (int i = 0; i < N; i++){
			double v = colorFunction(particles + i);
			colorFunctionArray[i] = v;
			if (min > v){
				min = v;
			}
			if (max < v){
				max = v;
			}
		}
		if (min == max){
			max += 1;
		}
	}
	for (int i = N - 1; i >= 0; i--){	
		if (colorParam2 != 0){
			particleColor = colorSelect((colorFunctionArray[i] - min)/(max - min));
		}
		else{
			particleColor = particleArray[particles[i].type];
		}
		DrawCircleV((Vector2){particles[i].x*factor, (Ly - particles[i].y)*factor}, particles[i].rad*factor, particleColor);
		
	}
	if (selected){
		DrawRing((Vector2){particleUnderClick->x*factor, (Ly - particleUnderClick->y)*factor}, 0.7* particleUnderClick->rad*factor,  particleUnderClick->rad*factor, 0, 360, 20, BLACK);  
	}
	if (addWallx){
		DrawLineEx((Vector2){0, 0}, (Vector2){0, Ly*factor} , 0.3*factor, BLACK);
		DrawLineEx((Vector2){Lx*factor, 0}, (Vector2){Lx*factor, Ly*factor} , 0.3*factor, BLACK);
	}
	if (addWally){
		DrawLineEx((Vector2){0, 0}, (Vector2){Lx*factor, 0} , 0.3*factor, BLACK);
		DrawLineEx((Vector2){0, Ly*factor}, (Vector2){Lx*factor, Ly*factor} , 0.3*factor, BLACK);
	}
	else if (addCircularWall){
		DrawRing((Vector2){Lx/2*factor, Lx/2*factor}, Lx/2*factor, (Lx/2 + 1)*factor, 0, 360, 360, BLACK);
	}
}

void draw(int argc, char *argv[]){

	
	BeginDrawing();

	ClearBackground(RAYWHITE);

	char name[200];	
	
	
	BeginMode2D(cam);
	
	drawParticlesAndBox();
	
	EndMode2D();

	DrawRectangle(start, 0, GetScreenWidth() - GetScreenHeight(), GetScreenHeight(), Fade((Color){220, 220, 220, 255}, 1.f));

	if (selected){
		int changes = 0;
		changes += doubleBox(&(particleUnderClick->vx), "vx", &editVx, (Rectangle){ start + 40*xGUI, 800*yGUI, 120*xGUI, 24*yGUI});
		changes += doubleBox(&(particleUnderClick->vy), "vy", &editVy, (Rectangle){ start + 40*xGUI, 824*yGUI, 120*xGUI, 24*yGUI});
		changes += doubleBox(&(particleUnderClick->m), "m", &editM, (Rectangle){ start + 40*xGUI, 848*yGUI, 120*xGUI, 24*yGUI});
		if (changes){
			particleUnderClick->coll++;

			removeEventFromQueue(eventList[N + particleUnderClick->num]);
			collisionEvent(particleUnderClick->num);
			removeEventFromQueue(eventList[particleUnderClick->num]);
			crossingEvent(particleUnderClick->num);
		}
	}

	if (colorEditing2) 
		GuiLock();

	int dirtyColorParam2 = colorParam2;
	if (GuiDropdownBox((Rectangle){start  + 100*xGUI, 740*yGUI, 120*xGUI, 24*yGUI }, "Uniform; Coll. Based; Vel. Based; Hex. Based; Square. Based", &colorParam2, colorEditing2)){
		
		colorEditing2 = !colorEditing2;
		if (colorParam2 == 0){
			particleColor = MAROON;
		}
		else if (colorParam2 == 1){
			colorFunction = &colorCollision;
		}
		else if (colorParam2 == 2){
			colorFunction = &colorVelocity;
		}
		else if (colorParam2 == 3){
			nSym = 6;
			colorFunction = &colorBOOP;
		}
		else if (colorParam2 == 4){
			nSym = 4;
			colorFunction = &colorBOOP;
		}
		if ((dirtyColorParam2 == 0) && (colorParam2 != 0)){
			colorFunctionArray = calloc(N, sizeof(double));
		}
		else if ((dirtyColorParam2 != 0) && (colorParam2 == 0)){
			free(colorFunctionArray);
		}
	}	

	if (leftClicked == 0)
		GuiUnlock();

	if (colorParam2){
		if (colorEditing) 
			GuiLock();

		if (GuiDropdownBox((Rectangle){start  + 220*xGUI, 740*yGUI, 120*xGUI, 24*yGUI }, "Plasma;Viridis;Copper", &colorParam, colorEditing)){
			colorEditing = !colorEditing;
			if (colorParam == 0){
				colorArray = Plasma;
			}
			else if (colorParam == 1){
				colorArray = Viridis;
			}
			else{
				colorArray = Copper;
			}
		}	

		if (leftClicked == 0)
			GuiUnlock();
	}		
		
	
	int dirtyNoise = noise;
	GuiToggleGroup((Rectangle){start  + 100*xGUI, 40*yGUI, 100*xGUI, 40*yGUI}, "No Thermostat;Langevin;Vel. Rescale", &noise); 
	
	if (dirtyNoise != noise){

		if (noise == 0){
			removeEventFromQueue(eventList[2*N + 2]);

		}
		if (dirtyNoise == 0){
			addEventNoise(t + dtnoise);
		}
	}	

	if (noise){
		sprintf(name, "%.3f", T);
		GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 80*yGUI, 505*xGUI, 40*yGUI}, "Temperature", name, &T, 0.0000000001f, 0.1f);
		if (noise == 1){
			sprintf(name, "%.3f", gamm);
			GuiSliderBarDouble((Rectangle){ start + 100*xGUI, 130*yGUI, 505*xGUI, 40*yGUI }, "Gamma", name, &gamm, 0.001f, 0.1f);
		}
	}
	

	bool dirtyWell = addWell;
	GuiCheckBox((Rectangle){ start  + 100*xGUI, 180*yGUI, 40*xGUI, 40*yGUI}, "Potential", &addWell);
	if (dirtyWell != addWell){
		for(int i = 0; i < N; i++){
			particle* p = particles + i;
			freeFly(p);
			if (addWell == 0){
				free(p->particlesInWell);
			}
			
		}
		for (int i = 0; i < Nxcells; i++){
			free(cellList[i]);
		}
		free(cellList);
		boxConstantHelper();
		cellListInit();
		if (addWell)
			pairsInit();
		
		for (int i = 0; i < N; i++){
			particle* p = particles + i;
			removeEventFromQueue(eventList[N + p->num]);
			collisionEvent(p->num);
			removeEventFromQueue(eventList[p->num]);
			crossingEvent(p->num);
		}
	}

	if (addWell){
		sprintf(name, "%.3f", U);
		GuiSliderBarDouble((Rectangle){ start + 100*xGUI, 220*yGUI, 505*xGUI, 40*yGUI }, "U", name, &U, -0.3f, 0.3f);
		sprintf(name, "%.3f", sig);
		float sigTemp = sig;
		GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 270*yGUI, 505*xGUI, 40*yGUI}, "Pot. rad.", name, &sig, 1.01f, 2.5f);
		if (sig != sigTemp){
			for (int i = 0; i < Nxcells; i++){
					free(cellList[i]);
			}
			free(cellList);
			boxConstantHelper();
			cellListInit();
			for(int i = 0; i < N; i++){
				particle *p = particles + i;
				freeFly(p);
				free(p->particlesInWell);
			} 
			pairsInit();
			
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				removeEventFromQueue(eventList[N + p->num]);
				collisionEvent(p->num);
				removeEventFromQueue(eventList[p->num]);
				crossingEvent(p->num);
				

			}
		}
	}

	bool dirtyField = addField;
	GuiCheckBox((Rectangle){ start  + 100*xGUI, 320*yGUI, 40*xGUI, 40*yGUI}, "Const. Field", &addField);
	if (dirtyField != addField){
		for(int i = 0; i < N; i++){			
			particle* p = particles + i;
			freeFly(p);
			removeEventFromQueue(eventList[N + p->num]);
			collisionEvent(p->num);
			removeEventFromQueue(eventList[p->num]);
			crossingEvent(p->num);
		}
	}

	if (addField){
		double fieldTemp = field;
		sprintf(name, "%.3f", field);
		GuiSliderBarDouble((Rectangle){ start + 100*xGUI, 360*yGUI, 505*xGUI, 40*yGUI}, "g", name, &field, -0.03, -0.0001);
		if (fieldTemp != field){
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				freeFly(p);
				removeEventFromQueue(eventList[N + p->num]);
				collisionEvent(p->num);
				removeEventFromQueue(eventList[p->num]);
				crossingEvent(p->num);

			}
		}
	}
	
	sprintf(name, "%.3f", res);
	GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 580*yGUI, 200*xGUI, 40*yGUI}, "Coeff of res.", name, &res, 0.3f, 1.f);

	float Ntemp = N;
	sprintf(name, "%d", N);
	GuiSliderBar((Rectangle){ start  + 100*xGUI, 615*yGUI, 200*xGUI, 25*yGUI}, "N. of particles", name, &Ntemp, 50.f, 9000.f);
	if ((int)Ntemp != N){
		if (structFactorActivated)
			free(positions);
		freeArrays();
		if (colorParam2)
			free(colorFunctionArray);
		N = (int)Ntemp;
		reset(argc, argv);
		if (colorParam2)
			colorFunctionArray = calloc(N, sizeof(double));
	}


	double sizeratioTemp = sizeratio;
	sprintf(name, "%.3lf", sizeratio);
	GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 640*yGUI, 200*xGUI, 25*yGUI}, "Sizeratio", name, &sizeratio, 0.25, 1);


	double fractionSmallNTemp = fractionSmallN;
	sprintf(name, "%.3lf", fractionSmallN);
	GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 665*yGUI, 200*xGUI, 25*yGUI}, "Frac. Of Small.", name, &fractionSmallN, 0, 1);

	

	double phiTemp = phi;
	sprintf(name, "%.3lf", phi);
	GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 690*yGUI, 200*xGUI, 25*yGUI}, "Packing fraction", name, &phi, 0.1, 0.86);
	if ((phiTemp != phi) || (fractionSmallNTemp != fractionSmallN) || (sizeratioTemp != sizeratio)){
	/*
		if ((phiTemp > phi) && (t > 1/vr)){
			double mult = sqrt(phiTemp/phi);
			Lx = Lx*mult;
			Ly = Lx;
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				freeFly(p);
				p->x = p->x*mult;
				p->y = p->y*mult;
			}
			printf("Lx = %lf\n", Lx);
			for (int i = 0; i < Nxcells; i++){
				free(cellList[i]);
			}
			free(cellList);
			boxConstantHelper();
			cellListInit();
			factor = screenHeight/Ly;
			if (addWell){
				for (int i = 0; i < N; i++){
					particle* p = particles + i;
					free(particles[i].particlesInWell);
				}
				pairsInit();
			}
			for (int i = 0; i < N; i++){
				particle* p = particles + i;
				p->coll++;
				removeEventFromQueue(eventList[N + p->num]);
				collisionEvent(p->num);
				removeEventFromQueue(eventList[p->num]);
				crossingEvent(p->num);
			}
		}
		else{
		*/
		if (structFactorActivated)
			free(positions);
		freeArrays();

		reset(argc, argv);
		//}
	}
	
	int dirtyWallParam = wallParam;
	if (GuiDropdownBox((Rectangle){start  + 100*xGUI, 400*yGUI, 120*xGUI, 24*yGUI }, "No Wall; Horiz. Wall; Vert. Wall; Square Walls; Circl. Wall", &wallParam, wallEditing)){
		wallEditing = !wallEditing;
		if (dirtyWallParam != wallParam){
			if (wallParam == 0){
				addWallx = 0;
				addWally = 0;
				addCircularWall = 0;
			}
			else if (wallParam == 1){
				addWallx = 0;
				addWally = 1;
				addCircularWall = 0;
			}
			else if (wallParam == 2){
				addWallx = 1;
				addWally = 0;
				addCircularWall = 0;
			}
			else if (wallParam == 3){
				addWallx = 1;
				addWally = 1;
				addCircularWall = 0;
			}
			else{
				addWallx = 0;
				addWally = 0;
				addCircularWall = 1;
			}
			if (wallParam != 0){
				if (structFactorActivated)
					free(positions);
				freeArrays();
				reset(argc, argv);
			}
			else{
				if (addWell){
					//To take care of the boundary conditions
					for (int i = 0; i < N; i++){
						particle* p = particles + i;
						freeFly(p);
						free(particles[i].particlesInWell);
					}
					pairsInit();
				}
				for (int i = 0; i < N; i++){
					particle* p = particles + i;
					freeFly(p);
					p->coll++;
					removeEventFromQueue(eventList[N + p->num]);
					collisionEvent(p->num);
					removeEventFromQueue(eventList[p->num]);
					crossingEvent(p->num);
				}
			}
		}
	}	

	
	bool tempStruct = structFactorActivated;
	GuiCheckBox((Rectangle){start  + 500*xGUI, 40*yGUI, 40*xGUI, 40*yGUI}, "Struct. Factor", &structFactorActivated);
	if (tempStruct != structFactorActivated){
		if (structFactorActivated){
			positions = calloc((N - (int)(N*fractionSmallN)), sizeof(position));
			#pragma omp parallel for
			for (int i = 0; i < N - (int)(N*fractionSmallN); i++){
				particle* p = particles + i;
				positions[i].x = p->x;
				positions[i].y = p->y;
			}
			asyncStructFactor();
		}
		else{
			free(positions);
		}
	}
	if (structFactorActivated){
		if ((counter%10 == 0) && (running)){
			//clock_t begin = clock();
			awaitStructFactor();
			//clock_t end = clock();
			//double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			//printf("time took await: %lf\n", time_spent);

			//begin = clock();
			for (int i = 0; i < N - (int)(N*fractionSmallN); i++){
				particle* p = particles + i;
				positions[i].x = p->x;
				positions[i].y = p->y;
			}
			//end = clock();
			//time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			//printf("time took FILL: %lf\n", time_spent);

			
			//begin = clock();
			normalizeStruct();
			//end = clock();
			//time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			//printf("time took NORMALIZE: %lf\n", time_spent);

			//begin = clock();
			UpdateTexture(texture, structFactor);
			//end = clock();
			//time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			//printf("time took UPDATE: %lf\n", time_spent);

			//begin = clock();
			asyncStructFactor();
			//end = clock();
			//time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			//printf("time took ASYNC: %lf\n", time_spent);
			counter = 0;
		}
		DrawTexturePro(texture,
            (Rectangle){ 0, 0, qN, qN },
            (Rectangle){ 400 + start*xGUI, 450*yGUI, SIZE*xGUI, SIZE*yGUI},
            (Vector2){ 0, 0 }, 0.0f, (Color){255, 255, 255, 255});
	}
	DrawFPS(GetScreenWidth() - 100, 10);
	
	EndDrawing();
	counter++;
}


double colorVelocity(particle* p){
	return sqrt(p->vx*p->vx + p->vy*p->vy); 
}

double colorCollision(particle* p){
	return (double)p->coll;
}

double colorBOOP(particle* p1){
	double treshold = 2;
	if ((addWell) && (U < 0)){
		treshold = 2*sig;
	}
	int count = 0;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	double re = 0;
	double im = 0;
	for (int j = -1; j <= 1; j++){
		for (int k = -1; k <= 1; k++){
			particle* p2 = cellList[PBCcell(X + j, 1)][PBCcell(Y + k, 0)];
			while ((p2 != NULL) && (p2->type == 1)){ //while there is a particle in the doubly linked list of the cellList do...
				if ((p1->num != p2->num) && (p1->type == 1)){
					double dy = p1->y - p2->y;
					double dx = p1->x - p2->x;
					PBC(&dx, &dy);
					if (dx*dx + dy*dy < treshold*pow(p1->rad + p2->rad, 2)){
						count++;
						double theta = atan2(dy, dx);
						re += cos(nSym*theta);
						im += sin(nSym*theta);
					}
				}
				p2 = p2->nxt;
			}
		}
	}
	if (count >= 1)
		return (re*re + im*im)/count;
	return 0;
}

Color colorSelect(double value){
	return colorArray[(int)(value*color_size)];
}

void* computeStructureFactor(void* arg){
	threadArg* threadArgument = (threadArg*)arg;
    int start = threadArgument->start;
    int end = threadArgument->end;
	position* p;
	for (int i = start; i < end; i++)
	{
		for (int j = 0; j < qN; j++)
		{
			double im = 0;
			double re = 0;

			for (int n = 0; n < N - (int)(N*fractionSmallN); n++)
			{
				p = positions + n;
				double qr = qx[i]*p->x + qx[j]*p->y;

				re += cos(qr);
				im += sin(qr);

			}
			structFactor[i][j] = (re*re + im*im) / N;
			structFactor[i][j] = log(structFactor[i][j]*structFactor[i][j] + 1);
			if ((i <= qN/2 + 3) && (i >= qN/2 - 3) && (j >= qN/2 - 3) && (j <= qN/2 + 3)){
				structFactor[i][j] = 0;
			}
		}
	}
	pthread_exit(NULL);
}

void normalizeStruct(){
	double max = 0;
	for (int i = 0; i < qN; i++){
		for (int j = 0; j < qN; j++){
			if (structFactor[i][j] > max){
				max = structFactor[i][j];
			}
		}
	}

	for (int i = 0; i < qN; i++){
		for (int j = 0; j < qN; j++){
			structFactor[i][j] = structFactor[i][j]/max;
		}
	}
}


void threadPoolInit(){

    int iterationsPerThread = qN/NUM_THREADS;
    int remainingIterations = qN%NUM_THREADS;
    
    for (int i = 0; i < NUM_THREADS; i++) {
        threadArgs[i].start = i*iterationsPerThread;
        threadArgs[i].end = (i + 1)*iterationsPerThread;
        
        // Distribute remaining iterations among threads
        if (i == NUM_THREADS - 1)
            threadArgs[i].end += remainingIterations;
	}
}

void asyncStructFactor(){
	for (int i = 0; i < NUM_THREADS; i++){
		pthread_create(&threads[i], NULL, computeStructureFactor, (void*)&threadArgs[i]);
	}
}

void awaitStructFactor(){
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }
}

void reset(int argc, char *argv[]){
	t = 0;
	ncol = 0;
	ncross = 0;
	paulTime = 0;
	actualPaulList = 0;
	selected = 0;
	leftClicked = 0;
	counter = 1;

	
	constantInit(argc, argv);
	particlesInit();
	cellListInit();
	//stupid allocation to be able to change sig and sigtemp in the graphical mode even though we don't need it yet!
	if (addWell)
		pairsInit();
	
	eventListInit();
	physicalQ();
	format();
	if (structFactorActivated){
		positions = calloc((N - (int)(N*fractionSmallN)), sizeof(position));

		/*
		Make the sim go laggy when changing N or phi continuously
		Not needed.
		for (int i = 0; i < N; i++){
			particle* p = particles + i;
			positions[i].x = p->x;
			positions[i].y = p->y;
		}
		asyncStructFactor();*/
	}
}

int doubleBox(double* ptr, char* text, bool* activate, Rectangle position){
	char textBuffer[128] = {0};
	char textBuffer2[128] = {0};
	sprintf(textBuffer, "%lf", *ptr);
	sprintf(textBuffer2, "%lf", *ptr);
	if (GuiTextBox(position, textBuffer, sizeof(textBuffer), *activate)){
		*activate = !(*activate);
	}
	Rectangle textBounds = {0};
	textBounds.width = (float)GetTextWidth(text);
	textBounds.height = (float)GuiGetStyle(DEFAULT, TEXT_SIZE);
	textBounds.x = position.x - textBounds.width - GuiGetStyle(SLIDER, TEXT_PADDING);
	textBounds.y = position.y + position.height/2 - GuiGetStyle(DEFAULT, TEXT_SIZE)/2;

	GuiDrawText(text, textBounds, TEXT_ALIGN_RIGHT, Fade(GetColor(GuiGetStyle(SLIDER, TEXT + (STATE_NORMAL*3))), guiAlpha));

	if (strcmp(textBuffer, textBuffer2) == 0){
		return 0;
	}
	int conversion = sscanf(textBuffer, "%lf", ptr);
	//awful dumb trick because sscanf "says" that the conversion is successful even when the number is -0

	if (*ptr == -0.0){
			//for some reason *ptr = 0 doesn't work :) the compiler optimize it out with -Ofast
            *ptr = 0.000000000000000001;
    }

	if (conversion == 1){
		return 1;
	}
	else{
		sscanf(textBuffer2, "%lf", ptr);
		return 0;
	}
}

void graphicalInit(){

	dtime = 1;
	res = 1;
	sig = 1.5;
	U = -0.15;
	gamm = 0.1;
	T = 0.01;
	Einit = 0.01;
	dtnoise = 1;
	load = 0;
	N = 2000;
	phi = 0.5;
	field = -0.001;
	firstScreen = 0;
	colorFunction = &colorCollision;


	cam.zoom = 1;

	for (int i = 0; i < qN; i++){
		qx[i] = -10 + 20.0*i/(double)qN;
	}
	
	threadPoolInit();

	SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT);
    InitWindow(screenWidth, screenHeight, "EDMD");
	SetTargetFPS(144);
	start = GetScreenHeight();
	image = GenImageColor(qN, qN, BLANK);
	image.data = structFactor;
	image.format = PIXELFORMAT_UNCOMPRESSED_R32;
	texture = LoadTextureFromImage(image);


}