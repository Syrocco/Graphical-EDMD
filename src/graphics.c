#define RAYGUI_IMPLEMENTATION
#include "raygui.h" 
#include "color.h"
#include "graphics.h"
#include "EDMD.h"
#include "tex_data.h"  
#include <pthread.h>

#define qN 180
#define SIZE 400
#define NUM_THREADS 11
int nSym = 4; 

double qxG[qN] = {0};
float structFactorGraphics[qN][qN] = {0};
bool structFactorActivated = 0;
position* positions;
pthread_t threads[NUM_THREADS];
threadArg threadArgs[NUM_THREADS];

Image image;
Texture2D texture;

Texture2D circleTexture;

#define MAX_BATCH_ELEMENTS 8192*4

int getWhatsUnderClick(window* screenWindow, state* screenState){
    double factor = screenWindow->factor;
    
	Vector2 pos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);
	if (GetMousePosition().x > screenWindow->start){
		return -1;
	}
	double x = pos.x/factor;
	double y = Ly - pos.y/factor; 

	if (t > 1/vr){
		if ((addWallx) && (x > Lx - 5) && (x < Lx + 5)){
			return 2;
		}
		if ((addWally) && (y > Ly - 5) && (y < Ly + 5)){
			return 3;
		}
		}
	if ((x < Lx) && (y < Ly)){
		int X = coordToCell(x, 1);
		int Y = coordToCell(y, 0);

		for (int j = -1; j <= 1; j++){
			for (int k = -1; k <= 1; k++){
				screenState->particleUnderClick =  cellList[PBCcellY(Y + j) * Nxcells + PBCcellX(X + k)];
				while (screenState->particleUnderClick != NULL){ //while there is a particle in the doubly linked list of the cellList do...
					if (pow((screenState->particleUnderClick->x - x), 2) + pow(screenState->particleUnderClick->y - y, 2) < screenState->particleUnderClick->rad*screenState->particleUnderClick->rad){
						return 1;
					}
					screenState->particleUnderClick = screenState->particleUnderClick->nxt;
				}
			}
		}
	}
	return 0;
}



void getInput(window* screenWindow, state* screenState){
	particle* particleUnderClick = screenState->particleUnderClick;

	 if (IsWindowResized() && !IsWindowFullscreen()){
            screenWindow->screenWidth = GetScreenWidth();
            screenWindow->screenHeight = GetScreenHeight();
			screenWindow->factor = GetScreenHeight()/Lx;
            screenWindow->start = GetScreenHeight();
			screenWindow->xGUI = GetScreenWidth()/1800.;
			screenWindow->yGUI = GetScreenHeight()/900.;
        }

        // check for alt + enter
 		if (IsKeyPressed(KEY_F))
 		{
            // see what display we are on right now
 			int display = GetCurrentMonitor();
 
            
            if (IsWindowFullscreen())
            {
                // if we are full screen, then go back to the windowed size
                SetWindowSize(screenWindow->screenWidth, screenWindow->screenHeight);
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
			screenWindow->factor = GetMonitorHeight(display)/Ly;
			screenWindow->start = GetMonitorHeight(display);
			screenWindow->xGUI = GetMonitorWidth(display)/1800.;
			screenWindow->yGUI = GetMonitorHeight(display)/900.;
			
			SetWindowSize(GetMonitorWidth(display), GetMonitorHeight(display));
			
			//ToggleFullscreen();
		}    
    }

	float wheel = GetMouseWheelMove();
	if (wheel != 0){
		
		Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);
		

		screenWindow->cam.offset = GetMousePosition();

		screenWindow->cam.target = mouseWorldPos;

		// evil trick
		double speed = 0.0625f*100/Ly;
		
		screenWindow->cam.zoom += wheel * speed*(screenWindow->cam.zoom*screenWindow->cam.zoom);
		if (screenWindow->cam.zoom < 0.2)
			screenWindow->cam.zoom = 0.2;
		if (screenWindow->cam.zoom > 10)
			screenWindow->cam.zoom = 10;
	}

	if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT))
		{
			Vector2 delta = GetMouseDelta();
			float scale = -1./screenWindow->cam.zoom;
			delta.x = delta.x*scale;
			delta.y = delta.y*scale;
			screenWindow->cam.target.x += delta.x;
			screenWindow->cam.target.y += delta.y;
		}

	if (IsKeyDown(KEY_R)){
		screenWindow->cam.zoom = 1;
		screenWindow->cam.offset = (Vector2){0, 0};
		screenWindow->cam.target = (Vector2){0, 0};
	}
	if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)){
		
		if (screenState->leftClicked > 0){
			freeFly(particleUnderClick);
			Vector2 pos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);

			particleUnderClick->coll++;
		
			double x = pos.x/screenWindow->factor;
			double y = Ly - pos.y/screenWindow->factor; 
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
			int resultat = getWhatsUnderClick(screenWindow, screenState);
			if (resultat == 2){
				screenState->wallMoving = true;
				Vector2 pos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);
				if (Lx < pos.x/screenWindow->factor)
					Lx = pos.x/screenWindow->factor;
			}
			if (resultat == 3){
				screenState->wallMoving = true;
				Vector2 pos = GetScreenToWorld2D(GetMousePosition(), screenWindow->cam);
				if (Ly < (Ly - pos.y/screenWindow->factor)){
					Ly = (Ly - pos.y/screenWindow->factor);
					screenWindow->cam.target.y += fabs(pos.y);
				}
			}
			if (resultat == 1){
				GuiLock();
				screenState->leftClicked = 1;
				screenState->selected = 1;
			}
			else if (resultat == 0){
				screenState->selected = 0;
			}
		}
	}
	else{
		GuiUnlock();
		screenState->leftClicked = 0;
		if (screenState->wallMoving){
			screenState->wallMoving = false;
			for(int i = 0; i < N; i++){
				particle* p = particles + i;
				freeFly(p);
				if (addWell){
					free(p->particlesInWell);
				}
				
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
	}

	if (IsKeyPressed(KEY_SPACE) && !screenState->spacePressed){
		screenState->spacePressed = true;
		screenState->running = !screenState->running;
	}
	
	if (IsKeyReleased(KEY_SPACE)){
		screenState->spacePressed = false;
	}
}

void GuiSliderBarDouble(Rectangle bounds, const char *textLeft, const char *textRight, double *value, double minValue, double maxValue){
	float floatValue = *value;
	GuiSliderBar(bounds, textLeft, textRight, &floatValue, minValue, maxValue);
	*value = (double)floatValue;
}

void drawParticlesAndBox(double factor, state* screenState, float zoom){
	double min = 1000000000000;
	double max = -1;
	particle* particleUnderClick = screenState->particleUnderClick;
	if (screenState->colorParam2 != 0){
		
		for (int i = 0; i < N; i++){
			double v = screenState->colorFunction(particles + i);
			screenState->colorFunctionArray[i] = v;
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
		if (screenState->colorParam2 != 0){
			screenState->particleColor = colorSelect((screenState->colorFunctionArray[i] - min)/(max - min), screenState->colorArray);
		}
		else{
			screenState->particleColor = screenState->particleArray[particles[i].type];
		}

		Rectangle sourceRect = {0.0f, 0.0f, circleTexture.width, circleTexture.height};

		// Calculate destination rectangle
		Rectangle destRect = {(particles[i].x - particles[i].rad)*factor, (Ly - particles[i].y - particles[i].rad)*factor, factor*particles[i].rad * 2, factor*particles[i].rad* 2};

		// Draw the circle using the texture
		DrawTexturePro(circleTexture, sourceRect, destRect, (Vector2){0, 0}, 0.0f, screenState->particleColor);

		//DrawCircleV((Vector2){particles[i].x*factor, (Ly - particles[i].y)*factor}, particles[i].rad*factor, screenState->particleColor);

	}
	if (screenState->selected){
		DrawRing((Vector2){particleUnderClick->x*factor, (Ly - particleUnderClick->y)*factor}, 0.7* particleUnderClick->rad*factor,  particleUnderClick->rad*factor, 0, 360, 20, BLACK);  
	}
	if (addWallx){
		DrawLineEx((Vector2){0, 0}, (Vector2){0, Ly*factor} , 0.3*factor/zoom, BLACK);
		DrawLineEx((Vector2){Lx*factor, 0}, (Vector2){Lx*factor, Ly*factor} , 0.3*factor/zoom, BLACK);
	}
	if (addWally){
		DrawLineEx((Vector2){0, 0}, (Vector2){Lx*factor, 0} , 0.3*factor/zoom, BLACK);
		DrawLineEx((Vector2){0, Ly*factor}, (Vector2){Lx*factor, Ly*factor} , 0.3*factor/zoom, BLACK);
	}
	else if (addCircularWall){
		DrawRing((Vector2){Lx/2*factor, Lx/2*factor}, Lx/2*factor, (Lx/2 + 1)*factor, 0, 360, 360, BLACK);
	}
}

void draw(int argc, char *argv[], window* screenWindow, state* screenState){
    int start = screenWindow->start;
	double xGUI = screenWindow->xGUI;
    double yGUI = screenWindow->yGUI;
    double factor = screenWindow->factor;
	particle* particleUnderClick = screenState->particleUnderClick;

	BeginDrawing();

	ClearBackground(RAYWHITE);

	char name[200];	
	
	BeginMode2D(screenWindow->cam);
	drawParticlesAndBox(factor, screenState, screenWindow->cam.zoom);
	
	EndMode2D();

	DrawRectangle(start, 0, GetScreenWidth() - GetScreenHeight(), GetScreenHeight(), Fade((Color){220, 220, 220, 255}, 1.f));

	if (screenState->selected){
		int changes = 0;
		changes += doubleBox(&(particleUnderClick->vx), "vx", &screenState->editVx, (Rectangle){ start + 40*xGUI, 800*yGUI, 120*xGUI, 24*yGUI});
		changes += doubleBox(&(particleUnderClick->vy), "vy", &screenState->editVy, (Rectangle){ start + 40*xGUI, 824*yGUI, 120*xGUI, 24*yGUI});
		changes += doubleBox(&(particleUnderClick->m), "m", &screenState->editM, (Rectangle){ start + 40*xGUI, 848*yGUI, 120*xGUI, 24*yGUI});
		if (changes){
			particleUnderClick->coll++;

			removeEventFromQueue(eventList[N + particleUnderClick->num]);
			collisionEvent(particleUnderClick->num);
			removeEventFromQueue(eventList[particleUnderClick->num]);
			crossingEvent(particleUnderClick->num);
		}
	}

	if (screenState->colorEditing2) 
		GuiLock();

	int dirtyColorParam2 = screenState->colorParam2;
	if (GuiDropdownBox((Rectangle){start  + 100*xGUI, 755*yGUI, 120*xGUI, 15*yGUI }, "Uniform; Coll. Based; Vel. Based; Hex. Based; Square. Based; Radius Based; Charge Based", &screenState->colorParam2, screenState->colorEditing2)){
		
		screenState->colorEditing2 = !screenState->colorEditing2;
		if (screenState->colorParam2 == 0){
			screenState->particleColor = MAROON;
		}
		else if (screenState->colorParam2 == 1){
			screenState->colorFunction = &colorCollision;
		}
		else if (screenState->colorParam2 == 2){
			screenState->colorFunction = &colorVelocity;
		}
		else if (screenState->colorParam2 == 3){
			nSym = 6;
			screenState->colorFunction = &colorBOOP;
		}
		else if (screenState->colorParam2 == 4){
			nSym = 4;
			screenState->colorFunction = &colorBOOP;
		}
		else if (screenState->colorParam2 == 5){
			nSym = 4;
			screenState->colorFunction = &colorRadius;
		}
		else if (screenState->colorParam2 == 6){
			screenState->colorFunction = &colorCharge;
		}
		if ((dirtyColorParam2 == 0) && (screenState->colorParam2 != 0)){
			screenState->colorFunctionArray = calloc(N, sizeof(double));
		}
		else if ((dirtyColorParam2 != 0) && (screenState->colorParam2 == 0)){
			free(screenState->colorFunctionArray);
		}
	}

	bool tempPoly = polydispersity;
	GuiCheckBox((Rectangle){start + 100*xGUI, 715*yGUI, 40*xGUI, 40*yGUI}, "Polydisp.", &polydispersity);
	if (polydispersity != tempPoly){
		if (polydispersity){
			//small hack to get the structure factor right!
			fractionSmallN = 0;
		}
		if (structFactorActivated)
			free(positions);
		freeArrays();

		reset(argc, argv, &screenWindow->factor, screenState);
		
	}

	if (screenState->leftClicked == 0)
		GuiUnlock();

	if (screenState->colorParam2){
		if (screenState->colorEditing) 
			GuiLock();

		if (GuiDropdownBox((Rectangle){start  + 220*xGUI, 755*yGUI, 120*xGUI, 15*yGUI }, "Plasma;Viridis;Copper", &screenState->colorParam, screenState->colorEditing)){
			screenState->colorEditing = !screenState->colorEditing;
			if (screenState->colorParam == 0){
				screenState->colorArray = Plasma;
			}
			else if (screenState->colorParam == 1){
				screenState->colorArray = Viridis;
			}
			else{
				screenState->colorArray = Copper;
			}
		}	

		if (screenState->leftClicked == 0)
			GuiUnlock();
	}		
	
	bool dirtyDamping = damping;
	GuiCheckBox((Rectangle){start + 406*xGUI, 40*yGUI, 40*xGUI, 40*yGUI}, "Damping", &damping);
	if (dirtyDamping != damping){
		for(int i = 0; i < N; i++){			
			particle* p = particles + i;
			freeFly(p);
			removeEventFromQueue(eventList[N + p->num]);
			collisionEvent(p->num);
			removeEventFromQueue(eventList[p->num]);
			crossingEvent(p->num);
		}
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

	if (noise != 0 || damping != 0){
		bool dirtyGamm = gamm;
		sprintf(name, "%.3f", gamm);
		GuiSliderBarDouble((Rectangle){ start + 100*xGUI, 130*yGUI, 505*xGUI, 40*yGUI }, "Gamma", name, &gamm, 0.001f, 0.1f);
		if (dirtyGamm != gamm){
			for(int i = 0; i < N; i++){			
				particle* p = particles + i;
				freeFly(p);
				removeEventFromQueue(eventList[N + p->num]);
				collisionEvent(p->num);
				removeEventFromQueue(eventList[p->num]);
				crossingEvent(p->num);
			}
		}
		if (noise == 1){
			sprintf(name, "%.3f", T);
			GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 80*yGUI, 505*xGUI, 40*yGUI}, "Temperature", name, &T, 0.0000000001f, 1.f);
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
		sprintf(name, "%.6f", U);
		GuiSliderBarDouble((Rectangle){ start + 100*xGUI, 220*yGUI, 505*xGUI, 40*yGUI }, "U", name, &U, -1.f, 1.f);
		sprintf(name, "%.3f", sig);
		float sigTemp = sig;
		GuiSliderBarDouble((Rectangle){ start  + 100*xGUI, 270*yGUI, 505*xGUI, 40*yGUI}, "Pot. rad.", name, &sig, 1.01f, 2.5f);
		if (sig != sigTemp){
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
		
		GuiCheckBox((Rectangle){ start  + 200*xGUI, 180*yGUI, 40*xGUI, 40*yGUI}, "Charged", &charged);
		if (charged){
			double dirtyPropP = proportionPositivelyCharged;
			sprintf(name, "%.2f", proportionPositivelyCharged);
			GuiSliderBarDouble((Rectangle){ start + 400*xGUI, 180*yGUI, 50*xGUI, 40*yGUI }, "%+", name, &proportionPositivelyCharged, 0.f, 0.5f);
			double dirtyPropN = proportionNeutralyCharged;
			sprintf(name, "%.2f", proportionNeutralyCharged);
			GuiSliderBarDouble((Rectangle){ start + 500*xGUI, 180*yGUI, 50*xGUI, 40*yGUI }, "%n", name, &proportionNeutralyCharged, 0.f, 0.5f);
			if ((dirtyPropP != proportionPositivelyCharged) || (dirtyPropN != proportionNeutralyCharged)){
				printf("Here!");
				for (int i = 0; i < N; i++){
					double randomTemp = drand(0, 1);
					if (randomTemp < proportionPositivelyCharged){
						particles[i].charge = 1;
					}
					else if (randomTemp < proportionNeutralyCharged + proportionPositivelyCharged){
						particles[i].charge = 0;
					}
					else{
						particles[i].charge = -1;
					}
				}
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

	GuiCheckBox((Rectangle){start + 250*xGUI, 320*yGUI, 40*xGUI, 40*yGUI}, "Delta", &addDelta);
	GuiCheckBox((Rectangle){start + 400*xGUI, 320*yGUI, 40*xGUI, 40*yGUI}, "Tangent Delta", &addDeltaTangent);
	if (addDelta || addDeltaTangent){
		sprintf(name, "%.3f", delta);
		GuiSliderBarDouble((Rectangle){ start + 100*xGUI, 360*yGUI, 505*xGUI, 20*yGUI}, "Delta", name, &delta, 0, 0.1);
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
	GuiSliderBarDouble((Rectangle){start + 100*xGUI, 580*yGUI, 200*xGUI, 40*yGUI}, "Coeff of res.", name, &res, 0.f, 1.f);

	float Ntemp = N;
	sprintf(name, "%d", N);
	GuiSliderBar((Rectangle){start + 100*xGUI, 615*yGUI, 200*xGUI, 25*yGUI}, "N. of particles", name, &Ntemp, 50.f, 50000.f);
	if ((int)Ntemp != N){
		if (structFactorActivated)
			free(positions);
		freeArrays();
		if (screenState->colorParam2)
			free(screenState->colorFunctionArray);
		N = (int)Ntemp;
		reset(argc, argv, &screenWindow->factor, screenState);
		if (screenState->colorParam2)
			screenState->colorFunctionArray = calloc(N, sizeof(double));
	}

	double sizeratioTemp = sizeratio;
	double fractionSmallNTemp = fractionSmallN;
	if (!polydispersity){
		
		sprintf(name, "%.3lf", sizeratio);
		GuiSliderBarDouble((Rectangle){start + 100*xGUI, 640*yGUI, 200*xGUI, 25*yGUI}, "Sizeratio", name, &sizeratio, 0.25, 1);


		sprintf(name, "%.3lf", fractionSmallN);
		GuiSliderBarDouble((Rectangle){start + 100*xGUI, 665*yGUI, 200*xGUI, 25*yGUI}, "Frac. Of Small.", name, &fractionSmallN, 0, 1);
	}
	

	double phiTemp = phi;
	sprintf(name, "%.3lf", phi);
	GuiSliderBarDouble((Rectangle){start + 100*xGUI, 690*yGUI, 200*xGUI, 25*yGUI}, "Packing fraction", name, &phi, 0.1, 0.88);
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

		reset(argc, argv, &screenWindow->factor, screenState);
        
		//}
	}
	



	int dirtyWallParam = screenState->wallParam;
	if (GuiDropdownBox((Rectangle){start + 100*xGUI, 400*yGUI, 120*xGUI, 24*yGUI }, "No Wall; Horiz. Wall; Vert. Wall; Square Walls; Circl. Wall", &screenState->wallParam, screenState->wallEditing)){
		screenState->wallEditing = !screenState->wallEditing;
		if (dirtyWallParam != screenState->wallParam){
			if (screenState->wallParam == 0){
				addWallx = 0;
				addWally = 0;
				addCircularWall = 0;
			}
			else if (screenState->wallParam == 1){
				addWallx = 0;
				addWally = 1;
				addCircularWall = 0;
			}
			else if (screenState->wallParam == 2){
				addWallx = 1;
				addWally = 0;
				addCircularWall = 0;
			}
			else if (screenState->wallParam == 3){
				addWallx = 1;
				addWally = 1;
				addCircularWall = 0;
			}
			else if (screenState->wallParam == 4){
				addWallx = 0;
				addWally = 0;
				addCircularWall = 1;
			}
			if (screenState->wallParam != 0){
				
				if (structFactorActivated)
					free(positions);
				freeArrays();
				reset(argc, argv, &screenWindow->factor, screenState);
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
	GuiCheckBox((Rectangle){start + 500*xGUI, 40*yGUI, 40*xGUI, 40*yGUI}, "Struct. Factor", &structFactorActivated);
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
		if (finishedStructComputation() && (screenState->running)){
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
			UpdateTexture(texture, structFactorGraphics);
			//end = clock();
			//time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			//printf("time took UPDATE: %lf\n", time_spent);

			//begin = clock();
			asyncStructFactor();
			//end = clock();
			//time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			//printf("time took ASYNC: %lf\n", time_spent);
			
		}
		DrawTexturePro(texture,
            (Rectangle){ 0, 0, qN, qN },
            (Rectangle){ 400 + start*xGUI, 450*yGUI, SIZE*xGUI, SIZE*yGUI},
            (Vector2){ 0, 0 }, 0.0f, (Color){255, 255, 255, 255});
	}
	DrawFPS(GetScreenWidth() - 100, 10);
	
	EndDrawing();
}


double colorVelocity(particle* p){
	return sqrt(p->vx*p->vx + p->vy*p->vy); 
}

double colorCollision(particle* p){
	return (double)p->coll;
}

double colorRadius(particle* p){
	return p->rad;
}


double colorCharge(particle* p){
	return p->charge;
}

double colorBOOP(particle* p1){
	double treshold = 2;
	if ((addWell) && (U > 0)){
		treshold *= sig;
	}
	int count = 0;
	int X = p1->cell[0];
	int Y = p1->cell[1];
	double re = 0;
	double im = 0;
	for (int j = -1; j <= 1; j++){
		for (int k = -1; k <= 1; k++){
			particle* p2 =  cellList[PBCcellY(Y + j) * Nxcells + PBCcellX(X + k)];
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

Color colorSelect(double value, Color* colorArray){
	return colorArray[(int)(value*color_size)];
}

void* computeStructureFactorGraphic(void* arg){
	threadArg* threadArgument = (threadArg*)arg;
    int start = threadArgument->start;
    int end = threadArgument->end;
	position* p;
	for (int i = start; i < end; i++)
	{
		for (int j = 0; j < qN; j++)
		{
			float im = 0;
			float re = 0;

			for (int n = 0; n < N - (int)(N*fractionSmallN); n++)
			{
				p = positions + n;
				float qr = qxG[i]*p->x + qxG[j]*p->y;

				re += cos(qr);
				im += sin(qr);

			}
			structFactorGraphics[i][j] = (re*re + im*im) / N;
			structFactorGraphics[i][j] = log(structFactorGraphics[i][j]*structFactorGraphics[i][j] + 1);
			if ((i <= qN/2 + 3) && (i >= qN/2 - 3) && (j >= qN/2 - 3) && (j <= qN/2 + 3)){
				structFactorGraphics[i][j] = 0;
			}
		}
	}
	threadArgument->resting = 1;
	pthread_exit(NULL);
}

void normalizeStruct(){
	float max = 0;
	for (int i = 0; i < qN; i++){
		for (int j = 0; j < qN; j++){
			if (structFactorGraphics[i][j] > max){
				max = structFactorGraphics[i][j];
			}
		}
	}

	for (int i = 0; i < qN; i++){
		for (int j = 0; j < qN; j++){
			structFactorGraphics[i][j] = structFactorGraphics[i][j]/max;
		}
	}
}


void threadPoolInit(){

    int iterationsPerThread = qN/NUM_THREADS;
    int remainingIterations = qN%NUM_THREADS;
    
    for (int i = 0; i < NUM_THREADS; i++) {
        threadArgs[i].start = i*iterationsPerThread;
        threadArgs[i].end = (i + 1)*iterationsPerThread;
        threadArgs[i].resting = 1;
        // Distribute remaining iterations among threads
        if (i == NUM_THREADS - 1)
            threadArgs[i].end += remainingIterations;
	}
}

bool finishedStructComputation(){
	int temp = 0;
	for (int i = 0; i < NUM_THREADS; i++){
		temp += threadArgs[i].resting;
	}
	return (temp == NUM_THREADS);
}

void asyncStructFactor(){
	for (int i = 0; i < NUM_THREADS; i++){
		threadArgs[i].resting = 0;
		pthread_create(&threads[i], NULL, computeStructureFactorGraphic, (void*)&threadArgs[i]);
	}
}

void awaitStructFactor(){
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }
}

void reset(int argc, char *argv[], double* factor, state* screenState){
	t = 0;
	ncol = 0;
	ncross = 0;
	paulTime = 0;
	actualPaulList = 0;
	screenState->selected = 0;
	screenState->leftClicked = 0;

	
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
	
	*factor = GetScreenHeight()/Ly;
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

window graphicalInit(){

	firstScreen = 0;
	load = 0;
	dtnoise = 1;
	dtime = 1;
	N = 2000;
	phi = 0.5;
	field = -0.001;
	res = 1;
	delta = 0.01;
	damping = 0.02;
	sig = 1.5;
	U = 0.15;
	gamm = 0.1;
	T = 0.01;
	Einit = 0.01;

	load = 0;
	S1 = 0;
	Hex = 0;
	coexistenceOld = 0;
	coexistence = 0;

	interfaceThermo = 0;
	interfaceThermoTemp = 0;

	clusterThermo = 0;
	dumpCluster = 0;
	killCluster = 0;
	strucThermo = 0;
	areaThermo = 0;
	boopThermo = 0;
	pcfThermo = 0;
	pcfBondOrderThermo = 0;
	pcfg6Thermo = 0;
	critical = 0;
	snapshotCritical = 0;

	aspectRatio = 1;
    Camera2D cam = {0};
    cam.zoom = 1;

	for (int i = 0; i < qN; i++){
		qxG[i] = -10 + 20.0*i/(double)qN;
	}
	
	threadPoolInit();
	SetConfigFlags(FLAG_WINDOW_RESIZABLE);// | FLAG_VSYNC_HINT);
	
    InitWindow(1800, 900, "EDMD");
	SetTargetFPS(144);
	image = GenImageColor(qN, qN, BLANK);
	UnloadImageColors(image.data);

	
	image.data = structFactorGraphics;
	image.format = PIXELFORMAT_UNCOMPRESSED_R32;
	texture = LoadTextureFromImage(image);
	
	// Load circle texture from embedded data
	Image circleImage = LoadImageFromMemory(".png", src_tex_png, src_tex_png_len);
	circleTexture = LoadTextureFromImage(circleImage);
	UnloadImage(circleImage);  // Free the image data after creating texture
	
    return (window){
        .screenWidth = 1800,
        .screenHeight = 900,
        .factor = 1,
        .start = GetScreenHeight(),
        .xGUI = 1,
        .yGUI = 1,
        .cam = cam,
    };
}

state GUIinit(){
	return (state){
		.leftClicked = 0,
		.selected = 0,
		.particleUnderClick = NULL,
		.particleArray = {GRAY, MAROON},
		.particleColor = BLACK,
		.colorEditing = false,
		.colorParam = 0,
		.colorEditing2 = false,
		.colorParam2 = 0,
		.colorArray = Plasma,
		.colorFunction = colorCollision,
		.colorFunctionArray = NULL,
		.running = true,
		.wallMoving = false,
		.spacePressed = false,
		.editVx = false,
		.editVy = false,
		.editM = false,
		.wallEditing = false,
		.wallParam = 0,
	};
}

void graphicsFree(state* screenState){
	if (structFactorActivated)
		free(positions);
	UnloadTexture(texture);     
	if (screenState->colorParam2)
		free(screenState->colorFunctionArray);
}
