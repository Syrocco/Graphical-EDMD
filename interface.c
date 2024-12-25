#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>

void saveDensityCoarse(particle *particles);
void initLocalDensity(double lx, double ly, int Dx, int Dy, int n, FILE* inter, int liquidliquid);
void computeLocalConcentration(particle *particles);
void freeInterface();

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

double** d;
double** dbis;
int Nx;
int Ny;
double dx;
double dy;
int N;
double Lx;
double Ly;
FILE* interface;

void initLocalDensity(double lx, double ly, int Dx, int Dy, int n, FILE* inter, int liquidliquid){
    Nx = Lx / Dx; 
    Ny = Ly / Dy;
    dx = Lx/Nx;
    dy = Ly/Ny; 
    Lx = lx;
    Ly = ly;
    N = n;
    d = malloc(Nx * sizeof(double*));
    if (liquidliquid){
       dbis = malloc(Nx * sizeof(double*)); 
    }
    interface = inter;

    for (int i = 0; i < Nx; i++) {
        d[i] = calloc(Ny, sizeof(double));
        if (liquidliquid){
            dbis[i] = calloc(Ny, sizeof(double));
        }
    }
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            d[i][j] = 0;
            if (liquidliquid){
                dbis[i][j] = 0;
            }
        }
    }
    if (liquidliquid == 0){
        for (int i = 0; i < Nx; i++){
            fprintf(interface, "%lf ", dx*(i+0.5));
        }
        fprintf(interface, "\n");
    }

    for (int i = 0; i < Ny; i++){
        fprintf(interface, "%lf ", dy*(i+0.5));
    }
    fprintf(interface, "\n");
}


void computeLocalConcentration(particle *particles){
    
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            d[i][j] = 0;
            if (dbis != NULL)
               dbis[i][j] = 0;
        }
    }
    
    for (int i = 0; i < N; i++) {
        int box_x = fmin((int)(particles[i].x / dx), Nx - 1);
        int box_y = fmin((int)(particles[i].y / dy), Ny - 1);
        
        if (particles[i].type == 1) {
            d[box_x][box_y] += M_PI/(dx*dy);
        }
        if (dbis != NULL){
            if (particles[i].type == 0) {
                dbis[box_x][box_y] += M_PI/(dx*dy);
            }
        }
    }
}

void saveDensityCoarse(particle *particles){
    computeLocalConcentration(particles);
    
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            fprintf(interface, "%lf ", d[i][j]);
        }
        fprintf(interface, "\n");
    }
}

void computeInterfacesPos(particle *particles, double tresh){
    computeLocalConcentration(particles);
    int J = Nx/4;
    for (int i = 0; i < Ny; i++){
        int count = 0;
        for (int j = J; j < Nx; j++){
            if ((d[j][i] < tresh) && (dbis[j][i] > tresh)){
                count++;
            }
            else{
                count = 0;
            }
            if (count == 2){
                fprintf(interface, "%lf ", (j - 1)*dy);
                break;
            }
            if (j == Nx - 1){
                fprintf(interface, "nan "); 
            } 
        } 
    }
    fprintf(interface, "\n");
    for (int i = 0; i < Ny; i++){
        int count = 0;
        for (int j = J; j > -Nx/4; j--){
            int k = j;
            if (j < 0){
                k =  Nx + j;
            }
            if ((d[k][i] < tresh) && (dbis[k][i] > tresh)){
                count++;
            }
            else{
                count = 0;
            }
            if (count == 2){
                fprintf(interface, "%lf ", (j + 2 + 0.5)*dy);
                break;
            }
            if (j == (int)(3*Nx/4)){
                fprintf(interface, "nan "); 
            } 
        } 
    }
    fprintf(interface, "\n");
}


void freeInterface(){
    for (int i = 0; i < Nx; i++) {
        free(d[i]);
        if (dbis != NULL) {
            free(dbis[i]);
        }
    }
    free(d);
    if (dbis != NULL) {
        free(dbis);
    }
}
