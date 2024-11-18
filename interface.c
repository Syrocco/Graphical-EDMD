void initLocalDensity(double lx, double ly, int Dx, int Dy, int n, FILE* inter);
void computeLocalConcentration(particle *particles);
#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>


#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

double** d;
int Nx;
int Ny;
double dx;
double dy;
int N;
double Lx;
double Ly;
FILE* interface;

void initLocalDensity(double lx, double ly, int Dx, int Dy, int n, FILE* inter){
    Nx = Lx / Dx; 
    Ny = Ly / Dy;
    dx = Lx/Nx;
    dy = Ly/Ny; 
    Lx = lx;
    Ly = ly;
    N = n;
    d = malloc(Nx * sizeof(double*));
    interface = inter;

    for (int i = 0; i < Nx; i++) {
        d[i] = calloc(Ny, sizeof(double));
    }
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            d[i][j] = 0;
        }
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
        }
    }
    
    for (int i = 0; i < N; i++) {
        int box_x = (int)(particles[i].x / dx);
        int box_y = (int)(particles[i].y / dy);
        
        if (particles[i].type == 1) {
            d[box_x][box_y] = d[box_x][box_y] + 1/(dx*dy)*M_PI;
        }
    }
    /*
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            d[i][j] = d[i][j]/(dx*dy)*M_PI;
             fprintf(interface, "%lf ", d[i][j]);
        }
        fprintf(interface, "\n");
    }*/
}

void computeInterfacesPos(particle *particles){
    computeLocalConcentration(particles);
    int J = Nx/4;
    for (int i = 0; i < Ny; i++){
        int count = 0;
        for (int j = J; J < Nx; j++){
            if (d[j][i] < 0.2){
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
            if (d[k][i] < 0.1){
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
