void init_local_density(double lx, double ly, int Dx, int Dy, int n, FILE* inter);
void compute_local_concentration(particle *particles);
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

void init_local_density(double lx, double ly, int Dx, int Dy, int n, FILE* inter){
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
    fprintf(interface, "%d %d\n", Nx, Ny);
}


void compute_local_concentration(particle *particles){
    
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            d[i][j] = 0;
        }
    }
    
    for (int i = 0; i < N; i++) {
        int box_x = (int)(particles[i].x / dx);
        int box_y = (int)(particles[i].y / dy);
        
        if (particles[i].type == 0) {
            d[box_x][box_y]++;
        }
    }
    
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            d[i][j] = d[i][j]/(dx*dy)*M_PI;
             fprintf(interface, "%lf ", d[i][j]);
        }
        fprintf(interface, "\n");
    }
}
