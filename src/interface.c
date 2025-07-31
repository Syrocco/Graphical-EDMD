#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>

void saveDensityCoarse(particle *particles);
#if THREE_D
void initLocalDensity(double lx, double ly, double lz, int Dx, int Dy, int Dz, int n, FILE* inter, int liquidliquid, int temperature);
#else
void initLocalDensity(double lx, double ly, int Dx, int Dy, int n, FILE* inter, int liquidliquid, int temperature);
#endif
void computeLocalConcentration(particle *particles);
void freeInterface();

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#if THREE_D
    double*** d;
    double*** dbis;
    double*** tp;
#else
    double** d;
    double** dbis;
    double** tp;
#endif

int Nx;
int Ny;
double dx;
double dy;
double Lx;
double Ly;
#if THREE_D
    double dz;
    double Lz;
    int Nz;
#endif
int N;
FILE* interface;

#if THREE_D
void initLocalDensity(double lx, double ly, double lz, int Dx, int Dy, int Dz, int n, FILE* inter, int liquidliquid, int temperature){
#else
void initLocalDensity(double lx, double ly, int Dx, int Dy, int n, FILE* inter, int liquidliquid, int temperature){
#endif
    Nx = lx / Dx; 
    Ny = ly / Dy;
    dx = lx/Nx;
    dy = ly/Ny; 
    Lx = lx;
    Ly = ly;
    #if THREE_D
        Nz = lz/Dz;    
        dz = lz/Nz;
        Lz = lz;
    #endif
    N = n;
    
    #if THREE_D
        // Allocate 3D arrays
        d = malloc(Nx * sizeof(double**));
        if (temperature){
            tp = malloc(Nx * sizeof(double**));
        }
        if (liquidliquid){
           dbis = malloc(Nx * sizeof(double**)); 
        }
        
        for (int i = 0; i < Nx; i++) {
            d[i] = malloc(Ny * sizeof(double*));
            if (liquidliquid){
                dbis[i] = malloc(Ny * sizeof(double*));
            }
            if (temperature){
                tp[i] = malloc(Ny * sizeof(double*));
            }
            
            for (int j = 0; j < Ny; j++) {
                d[i][j] = calloc(Nz, sizeof(double));
                if (liquidliquid){
                    dbis[i][j] = calloc(Nz, sizeof(double));
                }
                if (temperature){
                    tp[i][j] = calloc(Nz, sizeof(double));
                }
            }
        }
    #else
        // Original 2D allocation
        d = malloc(Nx * sizeof(double*));
        if (temperature){
            tp = malloc(Nx * sizeof(double*));
        }
        if (liquidliquid){
           dbis = malloc(Nx * sizeof(double*)); 
        }
        
        for (int i = 0; i < Nx; i++) {
            d[i] = calloc(Ny, sizeof(double));
            if (liquidliquid){
                dbis[i] = calloc(Ny, sizeof(double));
            }
            if (temperature){
                tp[i] = calloc(Ny, sizeof(double));
            }
        }
    #endif
    
    interface = inter;

    #if THREE_D
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nz; k++) {
                    d[i][j][k] = 0;
                    if (liquidliquid){
                        dbis[i][j][k] = 0;
                    }
                }
            }
        }
    #else
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                d[i][j] = 0;
                if (liquidliquid){
                    dbis[i][j] = 0;
                }
            }
        }
    #endif
    
    if (liquidliquid == 0){
        for (int i = 0; i < Nx; i++){
            fprintf(interface, "%lf ", dx*(i+0.5));
        }
        fprintf(interface, "\n");
    }

    if (temperature == 0){
        for (int i = 0; i < Ny; i++){
            fprintf(interface, "%lf ", dy*(i+0.5));
        }
        fprintf(interface, "\n");
        
        #if THREE_D
        for (int i = 0; i < Nz; i++){
            fprintf(interface, "%lf ", dz*(i+0.5));
        }
        fprintf(interface, "\n");
        #endif
    }
}


void computeLocalConcentration(particle *particles){
    
    #if THREE_D
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                d[i][j][k] = 0;
                if (dbis != NULL)
                   dbis[i][j][k] = 0;
                if (tp != NULL){
                    tp[i][j][k] = 0;
                }
            }
        }
    }
    
    for (int i = 0; i < N; i++) {
        int box_x = fmin((int)(particles[i].x / dx), Nx - 1);
        int box_y = fmin((int)(particles[i].y / dy), Ny - 1);
        int box_z = fmin((int)(particles[i].z / dz), Nz - 1);
        
        if (dbis == NULL){
            d[box_x][box_y] += 4./3.*particles[i].rad*particles[i].rad*particles[i].rad*M_PI/(dx*dy*dz);
        }
        
        else{
            if (particles[i].type == 1) {
                d[box_x][box_y] +=  4./3.*particles[i].rad*particles[i].rad*particles[i].rad*M_PI/(dx*dy*dz);
            }
            else{
                dbis[box_x][box_y] += 4./3.*particles[i].rad*particles[i].rad*particles[i].rad*M_PI/(dx*dy*dz);
            }
        }
        if (tp != NULL){
            tp[box_x][box_y][box_z] += 0.5*(particles[i].vx*particles[i].vx + 
                                           particles[i].vy*particles[i].vy + 
                                           particles[i].vz*particles[i].vz);
        }
    }
    #else
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            d[i][j] = 0;
            if (dbis != NULL)
               dbis[i][j] = 0;
            if (tp != NULL){
                tp[i][j] = 0;
            }
        }
    }
    
    for (int i = 0; i < N; i++) {
        int box_x = fmin((int)(particles[i].x / dx), Nx - 1);
        int box_y = fmin((int)(particles[i].y / dy), Ny - 1);
        
        if (dbis == NULL){
            d[box_x][box_y] += particles[i].rad*particles[i].rad*M_PI/(dx*dy);
        }
        
        else{
            if (particles[i].type == 1) {
                d[box_x][box_y] += particles[i].rad*particles[i].rad*M_PI/(dx*dy);
            }
            else{
                dbis[box_x][box_y] += particles[i].rad*particles[i].rad*M_PI/(dx*dy);
            }
        }
        if (tp != NULL){
            tp[box_x][box_y] += 0.5*(particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy);
        }
    }
    #endif
}

void saveDensityCoarse(particle *particles){
    computeLocalConcentration(particles);
    
    #if THREE_D
    if (tp == NULL){
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nz; k++) {
                    fprintf(interface, "%.6lf ", d[i][j][k]);
                }
                fprintf(interface, "\n");
            }
            //fprintf(interface, "\n");
        }
    }
    else{
        // Calculate average density across z-dimension
        for (int i = 0; i < Nx; i++) {
            double density_avg = 0;
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nz; k++) {
                    density_avg += d[i][j][k];
                }
            }
            fprintf(interface, "%.6lf ", density_avg/(Ny*Nz));
        }
        fprintf(interface, "\n");

        // Calculate average temperature across z-dimension
        for (int i = 0; i < Nx; i++) {
            double temp_avg = 0;
            double count = 0;
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nz; k++) {
                    if (d[i][j][k] > 0) {  // Avoid division by zero
                        temp_avg += tp[i][j][k]/(d[i][j][k]/((4.0/3.0)*M_PI/(dx*dy*dz)));
                        count++;
                    }
                }
            }
            fprintf(interface, "%.6lf ", count > 0 ? temp_avg/count : 0);
        }
        fprintf(interface, "\n");
    }
    #else
    if (tp == NULL){
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                fprintf(interface, "%.6lf ", d[i][j]);
            }
            fprintf(interface, "\n");
        }
    }
    else{
        for (int i = 0; i < Nx; i++) {
            double sss = 0;
            for (int j = 0; j < Ny; j++) {
                sss += d[i][j];
            }
            fprintf(interface, "%.6lf ", sss/Ny);
        }
        fprintf(interface, "\n");

        for (int i = 0; i < Nx; i++) {
            double sss = 0;
            for (int j = 0; j < Ny; j++) {
                sss += tp[i][j]/(d[i][j]/(M_PI/(dx*dy)));
            }
            fprintf(interface, "%.6lf ", sss/Ny);
        }
        fprintf(interface, "\n");
    }
    #endif
}

void computeInterfacesPos(particle *particles, double tresh){
    computeLocalConcentration(particles);
    
    #if THREE_D
    int J = Nx/4;
    // Calculate interface positions for each y-z plane
    for (int i = 0; i < Ny; i++){
        for (int k = 0; k < Nz; k++) {
            int count = 0;
            for (int j = J; j < Nx; j++){
                if ((d[j][i][k] < tresh) && (dbis[j][i][k] > tresh)){
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
    }
    //fprintf(interface, "\n");
    
    for (int i = 0; i < Ny; i++){
        for (int k = 0; k < Nz; k++) {
            int count = 0;
            for (int j = J; j > -Nx/4; j--){
                int m = j;
                if (j < 0){
                    m = Nx + j;
                }
                if ((d[m][i][k] < tresh) && (dbis[m][i][k] > tresh)){
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
    //fprintf(interface, "\n");
    #else
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
    #endif
}

void freeInterface(){
    #if THREE_D
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            free(d[i][j]);
            if (dbis != NULL) {
                free(dbis[i][j]);
            }
            if (tp != NULL) {
                free(tp[i][j]);
            }
        }
        free(d[i]);
        if (dbis != NULL) {
            free(dbis[i]);
        }
        if (tp != NULL) {
            free(tp[i]);
        }
    }
    #else
    for (int i = 0; i < Nx; i++) {
        free(d[i]);
        if (dbis != NULL) {
            free(dbis[i]);
        }
        if (tp != NULL) {
            free(tp[i]);
        }
    }
    #endif
    free(d);
    if (dbis != NULL) {
        free(dbis);
    }
    if (tp != NULL) {
        free(tp);
    }
}