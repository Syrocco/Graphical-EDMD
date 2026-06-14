#include "EDMD.h"
#include "osmosis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double DX;
double areaSlice;
int sliceNum;
osmosisSlice* data;
FILE* osmosis_file;

void init_osmosis(double lx, double ly, double sigma, FILE* inter){
    sliceNum = (int)(lx / sigma);
    DX = lx / sliceNum;


    areaSlice = ly*DX;
    data = calloc(sliceNum, sizeof(osmosisSlice));
    for (int i = 0; i < sliceNum; i++){
        data[i].x = (i + 0.5)*DX;
        data[i].densitySolvent = 0;
        data[i].densitySolute = 0;
        data[i].temperatureSolvent = 0;
        data[i].temperatureSolute = 0;
    }

    osmosis_file = inter;
    for (int i = 0; i < sliceNum; i++){
        fprintf(osmosis_file, "%.3lf ", data[i].x);
    }
    fprintf(osmosis_file, "\n");
    fflush(osmosis_file);
}

void save_osmosis(particle* particles, int N){
    for (int i = 0; i < sliceNum; i++){
        data[i].densitySolvent = 0;
        data[i].densitySolute = 0;
        data[i].temperatureSolvent = 0;
        data[i].temperatureSolute = 0;
    }
    
    double solventRad = 0; 
    double soluteRad = 0;

    for (int i = 0; i < N; i++){
        int box = fmin((int)(particles[i].x/DX), sliceNum - 1);
        if (particles[i].type == 0){
            solventRad = particles[i].rad;
            data[box].densitySolvent += 1;
            data[box].temperatureSolvent += 0.5*particles[i].m*(particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy);
        }
        else{
            soluteRad = particles[i].rad;
            data[box].densitySolute += 1;
            data[box].temperatureSolute += 0.5*particles[i].m*(particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy);
        }
    }
    
    for (int i = 0; i < sliceNum; i++){
        if (data[i].densitySolvent > 0){
            data[i].temperatureSolvent /= data[i].densitySolvent; 
        }
        if (data[i].densitySolute > 0){
            data[i].temperatureSolute /= data[i].densitySolute;
        }
        data[i].densitySolvent *= M_PI*solventRad*solventRad/areaSlice;
        data[i].densitySolute *= M_PI*soluteRad*soluteRad/areaSlice;
    }

    for (int i = 0; i < sliceNum; i++){
        fprintf(osmosis_file, "%.4lf ", data[i].densitySolvent);
    }
    fprintf(osmosis_file, "\n");
    for (int i = 0; i < sliceNum; i++){
        fprintf(osmosis_file, "%.4lf ", data[i].densitySolute);
    }
    fprintf(osmosis_file, "\n");
    for (int i = 0; i < sliceNum; i++){
        fprintf(osmosis_file, "%.4lf ", data[i].temperatureSolvent);
    }
    fprintf(osmosis_file, "\n");
    for (int i = 0; i < sliceNum; i++){
        fprintf(osmosis_file, "%.4lf ", data[i].temperatureSolute);
    }
    fprintf(osmosis_file, "\n");
    fflush(osmosis_file);
}