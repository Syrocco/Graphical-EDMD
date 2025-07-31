#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

int nqx, nqy, N;
double *qx, *qy;
double *structFactor;
double complex *structFactorComplex;
FILE *fileout;

void initStructureFactor(double qmax, double Lx, double Ly, int num, FILE* file, double doFQT){
    double xx = 2*M_PI/Lx;
    double yy = 2*M_PI/Ly;
    nqx = 2*qmax/xx + 1;
    nqy = 2*qmax/yy + 1;
    qx = calloc(nqx, sizeof(double));
    qy = calloc(nqy, sizeof(double));
    
    N = num;
    for (int i = 0; i < nqx; ++i){
        qx[i] = xx*(i - (nqx - 1)/2);
    }
    for (int i = 0; i < nqy; ++i){
        qy[i] = yy*(i - (nqy - 1)/2);
    }


    fileout = file;
    fprintf(file, "-%lf \n", doFQT);
    for (int i = 0; i < nqx; ++i)
        fprintf(fileout, "%g ", qx[i]);
    fprintf(fileout, "\n");

    for (int j = 0; j < nqy; ++j)
        fprintf(fileout, "%g ", qy[j]);
    fprintf(fileout, "\n");
    
    if (doFQT)
        structFactorComplex = calloc(nqy*nqx, sizeof(double complex));
    else
        structFactor = calloc(nqy*nqx, sizeof(double));
}

void computeStructureFactor(particle* particles){
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nqx; i++){
        for (int j = 0; j < nqy; j++){
            double im = 0;
            double re = 0;

            for (int n = 0; n < N; n++){
                particle* p = particles + n;
                double qr = qx[i]*p->x + qy[j]*p->y;

                re += cos(qr);
                im += sin(qr);
            }
            if (structFactorComplex != NULL)    
                structFactorComplex[i*nqy + j] = re + I*im;
            else
                structFactor[i*nqy + j] = (re*re + im*im)/N;
        }
    }
}

void saveStructureFactor(particle* particles){
    computeStructureFactor(particles);

    for (int i = 0; i < nqx; ++i){
        for (int j = 0; j < nqy; ++j){
            if (structFactorComplex != NULL)
                fprintf(fileout, "(%g+%gj) ", creal(structFactorComplex[i*nqy + j]), cimag(structFactorComplex[i*nqy + j]));
            else
                fprintf(fileout, "%g ", structFactor[i*nqy + j]);
        }
        fprintf(fileout, "\n");
    }
}