#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

double getCOMx(particle* particles, double N, double Lx);
void countConditions(particle* particles, int N, double Lx, double Ly, int* counts);

double getCOMx(particle* particles, double N, double Lx){
    double* xi = (double*)malloc(N * sizeof(double));
    double* zi = (double*)malloc(N * sizeof(double));

    // Constants
    double ri = Lx / (2 * M_PI);

    // Calculate thetai, xi, and zi
    for (int i = 0; i < N; ++i) {
        double thetai = particles[i].x/Lx*2*M_PI;
        xi[i] = ri*cos(thetai);
        zi[i] = ri*sin(thetai);
    }

    // Compute mean of xi and zi
    double mean_xi = 0.0;
    double mean_zi = 0.0;
    for (int i = 0; i < N; ++i) {
        mean_xi += xi[i];
        mean_zi += zi[i];
    }
    mean_xi /= N;
    mean_zi /= N;

    free(xi);
    free(zi);

    return Lx/(2*M_PI)*(atan2(-mean_zi, -mean_xi) + M_PI);
}

void countConditions(particle* particles, int N, double Lx, double Ly, int* counts){

    double COMx = getCOMx(particles, N, Lx);

    for (int i = 0; i < 4; ++i) {
        counts[i] = 0;
    }

    for (int i = 0; i < N; ++i) {
        double X = particles[i].x + (Lx / 4 - COMx);
        if (X < 0) {
            X += Lx;
        }
        else if (X > Lx){
            X -= Lx;
        }
        int firstx = ((X > (Lx/6)) && (X < (2*Lx/6)));
        int seconx = ((X > (4*Lx/6)) && (X < (5*Lx/6)));
        int firsty = (particles[i].y < (Ly/2));
        int secony = (particles[i].y > (Ly/2));

        if (firstx && firsty) {
            counts[0]++;
        }
        if (firstx && secony) {
            counts[1]++;
        }
        if (seconx && firsty) {
            counts[2]++;
        }
        if (seconx && secony) {
            counts[3]++;
        }
    }
}