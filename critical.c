#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>


#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

double getCOMx(particle* particles, double N, double Lx);
void countConditions(particle* particles, int N, double Lx, double Ly, int* counts, double* energyCounts, double dfrac);

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

void countConditions(particle* particles, int N, double Lx, double Ly, int* counts, double* energyCounts, double dfrac){

    double COMx = getCOMx(particles, N, Lx);

    for (int i = 0; i < 4; ++i) {
        counts[i] = 0;
        energyCounts[i] = 0;
    }
    particle* p;
    for (int i = 0; i < N; ++i) {
        p = particles + i;
        double X = particles[i].x + (Lx / 4 - COMx);
        if (X < 0) {
            X += Lx;
        }
        else if (X > Lx){
            X -= Lx;
        }

        double left = Lx/4;
        double right = 3*Lx/4;

        int firstx = ( (X > (left - dfrac)) && (X < (left + dfrac)) );
        int seconx = ( (X > (right - dfrac)) && (X < (right + dfrac)) );
        int firsty = (p->y < (Ly/2));
        int secony = (p->y > (Ly/2));

        if (firstx && firsty) {
            counts[0]++;
            energyCounts[0] += 0.5*p->m*(p->vx*p->vx + p->vy*p->vy);
        }
        if (firstx && secony) {
            counts[1]++;
            energyCounts[1] += 0.5*p->m*(p->vx*p->vx + p->vy*p->vy);
        }
        if (seconx && firsty) {
            counts[2]++;
            energyCounts[2] += 0.5*p->m*(p->vx*p->vx + p->vy*p->vy);
        }
        if (seconx && secony) {
            counts[3]++;
            energyCounts[3] += 0.5*p->m*(p->vx*p->vx + p->vy*p->vy);
        }
    }
    for (int i = 0; i < 4; i++){
        energyCounts[i] /= counts[i];
    }
}