#ifndef BOOP_H
#define BOOP_H

#include "EDMD.h"

typedef struct {
    double q6, q5, q7; // Order parameters
    double q6_arg; // Normalized order parameters
    int neighbors; // Number of neighbors
} boop_data;

boop_data* computeBOOPVoronoi(particle* particles, int N, double Lx, double Ly);
boop_data* computeBOOPCutoff(particle* particles, int N, double r, particle** cellList, int Nxcells);

#endif /* BOOP_H */
