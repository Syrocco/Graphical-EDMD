#include <stdio.h>
#include "EDMD.h"

typedef struct osmosisSlice osmosisSlice;
struct osmosisSlice {
	double x;
    double densitySolvent;
    double densitySolute;
    double temperatureSolute;
    double temperatureSolvent;
};

void init_osmosis(double lx, double ly, double sigma, FILE* inter);
void save_osmosis(particle* particles, int N);