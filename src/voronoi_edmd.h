#ifndef VORONOI_EDMD_H
#define VORONOI_EDMD_H

#include "jc_voronoi.h"
#include "EDMD.h"

jcv_diagram get_particle_voronoi(particle* particles, int N, double Lx, double Ly);

#endif // VORONOI_EDMD_H