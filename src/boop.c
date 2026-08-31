#include "EDMD.h"
#include <math.h>
#include <complex.h>
#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"
#include "voronoi_edmd.h"
#include <stdlib.h>
#include <string.h>
#include "boop.h"

extern int Nycells;
extern double cellxSize;
extern double cellySize;
#if THREE_D
extern int Nzcells;
extern double cellzSize;
#endif



/* BOOP uses a strict distance cutoff.  With cell width h, an index offset d
   can contribute iff (d - 1)*h < r_c, giving ceil(r_c/h) cells per axis. */
static int boopCellRange(double r_c, double cellSize) {
    return (int)ceil(fabs(r_c)/cellSize);
}

static int boopNeighborCellCount(int cellCount, int range) {
    return (range >= cellCount/2) ? cellCount : 2*range + 1;
}

static int boopNeighborCellCoordinate(int center, int cellCount,
                                      int range, int index) {
    if (range >= cellCount/2) {
        return index;
    }
    int coordinate = center - range + index;
    if (coordinate < 0) coordinate += cellCount;
    if (coordinate >= cellCount) coordinate -= cellCount;
    return coordinate;
}




boop_data* computeBOOPVoronoi(particle* particles, int N, double Lx, double Ly) {
    /* The Voronoi implementation is two-dimensional.  In a 3D build this
       observable intentionally describes the particles' XY projection. */
    boop_data* boop = calloc(N, sizeof(boop_data));
    jcv_diagram diagram = get_particle_voronoi(particles, N, Lx, Ly);
    const jcv_site* sites = jcv_diagram_get_sites(&diagram);
    for (int i = 0; i < diagram.numsites; i++) {
        int index = sites[i].index;
        if (index >= N) continue;
        double complex sum5 = 0.0 + 0.0*I;
        double complex sum6 = 0.0 + 0.0*I;
        double complex sum7 = 0.0 + 0.0*I;
        const jcv_graphedge* e = sites[i].edges;
        while (e) {
            if (e->neighbor) {
                double x2 = e->neighbor->p.x;
                double y2 = e->neighbor->p.y;
                double dx = x2 - particles[index].x;
                double dy = y2 - particles[index].y;

                double theta = atan2(dy, dx);
                
                
                sum5 += cexp(5 * theta * I);
                sum6 += cexp(6 * theta * I);
                sum7 += cexp(7 * theta * I);

                boop[index].neighbors++;
                }
            else {
                printf("no neighbor %d ???\n", i);
            }
            e = e->next;
        }
        if (boop[index].neighbors > 0) {
            boop[index].q5 = cabs(sum5) / boop[index].neighbors;
            boop[index].q6 = cabs(sum6) / boop[index].neighbors;
            boop[index].q7 = cabs(sum7) / boop[index].neighbors;
            boop[index].q6_arg = carg(sum6);
        } 
    }
    // Cleanup
    jcv_diagram_free(&diagram);
    
    return boop;
}

boop_data* computeBOOPCutoff(particle* particles, int N, double r_c, particle** cellList, int Nxcells) {
    boop_data* boop = calloc(N, sizeof(boop_data));
    const int rangeX = boopCellRange(r_c, cellxSize);
    const int rangeY = boopCellRange(r_c, cellySize);
    const int numberX = boopNeighborCellCount(Nxcells, rangeX);
    const int numberY = boopNeighborCellCount(Nycells, rangeY);
    #if THREE_D
    const int rangeZ = boopCellRange(r_c, cellzSize);
    const int numberZ = boopNeighborCellCount(Nzcells, rangeZ);
    #endif
    for (int i = 0; i < N; i++) {
        particle* p1 = &particles[i];
        int X = p1->cell[0];
        int Y = p1->cell[1];
        #if THREE_D
        int Z = p1->cell[2];
        #endif
        double complex sum5 = 0.0 + 0.0*I;
        double complex sum6 = 0.0 + 0.0*I;
        double complex sum7 = 0.0 + 0.0*I;
        int neighbors = 0;
        // Loop over neighboring cells (including the current cell)
        for (int yi = 0; yi < numberY; yi++) {
            const int neighborY = boopNeighborCellCoordinate(Y, Nycells, rangeY, yi);
            for (int xi = 0; xi < numberX; xi++) {
                const int neighborX = boopNeighborCellCoordinate(X, Nxcells, rangeX, xi);
                #if THREE_D
                for (int zi = 0; zi < numberZ; zi++) {
                    const int neighborZ = boopNeighborCellCoordinate(Z, Nzcells, rangeZ, zi);
                    particle* p2 = cellList[
                        neighborZ * Nycells * Nxcells +
                        neighborY * Nxcells + neighborX];
                #else
                particle* p2 = cellList[neighborY * Nxcells + neighborX];
                #endif
                while (p2 != NULL) {
                    if (p2->num != p1->num) {
                        double dx = p2->x - p1->x;
                        double dy = p2->y - p1->y;
                        #if THREE_D
                        double dz = p2->z - p1->z;
                        PBC(&dx, &dy, &dz);
                        double r2 = dx*dx + dy*dy + dz*dz;
                        #else
                        PBC(&dx, &dy);
                        double r2 = dx*dx + dy*dy;
                        #endif
                        if (r2 < r_c*r_c) {
                            /* q5/q6/q7 here are planar angular harmonics.
                               In 3D, select neighbors by their true 3D
                               separation but measure projected XY order. */
                            double theta = atan2(dy, dx);
                            sum5 += cexp(5 * theta * I);
                            sum6 += cexp(6 * theta * I);
                            sum7 += cexp(7 * theta * I);
                            neighbors++;
                        }
                    }
                    p2 = p2->nxt;
                }
                #if THREE_D
                }
                #endif
            }
        }
        boop[i].neighbors = neighbors;
        if (neighbors > 0) {
            boop[i].q5 = cabs(sum5) / neighbors;
            boop[i].q6 = cabs(sum6) / neighbors;
            boop[i].q7 = cabs(sum7) / neighbors;
            boop[i].q6_arg = carg(sum6);
        } else {
            boop[i].q5 = 0.0;
            boop[i].q6 = 0.0;
            boop[i].q7 = 0.0;
            boop[i].q6_arg = 0.0;
        }
    }
    return boop;
}
