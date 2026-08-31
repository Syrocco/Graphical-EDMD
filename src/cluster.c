#include "EDMD.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



typedef struct cluster {
    particle** particles;
    int size;
    int capacity;
} cluster;

int numClusters;
cluster* clusters;

extern int clusterMode;
extern int clusterNeighborThreshold;
extern double cellxSize;
extern double cellySize;
#if THREE_D
extern int Nzcells;
extern double cellzSize;
#endif

bool areParticlesClose(particle* p1, particle* p2, double r_c);
void addParticleToCluster(cluster* cl, particle* p);
void findClusters(particle* particles, int N, double r_c, particle** cellList, int Nxcells, int Nycells);
void findClustersOld(particle* particles, int N, double r_c, particle** cellList, int Nxcells, int Nycells);
void freeClusters();
int compareClusters(const void* a, const void* b);
void sortClustersBySize();



/* Cells whose indices differ by d can contain points separated by as little
   as (d - 1)*cellSize.  Clustering includes particles exactly at r_c, so the
   closed-cutoff stencil needs floor(r_c/cellSize) + 1 cells in each axis. */
static int clusterCellRange(double r_c, double cellSize) {
    return (int)floor(fabs(r_c)/cellSize) + 1;
}

/* If a periodic dimension is shorter than the requested stencil, traverse
   every cell once.  Re-visiting wrapped cells would inflate neighborCount. */
static int clusterNeighborCellCount(int cellCount, int range) {
    return (range >= cellCount/2) ? cellCount : 2*range + 1;
}

static int clusterNeighborCellCoordinate(int center, int cellCount,
                                         int range, int index) {
    if (range >= cellCount/2) {
        return index;
    }
    int coordinate = center - range + index;
    if (coordinate < 0) coordinate += cellCount;
    if (coordinate >= cellCount) coordinate -= cellCount;
    return coordinate;
}



bool areParticlesClose(particle* p1, particle* p2, double r_c) {
    double dx = p2->x - p1->x;
    double dy = p2->y - p1->y;
    #if THREE_D
    double dz = p2->z - p1->z;
    PBC(&dx, &dy, &dz);
    return (dx*dx + dy*dy + dz*dz) <= (r_c*r_c);
    #else
    PBC(&dx, &dy);
    return (dx*dx + dy*dy) <= (r_c*r_c);
    #endif
}

void addParticleToCluster(cluster* cl, particle* p) {
    if (cl->size == cl->capacity) {
        cl->capacity *= 2;
        cl->particles = realloc(cl->particles, cl->capacity * sizeof(particle*));
    }
    cl->particles[cl->size++] = p;
}

int compareClusters(const void* a, const void* b){
    cluster* cl1 = (cluster*)a;
    cluster* cl2 = (cluster*)b;
    return cl2->size - cl1->size; // Sort in descending order
}

void sortClustersBySize(){
    qsort(clusters, numClusters, sizeof(cluster), compareClusters);
}



void findClustersOld(particle* particles, int N, double r_c, particle** cellList, int Nxcells, __attribute__((unused)) int Nycells) {
    bool* visited = calloc(N, sizeof(bool));
    numClusters = 0;
    clusters = NULL;
    const int rangeX = clusterCellRange(r_c, cellxSize);
    const int rangeY = clusterCellRange(r_c, cellySize);
    const int numberX = clusterNeighborCellCount(Nxcells, rangeX);
    const int numberY = clusterNeighborCellCount(Nycells, rangeY);
    #if THREE_D
    const int rangeZ = clusterCellRange(r_c, cellzSize);
    const int numberZ = clusterNeighborCellCount(Nzcells, rangeZ);
    #endif
    for (int i = 0; i < N; i++) {
        if (!visited[i]) {
            cluster cl;
            cl.size = 0;
            cl.capacity = 10;
            cl.particles = malloc(cl.capacity * sizeof(particle*));
            addParticleToCluster(&cl, &particles[i]);
            visited[i] = true;

            for (int j = 0; j < cl.size; j++){
                particle* p1 = cl.particles[j];
                int X = p1->cell[0];
                int Y = p1->cell[1];
                #if THREE_D
                int Z = p1->cell[2];
                #endif
                for (int yi = 0; yi < numberY; yi++){
                    const int neighborY = clusterNeighborCellCoordinate(Y, Nycells, rangeY, yi);
                    for (int xi = 0; xi < numberX; xi++){
                        const int neighborX = clusterNeighborCellCoordinate(X, Nxcells, rangeX, xi);
                        #if THREE_D
                        for (int zi = 0; zi < numberZ; zi++){
                            const int neighborZ = clusterNeighborCellCoordinate(Z, Nzcells, rangeZ, zi);
                            particle* p2 = cellList[neighborZ * Nycells * Nxcells + neighborY * Nxcells + neighborX];
                        #else
                            particle* p2 = cellList[neighborY * Nxcells + neighborX];
                        #endif
                        while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
                            if (!visited[p2->num] && areParticlesClose(p1, p2, r_c)) {
                                addParticleToCluster(&cl, p2);
                                visited[p2->num] = true;
                            }
                            p2 = p2->nxt;
                        }
                    #if THREE_D
                    }
                    #endif
                    }
                }
            }
            clusters = realloc(clusters, (numClusters + 1) * sizeof(cluster));
            clusters[numClusters] = cl;
            numClusters++;
            
        }
    }
    sortClustersBySize();
    free(visited);
}

void findClusters(particle* particles, int N, double r_c, particle** cellList, int Nxcells, __attribute__((unused)) int Nycells) {
    const int neighborThreshold = clusterNeighborThreshold;
    bool* visited = calloc((unsigned int)N, sizeof(bool));
    bool* targetLike = calloc((unsigned int)N, sizeof(bool));
    numClusters = 0;
    clusters = NULL;
    const int rangeX = clusterCellRange(r_c, cellxSize);
    const int rangeY = clusterCellRange(r_c, cellySize);
    const int numberX = clusterNeighborCellCount(Nxcells, rangeX);
    const int numberY = clusterNeighborCellCount(Nycells, rangeY);
    #if THREE_D
    const int rangeZ = clusterCellRange(r_c, cellzSize);
    const int numberZ = clusterNeighborCellCount(Nzcells, rangeZ);
    #endif

    // First pass: identify liquidlike particles
    for (int i = 0; i < N; i++) {
        int neighborCount = 0;
        particle* p1 = &particles[i];
        int X = p1->cell[0];
        int Y = p1->cell[1];
        #if THREE_D
        int Z = p1->cell[2];
        #endif
        for (int yi = 0; yi < numberY; yi++) {
            const int neighborY = clusterNeighborCellCoordinate(Y, Nycells, rangeY, yi);
            for (int xi = 0; xi < numberX; xi++) {
                const int neighborX = clusterNeighborCellCoordinate(X, Nxcells, rangeX, xi);
                #if THREE_D
                for (int zi = 0; zi < numberZ; zi++) {
                    const int neighborZ = clusterNeighborCellCoordinate(Z, Nzcells, rangeZ, zi);
                    particle* p2 = cellList[neighborZ * Nycells * Nxcells + neighborY * Nxcells + neighborX];
                #else
                    particle* p2 = cellList[neighborY * Nxcells + neighborX];
                #endif
                while (p2 != NULL) {
                    if (p1->num != p2->num && areParticlesClose(p1, p2, r_c)) {
                        neighborCount++;
                    }
                    p2 = p2->nxt;
                }
                #if THREE_D
                }
                #endif
            }
        }
        if (clusterMode == 0) {
            if (neighborCount >= neighborThreshold) {
                targetLike[i] = true;
            }
        } else {
            if (neighborCount <= neighborThreshold) {
                targetLike[i] = true;
            }
        }
    }

    // Second pass: find clusters among liquidlike particles
    for (int i = 0; i < N; i++) {
        if (targetLike[i] && !visited[i]) {
            cluster cl;
            cl.size = 0;
            cl.capacity = 10;
            cl.particles = malloc(cl.capacity * sizeof(particle*));
            addParticleToCluster(&cl, &particles[i]);
            visited[i] = true;

            for (int j = 0; j < cl.size; j++) {
                particle* p1 = cl.particles[j];
                int X = p1->cell[0];
                int Y = p1->cell[1];
                #if THREE_D
                int Z = p1->cell[2];
                #endif
                for (int yi = 0; yi < numberY; yi++) {
                    const int neighborY = clusterNeighborCellCoordinate(Y, Nycells, rangeY, yi);
                    for (int xi = 0; xi < numberX; xi++) {
                        const int neighborX = clusterNeighborCellCoordinate(X, Nxcells, rangeX, xi);
                        #if THREE_D
                        for (int zi = 0; zi < numberZ; zi++) {
                            const int neighborZ = clusterNeighborCellCoordinate(Z, Nzcells, rangeZ, zi);
                            particle* p2 = cellList[neighborZ * Nycells * Nxcells + neighborY * Nxcells + neighborX];
                        #else
                            particle* p2 = cellList[neighborY * Nxcells + neighborX];
                        #endif
                        while (p2 != NULL) {
                            if (targetLike[p2->num] && !visited[p2->num] && areParticlesClose(p1, p2, r_c)) {
                                addParticleToCluster(&cl, p2);
                                visited[p2->num] = true;
                            }
                            p2 = p2->nxt;
                        }
                        #if THREE_D
                        }
                        #endif
                    }
                }
            }
            clusters = realloc(clusters, (numClusters + 1) * sizeof(cluster));
            clusters[numClusters] = cl;
            numClusters++;
        }
    }
    for (int i = 0; i < N; i++) {
        if (!targetLike[i]) {
            cluster cl;
            cl.size = 1;
            cl.capacity = 1;
            cl.particles = malloc(sizeof(particle*));
            cl.particles[0] = &particles[i];
            clusters = realloc(clusters, (numClusters + 1) * sizeof(cluster));
            clusters[numClusters] = cl;
            numClusters++;
        }
    }
    sortClustersBySize();
    free(visited);
    free(targetLike);
}

void freeClusters() {
    for (int i = 0; i < numClusters; i++) {
        free(clusters[i].particles);
    }
    free(clusters);
}


void saveHistrogramCluster(FILE* fileout){
    for (int i = 0; i < numClusters; i++) {
        if (clusters[i].size < 10){
            return;
        }
        fprintf(fileout, "%d ", clusters[i].size);
    }
}
