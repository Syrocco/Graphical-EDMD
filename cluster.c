#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>



typedef struct cluster {
    particle** particles;
    int size;
    int capacity;
} cluster;

int numClusters;
cluster* clusters;

bool areParticlesClose(particle* p1, particle* p2, double r_c);
void addParticleToCluster(cluster* cl, particle* p);
void findClusters(particle* particles, int N, double r_c, particle** cellList, int Nxcells);
void freeClusters();
int compareClusters(const void* a, const void* b);
void sortClustersBySize();



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

void findClusters(particle* particles, int N, double r_c, particle** cellList, int Nxcells) {
    bool* visited = calloc(N, sizeof(bool));
    numClusters = 0;
    clusters = NULL;

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

                for (int j = -1; j <= 1; j++){
                    for (int k = -1; k <= 1; k++){
                        #if THREE_D
                        for (int l = -1; l <= 1; l++){
                            particle* p2 = cellList[PBCcellZ(Z + k) * Nycells * Nxcells + PBCcellY(Y + j) * Nxcells + PBCcellX(X + l)];
                        #else
                            particle* p2 = cellList[PBCcellY(Y + j) * Nxcells + PBCcellX(X + k)];
                        #endif
                        while (p2 != NULL){ //while there is a particle in the doubly linked list of the cellList do...
                            if (p1->num != p2->num){
                                if (!visited[p2->num] && areParticlesClose(p1, p2, r_c)) {
                                    addParticleToCluster(&cl, p2);
                                    visited[p2->num] = true;
                                }
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

void freeClusters() {
    for (int i = 0; i < numClusters; i++) {
        free(clusters[i].particles);
    }
    free(clusters);
}

