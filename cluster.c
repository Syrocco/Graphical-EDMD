#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>


typedef struct cluster {
    particle** particles;
    int size;
    int capacity;
} cluster;

bool areParticlesClose(particle* p1, particle* p2, double r_c);
void addParticleToCluster(cluster* cl, particle* p);
void findClusters(particle* particles, int N, double r_c, cluster** clusters, int* numClusters);
void freeClusters(cluster* clusters, int numClusters);


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

void findClusters(particle* particles, int N, double r_c, cluster** clusters, int* numClusters) {
    bool* visited = calloc(N, sizeof(bool));
    *numClusters = 0;
    *clusters = NULL;

    for (int i = 0; i < N; i++) {
        if (!visited[i]) {
            cluster cl;
            cl.size = 0;
            cl.capacity = 10;
            cl.particles = malloc(cl.capacity * sizeof(particle*));
            addParticleToCluster(&cl, &particles[i]);
            visited[i] = true;

            for (int j = 0; j < cl.size; j++) {
                particle* p1 = cl.particles[j];
                for (int k = 0; k < N; k++) {
                    if (!visited[k] && areParticlesClose(p1, &particles[k], r_c)) {
                        addParticleToCluster(&cl, &particles[k]);
                        visited[k] = true;
                    }
                }
            }

            *clusters = realloc(*clusters, (*numClusters + 1) * sizeof(cluster));
            (*clusters)[*numClusters] = cl;
            (*numClusters)++;
        }
    }

    free(visited);
}

void freeClusters(cluster* clusters, int numClusters) {
    for (int i = 0; i < numClusters; i++) {
        free(clusters[i].particles);
    }
    free(clusters);
}