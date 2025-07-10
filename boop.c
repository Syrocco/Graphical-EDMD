#include "EDMD.h"
#include <math.h>
#include <complex.h>
#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"
#include <stdlib.h>
#include <string.h>
#include "boop.h"





boop_data* computeBOOPVoronoi(particle* particles, int N, double Lx, double Ly) {
    // First pass: count how many image points we need
    int image_count = 0;
    double boundary_distance = 6.0;
    
    for (int i = 0; i < N; i++) {
        double x = particles[i].x;
        double y = particles[i].y;
        
        // Check proximity to boundaries and count needed images
        if (x < boundary_distance) image_count++; // +Lx image
        if (x > Lx - boundary_distance) image_count++; // -Lx image
        if (y < boundary_distance) image_count++; // +Ly image
        if (y > Ly - boundary_distance) image_count++; // -Ly image
        
        // Corner images
        if (x < boundary_distance && y < boundary_distance) image_count++; // +Lx, +Ly
        if (x > Lx - boundary_distance && y < boundary_distance) image_count++; // -Lx, +Ly
        if (x < boundary_distance && y > Ly - boundary_distance) image_count++; // +Lx, -Ly
        if (x > Lx - boundary_distance && y > Ly - boundary_distance) image_count++; // -Lx, -Ly
    }
    
    int Nall = N + image_count;
    jcv_point* points = malloc(Nall * sizeof(jcv_point));
    
    // Add original points
    for (int i = 0; i < N; i++) {
        points[i].x = particles[i].x;
        points[i].y = particles[i].y;
    }
    
    // Add boundary images only when needed
    int point_index = N;
    for (int i = 0; i < N; i++) {
        double x = particles[i].x;
        double y = particles[i].y;
        
        // Add x-boundary images
        if (x < boundary_distance) {
            points[point_index].x = x + Lx;
            points[point_index].y = y;
            point_index++;
        }
        if (x > Lx - boundary_distance) {
            points[point_index].x = x - Lx;
            points[point_index].y = y;
            point_index++;
        }
        
        // Add y-boundary images
        if (y < boundary_distance) {
            points[point_index].x = x;
            points[point_index].y = y + Ly;
            point_index++;
        }
        if (y > Ly - boundary_distance) {
            points[point_index].x = x;
            points[point_index].y = y - Ly;
            point_index++;
        }
        
        // Add corner images
        if (x < boundary_distance && y < boundary_distance) {
            points[point_index].x = x + Lx;
            points[point_index].y = y + Ly;
            point_index++;
        }
        if (x > Lx - boundary_distance && y < boundary_distance) {
            points[point_index].x = x - Lx;
            points[point_index].y = y + Ly;
            point_index++;
        }
        if (x < boundary_distance && y > Ly - boundary_distance) {
            points[point_index].x = x + Lx;
            points[point_index].y = y - Ly;
            point_index++;
        }
        if (x > Lx - boundary_distance && y > Ly - boundary_distance) {
            points[point_index].x = x - Lx;
            points[point_index].y = y - Ly;
            point_index++;
        }
    }
    
    jcv_diagram diagram;
    memset(&diagram, 0, sizeof(jcv_diagram));
    jcv_diagram_generate(Nall, points, NULL, NULL, &diagram);
    boop_data* boop = calloc(N, sizeof(boop_data));
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
    free(points);
    return boop;
}

boop_data* computeBOOPCutoff(particle* particles, int N, double r_c, particle** cellList, int Nxcells) {
    boop_data* boop = calloc(N, sizeof(boop_data));
    for (int i = 0; i < N; i++) {
        particle* p1 = &particles[i];
        int X = p1->cell[0];
        int Y = p1->cell[1];
        double complex sum5 = 0.0 + 0.0*I;
        double complex sum6 = 0.0 + 0.0*I;
        double complex sum7 = 0.0 + 0.0*I;
        int neighbors = 0;
        // Loop over neighboring cells (including the current cell)
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                particle* p2 = cellList[PBCcellY(Y + j) * Nxcells + PBCcellX(X + k)];
                while (p2 != NULL) {
                    if (p2->num != p1->num) {
                        double dx = p2->x - p1->x;
                        double dy = p2->y - p1->y;
                        PBC(&dx, &dy);
                        double r2 = dx*dx + dy*dy;
                        if (r2 < r_c*r_c) {
                            double theta = atan2(dy, dx);
                            sum5 += cexp(5 * theta * I);
                            sum6 += cexp(6 * theta * I);
                            sum7 += cexp(7 * theta * I);
                            neighbors++;
                        }
                    }
                    p2 = p2->nxt;
                }
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

