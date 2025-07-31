#include "voronoi_edmd.h"
#include "jc_voronoi.h"
#include "EDMD.h"
#include <string.h>
#include <stdlib.h>



jcv_diagram get_particle_voronoi(particle* particles, int N, double Lx, double Ly) {
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
    free(points);
    return diagram;
}