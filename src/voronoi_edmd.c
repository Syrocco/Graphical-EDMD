#include "voronoi_edmd.h"
#include "jc_voronoi.h"
#include "EDMD.h"
#include <string.h>
#include <stdlib.h>


inline jcv_real jcv_cross(const jcv_point a, const jcv_point b){
    return a.x*b.y - a.y*b.x;
}

inline double jcv_perimeter(const jcv_site* a){
    double perimeter = 0;
    const jcv_graphedge* e = a->edges;
    while (e){
        perimeter += jcv_point_dist(&e->pos[0], &e->pos[1]);
        e = e->next;
        
    }
    return perimeter;
}

inline double jcv_area(const jcv_site* a){
    double area = 0;
    const jcv_graphedge* e = a->edges;
    while (e){
        area += jcv_cross(e->pos[0], e->pos[1]);
        e = e->next;
    }
    return 0.5*jcv_abs(area);
}

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

double* get_particle_voronoi_area(particle* particles, int N, double Lx, double Ly) {
    jcv_diagram diagram = get_particle_voronoi(particles, N, Lx, Ly);
    double* areas = malloc(N * sizeof(double));
    const jcv_site* sites = jcv_diagram_get_sites(&diagram);
    for (int i = 0; i < diagram.numsites; i++) {
        int index = sites[i].index;
        if (index >= N) continue;
        areas[index] = jcv_area(&sites[i]);
    }
    jcv_diagram_free(&diagram);
    return areas;
}


double* get_particle_voronoi_perimeter(particle* particles, int N, double Lx, double Ly) {
    jcv_diagram diagram = get_particle_voronoi(particles, N, Lx, Ly);
    double* perimeters = malloc(N * sizeof(double));
    const jcv_site* sites = jcv_diagram_get_sites(&diagram);
    for (int i = 0; i < diagram.numsites; i++) {
        int index = sites[i].index;
        if (index >= N) continue;
        perimeters[index] = jcv_perimeter(&sites[i]);
    }

    jcv_diagram_free(&diagram);
    return perimeters;
}