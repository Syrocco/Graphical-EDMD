#include "EDMD.h"
#include "intruder.h"
#include <math.h>
#include <stdio.h>




polygon square(double side_length){
    polygon p;
    p.n_vertices = 4;
    p.vertices[0][0] = -side_length/2;
    p.vertices[0][1] = -side_length/2;
    p.vertices[1][0] = side_length/2;
    p.vertices[1][1] = -side_length/2;
    p.vertices[2][0] = side_length/2;
    p.vertices[2][1] = side_length/2;
    p.vertices[3][0] = -side_length/2;
    p.vertices[3][1] = side_length/2;
    
    for(int i = 0; i < 4; i++){
        int j = (i+1)%4;
        double dx = p.vertices[j][0] - p.vertices[i][0];
        double dy = p.vertices[j][1] - p.vertices[i][1];
        double length = sqrt(dx*dx + dy*dy);
        p.lengths[i] = length;
        if(length > 0){
            p.normals[i][0] = -dy/length;
            p.normals[i][1] = dx/length;
            p.tangents[i][0] = dx/length;
            p.tangents[i][1] = dy/length;
        }
        else{
            p.normals[i][0] = 0;
            p.normals[i][1] = 0;
        }
    }
    
    return p;
}

polygon triangle(double base_length, double height){
    polygon p;
    p.n_vertices = 3;
    p.vertices[0][0] = -base_length/2;
    p.vertices[0][1] = -height/2;
    p.vertices[1][0] = base_length/2;
    p.vertices[1][1] = -height/2;
    p.vertices[2][0] = 0;
    p.vertices[2][1] = height/2;
    
    double x_com = 0;
    double y_com = 0;
    for(int i = 0; i < 3; i++){
        int j = (i+1)%3;
        double dx = p.vertices[j][0] - p.vertices[i][0];
        double dy = p.vertices[j][1] - p.vertices[i][1];
        double length = sqrt(dx*dx + dy*dy);
        p.lengths[i] = length;
        if(length > 0){
            p.normals[i][0] = -dy/length;
            p.normals[i][1] = dx/length;
            p.tangents[i][0] = dx/length;
            p.tangents[i][1] = dy/length;
        }
        else{
            p.normals[i][0] = 0;
            p.normals[i][1] = 0;
        }
        x_com += p.vertices[i][0];
        y_com += p.vertices[i][1];
    }
    for (int i = 0; i < 3; i++) {
        p.vertices[i][0] -= x_com / 3;
        p.vertices[i][1] -= y_com / 3;
    }
    
    return p;
}

void rotateTranslatePolygon(intruder* intr, double phi, double X, double Y){

    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    polygon* p_r = &intr->real_shape;
    const polygon* p_b = &intr->body_frame_shape;
    for(int i = 0; i < p_b->n_vertices; i++){
        double x = p_b->vertices[i][0];
        double y = p_b->vertices[i][1];
        p_r->vertices[i][0] = x * cos_phi - y * sin_phi + X;
        p_r->vertices[i][1] = x * sin_phi + y * cos_phi + Y;
        
        
        double nx = p_b->normals[i][0];
        double ny = p_b->normals[i][1];
        p_r->normals[i][0] = nx * cos_phi - ny * sin_phi;
        p_r->normals[i][1] = nx * sin_phi + ny * cos_phi;
        
        double tx = p_b->tangents[i][0];
        double ty = p_b->tangents[i][1];
        p_r->tangents[i][0] = tx * cos_phi - ty * sin_phi;
        p_r->tangents[i][1] = tx * sin_phi + ty * cos_phi;
    }

}

void intruderInit(intruder* intr, double x, double y, double vx, double vy, double angle, double omega, double M, double Im, int shape, double* info){
    intr->x = x;
    intr->y = y;
    intr->vx = vx;
    intr->vy = vy;
    intr->angle = angle;
    intr->omega = omega;
    intr->M = M;
    intr->Im = Im;
    intr->coll = 0;
    intr->t = 0;
    switch(shape){
        case SQUARE:
            intr->body_frame_shape = square(info[0]);
            break;
        case TRIANGLE:
            intr->body_frame_shape = triangle(info[0], info[1]);
            break;
        default:
            printf("Unknown intruder shape\n");
            break;
    }
    intr->real_shape = intr->body_frame_shape;
    rotateTranslatePolygon(intr, angle, x, y);
}

int countFakeParticlesForIntruder(intruder* intr, double rad){
    double dr = rad;
    int count = 0;
    polygon* p = &intr->real_shape;
    for(int i = 0; i < p->n_vertices; i++){
        double dist = 0;
        while (dist < p->lengths[i]){
            count++;
            dist += dr;
        }
    }
    return count;
}

void printIntruderShapeFromSmallParticles(FILE* fichier, intruder* intr, double rad, int N){
    double dr = rad;
    int N_count = N;
    polygon* p = &intr->real_shape;
    for(int i = 0; i < p->n_vertices; i++){
        double dist = 0;
        while (dist < p->lengths[i]){
            double x = p->vertices[i][0] + dist * p->tangents[i][0] + rad*p->normals[i][0];
            double y = p->vertices[i][1] + dist * p->tangents[i][1] + rad*p->normals[i][1];
            fprintf(fichier, "%d %d %lf %lf %d %d %lf %d %d\n", N_count, 2, x, y, 0, 0, rad, 1, 1);
            N_count++;
            dist += dr;
        }
    }
}

double collisionTimeIntruder(intruder* intr, particle* p) {

    const double tol = 1e-10;
    const int max_bisect_iters = 80;
    const double t_horizon = 10.0;
    const double min_dt = 1e-4;
    const double max_dt = 1e-1;

    double g0 = signedDistanceParticlePolygon(p, intr, 0.0);
    if (g0 <= 0.0){
        return 0.0;
    }

    double relvx = p->vx - intr->vx;
    double relvy = p->vy - intr->vy;
    double v_bound = sqrt(relvx*relvx + relvy*relvy) + fabs(intr->omega) * (2.0*p->rad + 1.0);
    if (v_bound <= 1e-14){
        return never;
    }

    double dt = 0.2 * p->rad / v_bound;
    if (dt < min_dt){
        dt = min_dt;
    }
    if (dt > max_dt){
        dt = max_dt;
    }

    double t_prev = 0.0;

    for (double tc = dt; tc <= t_horizon; tc += dt) {
        double g_now = signedDistanceParticlePolygon(p, intr, tc);
        if (g_now <= 0.0) {
            double a = t_prev;
            double b = tc;
            for (int iter = 0; iter < max_bisect_iters; ++iter) {
                double m = 0.5 * (a + b);
                double g_m = signedDistanceParticlePolygon(p, intr, m);
                if (g_m > 0.0) {
                    a = m;
                } else {
                    b = m;
                }
                if ((b - a) <= tol) {
                    break;
                }
            }
            return 0.5 * (a + b);
        }
        t_prev = tc;
    }
    return never;
}


double signedDistanceParticlePolygon(particle* p, intruder* intru, double t) {
    double px = p->x + p->vx * t;
    double py = p->y + p->vy * t;

    double cos_theta = cos(intru->angle + intru->omega * t);
    double sin_theta = sin(intru->angle + intru->omega * t);
    double min_dist = never;
    polygon* poly = &(intru->body_frame_shape);
    for (int i = 0; i < poly->n_vertices; ++i) {
        int j = (i + 1) % poly->n_vertices;

        // Rotate vertices based on omega and t
        double ax = poly->vertices[i][0] * cos_theta - poly->vertices[i][1] * sin_theta + intru->x + intru->vx * t;
        double ay = poly->vertices[i][0] * sin_theta + poly->vertices[i][1] * cos_theta + intru->y + intru->vy * t;
        double bx = poly->vertices[j][0] * cos_theta - poly->vertices[j][1] * sin_theta + intru->x + intru->vx * t;
        double by = poly->vertices[j][0] * sin_theta + poly->vertices[j][1] * cos_theta + intru->y + intru->vy * t;

        double dist = pointSegmentDistance(px, py, ax, ay, bx, by);
        if (dist < min_dist) {
            min_dist = dist;
        }
    }

    return min_dist - p->rad;
}

void computeCollisionNormal(particle* p, intruder* intru, double t, double* nx, double* ny) {
    double px = p->x + p->vx * t;
    double py = p->y + p->vy * t;

    double cos_theta = cos(intru->angle + intru->omega * t);
    double sin_theta = sin(intru->angle + intru->omega * t);
    int closest_edge = -1;
    double min_dist = never;
    polygon* poly = &(intru->body_frame_shape);
    for (int i = 0; i < poly->n_vertices; ++i) {
        int j = (i + 1) % poly->n_vertices;

        // Rotate vertices based on omega and t
        double ax = poly->vertices[i][0] * cos_theta - poly->vertices[i][1] * sin_theta + intru->x + intru->vx * t;
        double ay = poly->vertices[i][0] * sin_theta + poly->vertices[i][1] * cos_theta + intru->y + intru->vy * t;
        double bx = poly->vertices[j][0] * cos_theta - poly->vertices[j][1] * sin_theta + intru->x + intru->vx * t;
        double by = poly->vertices[j][0] * sin_theta + poly->vertices[j][1] * cos_theta + intru->y + intru->vy * t;

        double dist = pointSegmentDistance(px, py, ax, ay, bx, by);
        if (dist < min_dist) {
            min_dist = dist;
            closest_edge = i;
        }
    }

    if (closest_edge >= 0) {
        // Rotate normals based on omega and t
        *nx = poly->normals[closest_edge][0] * cos_theta - poly->normals[closest_edge][1] * sin_theta;
        *ny = poly->normals[closest_edge][0] * sin_theta + poly->normals[closest_edge][1] * cos_theta;
    } else {
        *nx = 0.0;
        *ny = 0.0;
    }
}


double pointSegmentDistance(double px, double py, double ax, double ay, double bx, double by) {
    double abx = bx - ax;
    double aby = by - ay;
    double apx = px - ax;
    double apy = py - ay;
    PBC(&abx, &aby);
    PBC(&apx, &apy);

    double ab_len2 = abx * abx + aby * aby;
    if (ab_len2 <= 0.0){
        return sqrt(apx * apx + apy * apy);
    }
    double proj = (apx * abx + apy * aby) / ab_len2;

    if (proj < 0.0) {
        return sqrt(apx * apx + apy * apy);
    } else if (proj > 1.0) {
        double bpx = px - bx;
        double bpy = py - by;
        PBC(&bpx, &bpy);
        return sqrt(bpx * bpx + bpy * bpy);
    } else {
        double cx = ax + proj * abx;
        double cy = ay + proj * aby;
        double cpx = px - cx;
        double cpy = py - cy;
        PBC(&cpx, &cpy);
        return sqrt(cpx * cpx + cpy * cpy);
    }
}

void freeFlyIntruder(intruder* intr){
    double dt = t - intr->t;
    intr->t = t;
    intr->angle += intr->omega*dt;
    intr->x += intr->vx*dt;
    intr->y += intr->vy*dt;
    PBCpostX(&intr->x);
    PBCpostY(&intr->y);
    rotateTranslatePolygon(intr, intr->angle, intr->x, intr->y);
}


int particleOutsideIntruder(particle* p, intruder* intru) {
    double dist = signedDistanceParticlePolygon(p, intru, 0.0);
    const double tol = 1e-8; // tolerance for touching
    if (dist < -tol) {
        return 0; // inside
    } 
    return 1; // outside or touching
}
