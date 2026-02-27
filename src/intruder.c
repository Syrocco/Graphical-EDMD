#include "EDMD.h"
#include "intruder.h"
#include <math.h>
#include <stdio.h>


void computeConstPol(polygon* p){
    p->max_extent = 0;
    double x_com = 0;
    double y_com = 0;
    for (int i = 0; i < p->n_vertices; i++){
        int j = (i+1) % p->n_vertices;
        x_com += p->vertices[i][0];
        y_com += p->vertices[i][1];
        double dist = sqrt(p->vertices[i][0]*p->vertices[i][0] + p->vertices[i][1]*p->vertices[i][1]);
        if (dist > p->max_extent){
            p->max_extent = dist;
        }
        double dx = p->vertices[j][0] - p->vertices[i][0];
        double dy = p->vertices[j][1] - p->vertices[i][1];
        double length = sqrt(dx*dx + dy*dy);
        p->lengths[i] = length;
        if(length > 0){
            p->normals[i][0] = -dy/length;
            p->normals[i][1] = dx/length;
            p->tangents[i][0] = dx/length;
            p->tangents[i][1] = dy/length;
        }
        else{
            p->normals[i][0] = 0;
            p->normals[i][1] = 0;
            p->tangents[i][0] = 0;
            p->tangents[i][1] = 0;
        }
    }
    x_com /= p->n_vertices;
    y_com /= p->n_vertices;
    for (int i = 0; i < p->n_vertices; i++){
        p->vertices[i][0] -= x_com;
        p->vertices[i][1] -= y_com;
    }
}

polygon wheel(double radius_intern, double radius_extern, int n_pikes){
    polygon p;
    p.n_vertices = 2*n_pikes;
    for(int i = 0; i < p.n_vertices; i++){
        double angle = i * (M_PI / n_pikes);
        double rad = (i % 2 == 0) ? radius_intern : radius_extern;
        p.vertices[i][0] = rad * cos(angle);
        p.vertices[i][1] = rad * sin(angle);
    }
    computeConstPol(&p);
    return p;
}

polygon chiral_weel(double radius_intern, double radius_extern, int n_pikes){
    polygon p;
    p.n_vertices = 2*n_pikes;
    for(int i = 0; i < n_pikes; i++){
        double angle = i *2*M_PI / n_pikes;
        p.vertices[2*i][0] = radius_extern * cos(angle);
        p.vertices[2*i][1] = radius_extern * sin(angle);
        p.vertices[2*i + 1][0] = radius_intern * cos(angle + M_PI/(2*n_pikes));
        p.vertices[2*i + 1][1] = radius_intern * sin(angle + M_PI/(2*n_pikes));
    }
    computeConstPol(&p);
    return p;
}



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
    
    computeConstPol(&p);
    
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
    
    computeConstPol(&p);
    
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
        case WHEEL:
            intr->body_frame_shape = wheel(info[0], info[1], (int)info[2]);
            break;
        case CHIRAL_WHEEL:
            intr->body_frame_shape = chiral_weel(info[0], info[1], (int)info[2]);
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
        return never;
    }
    else if (g0 > intr->real_shape.max_extent + 4*p->rad){
        return never;
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


// -----------------------------------------------------------------------------
// Signed distance helpers
//
// A "true" signed distance for a polygon is negative inside and positive
// outside, with magnitude equal to the distance to the boundary.
//
// Here we compute distances in the intruder *body frame* (non-periodic): we
// first bring the particle center in the same periodic image as the intruder
// center using PBC(), then rotate by -theta.
// -----------------------------------------------------------------------------

static inline double pointSegmentDistanceNoPBC(double px, double py,
                                               double ax, double ay,
                                               double bx, double by) {
    const double abx = bx - ax;
    const double aby = by - ay;
    const double apx = px - ax;
    const double apy = py - ay;

    const double ab_len2 = abx * abx + aby * aby;
    if (ab_len2 <= 0.0) {
        return sqrt(apx * apx + apy * apy);
    }

    const double proj = (apx * abx + apy * aby) / ab_len2;
    if (proj <= 0.0) {
        return sqrt(apx * apx + apy * apy);
    }
    if (proj >= 1.0) {
        const double bpx = px - bx;
        const double bpy = py - by;
        return sqrt(bpx * bpx + bpy * bpy);
    }

    const double cx = ax + proj * abx;
    const double cy = ay + proj * aby;
    const double dx = px - cx;
    const double dy = py - cy;
    return sqrt(dx * dx + dy * dy);
}

static inline int pointInPolygonEvenOdd(double px, double py, const polygon* poly) {
    // Even-odd (ray casting) rule. Works for any simple polygon.
    int inside = 0;
    for (int i = 0, j = poly->n_vertices - 1; i < poly->n_vertices; j = i++) {
        const double xi = poly->vertices[i][0];
        const double yi = poly->vertices[i][1];
        const double xj = poly->vertices[j][0];
        const double yj = poly->vertices[j][1];

        const int intersects = ((yi > py) != (yj > py)) &&
                               (px < (xj - xi) * (py - yi) / (yj - yi) + xi);
        if (intersects) {
            inside = !inside;
        }
    }
    return inside;
}


double signedDistanceParticlePolygon(particle* p, intruder* intru, double t) {
    // Particle center at time t (world frame)
    const double px = p->x + p->vx * t;
    const double py = p->y + p->vy * t;

    // Intruder center at time t (world frame)
    const double X = intru->x + intru->vx * t;
    const double Y = intru->y + intru->vy * t;

    // Relative vector (minimal image), then rotate into body frame
    double dx = px - X;
    double dy = py - Y;
    PBC(&dx, &dy);

    const double theta = intru->angle + intru->omega * t;
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);

    // World -> body (rotation by -theta)
    const double qx =  dx * cos_theta + dy * sin_theta;
    const double qy = -dx * sin_theta + dy * cos_theta;

    const polygon* poly = &(intru->body_frame_shape);
    double min_dist = never;
    for (int i = 0; i < poly->n_vertices; ++i) {
        const int j = (i + 1) % poly->n_vertices;
        const double ax = poly->vertices[i][0];
        const double ay = poly->vertices[i][1];
        const double bx = poly->vertices[j][0];
        const double by = poly->vertices[j][1];
        const double dist = pointSegmentDistanceNoPBC(qx, qy, ax, ay, bx, by);
        if (dist < min_dist) {
            min_dist = dist;
        }
    }


    // True signed distance for the point (inside => negative), then subtract r.
    const int inside = pointInPolygonEvenOdd(qx, qy, poly);
    const double sd_point = inside ? -min_dist : min_dist;
    return sd_point - p->rad;

}

void computeCollisionNormal(particle* p, intruder* intru, double t, double* nx, double* ny) {
    // Particle center at time t (world frame)
    const double px = p->x + p->vx * t;
    const double py = p->y + p->vy * t;

    // Intruder center at time t (world frame)
    const double X = intru->x + intru->vx * t;
    const double Y = intru->y + intru->vy * t;

    // Relative vector (minimal image), then rotate into body frame
    double dx = px - X;
    double dy = py - Y;
    PBC(&dx, &dy);

    const double theta = intru->angle + intru->omega * t;
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);

    // World -> body (rotation by -theta)
    const double qx =  dx * cos_theta + dy * sin_theta;
    const double qy = -dx * sin_theta + dy * cos_theta;

    const polygon* poly = &(intru->body_frame_shape);
    int closest_edge = -1;
    double min_dist = never;
    for (int i = 0; i < poly->n_vertices; ++i) {
        const int j = (i + 1) % poly->n_vertices;
        const double ax = poly->vertices[i][0];
        const double ay = poly->vertices[i][1];
        const double bx = poly->vertices[j][0];
        const double by = poly->vertices[j][1];
        const double dist = pointSegmentDistanceNoPBC(qx, qy, ax, ay, bx, by);
        if (dist < min_dist) {
            min_dist = dist;
            closest_edge = i;
        }
    }

    if (closest_edge >= 0) {
        // Body -> world rotation (by +theta)
        const double nbx = poly->normals[closest_edge][0];
        const double nby = poly->normals[closest_edge][1];
        *nx = nbx * cos_theta - nby * sin_theta;
        *ny = nbx * sin_theta + nby * cos_theta;
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
