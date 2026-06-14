#include "EDMD.h"
#include "intruder.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define fast_collision_time 1

double polygon_I_com(const polygon* p, double M){
    int N = p->n_vertices;

    double A2 = 0.0;     // 2A
    double Cx6A = 0.0;   // 6A*Cx
    double Cy6A = 0.0;   // 6A*Cy
    double J0  = 0.0;    // polar moment about origin, times density handled later

    for(int i=0;i<N;i++){
        int j = (i+1)%N;
        double xi = p->vertices[i][0], yi = p->vertices[i][1];
        double xj = p->vertices[j][0], yj = p->vertices[j][1];

        double c = xi*yj - xj*yi;      // cross product z-component
        A2  += c;
        Cx6A += (xi + xj) * c;
        Cy6A += (yi + yj) * c;

        double ri2 = xi*xi + yi*yi;
        double rj2 = xj*xj + yj*yj;
        double dot = xi*xj + yi*yj;

        J0 += c * (ri2 + dot + rj2);
    }

    double A = 0.5 * A2;
    double absA = (A >= 0.0) ? A : -A;
    if(absA == 0.0) return 0.0; // degenerate

    double Cx = Cx6A / (6.0 * A);
    double Cy = Cy6A / (6.0 * A);

    J0 *= (1.0/12.0);

    // areal density sigma = M/|A|
    double sigma = M / absA;

    // parallel-axis correction to COM
    double Icom = sigma * (J0 - A * (Cx*Cx + Cy*Cy));

    // Guard against tiny negative due to roundoff
    if(Icom < 0.0 && Icom > -1e-12) Icom = 0.0;
    return Icom;
}

void computeConstPol(polygon* p){
    double x_com = 0.0, y_com = 0.0;
    for (int i = 0; i < p->n_vertices; i++){
        x_com += p->vertices[i][0];
        y_com += p->vertices[i][1];
    }
    x_com /= p->n_vertices;
    y_com /= p->n_vertices;

    for (int i = 0; i < p->n_vertices; i++){
        p->vertices[i][0] -= x_com;
        p->vertices[i][1] -= y_com;
    }

    p->max_extent = 0.0;
    for (int i = 0; i < p->n_vertices; i++){
        int j = (i+1) % p->n_vertices;

        double dist = hypot(p->vertices[i][0], p->vertices[i][1]);
        if (dist > p->max_extent) p->max_extent = dist;

        double dx = p->vertices[j][0] - p->vertices[i][0];
        double dy = p->vertices[j][1] - p->vertices[i][1];
        double length = hypot(dx, dy);
        p->lengths[i] = length;

        if (length > 0.0){
            p->tangents[i][0] = dx/length;
            p->tangents[i][1] = dy/length;

            p->normals[i][0] = -dy/length;
            p->normals[i][1] = dx/length;
        } else {
            p->normals[i][0] = p->normals[i][1] = 0.0;
            p->tangents[i][0] = p->tangents[i][1] = 0.0;
        }
    }
}

polygon generalized_wheel(double radius_intern,
                          double radius_extern,
                          int n_spikes,
                          double chi_w)   // chirality knob: chi_w = 0 -> achiral
{
    polygon p;
    p.n_vertices = 2 * n_spikes;

    const double pitch = 2.0 * M_PI / n_spikes;

    /*  We define delta = pi/n + chi_w
        - chi_w = 0       -> inner vertices halfway between outer vertices (achiral)
        - chi_w != 0      -> chiral
        - mirror shape: chi_w -> -chi_w  (equivalently delta -> pitch - delta)
    */
    const double delta = M_PI / n_spikes + chi_w;

    for (int i = 0; i < n_spikes; i++) {
        double th = i * pitch;

        // outer tip
        p.vertices[2*i][0]     = radius_extern * cos(th);
        p.vertices[2*i][1]     = radius_extern * sin(th);

        // inner valley (angularly shifted by delta)
        p.vertices[2*i + 1][0] = radius_intern * cos(th + delta);
        p.vertices[2*i + 1][1] = radius_intern * sin(th + delta);
    }

    computeConstPol(&p);
    return p;
}

polygon poly(double radius, int n_vertices){
    polygon p;
    p.n_vertices = n_vertices;
    for(int i = 0; i < n_vertices; i++){
        double angle = (2.0 * M_PI * i) / (double)n_vertices;
        p.vertices[i][0] = radius * cos(angle);
        p.vertices[i][1] = radius * sin(angle);
    }
    computeConstPol(&p);
    return p;
}


polygon cross(double extension, double width){
    polygon p;
    p.n_vertices = 20;

    double w2 = width / 2.0;


    p.vertices[0][0] = w2; p.vertices[0][1] = w2;
    p.vertices[1][0] = w2 + extension; p.vertices[1][1] = w2;
    p.vertices[2][0] = w2 + extension; p.vertices[2][1] = w2 + extension;
    p.vertices[3][0] = w2 + extension + width; p.vertices[3][1] = w2 + extension;
    p.vertices[4][0] = w2 + extension + width; p.vertices[4][1] = -w2;
    p.vertices[5][0] = w2 ; p.vertices[5][1] = -w2;
    p.vertices[6][0] = w2 ; p.vertices[6][1] = -w2 - extension;
    p.vertices[7][0] = w2 + extension; p.vertices[7][1] = -w2 - extension;
    p.vertices[8][0] = w2 + extension; p.vertices[8][1] = -w2 - extension - width;
    p.vertices[9][0] = - w2; p.vertices[9][1] = -w2 - extension - width;
    p.vertices[10][0] = - w2; p.vertices[10][1] = -w2;
    p.vertices[11][0] = - w2 - extension; p.vertices[11][1] = -w2;
    p.vertices[12][0] = - w2 - extension; p.vertices[12][1] = -w2 - extension;
    p.vertices[13][0] = - w2 - extension - width; p.vertices[13][1] = -w2 - extension;
    p.vertices[14][0] = - w2 - extension - width; p.vertices[14][1] = w2;
    p.vertices[15][0] = - w2; p.vertices[15][1] = w2;
    p.vertices[16][0] = - w2; p.vertices[16][1] = w2 + extension;
    p.vertices[17][0] = - w2 - extension; p.vertices[17][1] = w2 + extension;
    p.vertices[18][0] = - w2 - extension; p.vertices[18][1] = w2 + extension + width;
    p.vertices[19][0] = w2; p.vertices[19][1] = w2 + extension + width;

    computeConstPol(&p);
    return p;
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

polygon stress_polygon_wheel(double radius_base, double spike_amp, double tri_amp, double phase){
    polygon p;
    const int n_vertices = 20;
    const double r_floor = 0.25 * radius_base;
    p.n_vertices = n_vertices;

    for (int i = 0; i < n_vertices; i++){
        double th = (2.0 * M_PI * i) / (double)n_vertices;

        // Hybrid radius profile: wheel-like spikes + triangular anisotropy + rough ripples.
        double wheel = spike_amp * pow(fabs(sin(6.0 * th + phase)), 1.6);
        double tri = tri_amp * cos(3.0 * th - 0.5 * phase);
        double ripple = 0.25 * spike_amp * sin(11.0 * th + 0.7 * phase);

        double r = radius_base + wheel + tri + ripple;
        if (r < r_floor) r = r_floor;

        p.vertices[i][0] = r * cos(th);
        p.vertices[i][1] = r * sin(th);
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

void intruderInit(intruder* intr, double x, double y, double vx, double vy, double angle, double omega, double M, double Im, int shape, double* info, int km_enabled, IntruderKM* km){
    intr->x = x;
    intr->y = y;
    intr->vx = vx;
    intr->vy = vy;
    intr->angle = angle;
    intr->omega = omega;
    intr->M = M;
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
        case GEN_WHEEL:
            intr->body_frame_shape = generalized_wheel(info[0], info[1], (int)info[2], info[3]);
            break;
        case COMPLEX_WHEEL:
            intr->body_frame_shape = stress_polygon_wheel(info[0], info[1], info[2], info[3]);
            break;
        case CROSS:
            intr->body_frame_shape = cross(info[0], info[1]);
            break;
        case POLYGON:
            intr->body_frame_shape = poly(info[0], (int)info[1]);
            intr->angle = M_PI/((int)info[1]);
            break;
        default:
            printf("Unknown intruder shape\n");
            break;
    }
    if (Im <= 0.0){
        intr->Im = polygon_I_com(&intr->body_frame_shape, M);
    }
    else{
        intr->Im = Im;
    }
    intr->real_shape = intr->body_frame_shape;
    rotateTranslatePolygon(intr, angle, x, y);

    km->enabled = km_enabled;  
	km->has_prev = 1;
	km->last_t = t;
	km->phi0 = intr->angle;
	km->vx0 = intr->vx;
	km->vy0 = intr->vy;
	km->om0 = intr->omega;
	km->nsamples = 0ULL;
	km->sum_dt = 0.0;
    km->sum_Etrans_dt = 0.0;
    km->sum_Erot_dt = 0.0;
	km->sum_dU[0] = km->sum_dU[1] = km->sum_dU[2] = 0.0;
	km->sum_intVb[0] = km->sum_intVb[1] = 0.0;
	km->sum_intOm_dt = 0.0;
	km->sum_intOmVb[0] = km->sum_intOmVb[1] = 0.0;
	for (int ii = 0; ii < 4; ++ii){
		for (int jj = 0; jj < 4; ++jj){
			km->XtWX[ii][jj] = 0.0;
			km->XTX[ii][jj] = 0.0;
		}
		for (int kk = 0; kk < 3; ++kk){
			km->XtWY[ii][kk] = 0.0;
			km->XTY[ii][kk] = 0.0;
		}
	}
	for (int kk = 0; kk < 3; ++kk){
		for (int ll = 0; ll < 3; ++ll){
			km->YTY[kk][ll] = 0.0;
		}
	}
}

int countFakeParticlesForIntruder(intruder* intr, double rad){
    if (rad == 0.0) return 0;
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
    if (rad == 0.0) return;
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

#if fast_collision_time
double collisionTimeIntruder(intruder* intr, particle* p) {

    const double tol = 1e-10;
    const int max_bisect_iters = 20;
    const double t_horizon = 10.0;
    const double min_dt = 1e-2;
    const double max_dt = 1e-1;

    const double extent = intr->body_frame_shape.max_extent;

    double relvx = p->vx - intr->vx;
    double relvy = p->vy - intr->vy;
    double L = extent + p->rad;
    double v_bound = hypot(relvx, relvy) + fabs(intr->omega) * L;
    if (v_bound <= 1e-14){
        return never;
    }

    // Cheap broad-phase reject in world frame before expensive signed distance.
    double dx0 = p->x - intr->x;
    double dy0 = p->y - intr->y;
    PBC(&dx0, &dy0);
    double r0 = hypot(dx0, dy0);
    if (r0 > (L + v_bound * t_horizon)){
        return never;
    }

    double g0 = signedDistanceParticlePolygon(p, intr, 0.0);
    if (g0 > extent + 4*p->rad){
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

    // If overlapping, walk forward to find exit
    if (g0 <= 0.0) {
        const double step = p->rad * 0.01 / v_bound;
        double tc_scan = 0.0;
        double g_scan = g0;
        for (int it = 0; it < 2000; ++it) {
            tc_scan += step;
            if (tc_scan > t_horizon) return never;
            g_scan = signedDistanceParticlePolygon(p, intr, tc_scan);
            if (g_scan > 1e-12) break;
        }
        if (g_scan <= 1e-12) return never;
        t_prev = tc_scan;
    }

    for (double tc = t_prev + dt; tc <= t_horizon; tc += dt) {
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
            return a;
        }
        t_prev = tc;
    }
    return never;
}
#else
double collisionTimeIntruder(intruder* intr, particle* p) {

    const double tol = 1e-10;
    const int max_bisect_iters = 80;
    const double t_horizon = 10.0;

    // "near contact" threshold: if we're closer than this, treat as immediate contact
    const double g_touch = 1e-12;

    // Cap only from above (do NOT force a minimum dt, that breaks conservativeness)
    const double max_dt = 1e-1;

    const double extent = intr->body_frame_shape.max_extent;

    // Conservative bound on relative speed of the particle vs the rotating intruder
    const double relvx = p->vx - intr->vx;
    const double relvy = p->vy - intr->vy;
    const double L = extent + p->rad;
    const double v_bound = hypot(relvx, relvy) + fabs(intr->omega) * L;

    if (v_bound <= 1e-14) {
        return never;
    }

    // Cheap broad-phase reject in world frame before expensive signed distance.
    double dx0 = p->x - intr->x;
    double dy0 = p->y - intr->y;
    PBC(&dx0, &dy0);
    double r0 = hypot(dx0, dy0);
    if (r0 > (L + v_bound * t_horizon)) {
        return never;
    }

    double g0 = signedDistanceParticlePolygon(p, intr, 0.0);

    // Early reject (keep your logic, but use the body-frame extent as the stable one)
    if (g0 > extent + 4.0 * p->rad) {
        return never;
    }

    double tc = 0.0;
    double g  = g0;

    // ---------------------------------------------------------------
    // FIX for non-convex polygons (e.g. chiral wheel):
    // If g0 <= 0 the particle numerically overlaps the polygon boundary
    // (can happen after bisection tolerance in a concave notch).
    // Instead of silently dropping the collision, walk forward in time
    // until the signed distance becomes clearly positive, then continue
    // with the normal conservative advancement.
    // ---------------------------------------------------------------
    if (g0 <= 0.0) {
        const double step = p->rad * 0.01 / v_bound;
        for (int it = 0; it < 2000; ++it) {
            tc += step;
            if (tc > t_horizon) return never;
            g = signedDistanceParticlePolygon(p, intr, tc);
            if (g > g_touch) break;
        }
        if (g <= g_touch) return never;   // could not exit overlap
        // fall through to conservative advancement below
    }

    for (int it = 0; it < 20000 && tc < t_horizon; ++it) {

        if (g <= g_touch) {
            return tc; // extremely close => schedule collision now
        }

        // Conservative advancement step: cannot "jump over" the boundary
        double dt = 0.5 * g / v_bound;
        if (dt > max_dt) dt = max_dt;

        double t_next = tc + dt;
        if (t_next > t_horizon) t_next = t_horizon;

        double g_next = signedDistanceParticlePolygon(p, intr, t_next);

        // Bracket found => bisection for accurate root
        if (g_next <= 0.0) {
            double a = tc, b = t_next;
            for (int k = 0; k < max_bisect_iters; ++k) {
                double m = 0.5 * (a + b);
                double g_m = signedDistanceParticlePolygon(p, intr, m);
                if (g_m > 0.0) a = m;
                else          b = m;
                if ((b - a) <= tol) break;
            }
            // Return the POSITIVE side of the bracket, not the midpoint.
            // This guarantees the particle stops just before the boundary
            // and is never placed in an overlapping state. Critical for
            // non-convex polygons where overlap with one edge can cause
            // the particle to be "inside" the polygon from another edge's
            // perspective, triggering g0<=0 → missed collisions.
            return a;
        }

        tc = t_next;
        g  = g_next;
    }

    return never;
}
#endif

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

    intr->angle -= 2.0*M_PI * floor(intr->angle / (2.0*M_PI));
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



inline void km_integrate_vbody(double vx, double vy, double phi0, double om, double dt,
                                      double* int_vbx, double* int_vby){
    // Computes ∫_0^dt V_body(t) dt, with V_body(t)=R(-(phi0+om t)) V_lab (vx,vy)
    double eps = 1e-14;
    if (fabs(om) < eps){
        const double c0 = cos(phi0);
        const double s0 = sin(phi0);
        const double vbx =  vx*c0 + vy*s0;
        const double vby = -vx*s0 + vy*c0;
        *int_vbx = dt * vbx;
        *int_vby = dt * vby;
        return;
    }
    const double phi1 = phi0 + om*dt;
    const double Icos = (sin(phi1) - sin(phi0)) / om;
    const double Isin = (cos(phi0) - cos(phi1)) / om;
    *int_vbx = vx*Icos + vy*Isin;
    *int_vby = -vx*Isin + vy*Icos;
}

int solve3x3(const double A_in[3][3], const double b_in[3], double x_out[3]){
    double A[3][3];
    double b[3];
    for (int i=0;i<3;++i){
        b[i]=b_in[i];
        for (int j=0;j<3;++j) A[i][j]=A_in[i][j];
    }
    for (int k=0;k<3;++k){
        int piv=k;
        double amax=fabs(A[k][k]);
        for (int i=k+1;i<3;++i){
            double v=fabs(A[i][k]);
            if (v>amax){amax=v;piv=i;}
        }
        if (amax<1e-18) return 0;
        if (piv!=k){
            for (int j=k;j<3;++j){ double tmp=A[k][j]; A[k][j]=A[piv][j]; A[piv][j]=tmp; }
            double tb=b[k]; b[k]=b[piv]; b[piv]=tb;
        }
        double inv = 1.0/A[k][k];
        for (int j=k;j<3;++j) A[k][j]*=inv;
        b[k]*=inv;
        for (int i=0;i<3;++i){
            if (i==k) continue;
            double f=A[i][k];
            if (f==0.0) continue;
            for (int j=k;j<3;++j) A[i][j]-=f*A[k][j];
            b[i]-=f*b[k];
        }
    }
    x_out[0]=b[0]; x_out[1]=b[1]; x_out[2]=b[2];
    return 1;
}

int solve4x4_3rhs(const double A_in[4][4], const double B_in[4][3], double X_out[4][3]){
    double A[4][4];
    double B[4][3];
    for (int i=0;i<4;++i){
        for (int j=0;j<4;++j) A[i][j]=A_in[i][j];
        for (int k=0;k<3;++k) B[i][k]=B_in[i][k];
    }
    for (int col=0; col<4; ++col){
        int piv=col;
        double amax=fabs(A[col][col]);
        for (int i=col+1;i<4;++i){
            double v=fabs(A[i][col]);
            if (v>amax){amax=v;piv=i;}
        }
        if (amax<1e-18) return 0;
        if (piv!=col){
            for (int j=col;j<4;++j){ double tmp=A[col][j]; A[col][j]=A[piv][j]; A[piv][j]=tmp; }
            for (int k=0;k<3;++k){ double tb=B[col][k]; B[col][k]=B[piv][k]; B[piv][k]=tb; }
        }
        double inv = 1.0/A[col][col];
        for (int j=col;j<4;++j) A[col][j]*=inv;
        for (int k=0;k<3;++k) B[col][k]*=inv;
        for (int i=0;i<4;++i){
            if (i==col) continue;
            double f=A[i][col];
            if (f==0.0) continue;
            for (int j=col;j<4;++j) A[i][j]-=f*A[col][j];
            for (int k=0;k<3;++k) B[i][k]-=f*B[col][k];
        }
    }
    for (int i=0;i<4;++i){
        for (int k=0;k<3;++k) X_out[i][k]=B[i][k];
    }
    return 1;
}

void km_finalize_and_dump(IntruderKM km, const char* intruderKMName, intruder intru){
    if (!km.enabled){
        return;
    }
    if (km.sum_dt <= 0.0 || km.nsamples < 10ULL){
        FILE* f = fopen(intruderKMName, "w");
        if (f){
            fprintf(f, "# intruderKM: insufficient samples (nsamples=%llu, sum_dt=%g)\n",
                    (unsigned long long)km.nsamples, km.sum_dt);
            fclose(f);
        }
        return;
    }

    // Solve weighted normal equations for B_fit (4x3): B = (XtWX)^{-1} XtWY
    double Bfit[4][3];
    if (!solve4x4_3rhs(km.XtWX, km.XtWY, Bfit)){
        fprintf(stderr, "intruderKM: solve4x4 failed (XtWX singular)\n");
        return;
    }

    // Extract L = (Bfit[1:4,:])^T  (3x3)
    double Lmat[3][3];
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            Lmat[i][j] = Bfit[1 + j][i];
        }
    }

    const double invT = 1.0 / km.sum_dt;
    double Ubar[3] = { km.sum_intVb[0]*invT, km.sum_intVb[1]*invT, km.sum_intOm_dt*invT };
    double abar[3] = { km.sum_dU[0]*invT, km.sum_dU[1]*invT, km.sum_dU[2]*invT };

    // F0_mean = <a_b> - L <U>
    double LU[3] = {0.0,0.0,0.0};
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            LU[i] += Lmat[i][j] * Ubar[j];
        }
    }
    double F0[3] = { abar[0] - LU[0], abar[1] - LU[1], abar[2] - LU[2] };

    // C_U = (-<Ω Vb_y>, <Ω Vb_x>, 0)
    double OmVbx = km.sum_intOmVb[0] * invT;
    double OmVby = km.sum_intOmVb[1] * invT;
    double CU[3] = { -OmVby, OmVbx, 0.0 };

    // U_ou_pred = -L^{-1} F0 ;  U_ou_equiv = <U> - L^{-1} C_U
    double tmp[3];
    if (!solve3x3(Lmat, F0, tmp)){
        fprintf(stderr, "intruderKM: solve3x3(L,F0) failed (L singular)\n");
        return;
    }
    double U_ou_pred[3] = { -tmp[0], -tmp[1], -tmp[2] };

    if (!solve3x3(Lmat, CU, tmp)){
        fprintf(stderr, "intruderKM: solve3x3(L,CU) failed (L singular)\n");
        return;
    }
    double U_ou_equiv[3] = { Ubar[0] - tmp[0], Ubar[1] - tmp[1], Ubar[2] - tmp[2] };

    // Diffusion from residual increments: D = (R^T R) / sum_dt, with R = Y - X Bfinal
    double Bfinal[4][3];
    for (int k=0;k<3;++k) Bfinal[0][k] = F0[k];
    for (int j=0;j<3;++j){
        for (int i=0;i<3;++i){
            Bfinal[1 + j][i] = Lmat[i][j];
        }
    }

    double term1[3][3] = {{0}}; // B^T XTY
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            double ss=0.0;
            for (int p=0;p<4;++p) ss += Bfinal[p][i] * km.XTY[p][j];
            term1[i][j]=ss;
        }
    }

    double M4x3[4][3]; // XTX * B
    for (int p=0;p<4;++p){
        for (int j=0;j<3;++j){
            double ss=0.0;
            for (int q=0;q<4;++q) ss += km.XTX[p][q] * Bfinal[q][j];
            M4x3[p][j]=ss;
        }
    }

    double term2[3][3] = {{0}}; // B^T XTX B
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            double ss=0.0;
            for (int p=0;p<4;++p) ss += Bfinal[p][i] * M4x3[p][j];
            term2[i][j]=ss;
        }
    }

    double Dmat[3][3];
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            double R = km.YTY[i][j] - term1[i][j] - term1[j][i] + term2[i][j];
            Dmat[i][j] = R * invT;
        }
    }
    for (int i=0;i<3;++i){
        for (int j=i+1;j<3;++j){
            double ss = 0.5*(Dmat[i][j] + Dmat[j][i]);
            Dmat[i][j]=ss; Dmat[j][i]=ss;
        }
    }

    FILE* f = fopen(intruderKMName, "w");
    if (!f){
        fprintf(stderr, "intruderKM: could not open output file %s\n", intruderKMName);
        return;
    }

    fprintf(f, "# intruderKM online estimates (collision-to-collision, lag=1)\n");
    fprintf(f, "# burn_time = %.10g\n", km.burn_time);
    fprintf(f, "# nsamples  = %llu\n", (unsigned long long)km.nsamples);
    fprintf(f, "# sum_dt    = %.17g\n", km.sum_dt);
    fprintf(f, "# <T_trans> = %.17g\n", km.sum_Etrans_dt * invT);
    fprintf(f, "# <T_rot>   = %.17g\n", 2*km.sum_Erot_dt * invT);

    fprintf(f, "\n# L (1/time)\n");
    for (int i=0;i<3;++i){
        fprintf(f, "%.17g %.17g %.17g\n", Lmat[i][0], Lmat[i][1], Lmat[i][2]);
    }

    fprintf(f, "\n# F0_mean (accel units): [ax0 ay0 aom0]\n");
    fprintf(f, "%.17g %.17g %.17g\n", F0[0], F0[1], F0[2]);

    fprintf(f, "\n# D (for U increments): <r r^T> = D dt\n");
    for (int i=0;i<3;++i){
        fprintf(f, "%.17g %.17g %.17g\n", Dmat[i][0], Dmat[i][1], Dmat[i][2]);
    }

    fprintf(f, "\n# Time-averaged U = <(Vb_x,Vb_y,Omega)>\n");
    fprintf(f, "%.17g %.17g %.17g\n", Ubar[0], Ubar[1], Ubar[2]);

    fprintf(f, "\n# C_U : (-<Omega Vb_y>, <Omega Vb_x>, 0)\n");
    fprintf(f, "%.17g %.17g %.17g\n", CU[0], CU[1], CU[2]);

    fprintf(f, "\n# U_ou_pred (OU mean ignoring rotation term)\n");
    fprintf(f, "%.17g %.17g %.17g\n", U_ou_pred[0], U_ou_pred[1], U_ou_pred[2]);

    fprintf(f, "\n# U_ou_equiv (mean-corrected for rotating frame)\n");
    fprintf(f, "%.17g %.17g %.17g\n", U_ou_equiv[0], U_ou_equiv[1], U_ou_equiv[2]);

    fprintf(f, "\n# Mean body force from <a_b>: [Fx Fy] = M*<a_b_trans>\n");
    fprintf(f, "%.17g %.17g\n", intru.M*abar[0], intru.M*abar[1]);
    fprintf(f, "# Mean body force from C: [Fx Fy] = [ -M<Omega Vb_y>, M<Omega Vb_x> ]\n");
    fprintf(f, "%.17g %.17g\n", -intru.M*OmVby, intru.M*OmVbx);

    fprintf(f, "\n# F0_reg (from regression intercept Bfit[0][:], should match F0_mean)\n");
    fprintf(f, "%.17g %.17g %.17g\n", Bfit[0][0], Bfit[0][1], Bfit[0][2]);

    fprintf(f, "\n# a_bar (mean body-frame accel from impulses = sum_dU/sum_dt)\n");
    fprintf(f, "%.17g %.17g %.17g\n", abar[0], abar[1], abar[2]);

    fclose(f);

}

void rethermalize_bath_particle_at_intruder_contact(particle* p, intruder* P,
                                                                  double contactX, double contactY,
                                                                  double rx, double ry,
                                                                  double nx, double ny,
                                                                  double tx, double ty,
                                                                  double T) {
    const double rethermalize_eps = 1e-6;
    double U1 = drand(0, 1);
    double U2 = drand(0, 1);
    double g1 = sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2);
    double g2 = sqrt(-2.0 * log(U1)) * sin(2.0 * M_PI * U2);                                                                
    const double sigma = sqrt(T / p->m);
    const double gn_out = fabs(3*sigma * g1);
    const double gt_out = sigma * g2;

    const double vcx = P->vx - P->omega * ry;
    const double vcy = P->vy + P->omega * rx;

    p->vx = vcx + gn_out * nx + gt_out * tx;
    p->vy = vcy + gn_out * ny + gt_out * ty;

    p->x = contactX + (p->rad + rethermalize_eps) * nx;
    p->y = contactY + (p->rad + rethermalize_eps) * ny;
    PBC(&(p->x), &(p->y));
}