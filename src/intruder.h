#ifndef INTRUDER_H
#define INTRUDER_H
#include "EDMD.h"

extern double t;
typedef struct polygon polygon;
struct polygon{
    double vertices[20][2];
    double lengths[20];
    double normals[20][2];
    double tangents[20][2];
    int n_vertices;
};

typedef struct intruder intruder;
struct intruder{
    double t;
    double x, y, vx, vy, angle, omega;
    double M, Im;
    polygon body_frame_shape;
    polygon real_shape;
    long unsigned int coll;
};

enum intruderShape{
    SQUARE,
    TRIANGLE
};

typedef struct collisionInfo collisionInfo;
struct collisionInfo{
    double time;
    double normal[2];
};

void intruderInit(intruder* intr, double x, double y, double vx, double vy, double angle, double omega, double M, double Im, int shape, double* info);
polygon square(double side_length);
polygon triangle(double base_length, double height);
void rotateTranslatePolygon(intruder* intr, double angle, double X, double Y);
void printIntruderShapeFromSmallParticles(FILE* fichier, intruder* intr, double rad, int N);
int countFakeParticlesForIntruder(intruder* intr, double rad);
double collisionTimeIntruder(intruder* intr, particle* p);

// Helper functions for collisionTimeIntruder
double signedDistanceParticlePolygon(particle* p, intruder* intru, double t);
void computeCollisionNormal(particle* p, intruder* intru, double t, double* nx, double* ny);
double pointSegmentDistance(double px, double py, double ax, double ay, double bx, double by);


void freeFlyIntruder(intruder* intru);
int particleOutsideIntruder(particle* p, intruder* intru);
#endif