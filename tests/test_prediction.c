#include "EDMD.h"
#include "intruder.h"

#include <assert.h>
#include <math.h>
#include <string.h>

extern double Lx, Ly, halfLx, halfLy;
extern double t;
#if THREE_D
extern double Lz, halfLz;
#endif

static void assertClose(double actual, double expected){
	const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
	assert(fabs(actual - expected) <= 1e-12*scale);
}

static particle sampleParticle(void){
	particle p;
	memset(&p, 0, sizeof(p));
	p.x = 3.0;
	p.y = 7.0;
	p.vx = 2.0;
	p.vy = -4.0;
	p.rad = 1.25;
	p.t = 1.0;
	p.vr = 0.125;
	p.m = 2.0;
	p.coll = 17;
	#if THREE_D
	p.z = 11.0;
	p.vz = 0.5;
	#endif
	return p;
}

static void testAllFlightPhysicsArePure(void){
	const particle p = sampleParticle();
	const particle before = p;
	const double targetTime = 3.5;
	const double dt = targetTime - p.t;

	particlePrediction prediction = predictParticle(
		&p, targetTime, BALLISTIC_FLIGHT, 0);
	assertClose(prediction.x, p.x + dt*p.vx);
	assertClose(prediction.y, p.y + dt*p.vy);
	assertClose(prediction.vx, p.vx);
	assertClose(prediction.vy, p.vy);
	assertClose(prediction.rad, p.rad);

	const double acceleration = -0.3;
	prediction = predictParticle(
		&p, targetTime, CONSTANT_FIELD_FLIGHT, acceleration);
	assertClose(prediction.x, p.x + dt*p.vx);
	assertClose(prediction.y, p.y + dt*p.vy + 0.5*acceleration*dt*dt);
	assertClose(prediction.vx, p.vx);
	assertClose(prediction.vy, p.vy + acceleration*dt);
	assertClose(prediction.rad, p.rad);

	const double gamma = 0.2;
	const double decay = exp(-gamma*dt);
	prediction = predictParticle(
		&p, targetTime, LINEAR_DAMPING_FLIGHT, gamma);
	assertClose(prediction.x, p.x + p.vx/gamma*(1 - decay));
	assertClose(prediction.y, p.y + p.vy/gamma*(1 - decay));
	assertClose(prediction.vx, p.vx*decay);
	assertClose(prediction.vy, p.vy*decay);
	assertClose(prediction.rad, p.rad);

	prediction = predictParticle(&p, targetTime, GROWTH_FLIGHT, 0);
	assertClose(prediction.x, p.x + dt*p.vx);
	assertClose(prediction.y, p.y + dt*p.vy);
	assertClose(prediction.vx, p.vx);
	assertClose(prediction.vy, p.vy);
	assertClose(prediction.rad, p.rad + dt*p.vr);

	#if THREE_D
	assertClose(prediction.z, p.z + dt*p.vz);
	assertClose(prediction.vz, p.vz);
	#endif
	assert(memcmp(&p, &before, sizeof(p)) == 0);
}

static void testPredictionMatchesCommittedFlight(void){
	Lx = Ly = 1000.0;
	halfLx = halfLy = 500.0;
	#if THREE_D
	Lz = 1000.0;
	halfLz = 500.0;
	#endif

	particle p = sampleParticle();
	particle before = p;
	t = 2.25;
	const particlePrediction normal = predictFreeFlyNormal(&p, t);
	assert(memcmp(&p, &before, sizeof(p)) == 0);
	freeFlyNormal(&p);
	assertClose(p.t, t);
	assertClose(p.x, normal.x);
	assertClose(p.y, normal.y);
	assertClose(p.vx, normal.vx);
	assertClose(p.vy, normal.vy);

	p = sampleParticle();
	before = p;
	const particlePrediction growth = predictFreeFlyGrow(&p, t);
	assert(memcmp(&p, &before, sizeof(p)) == 0);
	freeFlyGrow(&p);
	assertClose(p.t, t);
	assertClose(p.x, growth.x);
	assertClose(p.y, growth.y);
	assertClose(p.rad, growth.rad);
}

static void testCollisionQueriesDoNotAdvanceParticles(void){
	Lx = Ly = 1000.0;
	halfLx = halfLy = 500.0;
	t = 2.0;

	particle p1 = sampleParticle();
	particle p2 = sampleParticle();
	p1.x = 100.0;
	p1.y = 100.0;
	p1.vx = 0.0;
	p1.vy = 0.0;
	p1.rad = 1.0;
	p1.t = t;
	p2.x = 110.0;
	p2.y = 100.0;
	p2.vx = -1.0;
	p2.vy = 0.0;
	p2.rad = 1.0;
	p2.t = 0.0;
	const particle before1 = p1;
	const particle before2 = p2;

	(void)collisionTimeNormal(&p1, &p2);
	(void)sphereTime(&p1, &p2);
	(void)separateTime(&p1, &p2);
	(void)collisionTimeGrow(&p1, &p2);
	assert(memcmp(&p1, &before1, sizeof(p1)) == 0);
	assert(memcmp(&p2, &before2, sizeof(p2)) == 0);

	p2 = p1;
	const particle overlappingBefore = p2;
	assert(collisionTimeNormal(&p1, &p2) == -never);
	assert(memcmp(&p2, &overlappingBefore, sizeof(p2)) == 0);
}

static void testIntruderQueryUsesOffsetsWithoutMutation(void){
	Lx = Ly = 1000.0;
	halfLx = halfLy = 500.0;

	intruder body;
	memset(&body, 0, sizeof(body));
	body.x = 100.0;
	body.y = 100.0;
	body.vx = 0.25;
	body.angle = 0.1;
	body.omega = 0.05;
	body.t = 1.0;
	body.body_frame_shape = square(4.0);
	body.real_shape = body.body_frame_shape;

	particle p = sampleParticle();
	p.x = 110.0;
	p.y = 100.0;
	p.vx = -1.0;
	p.vy = 0.0;
	p.rad = 1.0;
	p.t = 0.5;
	const intruder bodyBefore = body;
	const particle particleBefore = p;

	const double predictedTime = collisionTimeIntruderAt(&body, &p, 2.0);
	assert(memcmp(&body, &bodyBefore, sizeof(body)) == 0);
	assert(memcmp(&p, &particleBefore, sizeof(p)) == 0);

	intruder committedBody = body;
	particle committedParticle = p;
	t = 2.0;
	freeFlyIntruder(&committedBody);
	freeFlyNormal(&committedParticle);
	const double committedTime = collisionTimeIntruder(
		&committedBody, &committedParticle);
	assertClose(predictedTime, committedTime);
}

int main(void){
	testAllFlightPhysicsArePure();
	testPredictionMatchesCommittedFlight();
	testCollisionQueriesDoNotAdvanceParticles();
	testIntruderQueryUsesOffsetsWithoutMutation();
	return 0;
}
