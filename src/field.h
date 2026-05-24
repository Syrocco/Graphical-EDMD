#ifndef FIELD_H
#define FIELD_H

typedef struct particle particle;

int pressure_field_init(double Lx, double Ly, double dr_pressure, const char* output_path);
void pressure_field_close_file(void);
void pressure_field_free(void);
void pressure_field_step(double t, double delta_time);

void pressure_field_add_kinetic(const particle* p);
void pressure_field_add_collision(
	const particle* pi,
	double dx,
	double dy,
	double coll_term_x,
	double coll_term_y,
	double coll_term_xy,
	double coll_term_yx
);

#endif
