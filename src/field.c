#include "field.h"

#include "EDMD.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	int bins_x;
	int bins_y;
	int discard_first;
	double dr;
	double lx;
	double ly;
	double* coll_x;
	double* coll_y;
	double* coll_xy;
	double* coll_yx;
	double* kinetic_x;
	double* kinetic_y;
	double* kinetic_xy;
	double* kinetic_yx;
	unsigned long int* coll_count;
	unsigned long int* kinetic_count;
	FILE* file;
} PressureFieldState;

static PressureFieldState g_pressure_field = {0};

static int pressure_field_bin_x(double x) {
	int bx = (int)(x / g_pressure_field.dr);
	if (bx < 0) {
		bx = 0;
	}
	if (bx >= g_pressure_field.bins_x) {
		bx = g_pressure_field.bins_x - 1;
	}
	return bx;
}

static int pressure_field_bin_y(double y) {
	int by = (int)(y / g_pressure_field.dr);
	if (by < 0) {
		by = 0;
	}
	if (by >= g_pressure_field.bins_y) {
		by = g_pressure_field.bins_y - 1;
	}
	return by;
}

static int pressure_field_bin_index(double x, double y) {
	if (g_pressure_field.bins_x <= 0 || g_pressure_field.bins_y <= 0) {
		return -1;
	}
	return pressure_field_bin_x(x) + g_pressure_field.bins_x * pressure_field_bin_y(y);
}

int pressure_field_init(double Lx, double Ly, double dr_pressure, const char* output_path) {
	g_pressure_field.lx = Lx;
	g_pressure_field.ly = Ly;
	g_pressure_field.dr = dr_pressure;
	g_pressure_field.discard_first = 1;
	g_pressure_field.bins_x = (int)ceil(Lx / dr_pressure);
	g_pressure_field.bins_y = (int)ceil(Ly / dr_pressure);

	int n_bins = g_pressure_field.bins_x * g_pressure_field.bins_y;
	g_pressure_field.coll_x = calloc((unsigned int)n_bins, sizeof(double));
	g_pressure_field.coll_y = calloc((unsigned int)n_bins, sizeof(double));
	g_pressure_field.coll_xy = calloc((unsigned int)n_bins, sizeof(double));
	g_pressure_field.coll_yx = calloc((unsigned int)n_bins, sizeof(double));
	g_pressure_field.kinetic_x = calloc((unsigned int)n_bins, sizeof(double));
	g_pressure_field.kinetic_y = calloc((unsigned int)n_bins, sizeof(double));
	g_pressure_field.kinetic_xy = calloc((unsigned int)n_bins, sizeof(double));
	g_pressure_field.kinetic_yx = calloc((unsigned int)n_bins, sizeof(double));
	g_pressure_field.coll_count = calloc((unsigned int)n_bins, sizeof(unsigned long int));
	g_pressure_field.kinetic_count = calloc((unsigned int)n_bins, sizeof(unsigned long int));

	if (
		g_pressure_field.coll_x == NULL || g_pressure_field.coll_y == NULL ||
		g_pressure_field.coll_xy == NULL || g_pressure_field.coll_yx == NULL ||
		g_pressure_field.kinetic_x == NULL || g_pressure_field.kinetic_y == NULL ||
		g_pressure_field.kinetic_xy == NULL || g_pressure_field.kinetic_yx == NULL ||
		g_pressure_field.coll_count == NULL || g_pressure_field.kinetic_count == NULL
	) {
		pressure_field_free();
		return 0;
	}

	g_pressure_field.file = fopen(output_path, "w");
	if (g_pressure_field.file == NULL) {
		pressure_field_free();
		return 0;
	}

	fprintf(g_pressure_field.file, "# t x y p pxx pyy pxy pyx p_coll p_kin coll_count kin_count\n");
	fprintf(g_pressure_field.file, "# dr_pressure %.8lf nx %d ny %d\n", dr_pressure, g_pressure_field.bins_x, g_pressure_field.bins_y);
	return 1;
}

void pressure_field_close_file(void) {
	if (g_pressure_field.file != NULL) {
		fclose(g_pressure_field.file);
		g_pressure_field.file = NULL;
	}
}

void pressure_field_free(void) {
	free(g_pressure_field.coll_x);
	free(g_pressure_field.coll_y);
	free(g_pressure_field.coll_xy);
	free(g_pressure_field.coll_yx);
	free(g_pressure_field.kinetic_x);
	free(g_pressure_field.kinetic_y);
	free(g_pressure_field.kinetic_xy);
	free(g_pressure_field.kinetic_yx);
	free(g_pressure_field.coll_count);
	free(g_pressure_field.kinetic_count);

	g_pressure_field.coll_x = NULL;
	g_pressure_field.coll_y = NULL;
	g_pressure_field.coll_xy = NULL;
	g_pressure_field.coll_yx = NULL;
	g_pressure_field.kinetic_x = NULL;
	g_pressure_field.kinetic_y = NULL;
	g_pressure_field.kinetic_xy = NULL;
	g_pressure_field.kinetic_yx = NULL;
	g_pressure_field.coll_count = NULL;
	g_pressure_field.kinetic_count = NULL;
	g_pressure_field.bins_x = 0;
	g_pressure_field.bins_y = 0;
	g_pressure_field.dr = 0;
	g_pressure_field.lx = 0;
	g_pressure_field.ly = 0;
	g_pressure_field.discard_first = 0;
}

void pressure_field_step(double t, double delta_time) {
	if (
		t != 0 &&
		g_pressure_field.file != NULL &&
		g_pressure_field.bins_x > 0 &&
		g_pressure_field.bins_y > 0 &&
		delta_time > 0 &&
		g_pressure_field.discard_first == 0
	) {
		for (int by = 0; by < g_pressure_field.bins_y; by++) {
			for (int bx = 0; bx < g_pressure_field.bins_x; bx++) {
				int idx = bx + g_pressure_field.bins_x * by;
				double x0 = bx * g_pressure_field.dr;
				double x1 = fmin((bx + 1) * g_pressure_field.dr, g_pressure_field.lx);
				double y0 = by * g_pressure_field.dr;
				double y1 = fmin((by + 1) * g_pressure_field.dr, g_pressure_field.ly);
				double bin_area = (x1 - x0) * (y1 - y0);
				double coll_px = -g_pressure_field.coll_x[idx] / (bin_area * delta_time);
				double coll_py = -g_pressure_field.coll_y[idx] / (bin_area * delta_time);
				double coll_pxy = -g_pressure_field.coll_xy[idx] / (bin_area * delta_time);
				double coll_pyx = -g_pressure_field.coll_yx[idx] / (bin_area * delta_time);
				double kin_px = g_pressure_field.kinetic_x[idx] / (bin_area * delta_time);
				double kin_py = g_pressure_field.kinetic_y[idx] / (bin_area * delta_time);
				double kin_pxy = g_pressure_field.kinetic_xy[idx] / (bin_area * delta_time);
				double kin_pyx = g_pressure_field.kinetic_yx[idx] / (bin_area * delta_time);
				double pxx = coll_px + kin_px;
				double pyy = coll_py + kin_py;
				double pxy = coll_pxy + kin_pxy;
				double pyx = coll_pyx + kin_pyx;
				double p = 0.5 * (pxx + pyy);
				double p_coll = 0.5 * (coll_px + coll_py);
				double p_kin = 0.5 * (kin_px + kin_py);
				double xc = x0 + 0.5 * (x1 - x0);
				double yc = y0 + 0.5 * (y1 - y0);
				fprintf(
					g_pressure_field.file,
					"%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %lu %lu\n",
					t,
					xc,
					yc,
					p,
					pxx,
					pyy,
					pxy,
					pyx,
					p_coll,
					p_kin,
					g_pressure_field.coll_count[idx],
					g_pressure_field.kinetic_count[idx]
				);
			}
		}
		fprintf(g_pressure_field.file, "\n");
		fflush(g_pressure_field.file);
	}

	if (t != 0 && g_pressure_field.discard_first == 1) {
		g_pressure_field.discard_first = 0;
	}

	if (g_pressure_field.bins_x > 0 && g_pressure_field.bins_y > 0) {
		int n_bins = g_pressure_field.bins_x * g_pressure_field.bins_y;
		memset(g_pressure_field.coll_x, 0, sizeof(double) * n_bins);
		memset(g_pressure_field.coll_y, 0, sizeof(double) * n_bins);
		memset(g_pressure_field.coll_xy, 0, sizeof(double) * n_bins);
		memset(g_pressure_field.coll_yx, 0, sizeof(double) * n_bins);
		memset(g_pressure_field.kinetic_x, 0, sizeof(double) * n_bins);
		memset(g_pressure_field.kinetic_y, 0, sizeof(double) * n_bins);
		memset(g_pressure_field.kinetic_xy, 0, sizeof(double) * n_bins);
		memset(g_pressure_field.kinetic_yx, 0, sizeof(double) * n_bins);
		memset(g_pressure_field.coll_count, 0, sizeof(unsigned long int) * n_bins);
		memset(g_pressure_field.kinetic_count, 0, sizeof(unsigned long int) * n_bins);
	}
}

void pressure_field_add_kinetic(const particle* p) {
	int idx = pressure_field_bin_index(p->x, p->y);
	if (idx < 0) {
		return;
	}
	g_pressure_field.kinetic_x[idx] += p->m * p->vx * p->vx;
	g_pressure_field.kinetic_y[idx] += p->m * p->vy * p->vy;
	g_pressure_field.kinetic_xy[idx] += p->m * p->vx * p->vy;
	g_pressure_field.kinetic_yx[idx] += p->m * p->vy * p->vx;
	g_pressure_field.kinetic_count[idx] += 1;
}

void pressure_field_add_collision(
	const particle* pi,
	double dx,
	double dy,
	double coll_term_x,
	double coll_term_y,
	double coll_term_xy,
	double coll_term_yx
) {
	double cx = pi->x + 0.5 * dx;
	double cy = pi->y + 0.5 * dy;
	PBCpostX(&cx);
	PBCpostY(&cy);
	int idx = pressure_field_bin_index(cx, cy);
	if (idx < 0) {
		return;
	}
	g_pressure_field.coll_x[idx] += coll_term_x;
	g_pressure_field.coll_y[idx] += coll_term_y;
	g_pressure_field.coll_xy[idx] += coll_term_xy;
	g_pressure_field.coll_yx[idx] += coll_term_yx;
	g_pressure_field.coll_count[idx] += 1;
}
