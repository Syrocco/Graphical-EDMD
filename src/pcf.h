#ifndef PCF_H
#define PCF_H

#include "EDMD.h"

typedef struct {
    double* g_r;
    double* r;
    int num_bins;
    double max_r;
    double bin_width;
} pcf_data;

typedef struct {
    double* g_r;
    double* g6_r;
    double* r;
    int num_bins;
    double max_r;
    double bin_width;
    double* k_vector;
} bond_order_pcf_data;

typedef struct {
    int num_bins;
    double max_r;
    double bin_width;
    double *r;
    double *g6_corr;
    int *counts;
} g6corr_data;

pcf_data* calculate_pcf(particle* particles, int N, double dr, double max_r, double Lx, double Ly);
bond_order_pcf_data* calculate_bond_order_pcf(particle* particles, int N, double dr, double max_r, double k_vector[2], double Lx, double Ly);
void free_pcf_data(pcf_data* data);
void free_bond_order_pcf_data(bond_order_pcf_data* data);
void save_pcf(const char* filename, particle* particles, int N, double dr, double max_r, double Lx, double Ly);
void save_pcf_boop(const char* pcf_filename, const char* boop_filename, particle* particles, int N, double dr, double max_r, double Lx, double Ly, double phi);
void find_max_structure_factor_bragg(particle* particles, int N, double Lx, double Ly, double expected_bragg, double k[2]);
void free_g6corr_data(g6corr_data* data);
g6corr_data* compute_g6_correlation(particle* particles, int N, double dr, double max_r, double Lx, double Ly);
void free_g6corr_data(g6corr_data* data);
void free_bond_order_pcf_data(bond_order_pcf_data* data);
void save_pcf_g6(const char* filename, particle* particles, int N, double dr, double max_r, double Lx, double Ly);
#endif /* PCF_H */
