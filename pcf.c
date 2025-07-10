#include "pcf.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "EDMD.h"
#include "boop.h"
#include <time.h>
#include <omp.h>
#include <complex.h>




pcf_data* calculate_pcf(particle* particles, int N, double dr, double max_r, double Lx, double Ly) {
    // Allocate memory for PCF data
    pcf_data* data = (pcf_data*)malloc(sizeof(pcf_data));
    if (data == NULL) return NULL;

    data->num_bins = (int)(max_r / dr);
    data->max_r = dr * data->num_bins;
    data->bin_width = dr;

    // Allocate memory for arrays
    data->g_r = (double*)calloc(data->num_bins, sizeof(double));
    data->r = (double*)malloc(data->num_bins * sizeof(double));
    if (data->g_r == NULL || data->r == NULL) {
        free_pcf_data(data);
        return NULL;
    }
    
    // Initialize r values at bin centers
    for (int i = 0; i < data->num_bins; i++) {
        data->r[i] = (i + 0.5) * data->bin_width;
    }
    
    // Calculate histogram of pair distances
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            
            // Apply periodic boundary conditions explicitly for large distances
            dx -= Lx * round(dx / Lx);
            dy -= Ly * round(dy / Ly);
            // Calculate distance
            double r = sqrt(dx*dx + dy*dy);
            
            if (r < max_r) {
                int bin = (int)(r / data->bin_width);
                if (bin < data->num_bins) {
                    data->g_r[bin] += 2.0; // Count each pair twice (i,j) and (j,i)
                }
            }
        }
    }
    
    // Normalize g(r)
    double volume = Lx * Ly;
    double density = N / volume;

    for (int i = 0; i < data->num_bins; i++) {
        // Calculate volume of shell: 2πr·dr in 2D
        double r = data->r[i];
        double shell_volume = 2 * M_PI * r * data->bin_width;
        
        // Normalize
        double norm = shell_volume * density * N;
        if (norm > 0) {
            data->g_r[i] /= norm;
        } else {
            data->g_r[i] = 0.0;
        }
    }
    
    return data;
}

bond_order_pcf_data* calculate_bond_order_pcf(particle* particles, int N, double dr, 
                                             double max_r, double k_vector[2], double Lx, double Ly) {
    bond_order_pcf_data* data = (bond_order_pcf_data*)malloc(sizeof(bond_order_pcf_data));
    if (data == NULL) return NULL;

    data->num_bins = (int)(max_r / dr);
    data->max_r = dr * data->num_bins;
    data->bin_width = dr;

    data->g_r = (double*)calloc(data->num_bins, sizeof(double));
    data->g6_r = (double*)calloc(data->num_bins, sizeof(double));
    data->r = (double*)malloc(data->num_bins * sizeof(double));
    data->k_vector = (double*)malloc(2 * sizeof(double));
    if (data->g_r == NULL || data->g6_r == NULL || data->r == NULL || data->k_vector == NULL) {
        free_bond_order_pcf_data(data);
        return NULL;
    }
    data->k_vector[0] = k_vector[0];
    data->k_vector[1] = k_vector[1];
    for (int i = 0; i < data->num_bins; i++) {
        data->r[i] = (i + 0.5) * data->bin_width;
    }

    // Parallelize both regular g(r) and bond order g6(r)
    int num_bins = data->num_bins;
    double *g_r_private = (double*)calloc(num_bins, sizeof(double));
    double *g6_r_private = (double*)calloc(num_bins, sizeof(double));
    #pragma omp parallel
    {
        double *g_r_local = (double*)calloc(num_bins, sizeof(double));
        double *g6_r_local = (double*)calloc(num_bins, sizeof(double));
        #pragma omp for nowait
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                PBC(&dx, &dy); // Apply periodic boundary conditions
                double r = sqrt(dx*dx + dy*dy);
                if (r < max_r) {
                    int bin = (int)(r / data->bin_width);
                    if (bin >= 0 && bin < num_bins) {
                        g_r_local[bin] += 1.0;
                        double k_dot_r_i = k_vector[0] * particles[i].x + k_vector[1] * particles[i].y;
                        double k_dot_r_j = k_vector[0] * particles[j].x + k_vector[1] * particles[j].y;
                        double cos_i = cos(k_dot_r_i);
                        double sin_i = sin(k_dot_r_i);
                        double cos_j = cos(k_dot_r_j);
                        double sin_j = sin(k_dot_r_j);
                        double bond_order = cos_i * cos_j + sin_i * sin_j;
                        g6_r_local[bin] += bond_order;
                    }
                }
            }
        }
        #pragma omp critical
        {
            for (int k = 0; k < num_bins; k++) {
                g_r_private[k] += g_r_local[k];
                g6_r_private[k] += g6_r_local[k];
            }
        }
        free(g_r_local);
        free(g6_r_local);
    }
    for (int k = 0; k < num_bins; k++) {
        data->g_r[k] = g_r_private[k];
        data->g6_r[k] = g6_r_private[k];
    }
    free(g_r_private);
    free(g6_r_private);

    // First, compute g6_r as the average bond order per pair in each bin
    for (int i = 0; i < data->num_bins; i++) {
        if (data->g_r[i] > 0) {
            data->g6_r[i] /= data->g_r[i];
        } else {
            data->g6_r[i] = 0.0;
        }
    }

    // Now normalize g(r) only
    double volume = Lx * Ly;
    double density = N / volume;
    for (int i = 0; i < data->num_bins; i++) {
        double r = data->r[i];
        double shell_volume = 2 * M_PI * r * data->bin_width;
        double expected_counts = shell_volume * density * N;
        if (expected_counts > 0) {
            data->g_r[i] /= expected_counts;
        } else {
            data->g_r[i] = 0.0;
        }
    }
    return data;
}

g6corr_data* compute_g6_correlation(particle* particles, int N, double dr, double max_r, double Lx, double Ly) {
    int num_bins = (int)(max_r / dr);
    g6corr_data* data = malloc(sizeof(g6corr_data));
    data->num_bins = num_bins;
    data->max_r = dr * num_bins;
    data->bin_width = dr;
    data->r = malloc(num_bins * sizeof(double));
    data->g6_corr = calloc(num_bins, sizeof(double));
    data->counts = calloc(num_bins, sizeof(int));
    for (int i = 0; i < num_bins; i++) {
        data->r[i] = (i + 0.5) * dr;
    }
    
    boop_data* boop = computeBOOPVoronoi(particles, N, Lx, Ly);
    double complex* psi6 = malloc(N * sizeof(double complex));
    for (int i = 0; i < N; i++) {
        psi6[i] = boop[i].q6 * cexp(I * boop[i].q6_arg);
    }
    // Accumulate correlation
    #pragma omp parallel
    {
        double *g6_corr_private = calloc(num_bins, sizeof(double));
        int *counts_private = calloc(num_bins, sizeof(int));
        #pragma omp for nowait
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                PBC(&dx, &dy);
                double r = sqrt(dx*dx + dy*dy);
                if (r < max_r) {
                    int bin = (int)(r / dr);
                    if (bin >= 0 && bin < num_bins) {
                        double corr = creal(conj(psi6[i]) * psi6[j]);
                        g6_corr_private[bin] += corr;
                        counts_private[bin]++;
                    }
                }
            }
        }
        #pragma omp critical
        {
            for (int k = 0; k < num_bins; k++) {
                data->g6_corr[k] += g6_corr_private[k];
                data->counts[k] += counts_private[k];
            }
        }
        free(g6_corr_private);
        free(counts_private);
    }
    free(psi6);
    // Normalize
    for (int i = 0; i < num_bins; i++) {
        if (data->counts[i] > 0) {
            data->g6_corr[i] /= data->counts[i];
        } else {
            data->g6_corr[i] = 0.0;
        }
    }
    return data;
}

void free_pcf_data(pcf_data* data) {
    if (data) {
        free(data->g_r);
        free(data->r);
        free(data);
    }
}

void free_bond_order_pcf_data(bond_order_pcf_data* data) {
    if (data) {
        free(data->g_r);
        free(data->g6_r);
        free(data->r);
        free(data->k_vector);
        free(data);
    }
}

void free_g6corr_data(g6corr_data* data) {
    if (data) {
        free(data->r);
        free(data->g6_corr);
        free(data->counts);
        free(data);
    }
}

void save_pcf(const char* filename, particle* particles, int N, double dr, double max_r, double Lx, double Ly) {
    pcf_data* data = calculate_pcf(particles, N, dr, max_r, Lx, Ly);
    
    // Check if the file already exists to determine if this is the first save
    FILE* check_file = fopen(filename, "r");
    int first_save = (check_file == NULL);
    if (check_file != NULL) {
        fclose(check_file);
    }
    
    FILE* file;
    if (first_save) {
        // First time saving - create new file with full data
        file = fopen(filename, "w");
        // Write both r and g(r) data
        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file, "%lf ", data->r[i]);
        }
        fprintf(file, "\n");
        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file, "%lf ", data->g_r[i]);
        }
        fprintf(file, "\n");
    } else {
       
        file = fopen(filename, "a");
        
        // Write only g(r) values (not the r values)
        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file, "%lf ", data->g_r[i]);
        }
        fprintf(file, "\n");
    }
    
    fclose(file);
    free_pcf_data(data);
}


void save_pcf_g6(const char* filename, particle* particles, int N, double dr, double max_r, double Lx, double Ly) {
   
    g6corr_data* data = compute_g6_correlation(particles, N, dr, max_r, Lx, Ly);
    
    FILE* check_file = fopen(filename, "r");
    int first_save = (check_file == NULL);
    if (check_file != NULL) {
        fclose(check_file);
    }
    
    FILE* file;
    if (first_save) {
        // First time saving - create new file with full data
        file = fopen(filename, "w");
        // Write both r and g(r) data
        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file, "%lf ", data->r[i]);
        }
        fprintf(file, "\n");
        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file, "%lf ", data->g6_corr[i]);
        }
        fprintf(file, "\n");
    } else {
       
        file = fopen(filename, "a");
        
        // Write only g(r) values (not the r values)
        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file, "%lf ", data->g6_corr[i]);
        }
        fprintf(file, "\n");
    }
    
    fclose(file);
    free_g6corr_data(data);
}


// Save both PCF and boop (q6) in two separate files
void save_pcf_boop(const char* pcf_filename, const char* boop_filename, particle* particles, int N, double dr, double max_r, double Lx, double Ly, double phi) {

    //double t0 = omp_get_wtime();
    double expected_bragg = sqrt(8*M_PI*phi/(sqrt(3)));
    double k_vector[2];
    find_max_structure_factor_bragg(particles, N, Lx, Ly, expected_bragg, k_vector);
    //double t1 = omp_get_wtime();

    //double dt_bragg = t1 - t0;
    //printf("Time to find max structure factor: %.6f s\n", dt_bragg);

    //double t2 = omp_get_wtime();
    bond_order_pcf_data* data = calculate_bond_order_pcf(particles, N, dr, max_r, k_vector, Lx, Ly);
    //double t3 = omp_get_wtime();
    if (data == NULL) return;
    //double dt_bondorder = t3 - t2;
    //printf("Time to calculate bond order PCF: %.6f s\n", dt_bondorder);

    FILE* check_file = fopen(pcf_filename, "r");
    int first_save = (check_file == NULL);
    if (check_file != NULL) {
        fclose(check_file);
    }

    FILE* file_pcf, *file_boop;
    if (first_save) {
        // First time saving - create new file with full data
        file_pcf = fopen(pcf_filename, "w");
        file_boop = fopen(boop_filename, "w");

        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file_boop, "%lf ", data->r[i]);
            fprintf(file_pcf, "%lf ", data->r[i]);
        }
        fprintf(file_pcf, "\n");
        fprintf(file_boop, "\n");
        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file_pcf, "%lf ", data->g_r[i]);
            fprintf(file_boop, "%lf ", data->g6_r[i]);
        }
        fprintf(file_pcf, "\n");
        fprintf(file_boop, "\n");
    } else {
        file_pcf = fopen(pcf_filename, "a");
        file_boop = fopen(boop_filename, "a");
        // Write only g(r) values (not the r values)
        for (int i = 0; i < data->num_bins; i++) {
            fprintf(file_pcf, "%lf ", data->g_r[i]);
            fprintf(file_boop, "%lf ", data->g6_r[i]);
        }
        fprintf(file_pcf, "\n");
        fprintf(file_boop, "\n");
    }

    fclose(file_pcf);
    fclose(file_boop);
    free_bond_order_pcf_data(data);
}

void find_max_structure_factor_bragg(particle* particles, int N, double Lx, double Ly, double expected_bragg, double k[2]) {
    double max_S = -1;
    double best_k[2] = {0.0, 0.0};
    double dkx = 2 * M_PI / Lx;
    double dky = 2 * M_PI / Ly;
    double kx_min = dkx;
    double kx_max = expected_bragg + 0.5;
    double ky_min = dky;
    double ky_max = expected_bragg + 0.5;

    #pragma omp parallel
    {
        double local_max_S = -1;
        double local_k[2] = {0.0, 0.0};
        #pragma omp for collapse(2) schedule(dynamic)
        for (int ikx = 0; ikx < (int)((kx_max - kx_min) / dkx); ikx++) {
            for (int iky = 0; iky < (int)((ky_max - ky_min) / dky); iky++) {
                double kx = kx_min + ikx * dkx;
                double ky = ky_min + iky * dky;
                double k_abs = sqrt(kx*kx + ky*ky);
                double theta = atan2(ky, kx);
                if (k_abs < expected_bragg - 0.5 || k_abs > expected_bragg + 0.5 || theta < 0 || theta > 2*M_PI/5){   
                    continue;
                }
                double re = 0.0, im = 0.0;
                for (int j = 0; j < N; j++) {
                    double phase = kx * particles[j].x + ky * particles[j].y;
                    re += cos(phase);
                    im += sin(phase);
                }
                double S = (re * re + im * im) / N;
                if (S > local_max_S) {
                    local_max_S = S;
                    local_k[0] = kx;
                    local_k[1] = ky;
                }
            }
        }
        #pragma omp critical
        {
            if (local_max_S > max_S) {
                max_S = local_max_S;
                best_k[0] = local_k[0];
                best_k[1] = local_k[1];
            }
        }
    }
    k[0] = best_k[0];
    k[1] = best_k[1];
}


