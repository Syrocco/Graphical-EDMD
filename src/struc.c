#include "EDMD.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

int nqx, nqy, N;
double *qx, *qy;
double *structFactor;
double complex *structFactorComplex;
FILE *fileout;


/*
 * Transverse velocity autocorrelation configuration.
 *
 * N_p samples, including both endpoints, cover the interval t_p.  When this
 * measurement is enabled, EDMD therefore samples thermodynamics every
 * t_p/(N_p - 1).  The default mode below is |k| = 3*2*pi/L for a square box;
 * all reciprocal-lattice vectors with the same physical |k| are averaged.
 */
#ifndef TVACF_ENABLED
#define TVACF_ENABLED 1
#endif
#ifndef TVACF_MODE_M
#define TVACF_MODE_M 3
#endif
#ifndef TVACF_MODE_N
#define TVACF_MODE_N 0
#endif
#ifndef TVACF_N_P
#define TVACF_N_P 101
#endif
#ifndef TVACF_T_P
#define TVACF_T_P 3000.0
#endif
#ifndef TVACF_CHECKPOINT_INTERVAL
#define TVACF_CHECKPOINT_INTERVAL 100
#endif

typedef struct transverse_mode {
    int m;
    int n;
    double kx;
    double ky;
    double tx;
    double ty;
} transverse_mode;

static transverse_mode *tvacfModes;
static int tvacfModeCount;
static int tvacfParticleCount;
static double complex *tvacfHistory;
static double complex *tvacfCorrelationSum;
static double *tvacfCorrelationM2;
static unsigned long long *tvacfOriginCount;
static unsigned long long tvacfSampleCount;
static int tvacfNextSlot;
static FILE *tvacfFile;
static char *tvacfFileName;
static double tvacfLx;
static double tvacfLy;

int transverseVelocityCorrelationIsEnabled(void){
    return TVACF_ENABLED;
}

double transverseVelocityCorrelationSampleInterval(void){
    if (TVACF_N_P < 2 || TVACF_T_P <= 0.0)
        return 0.0;
    return TVACF_T_P/(TVACF_N_P - 1);
}

void configureTransverseVelocityCorrelationSampling(double *thermoInterval,
                                                     int *useLogSpacing){
    if (!TVACF_ENABLED)
        return;

    double sampleInterval = transverseVelocityCorrelationSampleInterval();
    if (sampleInterval <= 0.0){
        fprintf(stderr,
                "ERROR: TVACF_N_P must be at least 2 and TVACF_T_P must be positive.\n");
        exit(EXIT_FAILURE);
    }

    if (*useLogSpacing){
        fprintf(stderr,
                "WARNING: TVACF requires uniform sampling; disabling logarithmic thermo spacing.\n");
        *useLogSpacing = 0;
    }
    if (fabs(*thermoInterval - sampleInterval) >
        1e-12*fmax(1.0, sampleInterval)){
        fprintf(stderr,
                "TVACF: setting dtimeThermo to %.17g so %d points span %.17g.\n",
                sampleInterval, TVACF_N_P, (double)TVACF_T_P);
        *thermoInterval = sampleInterval;
    }
}

static int transverseModeIsOnTargetShell(int m, int n, double targetK2,
                                         double Lx, double Ly){
    double kx = 2.0*M_PI*m/Lx;
    double ky = 2.0*M_PI*n/Ly;
    double k2 = kx*kx + ky*ky;
    return fabs(k2 - targetK2) <= 1e-10*targetK2;
}

static void writeTransverseVelocityCorrelation(void){
    if (tvacfFileName == NULL)
        return;

    if (tvacfFile == NULL){
        tvacfFile = fopen(tvacfFileName, "w");
    }
    else{
        tvacfFile = freopen(tvacfFileName, "w", tvacfFile);
    }
    if (tvacfFile == NULL){
        fprintf(stderr, "WARNING: could not write TVACF file %s.\n", tvacfFileName);
        return;
    }

    double targetKx = 2.0*M_PI*TVACF_MODE_M/tvacfLx;
    double targetKy = 2.0*M_PI*TVACF_MODE_N/tvacfLy;
    double targetK2 = targetKx*targetKx + targetKy*targetKy;
    fprintf(tvacfFile, "# Transverse velocity autocorrelation\n");
    fprintf(tvacfFile, "# target_mode_m %d target_mode_n %d\n",
            TVACF_MODE_M, TVACF_MODE_N);
    fprintf(tvacfFile, "# target_kx %.4lf target_ky %.4lf k_squared %.4lf\n",
            targetKx, targetKy, targetK2);
    fprintf(tvacfFile, "# sample_interval %.4lf requested_span %.4lf samples_seen %llu\n",
            transverseVelocityCorrelationSampleInterval(), (double)TVACF_T_P,
            tvacfSampleCount);
    fprintf(tvacfFile, "# isotropic_mode_count %d\n", tvacfModeCount);
    for (int q = 0; q < tvacfModeCount; q++){
        fprintf(tvacfFile, "# mode %d %d %.4lf %.4lf\n",
                tvacfModes[q].m, tvacfModes[q].n,
                tvacfModes[q].kx, tvacfModes[q].ky);
    }
    fprintf(tvacfFile,
            "# u_perp(k,t)=N^(-1/2) sum_j [v_j(t).e_perp(k)] exp[-i k.r_j(t)]\n");
    fprintf(tvacfFile,
            "# variance is the sample variance across time origins after isotropic averaging\n");
    fprintf(tvacfFile, "# columns: lag_time C_T_real variance\n");

    int lagCount = (tvacfSampleCount < (unsigned long long)TVACF_N_P)
                 ? (int)tvacfSampleCount : TVACF_N_P;
    double dt = transverseVelocityCorrelationSampleInterval();
    for (int lag = 0; lag < lagCount; lag++){
        double correlation = 0.0;
        double variance = 0.0;
        if (tvacfOriginCount[lag] > 0){
            correlation = creal(tvacfCorrelationSum[lag])/
                          ((double)tvacfOriginCount[lag]*tvacfModeCount);
        }
        if (tvacfOriginCount[lag] > 1){
            variance = tvacfCorrelationM2[lag]/(tvacfOriginCount[lag] - 1);
        }
        fprintf(tvacfFile, "%.4lf %.4lf %.4lf\n",
                lag*dt, correlation, variance);
    }
    fflush(tvacfFile);
}

void initTransverseVelocityCorrelation(double Lx, double Ly, int numParticles,
                                       const char *outputName){
    if (!TVACF_ENABLED)
        return;
    if (TVACF_MODE_M == 0 && TVACF_MODE_N == 0){
        fprintf(stderr, "ERROR: the TVACF wavevector cannot be k=(0,0).\n");
        exit(EXIT_FAILURE);
    }
    if (Lx <= 0.0 || Ly <= 0.0 || numParticles <= 0 || outputName == NULL){
        fprintf(stderr, "ERROR: invalid TVACF initialization parameters.\n");
        exit(EXIT_FAILURE);
    }

    tvacfLx = Lx;
    tvacfLy = Ly;
    tvacfParticleCount = numParticles;
    double targetKx = 2.0*M_PI*TVACF_MODE_M/Lx;
    double targetKy = 2.0*M_PI*TVACF_MODE_N/Ly;
    double targetK2 = targetKx*targetKx + targetKy*targetKy;
    double targetK = sqrt(targetK2);
    int maxM = (int)ceil(targetK*Lx/(2.0*M_PI) + 1e-12);
    int maxN = (int)ceil(targetK*Ly/(2.0*M_PI) + 1e-12);

    for (int m = -maxM; m <= maxM; m++){
        for (int n = -maxN; n <= maxN; n++){
            if ((m != 0 || n != 0) &&
                transverseModeIsOnTargetShell(m, n, targetK2, Lx, Ly)){
                tvacfModeCount++;
            }
        }
    }
    if (tvacfModeCount == 0){
        fprintf(stderr, "ERROR: no reciprocal-lattice vectors found for the TVACF shell.\n");
        exit(EXIT_FAILURE);
    }

    tvacfModes = calloc((size_t)tvacfModeCount, sizeof(*tvacfModes));
    tvacfHistory = calloc((size_t)TVACF_N_P*tvacfModeCount,
                          sizeof(*tvacfHistory));
    tvacfCorrelationSum = calloc((size_t)TVACF_N_P,
                                 sizeof(*tvacfCorrelationSum));
    tvacfCorrelationM2 = calloc((size_t)TVACF_N_P,
                                sizeof(*tvacfCorrelationM2));
    tvacfOriginCount = calloc((size_t)TVACF_N_P, sizeof(*tvacfOriginCount));
    tvacfFileName = malloc(strlen(outputName) + 1);
    if (tvacfModes == NULL || tvacfHistory == NULL ||
        tvacfCorrelationSum == NULL || tvacfCorrelationM2 == NULL ||
        tvacfOriginCount == NULL || tvacfFileName == NULL){
        fprintf(stderr, "ERROR: unable to allocate TVACF storage.\n");
        exit(EXIT_FAILURE);
    }
    strcpy(tvacfFileName, outputName);

    int q = 0;
    for (int m = -maxM; m <= maxM; m++){
        for (int n = -maxN; n <= maxN; n++){
            if ((m == 0 && n == 0) ||
                !transverseModeIsOnTargetShell(m, n, targetK2, Lx, Ly)){
                continue;
            }
            double kx = 2.0*M_PI*m/Lx;
            double ky = 2.0*M_PI*n/Ly;
            double k = sqrt(kx*kx + ky*ky);
            tvacfModes[q].m = m;
            tvacfModes[q].n = n;
            tvacfModes[q].kx = kx;
            tvacfModes[q].ky = ky;
            tvacfModes[q].tx = -ky/k;
            tvacfModes[q].ty = kx/k;
            q++;
        }
    }

    writeTransverseVelocityCorrelation();
}

void updateTransverseVelocityCorrelation(const particle *particles,
                                         double sampleTime){
    (void)sampleTime;
    if (!TVACF_ENABLED || tvacfHistory == NULL)
        return;

    int currentSlot = tvacfNextSlot;
    double normalization = 1.0/sqrt((double)tvacfParticleCount);
    for (int q = 0; q < tvacfModeCount; q++){
        double re = 0.0;
        double im = 0.0;
        const transverse_mode *mode = tvacfModes + q;
        for (int j = 0; j < tvacfParticleCount; j++){
            const particle *p = particles + j;
            double projectedVelocity = p->vx*mode->tx + p->vy*mode->ty;
            double phase = mode->kx*p->x + mode->ky*p->y;
            re += projectedVelocity*cos(phase);
            im -= projectedVelocity*sin(phase);
        }
        tvacfHistory[(size_t)currentSlot*tvacfModeCount + q] =
            normalization*(re + I*im);
    }

    int available = (tvacfSampleCount + 1 < (unsigned long long)TVACF_N_P)
                  ? (int)(tvacfSampleCount + 1) : TVACF_N_P;
    for (int lag = 0; lag < available; lag++){
        int pastSlot = currentSlot - lag;
        if (pastSlot < 0)
            pastSlot += TVACF_N_P;
        double complex shellCorrelation = 0.0;
        for (int q = 0; q < tvacfModeCount; q++){
            double complex current =
                tvacfHistory[(size_t)currentSlot*tvacfModeCount + q];
            double complex past =
                tvacfHistory[(size_t)pastSlot*tvacfModeCount + q];
            shellCorrelation += current*conj(past);
        }
        double sample = creal(shellCorrelation)/tvacfModeCount;
        double oldMean = 0.0;
        if (tvacfOriginCount[lag] > 0){
            oldMean = creal(tvacfCorrelationSum[lag])/
                      ((double)tvacfOriginCount[lag]*tvacfModeCount);
        }
        tvacfCorrelationSum[lag] += shellCorrelation;
        tvacfOriginCount[lag]++;
        double newMean = creal(tvacfCorrelationSum[lag])/
                         ((double)tvacfOriginCount[lag]*tvacfModeCount);
        tvacfCorrelationM2[lag] += (sample - oldMean)*(sample - newMean);
    }

    tvacfSampleCount++;
    tvacfNextSlot = (currentSlot + 1) % TVACF_N_P;
    if (TVACF_CHECKPOINT_INTERVAL > 0 &&
        tvacfSampleCount % TVACF_CHECKPOINT_INTERVAL == 0){
        writeTransverseVelocityCorrelation();
    }
}

void finalizeTransverseVelocityCorrelation(void){
    if (!TVACF_ENABLED)
        return;
    writeTransverseVelocityCorrelation();
    if (tvacfFile != NULL){
        fclose(tvacfFile);
        tvacfFile = NULL;
    }
    free(tvacfFileName);
    free(tvacfOriginCount);
    free(tvacfCorrelationM2);
    free(tvacfCorrelationSum);
    free(tvacfHistory);
    free(tvacfModes);
    tvacfFileName = NULL;
    tvacfOriginCount = NULL;
    tvacfCorrelationM2 = NULL;
    tvacfCorrelationSum = NULL;
    tvacfHistory = NULL;
    tvacfModes = NULL;
    tvacfModeCount = 0;
}


void initStructureFactor(double qmax, double Lx, double Ly, int num, FILE* file,  double doFQT){
    double xx = 2*M_PI/Lx;
    double yy = 2*M_PI/Ly;
    nqx = 2*qmax/xx + 1;
    nqy = 2*qmax/yy + 1;
    qx = calloc(nqx, sizeof(double));
    qy = calloc(nqy, sizeof(double));
    
    N = num;
    for (int i = 0; i < nqx; ++i){
        qx[i] = xx*(i - (nqx - 1)/2);
    }
    for (int i = 0; i < nqy; ++i){
        qy[i] = yy*(i - (nqy - 1)/2);
    }


    fileout = file;
    
    fprintf(file, "-%lf \n", doFQT);
    for (int i = 0; i < nqx; ++i)
        fprintf(fileout, "%g ", qx[i]);
    fprintf(fileout, "\n");

    for (int j = 0; j < nqy; ++j)
        fprintf(fileout, "%g ", qy[j]);
    fprintf(fileout, "\n");
    
    if (doFQT)
        structFactorComplex = calloc(nqy*nqx, sizeof(double complex));
    else
        structFactor = calloc(nqy*nqx, sizeof(double));
}

void computeStructureFactor(particle* particles){
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nqx; i++){
        for (int j = 0; j < nqy; j++){
            double im = 0;
            double re = 0;

            for (int n = 0; n < N; n++){
                particle* p = particles + n;
                double qr = qx[i]*p->x + qy[j]*p->y;

                re += cos(qr);
                im += sin(qr);
            }
            if (structFactorComplex != NULL)    
                structFactorComplex[i*nqy + j] = re + I*im;
            else
                structFactor[i*nqy + j] = (re*re + im*im)/N;
        }
    }
}

void computeVelocityStructureFactor(particle* particles){
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nqx; i++){
        for (int j = 0; j < nqy; j++){
            double im = 0;
            double re = 0;

            for (int n = 0; n < N; n++){
                particle* p = particles + n;
                double qr = qx[i]*p->x + qy[j]*p->y;

                re += p->vx*cos(qr) + p->vy*sin(qr);
                im += p->vx*sin(qr) - p->vy*cos(qr);
            }
            if (structFactorComplex != NULL)    
                structFactorComplex[i*nqy + j] = re + I*im;
            else
                structFactor[i*nqy + j] = (re*re + im*im)/N;
        }
    }
}
void saveStructureFactor(particle* particles, int mode){
    if (mode == 0)
        computeVelocityStructureFactor(particles);
    else
        computeStructureFactor(particles);

    for (int i = 0; i < nqx; ++i){
        for (int j = 0; j < nqy; ++j){
            if (structFactorComplex != NULL)
                fprintf(fileout, "(%g+%gj) ", creal(structFactorComplex[i*nqy + j]), cimag(structFactorComplex[i*nqy + j]));
            else
                fprintf(fileout, "%g ", structFactor[i*nqy + j]);
        }
        fprintf(fileout, "\n");
    }
}
