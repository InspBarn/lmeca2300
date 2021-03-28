/* cahn_hilliard.h */

#ifndef CAHN

#include <stdio.h>
// #include <stdlib.h>
// #include <unistd.h>
#include <time.h>
#include <math.h>

// #include "fft.h"
#include <complex.h>
#include <fftw3.h>

#include "plot.h"
#include "utils.h"

#define ind_spa(i,j) i+j*N
#define ind_fft(i,j) i+j*(N/2+1)

typedef struct {
    // Spectral derivation
    int N;
    double *laplace;

    double *in;
    fftw_complex *out;
    fftw_plan forward,backward;

    double *nonlin_in;
    fftw_complex *nonlin_out;
    fftw_plan nonlin_forward,nonlin_backward;

    // Simulation parameters
    double xMin, xMax, h;
    double tMax, dt, t;
    double a;
} cahn_hilliard;

cahn_hilliard* cahn_hilliard_init(const int N);
void cahn_hilliard_free(cahn_hilliard* problem);
void cahn_hilliard_solve(cahn_hilliard *problem);

#define CAHN
#endif
