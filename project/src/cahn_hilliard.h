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

#define ind_spa(i,j) i*N+j
#define ind_fft(i,j) i*(N/2+1)+j

typedef struct {
    // Spectral derivation
    int N;
    double *laplace;

    double *rval;
    fftw_complex *cval;
    fftw_plan forward,backward;

    // Simulation parameters
    double xMin, xMax, h;
    double tMax, dt, t;
    double a;
} cahn_hilliard;

cahn_hilliard* cahn_hilliard_init(const int N);
void cahn_hilliard_free(cahn_hilliard* problem);
void cahn_hilliard_solve(cahn_hilliard *problem, double *u);

#define CAHN
#endif
