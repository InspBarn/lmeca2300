/* cahn_hilliard.h */

#ifndef _cahn_hilliard_h_
#define _cahn_hilliard_h_

// Useful libraries
#include <stdio.h>
#include <string.h>
#include <time.h>
// Mathematics library
#include <math.h>
#include <complex.h>
#include <fftw3.h>
// My utilities
#include "utils.h"

#define ind_spa(i,j) i*N+j
#define ind_fft(i,j) i*(N/2+1)+j

typedef struct {
    // Spectral derivation
    int N;
    double *k;

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

#endif /* _cahn_hilliard_h_ */
