/* cahn_hilliard.h */

#ifndef CAHN

#include <stdio.h>
// #include <stdlib.h>
// #include <unistd.h>
#include <time.h>

#include "fft.h"
#include "plot.h"

typedef struct {
    // Spectral derivation
    int N;
    double complex *laplace;
    fft_plan *forward_plan;
    fft_plan *backward_plan;

    // Simulation parameters
    double xMin, xMax, h;
    double tMax, dt;
    double a;

    // Concentration
    double complex *u;
} cahn_hilliard;

cahn_hilliard* cahn_hilliard_init(const int N);
void cahn_hilliard_solve(cahn_hilliard *problem);

#define CAHN
#endif
