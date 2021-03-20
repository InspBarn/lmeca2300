/* fft.h */

#ifndef FFT

#include <stdlib.h>
#include <math.h>
#include <complex.h>

// #include "omp.h"

#define PI acos(-1)

#define FORWARD 1
#define BACKWARD -1

typedef struct{
    int N;
    double complex *twiddle;
    int *reindexing;
} fft_plan;

//
fft_plan* fft_init(int N, int direction);

// Standard FFTs
void fft(double complex *u, double complex *u_hat, fft_plan *plan);
void fft2(double complex *u, double complex *u_hat, fft_plan *plan);
void ifft(double complex *u, double complex *u_hat, fft_plan *plan);
void ifft2(double complex *u, double complex *u_hat, fft_plan *plan);

// Real FFTs
void rfft(double complex *u, double complex *u_hat, const int N);
void rfft2(double complex *u, double complex *u_hat, const int M, const int N);
void irfft(double complex *u, double complex *u_hat, const int N);
void irfft2(double complex *u, double complex *u_hat, const int M, const int N);

// Helper Routines
double* fftfreq(const int N);
void rfftfreq(double *k, const int N);

//
#define W(exp) cexp(-2.0*I*PI * exp)

#define FFT
#endif
