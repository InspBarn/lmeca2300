/* main.c */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

// #include "fft.h"
#include "cahn_hilliard.h"

int main(int argc, char *argv[])
{
    // Give a bit of entropy for the seed of rand()
    // or it will always be the same sequence
    int seed = (int) time(NULL);
    // srand(seed);

    // We print the seed so you can get the distribution of points back
    printf("seed=%d\n", seed);

    int N = 128;

    double *u = (double*) malloc(N*N*sizeof(double));
    for (int n=0; n<N*N; n++) {
        u[n] = (rand() / (double)RAND_MAX - 0.5)*2.0;
    }

    cahn_hilliard *problem = cahn_hilliard_init(N);
    cahn_hilliard_solve(problem, u);
    cahn_hilliard_free(problem);

    exit(0);
}
