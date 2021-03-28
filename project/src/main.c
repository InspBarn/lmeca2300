/* main.c */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

// #include "fft.h"
#include "cahn_hilliard.h"

int main(int argc, char *argv[]){
    clock_t start,end;

    int N = 8;

//    fft_plan *plan  = fft_init(N, FORWARD);
//    fft_plan *planb = fft_init(N, BACKWARD);
//
//    double complex *u = (double complex*)calloc(N, sizeof(double complex));
//    double complex *u_hat = (double complex*)calloc(N, sizeof(double complex));
//    for (int i=0; i<N; i++) {
//        u[i] = exp(sin(2*PI*(double)i/(double)N));
//        printf("%.2e,\t", creal(u[i]));
//    } printf("\n");
//
//    fft(u,u_hat, plan);
//    // ifft(u,u_hat,N);
//
//    for(int i=0; i<N; i++) {
//        printf("%.2e,\t", creal(u[i]));
//        // printf("%.2e + %.2ei\t", creal(u_hat[i]),cimag(u_hat[i]));
//    } printf("\n");
//
//    double complex *u = (double complex*)calloc(N*N, sizeof(double complex));
//    double complex *u_hat = (double complex*)calloc(N*N, sizeof(double complex));
//
//    for (int i=0; i<N; i++) {
//        for (int j=0; j<N; j++) {
//            u[i+j*N] = exp(sin(2*PI*(double)i/(double)N)) * cos(2*PI*(double)j/(double)N);
//        }
//    }
//
//    int max = 100;
//    double time = 0.0;
//    for (int i=0; i<max; i++) {
//        for (int j=0; j<N*N; j++)
//            u[j] = rand() / (double)RAND_MAX;
//        start = clock();
//        fft2(u,u_hat,plan);
//        end = clock();
//        time = time + (double)(end-start)/CLOCKS_PER_SEC*1e3/(double)max;
//    }
//    printf("Time = %7.3f [ms]\n", time);
//    ifft2(u,u_hat,N,N);
//    for(int j=0; j<N; j++) {
//        for (int i=0; i<N; i++) {
//            printf("%.2e + %.2ei\t", creal(u[i+j*N]),cimag(u[i+j*N]));
//            // printf("%.2e + %.2ei\t", creal(u_hat[i+j*N]),cimag(u_hat[i+j*N]));
//        }
//        printf("\n");
//    }
//
//    free(u); free(u_hat);
//    free(plan); free(planb);

    cahn_hilliard *problem = cahn_hilliard_init(N);
    cahn_hilliard_solve(problem);
    cahn_hilliard_free(problem);

    exit(0);
}
