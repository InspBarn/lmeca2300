/* fft.c */

#include "fft.h"

int* rearange(const int N)
{
    int *v = (int*)malloc(N*sizeof(int));
    int tgt=0, mask,temp;

    for (int pst=0; pst<N; pst++) {
        if (tgt>=pst) {
            v[tgt] = pst;
            v[pst] = tgt;
        }
        mask = N;
        while (tgt & (mask>>=1))
            tgt &= ~mask;
        tgt |= mask;
    }
    return v;
}

fft_plan* fft_init(int N, int direction)
{
    fft_plan *u = (fft_plan*)malloc(sizeof(fft_plan));
    u->N = N;

    u->twiddle = (double complex*)malloc((N/2)*sizeof(double complex));
    for (int k=0; k<N/2; k++)
        u->twiddle[k] = W((double)(direction*k)/(double)N);

    u->reindexing = rearange(N);
    return u;
}

/*
This function aims to rearange the indexes of a vector such that all
even indexes are then in the first half of the vector, and all odd indexes
are after. Here is an example for N = 8:
    0, 1, 2, 3, 4, 5, 6, 7
    ↓  ↓  ↓  ↓  ↓  ↓  ↓  ↓
    0, 4, 2, 6, 1, 5, 3, 7

    0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15
    ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓   ↓
    0,  8,  4, 12,  2, 10,  6, 14,  1,  9,  5, 13,  3, 11,  7, 15
*/
inline void rearange_fft(double complex *vec, double complex *sol, const int N)
{
    int tgt = 0, mask;
    double complex temp;
    for (int pst=0; pst<N; pst++) {
        if (tgt>=pst) {
            sol[tgt] = vec[pst];
            sol[pst] = vec[tgt];
        }
        mask = N;
        while (tgt & (mask>>=1))
            tgt &= ~mask;
        tgt |= mask;
    }
}

/*
The vector function has been rearanged according to the radix-2 fft algorithm. This function
aims to compute the fft based on this rearrangement.

Arguments
-------------------
u : function of length N on which we apply the fft algorithm
twiddle : helper, complex exponential of length N : exp(-2i*π * k/N), k = 0,...,(N-1)
direction : +/-1, specify if the fft is computed forward (+1) or backward (-1)
*/
inline void compute_fft(double complex *u, double complex *twiddle, const int N)
{
    int i,j,k,l,m,s; double complex E,O;

    // printf("\n");
    for (s=1; s<N; s<<=1) {
        m = s<<1;
        // printf("%d,%d\t", s,m);
        for (i=0; i<N; i+=m) {
            // printf("%d -- ", i);
            for (j=0; j<s; j++) {
                k = i+j; l = k+s;
                switch (j) {
                    case (0):
                        O = u[l];
                    default:
                        O = u[l]*twiddle[j*N/m];
                } E = u[k];
                u[k] = E + O;
                u[l] = E - O;
                // printf("(%d,%d,%d,%d), ", j,k,l,j*N/m);
            }
            // printf("\n\t");
        }
        // printf("\n");
    }
}

void fft(double complex *u, double complex *u_hat, fft_plan *plan)
{
    for (int k=0; k<plan->N; k++)
        u_hat[k] = u[plan->reindexing[k]];

    compute_fft(u_hat,plan->twiddle,plan->N);
}

void ifft(double complex *u, double complex *u_hat, fft_plan *plan)
{
    for (int k=0; k<plan->N; k++)
        u[k] = u_hat[plan->reindexing[k]];

    compute_fft(u,plan->twiddle,plan->N);

    for (int n=0; n<plan->N; n++)
        u[n] /= (double)plan->N;
}

void fft2(double complex *u, double complex *u_hat, fft_plan *plan)
{
    int i,j, N=plan->N;
    double complex *u_sub = (double complex*)malloc(N*sizeof(double complex));

    // #pragma omp for schedule(static)
    for (j=0; j<N; j++)
        fft(&u[j*N],&u_hat[j*N],plan);

    // #pragma omp parallel for // private(i,j,u_sub) default(shared)
    for (int i=0; i<N; i++) {
        // printf("We are now working on thread number %d\n", omp_get_thread_num());
        for (int j=0; j<N; j++)
            u_sub[j] = u_hat[i+plan->reindexing[j]*N];

        compute_fft(u_sub,plan->twiddle,N);

        for (int j=0; j<N; j++)
            u_hat[i+j*N] = u_sub[j];
    }

    free(u_sub);
}

void ifft2(double complex *u, double complex *u_hat, fft_plan *plan)
{
    int N = plan->N;
    double complex *u_sub = (double complex*)malloc(N*sizeof(double complex));

    // #pragma omp parallel for
    for (int j=0; j<N; j++)
        ifft(&u[j*N],&u_hat[j*N],plan);

    // #pragma omp parallel for
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++)
            u_sub[j] = u[i+plan->reindexing[j]*N];

        compute_fft(u_sub,plan->twiddle,N);

        for (int j=0; j<N; j++)
            u[i+j*N] = u_sub[j]/(double)N;
    } free(u_sub);
}

/* --------------------------------------------------------------------------------------------- */
// void compute_rfft(double complex *u, const int N, const int direction)
// {
//     int k,l,m,n; double complex E,O;
// 
//     // printf("\n");
//     for (int s=1; s<N; s<<=1) {
//         m = s<<1;
//         // printf("%d,%d\t", s,m);
//         for (int i=0; i<N/2+1; i+=m) {
//             // printf("%d -- ", i);
//             for (int j=0; j<s; j++) {
//                 k = i+j; l = k+s;
//                 switch (direction) {
//                     case (+1):
//                         E = u[k]; O = u[l]*W(j); //twiddle[j*N/m];
//                         break;
//                     case (-1):
//                         E = u[k]; O = u[l]/W(j); //twiddle[j*N/m];
//                         break;
//                     default:
//                         E = u[k]; O = u[l]*W(j); //twiddle[j*N/m];
//                 }
//                 u[k] = E + O;
//                 u[l] = E - O;
//                 // printf("(%d,%d,%d,%d), ", j,k,l,j*N/m);
//             }
//             // printf("\n\t");
//         }
//         // printf("\n");
//     }
// }
// 
// void rfft(double complex *u, double complex *u_hat, const int N)
// {
//     rearange_fft(u,u_hat,N);
//     compute_rfft(u_hat,N,+1);
// 
//     for (int i=1; i<N/2; i++)
//         u_hat[i] = u_hat[i] + I*u_hat[N/2];
// }

/* --------------------------------------------------------------------------------------------- */
double* fftfreq(const int N)
{
    double *freq = (double*)malloc(N*sizeof(double));
    for (int i=0; i<N/2; i++) {
        freq[i] = (double)i;
        freq[N/2+i] = (double)(i-N/2);
    }
    return freq;
}
