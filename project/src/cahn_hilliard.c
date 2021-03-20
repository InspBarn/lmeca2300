/* cahn_hilliard.c */

#include "cahn_hilliard.h"

void step_derivative(cahn_hilliard *problem, double t, double complex *u, double complex *du)
{
    const int N = problem->N;
    double complex *u_spat = (double complex*)malloc(N*N*sizeof(double complex));   
    double complex *u_cube = (double complex*)malloc(N*N*sizeof(double complex));

    ifft2(u_spat,u,problem->backward_plan);
    for (int k=0; k<N*N; k++) {
        u_spat[k] = cpow(u_spat[k], 3.0);
    } fft2(u_spat,u_cube,problem->forward_plan);

    double a_sqared = pow(problem->a,2.0);
    for (int k=0; k<N*N; k++) {
        du[k] = problem->laplace[k]*(u_cube[k] - u[k] + a_sqared*problem->laplace[k]*u[k]);
    }

    free(u_spat); free(u_cube);
}

void rk4(cahn_hilliard *problem, double t, double complex *u)
{
    const int N = problem->N;
    const double dt = problem->dt;

    double complex *k1 = (double complex*)malloc(N*N*sizeof(double complex));
    double complex *k2 = (double complex*)malloc(N*N*sizeof(double complex));
    double complex *k3 = (double complex*)malloc(N*N*sizeof(double complex));
    double complex *k4 = (double complex*)malloc(N*N*sizeof(double complex));
    double complex *temp = (double complex*)malloc(N*N*sizeof(double complex));

    step_derivative(problem,t,u,k1);

    for (int k=0; k<N*N; k++) {
        temp[k] = u[k] + dt/2.0 * k1[k];
    } step_derivative(problem,t+dt/2.0,temp,k2);

    for (int k=0; k<N*N; k++) {
        temp[k] = u[k] + dt/2.0 * k2[k];
    } step_derivative(problem,t+dt/2.0,temp,k3);

    for (int k=0; k<N*N; k++) {
        temp[k] = u[k] + dt * k3[k];
    } step_derivative(problem,t+dt,temp,k4);

    for (int k=0; k<N*N; k++) {
        u[k] += (k1[k] + k2[k] + k3[k] + k4[k]) * dt/6.0;
    }

    free(k1); free(k2);
    free(k3); free(k4);
    free(temp);
}

void cahn_hilliard_solve(cahn_hilliard *c)
{
    int N = c->N; int t=0;
    int Ntime = (int) (c->tMax / c->dt);
    double complex *u_hat = (double complex*)malloc(N*N*sizeof(double complex));
    fft2(c->u,u_hat,c->forward_plan);

    bov_window_t *window = bov_window_new(800, 800, "Cahn Hilliard Simulation");
    bov_window_set_color(window, (GLfloat[]){0.9f,0.85f,0.8f,1.0f});

    clock_t start,end;
    do {
        bov_window_update(window);

        // start = clock();
        // imshow(window, c->u, N,N);
        // end = clock();
        // printf("time : %.3f\n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i=0; i<4; i++) {
            rk4(c, t*c->dt, u_hat);
            t++;
        }

        printf("iter : %d -- time : %.6f\r", t,t*c->dt);

    } while(!bov_window_should_close(window));

    bov_window_delete(window);
    free(u_hat);
}

cahn_hilliard* cahn_hilliard_init(const int N)
{
    cahn_hilliard *problem = (cahn_hilliard*)malloc(sizeof(cahn_hilliard));
    problem->N = N;
    problem->xMin = 0.0;
    problem->xMax = 1.0;
    problem->tMax = 6e-3;
    problem->dt = 2.5e-7;
    problem->a = 0.01;

    problem->h = (problem->xMax - problem->xMin) / (double)N;

    problem->forward_plan = fft_init(problem->N,FORWARD);
    problem->backward_plan = fft_init(problem->N,BACKWARD);

    double xdfreq;
    double *freq = fftfreq(N);
    double complex constant = cpow(2.0*PI*I,2.0);
    problem->laplace = (double complex*)malloc(N*N*sizeof(double complex));
    for (int i=0; i<N; i++) {
        xdfreq = cpow(freq[i],2.0);
        for (int j=0; j<N; j++) {
            problem->laplace[i+j*N] = (xdfreq + cpow(freq[j],2.0)) * constant;
        }
    }

    problem->u = (double complex*)malloc(problem->N*problem->N*sizeof(double complex));
    for (int i=0; i<problem->N*problem->N; i++)
        problem->u[i] = (rand() / (double)RAND_MAX - 1.0)*2.0;
    return problem;
}
