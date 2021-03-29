/* cahn_hilliard.c */

#include "cahn_hilliard.h"

void normalize(double *u, const int N)
{
    for (int n=0; n<N*N; n++)
        u[n] /= (double)(N*N);
}

void step_derivative(cahn_hilliard *c_h, fftw_complex *u_hat, fftw_complex *u_dot_hat)
{
    const int N = c_h->N;

    memcpy(c_h->cval,u_hat, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);

    for (int n=0; n<N*N; n++)
        c_h->rval[n] = c_h->rval[n]*c_h->rval[n]*c_h->rval[n];
    fftw_execute(c_h->forward);

    double a_squared = c_h->a*c_h->a;
    for (int k=0; k<N*(N/2+1); k++)
        u_dot_hat[k] = c_h->laplace[k] * (c_h->cval[k] - u_hat[k] - a_squared*c_h->laplace[k]*u_hat[k]);
}

void rk4(cahn_hilliard *c_h)
{
    const int N = c_h->N;
    const double dt = c_h->dt;

    fftw_complex *k1 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *k2 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *k3 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *k4 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *u_hat = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));

    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));
    step_derivative(c_h,temp,k1);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt/2.0 * k1[k];
    step_derivative(c_h,u_hat,k2);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt/2.0 * k2[k];
    step_derivative(c_h,u_hat,k3);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt * k3[k];
    step_derivative(c_h,u_hat,k4);

    for (int k=0; k<N*(N/2+1); k++)
        c_h->cval[k] = temp[k] + (k1[k] + 2.0*k2[k] + 2.0*k3[k] + k4[k]) * dt/6.0;

    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);
    memcpy(c_h->cval,temp, N*(N/2+1)*sizeof(fftw_complex));

    fftw_free(k1); fftw_free(k2);
    fftw_free(k3); fftw_free(k4);
    fftw_free(temp); fftw_free(u_hat);
}

void cnicolson(cahn_hilliard *c_h)
{
    const int N = c_h->N;
    const double dt = c_h->dt;

    fftw_complex *progressive = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *regressive = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *u_hat = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));

    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));
    memcpy(u_hat,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));

    step_derivative(c_h,temp,progressive);
    for (int k=0; k<N*(N/2+1); k++)
        temp[k] += progressive[k]*dt/2.0;
    step_derivative(c_h,temp,regressive);

    for (int k=0; k<N*(N/2+1); k++)
        c_h->cval[k] = u_hat[k] + (progressive[k] + regressive[k]) * dt/2.0;

    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);
    memcpy(c_h->cval,temp, N*(N/2+1)*sizeof(fftw_complex));

    fftw_free(progressive);
    fftw_free(regressive);
    fftw_free(temp);
    fftw_free(u_hat);
}

void cahn_hilliard_solve(cahn_hilliard *c_h, double *u)
{
    c_h->dt /= 4.0;
    int N = c_h->N, t=0, Ntime = (int) (c_h->tMax/c_h->dt);

    memcpy(c_h->rval,u, N*N*sizeof(double));
    fftw_execute(c_h->forward);

    bov_window_t *window = bov_window_new(800, 800, "Cahn Hilliard Simulation");
    // bov_window_set_color(window, (GLfloat[]){0.9f,0.85f,0.8f,1.0f});
    window->param.zoom = 2.0/((double)N-1.0);
    window->param.translate[0] = -((double)N-1.0)/2.0;
    window->param.translate[1] = -((double)N-1.0)/2.0;

    clock_t start,end;
    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    do {
        printf("iteration : %5d / %d -- ", t,Ntime);
        bov_window_update(window);

        start = clock();
        imshow(window, c_h->rval, N,N);
        end = clock();
        printf("time for drawing : %.3f [ms] -- ", (double)(end-start)/CLOCKS_PER_SEC*1e3);

        start = clock();
        for (int i=0; i<40; i++) {
            c_h->t += c_h->dt;
            rk4(c_h);
            t++;
        } end = clock();

        printf("time for integration : %.3f [ms] -- ", (double)(end-start)/CLOCKS_PER_SEC*1e3);
        printf("time : %.6f\n", c_h->t);
        if (t>Ntime)
            break;

    } while(!bov_window_should_close(window));

    bov_window_delete(window);
}

cahn_hilliard* cahn_hilliard_init(const int N)
{
    cahn_hilliard *c_h = (cahn_hilliard*) malloc(sizeof(cahn_hilliard));
    c_h->N = N;
    c_h->xMin = 0.0;
    c_h->xMax = 1.0;
    c_h->tMax = 12e-3;
    c_h->dt = 1e-6;
    c_h->a = 0.01;

    c_h->h = (c_h->xMax - c_h->xMin) / (double)N;

    c_h->rval = (double*) fftw_malloc(N*N*sizeof(double));
    c_h->cval = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    c_h->forward = fftw_plan_dft_r2c_2d(N,N, c_h->rval, c_h->cval, FFTW_PATIENT);
    c_h->backward = fftw_plan_dft_c2r_2d(N,N, c_h->cval, c_h->rval, FFTW_PATIENT);

    double *k_x = (double*) malloc(N*(N/2+1)*sizeof(double));
    double *k_y = (double*) malloc(N*(N/2+1)*sizeof(double));
    for (int i=0; i<N; i++) {
        for (int j=0; j<(N/2+1); j++) {
            k_y[ind_fft(i,j)] = (double)j;
            if (i<N/2)
                k_x[ind_fft(i,j)] = (double)i;
            else
                k_x[ind_fft(i,j)] = (double)(i-N);
        }
    }

    double pi2_squared = (2.0*M_PI)*(2.0*M_PI)*(-1.0);
    c_h->laplace = (double*) fftw_malloc(N*(N/2+1)*sizeof(double));
    for (int kk=0; kk<N*(N/2+1); kk++) {
        c_h->laplace[kk] = pi2_squared * (k_x[kk]*k_x[kk] + k_y[kk]*k_y[kk]);
    }

    return c_h;
}

void cahn_hilliard_free(cahn_hilliard* c_h)
{
    free(c_h->laplace);

    fftw_destroy_plan(c_h->forward);
    fftw_destroy_plan(c_h->backward);
    fftw_free(c_h->rval);
    fftw_free(c_h->cval);

    fftw_cleanup();
}
