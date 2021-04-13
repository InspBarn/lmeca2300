/* cahn_hilliard.c */

#include "cahn_hilliard.h"

void normalize(double *u, const int N)
{
    for (int n=0; n<N*N; n++)
        u[n] /= (double)(N*N);
}

void slice_derivatives(cahn_hilliard *c_h)
{
    const int N = c_h->N;
    for (int k=0; k<N*(N/2+1); k++) {
        // c_h->cval_dot[3][k] = c_h->cval_dot[2][k];
        c_h->cval_dot[2][k] = c_h->cval_dot[1][k];
        c_h->cval_dot[1][k] = c_h->cval_dot[0][k];
    }
}

void f(cahn_hilliard *c_h, fftw_complex *u_hat, fftw_complex *u_dot_hat)
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

void df(cahn_hilliard *c_h, fftw_complex *u_hat, fftw_complex *du_dot_hat)
{
    const int N = c_h->N;

    memcpy(c_h->cval,u_hat, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);

    for (int n=0; n<N*N; n++)
        c_h->rval[n] = c_h->rval[n]*c_h->rval[n];
    fftw_execute(c_h->forward);

    double a_squared = c_h->a*c_h->a;
    for (int k=0; k<N*(N/2+1); k++)
        du_dot_hat[k] = c_h->laplace[k] * (c_h->cval[k]/(2.0*M_PI) - 1.0 - a_squared*c_h->laplace[k]);
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
    f(c_h,temp,k1);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt/2.0 * k1[k];
    f(c_h,u_hat,k2);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt/2.0 * k2[k];
    f(c_h,u_hat,k3);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt * k3[k];
    f(c_h,u_hat,k4);

    // slice_derivatives(c_h);
    for (int k=0; k<N*(N/2+1); k++) {
        c_h->cval_dot[0][k] = (k1[k] + 2.0*k2[k] + 2.0*k3[k] + k4[k]);
        c_h->cval[k] = temp[k] + c_h->cval_dot[0][k] * dt/6.0;
    }

    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);
    memcpy(c_h->cval,temp, N*(N/2+1)*sizeof(fftw_complex));

    fftw_free(k1); fftw_free(k2);
    fftw_free(k3); fftw_free(k4);
    fftw_free(temp); fftw_free(u_hat);
}

void euler_implicit(cahn_hilliard *c_h)
{
    const int N = c_h->N;
    const double dt = c_h->dt;

    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));

    int iter = 0;
    double error = 1.0, err;
    fftw_complex *c_hat = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *c_dot_hat = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *dc_dot_hat = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex dc_hat; memcpy(c_hat,temp, N*(N/2+1)*sizeof(fftw_complex));
    while ((error>1e-3) && (iter<100)) {
        error = 0.0;
        f(c_h,c_hat,c_dot_hat);
        df(c_h,c_hat,dc_dot_hat);
        for (int k=0; k<N*(N/2+1); k++) {
            dc_hat = (temp[k]+dt*c_dot_hat[k]-c_hat[k]) / (1-dt*dc_dot_hat[k]);
            c_hat[k] = c_hat[k] + dc_hat;
            err = pow(creal(dc_hat)*creal(dc_hat) + cimag(dc_hat)*cimag(dc_hat),0.5);
            if (err>error)
                error = err;
        } iter ++;
    }

    memcpy(c_h->cval,c_hat, N*(N/2+1)*sizeof(fftw_complex));
    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);
    memcpy(c_h->cval,temp, N*(N/2+1)*sizeof(fftw_complex));

    fftw_free(temp);
    fftw_free(c_hat);
    fftw_free(c_dot_hat);
    fftw_free(dc_dot_hat);
}

void ab4_am4(cahn_hilliard *c_h)
{
    const int N = c_h->N;
    const double dt = c_h->dt;

    fftw_complex *predictor = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *corrector = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *u_st = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));

    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));

    f(c_h,temp,predictor);
    for (int k=0; k<N*(N/2+1); k++)
        u_st[k] = temp[k] + dt/24.0*(55.0*predictor[k] - 59.0*c_h->cval_dot[0][k] \
                            + 37.0*c_h->cval_dot[1][k] - 9.0*c_h->cval_dot[2][k]);
    f(c_h,u_st,corrector);

    slice_derivatives(c_h);
    for (int k=0; k<N*(N/2+1); k++) {
        c_h->cval_dot[0][k] = predictor[k];
        c_h->cval[k] = temp[k] + dt/24.0*(9.0*corrector[k] + 19.0*c_h->cval_dot[0][k] \
                            - 5.0*c_h->cval_dot[1][k] + 1.0*c_h->cval_dot[2][k]);
    }

    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);
    memcpy(c_h->cval,temp, N*(N/2+1)*sizeof(fftw_complex));

    fftw_free(predictor);
    fftw_free(corrector);
    fftw_free(temp);
    fftw_free(u_st);
}

void cahn_hilliard_solve(cahn_hilliard *c_h, double *u)
{
    // c_h->dt /= 4.0;
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
        for (int i=0; i<100; i++) {
            c_h->t += c_h->dt;
            euler_implicit(c_h);
            // rk4(c_h);
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
    c_h->dt = 1e-7;
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

    c_h->cval_dot = fftw_malloc(N*(N/2+1)*sizeof(c_h->cval_dot[0]));

    return c_h;
}

void cahn_hilliard_free(cahn_hilliard* c_h)
{
    free(c_h->laplace);

    fftw_destroy_plan(c_h->forward);
    fftw_destroy_plan(c_h->backward);
    fftw_free(c_h->rval);
    fftw_free(c_h->cval);
    fftw_free(c_h->cval_dot);

    fftw_cleanup();
}
