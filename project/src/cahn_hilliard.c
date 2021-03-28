/* cahn_hilliard.c */

#include "cahn_hilliard.h"

void step_derivative(cahn_hilliard *c_h, fftw_complex *u_hat, fftw_complex *u_dot_hat)
{
    const int N = c_h->N;
    double *u = (double*) fftw_malloc(N*N*sizeof(double));
    fftw_complex *u_hat3 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));

    memcpy(temp,u_hat,N*(N/2+1)*sizeof(fftw_complex));

    fftw_execute_dft_c2r(c_h->backward,u_hat,u); cubic(u,N*N);
    fftw_execute_dft_r2c(c_h->forward,u,u_hat3);

    memcpy(u_hat,temp,N*(N/2+1)*sizeof(fftw_complex));

    double a_sqared = pow(c_h->a, 2.0);
    for (int k=0; k<N*(N/2+1); k++) {
        u_dot_hat[k] = c_h->laplace[k] * (u_hat3[k] - u_hat[k] - a_sqared*c_h->laplace[k]*u_hat[k])/(N*N);
        // u_dot_hat[k] = c_h->laplace[k] * (c_h->nonlin_out[k] - u_hat[k] + a_sqared*c_h->laplace[k]*u_hat[k])/(N*N);
    }

    free(u); free(u_hat3);
}

void rk4(cahn_hilliard *c_h)
{
    const int N = c_h->N;
    const double t = c_h->t;
    const double dt = c_h->dt;

    fftw_complex *k1 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *k2 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *k3 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *k4 = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *u_hat = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));

    // step_derivative(c_h,c_h->out,k1);
    memcpy(u_hat, c_h->out, N*(N/2+1)*sizeof(fftw_complex));
    step_derivative(c_h,u_hat,k1);

    for (int k=0; k<N*(N/2+1); k++)
        temp[k] = u_hat[k] + dt/2.0 * k1[k];
    step_derivative(c_h,temp,k2);

    for (int k=0; k<N*(N/2+1); k++)
        temp[k] = u_hat[k] + dt/2.0 * k2[k];
    step_derivative(c_h,temp,k3);

    for (int k=0; k<N*(N/2+1); k++)
        temp[k] = u_hat[k] + dt * k3[k];
    step_derivative(c_h,temp,k4);

    for (int k=0; k<N*(N/2+1); k++) {
        c_h->out[k] = u_hat[k] + (k1[k] + k2[k] + k3[k] + k4[k]) * dt/6.0;
        // printf("%f\n", creal(c_h->out[k]));
    } // printf("\n");

    free(k1); free(k2);
    free(k3); free(k4);
    free(temp); free(u_hat);
}

void cahn_hilliard_solve(cahn_hilliard *c_h)
{
    c_h->dt /= 4.0;
    int N = c_h->N, t=0, Ntime = (int) (c_h->tMax/c_h->dt);
    // printf("%d\n", Ntime);

    for (int n=0; n<N*N; n++)
        c_h->in[n] = (rand() / (double)RAND_MAX - 0.5)*2.0;

    fftw_execute(c_h->forward);
    for (int k=0; k<N*(N/2+1); k++)
        c_h->out[k] /= (N*N);

    for (int i=0; i<N; i++) {
        // printf("%.3f\n", c_h->in[i]);
        printf("%.3f + %.3fi\n", creal(c_h->out[i]),cimag(c_h->out[i]));
    } printf("\n");

    bov_window_t *window = bov_window_new(800, 800, "Cahn Hilliard Simulation");
    // bov_window_set_color(window, (GLfloat[]){0.9f,0.85f,0.8f,1.0f});
    window->param.zoom = 2.0/((double)N-1.0);
    window->param.translate[0] = -((double)N-1.0)/2.0;
    window->param.translate[1] = -((double)N-1.0)/2.0;

    clock_t start,end;
    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    do {
        bov_window_update(window);

        // memcpy(temp, c_h->out, N*(N/2+1)*sizeof(fftw_complex));
        fftw_execute(c_h->backward);
        // memcpy(c_h->out, temp, N*(N/2+1)*sizeof(fftw_complex));

        start = clock();
        imshow(window, c_h->in, N,N);
        end = clock();
        printf("time for drawing : %.3f -- ", (double)(end-start)/CLOCKS_PER_SEC);

        // fftw_execute(c_h->forward);

        for (int i=0; i<4; i++) {
            c_h->t += c_h->dt;
            rk4(c_h);
            t++;
        }

        printf("iter : %d -- time : %.6f\n", t,c_h->t);
        if (t>100)
            break;

    } while(!bov_window_should_close(window));

    for (int i=0; i<N; i++) {
        // printf("%.3f\n", c_h->in[i]);
        printf("%.3f + %.3fi\n", creal(c_h->out[i]),cimag(c_h->out[i]));
    } printf("\n");

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

    c_h->in = (double*) fftw_malloc(N*N*sizeof(double));
    c_h->out = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    c_h->forward = fftw_plan_dft_r2c_2d(N,N, c_h->in, c_h->out, FFTW_ESTIMATE);
    c_h->backward = fftw_plan_dft_c2r_2d(N,N, c_h->out, c_h->in, FFTW_ESTIMATE);

    c_h->nonlin_in = (double*) fftw_malloc(N*N*sizeof(double));
    c_h->nonlin_out = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    c_h->nonlin_forward = fftw_plan_dft_r2c_2d(N,N, c_h->nonlin_in, c_h->nonlin_out, FFTW_ESTIMATE);
    c_h->nonlin_backward = fftw_plan_dft_c2r_2d(N,N, c_h->nonlin_out, c_h->nonlin_in, FFTW_ESTIMATE);

    double pi2_squared = pow(2.0*I*M_PI,2.0);
    c_h->laplace = (double*) fftw_malloc(N*(N/2+1)*sizeof(double));
    for (int i=0; i<N; i++) {
        for (int j=0; j<(N/2+1); j++) {
            if (i<N/2)
                c_h->laplace[i+j*(N/2+1)] = pi2_squared * (pow(i,2.0) + pow(j,2.0));
            else
                c_h->laplace[i+j*(N/2+1)] = pi2_squared * (pow((i-N),2.0) + pow(j,2.0));
        }
    }
//    for (int i=0; i<(N/2+1); i++) {
//        for (int j=0; j<N; j++) {
//            if (j<N/2)
//                c_h->laplace[ind_fft(i,j)] = pi2_squared * (pow(i,2.0) + pow(j,2.0));
//            else
//                c_h->laplace[ind_fft(i,j)] = pi2_squared * (pow((i-N),2.0) + pow(j,2.0));
//        }
//    }

    return c_h;
}

void cahn_hilliard_free(cahn_hilliard* c_h)
{
    free(c_h->laplace);

    fftw_destroy_plan(c_h->forward);
    fftw_destroy_plan(c_h->backward);
    fftw_free(c_h->in);
    fftw_free(c_h->out);

    fftw_destroy_plan(c_h->nonlin_forward);
    fftw_destroy_plan(c_h->nonlin_backward);
    fftw_free(c_h->nonlin_in);
    fftw_free(c_h->nonlin_out);

    fftw_cleanup();
}
