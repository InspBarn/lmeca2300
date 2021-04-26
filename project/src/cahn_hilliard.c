/* cahn_hilliard.c */
#include "cahn_hilliard.h"

static const char* vertex_shader_text =
"#version 110\n"
"attribute vec3 vCol;\n"
"attribute vec2 vPos;\n"
"varying vec3 color;\n"
"void main()\n"
"{\n"
"    gl_Position = vec4(vPos, 0.0, 1.0);\n"
"    color = vCol;\n"
"}\n";

static const char* fragment_shader_text =
"#version 110\n"
"varying vec3 color;\n"
"void main()\n"
"{\n"
"    gl_FragColor = vec4(color, 1.0);\n"
"}\n";

static void normalize(double *u, const int N)
{
    for (int n=0; n<N*N; n++)
        u[n] /= (double)(N*N);
}

void f(cahn_hilliard *c_h, fftw_complex *u_hat, fftw_complex *u_dot_hat, double *k, double a)
{
    const int N = c_h->N;

    memcpy(c_h->cval,u_hat, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);

    for (int n=0; n<N*N; n++)
        c_h->rval[n] = c_h->rval[n]*c_h->rval[n]*c_h->rval[n];
    fftw_execute(c_h->forward);

    for (int i=0; i<N*(N/2+1); i++)
        u_dot_hat[i] = -k[i]*k[i] * (c_h->cval[i] - u_hat[i]) - (a*k[i]*k[i])*(a*k[i]*k[i])*u_hat[i];
}

void df(cahn_hilliard *c_h, fftw_complex *u_hat, fftw_complex *du_dot_hat, double *k, double a)
{
    const int N = c_h->N;

    memcpy(c_h->cval,u_hat, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);

    for (int n=0; n<N*N; n++)
        c_h->rval[n] = c_h->rval[n]*c_h->rval[n];
    fftw_execute(c_h->forward);

    for (int i=0; i<N*(N/2+1); i++)
        du_dot_hat[i] = -k[i]*k[i]*(c_h->cval[i]/(2.0*M_PI)-1.0) - (a*k[i]*k[i])*(a*k[i]*k[i]);
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
    f(c_h,temp,k1,c_h->k,c_h->a);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt/2.0 * k1[k];
    f(c_h,u_hat,k2,c_h->k,c_h->a);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt/2.0 * k2[k];
    f(c_h,u_hat,k3,c_h->k,c_h->a);

    for (int k=0; k<N*(N/2+1); k++)
        u_hat[k] = temp[k] + dt * k3[k];
    f(c_h,u_hat,k4,c_h->k,c_h->a);

    for (int k=0; k<N*(N/2+1); k++) {
        c_h->cval[k] = temp[k] + (k1[k] + 2.0*k2[k] + 2.0*k3[k] + k4[k]) * dt/6.0;
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
        f(c_h,c_hat,c_dot_hat,c_h->k,c_h->a);
        df(c_h,c_hat,dc_dot_hat,c_h->k,c_h->a);
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

void bdf_ab(cahn_hilliard *c_h, fftw_complex *c_prev, fftw_complex *c_prev_cub)
{
    const int N = c_h->N;
    const double dt = c_h->dt;

    fftw_complex *temp = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));

    fftw_execute(c_h->backward); normalize(c_h->rval,N);
    for (int n=0; n<N*N; n++)
        c_h->rval[n] = c_h->rval[n]*c_h->rval[n]*c_h->rval[n];
    fftw_execute(c_h->forward);

    fftw_complex *c_cub = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    memcpy(c_cub,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));

    double a2k4,k2; fftw_complex f_curr, f_prev;
    for (int i=0; i<N*(N/2+1); i++) {
        k2 = c_h->k[i]*c_h->k[i];
        a2k4 = (c_h->a*c_h->a) * (k2*k2);
        f_curr = c_cub[i]-temp[i];
        f_prev = c_prev_cub[i]-c_prev[i];
        if (c_h->t==0.0)
            c_h->cval[i] = (temp[i] - dt*k2*f_curr) / (1+dt*a2k4);
        else
            c_h->cval[i] = (4*temp[i] - c_prev[i] - 2.0*dt*k2*(2.0*f_curr - f_prev)) / (3.0+2.0*dt*a2k4);
    }

    memcpy(c_prev,temp, N*(N/2+1)*sizeof(fftw_complex));
    memcpy(c_prev_cub,c_cub, N*(N/2+1)*sizeof(fftw_complex));

    memcpy(temp,c_h->cval, N*(N/2+1)*sizeof(fftw_complex));
    fftw_execute(c_h->backward); normalize(c_h->rval,N);
    memcpy(c_h->cval,temp, N*(N/2+1)*sizeof(fftw_complex));

    fftw_free(temp); fftw_free(c_cub);
}

void cahn_hilliard_solve(cahn_hilliard *c_h, double *u)
{
    // c_h->dt /= 4.0;
    const int N = c_h->N, Ntime = (int) (c_h->tMax/c_h->dt);
    int t=0;

    memcpy(c_h->rval,u, N*N*sizeof(double));
    fftw_execute(c_h->forward);

    vertex_struct vertices[N*N];
    // int *indices = (int*)malloc(6*(N-1)*(N-1)*sizeof(int));
    int indices[6*(N-1)*(N-1)];

    float *color = malloc(3*sizeof(float));
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            vertices[i*N+j].x = 2.0 * ((float) j / (float) (N-1)) - 1.0;
            vertices[i*N+j].y = 2.0 * ((float) i / (float) (N-1)) - 1.0;
            vertices[i*N+j].r = 0.0;
            vertices[i*N+j].g = 0.0;
            vertices[i*N+j].b = 0.0;
        }
    }

    for (int i=0; i<(N-1); i++) {
        for (int j=0; j<(N-1); j++) {
            indices[6*(i*(N-1)+j)+0] = i*N+j;
            indices[6*(i*(N-1)+j)+1] = i*N+(j+1);
            indices[6*(i*(N-1)+j)+2] = (i+1)*N+(j+1);

            indices[6*(i*(N-1)+j)+3] = i*N+j;
            indices[6*(i*(N-1)+j)+4] = (i+1)*N+j;
            indices[6*(i*(N-1)+j)+5] = (i+1)*N+(j+1);
        }
    }

    GLFWwindow* window;
    GLuint vertex_array, element_buffer, vertex_buffer;
    GLuint vertex_shader, fragment_shader, program;
    GLint mvp_location, vpos_location, vcol_location;

    window = createWindow(800, 800, "Simple example");

    // NOTE: OpenGL error checks have been omitted for brevity
    glGenVertexArrays(1, &vertex_array);

    glGenBuffers(1, &vertex_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glGenBuffers(1, &element_buffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &vertex_shader_text, NULL);
    glCompileShader(vertex_shader);

    fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragment_shader_text, NULL);
    glCompileShader(fragment_shader);

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    // mvp_location = glGetUniformLocation(program, "MVP");
    vpos_location = glGetAttribLocation(program, "vPos");
    vcol_location = glGetAttribLocation(program, "vCol");

    glEnableVertexAttribArray(vpos_location);
    glVertexAttribPointer(vpos_location, 2, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (void*) 0);
    glEnableVertexAttribArray(vcol_location);
    glVertexAttribPointer(vcol_location, 3, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (void*) (sizeof(float)*2));

    clock_t start,end;
    fftw_complex *c_prev = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    fftw_complex *c_prev_cub = (fftw_complex*) fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
    do {
        printf("iteration : %5d / %d -- ", t,Ntime);

        start = clock();
        // glfwGetFramebufferSize(window, &width, &height);
        // ratio = width / (float) height;
        // glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT);

        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                jet(color, c_h->rval[i*N+j]);
                vertices[i*N+j].r = color[0];
                vertices[i*N+j].g = color[1];
                vertices[i*N+j].b = color[2];
            }
        }
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

        glUseProgram(program);
        glDrawElements(GL_TRIANGLES, 6*(N-1)*(N-1), GL_UNSIGNED_INT, 0);

        updateWindow(window);
        end = clock();
        printf("time for drawing : %.3f [ms] -- ", (double)(end-start)/CLOCKS_PER_SEC*1e3);

        start = clock();
        for (int i=0; i<5; i++) {
            c_h->t += c_h->dt;
            bdf_ab(c_h, c_prev, c_prev_cub);
            // euler_implicit(c_h);
            // rk4(c_h);
            t++;
        } end = clock();

        printf("time for integration : %.3f [ms] -- ", (double)(end-start)/CLOCKS_PER_SEC*1e3);
        printf("time : %.6f\n", c_h->t);
        if (t>Ntime)
            break;

    } while (!glfwWindowShouldClose(window));

    destroyWindow(window);
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

    double pi2 = (2.0*M_PI);
    c_h->k = (double*) fftw_malloc(N*(N/2+1)*sizeof(double));
    for (int kk=0; kk<N*(N/2+1); kk++) {
        c_h->k[kk] = pi2 * sqrt(k_x[kk]*k_x[kk] + k_y[kk]*k_y[kk]);
    }

    return c_h;
}

void cahn_hilliard_free(cahn_hilliard* c_h)
{
    free(c_h->k);

    fftw_destroy_plan(c_h->forward);
    fftw_destroy_plan(c_h->backward);
    fftw_free(c_h->rval);
    fftw_free(c_h->cval);

    fftw_cleanup();
}
