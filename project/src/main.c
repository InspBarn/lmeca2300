/* main.c */


// #include "utils.h"
// 
// #include <stdlib.h>
// #include <stdio.h>
// 
// static const struct
// {
//     float x, y;
//     float r, g, b;
// } vertices[4] =
// {
//     {   0.5f,  0.5f, 0.f, 0.f, 1.f },
//     {   0.5f, -0.5f, 1.f, 0.f, 0.f },
//     {  -0.5f, -0.5f, 0.f, 1.f, 0.f },
//     {  -0.5f,  0.5f, 0.f, 0.f, 1.f }
// };
// 
// unsigned int indices[] = {
//     0, 1, 2,
//     1, 2, 3
// };
// 
// static const char* vertex_shader_text =
// "#version 110\n"
// "uniform mat4 MVP;\n"
// "attribute vec3 vCol;\n"
// "attribute vec2 vPos;\n"
// "varying vec3 color;\n"
// "void main()\n"
// "{\n"
// "    gl_Position = MVP * vec4(vPos, 0.0, 1.0);\n"
// "    color = vCol;\n"
// "}\n";
// 
// static const char* fragment_shader_text =
// "#version 110\n"
// "varying vec3 color;\n"
// "void main()\n"
// "{\n"
// "    gl_FragColor = vec4(color, 1.0);\n"
// "}\n";
// 
// int main(void)
// {
//     GLFWwindow* window;
//     GLuint vertex_buffer, order_buffer, vertex_shader, fragment_shader, program;
//     GLint mvp_location, vpos_location, vcol_location;
// 
//     window = createWindow(640, 480, "Simple example");
// 
//     // NOTE: OpenGL error checks have been omitted for brevity
// 
//     glGenBuffers(1, &vertex_buffer);
//     glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
//     glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
// 
//     glGenBuffers(1, &order_buffer);
//     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, order_buffer);
//     glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
// 
//     vertex_shader = glCreateShader(GL_VERTEX_SHADER);
//     glShaderSource(vertex_shader, 1, &vertex_shader_text, NULL);
//     glCompileShader(vertex_shader);
// 
//     fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
//     glShaderSource(fragment_shader, 1, &fragment_shader_text, NULL);
//     glCompileShader(fragment_shader);
// 
//     program = glCreateProgram();
//     glAttachShader(program, vertex_shader);
//     glAttachShader(program, fragment_shader);
//     glLinkProgram(program);
// 
//     // mvp_location = glGetUniformLocation(program, "MVP");
//     // vpos_location = glGetAttribLocation(program, "vPos");
//     // vcol_location = glGetAttribLocation(program, "vCol");
// 
//     glEnableVertexAttribArray(vpos_location);
//     glVertexAttribPointer(vpos_location, 2, GL_FLOAT, GL_FALSE,
//                           sizeof(vertices[0]), (void*) 0);
//     glEnableVertexAttribArray(vcol_location);
//     glVertexAttribPointer(vcol_location, 3, GL_FLOAT, GL_FALSE,
//                           sizeof(vertices[0]), (void*) (sizeof(float) * 2));
// 
//     do {
//         float ratio;
//         int width, height;
//         mat4x4 m, p, mvp;
// 
//         glfwGetFramebufferSize(window, &width, &height);
//         ratio = width / (float) height;
// 
//         glViewport(0, 0, width, height);
//         glClear(GL_COLOR_BUFFER_BIT);
// 
//         mat4x4_identity(m);
//         mat4x4_rotate_Z(m, m, (float) glfwGetTime());
//         mat4x4_ortho(p, -ratio, ratio, -1.f, 1.f, 1.f, -1.f);
//         mat4x4_mul(mvp, p, m);
// 
//         glUseProgram(program);
//         glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (const GLfloat*) mvp);
//         glDrawArrays(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
// 
//         updateWindow(window);
//     } while (!glfwWindowShouldClose(window));
// 
//     destroyWindow(window);
//     exit(EXIT_SUCCESS);
// }
 

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

    exit(EXIT_SUCCESS);
}
