/* utils.c */

#include "utils.h"

void binary(GLfloat *color, GLfloat value)
{
    color[0] = (value+1.0)/2.0;
    color[1] = (value+1.0)/2.0;
    color[2] = (value+1.0)/2.0;
    color[3] = 1.0;
}

void hot_to_cold(GLfloat *color, GLfloat value)
{
    if (value<-0.5) {
        color[0] = 0.0;
        color[1] = (value+1.0)*2.0;
        color[2] = 1.0;
    } else if (value<0.0) {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 1.0 - (value+0.5)*2.0;
    } else if (value<0.5) {
        color[0] = value*2.0;
        color[1] = 1.0;
        color[2] = 0.0;
    } else {
        color[0] = 1.0;
        color[1] = 1.0 - (value-0.5)*2.0;
        color[2] = 0.0;
    } color[3] = 1.0;
}

void jet(GLfloat *color, GLfloat value)
{
    if (value<-0.75) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = (value+1.25)*2.0;
    } else if (value<-0.25) {
        color[0] = 0.0;
        color[1] = (value+0.75)*2.0;
        color[2] = 1.0;
    } else if (value<0.25) {
        color[0] = (value+0.25)*2.0;
        color[1] = 1.0;
        color[2] = 1.0 - (value+0.25)*2.0;
    } else if (value<0.75) {
        color[0] = 1.0;
        color[1] = 1.0 - (value-0.25)*2.0;
        color[2] = 0.0;
    } else {
        color[0] = 1.0 - (value-0.75)*2.0;
        color[1] = 0.0;
        color[2] = 0.0;
    } color[3] = 1.0;
}

static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

GLFWwindow* createWindow(int width, int height, const char *title)
{
	GLFWwindow *window;

    /* Initialize the library */
    if (!glfwInit())
        exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    /* Setting an error callback */
    glfwSetErrorCallback(error_callback);

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(width,height,title, NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    /* Set callbacks */
    glfwSetKeyCallback(window, key_callback);

    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    gladLoadGL();
    glfwSwapInterval(1);

    return window;
}

void destroyWindow(GLFWwindow *window)
{
    glfwDestroyWindow(window);
    glfwTerminate();
}

void updateWindow(GLFWwindow *window)
{
    glfwSwapBuffers(window);
    glfwPollEvents();
}

// bov_points_t* bov_points_new_with_value(GLfloat data[][3], GLsizei N, GLenum usage)
// {
//     GLfloat coord[N][2]; GLfloat value=0.0;
//     for (int i=0; i<N; i++) {
//         coord[i][0] = data[i][0];
//         coord[i][1] = data[i][1];
//         value += data[i][2]/(GLfloat)N;
//     }
//     bov_points_t *points = bov_points_new(coord,N,usage);
// 
//     GLfloat *color = (GLfloat*) malloc(4*sizeof(GLfloat));
//     binary(color,value);
// 
//     bov_points_set_color(points, color);
//     bov_points_set_width(points, 0);
//     bov_points_set_outline_color(points, color);
//     bov_points_set_outline_width(points, 0);
//     return points;
// }
// 
// void imshow(bov_window_t *window, double *z, int n1, int n2)
// {
//     GLfloat data[5][3];
// 
//     int ind;
//     for (float i=0.0; i<n1-1; i++) {
//         for (float j=0.0; j<n2-1; j++) {
//             ind = i*n2+j;
// 
//             data[0][0] = i;
//             data[0][1] = j;
//             data[0][2] = (float) z[ind];
// 
//             data[1][0] = (i+1);
//             data[1][1] = j;
//             data[1][2] = (float) z[ind+1];
// 
//             data[2][0] = (i+1);
//             data[2][1] = (j+1);
//             data[2][2] = (float) z[ind+n2+1];
// 
//             data[3][0] = i;
//             data[3][1] = (j+1);
//             data[3][2] = (float) z[ind+n2];
// 
//             data[4][0] = i;
//             data[4][1] = j;
//             data[4][2] = (float) z[ind];
// 
//             bov_points_t *points = bov_points_new_with_value(data, 5, GL_STREAM_DRAW);
//             bov_triangles_draw(window, points, 0, 3);
//             bov_triangles_draw(window, points, 2, 3);
//             bov_points_delete(points);
//         }
//     }
// }
