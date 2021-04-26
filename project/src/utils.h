/* utils.h */

// #include <glad/glad.h>
// #define GLFW_INCLUDE_NONE
// #include <GLFW/glfw3.h>
//  
// #include "linmath.h"
#include "BOV.h"
 
#include <stdlib.h>
#include <stdio.h>

void binary(float *color, float value);
void hot_to_cold(float *color, float value);
void jet(float *color, float value);

void imshow(bov_window_t *window, double *z, int n1, int n2);

// GLFWwindow* createWindow(int width, int height, const char *title);
// void destroyWindow(GLFWwindow *window);
// void updateWindow(GLFWwindow *window);
