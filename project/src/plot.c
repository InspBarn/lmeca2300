/* plot.c */

#include "plot.h"

bov_points_t* bov_points_new_with_value(GLfloat data[][3], GLsizei N, GLenum usage)
{
    GLfloat coord[N][2]; GLfloat value=0.0;
    for (int i=0; i<N; i++) {
        coord[i][0] = data[i][0];
        coord[i][1] = data[i][1];
        value += data[i][2]/(GLfloat)N;
    }
    bov_points_t *points = bov_points_new(coord,N,usage);

    GLfloat color[4];
    if (value>0.0) {
        color[0] = value;
        color[1] = 1.0-value;
        color[2] = 0.0;
    } else {
        color[0] = 0.0;
        color[1] = 1.0-value;
        color[2] = value;
    } color[3] = 1.0;
    bov_points_set_color(points, color);
    bov_points_set_width(points, 0);
    bov_points_set_outline_color(points, color);
    bov_points_set_outline_width(points, 0);

    return points;
}

void imshow(bov_window_t *window, double complex *z, int n1, int n2)
{
    float color[4] = {0.0,0.0,0.0,1.0};
    GLfloat data[5][3];

    int ind;
    for (float i=0.0; i<n1-1; i++) {
        for (float j=0.0; j<n2-1; j++) {
            ind = i*n2+j;

            data[0][0] = i/(double)n1;
            data[0][1] = j/(double)n2;
            data[0][2] = (float) creal(z[ind]);

            data[1][0] = (i+1)/(double)n1;
            data[1][1] = j/(double)n2;
            data[1][2] = (float) creal(z[ind+n2]);

            data[2][0] = (i+1)/(double)n1;
            data[2][1] = (j+1)/(double)n2;
            data[2][2] = (float) creal(z[ind+n2+1]);

            data[3][0] = i/(double)n1;
            data[3][1] = (j+1)/(double)n2;
            data[3][2] = (float) creal(z[ind+1]);

            data[4][0] = i/(double)n1;
            data[4][1] = j/(double)n2;
            data[4][2] = (float) creal(z[ind]);

            bov_points_t *points = bov_points_new_with_value(data, 5, GL_STREAM_DRAW);
            // bov_points_set_color(points, none);
            // bov_points_set_width(points, 0);
            // bov_points_set_outline_color(points, none);
            // bov_points_set_outline_width(points, 0);
            bov_triangles_draw(window, points, 0, 3);
            bov_triangles_draw(window, points, 2, 3);

            bov_points_delete(points);
        }
    }
}
