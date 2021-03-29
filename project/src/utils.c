/* utils.c */

#include "utils.h"

void cubic(double *u, const int N)
{
	for (int n=0; n<N; n++) {
		// printf("%.3f\t", u[n]);
		u[n] = u[n]*u[n]*u[n];
	} 
}
