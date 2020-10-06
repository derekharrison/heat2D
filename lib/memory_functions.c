/*
 * memory_functions.c
 *
 *  Created on: Oct 6, 2020
 *      Author: d-w-h
 */

#include <stdlib.h>

double* matrix1D(int np)
{
    double *a;

    a = (double *) calloc(np, sizeof(double));

    return a;
}

double** matrix2D( int nx, int ny)
{
    int i;
    double **m;

    m = (double **) calloc (nx, sizeof( double *));
    for ( i = 0; i < nx; i++)
        m[i] = (double *) calloc (ny, sizeof( double));

    return m;
}

void free2Df(double** m, int nx)
{
	int i;
    for (i = 0; i < nx; i++)
	    free(m[i]);

    free(m);
}
