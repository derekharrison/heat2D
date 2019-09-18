/*
 * memory_functions.c
 *
 *  Created on: Sep 18, 2019
 *      Author: d-w-h
 */
#include <stdio.h>
#include <stdlib.h>

#include "../inc/user_types.h"

/*-----------------------------------------------------------------------------------------------*/
double*** result_vector()
/*
 * Allocate memory for result vector
 *
 * return results
 */
{
	double*** results;
    results = malloc(sizeof(double**) * SIZE_RESULTS);

    return results;

}

/*-----------------------------------------------------------------------------------------------*/
double *matrix1D(int np)
/*
 * Allocate memory for vector of size np
 *
 * input np
 * return vec
 */
{
    double *vec;

    vec = (double *) calloc(np, sizeof(double));

    return vec;

}

/*-----------------------------------------------------------------------------------------------*/
void free_matrix1D(double *vec_ptr)
/*
 * Deallocate memory of vector
 *
 * input vec_ptr
 */
{
	free(vec_ptr);

}

/*-----------------------------------------------------------------------------------------------*/
double **matrix2D( int nm, int np)
/*
 * Allocate memory for matrix of size nm  by np
 *
 * input nm
 * input np
 * return mat
 */
{
   int i;
   double **mat;

   mat = (double **) calloc ( nm, sizeof( double *));
   for ( i = 0; i < nm; i++)
      mat[i] = (double *) calloc ( np, sizeof( double));

   return mat;

}

/*-----------------------------------------------------------------------------------------------*/
void free_matrix2D(double** m, int nm)
/*
 * Deallocate memory of matrix of size nm
 *
 * input m
 * input nm
 */
{
   int i;
   for ( i = 0; i < nm; i++)
      free(m[i]);

   free(m);

}

/*-----------------------------------------------------------------------------------------------*/
void allocate_solver_data_mem(solver_data_t* solver_data,
		                      grid_parameters_t grid_parameters,
		                      grid_coordinates_t* grid_coordinates)
/*
 * Allocate memory for solver data
 *
 * output   solver_data
 * input    grid_parameters
 * output   grid_coordinates
 */
{
	int nx, ny, nt;

	nx = grid_parameters.nx;
	ny = grid_parameters.ny;

	nt = nx * ny;

	solver_data->A = matrix2D(nt+1,3+1);
	solver_data->Astor = matrix2D(nt+1,3+1);
	solver_data->y = matrix1D(nt+1);
	solver_data->z = matrix1D(nt+1);
	solver_data->p = matrix1D(nt+1);
	solver_data->x = matrix1D(nt+1);
	solver_data->xo = matrix1D(nx*ny+1);
	solver_data->r = matrix1D(nt+1);
	solver_data->T = matrix2D(nx+1,ny+1);

	grid_coordinates->X = matrix2D(nx+1,ny+1);
	grid_coordinates->Y = matrix2D(nx+1,ny+1);

}

/*-----------------------------------------------------------------------------------------------*/
void deallocate_solver_data_mem(solver_data_t* solver_data,
		                        grid_parameters_t grid_parameters)
/*
 * Deallocate memory of solver data
 *
 * output   solver_data
 * input    grid_parameters
 */
{
	int nx, ny, nt;

	nx = grid_parameters.nx;
	ny = grid_parameters.ny;

	nt = nx * ny;

	free_matrix2D(solver_data->A, nt+1);
	free_matrix2D(solver_data->Astor, nt+1);
	free_matrix1D(solver_data->y);
	free_matrix1D(solver_data->z);
	free_matrix1D(solver_data->p);
	free_matrix1D(solver_data->x);
	free_matrix1D(solver_data->x);
	free_matrix1D(solver_data->r);

}
