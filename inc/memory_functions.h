/*
 * memory_functions.h
 *
 *  Created on: Sep 18, 2019
 *      Author: d-w-h
 */

#ifndef MEMORY_FUNCTIONS_H_
#define MEMORY_FUNCTIONS_H_

#include "user_types.h"

double*** result_vector();
void free_result_vector(double*** results);
double *matrix1D(int np);
void free_matrix1D(double *a);
double **matrix2D( int nm, int np);
void free_matrix2D(double** m, int nm);
double ***matrix3D( int nx, int ny, int nz);

void allocate_solver_data_mem(solver_data_t* solver_data,
                              grid_parameters_t grid_parameters,
                              grid_coordinates_t* grid_coordinates);
void deallocate_solver_data_mem(solver_data_t* solver_data,
                                grid_parameters_t grid_parameters);

#endif /* MEMORY_FUNCTIONS_H_ */
