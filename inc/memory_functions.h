/*
 * memory_functions.h
 *
 *  Created on: Sep 18, 2019
 *      Author: d-w-h
 */

#ifndef MEMORY_FUNCTIONS_H_
#define MEMORY_FUNCTIONS_H_

#include "user_types.h"

double *matrix1D(int n);
void free_matrix1D(double *vec_ptr);
double **matrix2D( int nx, int ny);
void free_matrix2D(double** mat_ptr, int nx);
void allocate_solver_data_mem(solver_data_t* solver_data,
                              grid_parameters_t grid_parameters,
                              grid_coordinates_t* grid_coordinates);
void deallocate_solver_data_mem(solver_data_t* solver_data,
                                grid_parameters_t grid_parameters,
                                grid_coordinates_t* grid_coordinates);

#endif /* MEMORY_FUNCTIONS_H_ */
