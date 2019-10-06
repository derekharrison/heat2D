/*
 * heat2D_support.h
 *
 *  Created on: Sep 18, 2019
 *      Author: d-w-h
 */

#ifndef HEAT2D_SUPPORT_H_
#define HEAT2D_SUPPORT_H_

#include "../inc/user_types.h"

void set_initial_temp(solver_data_t* solver_data,
                      grid_parameters_t grid_parameters);
void generate_grid_coordinates(grid_parameters_t grid_parameters,
                               grid_coordinates_t* grid_coordinates);
void calc_coefficient_matrix(grid_parameters_t grid_parameters,
                             grid_coordinates_t grid_coordinates,
                             time_parameters_t time_parameters,
                             physical_params_t physical_params,
                             solver_data_t* solver_data,
                             source_ptr q,
                             boundary_temperatures_t boundary_temperatures);
void store_coefficient_matrix(solver_data_t* solver_data,
                              grid_parameters_t grid_parameters);
void incomplete_cholesky_factorization(solver_data_t* solver_data,
                                       grid_parameters_t grid_parameters);
void solve_Mz_is_r(solver_data_t* solver_data,
                   grid_parameters_t grid_parameters);
void perform_cgm(solver_data_t* solver_data,
                 grid_parameters_t grid_parameters,
                 double* epsilon);
void set_solver_data(solver_data_t* solver_data,
                     grid_parameters_t grid_parameters);
void set_current_temperature(solver_data_t* solver_data,
                             grid_parameters_t grid_parameters);
void set_result_data(solver_data_t* solver_data,
                     grid_coordinates_t* grid_coordinates,
                     grid_parameters_t grid_parameters,
                     solver_results_t* solver_results);
#endif /* HEAT2D_SUPPORT_H_ */
