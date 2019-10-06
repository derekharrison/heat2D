/*
 * heat2D.c
 *
 *  Created on: Sep 18, 2019
 *      Author: d-w-h
 */

#include <stdio.h>
#include <stdlib.h>

#include "../inc/heat2D_support.h"
#include "../inc/memory_functions.h"
#include "../inc/user_types.h"
#include "../inc/vector_operations.h"

/*-----------------------------------------------------------------------------------------------*/
void heat2D(grid_parameters_t grid_parameters,
			time_parameters_t time_parameters,
			physical_params_t physical_params,
			boundary_temperatures_t boundary_temperatures,
			source_ptr source_equation,
			solver_results_t* solver_results)
/*
 * This function solves the 2D transient heat conduction equation,
 * gamma*div(grad(T))+q = rho*Cp*dT/dt
 * using the incomplete cholesky factorization conjugate gradient (ICCG) method
 *
 * input    grid_parameters
 * input    time_parameters
 * input    physical_params
 * input    boundary_temperatures
 * input    source_equation
 * return   results
 */
{
	int nx, ny, maxts;
	int imax, nt, it, countert;
	double error, dt, epsilon, to;
	grid_coordinates_t grid_coordinates = {0};
	solver_data_t solver_data = {0};

	nx = grid_parameters.nx;
	ny = grid_parameters.ny;

	to = time_parameters.to;
	maxts = time_parameters.maxts;

	dt = (time_parameters.tf - time_parameters.to) / time_parameters.maxts;

	nt = nx * ny;
	imax = 500;         //Maximum iterations ICCG
	error = 1e-30;      //Tolerance

	/* Allocate solver and result data */
	allocate_solver_data_mem(&solver_data,
			                 grid_parameters,
			                 &grid_coordinates);

	set_initial_temp(&solver_data,
			         grid_parameters);

	generate_grid_coordinates(grid_parameters,
			                  &grid_coordinates);

	/* Entering ICCG loop */
	countert = 0;
	time_parameters.t = to;
	do{
		calc_coefficient_matrix(grid_parameters,
								grid_coordinates,
								time_parameters,
								physical_params,
								&solver_data,
								source_equation,
								boundary_temperatures);

		store_coefficient_matrix(&solver_data,
				                 grid_parameters);

		incomplete_cholesky_factorization(&solver_data,
		                                  grid_parameters);

		epsilon = dot_product(solver_data.r,
				              solver_data.r,
				              nt);

		solve_Mz_is_r(&solver_data,
				      grid_parameters);

		set_solver_data(&solver_data,
				        grid_parameters);

		/*Perform solver iterations*/
		it = 0;
		do
		{
			perform_cgm(&solver_data,
					    grid_parameters,
					    &epsilon);

			it = it + 1;
		}while (it < imax && epsilon > error);

		set_current_temperature(&solver_data,
				                grid_parameters);

		++countert;
		time_parameters.t += dt;
	}while(countert < maxts);

	/* Processing results */
	set_temperature_result_data(&solver_data,
			                    grid_parameters);

	/* Deallocate solver data */
	deallocate_solver_data_mem(&solver_data, grid_parameters);

	/* Setting results */
	solver_results->grid_coordinates.X = grid_coordinates.X;
	solver_results->grid_coordinates.Y = grid_coordinates.Y;
	solver_results->T = solver_data.T;

}
