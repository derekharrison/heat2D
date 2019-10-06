/*
 * main.c
 *
 *  Created on: Oct 22, 2014
 *      Author: derek
 *
 *      This code solves the 2D transient heat conduction equation, gamma*div(grad(T))+q = rho*Cp*dT/dt,
 *      using the incomplete cholesky factorization conjugate gradient (ICCG) method. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../inc/user_types.h"
#include "../inc/heat2D.h"
#include "../inc/memory_functions.h"

double source_equation(double x,double y, double t);
double Tnfunc(double x, double t);
double Tsfunc(double x,double t);
double Twfunc(double y,double t);
double Tefunc(double y,double t);

/*-----------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[])
/*
 * Main function. Enter simulation parameters here.
 */
{
	clock_t begin, end;
	double time_spent;
	grid_parameters_t grid_parameters = {0};
	time_parameters_t time_parameters = {0};
	physical_params_t physical_params = {0};
	boundary_temperatures_t boundary_temperatures = {0};
	solver_results_t solver_results = {0};

	/* Adjustable parameters */
	grid_parameters.Lx = 4.0;           //Length of domain (m) along x coordinate
	grid_parameters.Ly = 7.0;			//length of domain (m) along y coordinate
	grid_parameters.nx = 8;             //Amount of nodes along x coordinate
	grid_parameters.ny = 8;				//amount of nodes along y coordinate

	time_parameters.maxts = 100;		//amount of timesteps
	time_parameters.to = 0.0;			//initial time (s). Start of simulation
	time_parameters.tf = 2.0;			//final time (s). End of simulation

	physical_params.gamma = 2.0;        //Heat conductivity (W/mK)
	physical_params.rho = 1.0;			//material density (kg/m3)
	physical_params.Cp = 2.0;			//heat capacity (J/kgK)
	/* End adjustable parameters */

	/* Setting boundary functions */
	boundary_temperatures.Tnfunc = Tnfunc;
	boundary_temperatures.Tsfunc = Tsfunc;
	boundary_temperatures.Tefunc = Tefunc;
	boundary_temperatures.Twfunc = Twfunc;

	/* Allocate result data */
	solver_results.T = matrix2D( grid_parameters.nx, grid_parameters.ny);
	solver_results.X = matrix2D( grid_parameters.nx, grid_parameters.ny);
	solver_results.Y = matrix2D( grid_parameters.nx, grid_parameters.ny);

	/* Executing solver */
	begin = clock();

	heat2D(grid_parameters,
		   time_parameters,
		   physical_params,
		   boundary_temperatures,
		   source_equation,
           &solver_results);

	end = clock();

	/* Print some results */
	int i, j;
	for (j = 1; j <= grid_parameters.ny; j++)
	{
	    for (i = 1; i <= grid_parameters.nx; i++)
	    {
	    	printf("T[%d][%d]: %f\t", i, j, solver_results.T[i][j]);
	    }
	    printf("\n");
	}

	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("\ntime spent: %f\n", time_spent);

	/* Deallocate result data */
	free_matrix2D(solver_results.T, grid_parameters.nx);
	free_matrix2D(solver_results.X, grid_parameters.nx);
	free_matrix2D(solver_results.Y, grid_parameters.nx);

	return 0;

}

/*-----------------------------------------------------------------------------------------------*/
double source_equation(double x, double y, double t)
/*
 * Source term in the equation:
 * gamma*div(grad(T))+q = rho*Cp*dT/dt
 *
 * input    x
 * input    y
 * input    t
 * return    Sq
 */
{
	double Sq = 10;

	return Sq;										//Enter equation for the source term here

}

/*-----------------------------------------------------------------------------------------------*/
double Tnfunc(double x, double t)
/*
 * North boundary temperature.
 *
 * input    x
 * input    t
 * return    T
 */
{
    double T = 1.0;

	return T + sin(2*t);								//Enter the equation for the north boundary condition here
}

/*-----------------------------------------------------------------------------------------------*/
double Tsfunc(double x, double t)
/*
 * South boundary temperature.
 *
 * input    x
 * input    t
 * return    T
 */
{
	double T = 1.0;

	return T + sin(2*t);								//Enter the equation for the south boundary condition here
}

/*-----------------------------------------------------------------------------------------------*/
double Twfunc(double y, double t)
/*
 * West boundary temperature.
 *
 * input    y
 * input    t
 */
{
	double T = 1.0;

	return T + sin(2*t);								//Enter the equation for the west boundary condition here
}

/*-----------------------------------------------------------------------------------------------*/
double Tefunc(double y, double t)
/*
 * East boundary temperature.
 *
 * input    y
 * input    t
 */
{
	double T = 4.0;

	return T + sin(2*t);								//Enter the equation for the east boundary condition here.
}


