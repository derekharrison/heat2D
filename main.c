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
#include "heat.h"
#include "memory_functions.h"
#include "parameters_and_boundaries.h"
#include "user_types.h"

int main(int argc, char *argv[])
{
    g_data grid_data;
    t_data time_data;
    p_params physical_params;
    s_data solver_data;
    boundaries_t boundaries;

    /*Adjustable parameters*/
    grid_data.Lx = 1.0;                       //Length of domain (m) along x coordinate
    grid_data.Ly = 1.0;                       //length of domain (m) along y coordinate
    grid_data.nx = 8;                         //Amount of nodes along x coordinate
    grid_data.ny = 8;                         //amount of nodes along y coordinate

    time_data.maxts = 20;                     //amount of timesteps
    time_data.to = 0.0;                       //initial time (s). Start of simulation
    time_data.tf = 0.01;                      //final time (s). End of simulation

    physical_params.gamma = 1.0;                    //Heat conductivity (W/mK)
    physical_params.rho = 1.0;                      //material density (kg/m3)
    physical_params.Cp = 1.0;                       //heat capacity (J/kgK)

    /*Allocate memory for solver data*/
    solver_data.X =  matrix2D(grid_data.nx+1,grid_data.ny+1);
    solver_data.Y =  matrix2D(grid_data.nx+1,grid_data.ny+1);
    solver_data.T =  matrix2D(grid_data.nx+1,grid_data.ny+1);

    /*Set boundaries and source term*/
    boundaries.source = q;
    boundaries.Tnfunc = Tnfunc;
    boundaries.Tsfunc = Tsfunc;
    boundaries.Tefunc = Tefunc;
    boundaries.Twfunc = Twfunc;

    /*Execute solver*/
    heat2D(grid_data, time_data, physical_params, boundaries, solver_data);

    /*Free allocated data*/
    free2Df(solver_data.X, grid_data.nx+1);
    free2Df(solver_data.Y, grid_data.nx+1);
    free2Df(solver_data.T, grid_data.nx+1);

    return 0;
}
