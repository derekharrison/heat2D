/*
 * heat2D_support.c
 *
 *  Created on: Sep 18, 2019
 *      Author: d-w-h
 */

#include <math.h>

#include "../inc/memory_functions.h"
#include "../inc/user_types.h"
#include "../inc/vector_operations.h"

/*-----------------------------------------------------------------------------------------------*/
void set_initial_temp(solver_data_t* solver_data, grid_parameters_t grid_parameters)
/*
 * Set initial temperature distribution
 *
 * output    solver_data
 * input    grid_parameters
 */
{
    int j;
    int nx, ny;
    double Ti;

    nx = grid_parameters.nx;
    ny = grid_parameters.ny;

    Ti =  0.0;

    for (j = 1; j <= nx*ny; j++)
        solver_data->xo[j] = Ti;
}

/*-----------------------------------------------------------------------------------------------*/
void generate_grid_coordinates(grid_parameters_t grid_parameters,
                               grid_coordinates_t* grid_coordinates)
/*
 * Set initial temperature distribution
 *
 * input    grid_parameters
 * output    grid_coordinates
 */
{
    int i, j;
    int nx, ny;
    double deltax, deltay;

    nx = grid_parameters.nx;
    ny = grid_parameters.ny;

    deltax = grid_parameters.Lx / nx;
    deltay = grid_parameters.Ly / ny;

    for (i = 1; i <= nx; i++)
        for (j = 1;j <= ny; j++)
        {
            grid_coordinates->X[i][j] = i*deltax-deltax/2;
            grid_coordinates->Y[i][j] = j*deltay-deltay/2;
        }
}

/*-----------------------------------------------------------------------------------------------*/
void calc_coefficient_matrix(grid_parameters_t grid_parameters,
                             grid_coordinates_t grid_coordinates,
                             time_parameters_t time_parameters,
                             physical_params_t physical_params,
                             solver_data_t* solver_data,
                             double (*q)(double, double, double),
                             boundary_temperatures_t boundary_temperatures)
/*
 * Determine coefficient matrix
 *
 * input    grid_parameters
 * input    grid_coordinates
 * input    time_parameters
 * output   solver_data
 * input    q
 * input    boundary_temperatures
 */
{
    int i,j, nn;
    int nx, ny;
    double** X, **Y;
    double deltax, deltay, dt, t;
    double c, a, b, d2, K;

    nx = grid_parameters.nx;
    ny = grid_parameters.ny;

    X = grid_coordinates.X;
    Y = grid_coordinates.Y;

    t = time_parameters.t;

    deltax = grid_parameters.Lx / grid_parameters.nx;
    deltay = grid_parameters.Ly / grid_parameters.ny;
    dt = (time_parameters.tf - time_parameters.to) / time_parameters.maxts;
    c = deltay / deltax;
    a = c * physical_params.gamma;
    b = physical_params.gamma / c;
    d2 = deltax * deltay;
    K = d2 * physical_params.rho * physical_params.Cp / dt;

    //Preliminary ICCG calculations
    //Generating central coefficients and source terms
    for (j = 2; j <= ny - 1; j++)
        for (i = 2; i<= nx - 1; i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = -b;
            solver_data->A[nn][2] = -a;
            solver_data->A[nn][3] = K + 2 * a + 2 * b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] - (solver_data->A[nn][3] * solver_data->xo[nn] -
                                a * solver_data->xo[nn-1] - a * solver_data->xo[nn+1] - b * solver_data->xo[nn+nx] - b*solver_data->xo[nn-nx]);
        }

    //Generating upper central coefficients and source terms
    for (j = ny;j <= ny; j++)
        for (i = 2; i <= nx - 1; i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = -b;
            solver_data->A[nn][2] = -a;
            solver_data->A[nn][3] = K + 2 * a + 3 * b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] + 2 * b * boundary_temperatures.Tnfunc(X[i][j],t) -
                               (solver_data->A[nn][3] * solver_data->xo[nn] - a * solver_data->xo[nn-1] - a * solver_data->xo[nn+1] -
                                b * solver_data->xo[nn-nx]);
        }

    //Generating upper left node coefficients and source terms
    for (j = ny; j <= ny; j++)
        for (i = 1; i <= 1; i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = -b;
            solver_data->A[nn][2] = 0;
            solver_data->A[nn][3] = K + 3 * a + 3 * b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] + 2 * b * boundary_temperatures.Tnfunc(X[i][j],t) +
                                2 * a * boundary_temperatures.Twfunc(Y[i][j],t) - (solver_data->A[nn][3] * solver_data->xo[nn] -
                                a * solver_data->xo[nn+1] - b * solver_data->xo[nn-nx]);
        }

    //Generating upper right node coefficients and source terms
    for (j = ny; j <= ny; j++)
        for (i = nx;i <= nx; i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = -b;
            solver_data->A[nn][2] = -a;
            solver_data->A[nn][3] = K + 3 * a + 3 * b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] + 2 * b * boundary_temperatures.Tnfunc(X[i][j],t) +
                                2 * a * boundary_temperatures.Tefunc(Y[i][j],t) - (solver_data->A[nn][3] * solver_data->xo[nn] -
                                a * solver_data->xo[nn-1] - b * solver_data->xo[nn-nx]);
        }

    //Generating left central node coefficients and source terms
    for (j = 2; j <= ny - 1; j++)
        for (i = 1; i <= 1; i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = -b;
            solver_data->A[nn][2] = 0;
            solver_data->A[nn][3] = K + 3 * a + 2 * b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] + 2 * a * boundary_temperatures.Twfunc(Y[i][j],t) -
                               (solver_data->A[nn][3] * solver_data->xo[nn] - a * solver_data->xo[nn+1] - b * solver_data->xo[nn+nx] -
                                b * solver_data->xo[nn-nx]);
        }

    //Generating lower left node coefficients and source terms
    for (j=1;j<=1;j++)
        for (i=1;i<=1;i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = 0;
            solver_data->A[nn][2] = 0;
            solver_data->A[nn][3] = K + 3 * a + 3 * b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] + 2 * a * boundary_temperatures.Twfunc(Y[i][j],t) +
                                2 * b * boundary_temperatures.Tsfunc(X[i][j],t) - (solver_data->A[nn][3] * solver_data->xo[nn] -
                                a * solver_data->xo[nn+1] - b * solver_data->xo[nn+nx]);
        }

    //Generating lower central node coefficients and source terms
    for (j = 1; j <= 1; j++)
        for (i = 2; i <= nx-1; i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = 0;
            solver_data->A[nn][2] = -a;
            solver_data->A[nn][3] = K + 2 * a + 3 * b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] + 2 * b * boundary_temperatures.Tsfunc(X[i][j],t) -
                               (solver_data->A[nn][3] * solver_data->xo[nn] - a * solver_data->xo[nn-1] - a * solver_data->xo[nn+1] -
                                b * solver_data->xo[nn+nx]);
        }

    //Generating lower right node coefficients and source terms
    for (j = 1; j <= 1; j++)
        for (i = nx; i <= nx; i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = 0;
            solver_data->A[nn][2] = -a;
            solver_data->A[nn][3] = K+3*a+3*b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] + 2 * a * boundary_temperatures.Tefunc(Y[i][j],t) +
                                2 * b * boundary_temperatures.Tsfunc(X[i][j],t)-(solver_data->A[nn][3] * solver_data->xo[nn] -
                                a * solver_data->xo[nn-1]- b * solver_data->xo[nn+nx]);
        }

    //Generating right central node coefficients and source terms
    for (j = 2; j <= ny-1; j++)
        for (i = nx;i <= nx; i++)
        {
            nn = i + (j - 1) * nx;
            solver_data->A[nn][1] = -b;
            solver_data->A[nn][2] = -a;
            solver_data->A[nn][3] = K + 3 * a + 2 * b;
            solver_data->r[nn] = q(X[i][j],Y[i][j],t) * d2 + K * solver_data->xo[nn] + 2 * a * boundary_temperatures.Tefunc(Y[i][j],t) -
                               (solver_data->A[nn][3] * solver_data->xo[nn] - a*solver_data->xo[nn-1] - b * solver_data->xo[nn+nx] -
                                b * solver_data->xo[nn-nx]);
        }

}

/*-----------------------------------------------------------------------------------------------*/
void store_coefficient_matrix(solver_data_t* solver_data,
                              grid_parameters_t grid_parameters)
/*
 * Store coefficient matrix
 *
 * output    solver_data
 * input     grid_parameters
 */
{
    int i, j;
    int nx, nt;

    nx = grid_parameters.nx;
    nt = nx * grid_parameters.ny;

    for (j = 1;j <= nt; j++)
        for (i = 1;i <= 3; i++)
            solver_data->Astor[j][i] = solver_data->A[j][i];
}

/*-----------------------------------------------------------------------------------------------*/
void incomplete_cholesky_factorization(solver_data_t* solver_data,
                                       grid_parameters_t grid_parameters)
/*
 * incomplete cholesky factorization
 *
 * output    solver_data
 * input     grid_parameters
 */
{
    int i, j;
    int nx, nt;

    nx = grid_parameters.nx;
    nt = nx * grid_parameters.ny;

    solver_data->A[1][1] = 0;
    solver_data->A[1][2] = 0;
    solver_data->A[1][3] = sqrt(solver_data->A[1][3]);

    for (j=2;j<=nt;j++)
        for (i=1;i<=3;i++)
            if (solver_data->A[j][i] != 0)
            {
                if (i==1 && j>=1+nx)
                    solver_data->A[j][i] = solver_data->A[j][i] / solver_data->A[j-nx][3];
                else if (i==2)
                    solver_data->A[j][i] = solver_data->A[j][i] / solver_data->A[j-1][3];
                else if (i==3)
                    solver_data->A[j][i]=sqrt(solver_data->A[j][i] - solver_data->A[j][i-2] * solver_data->A[j][i-2] -
                    solver_data->A[j][i-1] * solver_data->A[j][i-1]);
            }

}

/*-----------------------------------------------------------------------------------------------*/
void solve_Mz_is_r(solver_data_t* solver_data,
                   grid_parameters_t grid_parameters)
/*
 * Solve Mz = r
 *
 * output    solver_data
 * input     grid_parameters
 */
{
    int j, nu;
    int nx, nt;

    nx = grid_parameters.nx;
    nt = nx * grid_parameters.ny;

    solver_data->y[1] = solver_data->r[1] / solver_data->A[1][3];
    for (j = 2; j <=nx ; j++)
        solver_data->y[j] = (solver_data->r[j] - solver_data->A[j][2] * solver_data->y[j-1]) / solver_data->A[j][3];
    for (j=nx+1;j<=nt;j++)
        solver_data->y[j] = (solver_data->r[j] - solver_data->A[j][1] * solver_data->y[j-nx] -
                            solver_data->A[j][2] * solver_data->y[j-1]) / solver_data->A[j][3];

    //Solving L'z=y
    solver_data->z[nt] = solver_data->y[nt] / solver_data->A[nt][3];
    for (j = 2; j<=nx ; j++)
    {
        nu = nt + 1 - j;
        solver_data->z[nu] = (solver_data->y[nu] - solver_data->A[nu+1][2] * solver_data->z[nu+1]) / solver_data->A[nu][3];
    }
    for (j=nx+1;j<=nt;j++)
    {
        nu = nt + 1 - j;
        solver_data->z[nu] = (solver_data->y[nu] - solver_data->A[nu+1][2] * solver_data->z[nu+1] -
                             solver_data->A[nu+nx][1] * solver_data->z[nu+nx]) / solver_data->A[nu][3];
    }

}

/*-----------------------------------------------------------------------------------------------*/
void calc_Ap(solver_data_t* solver_data,
             grid_parameters_t grid_parameters,
             double* Astorp)
/*
 * Calculate A*p = Ap
 *
 * output    solver_data
 * input    grid_parameters
 * output    Astorp
 */
{
    int j;
    int nx, nt;

    nx = grid_parameters.nx;
    nt = nx * grid_parameters.ny;

    Astorp[1] = solver_data->Astor[1][3] * solver_data->p[1] + solver_data->Astor[2][2] * solver_data->p[2] +
                solver_data->Astor[1+nx][1] * solver_data->p[1+nx];

    for (j = 2; j <= nx; j++)
        Astorp[j] = solver_data->Astor[j][3] * solver_data->p[j] + solver_data->Astor[j][2] * solver_data->p[j-1] +
                    solver_data->Astor[j+1][2] * solver_data->p[j+1] + solver_data->Astor[j+nx][1] * solver_data->p[j+nx];

    for (j = nx + 1; j <= nt - nx; j++)
        Astorp[j] = solver_data->Astor[j][3] * solver_data->p[j] + solver_data->Astor[j][2] * solver_data->p[j-1] +
                    solver_data->Astor[j][1] * solver_data->p[j-nx] + solver_data->Astor[j+1][2] * solver_data->p[j+1] +
                    solver_data->Astor[j+nx][1] * solver_data->p[j+nx];

    for (j = nt - nx + 1; j <= nt - 1; j++)
        Astorp[j] = solver_data->Astor[j][3] * solver_data->p[j] + solver_data->Astor[j][2] * solver_data->p[j-1] +
                    solver_data->Astor[j][1] * solver_data->p[j-nx] + solver_data->Astor[j+1][2] * solver_data->p[j+1];

    Astorp[nt] = solver_data->Astor[nt][3] * solver_data->p[nt] + solver_data->Astor[nt][2] * solver_data->p[nt-1] +
                 solver_data->Astor[nt][1] * solver_data->p[nt-nx];

}

/*-----------------------------------------------------------------------------------------------*/
void perform_cgm(solver_data_t* solver_data,
                 grid_parameters_t grid_parameters,
                 double* epsilon)
/*
 * Perform conjugate gradient method
 *
 * output    solver_data
 * input    grid_parameters
 * output    epsilon
 */
{
    int nx, nt;
    double delold, *Astorp, pAstorp, alpha;
    double delnew, B;

    nx = grid_parameters.nx;
    nt = nx * grid_parameters.ny;

    Astorp = matrix1D(nt+1);

    calc_Ap(solver_data,
            grid_parameters,
            Astorp);

    delold = dot_product(solver_data->r,
                          solver_data->z,
                          nt);

    pAstorp = dot_product(solver_data->z,
                          Astorp,
                          nt);

    alpha = delold/(pAstorp);

    vector_sum(1.0,
               solver_data->x,
               alpha,
               solver_data->z,
               nt,
               solver_data->x);

    vector_sum(1.0,
               solver_data->r,
               -alpha,
               Astorp,
               nt,
               solver_data->r);

    solve_Mz_is_r(solver_data,
                  grid_parameters);

    delnew = dot_product(solver_data->r,
                          solver_data->z,
                          nt);

    B = delnew/delold;

    vector_sum(1.0,
               solver_data->z,
               B,
               solver_data->p,
               nt,
               solver_data->p);

    *epsilon = dot_product(solver_data->r,
                           solver_data->z,
                           nt);

    *epsilon = sqrt(*epsilon / nt);

    free_matrix1D(Astorp);

}

/*-----------------------------------------------------------------------------------------------*/
void set_solver_data(solver_data_t* solver_data,
                     grid_parameters_t grid_parameters)
/*
 * Set p and x solver data
 *
 * output    solver_data
 * input    grid_parameters
 */
{
    int j, nt, nx;

    nx = grid_parameters.nx;
    nt = nx * grid_parameters.ny;

    for (j = 1; j <= nt; j++)
    {
        solver_data->p[j] = solver_data->z[j];
        solver_data->x[j] = solver_data->xo[j];
    }

}

/*-----------------------------------------------------------------------------------------------*/
void set_current_temperature(solver_data_t* solver_data,
                             grid_parameters_t grid_parameters)
/*
 * Set current temperature
 *
 * output    solver_data
 * input    grid_parameters
 */
{
    int j, nt, nx;

    nx = grid_parameters.nx;
    nt = nx * grid_parameters.ny;

    for (j = 1; j <= nt; j++)
        solver_data->xo[j] = solver_data->x[j];

}

/*-----------------------------------------------------------------------------------------------*/
void set_result_data(solver_data_t* solver_data,
                     grid_coordinates_t* grid_coordinates,
                     grid_parameters_t grid_parameters,
                     solver_results_t* solver_results)
/*
 * Set result data
 *
 * input    solver_data
 * input    grid_coordinates
 * input    grid_parameters
 * output   solver_results
 */
{
    int i, j, nx, ny;

    nx = grid_parameters.nx;
    ny = grid_parameters.ny;

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            solver_results->T[i][j] = solver_data->x[i+(j-1)*nx];
            solver_results->grid_coordinates.X[i][j] = grid_coordinates->X[i][j];
            solver_results->grid_coordinates.Y[i][j] = grid_coordinates->Y[i][j];
        }

}
