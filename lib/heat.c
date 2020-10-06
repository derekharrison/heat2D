/*
 * heat.c
 *
 *  Created on: Oct 6, 2020
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../inc/memory_functions.h"
#include "../inc/user_types.h"

void heat2D(g_data grid_data, t_data time_data, p_params physical_params, boundaries_t boundaries, s_data solver_data)
{
    double Lx, Ly, gamma, tolerance, Ti;
    double rho, Cp, to, tf;
    int imax, nx, ny, maxts;

    /*Adjustable parameters*/
    Lx = grid_data.Lx;                       //Length of domain (m) along x coordinate
    Ly = grid_data.Ly;                       //length of domain (m) along y coordinate
    nx = grid_data.nx;                         //Amount of nodes along x coordinate
    ny = grid_data.ny;                         //amount of nodes along y coordinate

    maxts = time_data.maxts;                     //amount of timesteps
    to = time_data.to;                       //initial time (s). Start of simulation
    tf = time_data.tf;                      //final time (s). End of simulation

    gamma = physical_params.gamma;                    //Heat conductivity (W/mK)
    rho = physical_params.rho;                      //material density (kg/m3)
    Cp = physical_params.Cp;                       //heat capacity (J/kgK)

    imax = 500;                     //Maximum iterations ICCG
    tolerance = 1e-30;              //Tolerance
    Ti=0.0;                         //Initial temperature guess
    /*End adjustable parameters*/

    // Generating coefficient and source term matrices
    clock_t begin, end;
    double time_spent;

    begin = clock();

    double deltax, deltay, a, c, K, dt, d2, b;
    int nt;
    deltax = Lx / nx;
    deltay = Ly / ny;
    dt = (tf - to) / maxts;
    c = deltay / deltax;
    a=c * gamma;
    b = gamma / c;
    d2 = deltax * deltay;
    K = d2 * rho * Cp / dt;
    nt = nx * ny;

    double **A, **Astor, *Astorp, *y, *z, *p, *x, *r, *xo;
    A = matrix2D(nt+1,3+1);
    Astor = matrix2D(nt+1,3+1);
    y = matrix1D(nt+1);
    z = matrix1D(nt+1);
    p = matrix1D(nt+1);
    Astorp = matrix1D(nt+1);
    x = matrix1D(nt+1);
    xo = matrix1D(nx*ny+1);
    r = matrix1D(nt+1);

    int i, j;
    for (j = 1; j <= nx * ny; j++)
        xo[j]=Ti;

    for (i = 1; i<= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            solver_data.X[i][j] = i*deltax-deltax/2;
            solver_data.Y[i][j] = j*deltay-deltay/2;
        }

    /*Entering ICCG loop*/
    double epsilon, delold, delnew, pAstorp;
    double alpha, B, t;
    int nn, nu, it, countert;

    countert = 0;
    t = to;

    do{
        //Generating central coefficients and source terms
        for (j = 2; j <= ny - 1; j++)
            for (i = 2; i <= nx - 1; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = -b;
                A[nn][2] = -a;
                A[nn][3] = K+2*a+2*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]-(A[nn][3]*xo[nn]-a*xo[nn-1]-a*xo[nn+1]-b*xo[nn+nx]-b*xo[nn-nx]);
            }

        //Generating upper central coefficients and source terms
        for (j = ny; j <= ny; j++)
            for (i = 2; i <= nx - 1; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = -b;
                A[nn][2] = -a;
                A[nn][3] = K+2*a+3*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]+2*b*boundaries.Tnfunc(solver_data.X[i][j],t)-(A[nn][3]*xo[nn]-a*xo[nn-1]-a*xo[nn+1]-b*xo[nn-nx]);
            }

        //Generating upper left node coefficients and source terms
        for (j = ny; j <= ny; j++)
            for (i = 1; i<= 1; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = -b;
                A[nn][2] = 0;
                A[nn][3] = K+3*a+3*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]+2*b*boundaries.Tnfunc(solver_data.X[i][j],t)+2*a*boundaries.Twfunc(solver_data.Y[i][j],t)-(A[nn][3]*xo[nn]-a*xo[nn+1]-b*xo[nn-nx]);
            }

        //Generating upper right node coefficients and source terms
        for (j = ny; j<= ny; j++)
            for (i = nx; i <= nx; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = -b;
                A[nn][2] = -a;
                A[nn][3] = K+3*a+3*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]+2*b*boundaries.Tnfunc(solver_data.X[i][j],t)+2*a*boundaries.Tefunc(solver_data.Y[i][j],t)-(A[nn][3]*xo[nn]-a*xo[nn-1]-b*xo[nn-nx]);
            }

        //Generating left central node coefficients and source terms
        for (j = 2; j <= ny - 1; j++)
            for (i = 1; i <= 1; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = -b;
                A[nn][2] = 0;
                A[nn][3] = K+3*a+2*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]+2*a*boundaries.Twfunc(solver_data.Y[i][j],t)-(A[nn][3]*xo[nn]-a*xo[nn+1]-b*xo[nn+nx]-b*xo[nn-nx]);
            }

        //Generating lower left node coefficients and source terms
        for (j = 1; j <= 1; j++)
            for (i = 1; i<= 1; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = 0;
                A[nn][2] = 0;
                A[nn][3] = K+3*a+3*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]+2*a*boundaries.Twfunc(solver_data.Y[i][j],t)+2*b*boundaries.Tsfunc(solver_data.X[i][j],t)-(A[nn][3]*xo[nn]-a*xo[nn+1]-b*xo[nn+nx]);
            }

        //Generating lower central node coefficients and source terms
        for (j = 1; j <= 1; j++)
            for (i = 2; i <= nx - 1; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = 0;
                A[nn][2] = -a;
                A[nn][3] = K+2*a+3*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]+2*b*boundaries.Tsfunc(solver_data.X[i][j],t)-(A[nn][3]*xo[nn]-a*xo[nn-1]-a*xo[nn+1]-b*xo[nn+nx]);
            }

        //Generating lower right node coefficients and source terms
        for (j = 1;j <= 1; j++)
            for (i = nx; i <= nx; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = 0;
                A[nn][2] = -a;
                A[nn][3] = K+3*a+3*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]+2*a*boundaries.Tefunc(solver_data.Y[i][j],t)+2*b*boundaries.Tsfunc(solver_data.X[i][j],t)-(A[nn][3]*xo[nn]-a*xo[nn-1]-b*xo[nn+nx]);
            }

        //Generating right central node coefficients and source terms
        for (j = 2; j <= ny - 1; j++)
            for (i = nx; i <= nx; i++)
            {
                nn = i+(j-1)*nx;
                A[nn][1] = -b;
                A[nn][2] = -a;
                A[nn][3] = K+3*a+2*b;
                r[nn] = boundaries.source(solver_data.X[i][j],solver_data.Y[i][j],t)*d2+K*xo[nn]+2*a*boundaries.Tefunc(solver_data.Y[i][j],t)-(A[nn][3]*xo[nn]-a*xo[nn-1]-b*xo[nn+nx]-b*xo[nn-nx]);
            }

        //Storing coefficient matrix entries of A in Astor
        for (j = 1; j <= nt; j++)
            for (i = 1; i <= 3; i++)
                Astor[j][i] = A[j][i];


        // Incomplete Cholesky factorization of coefficient matrix A
        A[1][1]=0;
        A[1][2]=0;
        A[1][3]=sqrt(A[1][3]);

        for (j = 2; j <= nt; j++)
            for (i = 1; i <= 3; i++)
                if (A[j][i] != 0)
                {
                    if (i == 1 && j >= 1 + nx)
                        A[j][i]=A[j][i]/A[j-nx][3];
                    else if (i == 2)
                        A[j][i]=A[j][i]/A[j-1][3];
                    else if (i == 3)
                        A[j][i]=sqrt(A[j][i]-A[j][i-2]*A[j][i-2]-A[j][i-1]*A[j][i-1]);
                }

        /* ICCG main iterations */
        //epsilon = r'*r;
        epsilon = 0;
        for (j = 1;j <= nt; j++)
            epsilon = epsilon+r[j]*r[j];

        //Solving Mz=r
        //Solving Ly=r (y = L'z)
        y[1] = r[1]/A[1][3];
        for (j = 2; j <= nx; j++)
            y[j] = (r[j]-A[j][2]*y[j-1])/A[j][3];
        for (j = nx + 1; j <= nt; j++)
            y[j] = (r[j]-A[j][1]*y[j-nx]-A[j][2]*y[j-1])/A[j][3];

        //Solving L'z=y
        z[nt] = y[nt]/A[nt][3];
        for (j = 2; j <= nx; j++)
        {
            nu = nt+1-j;
            z[nu] = (y[nu]-A[nu+1][2]*z[nu+1])/A[nu][3];
        }
        for (j = nx + 1; j <= nt; j++)
        {
            nu = nt+1-j;
            z[nu] = (y[nu]-A[nu+1][2]*z[nu+1]-A[nu+nx][1]*z[nu+nx])/A[nu][3];
        }
        //End solving Mz=r

        for (j = 1;j <= nt; j++)
        {
            p[j] = z[j];
            x[j] = xo[j];
        }

        it = 0;                //ICCG iteration counter

        do
        {
            //Calculating Astor*p
            Astorp[1] = Astor[1][3]*p[1]+Astor[2][2]*p[2]+Astor[1+nx][1]*p[1+nx];
            for (j = 2; j <= nx; j++)
                Astorp[j] = Astor[j][3]*p[j]+Astor[j][2]*p[j-1]+Astor[j+1][2]*p[j+1]+Astor[j+nx][1]*p[j+nx];
            for (j = nx + 1; j <= nt - nx; j++)
                Astorp[j]=Astor[j][3]*p[j]+Astor[j][2]*p[j-1]+Astor[j][1]*p[j-nx]+Astor[j+1][2]*p[j+1]+Astor[j+nx][1]*p[j+nx];
            for (j = nt - nx + 1;j <= nt - 1; j++)
                Astorp[j]=Astor[j][3]*p[j]+Astor[j][2]*p[j-1]+Astor[j][1]*p[j-nx]+Astor[j+1][2]*p[j+1];
            Astorp[nt]=Astor[nt][3]*p[nt]+Astor[nt][2]*p[nt-1]+Astor[nt][1]*p[nt-nx];
            //End of calculation of Astor*p

            delold = 0.0;
            for (j = 1; j <= nt; j++)
                delold = delold+r[j]*z[j];

            //p'*Astorp
            pAstorp = 0.0;
            for (j = 1; j <= nt; j++)
                pAstorp = pAstorp+p[j]*Astorp[j];

            alpha = delold/(pAstorp);

            //x = x+alpha*p;
            for (j = 1; j <= nt; j++)
                x[j] = x[j]+alpha*p[j];

            //r = r-alpha*Astorp;
            for (j = 1; j <= nt; j++)
                r[j] = r[j]-alpha*Astorp[j];

            //Solving Mz=r
            //Solving Ly=r (y = L'z)
            y[1] = r[1]/A[1][3];
            for (j = 2; j <= nx; j++)
                y[j] = (r[j]-A[j][2]*y[j-1])/A[j][3];
            for (j = nx + 1; j <= nt; j++)
                y[j] = (r[j]-A[j][1]*y[j-nx]-A[j][2]*y[j-1])/A[j][3];

            //Solving L'z=y
            z[nt] = y[nt]/A[nt][3];
            for (j = 2; j <= nx; j++)
            {
                nu = nt+1-j;
                z[nu] = (y[nu]-A[nu+1][2]*z[nu+1])/A[nu][3];
            }
            for (j = nx + 1; j <= nt; j++)
            {
                nu = nt+1-j;
                z[nu] = (y[nu]-A[nu+1][2]*z[nu+1]-A[nu+nx][1]*z[nu+nx])/A[nu][3];
            }
            //End solving Mz=r

            delnew = 0.0;
            for (j = 1; j <= nt; j++)
                delnew = delnew+r[j]*z[j];

            B = delnew/delold;

            //p = z + B*p;
            for (j = 1; j <= nt; j++)
                p[j] = z[j]+B*p[j];

            //r'*r
            epsilon = 0.0;
            for (j = 1; j <= nt; j++)
                epsilon = epsilon+r[j]*r[j];

            epsilon = sqrt(epsilon/nt);
            it = it + 1;

        }while (it < imax && epsilon > tolerance);
        /*End of ICCG iteration loop. Advance timestep */

        for (j = 1; j <= nt; j++)
            xo[j] = x[j];

        ++countert;
        t += dt;
    }while(countert < maxts);

    /*Processing results*/
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        	solver_data.T[i][j] = x[i+(j-1)*nx];

    /*Exporting data*/
    FILE *file;
    file = fopen("numerical_temperature_data.txt","w");
    if (file != NULL)
    {
        for(j = 1;j <= ny; j++)
        {
            for(i = 1; i <= nx; i++)
            {
                fprintf(file,"%f %f %f",solver_data.X[i][j], solver_data.Y[i][j], solver_data.T[i][j]);
                fprintf(file,"\n");
            }
        }

        fclose(file);
    }
    else
    {
        printf("Could not open file");
    }

    file = fopen("grid_data.txt","w");
    if (file != NULL)
    {
        fprintf(file,"%i %i %f %f",nx, ny, Lx, Ly);
        fclose(file);
    }
    else
    {
        printf("Could not open file");
    }

    /*Freeing memory*/
    free(z);free(y);free(r);free2Df(A, nx+1);free2Df(Astor, nx+1);
    free(Astorp);free(x);free(xo);free(p);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("\niterations: %d\nrunning time: %f\ncountert: %d\ndone!\n",it,time_spent,countert);
}
