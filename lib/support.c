/*
 * support.c
 *
 *  Created on: Oct 6, 2020
 *      Author: d-w-h
 */

#include <stdio.h>
#include <stdlib.h>
#include "../inc/user_types.h"

void export_data(g_data grid_data, s_data solver_data)
{
    int i, j;
    /*Exporting data*/
    FILE *file;
    file = fopen("numerical_temperature_data.txt","w");
    if (file != NULL)
    {
        for(j = 1;j <= grid_data.ny; j++)
        {
            for(i = 1; i <= grid_data.nx; i++)
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
        fprintf(file,"%i %i %f %f",grid_data.nx, grid_data.ny, grid_data.Lx, grid_data.Ly);
        fclose(file);
    }
    else
    {
        printf("Could not open file");
    }
}
