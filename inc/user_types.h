/*
 * user_types.h
 *
 *  Created on: Sep 17, 2019
 *      Author: d-w-h
 */

#ifndef USER_TYPES_H_
#define USER_TYPES_H_

#define SIZE_RESULTS 3

typedef double (*source)(double x, double y, double t);
typedef double (*boundary_temperature)(double x, double t);

enum result_type {
    X_COORDINATES = 0,
    Y_COORDINATES = 1,
    TEMPERATURE = 2,
};

typedef struct grid_parameters {
    double Lx;
    double Ly;
    int nx;
    int ny;
} grid_parameters_t;

typedef struct time_parameters {
    double to;
    double tf;
    double t;
    int maxts;
} time_parameters_t;

typedef struct physical_params {
    double rho;
    double Cp;
    double gamma;
} physical_params_t;

typedef struct boundary_temperatures {
    boundary_temperature Tnfunc;
    boundary_temperature Tsfunc;;
    boundary_temperature Twfunc;
    boundary_temperature Tefunc;
} boundary_temperatures_t;

typedef struct grid_coordinates {
    double** X;
    double** Y;
} grid_coordinates_t;

typedef struct solver_data {
    double **A;
    double **Astor;
    double *y;
    double *z;
    double *p;
    double *x;
    double *xo;
    double *r;
    double **T;
} solver_data_t;

typedef struct solver_results {
    double** X;
    double** Y;
    double** T;
} solver_results_t;

#endif /* USER_TYPES_H_ */
