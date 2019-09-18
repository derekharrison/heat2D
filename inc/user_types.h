/*
 * user_types.h
 *
 *  Created on: Sep 17, 2019
 *      Author: d-w-h
 */

#ifndef USER_TYPES_H_
#define USER_TYPES_H_

#define SIZE_RESULTS 3

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
	double (*q)(double x, double y, double t);
	double (*Tnfunc)(double x, double t);
	double (*Tsfunc)(double x, double t);
	double (*Twfunc)(double y, double t);
	double (*Tefunc)(double y, double t);
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

#endif /* USER_TYPES_H_ */