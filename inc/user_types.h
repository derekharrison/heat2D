/*
 * user_types.h
 *
 *  Created on: Oct 6, 2020
 *      Author: d-w-h
 */

#ifndef USER_TYPES_H_
#define USER_TYPES_H_

typedef struct grid_data {
	double Lx;
	double Ly;
	int nx;
	int ny;
} g_data;

typedef struct time_data {
    double to;
    double tf;
    int maxts;
} t_data;

typedef struct physical_params {
	double gamma;
	double rho;
	double Cp;
} p_params;

typedef struct solver_data {
	double** X;
	double** Y;
	double** T;
} s_data;

typedef struct boundaries {
	double (*source)(double, double, double);
	double (*Tnfunc)(double, double);
	double (*Tsfunc)(double, double);
	double (*Tefunc)(double, double);
	double (*Twfunc)(double, double);
} boundaries_t;

#endif /* USER_TYPES_H_ */
