/*
 * memory_functions.h
 *
 *  Created on: Oct 6, 2020
 *      Author: d-w-h
 */

#ifndef MEMORY_FUNCTIONS_H_
#define MEMORY_FUNCTIONS_H_

double* matrix1D(int np);
double** matrix2D( int nx, int ny);
void free2Df(double** m, int nx);

#endif /* MEMORY_FUNCTIONS_H_ */
