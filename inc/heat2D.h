/*
 * heat2D.h
 *
 *  Created on: Sep 18, 2019
 *      Author: d-w-h
 */

#ifndef HEAT2D_H_
#define HEAT2D_H_

#include "user_types.h"

double*** heat2D(grid_parameters_t grid_parameters,
				 time_parameters_t time_parameters,
				 physical_params_t physical_params,
				 boundary_temperatures_t boundary_temperatures,
				 double (*source_equation) (double, double, double));

#endif /* HEAT2D_H_ */