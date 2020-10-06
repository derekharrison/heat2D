/*
 * heat.h
 *
 *  Created on: Oct 6, 2020
 *      Author: d-w-h
 */

#ifndef HEAT_H_
#define HEAT_H_

#include "user_types.h"

void heat2D(g_data grid_data, t_data time_data, p_params physical_params, boundaries_t boundaries, s_data* solver_data);

#endif /* HEAT_H_ */
