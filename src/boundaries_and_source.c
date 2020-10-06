/*
 * parameters_and_boundaries.c
 *
 *  Created on: Oct 6, 2020
 *      Author: d-w-h
 */

double q(double x, double y, double t)
{
    return 0.0;                                        //Enter equation for the source term here
}

double Tnfunc(double x, double t)
{
    return 1.0;                                //Enter the equation for the north boundary condition here
}
double Tsfunc(double x, double t)
{
    return 1.0;                                //Enter the equation for the south boundary condition here
}
double Twfunc(double y, double t)
{
    return 1.0;                                //Enter the equation for the west boundary condition here
}
double Tefunc(double y, double t)
{
    return 1.0;                                //Enter the equation for the east boundary condition here.
}
