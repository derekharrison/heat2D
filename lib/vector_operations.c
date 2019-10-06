/*
 * vector_operations.c
 *
 *  Created on: Sep 18, 2019
 *      Author: d-w-h
 */

/*-----------------------------------------------------------------------------------------------*/
double dot_product(double* x, double* y, int size)
/*
 * Determine the dot product dot_product = x * y
 *
 * input    x
 * input    y
 * input    size
 * return    dot_product
 */
{
	int j;
	double dot_product;

	dot_product = 0.0;
	for (j = 1; j <= size; j++)
		dot_product = dot_product + x[j] * y[j];

	return dot_product;
}

/*-----------------------------------------------------------------------------------------------*/
void vector_sum(double f1, double* vec_1, double f2, double* vec_2, int size, double* vec_sum)
/*
 * Calculate vec_sum = f1 * vec_1 + f2 * vec_2
 *
 * input    f1
 * input    vec_1
 * input    f2
 * input    vec_2
 * input    size
 * output    vec_sum
 */
{
	int j;

	for(j = 1; j <= size; ++j)
	{
		vec_sum[j] = f1 * vec_1[j] + f2 * vec_2[j];
	}

}
