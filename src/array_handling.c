/** \file array_handling.c
 *  Contains convient functions for dyanmically allocating 2D and 3D real or double complex arrays
 */
#include <stdlib.h>
#include "array_handling.h"

/** Allocate memeory for a real 2D array
*
* @param ni the number of rows
* @param nj the number of columns
* @return the pointer to the start of the memory block for the 2D array
*/
double** alloc_2D_double_array(int ni, int nj)
{
	double** array = (double**)calloc( ni, sizeof (double*));
	int i;
	for(i = 0; i < ni; i++){
		array[i] = (double*)calloc(nj , sizeof(double));
	}

	return array;
}

/** Free the memory associated with a real 2D array
*
* @param array the array to be freed
* @param ni the number of rows in the array
*/
void free_alloced_2D_double_array(double** array, int ni)
{
	int i;
	for (i = 0; i < ni; i++){
	   free(array[i]);
	}
	free(array);
}

/** Allocate memeory for a double complex 2D array
*
* @param ni the number of rows
* @param nj the number of columns
* @return the pointer to the start of the memory block for the 2D array
*/
double complex** alloc_2D_complex_array(int ni, int nj)
{
	double complex** array = (double complex**)calloc( ni, sizeof (double complex*));
	int i;
	for(i = 0; i < ni; i++){
		array[i] = (double complex*)calloc(nj , sizeof(double complex));
	}

	return array;
}

/** Free the memory associated with a double complex 2D array
*
* @param array the array to be freed
* @param ni the number of rows in the array
*/
void free_alloced_2D_complex_array(double complex** array, int ni)
{
	int i;
	for (i = 0; i < ni; i++){
	   free(array[i]);
	}
	free(array);
}


/** Allocate memeory for a real 3D array
*
* @param ni the number of rows
* @param nj the number of columns
* @param nk the depth of the array block
* @return the pointer to the start of the memory block for the 2D array
*/
double*** alloc_3D_double_array(int ni, int nj, int nk)
{
	double*** array = (double***)calloc( ni, sizeof (double**));
	int i,j;
	for(i = 0; i < ni; i++){
		array[i] = (double**)calloc(nj, sizeof(double*));
		for(j = 0; j < nj; j++){
			array[i][j] = (double*)calloc(nk , sizeof(double));
		}
	}

	return array;
}

/** Free the memory associated with a read 3D array
*
* @param array the array to be freed
* @param ni the number of rows in the array
* @param nj the number of columns in the array
*/
void free_alloced_3D_double_array(double*** array, int ni,  int nj)
{
	int i, j;
	for(i = 0; i < ni; i++){
		for(j = 0; j < nj; j++){
			free(array[i][j]);
		}
		free(array[i]);
	}

}

double**** alloc_4D_double_array(int n1, int n2, int n3, int n4)
{
	double**** array = (double****)calloc( n1, sizeof (double***));
	int i,j,k;
	for(i = 0; i < n1; i++){
		array[i] = (double***)calloc(n2, sizeof(double**));
		for(j = 0; j < n2; j++){
			array[i][j] = (double**)calloc(n3 , sizeof(double*));
			for(k=0; k < n3; k++){
				array[i][j][k] = (double*)calloc(n4, sizeof(double));
			}
		}
	}
				
	return array;
}

void free_alloced_4D_double_array(double**** array, int ni,  int nj, int nk)
{
	int i, j, k;
	for(i = 0; i < ni; i++){
		for(j = 0; j < nj; j++){
			for(k=0; k < nk; k++){
				free(array[i][j][k]);
			}
			free(array[i][j]);
		}
		free(array[i]);
	}

}

/** Allocate memeory for a double complex 3D array
*
* @param ni the number of rows
* @param nj the number of columns
* @param nk the depth of the array block
* @return the pointer to the start of the memory block for the 2D array
*/
double complex*** alloc_3D_complex_array(int ni, int nj, int nk)
{
	double complex*** array = (double complex***)calloc( ni, sizeof (double complex**));
	int i,j;
	for(i = 0; i < ni; i++){
		array[i] = (double complex**)calloc(nj, sizeof(double complex*));
		for(j = 0; j < nj; j++){
			array[i][j] = (double complex*)calloc(nk , sizeof(double complex));
		}
	}

	return array;
}

/** Free the memory associated with a double complex 3D array
*
* @param array the array to be freed
* @param ni the number of rows in the array
* @param nj the number of columns in the array
*/
void free_alloced_3D_complex_array(double complex*** array, int ni,  int nj)
{
	int i, j;
	for(i = 0; i < ni; i++){
		for(j = 0; j < nj; j++){
			free(array[i][j]);
		}
		free(array[i]);
	}

}
