#include <complex.h>

double** alloc_2D_double_array(int ni, int nj);
void free_alloced_2D_double_array(double** array, int ni);

double complex** alloc_2D_complex_array(int ni, int nj);
void free_alloced_2D_complex_array(double complex** array, int ni);

double*** alloc_3D_double_array(int ni, int nj, int nk);
void free_alloced_3D_double_array(double*** array, int ni,  int nj);

double complex*** alloc_3D_complex_array(int ni, int nj, int nk);
void free_alloced_3D_complex_array(double complex*** array, int ni,  int nj);

double**** alloc_4D_double_array(int n1, int n2, int n3, int n4);
void free_alloced_4D_double_array(double**** array, int ni,  int nj, int nk);


