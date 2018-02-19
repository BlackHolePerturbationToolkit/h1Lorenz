#include <complex.h>
#include <gsl/gsl_linalg.h>

int cisnan(double complex c);
double factorial(int x);
int sign(double num);
double normalized_associated_legendre_P(int l, int m, double theta);
double normalized_associated_legendre_P_theta_deriv(int l, int m, double theta);
double normalized_associated_legendre_P_theta_2nd_deriv(int l, int m, double theta);
double elliptic_integral_1st(double k);
double elliptic_integral_2nd(double k);
double legendre_P(int l, double x);
double legendre_Q(int l, double x);
double polygamma(int n, double x);
double lambertW(double x);
double clebsch_gordan_coefficient(int twice_j1, int twice_j2, int twice_j3, int twice_m1, int twice_m2, int twice_m3);
double real_matrix_condition_number(gsl_matrix *A, int m, int n);
double complex_matrix_condition_number(gsl_matrix_complex *A, int m, int n);
double complex_matrix_Frobenius_norm(gsl_matrix_complex *A, int m, int n);
