/** \file math_funcs.c
 *  Various useful mathematic functions. 
 *  All the function detailed here work to the accuracy set by DESIRED_FUNCTION_ACCURACY with the exception of legendre_P for x > 1.
 */

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_linalg.h>
#include "math_funcs.h"

extern const double DESIRED_FUNCTION_ACCURACY;


double complex_matrix_Frobenius_norm(gsl_matrix_complex *A, int m, int n)
{
	double sum = 0;
	int i,j;
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			gsl_complex value = gsl_matrix_complex_get(A, i, j);
			double complex a_ij = GSL_REAL(value) + I*GSL_IMAG(value);
			sum += pow(cabs(a_ij), 2.0);
		}
	}

	return sqrt(sum);
}


double complex_matrix_condition_number(gsl_matrix_complex *A, int m, int n)
{
	double norm_of_A = complex_matrix_Frobenius_norm(A, m, n);

	gsl_matrix_complex *A_copy = gsl_matrix_complex_calloc(m, n);
	gsl_matrix_complex_memcpy(A_copy, A);

	int s;

	gsl_permutation * p = gsl_permutation_alloc (m);
	gsl_linalg_complex_LU_decomp (A_copy, p, &s);

	gsl_matrix_complex *A_inverse = gsl_matrix_complex_calloc(m, n);

	gsl_linalg_complex_LU_invert (A_copy, p, A_inverse);

	double norm_of_A_inverse = complex_matrix_Frobenius_norm(A_inverse, m, n);

	gsl_permutation_free(p);
	gsl_matrix_complex_free(A_copy);
	gsl_matrix_complex_free(A_inverse);

	return norm_of_A*norm_of_A_inverse;
}

// This will probably give a different number from the complex method above due to different matrix norms
double real_matrix_condition_number(gsl_matrix *A, int m, int n)
{
	gsl_matrix *V = gsl_matrix_calloc(m, n);
	gsl_vector *S = gsl_vector_calloc(n);
	gsl_vector *work = gsl_vector_calloc(n);

	gsl_linalg_SV_decomp (A, V, S, work);

	double largest = -1e308;
	double smallest = 1e308;

	int i;
	for(i = 0; i < n; i++){
		double value = gsl_vector_get(S, i);
		if( value > largest) largest = value;
		if( value < smallest) smallest = value;
	}

	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);

	return largest/smallest;

}


int cisnan(double complex c)
{
	if(isnan(creal(c)) || isnan(cimag(c)) ) return 1;
	else return 0;
}

/** Computes the factorial of a given integer
 *
 * @param x the integer that you wish to know the factorial of
 * @return the factorial of the supplied integer
 */
double factorial(int x){
	if(x == 0 || x == 1) return 1.0;
	else{
		int i;
		double ans = 1.0;
		for(i = 2; i <= x; i++){
			ans *= i;
		}
		return ans;
	}
}

/** Returns the sign of a number (-1,0, or 1)
 * 
 * @param num the number you wish to know the sign of
 * @return the sign of the supplied number. 1 for positive, -1 for negative and 0 for zero.
 */
int sign(double num)
{
	if(num == 0) return 0;
	if(num > 0) return 1;
	else return -1;
}

/** Computes the value of normalised associated Legendre polynomial for l,m and theta;
 *
 * @param l the spherical harmonic l
 * @param m the spherical harmonic m
 * @param theta the theta angle on the sphere
 * @return the value of normalized associated legendre polynomial for the given values
 */
double normalized_associated_legendre_P(int l, int m, double theta)		//FIXME this function gets called a lot, would make more sense to calculate the answer once and store results in an array
{

	if( m > l ) return 0;
	if(theta == M_PI/2.0)
	{
		if((l+m)%2 == 1) return 0; 
	}

	gsl_sf_result result;
	int status = gsl_sf_legendre_sphPlm_e(l, abs(m), cos(theta), &result);
	if(status != GSL_SUCCESS){
		printf("\nError computer normaized associated LegendreP_lm\n");
		exit(0);
	}
	if( result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for normalized assoicated LegendreP_lm not reached.\n");
		printf("Error in normalized associated LegendreP_lm: %.20lf\n\n", result.err);
		exit(0);
	}

	if(m >= 0) return result.val;
	else return pow(-1.0, abs(m))*result.val;
}


double associated_legendre_P(int l, int m, double theta)		//FIXME this function gets called a lot, would make more sense to calculate the answer once and store results in an array
{

	if( abs(m) > l ) return 0;
	if(theta == M_PI/2.0)
	{
		if((l+m)%2 == 1) return 0; 
	}

	gsl_sf_result result;
	int status = gsl_sf_legendre_Plm_e(l, abs(m), cos(theta), &result);
	if(status != GSL_SUCCESS){
		printf("\nError computer normaized associated LegendreP_lm\n");
		exit(0);
	}
	if( result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for assoicated LegendreP_lm not reached.\n");
		printf("Error in associated LegendreP_lm: %.20lf\n\n", result.err);
		exit(0);
	}

	if(m >= 0) return result.val;
	else return pow(-1.0, abs(m))*result.val;
}


struct lm {
	int l;
	int m;
};


double associated_legendre_P_theta_deriv(int l, int m, double theta)
{
	double x = cos(theta);
	//return 0.5*(associated_legendre_P(l,m+1,theta) - (l+m)*(l-m+1)*associated_legendre_P(l,m-1,theta));


	//return associated_legendre_P(l, m+1, theta) + m * cos(theta)/sin(theta) * associated_legendre_P(l,m,theta);
	//return -sin(theta)*(sqrt(1.0-x*x) * associated_legendre_P(l, m+1, theta) + m * x * associated_legendre_P(l,m,theta))/(x*x-1.0);
	return -(-l*x*associated_legendre_P(l,m,theta) + (l+m)*associated_legendre_P(l-1,m,theta))/(1.0-x*x);			// Can calculate with coupling between l-modes too
}

double normalized_associated_legendre_P_theta_deriv(int l, int m, double theta)
{
	return sqrt( (2.0*l + 1.0)*factorial(l-m)/(4.0*M_PI*factorial(l+m)) ) * associated_legendre_P_theta_deriv(l,m,theta);
}

double normalized_associated_legendre_P_theta_2nd_deriv(int l, int m, double theta)
{
	double x = cos(theta);
	double normalization = sqrt( (2.0*l + 1.0)*factorial(l-m)/(4.0*M_PI*factorial(l+m)) );

	double P 		= associated_legendre_P(l,m,theta);
	double Pprime 	= associated_legendre_P_theta_deriv(l,m,theta)/(-sin(theta));
	double P1		= associated_legendre_P(l,m+1,theta);
	double P1prime	= associated_legendre_P_theta_deriv(l,m+1,theta)/(-sin(theta));

	return normalization * sin(theta) * ( (x*x-1.0)*(m*x*Pprime + sqrt(1.0-x*x)*P1prime) - (m+2.0*m*x*x)*P - 2.0*x*sqrt(1.0-x*x)*P1 )/pow(1.0-x*x,2.5);
}

/** Computes the value of the Elliptic Integral of the 1st kind for given k.
 *
 * @param k the elliptic function parameter
 * @return the value of the Elliptic Integral of the 1st kind for the given k
 */
double elliptic_integral_1st(double k)
{
	gsl_sf_result result;
	int status = gsl_sf_ellint_Kcomp_e(k, GSL_PREC_DOUBLE, &result);
	if(status != GSL_SUCCESS){
		printf("\nError computing Elliptic Integral of the 1st Kind\n");
		exit(0);
	}
	if(result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for Elliptic Integral of the 1st Kind not reached.\n");
		printf("Error in Elliptic Integral of the 1st Kind: %.20lf\n\n", result.err);
		exit(0);
	}
	return result.val;
}

/** Computes the value of the Elliptic Integral of the 2nd kind for given k.
 *
 * @param k the elliptic function parameter
 * @return the value of the Elliptic Integral of the 2nd kind for the given k
 */
double elliptic_integral_2nd(double k)
{
	gsl_sf_result result;
	int status = gsl_sf_ellint_Ecomp_e(k, GSL_PREC_DOUBLE, &result);
	if(status != GSL_SUCCESS){
		printf("\nError computing Elliptic Integral of the 2nd Kind\n");
		exit(0);
	}
	if(result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for Elliptic Integral of the 2nd Kind not reached.\n");
		printf("Error in Elliptic Integral of the 2nd Kind: %.20lf\n\n", result.err);
		exit(0);
	}
	return result.val;
}

/** Computes the value of the Elliptic Integral of the 3rd kind for given n and k
 *
 *
 *
 */
double elliptic_integral_3rd(double k, double n)
{
	gsl_sf_result result;
	int status = gsl_sf_ellint_Pcomp_e(k, n, GSL_PREC_DOUBLE, &result);
	if(status != GSL_SUCCESS){
		printf("\nError computing Elliptic Integral of the 3rd Kind\n");
		exit(0);
	}
	if(result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for Elliptic Integral of the 3rd Kind not reached.\n");
		printf("Error in Elliptic Integral of the 3rd Kind: %.20lf\n\n", result.err);
		exit(0);
	}
	return result.val;

}

/** Computes the value of the 1st Legendre function P_{l}(x) for given l and x
 *
 * @param l the spherical harmonic l parameter
 * @param x the x value
 * @return the value of the 1st Legendre function
 */
double legendre_P(int l, double x)
{
	if(l == 0) return 1.0;
 	if(l == 1) return x;
	// The gsl will only calculate the Legendre P for |x| <= 1
	if(fabs(x) <= 1.0){
		gsl_sf_result result;
		int status = gsl_sf_legendre_Pl_e(l, x, &result);
		if(status != GSL_SUCCESS){
			printf("Error computing LegendreP\n");
			exit(0);
		}
		if(result.err/result.val > DESIRED_FUNCTION_ACCURACY){
			printf("\nDesired function accuracy for LegendreP not reached.\n");
			printf("Error in LegendreP: %.20lf\n\n", result.err);
			exit(0);
		}
		return result.val;
	}else{
		double Plminus2 = 1.0;
		double Plminus1 = x;
		double Ptemp = 0;
		int i;
		for(i = 2; i <= l; i++){
			Ptemp = ((2.0*i-1.0)*x*Plminus1 - (i-1.0)*Plminus2)/i;
			Plminus2 = Plminus1;
			Plminus1 = Ptemp;
		}
		return Ptemp;
	}
}

/** Computes the value of the 2nd Legendre function P_{l}(x) for given l and x
 *
 * @param l the spherical harmonic l parameter
 * @param x the x value
 * @return the value of the 2nd Legendre function
 */
double legendre_Q(int l, double x)
{
	gsl_sf_result result;
	int status = gsl_sf_legendre_Ql_e (l, x, &result);
	if(status != GSL_SUCCESS){
		printf("\nError computing LegendreQ\n");
		exit(0);
	}
	if(result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for LegendreQ not reached.\n");
		printf("Error in LegendreQ: %.20lf\n\n", result.err);
		exit(0);
	}
	return result.val;
}

/** Computes the value of the Polygamma function Psi(n,x) for given n and x
 * 
 * @param n the polygamma n
 * @param x the polygamma x
 * @return the value of the Polygamma function
 */
double polygamma(int n, double x)
{
	gsl_sf_result result;
	int status = gsl_sf_psi_n_e(n, x, &result);
	int gsl_sf_psi_n_e (int n, double x, gsl_sf_result * result);
	if(status != GSL_SUCCESS){
		printf("\nError computing Polygamma\n");
		exit(0);
	}
	if(result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for polygamma not reached.\n");
		printf("Error in Polygamma: %.20lf\n\n", result.err);
		exit(0);
	}
	return result.val;
}

/** Computes the value of the LambertW for given x. 
 *
 * The LambertW function is the inverse of \f$f(x) = x*exp(x).\f$
 * @param x as defined above
 * @return the value of the LambertW function
 */
double lambertW(double x){
	gsl_sf_result result;
	int status = gsl_sf_lambert_W0_e(x, &result);
	if(status != GSL_SUCCESS){
		printf("\nError computing LambertW\n");
		exit(0);
	}
	if(result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for LambertW not reached.\n");
		printf("Error in LambertW: %.20lf\n\n", result.err);
		exit(0);
	}
	return result.val;
}

/** Computes the value of the Clebsch-Gordan coefficient
 *
 * Since the arguments of the standard coupling coefficient functions are integer or half-integer, the arguments of the following function is, by convention, integers equal to twice the actual spin value. They are \f$(-1)^{j_1-j_2+m_3}\sqrt{2j_3+1} \f$ times the Wigner 3j symbols.
 * @param twice_j1 twice the stand j1 coupling coefficient
 * @param twice_j2 twice the stand j2 coupling coefficient
 * @param twice_j3 twice the stand j3 coupling coefficient
 * @param twice_m1 twice the stand m1 coupling coefficient
 * @param twice_m2 twice the stand m2 coupling coefficient
 * @param twice_m3 twice the stand m3 coupling coefficient
 * @return the value of the Clebsch-Gordan coefficient
 */
double clebsch_gordan_coefficient(int twice_j1, int twice_j2, int twice_j3, int twice_m1, int twice_m2, int twice_m3)
{
	gsl_sf_result result;
	int status = gsl_sf_coupling_3j_e(twice_j1, twice_j2, twice_j3, twice_m1, twice_m2, -twice_m3, &result);
	if(status != GSL_SUCCESS){
		printf("\nError computing Clebsch-Gordan Coefficient\n");
		exit(0);
	}
	if(result.err/result.val > DESIRED_FUNCTION_ACCURACY){
		printf("\nDesired function accuracy for Clebsch-Gordan Coefficient not reached.\n");
		printf("Error in Clebsch-Gordan Coefficient: %.20lf\n\n", result.err);
		exit(0);
	}
	return pow(-1.0, twice_m3/2.0 + twice_j1/2.0 - twice_j2/2.0) * sqrt(twice_j3 + 1.0) * result.val;
}

















