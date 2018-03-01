/** \file scaling_coefficients.c
 * Calculates the coefficients required to correctly scale the homogeneous fields to get the inhomogeneous solution for that n mode. The rescaling is not done here but in the main file.
 */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include "scaling_coefficients.h"
#include "sources.h"
#include "Runge-Kutta.h"
#include "convergence_tests.h"
#include "error_estimate.h"

extern const double MATCHING_INTEGRATION_ACCURACY;
extern const int convergence_test_type;
extern const int abs_error_threshold_fix;

#define RE_IN_FIELD 	0
#define IM_IN_FIELD		1
#define RE_OUT_FIELD 	2
#define IM_OUT_FIELD	3

struct integrand_params
{
	int field_piece;
	int field_num;
	struct n_mode_data *n_mode;
	struct orbital_params *orbit;
	
};

double R_times_source(double chi, int field_piece, int field_num, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	struct coupled_set *cset = n_mode->cset;

	//Form the matrix system
	int i, j; 
	gsl_matrix_complex * m = gsl_matrix_complex_alloc (2*cset->num_coupled_fields, 2*cset->num_coupled_fields);

	for (i = 0; i < cset->num_coupled_fields; i++){
		for (j = 0; j < cset->num_coupled_fields; j++){
			gsl_complex minus_in_field 				= gsl_complex_rect(-creal(n_mode->R_in[i][j][0]), -cimag(n_mode->R_in[i][j][0]));
			gsl_complex minus_in_field_rs_deriv 	= gsl_complex_rect(-creal(n_mode->R_in_rs_deriv[i][j][0]), -cimag(n_mode->R_in_rs_deriv[i][j][0]));
			gsl_complex out_field 					= gsl_complex_rect(creal(n_mode->R_out[i][j][0]), cimag(n_mode->R_out[i][j][0]));
			gsl_complex out_field_rs_deriv 			= gsl_complex_rect(creal(n_mode->R_out_rs_deriv[i][j][0]), cimag(n_mode->R_out_rs_deriv[i][j][0]));

			gsl_matrix_complex_set (m, i, j, minus_in_field);
			gsl_matrix_complex_set (m, i, j + cset->num_coupled_fields, out_field);
			gsl_matrix_complex_set (m, i + cset->num_coupled_fields, j, minus_in_field_rs_deriv);
			gsl_matrix_complex_set (m, i + cset->num_coupled_fields, j + cset->num_coupled_fields, out_field_rs_deriv);
		}
	}

	// check the condition number
/*	if(chi > 3.14 && output_once == 0){
		double condition_number = complex_matrix_condition_number(m, 2*cset->num_coupled_fields, 2*cset->num_coupled_fields);
		printf("condition number = %e\t", condition_number);
		output_once = 1;
	}*/
	

	// Form the source vector
	gsl_vector_complex * J = gsl_vector_complex_alloc(2*cset->num_coupled_fields);
	for( i = 0; i < cset->num_coupled_fields; i++){
		gsl_complex zero 	= gsl_complex_rect(0,0);
		gsl_complex source	= gsl_complex_rect(creal(cset->source[i](chi, n_mode, orbit, SOURCE_TOTAL)), cimag(cset->source[i](chi, n_mode, orbit, SOURCE_TOTAL)));
		gsl_vector_complex_set(J, i, zero);
		gsl_vector_complex_set(J, i+cset->num_coupled_fields, source); 
	}



	// Invert the system to get Phi^{-1} * J = C
	gsl_vector_complex *scaling_coeffs = gsl_vector_complex_alloc (2*cset->num_coupled_fields);

	int s;

	gsl_permutation * p = gsl_permutation_alloc (2*cset->num_coupled_fields);
	gsl_linalg_complex_LU_decomp (m, p, &s);
	gsl_linalg_complex_LU_solve (m, p, J, scaling_coeffs);	



	double return_result = 0;

	if(field_piece == RE_IN_FIELD) return_result 		= GSL_REAL(gsl_vector_complex_get(scaling_coeffs, field_num));
	else if(field_piece == IM_IN_FIELD) return_result 	= GSL_IMAG(gsl_vector_complex_get(scaling_coeffs, field_num));
	else if(field_piece == RE_OUT_FIELD) return_result  = GSL_REAL(gsl_vector_complex_get(scaling_coeffs, cset->num_coupled_fields + field_num));
	else if(field_piece == IM_OUT_FIELD) return_result 	= GSL_IMAG(gsl_vector_complex_get(scaling_coeffs, cset->num_coupled_fields + field_num));

	//if(chi < 0.01) printf("%e %e\n", chi, return_result);

	gsl_vector_complex_free(J);
	gsl_permutation_free(p);
	gsl_vector_complex_free (scaling_coeffs);
	gsl_matrix_complex_free (m);

	return return_result;

}


double scaling_integrand(double chi, void *params){

	struct integrand_params* params2 = (struct integrand_params *)params;

	double f_p = f(r_of_chi(chi));

	return R_times_source(chi, params2->field_piece, params2->field_num, params2->n_mode, params2->orbit) * (1.0/(2.0*M_PI))  / f_p;
}



double integrate_to_find_scaling_coeffs(int field_piece, int field_num, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	struct integrand_params params = {field_piece, field_num, n_mode, orbit};

	return M_PI*scaling_integrand(M_PI/2.0, &params);

}


void calculate_scaling_coefficients(struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	int i;
	for(i = 0; i < n_mode->cset->num_coupled_fields; i++){
		n_mode->C_in[i] 	= integrate_to_find_scaling_coeffs(RE_IN_FIELD, i, n_mode, orbit) + I * integrate_to_find_scaling_coeffs(IM_IN_FIELD, i, n_mode, orbit);
		n_mode->C_out[i]	= integrate_to_find_scaling_coeffs(RE_OUT_FIELD, i, n_mode, orbit) + I * integrate_to_find_scaling_coeffs(IM_OUT_FIELD, i, n_mode, orbit);
	}
    

}

