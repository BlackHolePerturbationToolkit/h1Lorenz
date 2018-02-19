#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include <complex.h>
#include "GSF_ecc.h"
#include "sources.h"

#define RE_IN_FIELD 	0
#define IM_IN_FIELD		1
#define RE_OUT_FIELD 	2
#define IM_OUT_FIELD	3

double estimate_integrand(int field_piece, int field_num, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	struct coupled_set *cset = n_mode->cset;


	//Form the matrix system at periastron
	int i, j; 
	gsl_matrix_complex * m_periastron 			= gsl_matrix_complex_alloc (2*cset->num_coupled_fields, 2*cset->num_coupled_fields);
	gsl_matrix_complex * m_periastron_inverse 	= gsl_matrix_complex_alloc (2*cset->num_coupled_fields, 2*cset->num_coupled_fields);

	for (i = 0; i < cset->num_coupled_fields; i++){
		for (j = 0; j < cset->num_coupled_fields; j++){
			gsl_complex minus_in_field 				= gsl_complex_rect(-creal(n_mode->R_in[i][j][0]), -cimag(n_mode->R_in[i][j][0]));
			gsl_complex minus_in_field_rs_deriv 	= gsl_complex_rect(-creal(n_mode->R_in_rs_deriv[i][j][0]), -cimag(n_mode->R_in_rs_deriv[i][j][0]));
			gsl_complex out_field 					= gsl_complex_rect(creal(n_mode->R_out[i][j][0]), cimag(n_mode->R_out[i][j][0]));
			gsl_complex out_field_rs_deriv 			= gsl_complex_rect(creal(n_mode->R_out_rs_deriv[i][j][0]), cimag(n_mode->R_out_rs_deriv[i][j][0]));

			gsl_matrix_complex_set (m_periastron, i, j, minus_in_field);
			gsl_matrix_complex_set (m_periastron, i, j + cset->num_coupled_fields, out_field);
			gsl_matrix_complex_set (m_periastron, i + cset->num_coupled_fields, j, minus_in_field_rs_deriv);
			gsl_matrix_complex_set (m_periastron, i + cset->num_coupled_fields, j + cset->num_coupled_fields, out_field_rs_deriv);
		}
	}

	int s_periastron;

	gsl_permutation * p_periastron = gsl_permutation_alloc (2*cset->num_coupled_fields);
	gsl_linalg_complex_LU_decomp (m_periastron, p_periastron, &s_periastron);
	gsl_linalg_complex_LU_invert(m_periastron, p_periastron, m_periastron_inverse);	

	// Form the source vector
	gsl_vector_complex * J_periastron = gsl_vector_complex_alloc(2*cset->num_coupled_fields);
	for( i = 0; i < cset->num_coupled_fields; i++){
		gsl_complex zero 	= gsl_complex_rect(0,0);
		gsl_complex source	= gsl_complex_rect(creal(cset->source[i](orbit->chi_p[0], n_mode, orbit, SOURCE_COEFF)), cimag(cset->source[i](orbit->chi_p[0], n_mode, orbit, SOURCE_COEFF)));
		gsl_vector_complex_set(J_periastron, i, zero);
		gsl_vector_complex_set(J_periastron, i+cset->num_coupled_fields, source); 
	}

	//Form the matrix system at apastron
	gsl_matrix_complex * m_apastron 		= gsl_matrix_complex_alloc (2*cset->num_coupled_fields, 2*cset->num_coupled_fields);
	gsl_matrix_complex * m_apastron_inverse = gsl_matrix_complex_alloc (2*cset->num_coupled_fields, 2*cset->num_coupled_fields);

	for (i = 0; i < cset->num_coupled_fields; i++){
		for (j = 0; j < cset->num_coupled_fields; j++){
			gsl_complex minus_in_field 				= gsl_complex_rect(-creal(n_mode->R_in[i][j][NUM_CHI_VALUES]), -cimag(n_mode->R_in[i][j][NUM_CHI_VALUES]));
			gsl_complex minus_in_field_rs_deriv 	= gsl_complex_rect(-creal(n_mode->R_in_rs_deriv[i][j][NUM_CHI_VALUES]), -cimag(n_mode->R_in_rs_deriv[i][j][NUM_CHI_VALUES]));
			gsl_complex out_field 					= gsl_complex_rect(creal(n_mode->R_out[i][j][NUM_CHI_VALUES]), cimag(n_mode->R_out[i][j][NUM_CHI_VALUES]));
			gsl_complex out_field_rs_deriv 			= gsl_complex_rect(creal(n_mode->R_out_rs_deriv[i][j][NUM_CHI_VALUES]), cimag(n_mode->R_out_rs_deriv[i][j][NUM_CHI_VALUES]));

			gsl_matrix_complex_set (m_apastron, i, j, minus_in_field);
			gsl_matrix_complex_set (m_apastron, i, j + cset->num_coupled_fields, out_field);
			gsl_matrix_complex_set (m_apastron, i + cset->num_coupled_fields, j, minus_in_field_rs_deriv);
			gsl_matrix_complex_set (m_apastron, i + cset->num_coupled_fields, j + cset->num_coupled_fields, out_field_rs_deriv);
		}
	}

	// Form the source vector
	gsl_vector_complex * J_apastron = gsl_vector_complex_alloc(2*cset->num_coupled_fields);
	for( i = 0; i < cset->num_coupled_fields; i++){
		gsl_complex zero 	= gsl_complex_rect(0,0);
		gsl_complex source	= gsl_complex_rect(creal(cset->source[i](orbit->chi_p[NUM_CHI_VALUES], n_mode, orbit, SOURCE_COEFF)), cimag(cset->source[i](orbit->chi_p[NUM_CHI_VALUES], n_mode, orbit, SOURCE_COEFF)));
		gsl_vector_complex_set(J_apastron, i, zero);
		gsl_vector_complex_set(J_apastron, i+cset->num_coupled_fields, source); 
	}



	int s_apastron;

	gsl_permutation * p_apastron = gsl_permutation_alloc (2*cset->num_coupled_fields);
	gsl_linalg_complex_LU_decomp (m_apastron, p_apastron, &s_apastron);
	gsl_linalg_complex_LU_invert(m_apastron, p_apastron, m_apastron_inverse);	

	double apastron_field_contrib = 0;
	double periastron_field_contrib = 0;
	double fabs_deriv = 0;

	if(field_piece == RE_IN_FIELD){
		for(i=0; i < 2*cset->num_coupled_fields; i++){
			apastron_field_contrib 		+=	GSL_REAL(gsl_complex_mul(gsl_matrix_complex_get(m_apastron_inverse, field_num, i), gsl_vector_complex_get(J_apastron, i)));
			periastron_field_contrib 	+=	GSL_REAL(gsl_complex_mul(gsl_matrix_complex_get(m_periastron_inverse, field_num, i), gsl_vector_complex_get(J_periastron, i)));
		}
		fabs_deriv = fabs((apastron_field_contrib - periastron_field_contrib)/(orbit->rs_max - orbit->rs_min));
	}
	else if(field_piece == IM_IN_FIELD){
		for(i=0; i < 2*cset->num_coupled_fields; i++){
			apastron_field_contrib 		+=	GSL_IMAG(gsl_complex_mul(gsl_matrix_complex_get(m_apastron_inverse, field_num, i),gsl_vector_complex_get(J_apastron, i)));
			periastron_field_contrib 	+=	GSL_IMAG(gsl_complex_mul(gsl_matrix_complex_get(m_periastron_inverse, field_num, i), gsl_vector_complex_get(J_periastron, i)));
		}
		fabs_deriv = fabs((apastron_field_contrib - periastron_field_contrib)/(orbit->rs_max - orbit->rs_min));
	}
	else if(field_piece == RE_OUT_FIELD){
		for(i=0; i < 2*cset->num_coupled_fields; i++){
			apastron_field_contrib 		+=	GSL_REAL(gsl_complex_mul(gsl_matrix_complex_get(m_apastron_inverse, cset->num_coupled_fields + field_num, i), gsl_vector_complex_get(J_apastron,i)));
			periastron_field_contrib 	+=	GSL_REAL(gsl_complex_mul(gsl_matrix_complex_get(m_periastron_inverse, cset->num_coupled_fields + field_num, i), gsl_vector_complex_get(J_periastron, i)));
		}
		fabs_deriv = fabs((apastron_field_contrib - periastron_field_contrib)/(orbit->rs_max - orbit->rs_min));
	}
	else if(field_piece == IM_OUT_FIELD){
		for(i=0; i < 2*cset->num_coupled_fields; i++){
			apastron_field_contrib 		+=	GSL_IMAG(gsl_complex_mul(gsl_matrix_complex_get(m_apastron_inverse, cset->num_coupled_fields + field_num, i), gsl_vector_complex_get(J_apastron,i)));
			periastron_field_contrib 	+=	GSL_IMAG(gsl_complex_mul(gsl_matrix_complex_get(m_periastron_inverse, cset->num_coupled_fields + field_num, i), gsl_vector_complex_get(J_periastron, i)));
		}
		fabs_deriv = fabs((apastron_field_contrib - periastron_field_contrib)/(orbit->rs_max - orbit->rs_min));
	}

	gsl_permutation_free(p_apastron);
	gsl_permutation_free(p_periastron);
	gsl_matrix_complex_free (m_apastron);
	gsl_matrix_complex_free (m_periastron);
	gsl_matrix_complex_free (m_apastron_inverse);
	gsl_matrix_complex_free (m_periastron_inverse);
	gsl_vector_complex_free (J_periastron);
	gsl_vector_complex_free (J_apastron);

	
	if(n_mode->n !=0) return fabs_deriv/(n_mode->n * orbit->Omega_r);
	else return fabs_deriv;

}
