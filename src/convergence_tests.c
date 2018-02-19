#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "convergence_tests.h"
#include "GSF_ecc.h"

extern const double n_convergence_threshold;

extern const int convergence_test_type;

int inner_outer_diff_vs_jump(struct lm_mode_data *lm_mode, struct orbital_params *orbit, int n);
void threshold_outer_field_contributions(int *converged, struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit);

int test_index = 0;

void test_for_convergence(int *converged, struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit, int abs_n_max)
{
	if(n_mode->n >= abs_n_max){
		*converged = 1;
	}else{
		if(convergence_test_type == INNER_VS_OUTER){
			*converged = inner_outer_diff_vs_jump(lm_mode, orbit, n_mode->n);
		}else if(convergence_test_type == OUTER_CONTRIBUTIONS){
			threshold_outer_field_contributions(converged, lm_mode, n_mode, orbit);
		}
	}
}

/** The outer field converges much faster than the inner field, as doe the negative n modes. Therefore we can convergence test based on that.
* Of course this will only guarentee that the GSF calculated from the outer field is correct
*/
void threshold_outer_field_contributions(int *converged, struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	int k, d;

	double m = n_mode->m;
	double n = n_mode->n;

	double max_outer_field_contribution = 0;

	int test_index;
	for(test_index = 0; test_index < 2*NUM_CHI_VALUES+1; test_index++){

		int r_chi_index = test_index;
		if(test_index > NUM_CHI_VALUES) r_chi_index = 2*NUM_CHI_VALUES - test_index;

		for(k = 0; k < n_mode->cset->num_coupled_fields; k++){
			double complex outer_field = 0;
			for(d = 0; d < n_mode->cset->num_coupled_fields; d++){
				outer_field					 	+= n_fold_function(m, n, n_mode->C_out[d] * n_mode->R_out_rs_deriv[k][d][r_chi_index]*cexp(-I*n_mode->omega * orbit->t_p[test_index]));
			}
			if(cabs(outer_field) > max_outer_field_contribution) max_outer_field_contribution = cabs(outer_field);
		}
	}

	nprintf("\t outer field contributions: %e. Relative to field magnitude: %e", max_outer_field_contribution, cabs(lm_mode->h_out[0][0])*max_outer_field_contribution);

	// For m = 0 the inner and outer fields are complex conjugates of each other
	if(m == 0){
		if(max_outer_field_contribution < n_convergence_threshold && n_mode->n > 3) *converged = 1;
	}else{
		if(n < 0 && max_outer_field_contribution < n_convergence_threshold && n_mode->n < -3){
			*converged = 2;
			nprintf("\t-n terms have converged");	
		}
		if(n > 0 && *converged == N_MIN_CONVERGED && max_outer_field_contribution < n_convergence_threshold && n_mode->n > 3){
			*converged = 1;
		}
	}
}


/** Tests for convergence of the n sum by looking at the continuity of the fields and comparing the difference in their derivatives with the time domain jump from the delta function. There is also a maximum number of n modes set by n_max.
*
* @param lm_mode the lm mode being tested
* @param n the most recently calculated n mode index
*/
int inner_outer_diff_vs_jump(struct lm_mode_data *lm_mode, struct orbital_params *orbit, int n)
{
	double converged = 0;

	double E = orbit->E;
	double L = orbit->L;

	struct coupled_set *cset = lm_mode->cset;

	int d;
	nprintf("Convergence test:\t");
	for(d = 0; d < cset->num_coupled_fields; d++){
		nprintf("%.3e\t", fabs(creal(lm_mode->h_out[d][test_index] - lm_mode->h_in[d][test_index])) );
	}
	nprintf("|\t");

	for(d = 0; d < cset->num_gauge_fields; d++){
		nprintf("%.3e\t", fabs(creal(lm_mode->h_out[cset->num_coupled_fields+d][test_index] - lm_mode->h_in[cset->num_coupled_fields+d][test_index]) ));
	}

	int l = lm_mode->l;
	int m = lm_mode->m;

	int i;
	for(i = test_index; i <= test_index; i++){
		double r_p			= orbit->r_p[i];
		double phi_p		= orbit->phi_p[i];
		double contra_ur	= orbit_contra_ur(orbit->chi_p[i]);
		double contra_ut 	= orbit_contra_ut(r_p);
		double f_p			= f(r_p);

	if((l+m) % 2 != 0){

		double diff_R9 = (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[0][i] - lm_mode->h_in_rs_deriv[0][i]);
		double jump_R9 = cabs(8.0*M_PI*L/(E) * normalized_associated_legendre_P_theta_deriv(l, m, M_PI/2.0) * contra_ur / (r_p*r_p) * 4.0 * cexp(-I * m * phi_p));

		double diff_R10 = 1;
		double jump_R10 = 1;
		double R10_ratio = 0;

		double diff_R8 = (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[cset->num_coupled_fields][i] - lm_mode->h_in_rs_deriv[cset->num_coupled_fields][i]);
		double jump_R8 = cabs( -4.0 * 8.0 * M_PI * L / (r_p*r_p) *  normalized_associated_legendre_P_theta_deriv(l, m, M_PI/2.0) * cexp(-I*m*phi_p));

		double R8_ratio = fabs(1.0 - diff_R8/jump_R8);

		if( l > 1){
			diff_R10 	= (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[1][i] - lm_mode->h_in_rs_deriv[1][i]);
			jump_R10	= cabs( -8.0 * M_PI * I * m * pow(L,2.0)/ (pow(r_p,3.0) * E) * normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(-I*m*phi_p) * 4.0);

			if( m == 0){
				R10_ratio = cabs(lm_mode->h_out_rs_deriv[1][i] - lm_mode->h_in_rs_deriv[1][i]);
			}else{
				R10_ratio = fabs(1.0 - diff_R10/jump_R10);
			}
		}

		double R9_ratio = fabs(1.0-(diff_R9/jump_R9));

		if(isnan(R10_ratio) || isinf(R10_ratio)){
			R10_ratio = 0;
		}

		nprintf("\tchi = %lg ratio = %.3e %.3e %.3e", orbit->chi_p[i], R9_ratio, R10_ratio, R8_ratio);

		if( R9_ratio < n_convergence_threshold && R10_ratio < n_convergence_threshold && R8_ratio < n_convergence_threshold){
			converged = 1;
		}

	}else{

		double diff_R1 = fabs(1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[0][i] - lm_mode->h_in_rs_deriv[0][i]);
		double jump_R1 = cabs(-4.0 * 4.0 * M_PI /(E * pow(r_p,3.0)) * (2.0*pow(E,2.0)*r_p*r_p - f_p*r_p*r_p - pow(L,2.0)*f_p)*normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(-I *m *phi_p));

		double diff_R3 = fabs(1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[1][i] - lm_mode->h_in_rs_deriv[1][i]);
		double jump_R3 = cabs(-4.0 * 4.0 * M_PI / (E * pow(r_p,3.0)) * (r_p*r_p + pow(L,2.0)) * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(-I*m*phi_p));


		double diff_R5 = 1;
		double jump_R5 = cabs(4.0 * 8.0 * M_PI * I * contra_ur * m * L/(E*r_p*r_p) * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(-I*m*phi_p));

		double diff_R6 = 1; 
		double jump_R6 = cabs(-4.0 * 4.0 * M_PI * pow(L, 2.0) / (E*pow(r_p,3.0)) * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(-I*m*phi_p));

		double diff_R7 = 1;
		double jump_R7 = fabs((l*(l+1.0) - 2.0*m*m)*jump_R6);

		double diff_R2 = 1;
		double jump_R2 = cabs(4.0 * 8.0 * M_PI * contra_ur/r_p * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(-I*m*phi_p));

		double diff_R4 = 1;
		double jump_R4 = cabs(-4.0 * 8.0 * M_PI * I * m * L/(r_p*r_p) * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(-I*m*phi_p) ) ;

		double R1_ratio = fabs(1.0-diff_R1/jump_R1);
		double R2_ratio = 0;
		double R3_ratio = fabs(1.0-diff_R3/jump_R3);
		double R4_ratio = 0;
		double R5_ratio = 0;
		double R6_ratio = 0;
		double R7_ratio = 0;

		if( l == 0){
			diff_R6 = (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[2][i] - lm_mode->h_in_rs_deriv[2][i]);
			diff_R2 = (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[3][i] - lm_mode->h_in_rs_deriv[3][i]);
	
			R6_ratio = fabs(1.0-diff_R6/jump_R6);
			R2_ratio = fabs(1.0-diff_R2/jump_R2);
		}

		if( l == 1) {
			diff_R5 = (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[2][i] - lm_mode->h_in_rs_deriv[2][i]);
			diff_R6 = (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[3][i] - lm_mode->h_in_rs_deriv[3][i]);
			diff_R2 = (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[4][i] - lm_mode->h_in_rs_deriv[4][i]);
			diff_R4 = (1.0 - pow(contra_ur/(f_p*contra_ut),2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[5][i] - lm_mode->h_in_rs_deriv[5][i]);

			R5_ratio = fabs(1.0-diff_R5/jump_R5);
			R6_ratio = fabs(1.0-diff_R6/jump_R6);
			R2_ratio = fabs(1.0-diff_R2/jump_R2);
			R4_ratio = fabs(1.0-diff_R4/jump_R4);

			if(m == 0){
				R5_ratio = cabs(lm_mode->h_out_rs_deriv[2][i] - lm_mode->h_in_rs_deriv[2][i]);
				R4_ratio = cabs(lm_mode->h_out_rs_deriv[5][i] - lm_mode->h_in_rs_deriv[5][i]);
			}

		}

		if(l >= 2){				
			diff_R5 = (1.0 - pow(contra_ur/(f_p*contra_ut), 2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[2][i] - lm_mode->h_in_rs_deriv[2][i]);
			diff_R6 = (1.0 - pow(contra_ur/(f_p*contra_ut), 2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[3][i] - lm_mode->h_in_rs_deriv[3][i]);
			diff_R7 = (1.0 - pow(contra_ur/(f_p*contra_ut), 2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[4][i] - lm_mode->h_in_rs_deriv[4][i]);
			diff_R2 = (1.0 - pow(contra_ur/(f_p*contra_ut), 2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[5][i] - lm_mode->h_in_rs_deriv[5][i]);
			diff_R4 = (1.0 - pow(contra_ur/(f_p*contra_ut), 2.0)) * 1.0/f_p * cabs(lm_mode->h_out_rs_deriv[6][i] - lm_mode->h_in_rs_deriv[6][i]);

			R5_ratio = fabs(1.0-diff_R5/jump_R5);
			R6_ratio = fabs(1.0-diff_R6/jump_R6);
			R7_ratio = fabs(1.0-diff_R7/jump_R7);
			R2_ratio = fabs(1.0-diff_R2/jump_R2);
			R4_ratio = fabs(1.0-diff_R4/jump_R4);

			if(m == 0){
				R5_ratio = cabs(lm_mode->h_out_rs_deriv[2][i] - lm_mode->h_in_rs_deriv[2][i]);
				R4_ratio = cabs(lm_mode->h_out_rs_deriv[6][i] - lm_mode->h_in_rs_deriv[6][i]);
			}

		}

		if(isnan(R4_ratio) || isinf(R4_ratio)){
			R4_ratio = 0;
		}

		if(isnan(R5_ratio) || isinf(R5_ratio)){
			R5_ratio = 0;
		}

		if(isnan(R7_ratio) || isinf(R7_ratio)){
			R7_ratio = 0;
		}
		

		nprintf("chi = %lg ratio = %.3e %.3e %.3e %.3e %.3e %.3e %.3e ", orbit->chi_p[i], R1_ratio, R3_ratio, R5_ratio, R6_ratio, R7_ratio, R2_ratio, R4_ratio);

		if(R1_ratio < n_convergence_threshold && R3_ratio < n_convergence_threshold && R5_ratio < n_convergence_threshold && R6_ratio < n_convergence_threshold && R7_ratio < n_convergence_threshold && R2_ratio< n_convergence_threshold && R4_ratio < n_convergence_threshold){
			converged = 1;
		}

		}
	}

	return converged;
}
