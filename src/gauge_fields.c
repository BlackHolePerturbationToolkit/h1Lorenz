#include "gauge_fields.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <stdio.h>

#define IN_RE 0
#define IN_IM 1
#define OUT_RE 2
#define OUT_IM 3

#define MATCHING_INTEGRATION_ACCURACY 1e-10

struct static_integrand_params{

	struct n_mode_data *n_mode;
	int piece;
};

double R8_static_matching_integrand(double chi, void *params)
{
	struct static_integrand_params *params2 = (struct static_integrand_params *)params;
	
	int piece = params2->piece;

	int l				= params2->n_mode->l;

	double r_p 			= r_of_chi(chi);
	double f_p 			= f(r_p);

	double R8_out 				= R8_static_out_solution(l, r_p);
	double R8_out_deriv			= R8_static_out_rs_deriv_solution(l, r_p);
	double R8_in				= R8_static_in_solution(l, r_p);
	double R8_in_deriv			= R8_static_in_rs_deriv_solution(l, r_p);

	double wronskian 			= R8_in*R8_out_deriv - R8_out*R8_in_deriv ;

	if(piece == IN_RE) return R8_out*f_p/(wronskian*r_p*r_p)*(1.0/(2.0*M_PI));
	else  return R8_in*f_p/(wronskian*r_p*r_p)*(1.0/(2.0*M_PI));

}

double R8_integrate_static_component(struct orbital_params *orbit, struct n_mode_data *n_mode, int piece)
{
	struct static_integrand_params params = {n_mode, piece};

	double result, error = 0.0;

	gsl_function F;
	F.function = &R8_static_matching_integrand;
	F.params = &params;

	if(orbit->e == 0) return M_PI*R8_static_matching_integrand(M_PI/2.0, &params);

	gsl_integration_workspace *w = gsl_integration_workspace_alloc (100000);

	double return_code = gsl_integration_qag (&F, 0.0, M_PI, MATCHING_INTEGRATION_ACCURACY, MATCHING_INTEGRATION_ACCURACY, 100000, GSL_INTEG_GAUSS61, w, &result, &error); 
	if( return_code == GSL_EROUND){

	}

	gsl_integration_workspace_free (w);
	
	return result;

}


void construct_R8_gauge_field(struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit, int num_coupled_fields){

	double l 		= n_mode->l;
	double omega 	= n_mode->omega;

	// Variables used for the static case
	double c_out = 0;
	double c_in = 0;

	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES+1; i++){

		int r_chi_index = i;
		if(i > NUM_CHI_VALUES) r_chi_index = 2*NUM_CHI_VALUES - i;

		double t_p = orbit->t_p[i];
		double r_p = orbit->r_p[i];
		double f_p = f(r_p);
		double df = 2.0/(r_p*r_p);


		if( n_mode->m == 0 && n_mode->n == 0){
			double R8_coefficient 	= -4.0*16.0*M_PI*orbit->L * normalized_associated_legendre_P_theta_deriv(n_mode->l, n_mode->m, M_PI/2.0);
			
			// Only compute these once
			if(i==0){
				c_out 			= R8_coefficient*R8_integrate_static_component(orbit, n_mode, OUT_RE);
				c_in			= R8_coefficient*R8_integrate_static_component(orbit, n_mode, IN_RE);
			}

			lm_mode->h_out[num_coupled_fields][i] 			= c_out*R8_static_out_solution(l,r_p);
			lm_mode->h_out_rs_deriv[num_coupled_fields][i]	= c_out*R8_static_out_rs_deriv_solution(l,r_p);
			lm_mode->h_in[num_coupled_fields][i]			= c_in*R8_static_in_solution(l,r_p);
			lm_mode->h_in_rs_deriv[num_coupled_fields][i]	= c_in*R8_static_in_rs_deriv_solution(l,r_p);

			// The following two are zero for static modes
			lm_mode->h_out_t_deriv[num_coupled_fields][i] 	= 0;
			lm_mode->h_in_t_deriv[num_coupled_fields][i]	= 0;

			lm_mode->h_out[0][i] 			= 0;
			lm_mode->h_out_rs_deriv[0][i]	= 0;
			lm_mode->h_in[0][i]				= 0;
			lm_mode->h_in_rs_deriv[0][i]	= 0;
			lm_mode->h_out_t_deriv[0][i] 	= 0;
			lm_mode->h_in_t_deriv[0][i]		= 0;

			if(l>1){
				lm_mode->h_out[1][i] 			= 0;
				lm_mode->h_out_rs_deriv[1][i]	= 0;
				lm_mode->h_in[1][i]				= 0;
				lm_mode->h_in_rs_deriv[1][i]	= 0;
				lm_mode->h_out_t_deriv[1][i] 	= 0;
				lm_mode->h_in_t_deriv[1][i]		= 0;
			}
		}else{

			double complex R9_out 			= 0;
			double complex R9_in 			= 0;
			double complex R9_in_r_deriv 	= 0;
			double complex R9_out_r_deriv 	= 0;
			double complex R9_out_rr_deriv = 0;
			double complex R9_in_rr_deriv  = 0;


			double complex R10_out 		= 0;
			double complex R10_in 			= 0;
			double complex R10_out_r_deriv = 0;
			double complex R10_in_r_deriv 	= 0;

			double V_l 		= 0.25*f_p*(df/r_p + l*(l+1.0)/(r_p*r_p));

			int d;
			for(d = 0; d < num_coupled_fields; d++){
				R9_in				 		+= n_mode->C_in[d]  * n_mode->R_in[0][d][r_chi_index];
				R9_out						+= n_mode->C_out[d] * n_mode->R_out[0][d][r_chi_index];
				R9_in_r_deriv				+= 1.0/f_p * (n_mode->C_in[d] * n_mode->R_in_rs_deriv[0][d][r_chi_index]);
				R9_out_r_deriv				+= 1.0/f_p * (n_mode->C_out[d] * n_mode->R_out_rs_deriv[0][d][r_chi_index]);	

				if(num_coupled_fields != 1){
					R10_in				 	+= n_mode->C_in[d]  * n_mode->R_in[1][d][r_chi_index];
					R10_out					+= n_mode->C_out[d] * n_mode->R_out[1][d][r_chi_index];
					R10_in_r_deriv			+= 1.0/f_p * n_mode->C_in[d]  * n_mode->R_in_rs_deriv[1][d][r_chi_index];
					R10_out_r_deriv			+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[1][d][r_chi_index];
				}	

			}
			
			R9_in_rr_deriv				+= 1.0/(f_p*f_p) *( ( 4.0*V_l + 4.0*f_p/(r_p*r_p) * (1.0-4.5/r_p) - omega*omega) * R9_in - 2.0*f_p/(r_p*r_p)*(1.0-3.0/r_p)*R10_in - f_p*df*R9_in_r_deriv);
			R9_out_rr_deriv				+= 1.0/(f_p*f_p) *( ( 4.0*V_l + 4.0*f_p/(r_p*r_p) * (1.0-4.5/r_p) - omega*omega) * R9_out - 2.0*f_p/(r_p*r_p)*(1.0-3.0/r_p)*R10_out - f_p*df*R9_out_r_deriv);

			double m = n_mode->m;
			double n = n_mode->n;

			// Construct the h8 field from the gauge condition
			lm_mode->h_out[num_coupled_fields][i] 			+= n_fold_function(m, n, I*f_p/omega *  ( R9_out_r_deriv + 2.0 * R9_out/r_p - R10_out/r_p)*cexp(-I*omega*t_p)) ;
			lm_mode->h_in[num_coupled_fields][i]			+= n_fold_function(m, n, I*f_p/omega * ( R9_in_r_deriv + 2.0 * R9_in/r_p - R10_in/r_p)*cexp(-I*omega*t_p)) ;

			lm_mode->h_out_t_deriv[num_coupled_fields][i] 	+= n_fold_function(m, n, f_p * ( R9_out_r_deriv + 2.0 * R9_out/r_p - R10_out/r_p)*cexp(-I*omega*t_p)) ;
			lm_mode->h_in_t_deriv[num_coupled_fields][i]	+= n_fold_function(m, n, f_p * ( R9_in_r_deriv + 2.0 * R9_in/r_p - R10_in/r_p)*cexp(-I*omega*t_p)) ;

			lm_mode->h_out_rs_deriv[num_coupled_fields][i] 	+= n_fold_function(m, n, f_p * ( I/omega * ( f_p * R9_out_rr_deriv + df * R9_out_r_deriv + 2.0*f_p * (R9_out_r_deriv/r_p - R9_out/(r_p*r_p) ) + 2.0 * df * R9_out/r_p - f_p * (R10_out_r_deriv/r_p - R10_out/(r_p*r_p) ) - df * R10_out/r_p)*cexp(-I*omega*t_p) ));

			lm_mode->h_in_rs_deriv[num_coupled_fields][i] 	+= n_fold_function(m, n, f_p *(( I/omega * ( f_p * R9_in_rr_deriv + df * R9_in_r_deriv + 2.0*f_p * (R9_in_r_deriv/r_p - R9_in/(r_p*r_p) ) + 2.0 * df * R9_in/r_p - f_p * (R10_in_r_deriv/r_p - R10_in/(r_p*r_p) ) - df * R10_in/r_p))*cexp(-I*omega*t_p) ));

		}
	}


}

void construct_R2_R4_gauge_fields(struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit, int num_coupled_fields)
{
	double l 		= n_mode->l;
	double omega 	= n_mode->omega;

	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES+1; i++){

		int r_chi_index = i;
		if(i > NUM_CHI_VALUES) r_chi_index = 2*NUM_CHI_VALUES - i;

		double t_p = orbit->t_p[i];
		double r_p = orbit->r_p[i];
		double f_p = f(r_p);
		double df = 2.0/(r_p*r_p);

		double complex R1_out				= 0;
		double complex R1_out_r_deriv 		= 0;
		double complex R1_out_rr_deriv		= 0;
		double complex R1_in				= 0;
		double complex R1_in_r_deriv		= 0;
		double complex R1_in_rr_deriv		= 0;

		double complex R3_out			= 0;
		double complex R3_out_r_deriv 		= 0;
		double complex R3_out_rr_deriv		= 0;
		double complex R3_in			= 0;
		double complex R3_in_r_deriv		= 0;
		double complex R3_in_rr_deriv		= 0;

		double complex R5_out			= 0;
		double complex R5_out_r_deriv 		= 0;
		double complex R5_out_rr_deriv		= 0;
		double complex R5_in			= 0;
		double complex R5_in_r_deriv		= 0;
		double complex R5_in_rr_deriv		= 0;
		
		double complex R6_out			= 0;
		double complex R6_in			= 0;
		double complex R6_out_r_deriv		= 0;
		double complex R6_in_r_deriv		= 0;

		double complex R7_out			= 0;
		double complex R7_in			= 0;
		double complex R7_out_r_deriv		= 0;
		double complex R7_in_r_deriv		= 0;

		if(omega != 0){
			int d;
			for(d = 0; d < num_coupled_fields; d++){
				R1_in				 		+= n_mode->C_in[d]  * n_mode->R_in[0][d][r_chi_index];
				R1_out						+= n_mode->C_out[d] * n_mode->R_out[0][d][r_chi_index];
				R1_in_r_deriv				+= 1.0/f_p * n_mode->C_in[d] * n_mode->R_in_rs_deriv[0][d][r_chi_index];
				R1_out_r_deriv				+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[0][d][r_chi_index];	

				R3_in				 		+= n_mode->C_in[d]  * n_mode->R_in[1][d][r_chi_index];
				R3_out						+= n_mode->C_out[d] * n_mode->R_out[1][d][r_chi_index];
				R3_in_r_deriv				+= 1.0/f_p * n_mode->C_in[d]  * n_mode->R_in_rs_deriv[1][d][r_chi_index];
				R3_out_r_deriv				+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[1][d][r_chi_index];

				if(l == 0){
					R6_in				 	+= n_mode->C_in[d]  * n_mode->R_in[2][d][r_chi_index];
					R6_out					+= n_mode->C_out[d] * n_mode->R_out[2][d][r_chi_index];
					R6_in_r_deriv			+= 1.0/f_p * n_mode->C_in[d] * n_mode->R_in_rs_deriv[2][d][r_chi_index];
					R6_out_r_deriv			+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[2][d][r_chi_index];
				}else{
					R5_in				 	+= n_mode->C_in[d]  * n_mode->R_in[2][d][r_chi_index];
					R5_out					+= n_mode->C_out[d] * n_mode->R_out[2][d][r_chi_index];
					R5_in_r_deriv			+= 1.0/f_p * n_mode->C_in[d]  * n_mode->R_in_rs_deriv[2][d][r_chi_index];
					R5_out_r_deriv			+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[2][d][r_chi_index];

					R6_in				 	+= n_mode->C_in[d]  * n_mode->R_in[3][d][r_chi_index];
					R6_out					+= n_mode->C_out[d] * n_mode->R_out[3][d][r_chi_index];
					R6_in_r_deriv			+= 1.0/f_p * n_mode->C_in[d] * n_mode->R_in_rs_deriv[3][d][r_chi_index];
					R6_out_r_deriv			+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[3][d][r_chi_index];

					if(l >=2){
						R7_in				+= n_mode->C_in[d]  * n_mode->R_in[4][d][r_chi_index];
						R7_out				+= n_mode->C_out[d] * n_mode->R_out[4][d][r_chi_index];
						R7_in_r_deriv 		+= 1.0/f_p * n_mode->C_in[d] * n_mode->R_in_rs_deriv[4][d][r_chi_index];
						R7_out_r_deriv		+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[4][d][r_chi_index];
					}

				}
			}

			double V_l 				= 0.25*f_p*(df/r_p + l*(l+1.0)/(r_p*r_p)) - 0.25*omega*omega;

			R1_in_rr_deriv			= 	1.0/(f_p*f_p) * ( 4.0*V_l * R1_in + 4.0/(r_p*r_p) * f_p * f_p * R3_in_r_deriv + 2.0*f_p/(r_p*r_p) * (1.0 - 4.0/r_p)*(R1_in - R5_in - f_p*R3_in) - 2.0*f_p*f_p/(r_p*r_p)*(1.0-6.0/r_p)*R6_in  - f_p*df*R1_in_r_deriv);

			R1_out_rr_deriv			= 	1.0/(f_p*f_p) * ( 4.0*V_l * R1_out + 4.0/(r_p*r_p) * f_p * f_p * R3_out_r_deriv + 2.0*f_p/(r_p*r_p) * (1.0 - 4.0/r_p)*(R1_out - R5_out - f_p*R3_out) - 2.0*f_p*f_p/(r_p*r_p)*(1.0-6.0/r_p)*R6_out  - f_p*df*R1_out_r_deriv);

			R3_in_rr_deriv 			= 	1.0/(f_p*f_p) *  ( 4.0*V_l * R3_in  - 2.0*f_p/(r_p*r_p) * (R1_in  - R5_in  - (1.0-4.0/r_p)*(R3_in  + R6_in))   - f_p*df*R3_in_r_deriv);
			R3_out_rr_deriv 		= 	1.0/(f_p*f_p) *  ( 4.0*V_l * R3_out - 2.0*f_p/(r_p*r_p) * (R1_out - R5_out - (1.0-4.0/r_p)*(R3_out + R6_out))  - f_p*df*R3_out_r_deriv);

			R5_in_rr_deriv 			=	1.0/(f_p*f_p) * ( 4.0*V_l * R5_in + 4.0*f_p/(r_p*r_p)*((1.0-4.5/r_p)*R5_in - 0.5*l*(l+1.0)*(R1_in - f_p*R3_in) + 0.5*(1.0-3.0/r_p)*( l*(l+1.0)*R6_in - R7_in)) - f_p*df*R5_in_r_deriv);

			R5_out_rr_deriv 		=	1.0/(f_p*f_p) * ( 4.0*V_l * R5_out + 4.0*f_p/(r_p*r_p)*((1.0-4.5/r_p)*R5_out - 0.5*l*(l+1.0)*(R1_out - f_p*R3_out) + 0.5*(1.0-3.0/r_p)*( l*(l+1.0)*R6_out - R7_out)) - f_p*df*R5_out_r_deriv);

			double m = n_mode->m;
			double n = n_mode->n;

			// construct the h2 field
			lm_mode->h_out[num_coupled_fields][i] 			+= n_fold_function(m, n, I*f_p/omega * ( R1_out_r_deriv - f_p * R3_out_r_deriv + 1.0/r_p * (R1_out - R5_out - f_p*R3_out - 2.0*f_p*R6_out) )*cexp(-I*omega*t_p)) ;
			lm_mode->h_in[num_coupled_fields][i]			+= n_fold_function(m, n, I*f_p/omega * ( R1_in_r_deriv  - f_p * R3_in_r_deriv  + 1.0/r_p * (R1_in  - R5_in - f_p*R3_in  - 2.0*f_p*R6_in ) )*cexp(-I*omega*t_p)) ;

			lm_mode->h_out_t_deriv[num_coupled_fields][i]	+= n_fold_function(m, n, f_p * ( R1_out_r_deriv - f_p * R3_out_r_deriv + 1.0/r_p * (R1_out - R5_out - f_p*R3_out - 2.0*f_p*R6_out) )*cexp(-I*omega*t_p)) ;
			lm_mode->h_in_t_deriv[num_coupled_fields][i]	+= n_fold_function(m, n, f_p * ( R1_in_r_deriv  - f_p * R3_in_r_deriv  + 1.0/r_p * (R1_in  - R5_in - f_p*R3_in  - 2.0*f_p*R6_in ) )*cexp(-I*omega*t_p)) ;

			lm_mode->h_out_rs_deriv[num_coupled_fields][i] 	+= n_fold_function(m, n, f_p * (( I*df/omega * ( R1_out_r_deriv - f_p * R3_out_r_deriv + 1.0/r_p * (R1_out - R5_out - f_p*R3_out - 2.0*f_p*R6_out) ) + I*f_p/omega *(R1_out_rr_deriv - f_p*R3_out_rr_deriv - df*R3_out_r_deriv + R1_out_r_deriv/r_p - R1_out/(r_p*r_p) - R5_out_r_deriv/r_p + R5_out/(r_p*r_p) - f_p/r_p*R3_out_r_deriv - df/r_p*R3_out + f_p/(r_p*r_p)*R3_out - 2.0*f_p/r_p*R6_out_r_deriv - 2.0*df*R6_out/r_p + 2.0*f_p*R6_out/(r_p*r_p) ))*cexp(-I*omega*t_p) ));

			lm_mode->h_in_rs_deriv[num_coupled_fields][i] 	+= n_fold_function(m, n, f_p * (( I*df/omega * ( R1_in_r_deriv - f_p * R3_in_r_deriv + 1.0/r_p * (R1_in - R5_in - f_p*R3_in - 2.0*f_p*R6_in) ) + I*f_p/omega *(R1_in_rr_deriv - f_p*R3_in_rr_deriv - df*R3_in_r_deriv + R1_in_r_deriv/r_p - R1_in/(r_p*r_p) - R5_in_r_deriv/r_p + R5_in/(r_p*r_p) - f_p/r_p*R3_in_r_deriv - df/r_p*R3_in + f_p/(r_p*r_p)*R3_in - 2.0*f_p/r_p*R6_in_r_deriv - 2.0*df*R6_in/r_p + 2.0*f_p*R6_in/(r_p*r_p) ))*cexp(-I*omega*t_p) ));


			double complex R2_out			= (I*f_p/omega * ( R1_out_r_deriv - f_p * R3_out_r_deriv + 1.0/r_p * (R1_out - R5_out - f_p*R3_out - 2.0*f_p*R6_out) ));
			double complex R2_out_r_deriv 	= (I*df/omega * ( R1_out_r_deriv - f_p * R3_out_r_deriv + 1.0/r_p * (R1_out - R5_out - f_p*R3_out - 2.0*f_p*R6_out) ) + I*f_p/omega *(R1_out_rr_deriv - f_p*R3_out_rr_deriv - df*R3_out_r_deriv + R1_out_r_deriv/r_p - R1_out/(r_p*r_p) - R5_out_r_deriv/r_p + R5_out/(r_p*r_p) - f_p/r_p*R3_out_r_deriv - df/r_p*R3_out + f_p/(r_p*r_p)*R3_out - 2.0*f_p/r_p*R6_out_r_deriv - 2.0*df*R6_out/r_p + 2.0*f_p*R6_out/(r_p*r_p) ));
			double complex R4_out 			= (I*f_p/omega * (R5_out_r_deriv + 2.0*R5_out/r_p + l*(l+1.0) * R6_out/r_p - R7_out/r_p));

			double complex R2_in			= (I*f_p/omega * ( R1_in_r_deriv - f_p * R3_in_r_deriv + 1.0/r_p * (R1_in - R5_in - f_p*R3_in - 2.0*f_p*R6_in) ));
			double complex R2_in_r_deriv 	= (I*df/omega * ( R1_in_r_deriv - f_p * R3_in_r_deriv + 1.0/r_p * (R1_in - R5_in - f_p*R3_in - 2.0*f_p*R6_in) ) + I*f_p/omega *(R1_in_rr_deriv - f_p*R3_in_rr_deriv - df*R3_in_r_deriv + R1_in_r_deriv/r_p - R1_in/(r_p*r_p) - R5_in_r_deriv/r_p + R5_in/(r_p*r_p) - f_p/r_p*R3_in_r_deriv - df/r_p*R3_in + f_p/(r_p*r_p)*R3_in - 2.0*f_p/r_p*R6_in_r_deriv - 2.0*df*R6_in/r_p + 2.0*f_p*R6_in/(r_p*r_p) ));
			double complex R4_in 			= (I*f_p/omega * (R5_in_r_deriv + 2.0*R5_in/r_p + l*(l+1.0) * R6_in/r_p - R7_in/r_p));

			if(l >=1){
				//construct the h4 field
				lm_mode->h_out[num_coupled_fields+1][i] 			+= n_fold_function(m, n, I*f_p/omega * (R5_out_r_deriv + 2.0*R5_out/r_p + l*(l+1.0) * R6_out/r_p - R7_out/r_p)*cexp(-I*omega*t_p)) ;
				lm_mode->h_in[num_coupled_fields+1][i]				+= n_fold_function(m, n, I*f_p/omega * (R5_in_r_deriv  + 2.0*R5_in/r_p  + l*(l+1.0) * R6_in/r_p  - R7_in/r_p )*cexp(-I*omega*t_p)) ;	

				lm_mode->h_out_t_deriv[num_coupled_fields+1][i] 	+= n_fold_function(m, n, f_p * (R5_out_r_deriv + 2.0*R5_out/r_p + l*(l+1.0) * R6_out/r_p - R7_out/r_p)*cexp(-I*omega*t_p)) ;
				lm_mode->h_in_t_deriv[num_coupled_fields+1][i]		+= n_fold_function(m, n, f_p * (R5_in_r_deriv  + 2.0*R5_in/r_p  + l*(l+1.0) * R6_in/r_p  - R7_in/r_p )*cexp(-I*omega*t_p)) ;	
				

				lm_mode->h_out_rs_deriv[num_coupled_fields+1][i] 	+= n_fold_function(m, n, f_p * ( 	(I*df/omega * (R5_out_r_deriv + 2.0*R5_out/r_p + l*(l+1.0) * R6_out/r_p - R7_out/r_p)  + I*f_p/omega* (R5_out_rr_deriv + 2.0*R5_out_r_deriv/r_p - 2.0*R5_out/(r_p*r_p) + l*(l+1.0)*(R6_out_r_deriv/r_p - R6_out/(r_p*r_p)) - R7_out_r_deriv/r_p + R7_out/(r_p*r_p) ))*cexp(-I*omega*t_p) ));

				lm_mode->h_in_rs_deriv[num_coupled_fields+1][i] 	+= n_fold_function(m, n, f_p * ( ( 	I*df/omega * (R5_in_r_deriv + 2.0*R5_in/r_p + l*(l+1.0) * R6_in/r_p - R7_in/r_p)  + I*f_p/omega* (R5_in_rr_deriv + 2.0*R5_in_r_deriv/r_p - 2.0*R5_in/(r_p*r_p) + l*(l+1.0)*(R6_in_r_deriv/r_p - R6_in/(r_p*r_p)) - R7_in_r_deriv/r_p + R7_in/(r_p*r_p) ))*cexp(-I*omega*t_p) ));

			}

			//check gauge conditions
			// FIXME give worst value along the orbit
			if(i == 99) {
				double complex first_gauge_out		= I*omega*R1_out + f_p*(I*omega*R3_out 	+ R2_out_r_deriv 	+ R2_out/r_p 	- R4_out/r_p);
				double complex first_gauge_in		= I*omega*R1_in  + f_p*(I*omega*R3_in 	+ R2_in_r_deriv 	+ R2_in/r_p 	- R4_in/r_p);

				nprintf("1st gauge condition: %.3e %.3e ", cabs(first_gauge_out), cabs(first_gauge_in));

			}


		}


	}

	

}

void construct_R6_R7_R2_R4_gauge_fields(struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit, int num_coupled_fields)
{
	double l 		= n_mode->l;

	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES+1; i++){

		double r_p = orbit->r_p[i];
		double f_p = f(r_p);
		double df = 2.0/(r_p*r_p);

		double complex R1_out				= 0;
		double complex R1_out_r_deriv 		= 0;
		double complex R1_out_rr_deriv		= 0;
		double complex R1_in				= 0;
		double complex R1_in_r_deriv		= 0;
		double complex R1_in_rr_deriv		= 0;

		double complex R3_out				= 0;
		double complex R3_out_r_deriv 		= 0;
		double complex R3_out_rr_deriv		= 0;
		double complex R3_in				= 0;
		double complex R3_in_r_deriv		= 0;
		double complex R3_in_rr_deriv		= 0;

		double complex R5_out				= 0;
		double complex R5_out_r_deriv 		= 0;
		double complex R5_out_rr_deriv		= 0;
		double complex R5_in				= 0;
		double complex R5_in_r_deriv		= 0;
		double complex R5_in_rr_deriv		= 0;

		int r_chi_index = i;
		if(i > NUM_CHI_VALUES) r_chi_index = 2*NUM_CHI_VALUES - i;
		
		int d;
		for(d = 0; d < num_coupled_fields; d++){
			R1_in				 		+= n_mode->C_in[d]  * n_mode->R_in[0][d][r_chi_index];
			R1_out						+= n_mode->C_out[d] * n_mode->R_out[0][d][r_chi_index];
			R1_in_r_deriv				+= 1.0/f_p * n_mode->C_in[d] * n_mode->R_in_rs_deriv[0][d][r_chi_index];
			R1_out_r_deriv				+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[0][d][r_chi_index];	

			R3_in				 		+= n_mode->C_in[d]  * n_mode->R_in[1][d][r_chi_index];
			R3_out						+= n_mode->C_out[d] * n_mode->R_out[1][d][r_chi_index];
			R3_in_r_deriv				+= 1.0/f_p * n_mode->C_in[d]  * n_mode->R_in_rs_deriv[1][d][r_chi_index];
			R3_out_r_deriv				+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[1][d][r_chi_index];

			R5_in				 		+= n_mode->C_in[d]  * n_mode->R_in[2][d][r_chi_index];
			R5_out						+= n_mode->C_out[d] * n_mode->R_out[2][d][r_chi_index];
			R5_in_r_deriv				+= 1.0/f_p * n_mode->C_in[d]  * n_mode->R_in_rs_deriv[2][d][r_chi_index];
			R5_out_r_deriv				+= 1.0/f_p * n_mode->C_out[d] * n_mode->R_out_rs_deriv[2][d][r_chi_index];
		}

		double V_l 		= 0.25*f_p*(df/r_p + l*(l+1.0)/(r_p*r_p));

		// Construct the h6 field and its rs and t derivatives
		// No need to do n folding as this is only used in the n=0 mode. This mode is real.
		lm_mode->h_out[num_coupled_fields][i] 	+= 1.0/(2.0*f_p)*(r_p * R1_out_r_deriv - r_p*f_p*R3_out_r_deriv + R1_out - R5_out - f_p*R3_out);
		lm_mode->h_in[num_coupled_fields][i]	+= 1.0/(2.0*f_p)*(r_p * R1_in_r_deriv - r_p*f_p*R3_in_r_deriv + R1_in - R5_in - f_p*R3_in);

		double complex R6_out 	= lm_mode->h_out[num_coupled_fields][i];
		double complex R6_in	= lm_mode->h_in[num_coupled_fields][i];

		R1_in_rr_deriv		= 	1.0/(f_p*f_p) * ( 4.0*V_l * R1_in + 4.0/(r_p*r_p) * f_p * f_p * R3_in_r_deriv + 2.0*f_p/(r_p*r_p) * (1.0 - 4.0/r_p)*(R1_in - R5_in - f_p*R3_in) - 2.0*f_p*f_p/(r_p*r_p)*(1.0-6.0/r_p)*R6_in  - f_p*df*R1_in_r_deriv);

		R1_out_rr_deriv		= 	1.0/(f_p*f_p) * ( 4.0*V_l * R1_out + 4.0/(r_p*r_p) * f_p * f_p * R3_out_r_deriv + 2.0*f_p/(r_p*r_p) * (1.0 - 4.0/r_p)*(R1_out - R5_out - f_p*R3_out) - 2.0*f_p*f_p/(r_p*r_p)*(1.0-6.0/r_p)*R6_out  - f_p*df*R1_out_r_deriv);

		R3_in_rr_deriv 		= 	1.0/(f_p*f_p) *  ( 4.0*V_l * R3_in  - 2.0*f_p/(r_p*r_p) * (R1_in  - R5_in  - (1.0-4.0/r_p)*(R3_in  + R6_in))   - f_p*df*R3_in_r_deriv);
		R3_out_rr_deriv 	= 	1.0/(f_p*f_p) *  ( 4.0*V_l * R3_out - 2.0*f_p/(r_p*r_p) * (R1_out - R5_out - (1.0-4.0/r_p)*(R3_out + R6_out))  - f_p*df*R3_out_r_deriv);

		lm_mode->h_out_rs_deriv[num_coupled_fields][i] 	+=	-1.0/(2.0*f_p) * (R1_out*df - R5_out*df - 2.0*f_p*R1_out_r_deriv + r_p*df*R1_out_r_deriv + 2.0*f_p*f_p*R3_out_r_deriv + f_p*R5_out_r_deriv - r_p*f_p*R1_out_rr_deriv + r_p*f_p*f_p*R3_out_rr_deriv);
		lm_mode->h_in_rs_deriv[num_coupled_fields][i] 	+=	-1.0/(2.0*f_p) * (R1_in*df - R5_in*df - 2.0*f_p*R1_in_r_deriv + r_p*df*R1_in_r_deriv + 2.0*f_p*f_p*R3_in_r_deriv + f_p*R5_in_r_deriv - r_p*f_p*R1_in_rr_deriv + r_p*f_p*f_p*R3_in_rr_deriv);

		lm_mode->h_out_t_deriv[num_coupled_fields][i] 	= 0;
		lm_mode->h_in_t_deriv[num_coupled_fields][i]	= 0;

		//Construct the h7 field and its rs derivative
		lm_mode->h_out[num_coupled_fields+1][i]	+= r_p * R5_out_r_deriv + 2.0*R5_out + l*(l+1.0)*R6_out;
		lm_mode->h_in[num_coupled_fields+1][i]  += r_p * R5_in_r_deriv + 2.0*R5_in + l*(l+1.0)*R6_in;

		double complex R6_out_r_deriv 	= 1.0/f_p*lm_mode->h_out_rs_deriv[num_coupled_fields][i];
		double complex R6_in_r_deriv	= 1.0/f_p*lm_mode->h_in_rs_deriv[num_coupled_fields][i];

		double complex R7_out			= lm_mode->h_out[num_coupled_fields+1][i];
		double complex R7_in			= lm_mode->h_in[num_coupled_fields+1][i];

		R5_in_rr_deriv 			=	1.0/(f_p*f_p) * ( 4.0*V_l * R5_in + 4.0*f_p/(r_p*r_p)*((1.0-4.5/r_p)*R5_in - 0.5*l*(l+1.0)*(R1_in - f_p*R3_in) + 0.5*(1.0-3.0/r_p)*( l*(l+1.0)*R6_in - R7_in)) - f_p*df*R5_in_r_deriv);

		R5_out_rr_deriv 		=	1.0/(f_p*f_p) * ( 4.0*V_l * R5_out + 4.0*f_p/(r_p*r_p)*((1.0-4.5/r_p)*R5_out - 0.5*l*(l+1.0)*(R1_out - f_p*R3_out) + 0.5*(1.0-3.0/r_p)*( l*(l+1.0)*R6_out - R7_out)) - f_p*df*R5_out_r_deriv);


		lm_mode->h_out_rs_deriv[num_coupled_fields+1][i]	+= f_p*(r_p * R5_out_rr_deriv + R5_out_r_deriv + 2.0*R5_out_r_deriv + l*(l+1.0)*R6_out_r_deriv);
		lm_mode->h_in_rs_deriv[num_coupled_fields+1][i] 	+= f_p*(r_p * R5_in_rr_deriv + R5_in_r_deriv	+ 2.0*R5_in_r_deriv + l*(l+1.0)*R6_in_r_deriv);

		// The h2 and h4 fiels are zero
		lm_mode->h_out[num_coupled_fields+2][i]				= 0;
		lm_mode->h_in[num_coupled_fields+2][i] 				= 0;
		lm_mode->h_out_rs_deriv[num_coupled_fields+2][i]	= 0;
		lm_mode->h_in_rs_deriv[num_coupled_fields+2][i] 	= 0;
		lm_mode->h_out_t_deriv[num_coupled_fields+2][i]		= 0;
		lm_mode->h_in_t_deriv[num_coupled_fields+2][i] 		= 0;

		lm_mode->h_out[num_coupled_fields+3][i]				= 0;
		lm_mode->h_in[num_coupled_fields+3][i] 				= 0;
		lm_mode->h_out_rs_deriv[num_coupled_fields+3][i]	= 0;
		lm_mode->h_in_rs_deriv[num_coupled_fields+3][i] 	= 0;
		lm_mode->h_out_t_deriv[num_coupled_fields+3][i]		= 0;
		lm_mode->h_in_t_deriv[num_coupled_fields+3][i] 		= 0;

	}


}



