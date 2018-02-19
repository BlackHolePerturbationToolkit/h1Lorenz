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

double R8_static_matching_integrand_on_grid(double chi, void *params)
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

double R8_integrate_static_component_on_grid(struct orbital_params *orbit, struct n_mode_data *n_mode, int piece)
{
	struct static_integrand_params params = {n_mode, piece};

	double result, error = 0.0;

	gsl_function F;
	F.function = &R8_static_matching_integrand_on_grid;
	F.params = &params;

	if(orbit->e == 0) return M_PI*R8_static_matching_integrand_on_grid(M_PI/2.0, &params);

	gsl_integration_workspace *w = gsl_integration_workspace_alloc (100000);

	double return_code = gsl_integration_qag (&F, 0.0, M_PI, MATCHING_INTEGRATION_ACCURACY, MATCHING_INTEGRATION_ACCURACY, 100000, GSL_INTEG_GAUSS61, w, &result, &error); 
	if( return_code == GSL_EROUND){

	}

	gsl_integration_workspace_free (w);
	
	return result;

}

// The data input is not used here
void construct_R8_gauge_field_on_grid(double ***data, struct n_mode_data *n_mode, struct orbital_params *orbit, int num_coupled_fields){

	double l 		= n_mode->l;
	double m		= n_mode->m;
	double omega 	= n_mode->omega;

	// Variables used for the static case
	double c_out = 0;
	double c_in = 0;

	int i;
	for(i = 0; i < orbit->gridsize + 1; i++){

		double r_p;
		double f_p;
		double df;
		if(i == orbit->gridsize){
			r_p = orbit->grid[orbit->r0_grid_index];
			f_p = f(r_p);
			df = 2.0/(r_p*r_p);
		}else{
			r_p = orbit->grid[i];
			f_p = f(r_p);
			df = 2.0/(r_p*r_p);
		}


		if( n_mode->m == 0 && n_mode->n == 0){
			double R8_coefficient 	= -4.0*16.0*M_PI*orbit->L * normalized_associated_legendre_P_theta_deriv(n_mode->l, n_mode->m, M_PI/2.0);
			
			// Only compute these once and place the right-derivative at the particle at the end of the inhom array
			if(i==0){		
				c_out 			= R8_coefficient*R8_integrate_static_component_on_grid(orbit, n_mode, OUT_RE);
				c_in			= R8_coefficient*R8_integrate_static_component_on_grid(orbit, n_mode, IN_RE);
				
				//n_mode->inhom_data[num_coupled_fields][0][orbit->gridsize]	= c_out*R8_static_out_solution(l, orbit->grid[orbit->r0_grid_index]);
				//n_mode->inhom_data[num_coupled_fields][2][orbit->gridsize]	= c_out*R8_static_out_rs_deriv_solution(l,orbit->grid[orbit->r0_grid_index])/f(orbit->grid[orbit->r0_grid_index]);	
				
			}
			if(i <= orbit->r0_grid_index){
				n_mode->inhom_data[num_coupled_fields][0][i]	= c_in*R8_static_in_solution(l,r_p);
				n_mode->inhom_data[num_coupled_fields][2][i]	= c_in*R8_static_in_rs_deriv_solution(l,r_p)/f_p;
			}else{
				n_mode->inhom_data[num_coupled_fields][0][i]	= c_out*R8_static_out_solution(l,r_p);
				n_mode->inhom_data[num_coupled_fields][2][i]	= c_out*R8_static_out_rs_deriv_solution(l,r_p)/f_p;			
			}
				n_mode->inhom_data[num_coupled_fields][1][i]	= 0;			// The imaginary part is zero
				n_mode->inhom_data[num_coupled_fields][3][i]	= 0;

			/*if(i < orbit->rb_grid_index){
				n_mode->hom_data_in[0][0][0][i] = R8_static_in_solution(l,r_p);
				n_mode->hom_data_in[0][0][1][i] = 0;
				n_mode->hom_data_in[0][0][2][i] = R8_static_in_rs_deriv_solution(l,r_p)/f_p;
				n_mode->hom_data_in[0][0][3][i] = 0;
			}
			if(i >= orbit->ra_grid_index){
				//if(i==orbit->gridsize-orbit->r0_grid_index-1) printf("i = %d, j = %d\n", i, orbit->gridsize-i-1);
				n_mode->hom_data_out[0][0][0][orbit->gridsize-i-1] = R8_static_out_solution(l,r_p);
				n_mode->hom_data_out[0][0][1][orbit->gridsize-i-1] = 0;
				n_mode->hom_data_out[0][0][2][orbit->gridsize-i-1] = R8_static_out_rs_deriv_solution(l,r_p)/f_p;
				n_mode->hom_data_out[0][0][3][orbit->gridsize-i-1] = 0;
			}*/

			//if(i == orbit->rb_grid_index) printf("testing: %.14e\n", n_mode->hom_data_in[0][0][0][i]);
			
			// The arrays are calloc'ed so no need to zero
			n_mode->inhom_data[0][0][i] 		= 0;
			n_mode->inhom_data[0][1][i]			= 0;
			n_mode->inhom_data[0][2][i]			= 0;
			n_mode->inhom_data[0][3][i]			= 0;

			if(l>1){
				n_mode->inhom_data[1][0][i] 		= 0;
				n_mode->inhom_data[1][1][i]			= 0;
				n_mode->inhom_data[1][2][i]			= 0;
				n_mode->inhom_data[1][3][i]			= 0;
			}
		}else{

			double complex R9_out 			= 0;
			double complex R9_out_r_deriv 	= 0;
			double complex R9_out_rr_deriv = 0;

			double complex R10_out 		= 0;
			double complex R10_out_r_deriv = 0;

			double V_l 		= 0.25*f_p*(df/r_p + l*(l+1.0)/(r_p*r_p));

			R9_out						= data[0][0][i] + I*data[0][1][i];
			R9_out_r_deriv				= data[0][2][i] + I*data[0][3][i];

			if(num_coupled_fields != 1){
				R10_out					= data[1][0][i] + I*data[1][1][i];
				R10_out_r_deriv			= data[1][2][i] + I*data[1][3][i];
			}	

			R9_out_rr_deriv				= 1.0/(f_p*f_p) *( ( 4.0*V_l + 4.0*f_p/(r_p*r_p) * (1.0-4.5/r_p) - omega*omega) * R9_out - 2.0*f_p/(r_p*r_p)*(1.0-3.0/r_p)*R10_out - f_p*df*R9_out_r_deriv);

			if(data[0][0][i] == n_mode->hR[0][0][i]) R9_out_rr_deriv += n_mode->src[0][i];

			// Construct the h8 field from the gauge condition
			double complex R8 			= I*f_p/omega *  ( R9_out_r_deriv + 2.0 * R9_out/r_p - R10_out/r_p) ;

			double complex R8_r_deriv 	=  ( I/omega * ( f_p * R9_out_rr_deriv + df * R9_out_r_deriv + 2.0*f_p * (R9_out_r_deriv/r_p - R9_out/(r_p*r_p) ) + 2.0 * df * R9_out/r_p - f_p * (R10_out_r_deriv/r_p - R10_out/(r_p*r_p) ) - df * R10_out/r_p) );

			data[num_coupled_fields][0][i] = creal(R8);
			data[num_coupled_fields][1][i] = cimag(R8);
			data[num_coupled_fields][2][i] = creal(R8_r_deriv);
			data[num_coupled_fields][3][i] = cimag(R8_r_deriv);


		}
	}


}

void construct_R2_R4_gauge_fields_on_grid(double ***data, struct n_mode_data *n_mode, struct orbital_params *orbit, int num_coupled_fields)
{	
	double l 		= n_mode->l;
	double omega 	= n_mode->omega;

	int i;
	for(i = 0; i < orbit->gridsize + 1; i++){
		
		double r_p;
		double f_p;
		double df;
		if(i == orbit->gridsize){
			r_p = orbit->grid[orbit->r0_grid_index];
			f_p = f(r_p);
			df = 2.0/(r_p*r_p);
		}else{
			r_p = orbit->grid[i];
			f_p = f(r_p);
			df = 2.0/(r_p*r_p);
		}

		double complex R1_out				= 0;
		double complex R1_out_r_deriv 		= 0;
		double complex R1_out_rr_deriv		= 0;

		double complex R3_out				= 0;
		double complex R3_out_r_deriv 		= 0;
		double complex R3_out_rr_deriv		= 0;

		double complex R5_out				= 0;
		double complex R5_out_r_deriv 		= 0;
		double complex R5_out_rr_deriv		= 0;
		
		double complex R6_out				= 0;
		double complex R6_out_r_deriv		= 0;

		double complex R7_out				= 0;
		double complex R7_out_r_deriv		= 0;

		if(omega != 0){
			R1_out						+= data[0][0][i] + I*data[0][1][i];
			R1_out_r_deriv				+= data[0][2][i] + I*data[0][3][i];

			R3_out						+= data[1][0][i] + I*data[1][1][i];
			R3_out_r_deriv				+= data[1][2][i] + I*data[1][3][i];

			if(l == 0){
				R6_out					+= data[2][0][i] + I*data[2][1][i];
				R6_out_r_deriv			+= data[2][2][i] + I*data[2][3][i];
			}else{
				R5_out					+= data[2][0][i] + I*data[2][1][i];
				R5_out_r_deriv			+= data[2][2][i] + I*data[2][3][i];

				R6_out					+= data[3][0][i] + I*data[3][1][i];
				R6_out_r_deriv			+= data[3][2][i] + I*data[3][3][i];

				if(l >= 2){
					R7_out				+= data[4][0][i] + I*data[4][1][i];
					R7_out_r_deriv		+= data[4][2][i] + I*data[4][3][i];
				}

			}

			double V_l 				= 0.25*f_p*(df/r_p + l*(l+1.0)/(r_p*r_p)) - 0.25*omega*omega;


			R1_out_rr_deriv			= 	1.0/(f_p*f_p) * ( 4.0*V_l * R1_out + 4.0/(r_p*r_p) * f_p * f_p * R3_out_r_deriv + 2.0*f_p/(r_p*r_p) * (1.0 - 4.0/r_p)*(R1_out - R5_out - f_p*R3_out) - 2.0*f_p*f_p/(r_p*r_p)*(1.0-6.0/r_p)*R6_out  - f_p*df*R1_out_r_deriv);

			R3_out_rr_deriv 		= 	1.0/(f_p*f_p) *  ( 4.0*V_l * R3_out - 2.0*f_p/(r_p*r_p) * (R1_out - R5_out - (1.0-4.0/r_p)*(R3_out + R6_out))  - f_p*df*R3_out_r_deriv);

			R5_out_rr_deriv 		=	1.0/(f_p*f_p) * ( 4.0*V_l * R5_out + 4.0*f_p/(r_p*r_p)*((1.0-4.5/r_p)*R5_out - 0.5*l*(l+1.0)*(R1_out - f_p*R3_out) + 0.5*(1.0-3.0/r_p)*( l*(l+1.0)*R6_out - R7_out)) - f_p*df*R5_out_r_deriv);

			// Only needed for the 1st order calculation using the effective-source. Comment out for now.
			/*if(data[0][0][i] == n_mode->hR[0][0][i]){
				 R1_out_rr_deriv += n_mode->src[0][i];
				 R3_out_rr_deriv += n_mode->src[1][i];
				 R5_out_rr_deriv += n_mode->src[2][i];
			}*/
			
			//if(i == orbit->r0_grid_index) printf("R1 test: %.14e %.14e\n", creal(R5_out_rr_deriv), cimag(R5_out_rr_deriv));
			
			 
			// construct the h2 field
			double complex R2_out			= (I*f_p/omega * ( R1_out_r_deriv - f_p * R3_out_r_deriv + 1.0/r_p * (R1_out - R5_out - f_p*R3_out - 2.0*f_p*R6_out) ));
			double complex R2_out_r_deriv 	= (I*df/omega * ( R1_out_r_deriv - f_p * R3_out_r_deriv + 1.0/r_p * (R1_out - R5_out - f_p*R3_out - 2.0*f_p*R6_out) ) + I*f_p/omega *(R1_out_rr_deriv - f_p*R3_out_rr_deriv - df*R3_out_r_deriv + R1_out_r_deriv/r_p - R1_out/(r_p*r_p) - R5_out_r_deriv/r_p + R5_out/(r_p*r_p) - f_p/r_p*R3_out_r_deriv - df/r_p*R3_out + f_p/(r_p*r_p)*R3_out - 2.0*f_p/r_p*R6_out_r_deriv - 2.0*df*R6_out/r_p + 2.0*f_p*R6_out/(r_p*r_p) ));
			double complex R4_out 			= (I*f_p/omega * (R5_out_r_deriv + 2.0*R5_out/r_p + l*(l+1.0) * R6_out/r_p - R7_out/r_p));
			double complex R4_out_r_deriv 	= (I*df/omega * (R5_out_r_deriv + 2.0*R5_out/r_p + l*(l+1.0) * R6_out/r_p - R7_out/r_p)  + I*f_p/omega* (R5_out_rr_deriv + 2.0*R5_out_r_deriv/r_p - 2.0*R5_out/(r_p*r_p) + l*(l+1.0)*(R6_out_r_deriv/r_p - R6_out/(r_p*r_p)) - R7_out_r_deriv/r_p + R7_out/(r_p*r_p) ));

			data[num_coupled_fields][0][i] = creal(R2_out);
			data[num_coupled_fields][1][i] = cimag(R2_out);
			data[num_coupled_fields][2][i] = creal(R2_out_r_deriv);
			data[num_coupled_fields][3][i] = cimag(R2_out_r_deriv);

			if(l >=1){
				//construct the h4 field
				data[num_coupled_fields+1][0][i] = creal(R4_out);
				data[num_coupled_fields+1][1][i] = cimag(R4_out);
				data[num_coupled_fields+1][2][i] = creal(R4_out_r_deriv);
				data[num_coupled_fields+1][3][i] = cimag(R4_out_r_deriv);		
			}

			

			//check gauge conditions
			/*if(i == orbit->gridsize) {
				double complex first_gauge_out		= I*omega*R1_out + f_p*(I*omega*R3_out 	+ R2_out_r_deriv 	+ R2_out/r_p 	- R4_out/r_p);

				nprintf("1st gauge condition on grid: %.3e ", cabs(first_gauge_out));
			}*/


		}


	}

	

}

void construct_R6_R7_R2_R4_gauge_fields_on_grid(double ***data, struct n_mode_data *n_mode, struct orbital_params *orbit, int num_coupled_fields)
{
	double l 		= n_mode->l;
	
	int i;
	for(i = 0; i < orbit->gridsize + 1; i++){

		double r_p;
		double f_p;
		double df;
		if(i == orbit->gridsize){
			r_p = orbit->grid[orbit->r0_grid_index];
			f_p = f(r_p);
			df = 2.0/(r_p*r_p);
		}else{
			r_p = orbit->grid[i];
			f_p = f(r_p);
			df = 2.0/(r_p*r_p);
		}

		double complex R1_out_rr_deriv		= 0;
		double complex R3_out_rr_deriv		= 0;
		double complex R5_out_rr_deriv		= 0;

		double complex R1_out						= data[0][0][i] + I*data[0][1][i];
		double complex R1_out_r_deriv				= data[0][2][i] + I*data[0][3][i];

		double complex	R3_out						= data[1][0][i] + I*data[1][1][i];
		double complex R3_out_r_deriv				= data[1][2][i] + I*data[1][3][i];

		double complex R5_out						= data[2][0][i] + I*data[2][1][i];
		double complex R5_out_r_deriv				= data[2][2][i] + I*data[2][3][i];

		double V_l 		= 0.25*f_p*(df/r_p + l*(l+1.0)/(r_p*r_p));

		// Construct the h6 field and its rs and t derivatives
		// No need to do n folding as this is only used in the n=0 mode. This mode is real.
		double complex R6_out 	= 1.0/(2.0*f_p)*(r_p * R1_out_r_deriv - r_p*f_p*R3_out_r_deriv + R1_out - R5_out - f_p*R3_out);

		R1_out_rr_deriv		= 	1.0/(f_p*f_p) * ( 4.0*V_l * R1_out + 4.0/(r_p*r_p) * f_p * f_p * R3_out_r_deriv + 2.0*f_p/(r_p*r_p) * (1.0 - 4.0/r_p)*(R1_out - R5_out - f_p*R3_out) - 2.0*f_p*f_p/(r_p*r_p)*(1.0-6.0/r_p)*R6_out  - f_p*df*R1_out_r_deriv);

		R3_out_rr_deriv 	= 	1.0/(f_p*f_p) *  ( 4.0*V_l * R3_out - 2.0*f_p/(r_p*r_p) * (R1_out - R5_out - (1.0-4.0/r_p)*(R3_out + R6_out))  - f_p*df*R3_out_r_deriv);

		//Construct the h7 field and its rs derivative
		double complex R6_out_r_deriv 	= 1.0/f_p*(-1.0/(2.0*f_p) * (R1_out*df - R5_out*df - 2.0*f_p*R1_out_r_deriv + r_p*df*R1_out_r_deriv + 2.0*f_p*f_p*R3_out_r_deriv + f_p*R5_out_r_deriv - r_p*f_p*R1_out_rr_deriv + r_p*f_p*f_p*R3_out_rr_deriv));

		double complex R7_out			= r_p * R5_out_r_deriv + 2.0*R5_out + l*(l+1.0)*R6_out;

		R5_out_rr_deriv 		= 1.0/(f_p*f_p) * ( 4.0*V_l * R5_out + 4.0*f_p/(r_p*r_p)*((1.0-4.5/r_p)*R5_out - 0.5*l*(l+1.0)*(R1_out - f_p*R3_out) + 0.5*(1.0-3.0/r_p)*( l*(l+1.0)*R6_out - R7_out)) - f_p*df*R5_out_r_deriv);
		

		double complex R7_out_r_deriv	= (r_p * R5_out_rr_deriv + R5_out_r_deriv + 2.0*R5_out_r_deriv + l*(l+1.0)*R6_out_r_deriv);

		data[num_coupled_fields][0][i] = creal(R6_out);
		data[num_coupled_fields][1][i] = cimag(R6_out);
		data[num_coupled_fields][2][i] = creal(R6_out_r_deriv);
		data[num_coupled_fields][3][i] = cimag(R6_out_r_deriv);
		
		data[num_coupled_fields+1][0][i] = creal(R7_out);
		data[num_coupled_fields+1][1][i] = cimag(R7_out);
		data[num_coupled_fields+1][2][i] = creal(R7_out_r_deriv);
		data[num_coupled_fields+1][3][i] = cimag(R7_out_r_deriv);

		// The h2 and h4 fiels are zero
		data[num_coupled_fields+2][0][i] = 0;
		data[num_coupled_fields+2][1][i] = 0;
		data[num_coupled_fields+2][2][i] = 0;
		data[num_coupled_fields+2][3][i] = 0;

		data[num_coupled_fields+3][0][i] = 0;
		data[num_coupled_fields+3][1][i] = 0;
		data[num_coupled_fields+3][2][i] = 0;
		data[num_coupled_fields+3][3][i] = 0;

	}


}



