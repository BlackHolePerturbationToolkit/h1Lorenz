/** \file Runge-Kutta.c
 * Provides the functions to integrate the homogeneous fields from the domain boundaries. 
 * The boundary condition functions are set by the coupled set for the lmn mode. The inner and outer fields and their r* derivatives, between r_min and r_max, are stored in the n_mode_data structure.
 */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "Runge-Kutta.h"
#include "ode_systems.h"

extern const double RUNGE_KUTTA_ACCURACY;				//!< Sets the relative and absolute accuracy of the GSL Runge-Kutta algorithm

void integrate_out_field(struct orbital_params *orbit, struct n_mode_data *n_mode)
{
	struct coupled_set *cset = n_mode->cset;

	int d;
	for(d = 0; d < cset->num_coupled_fields; d++){		


		double* a0s = (double*)malloc(cset->num_coupled_fields * sizeof(double));
		int k;
		for(k = 0; k < cset->num_coupled_fields; k++){
			if(k == d) a0s[k] = 1;
			else a0s[k] = 0;
		}
		struct boundary_conditions bcs;
		bcs.R 					= (double complex*)malloc(cset->num_coupled_fields * sizeof(double complex));
		bcs.R_rs_deriv 			= (double complex*)malloc(cset->num_coupled_fields * sizeof(double complex));
		bcs.l 					= n_mode->l;
		bcs.omega				= n_mode->omega;

		cset->out_bc_func(a0s, &bcs);
		
		// Store the asymptotic ampltidues of the homogeneous fields
		for(k = 0; k < cset->num_coupled_fields; k++){
			n_mode->C_out_hom[k][d] = a0s[k];
		}

		const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;

		gsl_odeiv_step * s    	= gsl_odeiv_step_alloc (T, 4*cset->num_coupled_fields);
		gsl_odeiv_evolve * e  	= gsl_odeiv_evolve_alloc (4*cset->num_coupled_fields);
		gsl_odeiv_control * c  	= gsl_odeiv_control_y_new (0, RUNGE_KUTTA_ACCURACY);

		gsl_odeiv_system sys 	= {cset->ode_system, NULL, 4*cset->num_coupled_fields, n_mode};

		double t = bcs.rs_boundary;
		double h = -1e-10;

		//Boundary conditions
		double* y = (double*)malloc(4 * cset->num_coupled_fields * sizeof(double));


		// The static even mode bcs require scaling up as they are normally too small (ie. below the numerical thresholds in the code)
		double power_scaling = 0;
		if(n_mode->omega == 0){
			for(k = 0; k < cset->num_coupled_fields; k++){
				//printf("%e %e\n", creal(bcs.R[k]), -log10(fabs(creal(bcs.R[k]))));
				if( -log10(fabs(creal(bcs.R[k]))) > power_scaling) power_scaling = -log10(fabs(creal(bcs.R[k])));
			}
		}
		//printf("power_scaling=%e\n", power_scaling);

		for(k = 0; k < cset->num_coupled_fields; k++){
			y[0+4*k] = creal(bcs.R[k])*pow(10.0, power_scaling);
			y[1+4*k] = creal(bcs.R_rs_deriv[k])*pow(10.0, power_scaling);
			y[2+4*k] = cimag(bcs.R[k])*pow(10.0, power_scaling);
			y[3+4*k] = cimag(bcs.R_rs_deriv[k])*pow(10.0, power_scaling);
		}


		int i;
		for (i = orbit->gridsize-1; i >= orbit->ra_grid_index; i--) {
			double r 	= orbit->grid[i];
			double ef 	= f(r);
			double ti 	= rs_of_r(orbit->grid[i]);

			while (t > ti)
			{
				gsl_odeiv_evolve_apply (e, c, s, &sys, &t, ti, &h, y);
				//if(d == 1){ k = 0; printf("%.12lg %.12e %.12e %.12e %.12e %.12e\n", t, y[k+4*0], y[k+4*1], y[k+4*2], y[k+4*3], y[k+4*4]);}
			}
		//	if(i == orbit->gridsize-1 && bcs.l == 2) printf("test2: %e\n", y[0]);

			// store the field data
			for(k = 0; k < cset->num_coupled_fields; k++){
				n_mode->hom_data_out[k][d][0][orbit->gridsize-1 - i] = y[0 + 4*k];
				n_mode->hom_data_out[k][d][1][orbit->gridsize-1 - i] = y[2 + 4*k];
				n_mode->hom_data_out[k][d][2][orbit->gridsize-1 - i] = y[1 + 4*k]/ef;
				n_mode->hom_data_out[k][d][3][orbit->gridsize-1 - i] = y[3 + 4*k]/ef;
			}

			// store the field data at the particle
			if(i == orbit->r0_grid_index ){
				int s;
				for(s = 0; s <= NUM_CHI_VALUES; s++){
					for(k = 0; k < cset->num_coupled_fields; k++){
						n_mode->R_out[k][d][s] 				= y[0 + 4*k] + I*y[2 + 4*k];
						n_mode->R_out_rs_deriv[k][d][s] 	= y[1 + 4*k] + I*y[3 + 4*k];
					}
				}
			}

		}

		gsl_odeiv_evolve_free (e);
		gsl_odeiv_control_free (c);
		gsl_odeiv_step_free (s);
		free(y);
		free(a0s);
		free(bcs.R);
		free(bcs.R_rs_deriv);

	}
}

void integrate_in_field(struct orbital_params *orbit, struct n_mode_data *n_mode)
{
	struct coupled_set *cset = n_mode->cset;

	int d;
	for(d = 0; d < cset->num_coupled_fields; d++){

		double* a0s = (double*)malloc(cset->num_coupled_fields*sizeof(double));
		int k;
		for(k = 0; k < cset->num_coupled_fields; k++){
			if(k == d) a0s[k] = 1;
			else a0s[k] = 0;
		}

		struct boundary_conditions bcs;
		bcs.R 			= (double complex*)malloc(cset->num_coupled_fields * sizeof(double complex));
		bcs.R_rs_deriv 	= (double complex*)malloc(cset->num_coupled_fields * sizeof(double complex));
		bcs.l 			= n_mode->l;
		bcs.omega		= n_mode->omega;

		cset->in_bc_func(a0s, &bcs);
		
		// Store the asymptotic ampltidues of the homogeneous fields
		for(k = 0; k < cset->num_coupled_fields; k++){
			n_mode->C_in_hom[k][d] = a0s[k];
		}
		
		
		const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;

		gsl_odeiv_step * s     = gsl_odeiv_step_alloc (T, 4*cset->num_coupled_fields);
		gsl_odeiv_evolve * e   = gsl_odeiv_evolve_alloc (4*cset->num_coupled_fields);
		gsl_odeiv_control * c  = gsl_odeiv_control_y_new (RUNGE_KUTTA_ACCURACY, RUNGE_KUTTA_ACCURACY);

		gsl_odeiv_system sys = {cset->ode_system, NULL, 4*cset->num_coupled_fields, n_mode};

		double t = bcs.rs_boundary;
		double h = 1e-10;

		//Boundary conditions
		double* y = (double*)malloc(4*cset->num_coupled_fields*sizeof(double));

		for(k = 0; k < cset->num_coupled_fields; k++){
			y[0+4*k] = creal(bcs.R[k]);
			y[1+4*k] = creal(bcs.R_rs_deriv[k]);
			y[2+4*k] = cimag(bcs.R[k]);
			y[3+4*k] = cimag(bcs.R_rs_deriv[k]);
		}

		int i;

		for (i = 0; i <= orbit->rb_grid_index; i++) {
			double r 	= orbit->grid[i];
			double ef 	= f(r);
			double ti = rs_of_r(orbit->grid[i]);
	
			while (t < ti)
			{				
				gsl_odeiv_evolve_apply (e, c, s, &sys, &t, ti, &h, y);
			//	if(d == 0){ k = 0; printf("%.12lg %.12e %.12e %.12e %.12e %.12e\n", t, y[k+4*0], y[k+4*1], y[k+4*2], y[k+4*3], y[k+4*4]);}
			}

			// store the field data
			for(k = 0; k < cset->num_coupled_fields; k++){
				n_mode->hom_data_in[k][d][0][i] = y[0 + 4*k];
				n_mode->hom_data_in[k][d][1][i] = y[2 + 4*k];
				n_mode->hom_data_in[k][d][2][i] = y[1 + 4*k]/ef;
				n_mode->hom_data_in[k][d][3][i] = y[3 + 4*k]/ef;
			}

			// store the field data at the particle
			if(i == orbit->r0_grid_index){
				int s;
				for(s = 0; s <= NUM_CHI_VALUES; s++){
					for(k = 0; k < cset->num_coupled_fields; k++){
						n_mode->R_in[k][d][s] 				= y[0 + 4*k] + I*y[2 + 4*k];
						n_mode->R_in_rs_deriv[k][d][s] 		= y[1 + 4*k] + I*y[3 + 4*k];
					}
				}
			}

		}

		gsl_odeiv_evolve_free (e);
		gsl_odeiv_control_free (c);
		gsl_odeiv_step_free (s);
		free(y);
		free(a0s);
		free(bcs.R);
		free(bcs.R_rs_deriv);
	}
		//exit(0);
}

void integrate_field_equations(struct orbital_params *orbit, struct n_mode_data *n_mode)
{

	integrate_out_field(orbit, n_mode);
	integrate_in_field(orbit, n_mode);

}
