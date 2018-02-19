/** \file self_force.c
 *  Contains all the functions used to calculate the GSF from the lm modes
 */
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sf_hyperg.h>
#include "self_force.h"
#include "math_funcs.h"

double Anna_and_Barrys_D_r_circular_orbits(double l, double r1);
double Anna_and_Barrys_D_r(double l, double p, double e, double r1, double L);
double Anna_and_Barrys_D_t(double l, double p, double e, double r1, double ur1, double L);
double Anna_and_Barrys_D_H(double l, double r0, double E0, double L0);
double Anna_and_Barrys_F_H(double l, double r0, double E0, double L0);
double Anna_and_Barrys_F_r(double l, double p, double e, double r1, double ur1, double L);
double Anna_and_Barrys_F_t(double l, double p, double e, double r1, double ur1, double L);


/** Functions that map the coupled set orderings back to the standard ordering*/
void construct_monopole_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{
	fields->h1				= lm_mode->h_in[0][i];
	fields->h2				= lm_mode->h_in[3][i];
	fields->h3				= lm_mode->h_in[1][i];
	fields->h4				= 0;
	fields->h5				= 0;
	fields->h6				= lm_mode->h_in[2][i];
	fields->h7				= 0;
	fields->h8				= 0;
	fields->h9				= 0;
	fields->h10				= 0;

	fields->h1_t_deriv		= lm_mode->h_in_t_deriv[0][i];
	fields->h1_rs_deriv		= lm_mode->h_in_rs_deriv[0][i];
	fields->h2_t_deriv 		= lm_mode->h_in_t_deriv[3][i];
	fields->h2_rs_deriv		= lm_mode->h_in_rs_deriv[3][i];
	fields->h3_t_deriv		= lm_mode->h_in_t_deriv[1][i];
	fields->h3_rs_deriv		= lm_mode->h_in_rs_deriv[1][i];
	fields->h4_t_deriv 		= 0;
	fields->h4_rs_deriv		= 0;
	fields->h5_t_deriv 		= 0;
	fields->h5_rs_deriv		= 0;
	fields->h6_t_deriv		= lm_mode->h_in_t_deriv[2][i];
	fields->h6_rs_deriv		= lm_mode->h_in_rs_deriv[2][i];
	fields->h7_t_deriv 		= 0;
	fields->h7_rs_deriv		= 0;
	fields->h8_t_deriv 		= 0;
	fields->h8_rs_deriv		= 0;
	fields->h9_t_deriv 		= 0;
	fields->h9_rs_deriv		= 0;
	fields->h10_t_deriv 	= 0;
	fields->h10_rs_deriv	= 0;
}

void construct_odd_dipole_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{

	fields->h1				= 0;
	fields->h2				= 0;
	fields->h3				= 0;
	fields->h4				= 0;
	fields->h5				= 0;
	fields->h6				= 0;
	fields->h7				= 0;
	fields->h8				= lm_mode->h_in[1][i];
	fields->h9				= lm_mode->h_in[0][i];
	fields->h10				= 0;
	fields->h1_t_deriv		= 0;
	fields->h1_rs_deriv		= 0;
	fields->h2_t_deriv 		= 0;
	fields->h2_rs_deriv		= 0;
	fields->h3_t_deriv		= 0;
	fields->h3_rs_deriv		= 0;
	fields->h4_t_deriv 		= 0;
	fields->h4_rs_deriv		= 0;
	fields->h5_t_deriv 		= 0;
	fields->h5_rs_deriv		= 0;
	fields->h6_t_deriv		= 0;
	fields->h6_rs_deriv		= 0;
	fields->h7_t_deriv 		= 0;
	fields->h7_rs_deriv		= 0;
	fields->h8_t_deriv 		= lm_mode->h_in_t_deriv[1][i];
	fields->h8_rs_deriv		= lm_mode->h_in_rs_deriv[1][i];
	fields->h9_t_deriv 		= lm_mode->h_in_t_deriv[0][i];
	fields->h9_rs_deriv		= lm_mode->h_in_rs_deriv[0][i];
	fields->h10_t_deriv 	= 0;
	fields->h10_rs_deriv	= 0;
}

void construct_even_dipole_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{
	fields->h1				= lm_mode->h_in[0][i];
	fields->h2				= lm_mode->h_in[4][i];
	fields->h3				= lm_mode->h_in[1][i];
	fields->h4				= lm_mode->h_in[5][i];
	fields->h5				= lm_mode->h_in[2][i];
	fields->h6				= lm_mode->h_in[3][i];
	fields->h7				= 0;
	fields->h8				= 0;
	fields->h9				= 0;
	fields->h10				= 0;

	fields->h1_t_deriv		= lm_mode->h_in_t_deriv[0][i];
	fields->h1_rs_deriv		= lm_mode->h_in_rs_deriv[0][i];
	fields->h2_t_deriv 		= lm_mode->h_in_t_deriv[4][i];
	fields->h2_rs_deriv		= lm_mode->h_in_rs_deriv[4][i];
	fields->h3_t_deriv		= lm_mode->h_in_t_deriv[1][i];
	fields->h3_rs_deriv		= lm_mode->h_in_rs_deriv[1][i];
	fields->h4_t_deriv 		= lm_mode->h_in_t_deriv[5][i];
	fields->h4_rs_deriv		= lm_mode->h_in_rs_deriv[5][i];
	fields->h5_t_deriv 		= lm_mode->h_in_t_deriv[2][i];
	fields->h5_rs_deriv		= lm_mode->h_in_rs_deriv[2][i];
	fields->h6_t_deriv		= lm_mode->h_in_t_deriv[3][i];
	fields->h6_rs_deriv		= lm_mode->h_in_rs_deriv[3][i];
	fields->h7_t_deriv 		= 0;
	fields->h7_rs_deriv		= 0;
	fields->h8_t_deriv 		= 0;
	fields->h8_rs_deriv		= 0;
	fields->h9_t_deriv 		= 0;
	fields->h9_rs_deriv		= 0;
	fields->h10_t_deriv 	= 0;
	fields->h10_rs_deriv	= 0;
}

void construct_odd_modes_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{
	fields->h1				= 0;
	fields->h2				= 0;
	fields->h3				= 0;
	fields->h4				= 0;
	fields->h5				= 0;
	fields->h6				= 0;
	fields->h7				= 0;
	fields->h8				= lm_mode->h_in[2][i];
	fields->h9				= lm_mode->h_in[0][i];
	fields->h10				= lm_mode->h_in[1][i];

	fields->h1_t_deriv		= 0;
	fields->h1_rs_deriv		= 0;
	fields->h2_t_deriv 		= 0;
	fields->h2_rs_deriv		= 0;
	fields->h3_t_deriv		= 0;
	fields->h3_rs_deriv		= 0;
	fields->h4_t_deriv 		= 0;
	fields->h4_rs_deriv		= 0;
	fields->h5_t_deriv 		= 0;
	fields->h5_rs_deriv		= 0;
	fields->h6_t_deriv		= 0;
	fields->h6_rs_deriv		= 0;
	fields->h7_t_deriv 		= 0;
	fields->h7_rs_deriv		= 0;
	fields->h8_t_deriv 		= lm_mode->h_in_t_deriv[2][i];
	fields->h8_rs_deriv		= lm_mode->h_in_rs_deriv[2][i];
	fields->h9_t_deriv 		= lm_mode->h_in_t_deriv[0][i];
	fields->h9_rs_deriv		= lm_mode->h_in_rs_deriv[0][i];
	fields->h10_t_deriv 	= lm_mode->h_in_t_deriv[1][i];
	fields->h10_rs_deriv	= lm_mode->h_in_rs_deriv[1][i];
}

void construct_even_modes_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{
	fields->h1				= lm_mode->h_in[0][i];
	fields->h2				= lm_mode->h_in[5][i];
	fields->h3				= lm_mode->h_in[1][i];
	fields->h4				= lm_mode->h_in[6][i];
	fields->h5				= lm_mode->h_in[2][i];
	fields->h6				= lm_mode->h_in[3][i];
	fields->h7				= lm_mode->h_in[4][i];
	fields->h8				= 0;
	fields->h9				= 0;
	fields->h10				= 0;

	fields->h1_t_deriv		= lm_mode->h_in_t_deriv[0][i];
	fields->h1_rs_deriv		= lm_mode->h_in_rs_deriv[0][i];
	fields->h2_t_deriv 		= lm_mode->h_in_t_deriv[5][i];
	fields->h2_rs_deriv		= lm_mode->h_in_rs_deriv[5][i];
	fields->h3_t_deriv		= lm_mode->h_in_t_deriv[1][i];
	fields->h3_rs_deriv		= lm_mode->h_in_rs_deriv[1][i];
	fields->h4_t_deriv 		= lm_mode->h_in_t_deriv[6][i];
	fields->h4_rs_deriv		= lm_mode->h_in_rs_deriv[6][i];
	fields->h5_t_deriv 		= lm_mode->h_in_t_deriv[2][i];
	fields->h5_rs_deriv		= lm_mode->h_in_rs_deriv[2][i];
	fields->h6_t_deriv		= lm_mode->h_in_t_deriv[3][i];
	fields->h6_rs_deriv		= lm_mode->h_in_rs_deriv[3][i];
	fields->h7_t_deriv 		= lm_mode->h_in_t_deriv[4][i];
	fields->h7_rs_deriv		= lm_mode->h_in_rs_deriv[4][i];
	fields->h8_t_deriv 		= 0;
	fields->h8_rs_deriv		= 0;
	fields->h9_t_deriv 		= 0;
	fields->h9_rs_deriv		= 0;
	fields->h10_t_deriv 	= 0;
	fields->h10_rs_deriv	= 0;
}

void construct_monopole_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{
	fields->h1				= lm_mode->h_out[0][i];
	fields->h2				= lm_mode->h_out[3][i];
	fields->h3				= lm_mode->h_out[1][i];
	fields->h4				= 0;
	fields->h5				= 0;
	fields->h6				= lm_mode->h_out[2][i];
	fields->h7				= 0;
	fields->h8				= 0;
	fields->h9				= 0;
	fields->h10				= 0;

	fields->h1_t_deriv		= lm_mode->h_out_t_deriv[0][i];
	fields->h1_rs_deriv		= lm_mode->h_out_rs_deriv[0][i];
	fields->h2_t_deriv 		= lm_mode->h_out_t_deriv[3][i];
	fields->h2_rs_deriv		= lm_mode->h_out_rs_deriv[3][i];
	fields->h3_t_deriv		= lm_mode->h_out_t_deriv[1][i];
	fields->h3_rs_deriv		= lm_mode->h_out_rs_deriv[1][i];
	fields->h4_t_deriv 		= 0;
	fields->h4_rs_deriv		= 0;
	fields->h5_t_deriv 		= 0;
	fields->h5_rs_deriv		= 0;
	fields->h6_t_deriv		= lm_mode->h_out_t_deriv[2][i];
	fields->h6_rs_deriv		= lm_mode->h_out_rs_deriv[2][i];
	fields->h7_t_deriv 		= 0;
	fields->h7_rs_deriv		= 0;
	fields->h8_t_deriv 		= 0;
	fields->h8_rs_deriv		= 0;
	fields->h9_t_deriv 		= 0;
	fields->h9_rs_deriv		= 0;
	fields->h10_t_deriv 	= 0;
	fields->h10_rs_deriv	= 0;
}

void construct_odd_dipole_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{

	fields->h1				= 0;
	fields->h2				= 0;
	fields->h3				= 0;
	fields->h4				= 0;
	fields->h5				= 0;
	fields->h6				= 0;
	fields->h7				= 0;
	fields->h8				= lm_mode->h_out[1][i];
	fields->h9				= lm_mode->h_out[0][i];
	fields->h10				= 0;
	fields->h1_t_deriv		= 0;
	fields->h1_rs_deriv		= 0;
	fields->h2_t_deriv 		= 0;
	fields->h2_rs_deriv		= 0;
	fields->h3_t_deriv		= 0;
	fields->h3_rs_deriv		= 0;
	fields->h4_t_deriv 		= 0;
	fields->h4_rs_deriv		= 0;
	fields->h5_t_deriv 		= 0;
	fields->h5_rs_deriv		= 0;
	fields->h6_t_deriv		= 0;
	fields->h6_rs_deriv		= 0;
	fields->h7_t_deriv 		= 0;
	fields->h7_rs_deriv		= 0;
	fields->h8_t_deriv 		= lm_mode->h_out_t_deriv[1][i];
	fields->h8_rs_deriv		= lm_mode->h_out_rs_deriv[1][i];
	fields->h9_t_deriv 		= lm_mode->h_out_t_deriv[0][i];
	fields->h9_rs_deriv		= lm_mode->h_out_rs_deriv[0][i];
	fields->h10_t_deriv 	= 0;
	fields->h10_rs_deriv	= 0;
}

void construct_even_dipole_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{
	fields->h1				= lm_mode->h_out[0][i];
	fields->h2				= lm_mode->h_out[4][i];
	fields->h3				= lm_mode->h_out[1][i];
	fields->h4				= lm_mode->h_out[5][i];
	fields->h5				= lm_mode->h_out[2][i];
	fields->h6				= lm_mode->h_out[3][i];
	fields->h7				= 0;
	fields->h8				= 0;
	fields->h9				= 0;
	fields->h10				= 0;

	fields->h1_t_deriv		= lm_mode->h_out_t_deriv[0][i];
	fields->h1_rs_deriv		= lm_mode->h_out_rs_deriv[0][i];
	fields->h2_t_deriv 		= lm_mode->h_out_t_deriv[4][i];
	fields->h2_rs_deriv		= lm_mode->h_out_rs_deriv[4][i];
	fields->h3_t_deriv		= lm_mode->h_out_t_deriv[1][i];
	fields->h3_rs_deriv		= lm_mode->h_out_rs_deriv[1][i];
	fields->h4_t_deriv 		= lm_mode->h_out_t_deriv[5][i];
	fields->h4_rs_deriv		= lm_mode->h_out_rs_deriv[5][i];
	fields->h5_t_deriv 		= lm_mode->h_out_t_deriv[2][i];
	fields->h5_rs_deriv		= lm_mode->h_out_rs_deriv[2][i];
	fields->h6_t_deriv		= lm_mode->h_out_t_deriv[3][i];
	fields->h6_rs_deriv		= lm_mode->h_out_rs_deriv[3][i];
	fields->h7_t_deriv 		= 0;
	fields->h7_rs_deriv		= 0;
	fields->h8_t_deriv 		= 0;
	fields->h8_rs_deriv		= 0;
	fields->h9_t_deriv 		= 0;
	fields->h9_rs_deriv		= 0;
	fields->h10_t_deriv 	= 0;
	fields->h10_rs_deriv	= 0;
}

void construct_odd_modes_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{
	fields->h1				= 0;
	fields->h2				= 0;
	fields->h3				= 0;
	fields->h4				= 0;
	fields->h5				= 0;
	fields->h6				= 0;
	fields->h7				= 0;
	fields->h8				= lm_mode->h_out[2][i];
	fields->h9				= lm_mode->h_out[0][i];
	fields->h10				= lm_mode->h_out[1][i];

	fields->h1_t_deriv		= 0;
	fields->h1_rs_deriv		= 0;
	fields->h2_t_deriv 		= 0;
	fields->h2_rs_deriv		= 0;
	fields->h3_t_deriv		= 0;
	fields->h3_rs_deriv		= 0;
	fields->h4_t_deriv 		= 0;
	fields->h4_rs_deriv		= 0;
	fields->h5_t_deriv 		= 0;
	fields->h5_rs_deriv		= 0;
	fields->h6_t_deriv		= 0;
	fields->h6_rs_deriv		= 0;
	fields->h7_t_deriv 		= 0;
	fields->h7_rs_deriv		= 0;
	fields->h8_t_deriv 		= lm_mode->h_out_t_deriv[2][i];
	fields->h8_rs_deriv		= lm_mode->h_out_rs_deriv[2][i];
	fields->h9_t_deriv 		= lm_mode->h_out_t_deriv[0][i];
	fields->h9_rs_deriv		= lm_mode->h_out_rs_deriv[0][i];
	fields->h10_t_deriv 	= lm_mode->h_out_t_deriv[1][i];
	fields->h10_rs_deriv	= lm_mode->h_out_rs_deriv[1][i];
}

void construct_even_modes_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i)
{
	fields->h1				= lm_mode->h_out[0][i];
	fields->h2				= lm_mode->h_out[5][i];
	fields->h3				= lm_mode->h_out[1][i];
	fields->h4				= lm_mode->h_out[6][i];
	fields->h5				= lm_mode->h_out[2][i];
	fields->h6				= lm_mode->h_out[3][i];
	fields->h7				= lm_mode->h_out[4][i];
	fields->h8				= 0;
	fields->h9				= 0;
	fields->h10				= 0;

	fields->h1_t_deriv		= lm_mode->h_out_t_deriv[0][i];
	fields->h1_rs_deriv		= lm_mode->h_out_rs_deriv[0][i];
	fields->h2_t_deriv 		= lm_mode->h_out_t_deriv[5][i];
	fields->h2_rs_deriv		= lm_mode->h_out_rs_deriv[5][i];
	fields->h3_t_deriv		= lm_mode->h_out_t_deriv[1][i];
	fields->h3_rs_deriv		= lm_mode->h_out_rs_deriv[1][i];
	fields->h4_t_deriv 		= lm_mode->h_out_t_deriv[6][i];
	fields->h4_rs_deriv		= lm_mode->h_out_rs_deriv[6][i];
	fields->h5_t_deriv 		= lm_mode->h_out_t_deriv[2][i];
	fields->h5_rs_deriv		= lm_mode->h_out_rs_deriv[2][i];
	fields->h6_t_deriv		= lm_mode->h_out_t_deriv[3][i];
	fields->h6_rs_deriv		= lm_mode->h_out_rs_deriv[3][i];
	fields->h7_t_deriv 		= lm_mode->h_out_t_deriv[4][i];
	fields->h7_rs_deriv		= lm_mode->h_out_rs_deriv[4][i];
	fields->h8_t_deriv 		= 0;
	fields->h8_rs_deriv		= 0;
	fields->h9_t_deriv 		= 0;
	fields->h9_rs_deriv		= 0;
	fields->h10_t_deriv 	= 0;
	fields->h10_rs_deriv	= 0;
}

double C(double l, double m)
{
	return sqrt( (l*l - m*m)/( (2.0*l + 1.0)*(2.0*l - 1.0) ) );
}

double alpha2 (double l, double m)
{
	return -C(l+1.0, m) * C(l+2.0, m);
}

double alpha0(double l, double m)
{		
	return 1.0 - pow(C(l, m), 2.0) - pow(C(l+1.0, m), 2.0);
}

double alpham2(double l, double m)
{	
	return -C(l, m) * C(l-1.0, m);
}

double beta2(double l, double m)
{
	return l * C(l+1.0, m) * C(l+2.0, m);
}

double beta0(double l, double m)
{		
	return l * pow(C(l+1.0, m), 2.0) - (l+1.0) * pow(C(l, m), 2.0);
}

double betam2(double l, double m)
{
	return -(l+1.0) * C(l, m) * C(l-1.0, m);
}

double gamma2(double l, double m)
{	
	return l*l * C(l+1.0, m) * C(l+2.0, m);
}

double gamma0(double l, double m)
{
	return m*m - l*(l+1.0) + l*l * pow(C(l+1.0, m), 2.0) + pow(l+1.0, 2.0) * pow(C(l,m),2.0);
}

double gammam2(double l, double m)
{
	return pow(l+1.0, 2.0) * C(l, m) * C(l-1.0, m);
}

double delta1(double l, double m)
{
	return l * C(l+1.0, m);
}

double deltam1(double l, double m)
{
	return -(l+1.0) * C(l, m);
}

double epsilon1	(double l, double m)
{
	return (1.0 - l) * C(l+1.0, m);
}

double epsilonm1(double l, double m)
{
	return (l+2.0) * C(l, m);
}

double zeta3(double l, double m)
{	
	return -l * C(l+1.0, m) * C(l+2.0, m) * C(l+3.0, m);
}

double zeta1(double l, double m)
{
	return C(l+1.0, m) * (l * (1.0 - pow(C(l+1.0, m), 2.0) - pow(C(l+2.0, m), 2.0)) + (l+1.0) * pow(C(l, m), 2.0));
}

double zetam1(double l, double m)
{
	return -C(l,m) * ( (l+1.0)*(1.0 - pow(C(l-1.0, m), 2.0) - pow(C(l, m), 2.0)) + l * pow(C(l+1.0, m), 2.0) );
}

double zetam3(double l, double m)
{
	return (l+1.0) * C(l, m) * C(l-1.0, m) * C(l-2.0, m);
}

double xi3(double l, double m)
{
	return l*l * C(l+1.0, m) * C(l+2.0, m) * C(l+3.0, m);
}

double xi1(double l, double m)
{
	return C(l+1.0,m) * (m*m - l*(l+1.0) + l*l * pow(C(l+1.0, m),2.0) + pow(l+1.0, 2.0) * pow(C(l, m), 2.0) + l*l * pow(C(l+2.0, m), 2.0) );
}

double xim1	(double l, double m)
{	
	return C(l,m) * (m*m - l*(l+1.0) + l*l * pow(C(l+1.0, m), 2.0) + pow(l+1.0, 2.0)*pow(C(l, m), 2.0) + pow(l+1.0, 2.0)*pow(C(l-1.0, m), 2.0));
}

double xim3	(double l, double m)
{
	return pow(l+1.0, 2.0) * C(l,m) * C(l-1.0,m) * C(l-2.0,m);
}

double complex ft0_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h1 			= fields->h1;
	double complex h2			= fields->h2;
	double complex h3 			= fields->h3;
	double complex h4 			= fields->h4;
	double complex h5 			= fields->h5;
	double complex h6			= fields->h6;
	double complex h7 			= fields->h7;

	double complex h1_t_deriv	= fields->h1_t_deriv;
	double complex h1_rs_deriv	= fields->h1_rs_deriv;
	double complex h2_t_deriv	= fields->h2_t_deriv;
	double complex h2_rs_deriv	= fields->h2_rs_deriv;
	double complex h3_t_deriv	= fields->h3_t_deriv;
	double complex h3_rs_deriv	= fields->h3_rs_deriv;
	double complex h4_t_deriv	= fields->h4_t_deriv;
	double complex h4_rs_deriv	= fields->h4_rs_deriv;
	double complex h5_t_deriv	= fields->h5_t_deriv;
	double complex h5_rs_deriv	= fields->h5_rs_deriv;
	double complex h6_t_deriv	= fields->h6_t_deriv;
	double complex h6_rs_deriv	= fields->h6_rs_deriv;
	double complex h7_t_deriv	= fields->h7_t_deriv;
	double complex h7_rs_deriv	= fields->h7_rs_deriv;



	double complex line1	= E/(4.0*pow(eff,4.0)) * ( pow(contra_ur,3.0) + (E*E *(r+4.0)/r - 2.0*eff)*contra_ur - I*m*eff*Lhat*(pow(contra_ur,2.0) + E*E - 2.0*eff))*h1;

	double complex line2	= -r/(4.0*pow(eff,4.0)) * ((E*E+eff)*pow(contra_ur,2.0) + E*E*(E*E - eff))*h1_t_deriv - r/(4.0*pow(eff,4.0))*contra_ur*E*(pow(contra_ur,2.0) + E*E - 2.0*eff)*h1_rs_deriv;

	double complex line3	= 1.0/(2.0*pow(eff,4.0))*(pow(contra_ur,2.0)*(E*E*(r+1.0)/r - eff*(r-1.0)/r) +1.0/r*E*E*(E*E - eff) - I*m*eff*contra_ur*Lhat*(E*E-eff))*h2;

	double complex line4	= -r*contra_ur/(2.0*pow(eff,4.0))*(E*E*E * h2_t_deriv + contra_ur*(E*E-eff)*h2_rs_deriv) + E*Lhat*Lhat/(4.0*eff)*(contra_ur-I*m*Lhat)*h3 - r*Lhat*Lhat/(4.0*eff*eff)*((E*E + eff)*h3_t_deriv + contra_ur*E*h3_rs_deriv);

	double complex line5	= m/(2.0*l*(l+1.0)*pow(eff,3.0))*Lhat*(2.0*I*contra_ur*(E*E*(r-1.0)/r - eff*eff) + m*eff*Lhat*(E*E - eff))*h4;
	if(l == 0){
		line5 = 0;
	}

	double complex line6	= -I*m*r*Lhat/(2.0*l*(l+1.0)*pow(eff,3.0)) * (E*E*E*h4_t_deriv + contra_ur*(E*E - eff)*h4_rs_deriv) + m*E*Lhat/(2.0*l*(l+1.0)*pow(eff,3.0))*(I*pow(contra_ur,2.0)*(2.0-3.0/r) + I*E*E/r + m*eff*contra_ur*Lhat)*h5;
	if(l == 0){
		line6 = 0;
	}

	double complex line7	=  - I*m*r*contra_ur*Lhat/(2.0*l*(l+1.0)*pow(eff,3.0))*((E*E + eff)*h5_t_deriv + contra_ur*E*h5_rs_deriv) + E/(4.0*eff)*(-contra_ur + I*m*Lhat)*h6 + r/(4.0*eff*eff)*((E*E - eff)*h6_t_deriv + contra_ur*E*h6_rs_deriv);
	if(l == 0){
		line7 =  E/(4.0*eff)*(-contra_ur + I*m*Lhat)*h6 + r/(4.0*eff*eff)*((E*E - eff)*h6_t_deriv + contra_ur*E*h6_rs_deriv);
	}

	double complex line8	= - m*E*Lhat*Lhat/(4.0*l*(l+1.0)*lambda*eff) * (3.0*m*contra_ur - I*Lhat*(4.0 + m*m))*h7 + m*m*r*Lhat*Lhat/(4.0*l*(l+1.0)*lambda*eff*eff)*((E*E + eff)*h7_t_deriv + contra_ur*E*h7_rs_deriv);
	if(l == 0 || l == 1){
		line8 = 0;
	}

   

	return line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8;

}

double complex ft1_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h1 			= fields->h1;
	double complex h2			= fields->h2;
	double complex h3 			= fields->h3;
	double complex h5 			= fields->h5;
	double complex h6			= fields->h6;
	double complex h7 			= fields->h7;

	double complex h6_t_deriv	= fields->h6_t_deriv;
	double complex h6_rs_deriv	= fields->h6_rs_deriv;


	double complex line1 = -contra_ur*E*pow(Lhat,2.0)/(2.0*eff*eff) * h1 - pow(Lhat,2.0)/(2.0*eff*eff)*(E*E - eff)*h2 + contra_ur*E*pow(Lhat,2.0)/(2.0*eff)*h3 - I*m*E*pow(Lhat,3.0)/( 2.0*l*(l+1.0)*eff ) * h5;
	if(l == 0){
		line1 = -contra_ur*E*pow(Lhat,2.0)/(2.0*eff*eff) * h1 - pow(Lhat,2.0)/(2.0*eff*eff)*(E*E - eff)*h2 + contra_ur*E*pow(Lhat,2.0)/(2.0*eff)*h3;
	}

	double complex line2 = E*pow(Lhat,2.0)/(4.0*eff)*(3.0*contra_ur - I*m*Lhat)*h6 - r/(4.0*eff*eff)*Lhat*Lhat*(E*E + eff) * h6_t_deriv - r/(4.0*eff*eff)*contra_ur*E*Lhat*Lhat*h6_rs_deriv - I*m*E*pow(Lhat,3.0)/( l*(l+1.0)*lambda*eff)*h7;
	if(l == 0 || l == 1){
		line2 = E*pow(Lhat,2.0)/(4.0*eff)*(3.0*contra_ur - I*m*Lhat)*h6 - r/(4.0*eff*eff)*Lhat*Lhat*(E*E + eff) * h6_t_deriv - r/(4.0*eff*eff)*contra_ur*E*Lhat*Lhat*h6_rs_deriv;
	}

	return line1 + line2;
}

double complex ft2_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h4 			= fields->h4;
	double complex h5 			= fields->h5;
	double complex h7 			= fields->h7;

	double complex h7_t_deriv	= fields->h7_t_deriv;
	double complex h7_rs_deriv	= fields->h7_rs_deriv;


	double complex line1 = - Lhat*Lhat * (E*E - eff)/(2.0*l*(l+1.0)*eff*eff) * h4 - contra_ur*E*Lhat*Lhat/(2.0*l*(l+1.0)*eff*eff)*h5;
	if(l == 0){
		line1 = 0;
	}

	double complex line2 =   Lhat*Lhat/(4.0*l*(l+1.0)*lambda*eff*eff) * (eff*E*(3.0*contra_ur - 5.0*I*m*Lhat)*h7 - r*(E*E + eff)*h7_t_deriv - r*contra_ur*E*h7_rs_deriv);
	if(l == 0 || l == 1){
		line2 = 0;
	}

	return line1 + line2;
}

double complex ft3_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h7 			= fields->h7;

	double complex h7_t_deriv	= fields->h7_t_deriv;
	double complex h7_rs_deriv	= fields->h7_rs_deriv;

	double complex line1 = Lhat*Lhat/(4.0*l*(l+1.0)*lambda*eff*eff) * (-eff*E*(3.0*contra_ur - I*m*Lhat)*h7 + r*(E*E+eff)*h7_t_deriv + r*contra_ur*E*h7_rs_deriv );;
	if(l == 0 || l == 1){
		line1 = 0;		
	}

	return line1;
}

double complex ft4_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h8 			= fields->h8;
	double complex h9			= fields->h9;
	double complex h10			= fields->h10;

	double complex h10_t_deriv		= fields->h10_t_deriv;
	double complex h10_rs_deriv	= fields->h10_rs_deriv;

	double complex line1 	= -I*m*Lhat*Lhat/(2.0*l*(l+1.0)*eff*eff) * ( (E*E - eff)*h8 + contra_ur*E*h9 );
	if(l == 0){
		line1 		= 0;
	}

	double complex line2 = I*m*Lhat*Lhat/(2.0*l*(l+1.0)*lambda*eff*eff) * (eff*E * (3.0*contra_ur - 2.0*I*m*Lhat)*h10 - r*(E*E + eff)*h10_t_deriv - r*contra_ur*E*h10_rs_deriv); 
	if(l == 0 || l == 1){
		line2 		= 0;
	}

	return line1 + line2;
}

double complex ft5_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h8 			= fields->h8;
	double complex h9 			= fields->h9;
	double complex h10 			= fields->h10;

	double complex h8_t_deriv	= fields->h8_t_deriv;
	double complex h8_rs_deriv	= fields->h8_rs_deriv;
	double complex h9_t_deriv	= fields->h9_t_deriv;
	double complex h9_rs_deriv	= fields->h9_rs_deriv;

	double complex line1		= -contra_ur*Lhat/(l*(l+1.0)*pow(eff,3.0)) * (E*E*(r-1.0)/r - eff*eff)*h8 + r*E*E*E*Lhat/(2.0*l*(l+1.0)*pow(eff,3.0))*h8_t_deriv + r*contra_ur*Lhat/(2.0*l*(l+1.0)*pow(eff,3.0))*(E*E - eff)*h8_rs_deriv;
	if(l == 0){
		line1 			= 0;
	}

	double complex line2		= -E*Lhat/(2.0*l*(l+1.0)*pow(eff,3.0)) * ( pow(contra_ur,2.0) * (2.0*r - 3.0)/r + E*E/r )*h9 + r*contra_ur*Lhat/(2.0*l*(l+1.0)*pow(eff,3.0))*(E*E + eff)*h9_t_deriv + r*pow(contra_ur,2.0)*E*Lhat/(2.0*l*(l+1.0)*pow(eff,3.0))*h9_rs_deriv;
	if(l == 0){
		line2			= 0;
	}

	double complex line3		= (m*m - 1.0)*E*pow(Lhat,3.0)/(2.0*l*(l+1.0)*lambda*eff)*h10;
	if(l == 0 || l == 1){
		line3 			= 0;
	}

	return  line1 + line2 +  line3;
}

double complex ft6_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h9 			= fields->h9;
	double complex h10 		= fields->h10;

	double complex line1 		= E*pow(Lhat,3.0)/( 2.0*l*(l+1.0)*lambda*eff ) * (lambda * h9 + h10);

	if(l == 0) line1 	= 0;
	if(l == 1) line1 	= E*pow(Lhat,3.0)/( 2.0*l*(l+1.0)*eff ) * h9;

	return line1;

}

double complex ft7_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h10 		= fields->h10;

	double complex line1 	= E*pow(Lhat,3.0)/( 2.0*l*(l+1.0)*lambda*eff ) * h10;

	if(l == 0 || l == 1) line1 = 0;

	return line1;
}

double complex Ftm3(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double complex ft7			= ft7_func(l, fields, orbit);
	double complex ft6 		= ft6_func(l, fields, orbit);

	return	 zeta3(l,m)*ft6	+ xi3(l,m)*ft7;
}

double complex Ftm2(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double complex ft1			= ft1_func(l, fields, orbit);
	double complex	ft2			= ft2_func(l, fields, orbit);
	double complex ft3			= ft3_func(l, fields, orbit);

	return alpha2(l,m)*ft1 + beta2(l,m)*ft2 + gamma2(l,m)*ft3;
}

double complex Ftm1(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double complex ft4			= ft4_func(l, fields, orbit);
	double complex ft5			= ft5_func(l, fields, orbit);
	double complex ft6 		= ft6_func(l, fields, orbit);
	double complex ft7			= ft7_func(l, fields, orbit);

	 return epsilon1(l,m)*ft4 + delta1(l,m)*ft5 + zeta1(l,m)*ft6 + xi1(l,m)*ft7;
}

double complex Ft0	(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double complex ft0 		= ft0_func(l, fields, orbit);
	double complex ft1			= ft1_func(l, fields, orbit);
	double complex	ft2			= ft2_func(l, fields, orbit);
	double complex ft3			= ft3_func(l, fields, orbit);

	return ft0 + alpha0(l,m)*ft1 + beta0(l,m)*ft2 + gamma0(l,m)*ft3;
}

double complex Ft1(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double complex ft4			= ft4_func(l, fields, orbit);
	double complex ft5			= ft5_func(l, fields, orbit);
	double complex ft6 		= ft6_func(l, fields, orbit);
	double complex ft7			= ft7_func(l, fields, orbit);

	return epsilonm1(l,m)*ft4 + deltam1(l,m)*ft5 + zetam1(l,m)*ft6 + xim1(l,m)*ft7;
}

double complex Ft2(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double complex ft1			= ft1_func(l, fields, orbit);
	double complex	ft2			= ft2_func(l, fields, orbit);
	double complex ft3			= ft3_func(l, fields, orbit);

	return alpham2(l,m)*ft1 + betam2(l,m)*ft2 + gammam2(l,m)*ft3;
}

double complex Ft3	(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double complex ft7			= ft7_func(l, fields, orbit);
	double complex ft6 		= ft6_func(l, fields, orbit);

	return zetam3(l,m)*ft6 + xim3(l,m)*ft7;
}

double complex fr0_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h1 			= fields->h1;
	double complex h2			= fields->h2;
	double complex h3 			= fields->h3;
	double complex h4 			= fields->h4;
	double complex h5 			= fields->h5;
	double complex h6			= fields->h6;
	double complex h7 			= fields->h7;

	double complex h1_t_deriv	= fields->h1_t_deriv;
	double complex h1_rs_deriv	= fields->h1_rs_deriv;
	double complex h2_t_deriv	= fields->h2_t_deriv;
	double complex h2_rs_deriv	= fields->h2_rs_deriv;
	double complex h3_t_deriv	= fields->h3_t_deriv;
	double complex h3_rs_deriv	= fields->h3_rs_deriv;
	double complex h4_t_deriv	= fields->h4_t_deriv;
	double complex h4_rs_deriv	= fields->h4_rs_deriv;
	double complex h5_t_deriv	= fields->h5_t_deriv;
	double complex h5_rs_deriv	= fields->h5_rs_deriv;
	double complex h6_t_deriv	= fields->h6_t_deriv;
	double complex h6_rs_deriv	= fields->h6_rs_deriv;
	double complex h7_t_deriv	= fields->h7_t_deriv;
	double complex h7_rs_deriv	= fields->h7_rs_deriv;



	double complex line1	= 1.0/(4.0*pow(eff,3.0)) * ( pow(contra_ur,4.0) - I*m*eff*pow(contra_ur,3.0)*Lhat + pow(contra_ur,2.0)*( (r+4.0)/r*E*E + eff) - I*m*contra_ur*eff*Lhat*(E*E + 2.0*eff) - eff*E*E*(1.0-4.0/r) )*h1;

	double complex line2	= - r*contra_ur*E/(4.0*pow(eff,3.0)) * (pow(contra_ur,2.0) + E*E + 2.0*eff)*h1_t_deriv - r/(4.0*pow(eff,3.0)) * ( pow(contra_ur,4.0) + pow(contra_ur,2.0)*(E*E + eff) - eff*E*E )*h1_rs_deriv;

	double complex line3	= E/(2.0*pow(eff,3.0)) * (pow(contra_ur,3.0)*(1.0+1.0/r) - I*m*eff*pow(contra_ur,2.0)*Lhat + contra_ur/r*(E*E + 2.0*eff) - I*m*eff*eff*Lhat )*h2;

	double complex line4	= -r*E*E/(2.0*pow(eff,3.0))*( pow(contra_ur,2.0) + eff)*h2_t_deriv - r*pow(contra_ur,3.0)*E/(2.0*pow(eff,3.0))*h2_rs_deriv;

	double complex line5	= 1.0/(4.0*eff)*( -pow(contra_ur,4.0) + I*m*Lhat*pow(contra_ur,3.0) + E*E*pow(contra_ur,2.0) - (I*m*contra_ur*Lhat + eff)*(E*E - eff) )*h3;

	double complex line6	= r*contra_ur*E/(4.0*eff*eff) * ( pow(contra_ur,2.0) - E*E + eff)*h3_t_deriv + r/(4.0*eff*eff)*(pow(contra_ur,4.0) - E*E*pow(contra_ur,2.0) + eff*(E*E - eff))*h3_rs_deriv;
 
	double complex line7	=  I*m*contra_ur*E*Lhat/(2.0*l*(l+1.0)*eff*eff) * (2.0*contra_ur*(1.0-1.0/r) - I*m*eff*Lhat)*h4 - I*m*r*contra_ur*E*E*Lhat/(2.0*l*(l+1.0)*eff*eff)*h4_t_deriv - I*m*r*E*Lhat/(2.0*l*(l+1.0)*eff*eff)*(pow(contra_ur,2.0) - eff)*h4_rs_deriv;
	if(l == 0){
		line7 = 0;
	}

	double complex line8	= + I*m*Lhat/(2.0*l*(l+1.0)*eff*eff)*(pow(contra_ur,3.0)*(2.0-3.0/r) - I*m*eff*pow(contra_ur,2.0)*Lhat + contra_ur*(E*E/r + 2.0*eff*eff) - I*m*eff*eff*Lhat)*h5;
	if(l == 0){
		line8 = 0;
	}


	double complex line9	= - I*m*r*E*Lhat/(2.0*l*(l+1.0)*eff*eff)*(pow(contra_ur,2.0) + eff)*h5_t_deriv - I*m*r*pow(contra_ur,3.0)*Lhat/(2.0*l*(l+1.0)*eff*eff)*h5_rs_deriv - 0.25*(pow(contra_ur,2.0) - I*m*contra_ur*Lhat + eff)*h6 + r*contra_ur*E/(4.0*eff)*h6_t_deriv;
	if(l == 0){
		line9 = -0.25*(pow(contra_ur,2.0) - I*m*contra_ur*Lhat + eff)*h6 + r*contra_ur*E/(4.0*eff)*h6_t_deriv;
	}

	double complex line10	= r/(4.0*eff)*(pow(contra_ur,2.0) + eff)*h6_rs_deriv + m*Lhat*Lhat/(4.0*l*(l+1.0)*lambda)*(-3.0*m*pow(contra_ur,2.0) + I*contra_ur*Lhat*(4.0 + m*m) -m*eff)*h7 + m*m*r*contra_ur*E*Lhat*Lhat/(4.0*l*(l+1.0)*lambda*eff)*h7_t_deriv;
	if(l == 0 || l == 1){
		line10 = r/(4.0*eff)*(pow(contra_ur,2.0) + eff)*h6_rs_deriv;
	}

	double complex line11 	=  m*m*r*Lhat*Lhat/(4.0*l*(l+1.0)*lambda*eff)*(pow(contra_ur,2.0) - eff)*h7_rs_deriv; 
	if(l == 0 || l == 1){
		line11 = 0;
	}

   

	return line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8 + line9 + line10 + line11;

}

double complex fr1_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h1 			= fields->h1;
	double complex h2			= fields->h2;
	double complex h3 			= fields->h3;
	double complex h5 			= fields->h5;
	double complex h6			= fields->h6;
	double complex h7 			= fields->h7;

	double complex h6_t_deriv	= fields->h6_t_deriv;
	double complex h6_rs_deriv	= fields->h6_rs_deriv;

	double complex line1 = -Lhat*Lhat/(2.0*eff) * ( pow(contra_ur,2.0) + eff)*h1 - contra_ur*E*Lhat*Lhat/(2.0*eff)*h2 + Lhat*Lhat/2.0*( pow(contra_ur,2.0) + eff)*h3 - I*m*contra_ur*pow(Lhat,3.0)/(2.0*l*(l+1.0))*h5;
	if(l == 0){
		line1 = -Lhat*Lhat/(2.0*eff) * ( pow(contra_ur,2.0) + eff)*h1 - contra_ur*E*Lhat*Lhat/(2.0*eff)*h2 + Lhat*Lhat/2.0*( pow(contra_ur,2.0) + eff)*h3;
	}

	double complex line2 =  Lhat*Lhat/4.0 * (3.0*pow(contra_ur,2.0) - I*m*contra_ur*Lhat + eff)*h6 - r*contra_ur*E*Lhat*Lhat/(4.0*eff)*h6_t_deriv - r*Lhat*Lhat/(4.0*eff)*(pow(contra_ur,2.0) - eff)*h6_rs_deriv - I*m*contra_ur*pow(Lhat,3.0)/(l*(l+1.0)*lambda)*h7;
	if(l == 0 || l == 1){
		line2 = Lhat*Lhat/4.0 * (3.0*pow(contra_ur,2.0) - I*m*contra_ur*Lhat + eff)*h6 - r*contra_ur*E*Lhat*Lhat/(4.0*eff)*h6_t_deriv - r*Lhat*Lhat/(4.0*eff)*(pow(contra_ur,2.0) - eff)*h6_rs_deriv;
	}

	return line1 + line2;
}

double complex fr2_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h4 			= fields->h4;
	double complex h5 			= fields->h5;
	double complex h7 			= fields->h7;

	double complex h7_t_deriv	= fields->h7_t_deriv;
	double complex h7_rs_deriv	= fields->h7_rs_deriv;


 

	double complex line1 = - Lhat*Lhat/(2.0*l*(l+1.0)*eff)*(contra_ur*E*h4 + (pow(contra_ur,2.0) + eff)*h5);
	if(l == 0){
		line1 = 0;
	}

	double complex line2 =  Lhat*Lhat/(4.0*l*(l+1.0)*lambda*eff)*(eff*(3.0*pow(contra_ur,2.0) - 5.0*I*m*contra_ur*Lhat + eff)*h7 - r*contra_ur*E*h7_t_deriv - r*(pow(contra_ur,2.0) - eff)*h7_rs_deriv );
	if(l == 0 || l == 1){
		line2 = 0;
	}

	return line1 + line2;
}

double complex fr3_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);;
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h7 			= fields->h7;

	double complex h7_t_deriv	= fields->h7_t_deriv;
	double complex h7_rs_deriv	= fields->h7_rs_deriv;

	double complex line1 =  Lhat*Lhat/(4.0*l*(l+1.0)*lambda*eff) * (-eff* (3.0*pow(contra_ur,2.0) - I*m*contra_ur*Lhat + eff )*h7 + r*contra_ur*E*h7_t_deriv + r*(pow(contra_ur,2.0) - eff)*h7_rs_deriv);

	if(l == 0 || l == 1){
		line1 = 0;
	}

	return line1;
}

double complex fr4_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 			= (l+2.0)*(l-1.0);
	double contra_ur		= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 			= orbit->L/r;
	
	double complex h8 			= fields->h8;
	double complex h9			= fields->h9;
	double complex h10			= fields->h10;

	double complex h10_t_deriv		= fields->h10_t_deriv;
	double complex h10_rs_deriv		= fields->h10_rs_deriv;



 

	double complex line1 		= -I*m*Lhat*Lhat/(2.0*l*(l+1.0)*eff)*(contra_ur*E*h8 + (pow(contra_ur,2.0) + eff )*h9 );
	if(l == 0){
		line1 = 0;
	}

	double complex line2		= I*m*Lhat*Lhat/(2.0*l*(l+1.0)*lambda*eff)*(eff* (3.0*pow(contra_ur,2.0) - 2.0*I*m*contra_ur*Lhat + eff)*h10 - r*contra_ur*E*h10_t_deriv - r*(pow(contra_ur,2.0) - eff)*h10_rs_deriv );
	if(l == 0 || l == 1){
		line2 = 0;
	}

	return line1 + line2;
}

double complex fr5_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double m 			= fields->m;
	double r 			= fields->r;
	double eff			= fields->eff;
	double lambda 		= (l+2.0)*(l-1.0);
	double contra_ur	= fields->contra_ur;
	double E 			= orbit->E;

	double Lhat 		= orbit->L/r;
	
	double complex h8 			= fields->h8;
	double complex h9 			= fields->h9;
	double complex h10 		= fields->h10;

	double complex h8_t_deriv	= fields->h8_t_deriv;
	double complex h8_rs_deriv	= fields->h8_rs_deriv;
	double complex h9_t_deriv	= fields->h9_t_deriv;
	double complex h9_rs_deriv	= fields->h9_rs_deriv;

  

	double complex line1		= - pow(contra_ur,2.0)*E*Lhat/(l*(l+1.0)*eff*eff)*(1.0-1.0/r)*h8  + r*contra_ur*E*E*Lhat/(2.0*l*(l+1.0)*eff*eff)*h8_t_deriv + r*E*Lhat/(2.0*l*(l+1.0)*eff*eff)*(pow(contra_ur,2.0) - eff)*h8_rs_deriv;
	if(l == 0){
		line1 = 0;
	}

	double complex line2		= - contra_ur*Lhat/(2.0*l*(l+1.0)*eff*eff)*( (2.0*r-3.0)/r * pow(contra_ur,2.0) + E*E/r + 2.0*eff*eff )*h9 + r*E*Lhat/(2.0*l*(l+1.0)*eff*eff)*(pow(contra_ur,2.0) + eff)*h9_t_deriv;
	if(l == 0){
		line2 = 0;
	}

	double complex line3		=  r*pow(contra_ur,3.0)*Lhat/(2.0*l*(l+1.0)*eff*eff)*h9_rs_deriv + (m*m - 1.0)*contra_ur*pow(Lhat,3.0)/(2.0*l*(l+1.0)*lambda)*h10;
	if(l == 0){
		line3 = 0;
	}else if( l == 1){
		line3 =  r*pow(contra_ur,3.0)*Lhat/(2.0*l*(l+1.0)*eff*eff)*h9_rs_deriv;
	}

	return  line1 + line2 +  line3;
}

double complex fr6_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double r 			= fields->r;
	double lambda 		= (l+2.0)*(l-1.0);
	double contra_ur	= fields->contra_ur;

	double Lhat 		= orbit->L/r;
	
	double complex h9 			= fields->h9;
	double complex h10 		= fields->h10;

	double complex line1 = contra_ur * pow(Lhat,3.0)/(2.0*l*(l+1.0)*lambda)*(lambda*h9 + h10);
	if(l == 0 ){
		line1 = 0;
	}else if(l == 1){
		line1 = contra_ur*pow(Lhat,3.0)/(2.0*l*(l+1.0)) * h9;
	}

	return line1;

}

double complex fr7_func(double l, struct fields_data *fields, struct orbital_params *orbit)
{
	double r 			= fields->r;
	double lambda 		= (l+2.0)*(l-1.0);
	double contra_ur	= fields->contra_ur;

	double Lhat 		= orbit->L/r;
	
	double complex h10 		= fields->h10;

	double complex line1 = contra_ur * pow(Lhat,3.0)/(2.0*l*(l+1.0)*lambda)*h10;

	if(l == 0 || l == 1){
		line1 = 0;
	}

	return line1;
}

double complex Frm3(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{		
	double complex fr6			= fr6_func(l, fields, orbit);
	double complex fr7			= fr7_func(l, fields, orbit);

	return zeta3(l,m)*fr6 + xi3(l,m)*fr7;
}

double complex Frm2(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{		
	double complex fr1			= fr1_func(l, fields, orbit);
	double complex fr2			= fr2_func(l, fields, orbit);
	double complex fr3 		= fr3_func(l, fields, orbit);

	return alpha2(l,m)*fr1 + beta2(l,m)*fr2 + gamma2(l,m)*fr3;
}

double complex Frm1(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{		
	double complex fr4			= fr4_func(l, fields, orbit);
	double complex fr5 			= fr5_func(l, fields, orbit);
	double complex fr6			= fr6_func(l, fields, orbit);
	double complex fr7			= fr7_func(l, fields, orbit);

	return epsilon1(l,m)*fr4 + delta1(l,m)*fr5 + zeta1(l,m)*fr6 + xi1(l,m)*fr7;
}

double complex Fr0	(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{		
	double complex fr0 		= fr0_func(l, fields, orbit);
	double complex fr1			= fr1_func(l, fields, orbit);
	double complex fr2			= fr2_func(l, fields, orbit);
	double complex fr3 		= fr3_func(l, fields, orbit);

	return fr0 +alpha0(l,m)*fr1 + beta0(l,m)*fr2 + gamma0(l,m)*fr3;
}

double complex Fr1	(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{		
	double complex fr4			= fr4_func(l, fields, orbit);
	double complex fr5 		= fr5_func(l, fields, orbit);
	double complex fr6			= fr6_func(l, fields, orbit);
	double complex fr7			= fr7_func(l, fields, orbit);

	return epsilonm1(l,m)*fr4 + deltam1(l,m)*fr5 + zetam1(l,m)*fr6 + xim1(l,m)*fr7;
}

double complex Fr2	(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{		
	double complex fr1			= fr1_func(l, fields, orbit);
	double complex fr2			= fr2_func(l, fields, orbit);
	double complex fr3 		= fr3_func(l, fields, orbit);

	return alpham2(l,m)*fr1 + betam2(l,m)*fr2 + gammam2(l,m)*fr3;
}

double complex Fr3	(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{		
	double complex fr6			= fr6_func(l, fields, orbit);
	double complex fr7			= fr7_func(l, fields, orbit);

	return zetam3(l,m)*fr6 + xim3(l,m)*fr7;
}

double complex G2 (double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double r_p 	= fields->r;
	double uphi = fields->contra_uphi;
	double complex h3 	= fields->h3;
	double complex h7 	= fields->h7;

	double complex line1 = r_p*r_p*uphi*uphi * alpham2(l,m)*h3;
	double complex line2 = 0;
	if(l > 1){
		line2 = -r_p*r_p * uphi * uphi * factorial(l-2)/factorial(l+2) * (gammam2(l,m) - betam2(l,m)) * h7;
	}	

	return line1 + line2;
}

double complex G1 (double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double r_p 	= fields->r;
	double uphi = fields->contra_uphi;
	double ut 	= fields->contra_ut;
	double ur	= fields->contra_ur;
	double complex h10 = fields->h10;
	double complex h8	= fields->h8;
	double complex h9	= fields->h9;

	double complex line1 = 0;
	double complex line2 = 0;
	if(l > 0){
		line1 = -2.0*r_p*ut*uphi/(l*(l+1.0))*deltam1(l,m)*h8 - 2.0*r_p*ur*uphi/(f(r_p)*l*(l+1.0))*deltam1(l,m)*h9;
		if(l > 1){
			line2 = 2.0*I*m*r_p*r_p*uphi*uphi * factorial(l-2)/factorial(l+2) * epsilonm1(l,m) * h10;
		}
	}
	return line1 + line2;
}

double complex G0 (double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double r_p 	= fields->r;
	double uphi = fields->contra_uphi;
	double ut 	= fields->contra_ut;
	double ur	= fields->contra_ur;
	double eff	= fields->eff;
	double complex h1	= fields->h1;
	double complex h5	= fields->h5;
	double complex h6	= fields->h6;
	double complex h2 	= fields->h2;
	double complex h3	= fields->h3;
	double complex h4 	= fields->h4;
	double complex h7 	= fields->h7;

	double complex line1 =	(h1 + eff*h6)*ut*ut + 2.0/eff*h2*ut*ur + pow(eff,-2.0)*(h1 - eff*h6)*ur*ur + r_p*r_p*uphi*uphi*alpha0(l,m)*h3;
	double complex line2 = 0;
	double complex line3 = 0;
	if(l > 0){
		line2 = 2.0*I*m*r_p*h4/(l*(l+1.0))*ut*uphi + 2.0*I*m*r_p*h5/(eff*l*(l+1.0))*ur*uphi;
		if(l > 1){
			line3 = - factorial(l-2)/factorial(l+2)*(gamma0(l,m) - beta0(l,m) + m*m)*h7 * r_p*r_p*uphi*uphi;
		}
	}

	return line1 + line2 + line3;
}

double complex Gm1 (double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double r_p 	= fields->r;
	double uphi = fields->contra_uphi;
	double ut 	= fields->contra_ut;
	double ur	= fields->contra_ur;
	double eff	= fields->eff;
	double complex h8	= fields->h8;
	double complex h9	= fields->h9;
	double complex h10 = fields->h10;

	double complex line1 = 0;
	double complex line2 = 0;
	if(l > 0){
		line1 = - 2.0 * r_p*h8/(l*(l+1.0))*delta1(l,m)*ut*uphi - 2.0*r_p*h9/(eff*l*(l+1.0))*delta1(l,m)*ur*uphi;
		if(l > 1){
			line2 = 2.0*I*m*r_p*r_p * factorial(l-2)/factorial(l+2) * epsilon1(l,m)*h10*uphi*uphi;
		}
	}
	return line1 + line2;
}

double complex Gm2(double l, double m, struct fields_data *fields, struct orbital_params *orbit)
{
	double r_p 	= fields->r;
	double uphi = fields->contra_uphi;
	double complex h3 	= fields->h3;
	double complex h7 	= fields->h7;

	double complex line1 = r_p*r_p*uphi*uphi*alpha2(l,m)*h3;
	double complex line2 = 0;
	if(l > 1){
		line2 = -r_p*r_p*uphi*uphi*factorial(l-2)/factorial(l+2)*(gamma2(l,m) - beta2(l,m))*h7;
	}
	return line1 + line2;
}

void calculate_lm_mode_contribution_to_Fr(int l_max, struct lm_mode_data *lm_mode, double** Frl_full, struct orbital_params *orbit, int in_or_out)
{
	int l = (int)lm_mode->l;
	double m = lm_mode->m;

	struct coupled_set *cset = lm_mode->cset;

	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES + 1; i++){

		double r_p		= orbit->r_p[i];

		struct fields_data fields;
		if(in_or_out == INNER_FIELDS) cset->construct_inner_fields(lm_mode, &fields, i);
		else if(in_or_out == OUTER_FIELDS) cset->construct_outer_fields(lm_mode, &fields, i);
		fields.m 		= m;
		fields.r 		= r_p;
		fields.eff		= f(r_p);
		fields.contra_ur 	= orbit_contra_ur(orbit->chi_p[i]); 


		double factor 		= 2.0;
		if(m == 0) factor 	= 1.0;
	
		if(l <= l_max)	Frl_full[l][i]			+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Fr0(l,m, &fields, orbit);

		if(l+1 <= l_max) Frl_full[l+1][i]		+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l+1, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Frm1(l,m, &fields, orbit);
		if(l+2 <= l_max) Frl_full[l+2][i]		+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l+2, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Frm2(l,m, &fields, orbit);
		if(l+3 <= l_max) Frl_full[l+3][i]		+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l+3, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Frm3(l,m, &fields, orbit);

		if(l-1 >= 0 && l-1 <= l_max && l-1 >= m) Frl_full[l-1][i]	+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l-1, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Fr1(l,m, &fields, orbit);	
		if(l-2 >= 0 && l-2 <= l_max && l-2 >= m) Frl_full[l-2][i]	+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l-2, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Fr2(l,m, &fields, orbit);
		if(l-3 >= 0 && l-3 <= l_max && l-3 >= m) Frl_full[l-3][i]	+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l-3, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Fr3(l,m, &fields, orbit);	
	
	}

}

void calculate_lm_mode_contribution_to_Ft(int l_max, struct lm_mode_data *lm_mode, double** Ftl_full, struct orbital_params *orbit, int in_or_out)
{
	int l = (int)lm_mode->l;
	double m = lm_mode->m;

	struct coupled_set *cset = lm_mode->cset;

	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES + 1; i++){

		double r_p		= orbit->r_p[i];

		struct fields_data fields;
		if(in_or_out == INNER_FIELDS) cset->construct_inner_fields(lm_mode, &fields, i);
		else if(in_or_out == OUTER_FIELDS) cset->construct_outer_fields(lm_mode, &fields, i);
		fields.m 			= m;
		fields.r 			= r_p;
		fields.eff			= f(r_p);
		fields.contra_ur 	= orbit_contra_ur(orbit->chi_p[i]); 


		double factor 		= 2.0;
		if(m == 0) factor 	= 1.0;
	
		if(l <= l_max)	Ftl_full[l][i]			+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Ft0(l,m, &fields, orbit);

		if(l+1 <= l_max) Ftl_full[l+1][i]		+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l+1, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Ftm1(l,m, &fields, orbit);
		if(l+2 <= l_max) Ftl_full[l+2][i]		+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l+2, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Ftm2(l,m, &fields, orbit);
		if(l+3 <= l_max) Ftl_full[l+3][i]		+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l+3, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Ftm3(l,m, &fields, orbit);

		if(l-1 >= 0 && l-1 <= l_max && l-1 >= m) Ftl_full[l-1][i]	+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l-1, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Ft1(l,m, &fields, orbit);	
		if(l-2 >= 0 && l-2 <= l_max && l-2 >= m) Ftl_full[l-2][i]	+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l-2, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Ft2(l,m, &fields, orbit);
		if(l-3 >= 0 && l-3 <= l_max && l-3 >= m) Ftl_full[l-3][i]	+= factor*pow(r_p,-2.0) * normalized_associated_legendre_P(l-3, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Ft3(l,m, &fields, orbit);	

		//if(i == 0 ) printf("test: %d %.0lg %e\n", l, m, Ftl_full[l][i]);

	}

}

void calculate_lm_mode_contribution_to_huul(double l_max, double **huul_full, struct lm_mode_data *lm_mode, struct orbital_params *orbit, int in_or_out)
{
	int l = (int)lm_mode->l;
	double m = lm_mode->m;

	double factor 		= 2.0;
	if(m == 0) factor 	= 1.0;

	struct coupled_set *cset = lm_mode->cset;

	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES + 1; i++){

		double r_p		= orbit->r_p[i];

		struct fields_data fields;
		if(in_or_out == INNER_FIELDS) cset->construct_inner_fields(lm_mode, &fields, i);
		else if(in_or_out == OUTER_FIELDS) cset->construct_outer_fields(lm_mode, &fields, i);
		fields.m 			= m;
		fields.r 			= r_p;
		fields.eff			= f(r_p);
		fields.contra_ur 	= orbit_contra_ur(orbit->chi_p[i]); 
		fields.contra_ut	= orbit_contra_ut(r_p);
		fields.contra_uphi	= 1.0/(r_p*r_p) * orbit->L;

		if(l <= l_max)	huul_full[l][i]			+= factor/(2.0*r_p) * normalized_associated_legendre_P(l, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * G0(l,m, &fields, orbit);

		if(l+1 <= l_max) huul_full[l+1][i]		+= factor/(2.0*r_p) * normalized_associated_legendre_P(l+1, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Gm1(l,m, &fields, orbit);
		if(l+2 <= l_max) huul_full[l+2][i]		+= factor/(2.0*r_p) * normalized_associated_legendre_P(l+2, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * Gm2(l,m, &fields, orbit);

		if(l-1 >= 0 && l-1 <= l_max && l-1 >= m) huul_full[l-1][i]	+= factor/(2.0*r_p) * normalized_associated_legendre_P(l-1, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * G1(l,m, &fields, orbit);	
		if(l-2 >= 0 && l-2 <= l_max && l-2 >= m) huul_full[l-2][i]	+= factor/(2.0*r_p) * normalized_associated_legendre_P(l-2, m, M_PI/2.0) * cexp(I*m*orbit->phi_p[i]) * G2(l,m, &fields, orbit);

	}


}


// FIXME Below works for circular orbits only

///////////////////////////////////////////////////////////////////////
/// Beta from tensor modes	///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void calculate_tensor_lm_mode_contribution_to_beta(double l_max, double **beta_tensor_full, struct lm_mode_data *lm_mode, struct orbital_params *orbit, int in_or_out)
{
	int l 				= (int)lm_mode->l;
	double m 			= lm_mode->m;
	double omega 	= orbit->Omega_phi;
	double E 			= orbit->E;
	double L			= orbit->L;

	double factor 			= 2.0;
	if(m == 0) factor 	= 1.0;

	struct coupled_set *cset = lm_mode->cset;

	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES + 1; i++){
		double r0		= orbit->r_p[i];
		double f0		= f(r0);
		double df0 	= 2.0/(r0*r0);
		double ut		= orbit_contra_ut(r0);
		double uphi	= 1.0/(r0*r0) * L;
		double omega_p = 1.0/ut;

		struct fields_data fields;
		if(in_or_out == INNER_FIELDS) cset->construct_inner_fields(lm_mode, &fields, i);
		else if(in_or_out == OUTER_FIELDS) cset->construct_outer_fields(lm_mode, &fields, i);
		fields.m 					= m;
		fields.r 					= r0;
		fields.eff				= f0;
		fields.contra_ur 	= orbit_contra_ur(orbit->chi_p[i]); 
		fields.contra_ut	= ut;
		fields.contra_uphi	= uphi;
	
		double complex Y 	= normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		double complex YV1 = 0;
		double complex YV2 = 0;
		double complex YT1 = 0;
		double complex YT2 = 0;
		if(l>0) YV1 = 1.0/(l*(l+1.0)) * normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		if(l>0) YV2 = I*m/(l*(l+1.0)) * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		if(l>1) YT1 = factorial(l-2)/factorial(l+2)*(normalized_associated_legendre_P_theta_2nd_deriv(l,m,M_PI/2.0)* cexp(I*m*orbit->phi_p[i]) + m*m*normalized_associated_legendre_P(l,m, M_PI/2.0)* cexp(I*m*orbit->phi_p[i]));
		if(l>1) YT2 = 2.0*I*m*factorial(l-2)/factorial(l+2)*normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);

		double complex htt				= 1.0/(2.0*r0)*(fields.h1 + f0*fields.h6)*Y;
		double complex hrr 			= 1.0/(2.0*r0*f0*f0)*(fields.h1 - f0*fields.h6)*Y;
		double complex hrphi			= 1.0/(2.0*f0)*(fields.h5*YV2 - fields.h9*YV1);
		double complex htphi			= 1.0/2.0*(fields.h4*YV2 - fields.h8*YV1);
		double complex htr				=	1.0/(2.0*r0*f0)*fields.h2*Y;
		double complex hphiphi		= 1.0/2.0*r0*(fields.h3*Y - fields.h7*YT1 - fields.h10*YT2);

		double complex hrphidt		= -I*m*omega*hrphi;
		double complex htrdt			= -I*m*omega*htr;
		double complex htrdphi		= I*m*htr;
		double complex hrphidphi	= I*m*hrphi;

		double complex httdr			= -htt/r0 + 1.0/(2.0*r0)*(fields.h1_rs_deriv/f0 + fields.h6_rs_deriv + df0*fields.h6)*Y;
		double complex htphidr		= 1.0/(2.0*f0)*(fields.h4_rs_deriv*YV2 - fields.h8_rs_deriv*YV1);
		double complex hphiphidr	= 1.0/r0*hphiphi + 1.0/(2.0*f0)*r0*(fields.h3_rs_deriv*Y - fields.h7_rs_deriv*YT1 - fields.h10_rs_deriv*YT2);

		double complex Grut		= (2.0-r0)*E*hrr/pow(r0,3.0) + E*htrdt + (r0-2.0)*L/(2.0*pow(r0,3.0))*hrphidt -0.5*E*httdr + (2.0-r0)*L/(2.0*pow(r0,3.0))*htphidr + (r0-2.0)*L/(2.0*pow(r0,3.0))*htrdphi;

		double complex Gphiur 	= -E/pow(r0*(r0-2.0),2.0)*htphi - L/pow(r0,5.0)*hphiphi + E/(2.0*r0*r0-4.0*r0)*hrphidt + E/(2.0*r0*r0-4.0*r0)*htphidr + L/(2.0*pow(r0,4.0))*hphiphidr + E/(4.0*r0-2.0*r0*r0)*htrdphi;

		double complex Gruphi	= pow(r0-2.0,2.0)*L/pow(r0,3.0)*hrr + 0.5*E*hrphidt -0.5*E*htphidr + (2.0-r0)*L/(2.0*pow(r0,3.0))*hphiphidr + 0.5*E*htrdphi + (r0-2.0)*L/pow(r0,3.0)*hrphidphi;

		double complex hut			= ut*htt + uphi*htphi;

		double complex huphi		= ut*htphi + uphi*hphiphi;

		double complex F_r			= 0.5*(ut*ut*httdr + 2.0*ut*uphi*htphidr + uphi*uphi*hphiphidr);

		double beta			= creal(-1.0/(2.0*r0*omega)*(Gruphi - pow(omega_p,2.0)*r0*r0*Gphiur + r0*r0*omega/f0*Grut + omega/f0*hut +1.0/(r0*r0)*huphi));

		double beta_hat	= creal(beta - r0*omega_p/2.0*F_r);

	//	double E0 = E;
	//	double L0 = L;
	//	double omegap = omega_p;
	//	double Omega = omega;

	/*	double beta_hat2 = creal(-(E0*hrphidt - E0*htphidr + E0*htrdphi - (hphiphidr*L0*(-2 + r0))/Power(r0,3) + 
       (2*hrphidphi*L0*(-2 + r0))/Power(r0,3) + (2*hrr*L0*Power(-2 + r0,2))/Power(r0,3) + 
       (Omega*((hrphidt - htphidr + htrdphi)*L0*(-2 + r0) + 
            E0*(-2*hrr*(-2 + r0) + (2*htrdt - httdr)*Power(r0,3))))/(-2 + r0) + 
       (Power(omegap,2)*(2*hphiphi*L0*Power(-2 + r0,2) - 
            r0*(hphiphidr*L0*Power(-2 + r0,2) + 
               E0*Power(r0,2)*(-2*htphi + (hrphidt + htphidr - htrdphi)*(-2 + r0)*r0))))/
        (Power(-2 + r0,2)*Power(r0,3)) + (2*(hphiphi*uphi + htphi*ut))/Power(r0,2) + 
       (2*Omega*r0*(htphi*uphi + htt*ut))/(-2 + r0))/(4.*Omega*r0) - 
   (omegap*r0*(hphiphidr*Power(uphi,2) + ut*(2*htphidr*uphi + httdr*ut)))/4.);

		double beta_hat3 = creal((-(hphiphi*Power(-2 + r0,2)) + r0*(-4*hrphidphi + 12*hrr + hphiphidr*(-3 + r0)*(-2 + r0) + 4*htphi*Sqrt(r0) + 
        r0*(4*(hrphidphi - 4*hrr) + 2*(hrphidt - htphi + 3*htphidr - 2*htrdphi)*Sqrt(r0) - 
           (hrphidphi - 7*hrr + htt)*r0 - (hrphidt + 5*htphidr - 4*htrdphi)*Power(r0,1.5) - 
           (hrr - 2*htrdt)*Power(r0,2) + (htphidr - htrdphi)*Power(r0,2.5) - htrdt*Power(r0,3))))/
   (2.*Sqrt(-3 + r0)*(-2 + r0)*Power(r0,2.5)));

		double beta_hat4 = creal((-(hphiphi*Power(-2 + r0,2)) + 12*hrr*r0 + r0*(hphiphidr*(-3 + r0)*(-2 + r0) - double complex(0,1)*hrphi*m*(-3 + r0)*(-2 + r0) + 
        (-2*htphi + htphidr*(-3 + r0)*(-2 + r0) - double complex(0,1)*htr*m*Power(-2 + r0,2) - 
           (htt + hrr*(-7 + r0))*Sqrt(r0))*Power(r0,1.5) + double complex(0,1)*htr*m*(-2 + r0)*Power(r0,1.5) + 
        4*(htphi*Sqrt(r0) - 4*hrr*r0)))/(2.*Sqrt(-3 + r0)*(-2 + r0)*Power(r0,2.5)));*/

	//	if(i==0) printf("dpsi per lm-mode l=%d m=%.0lf : %.12e %.12e %.12e %.12e\n", l, m, beta_hat, beta_hat2, beta_hat3, beta_hat4);
	//	if(i==0) printf("TESTING: %d %.0e %.12e %.12e %.12e\n", l, m, creal(htt), creal(htphi), creal(hphiphi));

	//	if(i==0) printf("hpp = %.12e\nhrr = %.12e\nhrp,p = %.12e\nhpp,r = %.12e\nhtp = %.12e\nhrp,t = %.12e\nhtr,p = %.12e\nhtp,r = %.12e\nhtt = %.12e\nhtr,t = %.12e\n", creal(hphiphi), creal(hrr), creal(hrphidphi), creal(hphiphidr), creal(htphi), creal(hrphidt), creal(htrdphi), creal(htphidr), creal(htt), creal(htrdt));
	
		//if(i==0) printf("test: %.12e %.12e\n", creal(hrphi), creal(hrphidphi));

		if(l<= l_max) beta_tensor_full[l][i] 	+= factor * beta_hat;

	}



}



void calculate_tensor_lm_mode_contribution_to_huul(double l_max, double **huul_full, struct lm_mode_data *lm_mode, struct orbital_params *orbit, int in_or_out)
{
	int l = (int)lm_mode->l;
	double m = lm_mode->m;

	double factor 		= 2.0;
	if(m == 0) factor 	= 1.0;

	struct coupled_set *cset = lm_mode->cset;

	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES + 1; i++){
		double r_p		= orbit->r_p[i];

		struct fields_data fields;
		if(in_or_out == INNER_FIELDS) cset->construct_inner_fields(lm_mode, &fields, i);
		else if(in_or_out == OUTER_FIELDS) cset->construct_outer_fields(lm_mode, &fields, i);
		fields.m 			= m;
		fields.r 			= r_p;
		fields.eff			= f(r_p);
		fields.contra_ur 	= orbit_contra_ur(orbit->chi_p[i]); 
		fields.contra_ut	= orbit_contra_ut(r_p);
		fields.contra_uphi	= 1.0/(r_p*r_p) * orbit->L;
	
		double complex Y 	= normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		double complex YV1, YV2, YT1, YT2;
		if(l>0) YV1 = 1.0/(l*(l+1.0)) * normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		if(l>0) YV2 = I*m/(l*(l+1.0)) * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		if(l>1) YT1 = factorial(l-2)/factorial(l+2)*(normalized_associated_legendre_P_theta_2nd_deriv(l,m,M_PI/2.0)* cexp(I*m*orbit->phi_p[i]) + m*m*normalized_associated_legendre_P(l,m, M_PI/2.0)* cexp(I*m*orbit->phi_p[i]));
		if(l>1) YT2 = 2.0*I*m*factorial(l-2)/factorial(l+2)*normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);

		if(l<= l_max) huul_full[l][i] 	+= factor * creal(1.0/(2.0*r_p)*(fields.h1 + f(r_p)*fields.h6)*Y*pow(fields.contra_ut,2.0) + (fields.h4*YV2 - fields.h8*YV1) *fields.contra_ut * fields.contra_uphi + 0.5*r_p *(fields.h3*Y - fields.h7*YT1 - fields.h10 * YT2) * pow(fields.contra_uphi,2.0));

	//	if(l==2 && m == 0) printf("TESTING huu: %.16e\n", huul_full[l][i]);
		
	}
}


void calculate_tensor_lm_mode_contribution_to_hR_and_FrR(double l_max, double **huul_R, double **FrlR, struct n_mode_data *n_mode, struct orbital_params *orbit)
{	
	int l 		= n_mode->l;
	int m 		= n_mode->m;

	double factor 		= 2.0;
	if(m == 0) factor 	= 1.0;
	
	double r0 	= orbit->p;
	double phi0	= 0;
	
	double contra_ut 	= orbit_contra_ut(r0);
	double contra_uphi	= 1.0/(r0*r0) * orbit->L;
	
	double complex Y 	= normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(I*m*phi0);
	double complex YV1, YV2, YT1, YT2;
	if(l>0) YV1 = 1.0/(l*(l+1.0)) * normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(I*m*phi0);
	if(l>0) YV2 = I*m/(l*(l+1.0)) * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(I*m*phi0);
	if(l>1) YT1 = factorial(l-2)/factorial(l+2)*(normalized_associated_legendre_P_theta_2nd_deriv(l,m,M_PI/2.0)* cexp(I*m*phi0) + m*m*normalized_associated_legendre_P(l,m, M_PI/2.0)* cexp(I*m*phi0));
	if(l>1) YT2 = 2.0*I*m*factorial(l-2)/factorial(l+2)*normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(I*m*phi0);
	
	double complex hP_gauge1 	= 0;
	double complex hP_gauge2 	= 0;
	double complex dhP_gauge1 	= 0;
	double complex dhP_gauge2 	= 0;
	
	int i,j;
	for(i = 0; i < orbit->gridsize; i++){
		
		if(i==orbit->r0_grid_index){
			//printf("gauge test: %d %d %.14e %.14e\n", l, m, creal(orbit->hP_gauge[2*(l*(l+1)+2*m)]), creal(orbit->hP_gauge[8]));
			hP_gauge1 	= orbit->hP_gauge[2*l*(l+1) + 4*m];
			hP_gauge2 	= orbit->hP_gauge[2*l*(l+1) + 4*m + 1];
			dhP_gauge1 	= orbit->hP_gauge[2*l*(l+1) + 4*m + 2];
			dhP_gauge2 	= orbit->hP_gauge[2*l*(l+1) + 4*m + 3];
		}
		
		double complex *hR = calloc(10, sizeof(double complex));
		double complex *dhR = calloc(10, sizeof(double complex));
		for(j = 0; j < 10; j++) hR[j] = 0;
		if((l+m) % 2  == 0 ){
		 	if(m==0){
				hR[0] = n_mode->hR[0][0][i] + I*n_mode->hR[0][1][i];
				hR[2] = n_mode->hR[1][0][i] + I*n_mode->hR[1][1][i];
				hR[4] = n_mode->hR[2][0][i] + I*n_mode->hR[2][1][i];
				hR[5] = n_mode->hR[3][0][i] + I*n_mode->hR[3][1][i] + hP_gauge1;
				hR[6] = n_mode->hR[4][0][i] + I*n_mode->hR[4][1][i] + hP_gauge2 + l*(l+1)*hP_gauge1;
				
				dhR[0] = n_mode->hR[0][2][i] + I*n_mode->hR[0][3][i];
				dhR[2] = n_mode->hR[1][2][i] + I*n_mode->hR[1][3][i];
				dhR[4] = n_mode->hR[2][2][i] + I*n_mode->hR[2][3][i];
				dhR[5] = n_mode->hR[3][2][i] + I*n_mode->hR[3][3][i] + (24.0 - 16.0*r0 + 2.0*r0*r0)/(2.0*r0*pow(r0-2.0, 2.0))*hP_gauge1 + dhP_gauge1;
				dhR[6] = n_mode->hR[4][2][i] + I*n_mode->hR[4][3][i] + dhP_gauge2 + l*(l+1)*((24.0 - 16.0*r0 + 2.0*r0*r0)/(2.0*r0*pow(r0-2.0, 2.0))*hP_gauge1 +dhP_gauge1); // + gauge correction
			}else{
				hR[0] = n_mode->hR[0][0][i] + I*n_mode->hR[0][1][i];
				hR[2] = n_mode->hR[1][0][i] + I*n_mode->hR[1][1][i];
				hR[4] = n_mode->hR[2][0][i] + I*n_mode->hR[2][0][i];
				hR[5] = n_mode->hR[3][0][i] + I*n_mode->hR[3][1][i];
				
				dhR[0] = n_mode->hR[0][2][i] + I*n_mode->hR[0][3][i];
				dhR[2] = n_mode->hR[1][2][i] + I*n_mode->hR[1][3][i];
				if(l==1){
					hR[1] = n_mode->hR[4][0][i] + I*n_mode->hR[4][1][i] + hP_gauge1;
					hR[3] = n_mode->hR[5][0][i] + I*n_mode->hR[5][1][i] + hP_gauge2;
					
					dhR[1] = n_mode->hR[4][2][i] + I*n_mode->hR[4][3][i] + dhP_gauge1;
					dhR[3] = n_mode->hR[5][2][i] + I*n_mode->hR[5][3][i] + dhP_gauge2;
					dhR[5] = n_mode->hR[3][2][i] + I*n_mode->hR[3][3][i];
				}
				else{
					hR[1] = n_mode->hR[5][0][i] + I*n_mode->hR[5][1][i] + hP_gauge1;
					hR[3] = n_mode->hR[6][0][i] + I*n_mode->hR[6][1][i] + hP_gauge2;
					hR[6] = n_mode->hR[4][0][i] + I*n_mode->hR[4][1][i];
					
					dhR[1] = n_mode->hR[5][2][i] + I*n_mode->hR[5][3][i] + dhP_gauge1;
					dhR[3] = n_mode->hR[6][2][i] + I*n_mode->hR[6][3][i] + dhP_gauge2;
					dhR[5] = n_mode->hR[3][2][i] + I*n_mode->hR[3][3][i];
					dhR[6] = n_mode->hR[4][2][i] + I*n_mode->hR[4][3][i];
				}	
			}
		}else{
			if(m==0){
				//if(i== orbit->r0_grid_index) printf("new test: %.14e\n", n_mode->hR[0][0][i]);
				hR[7] = n_mode->hR[0][0][i] + I*n_mode->hR[0][1][i];
				
				dhR[7] = n_mode->hR[0][2][i] + I*n_mode->hR[0][3][i];
			}else{
				hR[7] = n_mode->hR[2][0][i] + I*n_mode->hR[2][1][i] + hP_gauge1;
				hR[8] = n_mode->hR[0][0][i] + I*n_mode->hR[0][1][i];
				hR[9] = n_mode->hR[1][0][i] + I*n_mode->hR[1][1][i]; 
				
				dhR[7] = n_mode->hR[2][2][i] + I*n_mode->hR[2][3][i] + dhP_gauge1; 
				dhR[8] = n_mode->hR[0][2][i] + I*n_mode->hR[0][3][i];
				dhR[9] = n_mode->hR[1][2][i] + I*n_mode->hR[1][3][i];
			}
		}
		
		
		if(l <= l_max){
			huul_R[l][i] += factor * creal(1.0/(2.0*r0)*(hR[0] + f(r0)*hR[5])*Y*pow(contra_ut,2.0) + (hR[3]*YV2 - hR[7]*YV1) *contra_ut * contra_uphi + 0.5*r0 *(hR[2]*Y - hR[6]*YT1 - hR[9] * YT2) * pow(contra_uphi,2.0));
			//if(i==orbit->r0_grid_index) printf("test: %.14e %.14e %.14e\n", creal(hR[3]), cimag(hR[3]), cimag(hP_gauge1));
			 
			double complex dhtt_R =  -1.0/(2.0*r0*r0*r0)*(r0*hR[0] + (r0-4.0)*hR[5] - r0*(r0*dhR[0]+(r0-2.0)*dhR[5]))*Y; 
			double complex dhpp_R =  1.0/2.0*(hR[2]*Y - hR[6]*YT1 - hR[9]*YT2) + r0/2.0*(dhR[2]*Y - dhR[6]*YT1 - dhR[9]*YT2);
			double complex dhtp_R =  1.0/2.0*(dhR[3]*YV2 - dhR[7]*YV1);
			
			/*if(i==orbit->r0_grid_index) {
				printf("test %.14e %.14e\n", creal(hR[2]), cimag(hR[2]));
				printf("test %.14e %.14e\n", creal(hR[5]), cimag(hR[5]));
				printf("test %.14e %.14e\n", creal(dhR[2]), cimag(dhR[2]));
				printf("test %.14e %.14e\n", creal(dhR[5]), cimag(dhR[5]));
				printf("test %.14e %.14e\n", creal(dhR[6]), cimag(dhR[6]));
				printf("test %.14e %.14e\n", creal(dhtt_R*pow(contra_ut, 2.0) + dhtp_R*contra_ut*contra_uphi + dhpp_R*pow(contra_uphi, 2.0)), cimag(dhtp_R));
			}*/
			
			FrlR[l][i] += 1.0/2.0*factor*creal(dhtt_R*pow(contra_ut, 2.0) + 2.0*dhtp_R*contra_ut*contra_uphi + dhpp_R*pow(contra_uphi, 2.0));
		 }
		 
		 free(hR);
		 free(dhR);
	 }
	
	
	
	
}

