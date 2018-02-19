#include "sources.h"
#include "math.h"
#include "stdio.h"


double complex source_sin(double chi, struct n_mode_data *n_mode)
{
	double omega 		= n_mode->omega;
	double m			= n_mode->m;
	
	double t_p 			= t_of_chi(chi);
	double phi_p 		= phi_of_chi(chi);

	return sin(omega*t_p - m*phi_p);
}

double complex source_cos(double chi, struct n_mode_data *n_mode)
{
	double omega 		= n_mode->omega;
	double m			= n_mode->m;
	
	double t_p 			= t_of_chi(chi);
	double phi_p 		= phi_of_chi(chi);

	return cos(omega*t_p - m*phi_p);
}

double complex R1_source_coeff(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	double r_p 			= r_of_chi(chi);
	double f_p			= f(r_p);

	double R1_source_coeff = - 4.0 * 8.0 * M_PI * normalized_associated_legendre_P(n_mode->l, n_mode->m, M_PI/2.0) ;

	return  R1_source_coeff * f_p*f_p * (2.0*pow(orbit->E,2.0)*r_p*r_p - f_p*r_p*r_p - pow(orbit->L,2.0) * f_p) / (orbit->E * pow(r_p, 3.0)) ;
}

double complex R1_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece)
{
	double complex return_value = 0;

	if(source_piece == SOURCE_TOTAL) 		return_value = R1_source_coeff(chi, n_mode, orbit) * source_cos(chi, n_mode);
	else if(source_piece == SOURCE_COEFF) 	return_value = R1_source_coeff(chi, n_mode, orbit);
	else if(source_piece == SOURCE_OSCIL) 	return_value = source_cos(chi, n_mode);

//	printf("test: %e\n", creal(return_value));

	return return_value;
}

double complex R3_source_coeff(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	double r_p 			= r_of_chi(chi);
	double f_p			= f(r_p);

	double R3_source_coeff = - 4.0 * 8.0 * M_PI * normalized_associated_legendre_P(n_mode->l, n_mode->m, M_PI/2.0);

	return R3_source_coeff * f_p*f_p * (r_p*r_p + pow(orbit->L, 2.0))/(orbit->E * pow(r_p, 3.0));

}

double complex R3_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece)
{
	double complex return_value = 0;

	if(source_piece == SOURCE_TOTAL) 		return_value = R3_source_coeff(chi, n_mode, orbit) * source_cos(chi, n_mode);
	else if(source_piece == SOURCE_COEFF) 	return_value = R3_source_coeff(chi, n_mode, orbit);
	else if(source_piece == SOURCE_OSCIL) 	return_value = source_cos(chi, n_mode);

	return return_value;
}

double complex R5_source_coeff(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	double m			= n_mode->m;
	
	double r_p 			= r_of_chi(chi);
	double contra_ur	= orbit_contra_ur(chi);
	double f_p			= f(r_p);

	double h5_source_coeff = - 4.0 * 16.0 * M_PI * m * orbit->L /(orbit->E) * normalized_associated_legendre_P(n_mode->l, n_mode->m, M_PI/2.0);

	return h5_source_coeff  * f_p*f_p * contra_ur / pow(r_p, 2.0);
}

double complex R5_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece)
{
	double complex return_value = 0;

	if(source_piece == SOURCE_TOTAL) 		return_value = R5_source_coeff(chi, n_mode, orbit) * source_sin(chi, n_mode);
	else if(source_piece == SOURCE_COEFF) 	return_value = R5_source_coeff(chi, n_mode, orbit);
	else if(source_piece == SOURCE_OSCIL) 	return_value = source_sin(chi, n_mode);

	return return_value;
}


double complex R6_source_coeff(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	double r_p 			= r_of_chi(chi);
	double f_p			= f(r_p);

	double h6_source_coeff = - 4.0 * 8.0 * M_PI * pow(orbit->L,2.0) /(orbit->E) * normalized_associated_legendre_P(n_mode->l, n_mode->m, M_PI/2.0);

	return h6_source_coeff * f_p*f_p / pow(r_p,3.0);

}

double complex R6_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece)
{
	double complex return_value = 0;

	if(source_piece == SOURCE_TOTAL) 		return_value = R6_source_coeff(chi, n_mode, orbit) * source_cos(chi, n_mode);
	else if(source_piece == SOURCE_COEFF) 	return_value = R6_source_coeff(chi, n_mode, orbit);
	else if(source_piece == SOURCE_OSCIL) 	return_value = source_cos(chi, n_mode);

	return return_value;
}

double complex R7_source_coeff(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	double l = n_mode->l;
	double m = n_mode->m;

	return ( l*(l+1.0) - 2.0*m*m) * R6_source_coeff(chi, n_mode, orbit);
}

double complex R7_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece)
{
	double l = n_mode->l;
	double m = n_mode->m;

	double complex return_value = 0;

	if(source_piece == SOURCE_TOTAL) 		return_value = ( l*(l+1.0) - 2.0*m*m) * R6_source_coeff(chi, n_mode, orbit)*source_cos(chi, n_mode);
	else if(source_piece == SOURCE_COEFF) 	return_value = ( l*(l+1.0) - 2.0*m*m) * R6_source_coeff(chi, n_mode, orbit);
	else if(source_piece == SOURCE_OSCIL) 	return_value = source_cos(chi, n_mode);

	return return_value;
}

///////////////////////// Odd sector sources below //////////////////////////////

double complex R9_source_coeff(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	double r_p 			= r_of_chi(chi);
	double contra_ur	= orbit_contra_ur(chi);
	double f_p			= f(r_p);

	double h9_source_coeff = 4.0 * 16.0 * M_PI * orbit->L / (orbit->E) * normalized_associated_legendre_P_theta_deriv(n_mode->l, n_mode->m, M_PI/2.0);

	return h9_source_coeff * I * f_p * f_p * contra_ur /(r_p*r_p);
}



double complex R9_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece)
{
	double complex return_value = 0;

	if(source_piece == SOURCE_TOTAL) 		return_value = R9_source_coeff(chi, n_mode, orbit) * source_sin(chi, n_mode);
	else if(source_piece == SOURCE_COEFF) 	return_value = R9_source_coeff(chi, n_mode, orbit);
	else if(source_piece == SOURCE_OSCIL) 	return_value = source_sin(chi, n_mode);

	return return_value;
}

double complex R10_source_coeff(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	double r_p 			= r_of_chi(chi);
	double f_p			= f(r_p);

	double h10_source_coeff = -4.0 * 16.0 * M_PI * n_mode->m * pow(orbit->L, 2.0) / (orbit->E) * normalized_associated_legendre_P_theta_deriv(n_mode->l, n_mode->m, M_PI/2.0);

	return h10_source_coeff * I * f_p*f_p/(r_p*r_p*r_p);
}

double complex R10_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece)
{
	double complex return_value = 0;

	if(source_piece == SOURCE_TOTAL) 		return_value = R10_source_coeff(chi, n_mode, orbit) * source_cos(chi, n_mode);
	else if(source_piece == SOURCE_COEFF) 	return_value = R10_source_coeff(chi, n_mode, orbit);
	else if(source_piece == SOURCE_OSCIL) 	return_value = source_cos(chi, n_mode);

	return return_value;
}
