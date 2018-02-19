#include <gsl/gsl_errno.h>
#include "GSF_ecc.h"

int monopole_ode_system(double rs, const double y[], double g[], void *params)
{
	struct n_mode_data *n_mode = (struct n_mode_data *)params;

	double omega 	= n_mode->omega; 

	double r 		= r_of_rs(rs);
	double eff 		= f(r);
	double df		= 2.0/(r*r);

	double V_l 		= 0.25*eff*(df/r);

	//R1 ode
	g[0] = y[1];
	g[1] = (4.0 * V_l  - omega*omega) * y[0] + 2.0*eff*df*y[5] + 2.0*eff/(r*r)*(1.0-4.0/r)*(y[0] - eff*y[4]) - 2.0*eff*eff/(r*r)*(1.0- 6.0/r) * y[8] ;
	g[2] = y[3];
	g[3] = (4.0 * V_l  - omega*omega) * y[2] + 2.0*eff*df*y[7] + 2.0*eff/(r*r)*(1.0-4.0/r)*(y[2] - eff*y[6]) - 2.0*eff*eff/(r*r)*(1.0- 6.0/r) * y[10] ;

	//R3 ode
	g[4] = y[5];
	g[5] = (4.0 * V_l  - omega*omega) * y[4] - 2.0*eff/(r*r) * ( y[0] - (1.0 - 4.0/r)*(y[4] + y[8]) );
	g[6] = y[7];
	g[7] = (4.0 * V_l  - omega*omega) * y[6] - 2.0*eff/(r*r) * ( y[2] - (1.0 - 4.0/r)*(y[6] + y[10]) );

	//R6 ode
	g[8] = y[9];
	g[9] = (4.0 * V_l - omega*omega) * y[8] - 2.0*eff/(r*r) * (y[0] - (1.0 - 4.0/r)*(y[4] + y[8]));
	g[10] = y[11];
	g[11] = (4.0 * V_l - omega*omega) * y[10] - 2.0*eff/(r*r) * (y[2]- (1.0 - 4.0/r)*(y[6] + y[10]));


	return GSL_SUCCESS;

}


int even_dipole_ode_system(double rs, const double y[], double g[], void *params)
{

	struct n_mode_data *n_mode = (struct n_mode_data *)params;

	double l 		= n_mode->l;
	double omega 	= n_mode->omega; 

	double r 		= r_of_rs(rs);
	double eff 		= f(r);
	double df		= 2.0/(r*r);

	double V_l 		= 0.25*eff*(df/r + l*(l+1.0)/(r*r));

	//R1 ode
	g[0] = y[1];
	g[1] = (4.0 * V_l  - omega*omega) * y[0] + 2.0*eff*df*y[5] + 2.0*eff/(r*r)*(1.0-4.0/r)*(y[0] - y[8]  - eff*y[4]) - 2.0*eff*eff/(r*r)*(1.0- 6.0/r) * y[12] ;
	g[2] = y[3];
	g[3] = (4.0 * V_l  - omega*omega) * y[2] + 2.0*eff*df*y[7] + 2.0*eff/(r*r)*(1.0-4.0/r)*(y[2] - y[10] - eff*y[6]) - 2.0*eff*eff/(r*r)*(1.0- 6.0/r) * y[14] ;

	//R3 ode
	g[4] = y[5];
	g[5] = (4.0 * V_l  - omega*omega) * y[4] - 2.0*eff/(r*r) * ( y[0] - y[8]  - (1.0 - 4.0/r)*(y[4] + y[12]) );
	g[6] = y[7];
	g[7] = (4.0 * V_l  - omega*omega) * y[6] - 2.0*eff/(r*r) * ( y[2] - y[10] - (1.0 - 4.0/r)*(y[6] + y[14]) );

	//R5 ode
	g[8]  = y[9];
	g[9]  = (4.0 * V_l  - omega*omega) * y[8]  + 4.0*eff/(r*r) * ( (1.0-4.5/r)*y[8]  - 0.5*l*(l+1.0)*(y[0] - eff*y[4]) + 0.5*(1.0 - 3.0/r)*(l*(l+1.0)*y[12]));
	g[10] = y[11];
	g[11] = (4.0 * V_l  - omega*omega) * y[10] + 4.0*eff/(r*r) * ( (1.0-4.5/r)*y[10] - 0.5*l*(l+1.0)*(y[2] - eff*y[6]) + 0.5*(1.0 - 3.0/r)*(l*(l+1.0)*y[14] ));

	//R6 ode
	g[12] = y[13];
	g[13] = (4.0 * V_l - omega*omega) * y[12] - 2.0*eff/(r*r) * (y[0] - y[8]  - (1.0 - 4.0/r)*(y[4] + y[12]));
	g[14] = y[15];
	g[15] = (4.0 * V_l - omega*omega) * y[14] - 2.0*eff/(r*r) * (y[2] - y[10] - (1.0 - 4.0/r)*(y[6] + y[14]));



	return GSL_SUCCESS;


}

int even_ode_system(double rs, const double y[], double g[], void *params)
{
	struct n_mode_data *n_mode = (struct n_mode_data *)params;

	double l 		= n_mode->l;
	double omega 	= n_mode->omega; 

	double r 		= r_of_rs(rs);
	double eff 		= f(r);
	double df		= 2.0/(r*r);

	double V_l 		= 0.25*eff*(df/r + l*(l+1.0)/(r*r));

	//R1 ode
	g[0] = y[1];
	g[1] = (4.0 * V_l  - omega*omega) * y[0] + 2.0*eff*df*y[5] + 2.0*eff/(r*r)*(1.0-4.0/r)*(y[0] - y[8] - eff*y[4]) - 2.0*eff*eff/(r*r)*(1.0- 6.0/r) * y[12] ;
	g[2] = y[3];
	g[3] = (4.0 * V_l  - omega*omega) * y[2] + 2.0*eff*df*y[7] + 2.0*eff/(r*r)*(1.0-4.0/r)*(y[2] - y[10] - eff*y[6]) - 2.0*eff*eff/(r*r)*(1.0- 6.0/r) * y[14] ;

	//R3 ode
	g[4] = y[5];
	g[5] = (4.0 * V_l  - omega*omega) * y[4] - 2.0*eff/(r*r) * ( y[0] - y[8] - (1.0 - 4.0/r)*(y[4] + y[12]) );
	g[6] = y[7];
	g[7] = (4.0 * V_l  - omega*omega) * y[6] - 2.0*eff/(r*r) * ( y[2] - y[10] - (1.0 - 4.0/r)*(y[6] + y[14]) );

	//R5 ode
	g[8] = y[9];
	g[9] = (4.0 * V_l  - omega*omega) * y[8] + 4.0*eff/(r*r) * ( (1.0-4.5/r)*y[8] - 0.5*l*(l+1.0)*(y[0] - eff*y[4]) + 0.5*(1.0 - 3.0/r)*(l*(l+1.0)*y[12] - y[16]));
	g[10] = y[11];
	g[11] = (4.0 * V_l  - omega*omega) * y[10] + 4.0*eff/(r*r) * ( (1.0-4.5/r)*y[10] - 0.5*l*(l+1.0)*(y[2] - eff*y[6]) + 0.5*(1.0 - 3.0/r)*(l*(l+1.0)*y[14] - y[18]));

	//R6 ode
	g[12] = y[13];
	g[13] = (4.0 * V_l - omega*omega) * y[12] - 2.0*eff/(r*r) * (y[0] - y[8] - (1.0 - 4.0/r)*(y[4] + y[12]));
	g[14] = y[15];
	g[15] = (4.0 * V_l - omega*omega) * y[14] - 2.0*eff/(r*r) * (y[2] - y[10] - (1.0 - 4.0/r)*(y[6] + y[14]));

	//R7 ode
	g[16] = y[17];
	g[17] = (4.0 * V_l - omega*omega) * y[16] - 2.0*eff/(r*r) * (y[16] + (l+2.0)*(l-1.0) * y[8]);
	g[18] = y[19];
	g[19] = (4.0 * V_l - omega*omega) * y[18] - 2.0*eff/(r*r) * (y[18] + (l+2.0)*(l-1.0) * y[10]);

	return GSL_SUCCESS;


}

int even_static_ode_system(double rs, const double y[], double g[], void *params)
{
	struct n_mode_data *n_mode = (struct n_mode_data *)params;

	double l 		= n_mode->l;
	double omega 	= n_mode->omega; 

	double r 		= r_of_rs(rs);
	double eff 		= f(r);
	double df		= 2.0/(r*r);

	double V_l 		= 0.25*eff*(df/r + l*(l+1.0)/(r*r));

	//R1 ode
	g[0] = y[1];
	g[1] = (4.0 * V_l  - omega*omega) * y[0] + 2.0*eff*df*y[5] + 2.0*eff/(r*r)*(1.0-4.0/r)*(y[0] - y[8] - eff*y[4]) - 2.0*eff*eff/(r*r)*(1.0- 6.0/r) * 1.0/(2.0*eff)*(r/eff*y[1] - r*y[5] + y[0] - y[8]  - eff*y[4]) ;
	g[2] = y[3];
	g[3] = (4.0 * V_l  - omega*omega) * y[2] + 2.0*eff*df*y[7] + 2.0*eff/(r*r)*(1.0-4.0/r)*(y[2] - y[10] - eff*y[6]) - 2.0*eff*eff/(r*r)*(1.0- 6.0/r) * 1.0/(2.0*eff)*(r/eff*y[3] - r*y[7] + y[2] - y[10] - eff*y[6]) ;

	//R3 ode
	g[4] = y[5];
	g[5] = (4.0 * V_l  - omega*omega) * y[4] - 2.0*eff/(r*r) * ( y[0] - y[8]  - (1.0 - 4.0/r)*(y[4] + 1.0/(2.0*eff)*(r/eff*y[1] - r*y[5] + y[0] - y[8]  - eff*y[4])) );
	g[6] = y[7];
	g[7] = (4.0 * V_l  - omega*omega) * y[6] - 2.0*eff/(r*r) * ( y[2] - y[10] - (1.0 - 4.0/r)*(y[6] + 1.0/(2.0*eff)*(r/eff*y[3] - r*y[7] + y[2] - y[10] - eff*y[6])) );

	//R5 ode
	g[8] = y[9];
	g[9] = (4.0 * V_l  - omega*omega) * y[8] + 4.0*eff/(r*r) * ( (1.0-4.5/r)*y[8] - 0.5*l*(l+1.0)*(y[0] - eff*y[4]) - 0.5*(1.0 - 3.0/r)*(r/eff*y[9] + 2.0*y[8]));
	g[10] = y[11];
	g[11] = (4.0 * V_l  - omega*omega) * y[10] + 4.0*eff/(r*r) * ( (1.0-4.5/r)*y[10] - 0.5*l*(l+1.0)*(y[2] - eff*y[6]) - 0.5*(1.0 - 3.0/r)*(r/eff*y[11] + 2.0*y[10]));

	return GSL_SUCCESS;

}

int odd_dipole_ode_system (double rs, const double y[], double g[], void *params)
{

	struct n_mode_data *n_mode = (struct n_mode_data *)params;

	double l 		= n_mode->l;
	double omega 	= n_mode->omega; 

	double r 		= r_of_rs(rs);
	double eff 		= f(r);

	double df		= 2.0/(r*r);

	double V_l 		= 0.25*eff*(df/r + l*(l+1.0)/(r*r));


	// R9 ode
	g[0] = y[1];
	g[1] = (4.0*V_l + 4.0*eff/(r*r)*(1.0-4.5/r) - omega*omega) * y[0];
	g[2] = y[3];
	g[3] = (4.0*V_l + 4.0*eff/(r*r)*(1.0-4.5/r) - omega*omega) * y[2];

	return GSL_SUCCESS;

}

int odd_ode_system (double rs, const double y[], double g[], void *params)
{

	struct n_mode_data *n_mode = (struct n_mode_data *)params;

	double l 		= n_mode->l;
	double omega 	= n_mode->omega; 

	double r 		= r_of_rs(rs);
	double eff 		= f(r);

	double df		= 2.0/(r*r);

	double V_l 		= 0.25*eff*(df/r + l*(l+1.0)/(r*r));

	double lambda = (l+2.0)*(l-1.0);

	// R9 ode
	g[0] = y[1];
	g[1] = (4.0*V_l + 4.0*eff/(r*r)*(1.0-9.0/(2.0*r)) - omega*omega) * y[0] - 2.0*eff/(r*r)*(1.0-3.0/r)*y[4];
	g[2] = y[3];
	g[3] = (4.0*V_l + 4.0*eff/(r*r)*(1.0-9.0/(2.0*r)) - omega*omega) * y[2] - 2.0*eff/(r*r)*(1.0-3.0/r)*y[6];

	// R10 ode
	g[4] = y[5];
	g[5] = (4.0*V_l - omega*omega - 2.0*eff/(r*r))*y[4] - 2.0*eff/(r*r)*lambda*y[0];
	g[6] = y[7];
	g[7] = (4.0*V_l - omega*omega - 2.0*eff/(r*r))*y[6] - 2.0*eff/(r*r)*lambda*y[2];

	return GSL_SUCCESS;

}
