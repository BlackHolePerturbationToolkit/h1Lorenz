/** \file libSch.c
* Library for calculating various properties of orbits in Schwarzschild spacetime
*/

#include "libSch.h"
#include <stdio.h>
#include "math_funcs.h"
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_lambert.h>

extern const double DESIRED_INVERT_ACCURACY;

extern double e, p;
extern double M;

double R8_out_odd_static_n_zero(double r_chi, int l);
double R8_out_odd_static_n_zero_r_deriv(double r_chi, int l);

void	set_primary_orbital_params(double p_new, double e_new)
{
	p = p_new;
	e = e_new;
}
void 	get_orbit_params(struct orbital_params *orbit, int num_chi_values)
{
	orbit->e 			= e;
	orbit->p			= p;

	orbit->E 			= orbit_energy();
	orbit->L 			= orbit_ang_mom();
	orbit->r_max 		= r_max();
	orbit->r_min 		= r_min();
	orbit->T_r 			= orbit_radial_period();
	orbit->Omega_phi	= orbit_phi_freq();
	orbit->Omega_r		= orbit_radial_freq();
	orbit->rs_max 		= rs_of_r(r_max());
	orbit->rs_min		= rs_of_r(r_min());

	orbit->chi_p 		= (double *)malloc((2*num_chi_values + 1)*sizeof(double));
	orbit->r_p			= (double *)malloc((2*num_chi_values + 1)*sizeof(double));
	orbit->rs_p			= (double *)malloc((2*num_chi_values + 1)*sizeof(double));
	orbit->phi_p		= (double *)malloc((2*num_chi_values + 1)*sizeof(double));
	orbit->t_p			= (double *)malloc((2*num_chi_values + 1)*sizeof(double));

	int i;
	for(i = 0; i <= 2.0*num_chi_values; i++){	
		double chi_p 	= i*M_PI/num_chi_values;
		double r_p		= r_of_chi(chi_p);
		orbit->chi_p[i] = chi_p;
		orbit->r_p[i]	= r_p;
		orbit->rs_p[i] 	= rs_of_r(r_p);
		orbit->phi_p[i] = 0;//phi_of_chi(chi_p);
		orbit->t_p[i]	= 0;//t_of_chi(chi_p);
	}

}


void 	print_orbital_params(struct orbital_params *orbit)
{
	printf("Orbital parameters:\n");
		printf("\tp = %lg\n", orbit->p);
		printf("\te = %lg\n", orbit->e);
		printf("\tE = %lg\n", orbit->E);
		printf("\tL = %lg\n", orbit->L);
		printf("\tr_max = %lg\n", orbit->r_max);
		printf("\tr_min = %lg\n", orbit->r_min);
		printf("\trs_max = %lg\n", orbit->rs_max);
		printf("\trs_min = %lg\n", orbit->rs_min);
		printf("\tT_r = %lg\n", orbit->T_r);
		printf("\tOmega_r = %lg\n", orbit->Omega_r);
		printf("\tOmega_phi = %lg\n", orbit->Omega_phi);
		printf("\tOmega_r/Omega_phi = %lg\n", orbit->Omega_r/orbit->Omega_phi);
}

double 	get_orbital_p()
{
	return p;
}

double 	get_orbital_e()
{
	return e;
}

/** Calculates \f$ f= \left(1-\frac{2M}{r}\right) \f$
*
* @param r the Schwarzschild radius
* @return f(r)
*/
double f(double r)
{
	return 1.0 - 2.0 * M/r;
}



double orbit_contra_ut(double r)
{
	return 1.0/f(r)*orbit_energy();
}

double orbit_covariant_ut(double r)
{
	return -orbit_energy();
}

/**
* Returns contravariant r-component of the 4-velocity for a given chi
*/
double orbit_contra_ur(double chi)
{
	double r = r_of_chi(chi);
	if(r == r_max() || r == r_min() ) return 0;

	double inside_sqrt = orbit_energy()*orbit_energy() - f(r)*(orbit_ang_mom()*orbit_ang_mom()/(r*r) + 1.0);

	if(inside_sqrt < 0){
		if ( fabs(inside_sqrt) < 1e-10) return 0;
		else printf("Something is very wrong with the calculation of contra_ur (sqrt of a negative number)\n");
	}

	double contra_ur  = sqrt( inside_sqrt);

	if(chi > M_PI)  contra_ur = -contra_ur;

	return contra_ur;
}

double orbit_covariant_ur(double chi)
{
	return 1.0/f(r_of_chi(chi)) * orbit_contra_ur(chi);
}

double dt_dtau(double chi)
{
	return orbit_energy()/ f( r_of_chi(chi) );
}

/** Calculates \f$d\chi/dt\f$ for eccentric orbit about Schwarzschild black hole
 *
 * @param chi the chi value at which to determine \f$d\chi/dt\f$
 * @return \f$d\chi/dt\f$
 */
double dt_dchi(double chi)
{
	 return p*p/( pow((1.0+e*cos(chi)), 2.0)*(p-2.0-2.0*e*cos(chi)) )*sqrt( (pow(p - 2.0, 2.0) - 4.0*e*e)/ (p - 6.0 - 2.0*e*cos(chi)) );
}


double dphi_dchi(double chi)
{
	return sqrt( p/(p-6.0-2.0*e*cos(chi)) );
}

// The integrand required to compute T_r
 
 
static double t_integrand(double chi, void *params)
{
	return dt_dchi(chi);
}

// The integrand required to compute the phi orbital frequency
 
 
static double phi_integrand(double chi, void *params)
{

	return dphi_dchi(chi);
}

double t_of_chi(double chi)
{
	if(e == 0) return 0;

	double result, error;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);

	gsl_function F;
	F.function = &t_integrand;
	F.params = NULL;

	gsl_integration_qag (&F, 0.0, chi, 0, 1e-13, 100000, GSL_INTEG_GAUSS61, w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;
}

double phi_of_chi(double chi)
{
	if(e == 0) return 0;

	double result, error;


	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);

	gsl_function F;
	F.function = &phi_integrand;
	F.params = NULL;

	gsl_integration_qag (&F, 0.0, chi, 0, 1e-13, 100000, GSL_INTEG_GAUSS61, w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;
}

// Computes the radial orbital period

double orbit_radial_period()
{
	if(e==0) return 0;
	double result;
	if(p == 6.0+2.0*e) result = INFINITY;
	else result = t_of_chi(2.0*M_PI);

	return result;
}

/** Computes the radial orbital frequency \f$ \Omega_r = 2\pi/T_r \f$
*
* @return the radial orbital frequency \f$ \Omega_r \f$
*/
double orbit_radial_freq()
{
	if(e==0) return 0;
	if(e == 0 && p <= 6) return 0;
	double result;
	if(p == 6.0+2.0*e) result = 0;
	else result = 2.0*M_PI/orbit_radial_period();

	return result;
}

/** Computes the azimuthal(phi) orbital frequency \f$  \Omega_\phi = \phi(\chi=2\pi)/T_r = 2\phi(\chi=\pi)/T_r \f$
*
* @return the azimuthal orbital frequency \f$ \Omega_\phi \f$
*/
double orbit_phi_freq()
{
	double result;

	if(e==0) return pow(p,-1.5);


	if(p == 6.0+2.0*e) result = pow( (1.0+e)/(6.0+2.0*e), 1.5);
	else result = phi_of_chi(2.0*M_PI)/orbit_radial_period();;

	return result;
}

/** Calculates the minimum radius attained by the orbit  \f$ =\left( \frac{p}{1+e}\right) \f$
*
* @return the minimum radius of the orbit
*/
double r_min()
{
	return p/(1.0+e);
}

/** Calculates the maximum radius attained by the orbit \f$ =\left( \frac{p}{1-e}\right) \f$
*
* @return the maximum radius of the orbit
*/
double r_max()
{
	return p/(1.0-e);
}

/** The value of r for a given \f$\chi\f$
*
* @param chi \f$\chi\f$
* @return the radius for the supplied value of \f$\chi\f$
*/
double r_of_chi(double chi)
{
	return p/(1.0+e*cos(chi));
}

/** Calculates the orbital energy \f$ = \left( \sqrt{\frac{(p-2-2e)(p-2+2e)}{p(p-3-e^2)}}\right) \f$
*
* @return the orbital energy
*/
double orbit_energy()
{
	return sqrt( (p-2.0-2.0*e)*(p-2.0+2.0*e)/( p*(p-3.0-e*e) ) );
}

/** Calculates the orbital angular momentum \f$ = \left( \frac{p}{\sqrt{p-3-e^2}} \right) \f$
*
* @return the orbital angular momentum
*/
double orbit_ang_mom()
{
	return p/sqrt( (p-3.0-e*e) );
}

/** Calculates the tortoise (\f$r_*\f$) coordinate for a given r
*
* @param  r the Schwarzschild radius
* @return the tortoise (\f$r_*\f$) coordinate
*/
double rs_of_r(double r){
	return r + 2.0*M*log(r/(2.0*M) - 1.0);
}

/*static double invert_formula(double r, void *params)
{
	double r_star = *(double *)params;
	return rs_of_r(r) - r_star;
}*/


// Converts the r_star coordinate to the Boyer-Lindquist radial coordinate.
 
// @param rs the r_star coordinate.

// Below is not a good solution as exp(rs) gets very large quickly
/*double r_of_rs(double rs)
{
	printf("%e %e %e\n", rs, rs/(2.0*M)-1.0, exp(rs/(2.0*M)-1.0));
	return 2.0*(M+M*gsl_sf_lambert_W0(exp(rs/(2.0*M)-1.0)));
}*/

double r_of_rs ( double r_st )
{
    
    double r = 5.0;  
    
	if(r_st < -50.0) return 2.0;
    
    double eps = 1e-15;
    double tol = 0.0000001;
     
// the following retarded if &else statements are needed because for values of R_star<0, we are close to r=2 singularity
// thus we must pick the initial guess values for r_0 very carefully. otherwise Newton-Raphson will give the first tangent
// line with root at r<2.
// for some reason, the tolerance also needed to be modified for r very close to r=2. even with this drawback (tol upto 0.1) 
// the residual values are shown to equal zero up to one part in 10^-14.

  if (r_st >= 0)
  {
  r = r_st+3;  // initialize r
  }
  else
 {
  if ( r_st >=-15 )
  { r = 2 + 1e-5; }
  else
  {
  if ( r_st >= -29 )
  {  r = 2 + 1e-7; }
  else
  {
     if (r_st >= -38 )
     { r = 2 + 1e-8;
       tol = 1e-6; 
        }
     else 
     { if (r_st >= -42 )
         { r = 2 + 1e-9;
           tol = 1e-4; 
         }
       else
       {  if (r_st >= -47 )
          { r = 2+1e-10;
            tol = 1e-4;
          }
          else 
          {  if (r_st >=-52 )
             { r = 2+1e-12;
               tol = 1e-3;
              }
       //    else
         //    {  if (r_st >=-54 )
          //      { r = 2+1e-12;
           //       tol = 1e-3;
           //     }
              //  else
             //   {  if (r_st >= -56 )
              //     { r = 2+1e-12;
              //      tol = 1e-3;
              //     }  
                   else
                   {   if (r_st >= -60 )
                       { r = 2 + 1e-14;
                         tol = 1e-2;
                       }
                       else 
                       {   if (r_st >= -65 )
                           { r = 2+1e-14;
                             tol = 0.1;
                           }
                       }
                    }
                  
              //  }
          //   }
          
            }
        }
     }
   }
  }
 }


//********************************************


  double g = r+2*log(r/2-1)-r_st ;

// the while loop evaluates the RHS until the residue is less than tolerance (but greater than eps). This is why eps < tol
// so that we can make the while loop break when we meet our criterion for the root.

   while (fabs(g) > eps )
    {
         
          double g =r+2.0*log(r/2.0-1.0)-r_st ;
          double dg_dr = r/(r-2.0);
          double y = g;
          double m = dg_dr;
          r = -y/m + r;
         
        if ( fabs(g) < tol)
               break;
      
  
    }
  
  //printf("R_star=%d Residue=%.12e r=%.16e number of steps=%d\n",R_star[j], y, r, count);
 
    return r;
}

/*static double invert_formula(double r, void *params)
{
	double *r_star = (double *)params;
	return rs_of_r(r) - *r_star;
}

double r_of_rs(double rs)
{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *solver_type;
	gsl_root_fsolver *solver;
	double root = 0;
	double x_lo = 2.0+1E-14, x_hi = rs + 1000;
	gsl_function F;

	F.function = &invert_formula;
	F.params = &rs;

	solver_type = gsl_root_fsolver_brent;
	solver = gsl_root_fsolver_alloc (solver_type);
	gsl_root_fsolver_set (solver, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (solver);
		root = gsl_root_fsolver_root (solver);
		x_lo = gsl_root_fsolver_x_lower (solver);
		x_hi = gsl_root_fsolver_x_upper (solver);
		status = gsl_root_test_interval (x_lo, x_hi, 0, DESIRED_INVERT_ACCURACY);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (solver);

	return root;
}*/

double locate_effective_horizon()
{
	return ceil(rs_of_r(2.0 + 1E-14));
}



////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Circular orbit stuff ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

double circular_orbit_phi_freq()
{
	return	sqrt(M)/pow(p, 1.5);
}

double circular_orbit_ang_mom()
{
	return pow(M*p,0.5)*pow(1.0 - 3.0*M/p, -0.5);
}

double circular_orbit_energy()
{
	return (1.0 - 2.0*M/p)*pow(1.0 - 3.0*M/p,-0.5);
}






/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////	R8 static mode			/////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


double calc_alk(int l, int k)
{

	double f1 = factorial(l + k - 1.0);
	double f2 = factorial(l - k + 1);
	double f3 = factorial(k	+ 1);
	double f4 = factorial(k);

	return	l*(l+1.0)*f1/( f2*f3*f4 );
}

double calc_blk(int l, int k)
{
	double sum = 0;
	int t;
	for(t = 0; t <= l-k+1; t++){
		sum += pow(-1.0, t) * calc_alk(l, k+t)/(t+1.0);
	}
	return sum;
}


double R8_static_out_solution(int l, double r_p)
{
	double phi_out = 0;


	if(l == 1) phi_out = 1.0/r_p;
	if(l >= 2){
		/*int k;
		double sum = 0;
		double x = r_p/2.0 - 1.0;
		for(k = 0; k <= l+1; k++){
			double blk = calc_blk( l,  k);
			sum += 1.0/(1.0+x) * blk * pow(x,k);

		}*/
		phi_out =  R8_out_odd_static_n_zero(r_p, l);//static_in_solution(l, r_p) * log(f_p) + sum;
	}

	return phi_out;
}

double R8_static_out_rs_deriv_solution(int l, double r_p)
{
	double phi_out_rs_deriv = 0;

	double f_p = f(r_p);

	if(l == 1) phi_out_rs_deriv = -1.0/(r_p*r_p)*f_p;
	if(l >= 2){
		/*int k;
		double sum = 0;
		double x = r_p/2.0 - 1.0;
		for(k = 0; k <= l+1; k++){
			double blk = calc_blk( l,  k);
			sum += 1.0/pow(x+1.0,2.0) * ( (1.0 + x)* blk*k*pow(x,k-1.0) - blk*pow(x,k) );
		}*/

		phi_out_rs_deriv = f_p*R8_out_odd_static_n_zero_r_deriv(r_p,l);//static_in_rs_deriv_solution(l, r_p)*log(f_p) + static_in_solution(l,r_p)*2.0/(r_p*r_p) + f_p/2.0 * sum;
	}



	return phi_out_rs_deriv;
}

double R8_static_in_solution(int l, double r_p)
{
	double phi_in = 0;

	if(l == 1) phi_in = r_p*r_p;
	if(l >= 2){
		double x 	= r_p/2.0 - 1.0;
		int k;
		double sum = 0;
		for(k = 0; k <= l+1; k++){
			double alk 	= calc_alk(l,k);
			sum 		+= x / ( 1.0 + x ) * alk*pow(x,k);
		}
		phi_in = sum;
	}

	return phi_in;
}

double R8_static_in_rs_deriv_solution(int l, double r_p)
{
	double phi_in_rs_deriv = 0;

	double f_p = f(r_p);

	if(l == 1) phi_in_rs_deriv = 2.0*r_p*f_p;
	if(l >= 2){
		double x 	= r_p/2.0 - 1.0;
		int k;
		double sum = 0;
		for(k = 0; k <= l+1; k++){
			double alk 	= calc_alk(l,k);
			sum 		+= 1.0/(1.0 + x) * alk*( k + 1.0/(1.0 + x))*pow(x,k);
		}
		phi_in_rs_deriv = f_p/2.0*(sum);
	}

	return phi_in_rs_deriv;
}

double R8_out_odd_static_n_zero(double r_chi, int l)
{
	double r_0 = r_chi;
	double sum8, Dsum8;
	double M = 1.0;
	int k; 
	int k_max = 1;  
	double c[2];
	double del;

	// we have a_l=1, smaller terms are zero
	c[0] 	= 1.0; 

	sum8 	= c[0]/pow( r_0, l ); 
	Dsum8 	= -l*c[0]/pow( r_0, l+1.0 );
	//printf("0, %e\t", Dsum8);
	for (k=1; k<=k_max; k++)
	{
		c[1] 	= 2*M*c[0]*( ( (k-1) + l)*( (k-1) + l + 1.) - 2. ) /( (k )*(k + 2.*l + 1.) );

		sum8 	+= c[1]/pow( r_0, l+k );
		Dsum8 	+= -(l+k)*c[1]/pow( r_0, l+k+1. );

		del	=  c[1]/pow( r_0, l+k );

		//printf("%d, %e\t", k, Dsum8);

		if( fabs(del/sum8 ) > 1e-16 )	
		{	k_max = k_max + 1;	}

		c[0] = c[1];
		c[1] = 0.0; // freeing c[1]
		
	}

	return sum8;

}

double R8_out_odd_static_n_zero_r_deriv(double r_chi, int l)
{
	double r_0 = r_chi;
	double sum8, Dsum8;
	double M = 1.0;
	int k; 
	int k_max = 1;  
	double c[2];
	double del;

	// we have a_l=1, smaller terms are zero
	c[0] 	= 1.0; 

	sum8 	= c[0]/pow( r_0, l ); 
	Dsum8 	= -l*c[0]/pow( r_0, l+1.0 );
	//printf("0, %e\t", Dsum8);
	for (k=1; k<=k_max; k++)
	{
		c[1] 	= 2*M*c[0]*( ( (k-1) + l)*( (k-1) + l + 1.) - 2. ) /( (k )*(k + 2.*l + 1.) );

		sum8 	+= c[1]/pow( r_0, l+k );
		Dsum8 	+= -(l+k)*c[1]/pow( r_0, l+k+1. );

		del	=  c[1]/pow( r_0, l+k );

		//printf("%d, %e\t", k, Dsum8);

		if( fabs(del/sum8 ) > 1e-16 )	
		{	k_max = k_max + 1;	}

		c[0] = c[1];
		c[1] = 0.0; // freeing c[1]
		
	}

	return Dsum8;

}
