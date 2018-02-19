/*
	This file uses {p,e} as input to compute the 4 unknown coefficients a1,a2,b1,b2 of Golbourn pg. 139
	to construct the inhomogeneous solutions for h^1, h^3 and h^6 for the static monopole (m=n=0) piece.

	It then needs r_p as input to output the value of the fields (and their r-derivatives) at r=r_p .

	The degree of freedom is 4 because we have 3 field and 1 gauge equation, this (3-1)*2= 4 d.o.f.

	The first half of the code performs the 4 integrals necessary to compute the coefficients. This is done by inverting
	a 4 by 4 Real matrix of homogeneous solutions construted from h^1, h^3 and their r-derivatives. One could instead
	have used {h^1,h^6} or {h^3,h^6} pairings, but I stuck with {h^1,h^3} which also turned out to be what Golbourn did.

	Since the matrix elements are analytic functions of r, I did the matrix inversion using Maple since C only inverts numerically.
	As this matrix, Phi_inv, is multiplied by a column vector J= (0,0,J^1,J^3)^T, we only need to know the right hand half
	(last 2 columns) of the matrix. After Phi_inv*J, we do 4 integrals, which give us a1,a2,b1,b2.

	Second half deals with constructing the imhomogeneous solutions. Again refer to Golbourn pg. 139 for details but
	basically we simply multiply the coefficients by homog. sol. and add them. There is a fine detail about adding an inner	
	homogeneous piece to the outer homog. solution to make sure it has the right 1/r drop-off as r-> \infty. This is 
	explained in Golbourn pg. 138 and Barack's manuscript as well.


	CODE details:
	- Currently I use 2 of my existing libraries for things like E_p, L_p, dt_dchi etc., you will need to change these to your own libraries.
	- Since I used quite a bit of CForm stuff from mathematica, I added some new definitions.
	- As it is the code simply prints out h^1,h^3,h^6 and their r-derivatives at r=r_p. You will need to modify it however you like to make it
	  a function that spits out these fields. 
	- Because we still want to know the value of the fields at a given point, I artificially inserted r_p (also used R) in the 2nd half.
	  That obviously needs to be turned into an argument that goes into your function that computes these fields.

*/

#include "GSF_ecc.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

double mu = 1.0;   // small particle mass.

#define Power(x, y)	(pow((double)(x), (double)(y)))
#define Sqrt(x)		(sqrt((double)(x)))

#define Abs(x)		(fabs((double)(x)))

#define Exp(x)		(exp((double)(x)))
#define Log(x)		(log((double)(x)))
#define ln(x)		(log((double)(x)))
#define Cos(x)		(cos((double)(x)))

#define Pi		M_PI

// Integrand for a1_hat
double integ1( double chi, void *params)
{
	struct orbital_params *orbit = (struct orbital_params *)params;
	double E_p = orbit->E;
	double L_p = orbit->L;
	
	double r = r_of_chi(chi);

	double Phi_inv_13 = (mu*(2.*M - 1.*r)*r*(6.4e-16*Power(M,7) - 
	   6.4e-16*Power(M,6)*r - 3.e-16*Power(M,4)*Power(r,3) + 
	   1.e-16*Power(M,3)*Power(r,4) - 1152.*M*Power(r,6) + 
	   576.*Power(r,7)))/
   (M*Sqrt(Pi)*(5.12e-14*Power(M,8) + 2.56e-14*Power(M,7)*r - 
	   5.12e-15*Power(M,6)*Power(r,2) + 
	   4.e-15*Power(M,5)*Power(r,3) - 
	   4.e-15*Power(M,4)*Power(r,4) - 
	   4.e-15*Power(M,3)*Power(r,5) - 
	   18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
	   4608.*Power(r,8)));

	double Phi_inv_14 = -((mu*(3.84e-15*Power(M,9) + 4.48e-15*Power(M,8)*r - 
		 1.92e-15*Power(M,7)*Power(r,2) - 
		 1.e-15*Power(M,6)*Power(r,3) + 
		 1.e-15*Power(M,5)*Power(r,4) - 
		 1.e-15*Power(M,4)*Power(r,5) - 
		 4608.*Power(M,3)*Power(r,6) + 
		 6912.*Power(M,2)*Power(r,7) - 3456.*M*Power(r,8) + 
		 576.*Power(r,9)))/
	 (M*Sqrt(Pi)*(5.12e-14*Power(M,8) + 2.56e-14*Power(M,7)*r - 
		 5.12e-15*Power(M,6)*Power(r,2) + 
		 4.e-15*Power(M,5)*Power(r,3) - 
		 4.e-15*Power(M,4)*Power(r,4) - 
		 4.e-15*Power(M,3)*Power(r,5) - 
		 18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
		 4608.*Power(r,8))));

	double J1 = -16.0*sqrt(Pi)/(E_p*pow(r,3.0))*(2.0*E_p*E_p*r*r-(1-2/r)*r*r-(1-2/r)*L_p*L_p);
	double J3 = -16.0*sqrt(Pi)/(E_p*pow(r,3.0))*(r*r+L_p*L_p);

	return (Phi_inv_13*J1 + Phi_inv_14*J3)*(1.0/(2.0*M_PI));
}

// Integrand for a2_hat
double integ2( double chi, void *params)
{
	struct orbital_params *orbit = (struct orbital_params *)params;
	
	double E_p = orbit->E;
	double L_p = orbit->L;
	
	double r = r_of_chi(chi);
	double Phi_inv_23 = (mu*Power(r,7)*(-2.*M + 1.*r)*
	 (-576.*M + 288.*r + 384.*M*Log((-2*M + r)/r) - 
	   96.*r*Log((-2*M + r)/r)))/
   (M*Sqrt(Pi)*(5.12e-14*Power(M,8) + 2.56e-14*Power(M,7)*r - 
	   5.12e-15*Power(M,6)*Power(r,2) + 
	   4.e-15*Power(M,5)*Power(r,3) - 
	   4.e-15*Power(M,4)*Power(r,4) - 
	   4.e-15*Power(M,3)*Power(r,5) - 
	   18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
	   4608.*Power(r,8)));

	double Phi_inv_24 =(-4*mu*Power(r,4)*(1.e-16*Power(M,5) + 
	   384.*Power(M,3)*Power(r,2) - 
	   672.*Power(M,2)*Power(r,3) + 384.*M*Power(r,4) - 
	   72.*Power(r,5) - 192.*Power(M,3)*Power(r,2)*
		Log((-2*M + r)/r) + 
	   288.*Power(M,2)*Power(r,3)*Log((-2*M + r)/r) - 
	   144.*M*Power(r,4)*Log((-2*M + r)/r) + 
	   24.*Power(r,5)*Log((-2*M + r)/r)))/
   (M*Sqrt(Pi)*(5.12e-14*Power(M,8) + 2.56e-14*Power(M,7)*r - 
	   5.12e-15*Power(M,6)*Power(r,2) + 
	   4.e-15*Power(M,5)*Power(r,3) - 
	   4.e-15*Power(M,4)*Power(r,4) - 
	   4.e-15*Power(M,3)*Power(r,5) - 
	   18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
	   4608.*Power(r,8)));

	double J1 = -16.0*sqrt(Pi)/(E_p*pow(r,3.0))*(2.0*E_p*E_p*r*r-(1-2/r)*r*r-(1-2/r)*L_p*L_p);
	double J3= -16.0*sqrt(Pi)/(E_p*pow(r,3.0))*(r*r+L_p*L_p);

	return (Phi_inv_23*J1 + Phi_inv_24*J3)*(1.0/(2.0*M_PI));
}

// Integrand for b1
double integ3( double chi, void *params)
{
	struct orbital_params *orbit = (struct orbital_params *)params;	
	
	double E_p = orbit->E;
	double L_p = orbit->L;
	
	double r = r_of_chi(chi);
	double Phi_inv_33 = (mu*(2.*M - 1.*r)*Power(r,4)*
	 (4.e-16*Power(M,6) - 8.e-16*Power(M,5)*r + 
	   3.2e-16*Power(M,4)*Power(r,2) + 
	   2048.*Power(M,3)*Power(r,3) - 
	   1280.*Power(M,2)*Power(r,4) + 192.*Power(r,6) - 
	   6144.*Power(M,3)*Power(r,3)*Log(r/M) + 
	   1536.*Power(M,2)*Power(r,4)*Log(r/M)))/
   (Power(M,3)*Sqrt(Pi)*(5.12e-14*Power(M,8) + 
	   2.56e-14*Power(M,7)*r - 5.12e-15*Power(M,6)*Power(r,2) + 
	   4.e-15*Power(M,5)*Power(r,3) - 
	   4.e-15*Power(M,4)*Power(r,4) - 
	   4.e-15*Power(M,3)*Power(r,5) - 
	   18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
	   4608.*Power(r,8)));

	double Phi_inv_34 = (-4*mu*(2.*M - 1.*r)*Power(r,4)*
	 (1.e-15*Power(M,6) - 1.e-15*Power(M,5)*r - 
	   256.*Power(M,4)*Power(r,2) + 
	   256.0000000000000001*Power(M,3)*Power(r,3) + 
	   128.*Power(M,2)*Power(r,4) - 192.*M*Power(r,5) + 
	   48.*Power(r,6) - 1536.*Power(M,4)*Power(r,2)*Log(r/M) + 
	   1536.*Power(M,3)*Power(r,3)*Log(r/M) - 
	   384.*Power(M,2)*Power(r,4)*Log(r/M)))/
   (Power(M,3)*Sqrt(Pi)*(5.12e-14*Power(M,8) + 
	   2.56e-14*Power(M,7)*r - 5.12e-15*Power(M,6)*Power(r,2) + 
	   4.e-15*Power(M,5)*Power(r,3) - 
	   4.e-15*Power(M,4)*Power(r,4) - 
	   4.e-15*Power(M,3)*Power(r,5) - 
	   18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
	   4608.*Power(r,8)));

	double J1 = -16.0*sqrt(Pi)/(E_p*pow(r,3.0))*(2.0*E_p*E_p*r*r-(1-2/r)*r*r-(1-2/r)*L_p*L_p);
	double J3= -16.0*sqrt(Pi)/(E_p*pow(r,3.0))*(r*r+L_p*L_p);

	return (Phi_inv_33*J1 + Phi_inv_34*J3)*(1.0/(2.0*M_PI));
}

// Integrand for b2
double integ4( double chi, void *params)
{
	struct orbital_params *orbit = (struct orbital_params *)params;	
	
	double E_p = orbit->E;
	double L_p = orbit->L;
	
	double r = r_of_chi(chi);
	double Phi_inv_43 = (192*mu*(2.*M - 1.*r)*(4.*M - 1.*r)*Power(r,7))/
 (M*Sqrt(Pi)*(5.12e-14*Power(M,8) + 2.56e-14*Power(M,7)*r - 
	   5.12e-15*Power(M,6)*Power(r,2) + 
	   4.e-15*Power(M,5)*Power(r,3) - 
	   4.e-15*Power(M,4)*Power(r,4) - 
	   4.e-15*Power(M,3)*Power(r,5) - 
	   18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
	   4608.*Power(r,8)));

	double Phi_inv_44 = (-192*mu*(2.*M - 1.*r)*Power(r,6)*
	 (4.*Power(M,2) - 4.*M*r + Power(r,2)))/
   (M*Sqrt(Pi)*(5.12e-14*Power(M,8) + 2.56e-14*Power(M,7)*r - 
	   5.12e-15*Power(M,6)*Power(r,2) + 
	   4.e-15*Power(M,5)*Power(r,3) - 
	   4.e-15*Power(M,4)*Power(r,4) - 
	   4.e-15*Power(M,3)*Power(r,5) - 
	   18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
	   4608.*Power(r,8)));

	double J1 = -16.0*sqrt(Pi)/(E_p*pow(r,3.0))*(2.0*E_p*E_p*r*r-(1-2/r)*r*r-(1-2/r)*L_p*L_p);
	double J3= -16.0*sqrt(Pi)/(E_p*pow(r,3.0))*(r*r+L_p*L_p);


	return (Phi_inv_43*J1 + Phi_inv_44*J3)*(1.0/(2.0*M_PI));
}


void construct_static_monopole (struct lm_mode_data *lm_mode, struct orbital_params *orbit )
{
		printf("Calculating l=0 m=0 n=0 Analytically\n");



	double a1_hat, a2_hat, b1, b2, a1, a2;  // the unknown coefficients we are solving for. a1_hat = a1 - 2E_p, a2_hat = a2 + E_p.
	double R;   // Declaring the variables. r and chi are used for the integral. r_p and R are for the location of particle.
	double P, Q;//, W, K, L;  // These are all functions of r used to construct the homo. sol. See Barack's notes.



	double E_p, L_p;//, T_r; // Energy, ang. momentum and radial period.
	//double J1, J3;  // Sources for fields h^1 and h^3 for l=0,m=0.


	// Computing E_p, L_p, T_r. 
	E_p = orbit->E;
	L_p = orbit->L;
	//T_r = orbit->T_r;

	//printf(" %.14f %.14f %.14f\n", T_r, E_p, L_p);


//*********************  Inverse Matrix *****************************************
	// Since C can't invert analytically, we use the analytical expressions obtained by Maple (or Mathematica) for the
	// matrix inversion. And since the source vector's first 2 entries are zero, we only need the last 2 columns of this
	// inverted matrix, which we call Phi_inv. So we need matrix elements: 13,14,23,24,33,34,43,44.

	// Declaring the 8 matrix elements.
	//double Phi_inv_13, Phi_inv_14, Phi_inv_23, Phi_inv_24, Phi_inv_33, Phi_inv_34, Phi_inv_43, Phi_inv_44; 

	// Writing in the 4 integrals used to compute a1_hat, a2_hat, b1, b2.
	// I used Maple to analytically invert the matrix Phi so that I could retain the functional dependence on r.
	// For example: the integral for a1 is given by the first row of the column vector C = Phi_inv*J so this is given by:
	// 			Int[ Phi_inv_13*J1 + Phi_inv_14*J3) dr
	// For sake of clarity I pull the factor of dt_dchi out of J1, J3 and multiply both terms by it inside the integral.






//************* doing the 4 integrals all back-to-back.*************


	if(orbit->e==0){
		a1_hat 	= Pi*integ1(M_PI/2.0, orbit);
		a2_hat	= Pi*integ2(M_PI/2.0, orbit);
		b1		= Pi*integ3(M_PI/2.0, orbit);
		b2		= Pi*integ4(M_PI/2.0, orbit);
	}else{
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
		gsl_integration_workspace * x = gsl_integration_workspace_alloc (100000);
		gsl_integration_workspace * y = gsl_integration_workspace_alloc (100000);
		gsl_integration_workspace * z = gsl_integration_workspace_alloc (100000);
	   
	  	double error;
	 
	  	gsl_function F1, F2, F3, F4;
		F1.function = &integ1;
		F2.function = &integ2;
	 	F3.function = &integ3;
		F4.function = &integ4;

		gsl_integration_qag (&F1, 0.0, Pi, 0, 1e-13, 100000, GSL_INTEG_GAUSS61, w, &a1_hat, &error);
		gsl_integration_qag (&F2, 0.0, Pi, 0, 1e-13, 100000, GSL_INTEG_GAUSS61, x, &a2_hat, &error); 
		gsl_integration_qag (&F3, 0.0, Pi, 0, 1e-13, 100000, GSL_INTEG_GAUSS61, y, &b1, &error); 
		gsl_integration_qag (&F4, 0.0, Pi, 0, 1e-13, 100000, GSL_INTEG_GAUSS61, z, &b2, &error); 
		  // printf ("result  = % .10f %.10f %.10f %.10f\n", a1_hat, a2_hat, b1,b2);
		gsl_integration_workspace_free (w); 
		gsl_integration_workspace_free (x); 
		gsl_integration_workspace_free (y); 
		gsl_integration_workspace_free (z); 
	}
// ****** Done with all 4 integrals ************

	// Integral gives us a1_hat, a2_hat. We want a1, a2. See Golbourn pg. 139 for details.
	a1 = a1_hat + 2.0*E_p;
	a2 = a2_hat - E_p;

	//printf("a1=%.12f \t a2=%.12f \t b1=%.10f \t b2=%.10f\n", a1,a2,b1, b2);


	int i;
	for(i = 0; i < 2*NUM_CHI_VALUES +1; i++){

		double r_p = orbit->r_p[i];
		double f_p = f(r_p);


	//**************** Construct Inhomogeneous Solutions ******************
		//r_p =12.5;
		double En = E_p; // usual BS because mathematica does its own shit with E_p
		R =r_p;   // Again because mathematica does not like r_p

		P=r_p*r_p+2*r_p*M+4*M*M; 
		Q=r_p*r_p*r_p-r_p*r_p*M-2*r_p*M*M+12*M*M*M; 
		//W=3*r_p*r_p*r_p-r_p*r_p*M-4*r_p*M*M-28*M*M*M/3; 
		//K=r_p*r_p*r_p*M-5*r_p*r_p*M*M-20*r_p*M*M*M/3+28*pow(M,4); 
		//L=pow(r_p,4)-3*pow(r_p,3)*M+16*r_p*M*M*M-24*pow(M,4);

		double h1[3] = {-f_p, 1/f_p, 1};
		double h2[3] = {-2*f_p*M/pow(r_p,3)*P, 2/f_p/pow(r_p,3)*Q, 2*f_p*P/pow(r_p,2)};
		double h3[3] = {-pow(M/r_p,4), M*M*M*(3*M - 2*r_p)/(pow(r_p*r_p*f_p,2)), pow(M/r_p,3)};
		double h4[3] = {-(M*(28*Power(M,3) + 12*Power(M,2)*r_p + 3*M*Power(r_p,2) - 
		     9*Power(r_p,3) + 3*(8*Power(M,3) - Power(r_p,3))*Log(f_p) + 
		     24*Power(M,3)*Log(r_p/M)))/(3.*Power(r_p,4)),
	   -(-84*Power(M,4) + 20*Power(M,3)*r_p + 
		   15*Power(M,2)*Power(r_p,2) - 3*M*Power(r_p,3) + 
		   (-72*Power(M,4) + 48*Power(M,3)*r_p - 9*M*Power(r_p,3) + 
		      3*Power(r_p,4))*Log(f_p) - 
		   (72*Power(M,4) - 48*Power(M,3)*r_p)*Log(r_p/M))/
		(3.*Power(f_p,2)*Power(r_p,4)),
	   -(-28*Power(M,3) - 12*Power(M,2)*r_p - 3*M*Power(r_p,2) + 
		   (-24*Power(M,3) + 3*Power(r_p,3))*Log(f_p) - 
		   24*Power(M,3)*Log(r_p/M))/(3.*Power(r_p,3))};

		double h5[3];

		int j;
		for (j=0; j<=2; j++)
		{
			h5[j] = h4[j] - 3.0*h1[j] + 1.5*h2[j];
		}
	
		double H1_in=2*sqrt(Pi)*r_p/mu*( a1*(h1[0]+f_p*f_p*h1[1]) + a2*(h2[0]+f_p*f_p*h2[1]) ); 
		double H1_out=2*sqrt(Pi)*r_p/mu*( b1*(h3[0]+f_p*f_p*h3[1]) + b2*(h5[0]+f_p*f_p*h5[1]) + E_p*( 2.0*(h1[0]+f_p*f_p*h1[1]) - (h2[0]+f_p*f_p*h2[1]) ) ); 


		double H3_in=4*sqrt(Pi)*r_p/mu*( a1*h1[2] + a2*h2[2] );
		double H3_out=4*sqrt(Pi)*r_p/mu*( b1*(h3[2]) + b2*(h5[2]) + E_p*( 2.0*h1[2] - h2[2] ) );
		double H6_in=2*sqrt(Pi)*r_p/(mu*f_p)*( a1*(h1[0]-f_p*f_p*h1[1]) + a2*(h2[0]-f_p*f_p*h2[1]) ); 
		double H6_out=2*sqrt(Pi)*r_p/(mu*f_p)*( b1*(h3[0]-f_p*f_p*h3[1]) + b2*(h5[0]-f_p*f_p*h5[1]) + E_p*( 2.0*(h1[0]-f_p*f_p*h1[1]) - (h2[0]-f_p*f_p*h2[1]) ) ); 

		double dH1_in = (2*a2*Sqrt(Pi)*R*((-2*M*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,3) - (4*Power(M,2)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,5) + 
		    (6*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) + 
		    (2*(1 - (2*M)/R)*(-2*Power(M,2) - 2*M*R + 3*Power(R,2)))/Power(R,3) + 
		    (4*M*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,5) - 
		    (6*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,4)))/mu + 
	   (2*a2*Sqrt(Pi)*((-2*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) + 
		    (2*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,3)))/mu;

		double dH3_in = (4*Sqrt(Pi)*R*((2*a2*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,2) + (4*a2*M*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) - 
		    (4*a2*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3)))/mu + 
	   (4*Sqrt(Pi)*(a1 + (2*a2*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,2)))/mu;
	

		double dH1_out = (2*Sqrt(Pi)*(b1*(-(Power(M,4)/Power(R,4)) + (Power(M,3)*(3*M - 2*R))/Power(R,4)) + 
		    En*((2*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) - 
		       (2*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,3)) + 
		    b2*(-3*(-1 + (2*M)/R) - (3.*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) - 
		       (M*(28*Power(M,3) + 12*Power(M,2)*R + 3*M*Power(R,2) - 9*Power(R,3) + 3*(8*Power(M,3) - Power(R,3))*Log(1 - (2*M)/R) + 
		            24*Power(M,3)*Log(R/M)))/(3.*Power(R,4)) + 
		       Power(1 - (2*M)/R,2)*(-3/(1 - (2*M)/R) + (3.*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/
		           ((1 - (2*M)/R)*Power(R,3)) - (-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		             (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		             (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M))/(3.*Power(1 - (2*M)/R,2)*Power(R,4))))))/mu + 
	   (2*Sqrt(Pi)*R*(b1*((4*Power(M,4))/Power(R,5) - (4*Power(M,3)*(3*M - 2*R))/Power(R,5) - (2*Power(M,3))/Power(R,4)) + 
		    En*((2*M*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,3) + (4*Power(M,2)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,5) - 
		       (6*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) - 
		       (2*(1 - (2*M)/R)*(-2*Power(M,2) - 2*M*R + 3*Power(R,2)))/Power(R,3) - 
		       (4*M*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,5) + 
		       (6*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,4)) + 
		    b2*((6*M)/Power(R,2) - (3.*M*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,3) - (6.*Power(M,2)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,5) + 
		       (9.*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) - 
		       (M*(12*Power(M,2) + (24*Power(M,3))/R + 6*M*R - 27*Power(R,2) + (6*M*(8*Power(M,3) - Power(R,3)))/((1 - (2*M)/R)*Power(R,2)) - 
		            9*Power(R,2)*Log(1 - (2*M)/R)))/(3.*Power(R,4)) + 
		       (4*M*(28*Power(M,3) + 12*Power(M,2)*R + 3*M*Power(R,2) - 9*Power(R,3) + 3*(8*Power(M,3) - Power(R,3))*Log(1 - (2*M)/R) + 
		            24*Power(M,3)*Log(R/M)))/(3.*Power(R,5)) + 
		       Power(1 - (2*M)/R,2)*((6*M)/(Power(1 - (2*M)/R,2)*Power(R,2)) + 
		          (3.*(-2*Power(M,2) - 2*M*R + 3*Power(R,2)))/((1 - (2*M)/R)*Power(R,3)) - 
		          (6.*M*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/(Power(1 - (2*M)/R,2)*Power(R,5)) - 
		          (9.*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/((1 - (2*M)/R)*Power(R,4)) - 
		          (20*Power(M,3) + 30*Power(M,2)*R - 9*M*Power(R,2) - (72*Power(M,4) - 48*Power(M,3)*R)/R + 
		             (2*M*(-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4)))/((1 - (2*M)/R)*Power(R,2)) + 
		             (48*Power(M,3) - 27*M*Power(R,2) + 12*Power(R,3))*Log(1 - (2*M)/R) + 48*Power(M,3)*Log(R/M))/
		           (3.*Power(1 - (2*M)/R,2)*Power(R,4)) + (4*M*
		             (-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		               (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		               (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M)))/(3.*Power(1 - (2*M)/R,3)*Power(R,6)) + 
		          (4*(-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		               (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		               (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M)))/(3.*Power(1 - (2*M)/R,2)*Power(R,5))) + 
		       (4*M*(1 - (2*M)/R)*(-3/(1 - (2*M)/R) + (3.*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/
		             ((1 - (2*M)/R)*Power(R,3)) - (-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		               (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		               (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M))/(3.*Power(1 - (2*M)/R,2)*Power(R,4))))/Power(R,2))))/mu;



		double dH3_out = (4*Sqrt(Pi)*R*((-3*b1*Power(M,3))/Power(R,4) + En*((-2*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,2) - 
		       (4*M*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) + (4*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3)) + 
		    b2*((3.*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,2) + (6.*M*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) - 
		       (6.*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) - 
		       (-12*Power(M,2) - (24*Power(M,3))/R - 6*M*R + (2*M*(-24*Power(M,3) + 3*Power(R,3)))/((1 - (2*M)/R)*Power(R,2)) + 
		          9*Power(R,2)*Log(1 - (2*M)/R))/(3.*Power(R,3)) + 
		       (-28*Power(M,3) - 12*Power(M,2)*R - 3*M*Power(R,2) + (-24*Power(M,3) + 3*Power(R,3))*Log(1 - (2*M)/R) - 24*Power(M,3)*Log(R/M))/
		        Power(R,4))))/mu + (4*Sqrt(Pi)*((b1*Power(M,3))/Power(R,3) + 
		    En*(2. - (2*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,2)) + 
		    b2*(-3 + (3.*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,2) - 
		       (-28*Power(M,3) - 12*Power(M,2)*R - 3*M*Power(R,2) + (-24*Power(M,3) + 3*Power(R,3))*Log(1 - (2*M)/R) - 24*Power(M,3)*Log(R/M))/
		        (3.*Power(R,3)))))/mu;

		double dH6_in = (2*Sqrt(Pi)*R*((-4*a1*M)/Power(R,2) + a2*((-2*M*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,3) - 
		       (4*Power(M,2)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,5) + (6*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) - 
		       (2*(1 - (2*M)/R)*(-2*Power(M,2) - 2*M*R + 3*Power(R,2)))/Power(R,3) - 
		       (4*M*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,5) + 
		       (6*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,4))))/(mu*(1 - (2*M)/R)) + 
	   (2*Sqrt(Pi)*(a1*(-2 + (4*M)/R) + a2*((-2*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) - 
		       (2*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,3))))/(mu*(1 - (2*M)/R)) - 
	   (4*M*Sqrt(Pi)*(a1*(-2 + (4*M)/R) + a2*((-2*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) - 
		       (2*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,3))))/(mu*Power(1 - (2*M)/R,2)*R);

		double dH6_out = (2*Sqrt(Pi)*(b1*(-(Power(M,4)/Power(R,4)) - (Power(M,3)*(3*M - 2*R))/Power(R,4)) + 
		    En*(2.*(-2 + (4*M)/R) + (2*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) + 
		       (2*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,3)) + 
		    b2*(-3*(-1 + (2*M)/R) - (3.*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) - 
		       (M*(28*Power(M,3) + 12*Power(M,2)*R + 3*M*Power(R,2) - 9*Power(R,3) + 3*(8*Power(M,3) - Power(R,3))*Log(1 - (2*M)/R) + 
		            24*Power(M,3)*Log(R/M)))/(3.*Power(R,4)) - 
		       Power(1 - (2*M)/R,2)*(-3/(1 - (2*M)/R) + (3.*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/
		           ((1 - (2*M)/R)*Power(R,3)) - (-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		             (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		             (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M))/(3.*Power(1 - (2*M)/R,2)*Power(R,4))))))/(mu*(1 - (2*M)/R)) - 
	   (4*M*Sqrt(Pi)*(b1*(-(Power(M,4)/Power(R,4)) - (Power(M,3)*(3*M - 2*R))/Power(R,4)) + 
		    En*(2.*(-2 + (4*M)/R) + (2*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) + 
		       (2*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,3)) + 
		    b2*(-3*(-1 + (2*M)/R) - (3.*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,3) - 
		       (M*(28*Power(M,3) + 12*Power(M,2)*R + 3*M*Power(R,2) - 9*Power(R,3) + 3*(8*Power(M,3) - Power(R,3))*Log(1 - (2*M)/R) + 
		            24*Power(M,3)*Log(R/M)))/(3.*Power(R,4)) - 
		       Power(1 - (2*M)/R,2)*(-3/(1 - (2*M)/R) + (3.*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/
		           ((1 - (2*M)/R)*Power(R,3)) - (-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		             (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		             (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M))/(3.*Power(1 - (2*M)/R,2)*Power(R,4))))))/(mu*Power(1 - (2*M)/R,2)*R) + 
	   (2*Sqrt(Pi)*R*(b1*((4*Power(M,4))/Power(R,5) + (4*Power(M,3)*(3*M - 2*R))/Power(R,5) + (2*Power(M,3))/Power(R,4)) + 
		    En*((-8.*M)/Power(R,2) + (2*M*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,3) + (4*Power(M,2)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,5) - 
		       (6*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) + 
		       (2*(1 - (2*M)/R)*(-2*Power(M,2) - 2*M*R + 3*Power(R,2)))/Power(R,3) + 
		       (4*M*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,5) - 
		       (6*(1 - (2*M)/R)*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/Power(R,4)) + 
		    b2*((6*M)/Power(R,2) - (3.*M*(1 - (2*M)/R)*(2*M + 2*R))/Power(R,3) - (6.*Power(M,2)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,5) + 
		       (9.*M*(1 - (2*M)/R)*(4*Power(M,2) + 2*M*R + Power(R,2)))/Power(R,4) - 
		       (M*(12*Power(M,2) + (24*Power(M,3))/R + 6*M*R - 27*Power(R,2) + (6*M*(8*Power(M,3) - Power(R,3)))/((1 - (2*M)/R)*Power(R,2)) - 
		            9*Power(R,2)*Log(1 - (2*M)/R)))/(3.*Power(R,4)) + 
		       (4*M*(28*Power(M,3) + 12*Power(M,2)*R + 3*M*Power(R,2) - 9*Power(R,3) + 3*(8*Power(M,3) - Power(R,3))*Log(1 - (2*M)/R) + 
		            24*Power(M,3)*Log(R/M)))/(3.*Power(R,5)) - 
		       Power(1 - (2*M)/R,2)*((6*M)/(Power(1 - (2*M)/R,2)*Power(R,2)) + 
		          (3.*(-2*Power(M,2) - 2*M*R + 3*Power(R,2)))/((1 - (2*M)/R)*Power(R,3)) - 
		          (6.*M*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/(Power(1 - (2*M)/R,2)*Power(R,5)) - 
		          (9.*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/((1 - (2*M)/R)*Power(R,4)) - 
		          (20*Power(M,3) + 30*Power(M,2)*R - 9*M*Power(R,2) - (72*Power(M,4) - 48*Power(M,3)*R)/R + 
		             (2*M*(-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4)))/((1 - (2*M)/R)*Power(R,2)) + 
		             (48*Power(M,3) - 27*M*Power(R,2) + 12*Power(R,3))*Log(1 - (2*M)/R) + 48*Power(M,3)*Log(R/M))/
		           (3.*Power(1 - (2*M)/R,2)*Power(R,4)) + (4*M*
		             (-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		               (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		               (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M)))/(3.*Power(1 - (2*M)/R,3)*Power(R,6)) + 
		          (4*(-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		               (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		               (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M)))/(3.*Power(1 - (2*M)/R,2)*Power(R,5))) - 
		       (4*M*(1 - (2*M)/R)*(-3/(1 - (2*M)/R) + (3.*(12*Power(M,3) - 2*Power(M,2)*R - M*Power(R,2) + Power(R,3)))/
		             ((1 - (2*M)/R)*Power(R,3)) - (-84*Power(M,4) + 20*Power(M,3)*R + 15*Power(M,2)*Power(R,2) - 3*M*Power(R,3) + 
		               (-72*Power(M,4) + 48*Power(M,3)*R - 9*M*Power(R,3) + 3*Power(R,4))*Log(1 - (2*M)/R) - 
		               (72*Power(M,4) - 48*Power(M,3)*R)*Log(R/M))/(3.*Power(1 - (2*M)/R,2)*Power(R,4))))/Power(R,2))))/(mu*(1 - (2*M)/R));

		lm_mode->h_in[0][i] = H1_in;
		lm_mode->h_in[1][i] = H3_in;
		lm_mode->h_in[2][i] = H6_in;

		lm_mode->h_in_rs_deriv[0][i] = f_p*dH1_in;
		lm_mode->h_in_rs_deriv[1][i] = f_p*dH3_in;
		lm_mode->h_in_rs_deriv[2][i] = f_p*dH6_in;
	
		lm_mode->h_out[0][i] = H1_out;
		lm_mode->h_out[1][i] = H3_out;
		lm_mode->h_out[2][i] = H6_out;

		lm_mode->h_out_rs_deriv[0][i] = f_p*dH1_out;
		lm_mode->h_out_rs_deriv[1][i] = f_p*dH3_out;
		lm_mode->h_out_rs_deriv[2][i] = f_p*dH6_out;

		if(i==1){
			printf("at chi = %lg\n", orbit->chi_p[i]);
			printf(" h1_in=%.14f\th3_in=%.14f\th6_in=%.14f\n", H1_in, H3_in, H6_in);
			printf(" h1_out=%.14f\th3_out=%.14f\th6_out=%.14f\n", H1_out, H3_out, H6_out);
			printf(" dh1_in=%.14f\tdh3_in=%.14f\tdh6_in=%.14f\n", dH1_in, dH3_in, dH6_in);
			printf(" dh1_out=%.14f\tdh3_out=%.14f\tdh6_out=%.14f\n", dH1_out, dH3_out, dH6_out);
		}

	}

}
