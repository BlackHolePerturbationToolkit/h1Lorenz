#include "boundary_conditions.h"
#include "libSch.h"
#include "math_funcs.h"
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#define Power(x, y)	(pow((double)(x), (double)(y)))


void monopole_out_bcs(double* a0s, struct boundary_conditions *bcs)
{
	double rs_out 	= 110000.0;

	int converged = 0;

	while(!converged){
		double r_out 	= r_of_rs(rs_out);

		double l 		= bcs->l;
		double omega	= bcs->omega;

	   	int j;    // j = k+3 where k is the actual index used in analytic work.
		int j_max = 5;
		double L = l*(l+1);
		double lambda = (l+2)*(l-1);
		double w = bcs->omega;

		double complex a1[5], a3[5], a6[5], a5[5], a7[5];
		double delta1, delta3, delta6;

		//double complex a1[200], a5[200], a6[200], a7[200], a3[200];
		//double delta1[200], delta3[200], delta5[200], delta6[200], delta7[200];
		
		double complex C1, C5, C7, D1, D3, D5, D7, E1, E3, E5, E7, F3;
	
	
		a1[0] = 0.0;	a6[0] = 0.0; 	
		a1[1] = 0.0; 	a6[1] = 0.0;
		a1[2] = 0.0;	a6[2] = 0.0;
		a5[0] = 0.0;	a7[0] = 0.0;
		a5[1] = 0.0; 	a7[1] = 0.0;
		a5[2] = 0.0;	a7[2] = 0.0;
		a3[0] = 0.0; 	a3[1] = 0.0; a3[2] = 0.0;

		a1[3] = a0s[0];
		a3[3] = a0s[1];  
		a5[3] = 0;
		a6[3] = a0s[2];
		a7[3] = 0;

		double complex sum1 = a1[3]; double complex sum3 = a3[3]; double complex sum5 = a5[3]; double complex sum6 = a6[3]; double complex sum7 = a7[3];
	
		double complex diff_sum1 = 1.0*I*w*a1[3]; double complex diff_sum3 = 1.0*I*w*a3[3]; double complex diff_sum5 = 1.0*I*w*a5[3]; double complex diff_sum6 = 1.0*I*w*a6[3];
		double complex diff_sum7 = 1.0*I*w*a7[3];


		for ( j = 4; j < j_max; j++ )
		{ 
			C1 = (j-4)*(j-3)+4.0*I*w*(j-4)-2.0-L;
			D1 = 2.0*(5.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			D3 = 2.0*( 2.0*(j-5) -8.0 + 4.0*I*w );
			E1 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) -3.0);
			E3 = 8.0*(5.0-2.0*(j-6) );
			F3 = 16.0*(j-7-2.0);
			C5 = (j-4)*(j-3)+4.0*I*w*(j-4)-4.0-L;
			D5 = 2.0*(12.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			E5 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) -8.0);
			C7 = (j-4)*(j-3)+4.0*I*w*(j-4)+2.0-L;
			D7 = 2.0*(-3.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			E7 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) +1.0);
	
			a1[4] = ( C1*a1[3] + (2.0-4.0*I*w)*a3[3] +2.0*a5[3] + 2.0*a6[3] + D1*a1[2] + D3*a3[2] -12.0*a5[2] -20.0*a6[2] + E1*a1[1] + E3*a3[1] + 16.0*a5[1] + 56.0*a6[1] + F3*a3[j-4] - 48.0*a6[0] ) /( 2.0*I*w*(j-3) );
			a3[4] = ( C1*a3[3] + 2.0*(a1[3]- a5[3] -a6[3]) + D1*a3[2] + 4.0*(-a1[2]+ a5[2] +3.0*a6[2]) + E1*a3[1] - 16.0*a6[1] ) /( 2.0*I*w*(j-3) );
			a6[4] = ( C1*a6[3] + 2.0*(a1[3]- a5[3] -a3[3]) + D1*a6[2] + 4.0*(-a1[2]+ a5[2] +3.0*a3[2]) + E1*a6[1] - 16.0*a3[1] ) /( 2.0*I*w*(j-3) );
			a5[4] = ( C5*a5[3] + 2.0*L*(a1[3]-a3[3]-a6[3]) + 2.0*a7[3] + D5*a5[2] + 2.0*L*(-2.0*a1[2] + 4.0*a3[2] + 5.0*a6[2]) -10.0*a7[2] + E5*a5[1] + 4.0*(-2.0*L*a3[1] -3.0*L*a6[1] + 3.0*a7[1]) ) /( 2.0*I*w*(j-3) );
			a7[4] = ( C7*a7[3] +2.0*lambda*a5[3] + D7*a7[2] -4.0*lambda*a5[2] + E7*a7[1] ) /( 2.0*I*w*(j-3) );
	
	
			sum1 +=  a1[4]/pow( r_out, j-3.0);
			sum3 +=  a3[4]/pow( r_out, j-3.0);
			sum6 +=  a6[4]/pow( r_out, j-3.0);
			sum5 +=  a5[4]/pow( r_out, j-3.0);
			sum7 +=  a7[4]/pow( r_out, j-3.0);
			
			diff_sum1 +=  ( I*w*a1[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a1[4]/pow(r_out,j-2.0) );
			diff_sum3 +=  ( I*w*a3[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a3[4]/pow(r_out,j-2.0) );
			diff_sum6 +=  ( I*w*a6[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a6[4]/pow(r_out,j-2.0) );
			diff_sum5 +=  ( I*w*a5[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a5[4]/pow(r_out,j-2.0) );
			diff_sum7 +=  ( I*w*a7[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a7[4]/pow(r_out,j-2.0) );


			delta1 = cabs( a1[4]/pow( r_out, j-3.0) ); delta3 = cabs( a3[4]/pow( r_out, j-3.0) ); delta6 = cabs( a6[4]/pow( r_out, j-3.0) ); 
			//printf("k=%d\td1=%e\td3=%e\td6=%e\n", j-3, delta1[j], delta3[j], delta6[j] );
	
			if ( delta1 >= 1e-14 || delta3 >= 1e-14 || delta6 >= 1e-14 )
			{ j_max = j_max+1; }

			a1[0] = a1[1];  
			a1[1] = a1[2];
			a1[2] = a1[3];
			a1[3] = a1[4];

			a3[0] = a3[1];  
			a3[1] = a3[2];
			a3[2] = a3[3];
			a3[3] = a3[4];

			a5[0] = a5[1];  
			a5[1] = a5[2];
			a5[2] = a5[3];
			a5[3] = a5[4];	

			a6[0] = a6[1];  
			a6[1] = a6[2];
			a6[2] = a6[3];
			a6[3] = a6[4];	

			a7[0] = a7[1];  
			a7[1] = a7[2];
			a7[2] = a7[3];
			a7[3] = a7[4];	

		 }

		if(bcs->omega == 0) break;

		if( !cisnan(sum1) && !cisnan(sum3) && !cisnan(sum5) && !cisnan(sum6) && !cisnan(sum7) && !cisnan(diff_sum1) && !cisnan(diff_sum3) && !cisnan(diff_sum5) && !cisnan(diff_sum6) && !cisnan(diff_sum7) ){
		bcs->R[0] = sum1 * cexp(I*omega*rs_out);
		bcs->R[1] = sum3 * cexp(I*omega*rs_out);
		bcs->R[2] = sum6 * cexp(I*omega*rs_out);

		bcs->R_rs_deriv[0] = diff_sum1 * cexp(I*omega*rs_out);
		bcs->R_rs_deriv[1] = diff_sum3 * cexp(I*omega*rs_out);
		bcs->R_rs_deriv[2] = diff_sum6 * cexp(I*omega*rs_out);

			converged = 1;
		}else{
			rs_out += 5000.0;	//printf("BOUNDARY INCREASES\n\n");
		}
	}

	bcs->rs_boundary 	= rs_out; 

}

void monopole_in_bcs(double* a0s, struct boundary_conditions *bcs)
{

	double rs_in = -50.0;//-10.0/fabs(bcs->omega);;

	double omega = bcs->omega;

	double r_in = r_of_rs(rs_in);

	double epsilon = r_in - 2.0;    // expansion term for inner BC sum
	int j;    // j = k+3 where k is the actual index used in analytic work.
	int j_max = 5;
	double l = bcs->l;
	double L = l*(l+1);
	double lambda = (l+2.0)*(l-1.0);
	double w = bcs->omega;

	double complex b1[5], b3[5], b6[5], b5[5], b7[5];
	double delta1, delta3, delta6;	

	double complex C1, C3, D1, D3, E1, F1, F3, G3, H3, I5, J5, K7;
	
	
  	b1[0] = 0.0;	b6[0] = 0.0; 	
	b1[1] = 0.0; 	b6[1] = 0.0;
	b1[2] = 0.0;	b6[2] = 0.0;
	b5[0] = 0.0;	b7[0] = 0.0;
	b5[1] = 0.0; 	b7[1] = 0.0;
	b5[2] = 0.0;	b7[2] = 0.0;
	b3[0] =0.0; b3[1] = 0.0; b3[2] = 0.0;


	b1[3] = a0s[0]; 
	b3[3] = a0s[1];
	b5[3] = 0;
	b6[3] = a0s[2];
	b7[3] = 0;

	double complex sum1 = b1[3]; double complex sum3 = b3[3]; double complex sum5 = b5[3]; double complex sum6 = b6[3]; double complex sum7 = b7[3];
	
	double complex diff_sum1 = -1.0*I*w*b1[3]; double complex diff_sum3 = -1.0*I*w*b3[3]; double complex diff_sum5 = -1.0*I*w*b5[3]; double complex diff_sum6 = -1.0*I*w*b6[3];
	double complex diff_sum7 = -1.0*I*w*b7[3];
	
	
     	for ( j = 4; j < j_max; j++ )
     	{ 
		C1 = 4.0*(1.0 - L + 3.0*(j-4)*(j-4)-16.0*I*w*(j-4) -(j-4) );
		C3 = 8.0*(2.0*I*w- (j-4) );
		D1 = 2.0*( (j-5)*(3.0*(j-5)-2.0) -2.0*L-1.0 -24.0*I*w*(j-5) );
		D3 = 4.0*(4.0*I*w- (j-5) -1.0 );
		E1 = (j-6)*(j-7)-16.0*I*w*(j-6) -2.0 - L;
		F1 = 2.0*I*w*(j-7);
		F3 = 2.0*I*w*(j-6);
		G3 = 2.0*( (j-4)*( 2.0*(j-4) -1.0 ) - L + 1.0 -12.0*I*w*(j-4) );
		H3 = (j-5)*(j-6)- L - 2.0 -12.0*I*w*(j-5);
		I5 = 2.0*( (j-4)*( 2.0*(j-4) -1.0 ) - L + 4.0 -12.0*I*w*(j-4) );
		J5 = (j-5)*(j-6)- L - 4.0 -12.0*I*w*(j-5);
		K7 = (j-5)*(j-6)- L + 2.0 -12.0*I*w*(j-5);

		b1[4] = ( C1*b1[3]+ C3*b3[3] -8.0*b5[3] + D1*b1[2] + D3*b3[2] -8.0*b6[2] + E1*b1[1] + (2.0+4.0*I*w)*b3[1] +2.0*b5[1] +2.0*b6[1] - F1*b1[0] ) / ( (j-3.0)*8.0*( 4.0*I*w-(j-3) ) ); 
		b3[4] = ( G3*b3[3] + 4.0*(b1[3] -b5[3] +b6[3]) + H3*b3[2] + 2.0*(b1[2] -b5[2] -b6[2]) - F3*b3[1] ) / ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b6[4] = ( G3*b6[3] + 4.0*(b1[3] -b5[3] +b3[3]) + H3*b6[2] + 2.0*(b1[2] -b5[2] -b3[2]) - F3*b6[1] ) / ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b5[4] = ( I5*b5[3] + 2.0*L*(2.0*b1[3] + b6[3]) -2.0*b7[3] + J5*b5[2] + 2.0*L*(b1[2] - b3[2] -b6[2]) +2.0*b7[2] -F3*b5[1] )/ ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b7[4] = ( G3*b7[3] + 4.0*lambda*b5[3] + K7*b7[2] + 2.0*lambda*b5[2] - F3*b7[1] )/ ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) );

		//printf("k=%d   b1= %.14f + I*(%.14f)\tb3 = %.14f +I*(%.14f)\tb6 = %.14f +I*(%.14f)\n", k, creal(b1[4]), cimag(b1[4]), creal(b3[4]), cimag(b3[4]), creal(b6[4]), cimag(b6[4]) );
	
		sum1 = sum1 + b1[4]*pow( epsilon, j-3.0);
		sum3 = sum3 + b3[4]*pow( epsilon, j-3.0);
		sum6 = sum6 + b6[4]*pow( epsilon, j-3.0);
		sum5 = sum5 + b5[4]*pow( epsilon, j-3.0);
		sum7 = sum7 + b7[4]*pow( epsilon, j-3.0);
	
		diff_sum1 = diff_sum1 - 1.0*I*w*b1[4]*pow( epsilon, j-3.0) + (j-3.0)*b1[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum3 = diff_sum3 - 1.0*I*w*b3[4]*pow( epsilon, j-3.0) + (j-3.0)*b3[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum6 = diff_sum6 - 1.0*I*w*b6[4]*pow( epsilon, j-3.0) + (j-3.0)*b6[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum5 = diff_sum5 - 1.0*I*w*b5[4]*pow( epsilon, j-3.0) + (j-3.0)*b5[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum7 = diff_sum7 - 1.0*I*w*b7[4]*pow( epsilon, j-3.0) + (j-3.0)*b7[4]*pow( epsilon, j-3.0)/ r_in;
	

		delta1 = cabs( b1[4]*pow( epsilon, j-3.0) ); delta3 = cabs( b3[4]*pow( epsilon, j-3.0) );  delta6 = cabs( b6[4]*pow( epsilon, j-3.0) ); 
		//printf("%e\t%e\t%e\n", delta1, delta3, delta6);

	        if ( delta1 >= 1e-14 || delta3 >= 1e-14 || delta6 >= 1e-14  )
	        { j_max = j_max+1; }

		b1[0] = b1[1];  
		b1[1] = b1[2];
		b1[2] = b1[3];
		b1[3] = b1[4];

		b3[0] = b3[1];  
		b3[1] = b3[2];
		b3[2] = b3[3];
		b3[3] = b3[4];

		b5[0] = b5[1];  
		b5[1] = b5[2];
		b5[2] = b5[3];
		b5[3] = b5[4];

		b6[0] = b6[1];  
		b6[1] = b6[2];
		b6[2] = b6[3];
		b6[3] = b6[4];

		b7[0] = b7[1];  
		b7[1] = b7[2];
		b7[2] = b7[3];
		b7[3] = b7[4];
    
	//printf("j=%d, k=%d   b1= %.14f + (%.14f) i \t b7 = %.14f + (%.14f) i \n", j, j-3, creal(b1[j]), cimag(b1[j]), creal(b7[j]), cimag(b7[j]) );
	//printf("j=%d, k=%d   b3= %.14f + (%.14f) i \t b5 = %.14f + (%.14f) i \n", j, j-3, creal(b3[j]), cimag(b3[j]), creal(b5[j]), cimag(b5[j]) );
	//printf("b6= %.14f + (%.14f) i \t delta1 = %.14f  \n",  creal(b6[j]), cimag(b6[j]),delta1[j] );

     	}

	bcs->R[0] = sum1 * cexp(-I*omega*rs_in);
	bcs->R[1] = sum3 * cexp(-I*omega*rs_in);
	bcs->R[2] = sum6 * cexp(-I*omega*rs_in);

	bcs->R_rs_deriv[0] = diff_sum1 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[1] = diff_sum3 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[2] = diff_sum6 * cexp(-I*omega*rs_in);


	bcs->rs_boundary 	= rs_in; 
}

void even_dipole_out_bcs(double* a0s, struct boundary_conditions *bcs)
{

	double rs_out 	= 110000.0;
	int converged = 0;
	while(!converged){
	
		double r_out 	= r_of_rs(rs_out);
		int j_max 		= 5;
		double l 		= bcs->l;
		double omega	= bcs->omega;

	   	int j;    // j = k+3 where k is the actual index used in analytic work.

		double L = l*(l+1.0);
		double lambda = (l+2.0)*(l-1.0);
		double w = omega;

 	    	double complex a1[5], a3[5], a5[5], a6[5], a7[5];
		double delta1, delta3, delta5, delta6, delta7;
		
		double complex C1, C5, C7, D1, D3, D5, D7, E1, E3, E5, E7, F3;
	
	
		a1[0] = 0.0;	a6[0] = 0.0; 	
		a1[1] = 0.0; 	a6[1] = 0.0;
		a1[2] = 0.0;	a6[2] = 0.0;
		a5[0] = 0.0;	a7[0] = 0.0;
		a5[1] = 0.0; 	a7[1] = 0.0;
		a5[2] = 0.0;	a7[2] = 0.0;
		a3[0] = 0.0; 	a3[1] = 0.0; a3[2] = 0.0;

		double p = get_orbital_p();

		if(omega < 1e-2){
			if(a0s[0] == 1){
				a1[3] = -1.0;
				a3[3] = 1.0;
				a5[3] = 2.0;
				a6[3] = -4.0;
			}else if(a0s[1] == 1){
				a1[3] = 1.0;
				a3[3] = -1.0;
				a5[3] = 3.0;
				a6[3] = -1.0;
			}else if(a0s[2] == 1){
				a1[3] = 3.0;
				a3[3] = 2.0;
				a5[3] = -1.0;
				a6[3] = 2.0;
			}else if(a0s[3] == 1){
				a1[3] = pow(omega*p, 2.0)*1.0;
				a3[3] = pow(omega*p, 2.0)*-1.0;
				a5[3] = pow(omega*p, 2.0)*-2.0;
				a6[3] = pow(omega*p, 2.0)*-1.0;
			}
		}else{
			a1[3] = a0s[0];
			a3[3] = a0s[1];  
			a5[3] = a0s[2];
			a6[3] = a0s[3];
		}
		a7[3] = 0;



		double complex sum1 = a1[3]; double complex sum3 = a3[3]; double complex sum5 = a5[3]; double complex sum6 = a6[3]; double complex sum7 = a7[3];
	
		double complex diff_sum1 = 1.0*I*w*a1[3]; double complex diff_sum3 = 1.0*I*w*a3[3]; double complex diff_sum5 = 1.0*I*w*a5[3]; double complex diff_sum6 = 1.0*I*w*a6[3];
		double complex diff_sum7 = 1.0*I*w*a7[3];
	
	
		for ( j = 4; j < j_max; j++ )
		{ 
			C1 = (j-4)*(j-3)+4.0*I*w*(j-4)-2.0-L;
			D1 = 2.0*(5.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			D3 = 2.0*( 2.0*(j-5) -8.0 + 4.0*I*w );
			E1 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) -3.0);
			E3 = 8.0*(5.0-2.0*(j-6) );
			F3 = 16.0*(j-7-2.0);
			C5 = (j-4)*(j-3)+4.0*I*w*(j-4)-4.0-L;
			D5 = 2.0*(12.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			E5 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) -8.0);
			C7 = (j-4)*(j-3)+4.0*I*w*(j-4)+2.0-L;
			D7 = 2.0*(-3.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			E7 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) +1.0);

			a1[4] = ( C1*a1[3] + (2.0-4.0*I*w)*a3[3] +2.0*a5[3] + 2.0*a6[3] + D1*a1[2] + D3*a3[2] -12.0*a5[2] -20.0*a6[2] + E1*a1[1] + E3*a3[1] + 16.0*a5[1] + 56.0*a6[1] + F3*a3[0] - 48.0*a6[0] ) /( 2.0*I*w*(j-3) );
			a3[4] = ( C1*a3[3] + 2.0*(a1[3]- a5[3] -a6[3]) + D1*a3[2] + 4.0*(-a1[2]+ a5[2] +3.0*a6[2]) + E1*a3[1] - 16.0*a6[1] ) /( 2.0*I*w*(j-3) );
			a6[4] = ( C1*a6[3] + 2.0*(a1[3]- a5[3] -a3[3]) + D1*a6[2] + 4.0*(-a1[2]+ a5[2] +3.0*a3[2]) + E1*a6[1] - 16.0*a3[1] ) /( 2.0*I*w*(j-3) );
			a5[4] = ( C5*a5[3] + 2.0*L*(a1[3]-a3[3]-a6[3]) + 2.0*a7[3] + D5*a5[2] + 2.0*L*(-2.0*a1[2] + 4.0*a3[2] + 5.0*a6[2]) -10.0*a7[2] + E5*a5[1] + 4.0*(-2.0*L*a3[1] -3.0*L*a6[1] + 3.0*a7[1]) ) /( 2.0*I*w*(j-3) );
			a7[4] = ( C7*a7[3] +2.0*lambda*a5[3] + D7*a7[2] -4.0*lambda*a5[2] + E7*a7[1] ) /( 2.0*I*w*(j-3) );
	
			sum1 = sum1 + a1[4]/pow( r_out, j-3.0);
			sum3 = sum3 + a3[4]/pow( r_out, j-3.0);
			sum5 = sum5 + a5[4]/pow( r_out, j-3.0);
			sum6 = sum6 + a6[4]/pow( r_out, j-3.0);
			sum7 = sum7 + a7[4]/pow( r_out, j-3.0);
	
			diff_sum1 = diff_sum1 + ( I*w*a1[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a1[4]/pow(r_out,j-2.0) );
			diff_sum3 = diff_sum3 + ( I*w*a3[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a3[4]/pow(r_out,j-2.0) );
			diff_sum5 = diff_sum5 + ( I*w*a5[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a5[4]/pow(r_out,j-2.0) );
			diff_sum6 = diff_sum6 + ( I*w*a6[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a6[4]/pow(r_out,j-2.0) );
			diff_sum7 = diff_sum7 + ( I*w*a7[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a7[4]/pow(r_out,j-2.0) );

			delta1 = cabs( a1[4]/pow( r_out, j-3.0) ); delta3 = cabs( a3[4]/pow( r_out, j-3.0) ); delta5 = cabs( a5[4]/pow( r_out, j-3.0) ); delta6 = cabs( a6[4]/pow( r_out, j-3.0) ); delta7 = cabs( a7[4]/pow( r_out, j-3.0) );

			if ( delta1 >= 1e-14 || delta3 >= 1e-14 || delta5 >= 1e-14 || delta6 >= 1e-14 || delta7 >= 1e-14  )
			{ j_max = j_max+1; }

			a1[0] = a1[1];  
			a1[1] = a1[2];
			a1[2] = a1[3];
			a1[3] = a1[4];

			a3[0] = a3[1];  
			a3[1] = a3[2];
			a3[2] = a3[3];
			a3[3] = a3[4];

			a6[0] = a6[1];  
			a6[1] = a6[2];
			a6[2] = a6[3];
			a6[3] = a6[4];

			a5[0] = a5[1];  
			a5[1] = a5[2];
			a5[2] = a5[3];
			a5[3] = a5[4];

			a7[0] = a7[1];  
			a7[1] = a7[2];
			a7[2] = a7[3];
			a7[3] = a7[4];
		

		 }

			if( !cisnan(sum1) && !cisnan(sum3) && !cisnan(sum5) && !cisnan(sum6) && !cisnan(sum7) && !cisnan(diff_sum1) && !cisnan(diff_sum3) && !cisnan(diff_sum5) && !cisnan(diff_sum6) && !cisnan(diff_sum7) ){
				bcs->R[0] = sum1 * cexp(I*omega*rs_out);
				bcs->R[1] = sum3 * cexp(I*omega*rs_out);
				bcs->R[2] = sum5 * cexp(I*omega*rs_out);
				bcs->R[3] = sum6 * cexp(I*omega*rs_out);

				bcs->R_rs_deriv[0] = diff_sum1 * cexp(I*omega*rs_out);
				bcs->R_rs_deriv[1] = diff_sum3 * cexp(I*omega*rs_out);
				bcs->R_rs_deriv[2] = diff_sum5 * cexp(I*omega*rs_out);
				bcs->R_rs_deriv[3] = diff_sum6 * cexp(I*omega*rs_out);

				converged = 1;
			}else{
				rs_out += 5000.0;	//printf("BOUNDARY INCREASES to %lf\n\n", rs_out);
			}
	}

	bcs->rs_boundary 	= rs_out; 

}

void even_dipole_in_bcs(double* a0s, struct boundary_conditions *bcs)
{
	double rs_in = -50.00;//-10.0/fabs(bcs->omega);

	double omega = bcs->omega;

	double r_in = r_of_rs(rs_in);

	double epsilon = r_in - 2.0;    // expansion term for inner BC sum
	int j;    // j = k+3 where k is the actual index used in analytic work.
	int j_max = 5;
	double l = bcs->l;
	double L = l*(l+1);
	double lambda = (l+2.0)*(l-1.0);
	double w = bcs->omega;

     	double complex b1[5], b3[5], b5[5], b6[5], b7[5];
	double delta1, delta3, delta5, delta6, delta7;
        
	double complex C1, C3, D1, D3, E1, F1, F3, G3, H3, I5, J5, K7;
	
	
   	b1[0] = 0.0;	b6[0] = 0.0; 	
	b1[1] = 0.0; 	b6[1] = 0.0;
	b1[2] = 0.0;	b6[2] = 0.0;
	b5[0] = 0.0;	b7[0] = 0.0;
	b5[1] = 0.0; 	b7[1] = 0.0;
	b5[2] = 0.0;	b7[2] = 0.0;
	b3[0] =0.0; b3[1] = 0.0; b3[2] = 0.0;

	
	// initial value for b3_{k=0} i.e. b3[3] comes from the gauge condition


	b1[3] = a0s[0]; 
	b3[3] = a0s[1];
	b5[3] = a0s[2];
	b6[3] = a0s[3];
	b7[3] = 0;

	double complex sum1 = b1[3]; double complex sum3 = b3[3]; double complex sum5 = b5[3]; double complex sum6 = b6[3]; double complex sum7 = b7[3];
	
	double complex diff_sum1 = -1.0*I*w*b1[3]; double complex diff_sum3 = -1.0*I*w*b3[3]; double complex diff_sum5 = -1.0*I*w*b5[3]; double complex diff_sum6 = -1.0*I*w*b6[3];
	double complex diff_sum7 = -1.0*I*w*b7[3];
	
	
     	for ( j = 4; j < j_max; j++ )
  	{ 
		C1 = 4.0*(1.0 - L + 3.0*(j-4)*(j-4)-16.0*I*w*(j-4) -(j-4) );
		C3 = 8.0*(2.0*I*w- (j-4) );
		D1 = 2.0*( (j-5)*(3.0*(j-5)-2.0) -2.0*L-1.0 -24.0*I*w*(j-5) );
		D3 = 4.0*(4.0*I*w- (j-5) -1.0 );
		E1 = (j-6)*(j-7)-16.0*I*w*(j-6) -2.0 - L;
		F1 = 2.0*I*w*(j-7);
		F3 = 2.0*I*w*(j-6);
		G3 = 2.0*( (j-4)*( 2.0*(j-4) -1.0 ) - L + 1.0 -12.0*I*w*(j-4) );
		H3 = (j-5)*(j-6)- L - 2.0 -12.0*I*w*(j-5);
		I5 = 2.0*( (j-4)*( 2.0*(j-4) -1.0 ) - L + 4.0 -12.0*I*w*(j-4) );
		J5 = (j-5)*(j-6)- L - 4.0 -12.0*I*w*(j-5);
		K7 = (j-5)*(j-6)- L + 2.0 -12.0*I*w*(j-5);

 
	
		b1[4] = ( C1*b1[3]+ C3*b3[3] -8.0*b5[3] + D1*b1[2] + D3*b3[2] -8.0*b6[2] + E1*b1[1] + (2.0+4.0*I*w)*b3[1] +2.0*b5[1] +2.0*b6[1] - F1*b1[0] ) / ( (j-3.0)*8.0*( 4.0*I*w-(j-3) ) ); 
		b3[4] = ( G3*b3[3] + 4.0*(b1[3] -b5[3] +b6[3]) + H3*b3[2] + 2.0*(b1[2] -b5[2] -b6[2]) - F3*b3[1] ) / ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b6[4] = ( G3*b6[3] + 4.0*(b1[3] -b5[3] +b3[3]) + H3*b6[2] + 2.0*(b1[2] -b5[2] -b3[2]) - F3*b6[1] ) / ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b5[4] = ( I5*b5[3] + 2.0*L*(2.0*b1[3] + b6[3]) -2.0*b7[3] + J5*b5[2] + 2.0*L*(b1[2] - b3[2] -b6[2]) +2.0*b7[2] -F3*b5[1] )/ ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b7[4] = ( G3*b7[3] + 4.0*lambda*b5[3] + K7*b7[2] + 2.0*lambda*b5[2] - F3*b7[1] )/ ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) );
	
		sum1 = sum1 + b1[4]*pow( epsilon, j-3.0);
		sum3 = sum3 + b3[4]*pow( epsilon, j-3.0);
		sum5 = sum5 + b5[4]*pow( epsilon, j-3.0);
		sum6 = sum6 + b6[4]*pow( epsilon, j-3.0);
		sum7 = sum7 + b7[4]*pow( epsilon, j-3.0);
	
		diff_sum1 = diff_sum1 - 1.0*I*w*b1[4]*pow( epsilon, j-3.0) + (j-3.0)*b1[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum3 = diff_sum3 - 1.0*I*w*b3[4]*pow( epsilon, j-3.0) + (j-3.0)*b3[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum5 = diff_sum5 - 1.0*I*w*b5[4]*pow( epsilon, j-3.0) + (j-3.0)*b5[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum6 = diff_sum6 - 1.0*I*w*b6[4]*pow( epsilon, j-3.0) + (j-3.0)*b6[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum7 = diff_sum7 - 1.0*I*w*b7[4]*pow( epsilon, j-3.0) + (j-3.0)*b7[4]*pow( epsilon, j-3.0)/ r_in;
	
		delta1 = cabs( b1[4]*pow( epsilon, j-3.0) ); delta3 = cabs( b3[4]*pow( epsilon, j-3.0) ); delta5 = cabs( b5[4]*pow( epsilon, j-3.0) ); delta6 = cabs( b6[4]*pow( epsilon, j-3.0) ); delta7 = cabs( b7[4]*pow( epsilon, j-3.0) ); 

		//printf("k=%d   b1= %.14f + I*(%.14f)\tb3 = %.14f +I*(%.14f)\tb5 = %.14f +I*(%.14f)\tb6 = %.14f +I*(%.14f)\tb7 = %.14f +I*(%.14f)\n", j-3, creal(b1[4]), cimag(b1[4]), creal(b3[4]), cimag(b3[4]), creal(b5[4]), cimag(b5[4]), creal(b6[4]), cimag(b6[4]), creal(b7[4]), cimag(b7[4]) );

	       if ( delta1 >= 1e-14 || delta3 >= 1e-14 || delta5 >= 1e-14 || delta6 >= 1e-14 || delta7 >= 1e-14  )
	       { j_max = j_max+1; }

		b1[0] = b1[1];  
		b1[1] = b1[2];
		b1[2] = b1[3];
		b1[3] = b1[4];

		b3[0] = b3[1];  
		b3[1] = b3[2];
		b3[2] = b3[3];
		b3[3] = b3[4];

		b6[0] = b6[1];  
		b6[1] = b6[2];
		b6[2] = b6[3];
		b6[3] = b6[4];

		b5[0] = b5[1];  
		b5[1] = b5[2];
		b5[2] = b5[3];
		b5[3] = b5[4];

		b7[0] = b7[1];  
		b7[1] = b7[2];
		b7[2] = b7[3];
		b7[3] = b7[4];
  
		//printf("j=%d, k=%d   b1= %.14f + (%.14f) i \t b7 = %.14f + (%.14f) i \n", j, j-3, creal(b1[j]), cimag(b1[j]), creal(b7[j]), cimag(b7[j]) );
		//printf("j=%d, k=%d   b3= %.14f + (%.14f) i \t b5 = %.14f + (%.14f) i \n", j, j-3, creal(b3[j]), cimag(b3[j]), creal(b5[j]), cimag(b5[j]) );
		//printf("b6= %.14f + (%.14f) i \t delta1 = %.14f  \n",  creal(b6[j]), cimag(b6[j]),delta1[j] );

     }

	bcs->R[0] = sum1 * cexp(-I*omega*rs_in);
	bcs->R[1] = sum3 * cexp(-I*omega*rs_in);
	bcs->R[2] = sum5 * cexp(-I*omega*rs_in);
	bcs->R[3] = sum6 * cexp(-I*omega*rs_in);

	bcs->R_rs_deriv[0] = diff_sum1 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[1] = diff_sum3 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[2] = diff_sum5 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[3] = diff_sum6 * cexp(-I*omega*rs_in);

	bcs->rs_boundary 	= rs_in; 

	/*printf("\n%.12e+I%12e %.12e+I%12e %.12e+I%12e %.12e+I%12e %.12e+I%12e %.12e+I%12e %.12e+I%12e %.12e+I%12e %e\n", 
		creal(bcs->R[0]), cimag(bcs->R[0]), 
		creal(bcs->R[1]), cimag(bcs->R[1]), 
		creal(bcs->R[2]), cimag(bcs->R[2]), 
		creal(bcs->R[3]), cimag(bcs->R[3]), 
		creal(bcs->R_rs_deriv[0]), cimag(bcs->R_rs_deriv[0]), 
		creal(bcs->R_rs_deriv[1]), cimag(bcs->R_rs_deriv[1]), 
		creal(bcs->R_rs_deriv[2]), cimag(bcs->R_rs_deriv[2]), 
		creal(bcs->R_rs_deriv[3]), cimag(bcs->R_rs_deriv[3]),
		rs_in);*/
}


void even_out_bcs(double* a0s, struct boundary_conditions *bcs)
{
	double omega = bcs->omega;
	
	double l 		= bcs->l;

	double rs_out 	= 110000.0;//10.0/omega;

	int converged = 0;
	int increase_level = 0;
	
	while(!converged){
		double r_out 	= r_of_rs(rs_out);



	   	int j;    // j = k+3 where k is the actual index used in analytic work.
		int j_max = 5;
		double L = l*(l+1.0);
		double lambda = (l+2.0)*(l-1.0);
		double w = bcs->omega;

     		double complex a1[5], a3[5], a5[5], a6[5], a7[5];
		double delta1, delta3, delta5, delta6, delta7;
		
		double complex C1, C5, C7, D1, D3, D5, D7, E1, E3, E5, E7, F3;

		double p = get_orbital_p();	
	
		a1[0] = 0.0;	a6[0] = 0.0; 	
		a1[1] = 0.0; 	a6[1] = 0.0;
		a1[2] = 0.0;	a6[2] = 0.0;
		a5[0] = 0.0;	a7[0] = 0.0;
		a5[1] = 0.0; 	a7[1] = 0.0;
		a5[2] = 0.0;	a7[2] = 0.0;
		a3[0] = 0.0; 	a3[1] = 0.0; a3[2] = 0.0;


		//if(omega < 1e-2){

			double a_data[] = { -(2.0+L), 2.00, 2.00, 2.00, 0.00,
						       2.0, -(2.0+L), -2.0, -2.0, 0.0,
						       2.0*L, -2.0*L, -(4.0+L), -2.0*L, 2.0,
						       2.0, -2.0, -2.0, -(2.0+L), 0.0,
							   0.00, 0.00, 2.0*lambda, 0.00, -(L-2.0) };

			gsl_matrix_view M = gsl_matrix_view_array (a_data, 5, 5);

			gsl_matrix_transpose(&M.matrix);

		    gsl_vector_complex *eigenvalues 	= gsl_vector_complex_alloc (5);
		    gsl_matrix_complex *eigenvectors 	= gsl_matrix_complex_alloc (5, 5);
		    gsl_eigen_nonsymmv_workspace * ws 	= gsl_eigen_nonsymmv_alloc (5);
		    gsl_eigen_nonsymmv (&M.matrix, eigenvalues, eigenvectors, ws);
		    gsl_eigen_nonsymmv_free (ws);
		    gsl_eigen_nonsymmv_sort (eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);

			int highest_power	= (int)round((GSL_REAL(gsl_vector_complex_get (eigenvalues, 0))));
			int middle_power 	= (int)round((GSL_REAL(gsl_vector_complex_get (eigenvalues, 1))));
			int lowest_power 	= (int)round((GSL_REAL(gsl_vector_complex_get (eigenvalues, 4))));


			int highest_l 		= (int)round((-1.0 + sqrt(1.0 + 4.0*abs(highest_power)))/2.0);
			int middle_l		= (int)round((-1.0 + sqrt(1.0 + 4.0*abs(middle_power)))/2.0);
			int lowest_l		= (int)round((-1.0 + sqrt(1.0 + 4.0*abs(lowest_power)))/2.0);


			gsl_matrix_complex_transpose(eigenvectors);


			int s;

			gsl_permutation * perm = gsl_permutation_alloc (5);
			gsl_linalg_complex_LU_decomp (eigenvectors, perm, &s);

			gsl_matrix_complex *inverse_eigen = gsl_matrix_complex_alloc(5,5);

			gsl_linalg_complex_LU_invert(eigenvectors, perm, inverse_eigen);

			if(a0s[0] == 1){
				a1[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 0, 0))*pow(omega*p, (highest_l - middle_l));
				a3[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 1, 0))*pow(omega*p, (highest_l - middle_l));
				a5[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 2, 0))*pow(omega*p, (highest_l - middle_l));
				a6[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 3, 0))*pow(omega*p, (highest_l - middle_l));
				a7[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 4, 0))*pow(omega*p, (highest_l - middle_l));
			}else if(a0s[1] == 1){
				a1[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 0, 1));
				a3[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 1, 1));
				a5[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 2, 1));
				a6[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 3, 1));
				a7[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 4, 1));
			}else if(a0s[2] == 1){
				a1[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 0, 2));
				a3[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 1, 2));
				a5[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 2, 2));
				a6[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 3, 2));
				a7[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 4, 2));
			}else if(a0s[3] == 1){
				a1[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 0, 3));						
				a3[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 1, 3));
				a5[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 2, 3));
				a6[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 3, 3));
				a7[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 4, 3));
			}else if(a0s[4] == 1){
				a1[3] =	GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 0, 4))*pow(omega*p, (lowest_l - middle_l));										
				a3[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 1, 4))*pow(omega*p, (lowest_l - middle_l));
				a5[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 2, 4))*pow(omega*p, (lowest_l - middle_l));
				a6[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 3, 4))*pow(omega*p, (lowest_l - middle_l));
				a7[3] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 4, 4))*pow(omega*p, (lowest_l - middle_l));
			}		

			gsl_permutation_free(perm);
			gsl_vector_complex_free(eigenvalues);
			gsl_matrix_complex_free(eigenvectors);
			gsl_matrix_complex_free(inverse_eigen);
			//}else{
			// a1[3] = a0s[0];
// 			a3[3] = a0s[1];
// 			a5[3] = a0s[2];
// 			a6[3] = a0s[3];
// 			a7[3] = a0s[4];
			//}
			
		// Store the rotated basis values in a0s
			a0s[0] = a1[3];
			a0s[1] = a3[3];
			a0s[2] = a5[3];
			a0s[3] = a6[3];
			a0s[4] = a7[3];

		double complex sum1 = a1[3]; double complex sum3 = a3[3]; double complex sum5 = a5[3]; double complex sum6 = a6[3]; double complex sum7 = a7[3];
	
		double complex diff_sum1 = 1.0*I*w*a1[3]; double complex diff_sum3 = 1.0*I*w*a3[3]; double complex diff_sum5 = 1.0*I*w*a5[3]; double complex diff_sum6 = 1.0*I*w*a6[3];
		double complex diff_sum7 = 1.0*I*w*a7[3];
	
	
		for ( j = 4; j < j_max; j++ )
		{ 
			C1 = (j-4)*(j-3)+4.0*I*w*(j-4)-2.0-L;
			D1 = 2.0*(5.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			D3 = 2.0*( 2.0*(j-5) -8.0 + 4.0*I*w );
			E1 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) -3.0);
			E3 = 8.0*(5.0-2.0*(j-6) );
			F3 = 16.0*(j-7-2.0);
			C5 = (j-4)*(j-3)+4.0*I*w*(j-4)-4.0-L;
			D5 = 2.0*(12.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			E5 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) -8.0);
			C7 = (j-4)*(j-3)+4.0*I*w*(j-4)+2.0-L;
			D7 = 2.0*(-3.0 + L -2.0*(j-5)*(j-5)-3.0*(j-5) );
			E7 = 4.0*( (j-6)*(j-6) + 2.0*(j-6) +1.0);

			a1[4] = ( C1*a1[3] + (2.0-4.0*I*w)*a3[3] +2.0*a5[3] + 2.0*a6[3] + D1*a1[2] + D3*a3[2] -12.0*a5[2] -20.0*a6[2] + E1*a1[1] + E3*a3[1] + 16.0*a5[1] + 56.0*a6[1] + F3*a3[0] - 48.0*a6[0] ) /( 2.0*I*w*(j-3) );
			a3[4] = ( C1*a3[3] + 2.0*(a1[3]- a5[3] -a6[3]) + D1*a3[2] + 4.0*(-a1[2]+ a5[2] +3.0*a6[2]) + E1*a3[1] - 16.0*a6[1] ) /( 2.0*I*w*(j-3) );
			a6[4] = ( C1*a6[3] + 2.0*(a1[3]- a5[3] -a3[3]) + D1*a6[2] + 4.0*(-a1[2]+ a5[2] +3.0*a3[2]) + E1*a6[1] - 16.0*a3[1] ) /( 2.0*I*w*(j-3) );
			a5[4] = ( C5*a5[3] + 2.0*L*(a1[3]-a3[3]-a6[3]) + 2.0*a7[3] + D5*a5[2] + 2.0*L*(-2.0*a1[2] + 4.0*a3[2] + 5.0*a6[2]) -10.0*a7[2] + E5*a5[1] + 4.0*(-2.0*L*a3[1] -3.0*L*a6[1] + 3.0*a7[1]) ) /( 2.0*I*w*(j-3) );
			a7[4] = ( C7*a7[3] +2.0*lambda*a5[3] + D7*a7[2] -4.0*lambda*a5[2] + E7*a7[1] ) /( 2.0*I*w*(j-3) );
	
			sum1 = sum1 + a1[4]/pow( r_out, j-3.0);
			sum3 = sum3 + a3[4]/pow( r_out, j-3.0);
			sum5 = sum5 + a5[4]/pow( r_out, j-3.0);
			sum6 = sum6 + a6[4]/pow( r_out, j-3.0);
			sum7 = sum7 + a7[4]/pow( r_out, j-3.0);
	
			diff_sum1 = diff_sum1 + ( I*w*a1[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a1[4]/pow(r_out,j-2.0) );
			diff_sum3 = diff_sum3 + ( I*w*a3[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a3[4]/pow(r_out,j-2.0) );
			diff_sum5 = diff_sum5 + ( I*w*a5[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a5[4]/pow(r_out,j-2.0) );
			diff_sum6 = diff_sum6 + ( I*w*a6[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a6[4]/pow(r_out,j-2.0) );
			diff_sum7 = diff_sum7 + ( I*w*a7[4]/pow(r_out, j-3.0) - (1.0-2.0/r_out)*(j-3.0)*a7[4]/pow(r_out,j-2.0) );


			delta1 = cabs( a1[4]/pow( r_out, j-3.0) ); delta3 = cabs( a3[4]/pow( r_out, j-3.0) ); delta5 = cabs( a5[4]/pow( r_out, j-3.0) ); delta6 = cabs( a6[4]/pow( r_out, j-3.0) ); delta7 = cabs( a7[4]/pow( r_out, j-3.0) );

		       if ( delta1 >= 1e-14 || delta3 >= 1e-14 || delta5 >= 1e-14 || delta6 >= 1e-14 || delta7 >= 1e-14  )
		       { j_max = j_max+1; }	

			a1[0] = a1[1];  
			a1[1] = a1[2];
			a1[2] = a1[3];
			a1[3] = a1[4];

			a3[0] = a3[1];  
			a3[1] = a3[2];
			a3[2] = a3[3];
			a3[3] = a3[4];

			a6[0] = a6[1];  
			a6[1] = a6[2];
			a6[2] = a6[3];
			a6[3] = a6[4];

			a5[0] = a5[1];  
			a5[1] = a5[2];
			a5[2] = a5[3];
			a5[3] = a5[4];

			a7[0] = a7[1];  
			a7[1] = a7[2];
			a7[2] = a7[3];
			a7[3] = a7[4];

		}

		if(omega == 0) break;

		if( !cisnan(sum1) && !cisnan(sum3) && !cisnan(sum5) && !cisnan(sum6) && !cisnan(sum7) && !cisnan(diff_sum1) && !cisnan(diff_sum3) && !cisnan(diff_sum5) && !cisnan(diff_sum6) && !cisnan(diff_sum7) ){
			bcs->R[0] = sum1 * cexp(I*omega*rs_out);
			bcs->R[1] = sum3 * cexp(I*omega*rs_out);
			bcs->R[2] = sum5 * cexp(I*omega*rs_out);
			bcs->R[3] = sum6 * cexp(I*omega*rs_out);
			bcs->R[4] = sum7 * cexp(I*omega*rs_out);

			bcs->R_rs_deriv[0] = diff_sum1 * cexp(I*omega*rs_out);
			bcs->R_rs_deriv[1] = diff_sum3 * cexp(I*omega*rs_out);
			bcs->R_rs_deriv[2] = diff_sum5 * cexp(I*omega*rs_out);
			bcs->R_rs_deriv[3] = diff_sum6 * cexp(I*omega*rs_out);
			bcs->R_rs_deriv[4] = diff_sum7 * cexp(I*omega*rs_out);

			converged = 1;
		}else{
			if(increase_level == 0){
				rs_out += 500.0;
				increase_level = 1;	
			}
			else if(increase_level == 1){
				rs_out += 1000.0;
				increase_level = 2;
			}else{
				rs_out += 5000.0;
			}
			//printf("BOUNDARY INCREASES\n\n");
		}
	}


	bcs->rs_boundary 	= rs_out; 

}

void even_in_bcs(double* a0s, struct boundary_conditions *bcs)
{

	double rs_in = -50.0;//-10.0/fabs(bcs->omega);;

	double omega = bcs->omega;

	double r_in = r_of_rs(rs_in);
	double epsilon = r_in - 2.0;    // expansion term for inner BC sum
	int j;    // j = k+3 where k is the actual index used in analytic work.
	int j_max = 5;
	double l = bcs->l;
	double L = l*(l+1);
	double lambda = (l+2.0)*(l-1.0);
	double w = bcs->omega;

    	double complex b1[5], b3[5], b5[5], b6[5], b7[5];
	double delta1, delta3, delta5, delta6, delta7;
    
	double complex C1, C3, D1, D3, E1, F1, F3, G3, H3, I5, J5, K7;
	
	
   	b1[0] = 0.0;	b6[0] = 0.0; 	
	b1[1] = 0.0; 	b6[1] = 0.0;
	b1[2] = 0.0;	b6[2] = 0.0;
	b5[0] = 0.0;	b7[0] = 0.0;
	b5[1] = 0.0; 	b7[1] = 0.0;
	b5[2] = 0.0;	b7[2] = 0.0;
	b3[0] =0.0; b3[1] = 0.0; b3[2] = 0.0;

	
	// initial value for b3_{k=0} i.e. b3[3] comes from the gauge condition
	b1[3] = a0s[0]; 
	b3[3] = a0s[1];
	b5[3] = a0s[2];
	b6[3] = a0s[3];
	b7[3] = a0s[4];

	double complex sum1 = b1[3]; double complex sum3 = b3[3]; double complex sum5 = b5[3]; double complex sum6 = b6[3]; double complex sum7 = b7[3];
	
	double complex diff_sum1 = -1.0*I*w*b1[3]; double complex diff_sum3 = -1.0*I*w*b3[3]; double complex diff_sum5 = -1.0*I*w*b5[3]; double complex diff_sum6 = -1.0*I*w*b6[3];
	double complex diff_sum7 = -1.0*I*w*b7[3];
	
	
     	for ( j = 4; j < j_max; j++ )
     	{ 
		C1 = 4.0*(1.0 - L + 3.0*(j-4)*(j-4)-16.0*I*w*(j-4) -(j-4) );
		C3 = 8.0*(2.0*I*w- (j-4) );
		D1 = 2.0*( (j-5)*(3.0*(j-5)-2.0) -2.0*L-1.0 -24.0*I*w*(j-5) );
		D3 = 4.0*(4.0*I*w- (j-5) -1.0 );
		E1 = (j-6)*(j-7)-16.0*I*w*(j-6) -2.0 - L;
		F1 = 2.0*I*w*(j-7);
		F3 = 2.0*I*w*(j-6);
		G3 = 2.0*( (j-4)*( 2.0*(j-4) -1.0 ) - L + 1.0 -12.0*I*w*(j-4) );
		H3 = (j-5)*(j-6)- L - 2.0 -12.0*I*w*(j-5);
		I5 = 2.0*( (j-4)*( 2.0*(j-4) -1.0 ) - L + 4.0 -12.0*I*w*(j-4) );
		J5 = (j-5)*(j-6)- L - 4.0 -12.0*I*w*(j-5);
		K7 = (j-5)*(j-6)- L + 2.0 -12.0*I*w*(j-5);
	
		b1[4] = ( C1*b1[3]+ C3*b3[3] -8.0*b5[3] + D1*b1[2] + D3*b3[2] -8.0*b6[2] + E1*b1[1] + (2.0+4.0*I*w)*b3[1] +2.0*b5[1] +2.0*b6[1] - F1*b1[0] ) / ( (j-3.0)*8.0*( 4.0*I*w-(j-3) ) ); 
		b3[4] = ( G3*b3[3] + 4.0*(b1[3] -b5[3] +b6[3]) + H3*b3[2] + 2.0*(b1[2] -b5[2] -b6[2]) - F3*b3[1] ) / ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b6[4] = ( G3*b6[3] + 4.0*(b1[3] -b5[3] +b3[3]) + H3*b6[2] + 2.0*(b1[2] -b5[2] -b3[2]) - F3*b6[1] ) / ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b5[4] = ( I5*b5[3] + 2.0*L*(2.0*b1[3] + b6[3]) -2.0*b7[3] + J5*b5[2] + 2.0*L*(b1[2] - b3[2] -b6[2]) +2.0*b7[2] -F3*b5[1] )/ ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) ); 
		b7[4] = ( G3*b7[3] + 4.0*lambda*b5[3] + K7*b7[2] + 2.0*lambda*b5[2] - F3*b7[1] )/ ( (j-3.0)*4.0*( 4.0*I*w-(j-3) ) );
	
		sum1 = sum1 + b1[4]*pow( epsilon, j-3.0);
		sum3 = sum3 + b3[4]*pow( epsilon, j-3.0);
		sum5 = sum5 + b5[4]*pow( epsilon, j-3.0);
		sum6 = sum6 + b6[4]*pow( epsilon, j-3.0);
		sum7 = sum7 + b7[4]*pow( epsilon, j-3.0);
	
		diff_sum1 = diff_sum1 - 1.0*I*w*b1[4]*pow( epsilon, j-3.0) + (j-3.0)*b1[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum3 = diff_sum3 - 1.0*I*w*b3[4]*pow( epsilon, j-3.0) + (j-3.0)*b3[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum5 = diff_sum5 - 1.0*I*w*b5[4]*pow( epsilon, j-3.0) + (j-3.0)*b5[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum6 = diff_sum6 - 1.0*I*w*b6[4]*pow( epsilon, j-3.0) + (j-3.0)*b6[4]*pow( epsilon, j-3.0)/ r_in;
		diff_sum7 = diff_sum7 - 1.0*I*w*b7[4]*pow( epsilon, j-3.0) + (j-3.0)*b7[4]*pow( epsilon, j-3.0)/ r_in;

		delta1 = cabs( b1[4]*pow( epsilon, j-3.0) ); delta3 = cabs( b3[4]*pow( epsilon, j-3.0) ); delta5 = cabs( b5[4]*pow( epsilon, j-3.0) ); delta6 = cabs( b6[4]*pow( epsilon, j-3.0) ); delta7 = cabs( b7[4]*pow( epsilon, j-3.0) ); 

		//printf("k=%d   b1= %.14f + I*(%.14f)\tb3 = %.14f +I*(%.14f)\tb5 = %.14f +I*(%.14f)\tb6 = %.14f +I*(%.14f)\tb7 = %.14f +I*(%.14f)\n", j-3, creal(b1[4]), cimag(b1[4]), creal(b3[4]), cimag(b3[4]), creal(b5[4]), cimag(b5[4]), creal(b6[4]), cimag(b6[4]), creal(b7[4]), cimag(b7[4]) );

	       if ( delta1 >= 1e-14 || delta3 >= 1e-14 || delta5 >= 1e-14 || delta6 >= 1e-14 || delta7 >= 1e-14  )
	       { j_max = j_max+1; }

		b1[0] = b1[1];  
		b1[1] = b1[2];
		b1[2] = b1[3];
		b1[3] = b1[4];

		b3[0] = b3[1];  
		b3[1] = b3[2];
		b3[2] = b3[3];
		b3[3] = b3[4];

		b6[0] = b6[1];  
		b6[1] = b6[2];
		b6[2] = b6[3];
		b6[3] = b6[4];

		b5[0] = b5[1];  
		b5[1] = b5[2];
		b5[2] = b5[3];
		b5[3] = b5[4];

		b7[0] = b7[1];  
		b7[1] = b7[2];
		b7[2] = b7[3];
		b7[3] = b7[4];

     }

	bcs->R[0] = sum1 * cexp(-I*omega*rs_in);
	bcs->R[1] = sum3 * cexp(-I*omega*rs_in);
	bcs->R[2] = sum5 * cexp(-I*omega*rs_in);
	bcs->R[3] = sum6 * cexp(-I*omega*rs_in);
	bcs->R[4] = sum7 * cexp(-I*omega*rs_in);

	bcs->R_rs_deriv[0] = diff_sum1 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[1] = diff_sum3 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[2] = diff_sum5 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[3] = diff_sum6 * cexp(-I*omega*rs_in);
	bcs->R_rs_deriv[4] = diff_sum7 * cexp(-I*omega*rs_in);

	bcs->rs_boundary 	= rs_in; 
}

void odd_out_bcs(double* a0s, struct boundary_conditions *bcs)
{
	double l = bcs->l;
	double w = bcs->omega;
	double omega = bcs->omega;

	double rs_out = 110000.0;
	int converged = 0;

	while(!converged){

		double r_out = r_of_rs(rs_out);

		int j, k;    // j = k+2 where k is the actual index used in analytic work.
		int j_max = 5;
		double L = l*(l+1.0);
		double lambda = (l+2.0)*(l-1.0);

	     	double complex a9[4];
		double complex a10[4];  
		double delta9, delta10; 

		double complex C, D, E, Ii, J, K;
		double F, G, H, N, P;

		F = 2.0;
		G = -10.0;
		H = 12.0;
		N = 2.0*lambda;
		P = -4.0*lambda;

		a9[0] = 0.0;    a10[0] = 0.0;
		a9[1] = 0.0;     a10[1] = 0.0;


		a9[2] = a0s[0];  // this a9_0 = 1.0 for k=0
		a10[2] = a0s[1];

		if(omega != 0){
			double a_data[] = { -(L+4.0), 2.0,
								2.0*lambda, -(L - 2.0)};

			gsl_matrix_view M = gsl_matrix_view_array (a_data, 2, 2);

			gsl_matrix_transpose(&M.matrix);

		    gsl_vector_complex *eigenvalues 	= gsl_vector_complex_alloc (2);
		    gsl_matrix_complex *eigenvectors 	= gsl_matrix_complex_alloc (2, 2);
		    gsl_eigen_nonsymmv_workspace * ws 	= gsl_eigen_nonsymmv_alloc (2);
		    gsl_eigen_nonsymmv (&M.matrix, eigenvalues, eigenvectors, ws);
		    gsl_eigen_nonsymmv_free (ws);
		    gsl_eigen_nonsymmv_sort (eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);
			//printf("\n");

			int highest_power	= (int)round((GSL_REAL(gsl_vector_complex_get (eigenvalues, 0))));
			int lowest_power 	= (int)round((GSL_REAL(gsl_vector_complex_get (eigenvalues, 1))));

	
			int highest_l 		= (int)round((-1.0 + sqrt(1.0 + 4.0*abs(highest_power)))/2.0);
			int lowest_l		= (int)round((-1.0 + sqrt(1.0 + 4.0*abs(lowest_power)))/2.0);


			gsl_matrix_complex_transpose(eigenvectors);

			int s;

			gsl_permutation * perm = gsl_permutation_alloc (2);
			gsl_linalg_complex_LU_decomp (eigenvectors, perm, &s);

			gsl_matrix_complex *inverse_eigen = gsl_matrix_complex_alloc(2,2);

			gsl_linalg_complex_LU_invert(eigenvectors, perm, inverse_eigen);

			double p = get_orbital_p();

			if(a0s[0] == 1){
				a9[2] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 0, 0))*pow(omega*p, highest_l - lowest_l);
				a10[2] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 1, 0))*pow(omega*p, highest_l - lowest_l);
			}else if(a0s[1] == 1){
				a9[2] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 0, 1));
				a10[2] = GSL_REAL(gsl_matrix_complex_get(inverse_eigen, 1, 1));
			}

			gsl_permutation_free(perm);
			gsl_vector_complex_free(eigenvalues);
			gsl_matrix_complex_free(eigenvectors);
			gsl_matrix_complex_free(inverse_eigen);

		}
		 
		// Store the complex amplitude of the homogeneous fields
		a0s[0] = a9[2];
		a0s[1] = a10[2];
     
		double complex sum9 = a9[2];
		double complex sum10 = a10[2];
		double complex diff_sum9 = I*w*a9[2];
		double complex diff_sum10 = I*w*a10[2];


		for ( j = 3; j <= j_max; j++ )
		{
			k = j - 2;
			C = ( 4.0*I*w*(k-1) + (k-1)*k - L - 4.0 );
			D = -6.0*(k-2)-4.0*(k-2)*(k-2)+24.0+2.0*L;
			E = 4.0*( (k-3)*(k-3)+2.0*(k-3)-8.0 );
			Ii = 4.0*I*w*(k-1) + (k-1)*k - L + 2.0;
			J = -6.0*(k-2) - 4.0*(k-2)*(k-2) - 6.0 + 2.0*L;
			K = 4.0*( (k-3)*(k-3) + 2.0*(k-3) + 1.0 );

			a9[3] = ( C*a9[2]+D*a9[1]+E*a9[0]+F*a10[2]+G*a10[1] + H*a10[0] ) / (2.0*I*w*k ); // recursion relation for the b9_k's
			a10[3] = ( Ii*a10[2]+J*a10[1]+K*a10[0]+N*a9[2]+P*a9[1] ) / (2.0*I*w*k );

			sum9 += a9[3] / pow( r_out, j-2.0);
			diff_sum9 += (I*w*a9[3]/pow(r_out, j-2.0) - (1.0-2.0/r_out)*(j-2.0)*a9[3]/pow(r_out,j-1.0));
			sum10 += a10[3] / pow( r_out, j-2.0);
			diff_sum10 += (I*w*a10[3]/pow(r_out, j-2.0) - (1.0-2.0/r_out)*(j-2.0)*a10[3]/pow(r_out,j-1.0));

			delta9 = cabs( a9[3] / pow( r_out, k) ); 
			delta10 = cabs( a9[3] / pow( r_out, k) ); 


			if ( delta9 > 1e-14 || delta10 > 1e-14)
			{ j_max=j_max+1;}

			a9[0] = a9[1];  
			a9[1] = a9[2];
			a9[2] = a9[3];

			a10[0] = a10[1];  
			a10[1] = a10[2];
			a10[2] = a10[3];

	 
		 }
		//printf("k_max =%d  sum9=%g + I*%g , sum10=%g+ I*%g, Dsum9=%g+I*%g,  Dsum10=%g + I*%g\n",j_max-3, creal(sum9), cimag(sum9), creal(sum10), cimag(sum10),creal(diff_sum9), cimag(diff_sum9),creal(diff_sum10), cimag(diff_sum10));

		if(w == 0) break;

		if(!cisnan(sum9) && !cisnan(sum10) && !cisnan(diff_sum9) && !cisnan(diff_sum10)  ) {
			bcs->R[0] 		= sum9 * cexp(I*omega*rs_out);
			bcs->R_rs_deriv[0] 	= diff_sum9 * cexp(I*omega*rs_out);
			bcs->R[1] 		= sum10 * cexp(I*omega*rs_out);
			bcs->R_rs_deriv[1] 	= diff_sum10 * cexp(I*omega*rs_out);
			converged = 1;
		}else{
			rs_out += 5000.0;		//printf("BOUNDARY INCREASES\n\n");
		}

	}

	bcs->rs_boundary 	= rs_out; 

}


void odd_in_bcs(double* a0s, struct boundary_conditions *bcs)
{

	double rs_in;
	//if(bcs->omega != 0){
	//	rs_in = -10.0/fabs(bcs->omega);
	//}else{
		rs_in = -50.0;
	//}

	double r_in = r_of_rs(rs_in);

	double l = bcs->l;
	double w = bcs->omega;
	double omega = w;

	    double epsilon = r_in - 2.0;    // expansion term for inner BC sum
	    int j, k;    // j = k+3 where k is the actual index used in analytic work.
	    int j_max = 5;
	    double L = l*(l+1.0);
	    double lambda = (l+2.0)*(l-1.0);

	    double complex b9[4], b10[4];  
	    double delta9, delta10; 
	    double complex C, D, E, F, G, H, J, N, P;
	   
	   
	    b9[0] = 0.0;    b10[0] = 0.0;
	    b9[1] = 0.0;     b10[1] = 0.0;


	    b9[2] = a0s[0];  // this b9_0 = 1.0 for k=0
	    b10[2] = a0s[1];

	double complex sum9 = b9[2];
	double complex sum10 = b10[2];
	double complex diff_sum9 = -1.0*I*w*b9[2];
	double complex diff_sum10 = -1.0*I*w*b10[2];
   
     	for ( j = 3; j < j_max; j++ )
     	{
		k = j - 2;
		C = 2.0*((k-1)+ 12.0*I*w*(k-1)-2.0*(k-1)*(k-1) + L - 4.0);
		D = (4.0 + 12*I*w*(k-2) + L -(k-2)*(k-3.0));
		E = 2.0*I*w*(k-3);
		F = 2.0;
		G = -2.0;
		H = 2*(k-1.0+ 12.0*I*w*(k-1)-2.0*(k-1)*(k-1) + L - 1.0);
		J = (-2.0 + 12*I*w*(k-2) + L -(k-2)*(k-3.0));
		N = -4.0*lambda;
		P = -2.0*lambda;
	

		b9[3] = ( C*b9[2]+D*b9[1]+E*b9[0]+F*b10[2]+G*b10[1] ) / ( k*4.0*(k-4.0*I*w) ); // recursion relation for the b9_k's
		b10[3] = ( H*b10[2]+J*b10[1]+E*b10[0]+N*b9[2]+P*b9[1] ) / ( 4.0*k*(k-4.0*I*w) );
	
		sum9 +=  b9[3]*pow( epsilon, j-2.0);
		sum10 +=  + b10[3]*pow( epsilon, j-2.0);
		diff_sum9 +=  - 1.0*I*w*b9[3]*pow( epsilon, j-2.0) + (j-2.0)*b9[3]*pow( epsilon, j-2.0)/ r_in;
		diff_sum10 +=  - 1.0*I*w*b10[3]*pow( epsilon, j-2.0) + (j-2.0)*b10[3]*pow( epsilon, j-2.0)/ r_in;

		delta9 = cabs( b9[3]*pow( epsilon, j-2.0) );   
		delta10 = cabs( b10[3]*pow( epsilon, j-2.0) );

	 	if ( delta9 > 1e-14 || delta10 > 1e-14)
		{ j_max=j_max+1;}

		b9[0] = b9[1];  
		b9[1] = b9[2];
		b9[2] = b9[3];

		b10[0] = b10[1];  
		b10[1] = b10[2];
		b10[2] = b10[3];
   
     }
   
    bcs->R[0] 		= sum9 * cexp(-I*omega*rs_in);
    bcs->R_rs_deriv[0] 	= diff_sum9 * cexp(-I*omega*rs_in);
    bcs->R[1] 		= sum10 * cexp(-I*omega*rs_in);
    bcs->R_rs_deriv[1] 	= diff_sum10 * cexp(-I*omega*rs_in);


	bcs->rs_boundary 	= rs_in; 
	

}


void odd_dipole_out_bcs(double* a0s, struct boundary_conditions *bcs)
{
	double l = bcs->l;
	double w = bcs->omega;
	double omega = w;

	double rs_out = 110000.0;

	int converged = 0;

	while(!converged){


		double r_out = r_of_rs(rs_out);

		int j, k;    // j = k+3 where k is the actual index used in analytic work.
		int j_max = 5;
		double L = l*(l+1.0);
		double lambda = (l+2.0)*(l-1.0);

	     	double complex a9[4];
		double complex a10[4];  
		double delta9, delta10; 

		double complex C, D, E, Ii, J, K;
		double F, G, H, N, P;

		F = 2.0;
		G = -10.0;
		H = 12.0;
		N = 2.0*lambda;
		P = -4.0*lambda;

		a9[0] = 0.0;    a10[0] = 0.0;
		a9[1] = 0.0;     a10[1] = 0.0;

		a9[2] = a0s[0];  // this a9_0 = 1.0 for k=0
		a10[2] = 0;
 	
		double complex sum9 = a9[2];
		double complex sum10 = a10[2];
		double complex diff_sum9 = I*w*a9[2];
		double complex diff_sum10 = I*w*a10[2];


		for ( j = 3; j <= j_max; j++ )
		{
			k = j - 2;
			C = ( 4.0*I*w*(k-1) + (k-1)*k - L - 4.0 );
			D = -6.0*(k-2)-4.0*(k-2)*(k-2)+24.0+2.0*L;
			E = 4.0*( (k-3)*(k-3)+2.0*(k-3)-8.0 );
			Ii = 4.0*I*w*(k-1) + (k-1)*k - L + 2.0;
			J = -6.0*(k-2) - 4.0*(k-2)*(k-2) - 6.0 + 2.0*L;
			K = 4.0*( (k-3)*(k-3) + 2.0*(k-3) + 1.0 );

			a9[3] = ( C*a9[2]+D*a9[1]+E*a9[0]+F*a10[2]+G*a10[1] + H*a10[0] ) / (2.0*I*w*k ); // recursion relation for the b9_k's
			a10[3] = ( Ii*a10[2]+J*a10[1]+K*a10[0]+N*a9[2]+P*a9[1] ) / (2.0*I*w*k );

			sum9 += a9[3] / pow( r_out, j-2.0);
			diff_sum9 += (I*w*a9[3]/pow(r_out, j-2.0) - (1.0-2.0/r_out)*(j-2.0)*a9[3]/pow(r_out,j-1.0));
			sum10 += a10[3] / pow( r_out, j-2.0);
			diff_sum10 += (I*w*a10[3]/pow(r_out, j-2.0) - (1.0-2.0/r_out)*(j-2.0)*a10[3]/pow(r_out,j-1.0));

			delta9 = cabs( a9[3] / pow( r_out, k) ); 
			delta10 = cabs( a9[3] / pow( r_out, k) ); 


			if ( delta9 > 1e-14 || delta10 > 1e-14)
			{ j_max=j_max+1;}

			a9[0] = a9[1];  
			a9[1] = a9[2];
			a9[2] = a9[3];

			a10[0] = a10[1];  
			a10[1] = a10[2];
			a10[2] = a10[3];

	 
		 }

		if(w == 0) break;

		if(!cisnan(sum9) && !cisnan(diff_sum9)){
			bcs->R[0] 		= sum9 * cexp(I*omega*rs_out);;
			bcs->R_rs_deriv[0] 	= diff_sum9 * cexp(I*omega*rs_out);;
			converged = 1;
		}else{
			rs_out += 5000.0;		//printf("BOUNDARY INCREASES\n\n");
		}

	}

	bcs->rs_boundary 	= rs_out; 

}

void odd_dipole_in_bcs(double* a0s, struct boundary_conditions *bcs)
{
	double rs_in;
//	if(bcs->omega != 0){
	//	rs_in = -10.0/fabs(bcs->omega);
	//}else{
		rs_in = -50.0;
	//}
	double r_in = r_of_rs(rs_in);


	double l = bcs->l;
	double w = bcs->omega;
	double omega = w;

    double epsilon = r_in - 2.0;    // expansion term for inner BC sum
    int j, k;    // j = k+3 where k is the actual index used in analytic work.
    int j_max = 5;
    double L = l*(l+1.0);
    double lambda = (l+2.0)*(l-1.0);

	    double complex b9[4], b10[4];  
	    double delta9, delta10; 
	    double complex C, D, E, F, G, H, J, N, P;
	   
	   
	b9[0] = 0.0;    b10[0] = 0.0;
	b9[1] = 0.0;     b10[1] = 0.0;

	b9[2] = a0s[0];  // this b9_0 = 1.0 for k=0
	b10[2] = 0.0;

	double complex sum9 = b9[2];
	double complex sum10 = b10[2];
	double complex diff_sum9 = -1.0*I*w*b9[2];
	double complex diff_sum10 = -1.0*I*w*b10[2];
   
     	for ( j = 3; j < j_max; j++ )
     	{
		k = j - 2;
		C = 2.0*((k-1)+ 12.0*I*w*(k-1)-2.0*(k-1)*(k-1) + L - 4.0);
		D = (4.0 + 12*I*w*(k-2) + L -(k-2)*(k-3.0));
		E = 2.0*I*w*(k-3);
		F = 2.0;
		G = -2.0;
		H = 2*(k-1.0+ 12.0*I*w*(k-1)-2.0*(k-1)*(k-1) + L - 1.0);
		J = (-2.0 + 12*I*w*(k-2) + L -(k-2)*(k-3.0));
		N = -4.0*lambda;
		P = -2.0*lambda;
	

		b9[3] = ( C*b9[2]+D*b9[1]+E*b9[0]+F*b10[2]+G*b10[1] ) / ( k*4.0*(k-4.0*I*w) ); // recursion relation for the b9_k's
		b10[3] = ( H*b10[2]+J*b10[1]+E*b10[0]+N*b9[2]+P*b9[1] ) / ( 4.0*k*(k-4.0*I*w) );
	
		sum9 +=  b9[3]*pow( epsilon, j-2.0);
		sum10 +=  + b10[3]*pow( epsilon, j-2.0);
		diff_sum9 +=  - 1.0*I*w*b9[3]*pow( epsilon, j-2.0) + (j-2.0)*b9[3]*pow( epsilon, j-2.0)/ r_in;
		diff_sum10 +=  - 1.0*I*w*b10[3]*pow( epsilon, j-2.0) + (j-2.0)*b10[3]*pow( epsilon, j-2.0)/ r_in;

		delta9 = cabs( b9[3]*pow( epsilon, j-2.0) );   
		delta10 = cabs( b10[3]*pow( epsilon, j-2.0) );

	 	if ( delta9 > 1e-14 || delta10 > 1e-14)
		{ j_max=j_max+1;}

		b9[0] = b9[1];  
		b9[1] = b9[2];
		b9[2] = b9[3];

		b10[0] = b10[1];  
		b10[1] = b10[2];
		b10[2] = b10[3];
   
     }

   
    	bcs->R[0] 		= sum9 * cexp(-I*omega*rs_in);
    	bcs->R_rs_deriv[0] 	= diff_sum9 * cexp(-I*omega*rs_in);


	bcs->rs_boundary 	= rs_in; 

}

void even_static_out_bcs(double* a0s, struct boundary_conditions *bcs)
{
	double rs_out = 110000.0;//r_max() + 200.0;
	double r_out = r_of_rs(rs_out);

	int l = bcs->l;

	double sum1, sum3, sum5, diff_sum1, diff_sum3, diff_sum5;

	sum1 = 0.0;		sum3 = 0.0;		sum5 = 0.0;
	diff_sum1 = 0.0;        diff_sum3 = 0.0;        diff_sum5 = 0.0;

	int k;
	double M = 1.0;

	double a1[4], a3[4], a5[4], a1l[4], a3l[4], a5l[4];
	double c1, c2, c3;

	c1 = a0s[0];
	c2 = a0s[1];
	c3 = a0s[2];


	// The analytic expression below for the coefficients a1,3,5[1,2] and a1l,3l,5l[1,2] are derived using the recursion relations in Mathematica
	// See Mathematica files static_even_outer_BC_Jan11_2012.nb and static_even_outer_BC_el_plus_1_case_Jan9_2012.nb
	//*********************************************************

		a3[0]  = c1;
		a5[0]  = c2;
		a1[0]  = a3[0] + a5[0]/(l+1.0);

		a1l[0] = 0.0;
		a3l[0] = 0.0;
		a5l[0] = 0.0;

		a1[1] = ((-1 + l)*(c2 + c1*l))/l;
		a3[1] = c1 + c2/l + c1*l;
		a5[1] = (2*c1*l*(1 + l) + c2*(-2 - Power(l,2) + Power(l,3)))/(l*(1 + l));

		a1l[1] = 0.0;
		a3l[1] = 0.0;
		a5l[1] = 0.0;

		a1[2] = -c3/(2.*l) + (2*c1*l*(-6 - 2*l + 12*Power(l,2) + 9*Power(l,3) +  4*Power(l,4) + 5*Power(l,5) + 2*Power(l,6)) + c2*(12 + 8*l + 6*Power(l,2) - 15*Power(l,3) - 8*Power(l,4) + 15*Power(l,5) + 6*Power(l,6)))/ (2.*Power(l,2)*(1 + l)*(1 + 2*l)*(3 + 2*l));
		a3[2] = c3/(2.*l) + (c2*(-12 - 8*l + 6*Power(l,2) + 51*Power(l,3) + 42*Power(l,4) + 7*Power(l,5) - 2*Power(l,6)) + 2*c1*l*(6 + 6*l + 6*Power(l,2) + 23*Power(l,3) + 28*Power(l,4) + 13*Power(l,5) + 2*Power(l,6)))/ (2.*Power(l,2)*(1 + l)*(1 + 2*l)*(3 + 2*l));
		a5[2] = c3;

		a1l[2] = (-2*(2 + l)*(-(c2*(-3 + l)) + 2*c1*(1 + l)))/((1 + l)*(1 + 2*l)*(3 + 2*l));
		a3l[2] = -1.0*a1l[2];
		a5l[2] = 2.0*l*a3l[2];


		for( k=0; k<=2; k++ )
		{

			sum1 +=	 ( a1[k] + a1l[k]*log(r_out) )*pow( r_out, -(l+k) );
			sum3 +=  ( a3[k] + a3l[k]*log(r_out) )*pow( r_out, -(l+k) );
			sum5 +=  ( a5[k] + a5l[k]*log(r_out) )*pow( r_out, -(l+k) );

			diff_sum1 += f(r_out)*( -(l+k)*a1[k] + ( 1.0-(l+k)*log(r_out) )*a1l[k] )/pow(r_out, l+k+1.0);
			diff_sum3 += f(r_out)*( -(l+k)*a3[k] + ( 1.0-(l+k)*log(r_out) )*a3l[k] )/pow(r_out, l+k+1.0);
			diff_sum5 += f(r_out)*( -(l+k)*a5[k] + ( 1.0-(l+k)*log(r_out) )*a5l[k] )/pow(r_out, l+k+1.0);
		}


		for( k=l+3; k<=l+10; k++ )
		{

			double C1 = ( l*(l+1) + 1.0 - k*k );
			double C5 = ( l*(l+1) + k*(1.0-k) );
			double D5 = 2*k - 1.0;

			int j = k - 1;
			double D1 = j*(j - 1.0);
			double D3 = 2*(j + 1.0);
			double E1 = 1 - 2.0*j;
			double G3 = ( l*(l+1) + 2.0 - 2*j*j );
			double G1 = -2.0*j;
			double H3 = 4.0*j;
			double E5 = ( j*(1-j) + 2.0 );
			double D52 = 2*j - 1.0;

			int n = k - 2;
			double F3 = n + 1.0;
			double I3 = n*n - 1.0;
			double J3 = -2.0*n;


			double L = l*(l + 1.0);
			
			double Mat[] = { 	C1, 2.0*k, -(k + 1.0), 1, -1, 0, 
						0, C1, 0, -(k + 1.0), 0, -1, 
						-(k + 1.0), 1, C1, 2.0*k, 1, 0, 
						0, -(k + 1.0), 0, C1, 0, 1, 
						-2*L, 0, 2*L, 0, C5, D5, 
						0, -2*L, 0, 2*L, 0, C5	
					};

								
			double A = -2*M*( D1*a1[2] + D3*a3[2] + a5[2] + E1*a1l[2] - 2*a3l[2] ) + 
					4.0*M*M*( F3*a3[1] - a3l[1] );
			double B = -2*M*( D1*a1l[2] + D3*a3l[2] + a5l[2] ) + 
			  		 4*M*M*F3*a3l[1];
			double C = 2*M*( G3*a3[2] + G1*a1[2] + H3*a3l[2] + 2*a1l[2] ) + 
					   4*M*M*( I3*a3[1] + J3*a3l[1] );
			double D = 2*M*( G3*a3l[2] + G1*a1l[2] ) + 
					   4*M*M*I3*a3l[1];
			double E = 2*M*( E5*a5[2] + 2*L*a3[2] + D52*a5l[2] );
			double Ff = 2*M*( E5*a5l[2] + 2*L*a3l[2] );

			double Vec[] = { A, B, C, D, E, Ff };

		       gsl_matrix_view m 
			 = gsl_matrix_view_array (Mat, 6, 6);
	     
		       gsl_vector_view b
			 = gsl_vector_view_array (Vec, 6);
		     
		       gsl_vector *x = gsl_vector_alloc (6);
			       
		       int s;
		     
		       gsl_permutation * p = gsl_permutation_alloc (6);
		       gsl_linalg_LU_decomp (&m.matrix, p, &s);
		       gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

			int i;
			double row[6]; 

			for (i = 0; i < 6; i++)
			{
				row[i]= gsl_vector_get(x,i);
				//printf("%.14f\t", row[i] );
			}
 
			a1l[3] = row[1];
			a3l[3] = row[3];
			a5l[3] = row[5];

			a1[3] = row[0];
			a3[3] = row[2];
			a5[3] = row[4];
		     
		       gsl_permutation_free (p);
		       gsl_vector_free (x);
				
			//printf("k = %d:\ta1=%f\t\ta3=%f\t\ta5=%f\t\t", k, a1[3], a3[3], a5[3]);
			//printf("a1l=%f\t  a3l=%f\t  a5l=%f\n",a1l[3],a3l[3],a5l[3]);
		
			sum1 += ( a1[3] + a1l[3]*log(r_out) )/pow(r_out,k); 
			sum3 += ( a3[3] + a3l[3]*log(r_out) )/pow(r_out,k);  
			sum5 += ( a5[3] + a5l[3]*log(r_out) )/pow(r_out,k); 
			
			diff_sum1 += f(r_out)*( -k*a1[3]+(1.0-k*log(r_out))*a1l[3] )/pow(r_out,k+1.0);
			diff_sum3 += f(r_out)*( -k*a3[3]+(1.0-k*log(r_out))*a3l[3] )/pow(r_out,k+1.0);
			diff_sum5 += f(r_out)*( -k*a5[3]+(1.0-k*log(r_out))*a5l[3] )/pow(r_out,k+1.0); 


			a1[0] = a1[1];  
			a1[1] = a1[2];
			a1[2] = a1[3];

			a3[0] = a3[1];  
			a3[1] = a3[2];
			a3[2] = a3[3];

			a5[0] = a5[1];  
			a5[1] = a5[2];
			a5[2] = a5[3];

			a1l[0] = a1l[1];  
			a1l[1] = a1l[2];
			a1l[2] = a1l[3];

			a3l[0] = a3l[1];  
			a3l[1] = a3l[2];
			a3l[2] = a3l[3];

			a5l[0] = a5l[1];  
			a5l[1] = a5l[2];
			a5l[2] = a5l[3];


		} // end k-loop
	
		//***** Beyond l = 18, the BCs are smaller than ~ 1e-40, which causes the ODE solver to crap out
		//      I had determined that its tolerance for delta_abs is about 1e-42
		//   	So for l>=18 modes, we multiply the BCs by r_out^l
		//	There is no need to divide by r_out^l later because the way the scaling coefficients are computed
		//	is scale invariant.

		if ( l >=18 )
		{
			sum1 	= sum1*pow( r_out, l);
			sum3 	= sum3*pow( r_out, l);
			sum5 	= sum5*pow( r_out, l);
			diff_sum1 = diff_sum1*pow( r_out, l);
			diff_sum3 = diff_sum3*pow( r_out, l);
			diff_sum5 = diff_sum5*pow( r_out, l);

		}
	



	
	// No exponential as omega is zero
	bcs->R[0] = sum1;
	bcs->R[1] = sum3;
	bcs->R[2] = sum5;

	bcs->R_rs_deriv[0] = diff_sum1;
	bcs->R_rs_deriv[1] = diff_sum3;
	bcs->R_rs_deriv[2] = diff_sum5;

	bcs->rs_boundary 	= rs_out; 
}

void even_static_in_bcs(double* a0s, struct boundary_conditions *bcs)
{
	double rs_in = -30.0;
	double r_in = r_of_rs(rs_in);

	double l = bcs->l;

	double sum1, sum3, sum5, diff_sum1, diff_sum3, diff_sum5;
	double eps = r_in -2.0;

	int k, k_max=6;
	double L = l*(l+1.0);
	double b1[4], b3[4], b5[4];
	double delta1, delta3, delta5;
	double C5, D5, C3, D3,D1, E3, E1, C1, F1, G1, G3;

	double c1, c2, c3;

	c1 = a0s[0];
	c2 = a0s[1];
	c3 = a0s[2];

	//printf("c1=%g c2=%g c3=%g\n",c1,c2,c3);

	b1[0] = 0.0; b1[1] = 0.0; b5[0] = 0.0;
	b3[0] = c1;  b3[1] = c2;  b5[1] = c3;   // free parameters
	b5[2] = 0.25*(-4.0*b5[1] + l*(b3[0]+b5[1]) + l*l*(b3[0]+b5[1]) ); 
	b1[2] = ( l*(l+1.0)*b3[0]+ 2.0*b5[1])/8.0;
	

	sum3 = c1 + c2*eps; 
	diff_sum3 = b3[1]*pow(eps, 1.0)/r_in; 
	sum1 = 0.0; diff_sum1 = 0.0;
	sum5 = 0.0; diff_sum5 = 0.0;

	for ( k = 0; k <= 2; k++)
	{
		sum1 += b1[k]*pow(eps, k); 
		sum5 += b5[k]*pow(eps, k);
		diff_sum1 += k*b1[k]*pow(eps,k)/r_in;
		diff_sum5 += k*b5[k]*pow(eps,k)/r_in;

	}

	for ( k=3 ; k <= k_max; k++)
	{
		C5=2.0*L-4.0*((k-1.0)*(k-1.0)+1.0); D5=L-(k-2.0)*(k-1.0); C3=4.0*(k-1.0)*(2.0-k); C1=-4.0*k; 
		D3=2.0*(L+(k-2.0)*(1.0-2.0*(k-2.0))); D1=2.0*(k-4.0); E3=L+1.0-(k-3.0)*(k-3.0); E1=k-4.0; 
		F1=4.0*(L+1.0+(k-1.0)*(4.0-3.0*(k-1.0))); G1=2.0*(2.0*L+2.0+(k-2.0)*(2.0-3.0*(k-2.0))); G3=2.0*(k-2.0);

		b5[3]=(C5*b5[2]-4.0*L*b1[2]+D5*b5[1]+2.0*L*(b3[1]-b1[1]))/(4.0*k*(k-1.0));
		b1[3]=(F1*b1[2]+G1*b1[1]+G3*b3[1]-2.0*b5[1]+E3*b1[0]+E1*b3[0]-b5[0])/(8.0*k*(k-2.0));
		b3[2]=(8.0*k*b1[3]-C1*b1[2]-4.0*b5[2]-D3*b3[1]-D1*b1[1]-4.0*b5[1]-E3*b3[0]-E1*b1[0]-b5[0] )/C3;	
		//printf("k=%d, b1_k=%.14f  b3_k=%.14f  b5_k=%.14f\n",k, b1[3], b3[2],b5[3]);

		delta1 = fabs( b1[3]*pow( eps, k) ); delta3 = fabs( b3[2]*pow( eps,k-1) ); delta5 = fabs( b5[3]*pow( eps, k) ); 
	
      		 if ( delta1 >= 1e-14 || delta3 >= 1e-14 || delta5 >= 1e-14  )
        	{ k_max = k_max+1; }

		sum1 = sum1 + b1[3]*pow(eps, k); 
		sum3 = sum3 + b3[2]*pow(eps, k-1); 
		sum5 = sum5 + b5[3]*pow(eps, k);
				
		diff_sum1 = diff_sum1 + k*b1[3]*pow(eps,k)/r_in;
		diff_sum3 = diff_sum3 + (k-1)*b3[2]*pow(eps,k-1)/r_in;
		diff_sum5 = diff_sum5 + k*b5[3]*pow(eps,k)/r_in;

		b1[0] = b1[1];  
		b1[1] = b1[2];
		b1[2] = b1[3];

		b3[0] = b3[1];  
		b3[1] = b3[2];

		b5[0] = b5[1];  
		b5[1] = b5[2];
		b5[2] = b5[3];

	}

	// No exponential as omega is zero
	bcs->R[0] = sum1;
	bcs->R[1] = sum3;
	bcs->R[2] = sum5;

	bcs->R_rs_deriv[0] = diff_sum1;
	bcs->R_rs_deriv[1] = diff_sum3;
	bcs->R_rs_deriv[2] = diff_sum5;

	bcs->rs_boundary 	= rs_in; 
	
	//printf("BC test: %.16e %.16e %.16e %.16e\n", r_in, creal(bcs->R[0]), creal(bcs->R[1]), creal(bcs->R[2]));
	//exit(0);
}



