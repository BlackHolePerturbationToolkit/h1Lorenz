#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sf_gamma.h>
//#include "SF_data.h"
//#include "libSch.h"
#include "math_funcs.h"
#include "reg_params.h"
#include "file_output.h"

#define Power(x, y)    				 				(pow((double)(x), (double)(y)))
#define Sqrt(x)        				 				(sqrt((double)(x)))
#define EllipticK(x)									(elliptic_integral_1st((double)(sqrt(x))))
#define EllipticE(x)									(elliptic_integral_2nd((double)(sqrt(x))))
#define Pi														M_PI

int l_max;
int tail_l_max = 0;
int num_chi;
double p;
double e;
struct orbital_params orbit;

int verbose = 1;

const double DESIRED_FUNCTION_ACCURACY 		= 1e-10;		//!< Sets the desired function accuracy for the wrappers around various GSL functions

void calculate_GSF_from_full_fields(int l_max, int num_chi, struct GSF_data *GSF, struct orbital_params *orbit);
void calculate_tail(int l_max, int num_chi, struct GSF_component_data_piece *GSF_comp_piece, int n, int Nmin, int Nmax);
double ultimate_tail(double *data, int lmax, int n, int Nmin, int Nmax);
void calculate_hat_beta(struct GSF_data *GSF, struct orbital_params *orbit);
void	construct_beta_full(int l_max, int num_chi, struct GSF_data *GSF, struct orbital_params *orbit);
void calculate_hat_beta_from_beta(struct GSF_data *GSF, struct orbital_params *orbit);

int main(int argc, char *argv[])
{

	if( argv[1] != NULL && argv[2] != NULL){
			    e = (double)strtod(argv[1], NULL);
			    p = (double)strtod(argv[2], NULL);
		if(e < 0.0 || e >= 1.0){
			printf("\ne must be great than or equal to 0 and less than 1\n\n");
			exit(0);
		}
	}else{
			    printf("\n e and p values required\n\n");
			    exit(0);
	}	

	// Assertain lmax and num_chi
	FILE* l_mode_file;
	char line[1000];
	int line_count = 0;
	char file_location[100];


	int l;
	for(l = 0; l < 100; l++){
		sprintf(file_location, "data/p%lf/e%lf/modes/l%d", p, e, l);
		//	printf("Opening file: %s\n", file_location);
		l_mode_file=fopen(file_location, "r");

		if(l_mode_file != NULL){
			if(l == 0){
			// count the number of lines
				while ( fgets(line, 500, l_mode_file) != NULL){ 
					if(line[0] == '#') continue;
					line_count++;
				}
				fclose(l_mode_file);
			}
		}else break;

	}
	

	l_max = l-1;
	num_chi = line_count;

	if(verbose) printf("l_max = %d, num_chi = %d\n", l_max, num_chi);

	set_primary_orbital_params(p, e);
	get_orbit_params(&orbit, (num_chi-1)/2);

	// Display the orbital parameters
	if(verbose)print_orbital_params(&orbit);

	// Prepare the structures
	struct GSF_data GSF;
	setup_GSF_data_structure(&GSF, l_max, num_chi);

	// Load the data
	for(l = 0; l <= l_max; l++){
		int i = 0;
		sprintf(file_location, "data/p%lf/e%lf/modes/l%d", p, e, l);
		//	printf("Opening file: %s\n", file_location);
		l_mode_file=fopen(file_location, "r");

		if(l_mode_file != NULL){
			// count the number of lines
			while ( fgets(line, 1000, l_mode_file) != NULL){ 
				if(line[0] == '#') continue;
					strtok(line, " ");
					char *str2 = strtok(NULL, " ");
					char *str3 = strtok(NULL, " ");
					char *str4 = strtok(NULL, " ");
					char *str5 = strtok(NULL, " ");
					char *str6 = strtok(NULL, " ");
					char *str7 = strtok(NULL, " ");
					char *str8 = strtok(NULL, " ");
					char *str9 = strtok(NULL, " ");
					char *str10 = strtok(NULL, " ");
					char *str11 = strtok(NULL, " ");
					char *str12 = strtok(NULL, " ");
					char *str13 = strtok(NULL, " ");
					char *str14 = strtok(NULL, " ");
					char *str15 = strtok(NULL, " ");
					char *str16 = strtok(NULL, " ");
					char *str17 = strtok(NULL, " ");
					char *str18 = strtok(NULL, " ");
					char *str19 = strtok(NULL, " ");

					GSF.Fr_in->total->l_full[l][i] = strtod(str2, NULL);
					GSF.Fr_out->total->l_full[l][i] = strtod(str3, NULL);
					GSF.Ft_in->total->l_full[l][i] = strtod(str4, NULL);
					GSF.Ft_out->total->l_full[l][i] = strtod(str5, NULL);
					GSF.huu_in->total->l_full[l][i] = strtod(str6, NULL);
					GSF.huu_out->total->l_full[l][i] = strtod(str7, NULL);


					GSF.hut_in->total->l_full[l][i] = strtod(str8, NULL);
					GSF.hut_out->total->l_full[l][i] = strtod(str9, NULL);

					GSF.huphi_in->total->l_full[l][i] = strtod(str10, NULL);
					GSF.huphi_out->total->l_full[l][i] = strtod(str11, NULL);

					GSF.Grut_in->total->l_full[l][i] = strtod(str12, NULL);
					GSF.Grut_out->total->l_full[l][i] = strtod(str13, NULL);

					GSF.Gruphi_in->total->l_full[l][i] = strtod(str14, NULL);
					GSF.Gruphi_out->total->l_full[l][i] = strtod(str15, NULL);

					GSF.Gphiur_in->total->l_full[l][i] = strtod(str16, NULL);
					GSF.Gphiur_out->total->l_full[l][i] = strtod(str17, NULL);

					GSF.beta_tensor_in->total->l_full[l][i] = strtod(str18, NULL);
					GSF.beta_tensor_out->total->l_full[l][i] = strtod(str19, NULL);


					i++;
			}
		}else break;

	}

	// Construct beta_full from hut, huphi, Grut, Gruphi, Gphiur
	construct_beta_full(l_max, num_chi, &GSF, &orbit);

	// Regularize the data
	calculate_GSF_from_full_fields(l_max, num_chi, &GSF, &orbit);

	// Output the data
	output_reg_data_to_file(l_max, &GSF, &orbit);

	// calculate hat_beta
//	calculate_hat_beta(&GSF, &orbit);
	calculate_hat_beta_from_beta(&GSF, &orbit);

	int k;
	for(k=0;k<=l_max; k++){
	if(verbose)	printf("%d %.12e %.12e %.12e %.12e %.12e\n", k, GSF.beta_tensor_in->total->l_full[k][0], GSF.beta_tensor_out->total->l_full[k][0], GSF.beta_in->total->l_full[k][0], GSF.beta_out->total->l_full[k][0],  GSF.huu_in->total->l_full[k][0]);
	}

	return 0;
}

void calculate_GSF_from_full_fields(int l_max, int num_chi, struct GSF_data *GSF, struct orbital_params *orbit)
{
		int i;

		double E = orbit->E;
		double L = orbit->L;


		// Regularize
		for(i=0; i < num_chi; i++){
			double r_p = orbit->r_p[i];
			double f_p = f(r_p);
			double contra_ur = orbit_contra_ur(orbit->chi_p[i]);
			int l;

			for(l = 0; l <= l_max; l++){

				// Calculate the conservative and dissipative pieces of the GSF
				GSF->Fr_in->cons->l_full[l][i] 		= 0.5 * ( GSF->Fr_in->total->l_full[l][i] + GSF->Fr_in->total->l_full[l][num_chi -1 - i]);
				GSF->Fr_out->cons->l_full[l][i]		= 0.5 * ( GSF->Fr_out->total->l_full[l][i] + GSF->Fr_out->total->l_full[l][num_chi -1 - i]);

				GSF->Ft_in->cons->l_full[l][i]		= 0.5 * ( GSF->Ft_in->total->l_full[l][i] - GSF->Ft_in->total->l_full[l][num_chi -1 - i]);
				GSF->Ft_out->cons->l_full[l][i] 	= 0.5 * ( GSF->Ft_out->total->l_full[l][i] - GSF->Ft_out->total->l_full[l][num_chi -1 - i]);

				GSF->Fr_in->l_full_diss[l][i] 		= 0.5 * ( GSF->Fr_in->total->l_full[l][i] - GSF->Fr_in->total->l_full[l][num_chi -1 - i]);
				GSF->Fr_out->l_full_diss[l][i]		= 0.5 * ( GSF->Fr_out->total->l_full[l][i] - GSF->Fr_out->total->l_full[l][num_chi -1 - i]);

				GSF->Ft_in->l_full_diss[l][i]			= 0.5 * ( GSF->Ft_in->total->l_full[l][i] + GSF->Ft_in->total->l_full[l][num_chi -1 - i]);
				GSF->Ft_out->l_full_diss[l][i] 		= 0.5 * ( GSF->Ft_out->total->l_full[l][i] + GSF->Ft_out->total->l_full[l][num_chi -1 - i]);

				// Calculate the regularization parameters
				double U 	= 1.0 + pow(L,2.0)/pow(r_p,2.0);
				double w	= pow(L, 2.0)/(r_p*r_p + pow(L, 2.0));

				double Khat	= elliptic_integral_1st(sqrt(w));
				double Ehat	= elliptic_integral_2nd(sqrt(w));

				double Ar	= E/(pow(r_p, 2.0)*U);
				double Br	= -1.0/(M_PI*r_p*r_p*pow(U,1.5)) * ((E*E + f_p*U)*Khat - (2.0*E*E*(1.0-U) - f_p*U*(1.0-2.0*U))*Ehat );
				double Drl	= Anna_and_Barrys_D_r(l, orbit->p, orbit->e, r_p, orbit->L);
				double Frl  = Anna_and_Barrys_F_r(l, orbit->p, orbit->e, r_p, contra_ur, orbit->L);

				double At	= contra_ur/(r_p*r_p*f_p*U);
				double Bt	= E*contra_ur/(M_PI*r_p*r_p*f_p*pow(U,1.5))*(-Khat + 2.0*(1.0-U)*Ehat );
				double Dtl  = Anna_and_Barrys_D_t(l, orbit->p, orbit->e, r_p, contra_ur, orbit->L);
				double Ftl  = Anna_and_Barrys_F_t(l, orbit->p, orbit->e, r_p, contra_ur, orbit->L);

				// The metric perturbation regularization parameters
				double BH 	= 2.0/(M_PI*sqrt(r_p*r_p + L*L)) * Khat;
				double DH	= Anna_and_Barrys_D_H(l, r_p, orbit->E, orbit->L);
				double FH	= Anna_and_Barrys_F_H(l, r_p, orbit->E, orbit->L);

				// The metric perturbation required for spin prececession
				double Bhut 	= 4.0 * E * Khat/(M_PI * sqrt(L*L + r_p*r_p));
				double Bhuphi = - 4.0 * L * Khat/(M_PI * sqrt(L*L + r_p*r_p));

				double AGrutl	= (r_p-2.0)*(2.0*l+1.0)/(2.0*pow(r_p,3.0)); 
				double BGrutl = E*((L*L*(8*M - 4.0*r_p) + r_p*r_p*(2.0*M-r_p+4.0*E*E*r_p))*Ehat + 2.0*(2.0*M-r_p)*r_p*r_p*Khat)/(M_PI*pow(r_p,4.0)*sqrt(L*L+r_p*r_p));

				double AGruphil = -E*L*(2.0*l+1.0)/(2.0*(L*L+r_p*r_p));
				double BGruphil = ((-(Power(L,6)*(8 - 4*r_p)) - Power(r_p,6)*(-2 + r_p + 3*Power(E,2)*r_p) - 2*Power(L,4)*Power(r_p,2)*(6 + (-3 + 2*Power(E,2))*r_p) - Power(L,2)*Power(r_p,4)*(2 + (-1 + 5*Power(E,2))*r_p))*Ehat + Power(r_p,4)*(Power(L,2)*(-2 + r_p + 2*Power(E,2)*r_p) + Power(r_p,2)*(-2 + r_p + 3*Power(E,2)*r_p))*Khat)/(L*Pi*Power(r_p,4)*Power(Power(L,2) + Power(r_p,2),1.5));

				double AGphiurl = E*L*(2.0*l+1.0)/(2*r_p*(r_p-2.0)*(L*L+r_p*r_p));
				double BGphiurl = ((Power(L,6)*(-8 + 4*r_p) + Power(r_p,6)*(-2 + r_p + 3*Power(E,2)*r_p) + Power(L,2)*Power(r_p,4)*(-14 + (7 - 3*Power(E,2))*r_p) -         2*Power(L,4)*Power(r_p,2)*(10 + (-5 + 2*Power(E,2))*r_p))*Ehat - Power(r_p,4)*(Power(L,2)*(-2 + r_p + 2*Power(E,2)*r_p) + Power(r_p,2)*(-2 + r_p + 3*Power(E,2)*r_p))*Khat)/(L*Pi*(-2 + r_p)*Power(r_p,5)*Power(Power(L,2) + Power(r_p,2),1.5));

				double Abeta = ((1 + 2*l)*(3*M - r_p))/(4.*(2*M - r_p)*r_p);
				double Bbeta = ((3*M - r_p)*((17*M - 8*r_p)*EllipticE(M/(-2*M + r_p)) + 2*(-9*M + 4*r_p)*EllipticK(M/(-2*M + r_p))))/(2.*M*Pi*Power(r_p,1.5)*Sqrt(-2*M + r_p));

				if(i == 0 && l == 0) printf("A = %.12e, B = %.12e\n", Abeta/(2.0*l+1.0), Bbeta);

				double Dbeta = -((4652*Power(M,4) + 16*Power(r_p,3)*(2 + 3*r_p) - M*Power(r_p,2)*(192 + 727*r_p) + 4*Power(M,2)*r_p*(88 + 863*r_p) - 
        3*Power(M,3)*(64 + 2235*r_p))*EllipticE(M/(-2*M + r_p)) - 
     (3*M - r_p)*(2021*Power(M,3) - 16*Power(r_p,2)*(2 + 3*r_p) - 8*Power(M,2)*(12 + 253*r_p) + M*r_p*(128 + 607*r_p))*EllipticK(M/(-2*M + r_p)))/
   (4.*(-1 + 2*l)*(3 + 2*l)*M*Pi*(3*M - r_p)*Power(r_p,2.5)*Sqrt(-2*M + r_p));

				double Fbeta = -(-3*((10407790*Power(M,6) + 2*Power(r_p,5)*(2080 + 12747*r_p) - M*Power(r_p,4)*(52480 + 422477*r_p) + 2*Power(M,2)*Power(r_p,3)*(153280 + 1460673*r_p) - 
          5*Power(M,5)*(133760 + 4713413*r_p) - 2*Power(M,3)*Power(r_p,2)*(455360 + 5359257*r_p) + 2*Power(M,4)*r_p*(641440 + 10952593*r_p))*EllipticE(M/(-2*M + r_p))\
        - 2*(3*M - r_p)*(2132870*Power(M,5) - Power(r_p,4)*(1600 + 12507*r_p) - 80*Power(M,4)*(1336 + 48531*r_p) + M*Power(r_p,3)*(22720 + 175621*r_p) + 
          5*Power(M,3)*r_p*(36736 + 557813*r_p) - Power(M,2)*Power(r_p,2)*(103680 + 991649*r_p))*EllipticK(M/(-2*M + r_p))))/
   (80.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(3*M - r_p,3)*Power(r_p,3.5)*Sqrt(-2*M + r_p));

				double Dhut = (2*Sqrt(1 - (2*M)/r_p)*((-50*Power(M,2) + 23*M*r_p + Power(r_p,2))*EllipticE(M/(-2*M + r_p)) + 
       (72*Power(M,2) - 29*M*r_p - Power(r_p,2))*EllipticK(M/(-2*M + r_p))))/((-3 + 4*l*(1 + l))*Pi*Power(r_p,2)*(-3*M + r_p));

				double Dhuphi = (2*(2*M - r_p)*(49*Power(M,2) - 31*M*r_p + 8*Power(r_p,2))*EllipticE(M/(-2*M + r_p)) + 
     2*(-144*Power(M,3) + 149*Power(M,2)*r_p - 55*M*Power(r_p,2) + 8*Power(r_p,3))*EllipticK(M/(-2*M + r_p)))/
   ((-3 + 4*l*(1 + l))*Pi*r_p*(-3*M + r_p)*Sqrt(M*(-2*M + r_p)));

			double DGruphi = ((2*M - r_p)*((2652*Power(M,4) - 3593*Power(M,3)*r_p + 1724*Power(M,2)*Power(r_p,2) - 343*M*Power(r_p,3) + 24*Power(r_p,4))*
        EllipticE(M/(-2*M + r_p)) - (3*M - r_p)*(957*Power(M,3) - 924*Power(M,2)*r_p + 283*M*Power(r_p,2) - 24*Power(r_p,3))*EllipticK(M/(-2*M + r_p)))
     )/(2.*(-3 + 4*l*(1 + l))*Pi*Power(r_p,3)*Power(-3*M + r_p,2)*Sqrt(M*(-2*M + r_p)));

		double DGrut = (Power(-2*M + r_p,1.5)*((-60*Power(M,3) - Power(M,2)*r_p + 12*M*Power(r_p,2) + Power(r_p,3))*EllipticE(M/(-2*M + r_p)) - 
       (3*M - r_p)*(81*Power(M,2) + 4*M*r_p - Power(r_p,2))*EllipticK(M/(-2*M + r_p))))/(2.*(-1 + 2*l)*(3 + 2*l)*Pi*Power(r_p,4.5)*Power(-3*M + r_p,2));

		double DGphiur = ((2844*Power(M,4) - 4009*Power(M,3)*r_p + 2012*Power(M,2)*Power(r_p,2) - 407*M*Power(r_p,3) + 24*Power(r_p,4))*EllipticE(M/(-2*M + r_p)) - 
     (3*M - r_p)*(1437*Power(M,3) - 1308*Power(M,2)*r_p + 347*M*Power(r_p,2) - 24*Power(r_p,3))*EllipticK(M/(-2*M + r_p)))/
   (2.*(-1 + 2*l)*(3 + 2*l)*Pi*Power(r_p,4)*Power(-3*M + r_p,2)*Sqrt(M*(-2*M + r_p)));

			double Fhut = (3*Sqrt(r_p*(-2*M + r_p))*(-2*(2*M - r_p)*(29585*Power(M,5) - 47115*Power(M,4)*r_p + 29578*Power(M,3)*Power(r_p,2) - 
          8470*Power(M,2)*Power(r_p,3) + 949*M*Power(r_p,4) - 15*Power(r_p,5))*EllipticE(M/(-2*M + r_p)) + 
       (3*M - r_p)*(107090*Power(M,5) - 150065*Power(M,4)*r_p + 80015*Power(M,3)*Power(r_p,2) - 19381*Power(M,2)*Power(r_p,3) + 
          1883*M*Power(r_p,4) - 30*Power(r_p,5))*EllipticK(M/(-2*M + r_p))))/
   (10.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r_p,4)*Power(-3*M + r_p,4));

			double Fhuphi = (6*M*(2*M - r_p)*(154985*Power(M,5) - 266755*Power(M,4)*r_p + 170698*Power(M,3)*Power(r_p,2) - 52310*Power(M,2)*Power(r_p,3) + 
        8029*M*Power(r_p,4) - 535*Power(r_p,5))*EllipticE(M/(-2*M + r_p)) - 
     3*M*(3*M - r_p)*(267410*Power(M,5) - 479025*Power(M,4)*r_p + 327375*Power(M,3)*Power(r_p,2) - 105301*Power(M,2)*Power(r_p,3) + 
        15643*M*Power(r_p,4) - 830*Power(r_p,5))*EllipticK(M/(-2*M + r_p)))/
   (10.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r_p,2)*Power(-3*M + r_p,4)*Sqrt(M*(-2*M + r_p)));

			double FGruphi = (Sqrt(M*(-2*M + r_p))*((15845850*Power(M,6) - 40595915*Power(M,5)*r_p + 42696770*Power(M,4)*Power(r_p,2) - 23431254*Power(M,3)*Power(r_p,3) + 
          7019182*Power(M,2)*Power(r_p,4) - 1077127*M*Power(r_p,5) + 65166*Power(r_p,6))*EllipticE(M/(-2*M + r_p)) - 
       2*(3*M - r_p)*(3552570*Power(M,5) - 7094080*Power(M,4)*r_p + 5654325*Power(M,3)*Power(r_p,2) - 2227245*Power(M,2)*Power(r_p,3) + 
          429701*M*Power(r_p,4) - 32223*Power(r_p,5))*EllipticK(M/(-2*M + r_p))))/
   (40.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r_p,4)*Power(-3*M + r_p,4));

			double FGrut = (Power(-2*M + r_p,1.5)*((-178290*Power(M,6) - 300025*Power(M,5)*r_p + 338642*Power(M,4)*Power(r_p,2) + 29614*Power(M,3)*Power(r_p,3) - 
          84254*Power(M,2)*Power(r_p,4) + 16947*M*Power(r_p,5) - 90*Power(r_p,6))*EllipticE(M/(-2*M + r_p)) + 
       2*(3*M - r_p)*(403830*Power(M,5) - 422590*Power(M,4)*r_p + 205765*Power(M,3)*Power(r_p,2) - 61297*Power(M,2)*Power(r_p,3) + 
          8361*M*Power(r_p,4) - 45*Power(r_p,5))*EllipticK(M/(-2*M + r_p))))/
   (40.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r_p,5.5)*Power(-3*M + r_p,4));

			double Fphiur = (M*(-20778330*Power(M,6) + 41696395*Power(M,5)*r_p - 33061570*Power(M,4)*Power(r_p,2) + 13093974*Power(M,3)*Power(r_p,3) - 
        2694062*Power(M,2)*Power(r_p,4) + 272327*M*Power(r_p,5) - 11406*Power(r_p,6))*EllipticE(M/(-2*M + r_p)) + 
     2*M*(3*M - r_p)*(4299450*Power(M,5) - 6849920*Power(M,4)*r_p + 4110165*Power(M,3)*Power(r_p,2) - 1133805*Power(M,2)*Power(r_p,3) + 
        137701*M*Power(r_p,4) - 5343*Power(r_p,5))*EllipticK(M/(-2*M + r_p)))/
   (40.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r_p,5)*Power(-3*M + r_p,4)*Sqrt(M*(-2*M + r_p)));

				//Fbeta = 0;

				// Regularize the various components of the GSF and MP
				GSF->Fr_in->total->l_reg[l][i]		= GSF->Fr_in->total->l_full[l][i] - Ar*(l+0.5) - Br - Drl - Frl;
				GSF->Fr_out->total->l_reg[l][i]		= GSF->Fr_out->total->l_full[l][i] + Ar*(l+0.5) - Br - Drl - Frl;

				GSF->Ft_in->total->l_reg[l][i]		= GSF->Ft_in->total->l_full[l][i] - At*(l+0.5) - Bt - Dtl - Ftl;
				GSF->Ft_out->total->l_reg[l][i]		= GSF->Ft_out->total->l_full[l][i] + At*(l+0.5) - Bt - Dtl - Ftl;

				GSF->huu_in->total->l_reg[l][i]		= GSF->huu_in->total->l_full[l][i] - 2.0*BH - 2.0*DH - 2.0*FH;
				GSF->huu_out->total->l_reg[l][i]	= GSF->huu_out->total->l_full[l][i] - 2.0*BH - 2.0*DH - 2.0*FH;

				GSF->hut_in->total->l_reg[l][i]		= GSF->hut_in->total->l_full[l][i] - Bhut - Dhut - Fhut;
				GSF->hut_out->total->l_reg[l][i]	= GSF->hut_out->total->l_full[l][i] - Bhut - Dhut - Fhut;

				GSF->huphi_in->total->l_reg[l][i]		= GSF->huphi_in->total->l_full[l][i] - Bhuphi - Dhuphi - Fhuphi;
				GSF->huphi_out->total->l_reg[l][i]	= GSF->huphi_out->total->l_full[l][i] - Bhuphi - Dhuphi - Fhuphi;

				GSF->Grut_in->total->l_reg[l][i]		= GSF->Grut_in->total->l_full[l][i] + AGrutl - BGrutl - DGrut - FGrut;
				GSF->Grut_out->total->l_reg[l][i]		= GSF->Grut_out->total->l_full[l][i] - AGrutl - BGrutl - DGrut - FGrut;

				GSF->Gruphi_in->total->l_reg[l][i]		= GSF->Gruphi_in->total->l_full[l][i] + AGruphil - BGruphil - DGruphi - FGruphi;
				GSF->Gruphi_out->total->l_reg[l][i]		= GSF->Gruphi_out->total->l_full[l][i] - AGruphil - BGruphil - DGruphi - FGruphi;

				GSF->Gphiur_in->total->l_reg[l][i]		= GSF->Gphiur_in->total->l_full[l][i] + AGphiurl - BGphiurl - DGphiur - Fphiur;
				GSF->Gphiur_out->total->l_reg[l][i]		= GSF->Gphiur_out->total->l_full[l][i] - AGphiurl - BGphiurl - DGphiur - Fphiur;

				GSF->beta_in->total->l_reg[l][i]		= GSF->beta_in->total->l_full[l][i] + Abeta - Bbeta - Dbeta - Fbeta;
				GSF->beta_out->total->l_reg[l][i]		= GSF->beta_out->total->l_full[l][i] - Abeta - Bbeta - Dbeta - Fbeta;

				// The conservative piece of the GSF is regularized using the same regularization paramters. 
				GSF->Fr_in->cons->l_reg[l][i]		= GSF->Fr_in->cons->l_full[l][i] - Ar*(l+0.5) - Br - Drl - Frl;
				GSF->Fr_out->cons->l_reg[l][i]		= GSF->Fr_out->cons->l_full[l][i] + Ar*(l+0.5) - Br - Drl - Frl;

				GSF->Ft_in->cons->l_reg[l][i]		= GSF->Ft_in->cons->l_full[l][i] - At*(l+0.5) - Bt - Dtl - Ftl;
				GSF->Ft_out->cons->l_reg[l][i]		= GSF->Ft_out->cons->l_full[l][i] + At*(l+0.5) - Bt - Dtl - Ftl;
				
				// Add the regularized components to the numerical sum
				GSF->Fr_in->total->reg_numerical_sum[i] 		+= GSF->Fr_in->total->l_reg[l][i];
				GSF->Fr_out->total->reg_numerical_sum[i]		+= GSF->Fr_out->total->l_reg[l][i];

				GSF->Ft_in->total->reg_numerical_sum[i] 		+= GSF->Ft_in->total->l_reg[l][i];
				GSF->Ft_out->total->reg_numerical_sum[i] 		+= GSF->Ft_out->total->l_reg[l][i];

				GSF->Fr_in->cons->reg_numerical_sum[i] 		+= GSF->Fr_in->cons->l_reg[l][i];
				GSF->Fr_out->cons->reg_numerical_sum[i]		+= GSF->Fr_out->cons->l_reg[l][i];

				GSF->Ft_in->cons->reg_numerical_sum[i] 		+= GSF->Ft_in->cons->l_reg[l][i];
				GSF->Ft_out->cons->reg_numerical_sum[i] 	+= GSF->Ft_out->cons->l_reg[l][i];

				GSF->huu_in->total->reg_numerical_sum[i]	+= GSF->huu_in->total->l_reg[l][i];
				GSF->huu_out->total->reg_numerical_sum[i]	+= GSF->huu_out->total->l_reg[l][i];

				GSF->hut_in->total->reg_numerical_sum[i]	+= GSF->hut_in->total->l_reg[l][i];
				GSF->hut_out->total->reg_numerical_sum[i]	+= GSF->hut_out->total->l_reg[l][i];

				GSF->huphi_in->total->reg_numerical_sum[i]	+= GSF->huphi_in->total->l_reg[l][i];
				GSF->huphi_out->total->reg_numerical_sum[i]	+= GSF->huphi_out->total->l_reg[l][i];

				GSF->Grut_in->total->reg_numerical_sum[i]	+= GSF->Grut_in->total->l_reg[l][i];
				GSF->Grut_out->total->reg_numerical_sum[i]	+= GSF->Grut_out->total->l_reg[l][i];

				GSF->Gruphi_in->total->reg_numerical_sum[i]	+= GSF->Gruphi_in->total->l_reg[l][i];
				GSF->Gruphi_out->total->reg_numerical_sum[i]	+= GSF->Gruphi_out->total->l_reg[l][i];

				GSF->Gphiur_in->total->reg_numerical_sum[i]	+= GSF->Gphiur_in->total->l_reg[l][i];
				GSF->Gphiur_out->total->reg_numerical_sum[i]	+= GSF->Gphiur_out->total->l_reg[l][i];

				if( fabs(GSF->beta_out->total->l_reg[l][i]) > 3e-12){
					tail_l_max = l;
					GSF->beta_in->total->reg_numerical_sum[i]	+= GSF->beta_in->total->l_reg[l][i];
					GSF->beta_out->total->reg_numerical_sum[i]	+= GSF->beta_out->total->l_reg[l][i];
				}

				// The dissipative piece does not need regularizing so just add it to the total
				GSF->Fr_in->reg_diss_total[i]				+= GSF->Fr_in->l_full_diss[l][i];
				GSF->Fr_out->reg_diss_total[i]				+= GSF->Fr_out->l_full_diss[l][i];

				GSF->Ft_in->reg_diss_total[i]				+= GSF->Ft_in->l_full_diss[l][i];
				GSF->Ft_out->reg_diss_total[i]				+= GSF->Ft_out->l_full_diss[l][i];
			}

		}

		//calculate the tail for the conservative piece of the GSF

		calculate_tail(l_max, num_chi, GSF->Fr_in->cons, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->Fr_out->cons, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->Ft_in->cons, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->Ft_out->cons, 5, 2, 6);

		calculate_tail(l_max, num_chi, GSF->huu_in->total, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->huu_out->total, 5, 2, 6);

		calculate_tail(l_max, num_chi, GSF->hut_in->total, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->hut_out->total, 5, 2, 6);

		calculate_tail(l_max, num_chi, GSF->huphi_in->total, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->huphi_out->total, 5, 2, 6);

		calculate_tail(l_max, num_chi, GSF->Grut_in->total, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->Grut_out->total, 5, 2, 6);

		calculate_tail(l_max, num_chi, GSF->Gruphi_in->total, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->Gruphi_out->total, 5, 2, 6);

		calculate_tail(l_max, num_chi, GSF->Gphiur_in->total, 5, 2, 6);
		calculate_tail(l_max, num_chi, GSF->Gphiur_out->total, 5, 2, 6);

		if(verbose)	printf("tail_l_max: %d\n", tail_l_max);

		calculate_tail(tail_l_max, num_chi, GSF->beta_in->total, 5, 2, 6);
		calculate_tail(tail_l_max, num_chi, GSF->beta_out->total, 5, 2, 6);

		// The total SF is the sum of the conservative and dissipative pieces
		for(i=0; i < 2*num_chi + 1; i++){
			GSF->Fr_in->total->reg_total[i] = GSF->Fr_in->cons->reg_total[i] + GSF->Fr_in->reg_diss_total[i];
			GSF->Fr_out->total->reg_total[i] = GSF->Fr_out->cons->reg_total[i] + GSF->Fr_out->reg_diss_total[i];

			GSF->Ft_in->total->reg_total[i] = GSF->Ft_in->cons->reg_total[i] + GSF->Ft_in->reg_diss_total[i];
			GSF->Ft_out->total->reg_total[i] = GSF->Ft_out->cons->reg_total[i] + GSF->Ft_out->reg_diss_total[i];
		}

}

void calculate_tail(int l_max, int num_chi, struct GSF_component_data_piece *GSF_comp_piece, int n, int Nmin, int Nmax)
{
	int i;
	for(i = 0; i <= num_chi; i++){

	//	double* data = (double*)calloc(l_max+1, sizeof(double));
		double data[l_max+1];
		int l;
		for(l = 0; l <= l_max; l++){
			data[l] = GSF_comp_piece->l_reg[l][i];
		}

		double tail = ultimate_tail( data, l_max, n, Nmin, Nmax);

		//FIXME add error estimates

		GSF_comp_piece->reg_tail[i] = tail;
		GSF_comp_piece->reg_total[i] 	= GSF_comp_piece->reg_numerical_sum[i] + tail;

	//	free(data);

	}

}

/* Nmin = 0 with only AB used, 1 when ABD used etc*/
double ultimate_tail(double *data, int lmax, int n, int Nmin, int Nmax) 
{
	double chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *c;
	 
	X = gsl_matrix_alloc (n+1, Nmax-Nmin);
	y = gsl_vector_alloc (n+1);

	c = gsl_vector_alloc (Nmax-Nmin);
	cov = gsl_matrix_alloc (Nmax-Nmin, Nmax-Nmin);

	// Fit the data to the `Detweiler polynomials'
	int l;
	for (l = 0; l <= n; l++)
	{
		double li = lmax - n + l;

		int i;
		for(i = 0; i < Nmax-Nmin; i++){
			gsl_matrix_set (X, l, i, pow(-4.0,-(i+Nmin+1))/(gsl_sf_poch(0.5-li, i+Nmin+1)*gsl_sf_poch(1.5+li,i+Nmin+1)) );
		}
		gsl_vector_set (y, l, data[lmax - n + l]);
	}
	 
	gsl_multifit_linear_workspace * fit_workspace = gsl_multifit_linear_alloc (n+1, Nmax-Nmin);
	gsl_multifit_linear (X, y, c, cov, &chisq, fit_workspace);
	gsl_multifit_linear_free (fit_workspace); 

	// Now to estimate the tail
	double SF_tail = 0;

	int en = Nmin+1;
	for(en=Nmin+1; en <= Nmax-Nmin; en++){
		SF_tail += gsl_vector_get(c, en-Nmin-1)*( -pow(-1.0, -en)*pow(4.0, 1-en)*(lmax+1)/( (1.0+2.0*lmax)*(3.0+2.0*lmax)*(2.0*en-1.0)*gsl_sf_poch(0.5-lmax,en-1.0)*gsl_sf_poch(2.5+lmax,en-1.0) ) );
	}

	gsl_matrix_free(cov);
	gsl_vector_free(c);
	gsl_matrix_free(X);
	gsl_vector_free(y);

	return SF_tail;
}

void calculate_hat_beta(struct GSF_data *GSF, struct orbital_params *orbit)
{
	int i = 10;
	double r0 = orbit->r_p[i];
	double f0 = 1.0 - 2.0/r0;
	
	double Omega 		= pow(r0, -1.5);
	double omega_p 	= pow(1.0 - 3.0/r0, 0.5);

	double alpha = -pow(r0*(r0-3.0), -0.5);

	double F_r_out = GSF->Fr_out->total->reg_total[i]/f0;
	double hut_out = GSF->hut_out->total->reg_total[i];
	double huphi_out = GSF->huphi_out->total->reg_total[i];
	double Gruphi_out = GSF->Gruphi_out->total->reg_total[i];
	double Gphiur_out = GSF->Gphiur_out->total->reg_total[i];
	double Grut_out = GSF->Grut_out->total->reg_total[i];

	double beta_out = -1.0/(2.0*r0*Omega)*(Gruphi_out - pow(omega_p,2.0)*r0*r0*Gphiur_out + r0*r0*Omega/f0*Grut_out + Omega/f0*hut_out + 1.0/(r0*r0)*huphi_out);
	double beta_hat_out = beta_out - r0*omega_p/2.0 * F_r_out + alpha/(r0*omega_p);

	double F_r_in = GSF->Fr_in->total->reg_total[i]/f0;
	double hut_in = GSF->hut_in->total->reg_total[i];
	double huphi_in = GSF->huphi_in->total->reg_total[i];
	double Gruphi_in = GSF->Gruphi_in->total->reg_total[i];
	double Gphiur_in = GSF->Gphiur_in->total->reg_total[i];
	double Grut_in = GSF->Grut_in->total->reg_total[i];

	double beta_in = -1.0/(2.0*r0*Omega)*(Gruphi_in - pow(omega_p,2.0)*r0*r0*Gphiur_in + r0*r0*Omega/f0*Grut_in + Omega/f0*hut_in + 1.0/(r0*r0)*huphi_in);
	double beta_hat_in = beta_in - r0*omega_p/2.0 * F_r_in + alpha/(r0*omega_p);

	if(verbose) printf("hat{beta} (in/out) = %.16e / %.16e. Relative difference = %.3e\n", beta_hat_in, beta_hat_out, 1.0 - fabs(beta_hat_in/beta_hat_out));
	if(!verbose) printf("%.3lf\t%.16e\t%.3e\n", r0, beta_hat_out, 1.0 - fabs(beta_hat_in/beta_hat_out));

}

void	construct_beta_full(int l_max, int num_chi, struct GSF_data *GSF, struct orbital_params *orbit)
{
	int i;
	for(i=0; i < num_chi; i++){
		double r0 = orbit->r_p[i];
		double f0 = 1.0 - 2.0/r0;

		double Omega 		= pow(r0, -1.5);
		double omega_p 	= pow(1.0 - 3.0/r0, 0.5);
		int l;
		for(l = 0; l <= l_max; l++){

			double hut_out = GSF->hut_out->total->l_full[l][i];
			double huphi_out = GSF->huphi_out->total->l_full[l][i];
			double Gruphi_out = GSF->Gruphi_out->total->l_full[l][i];
			double Gphiur_out = GSF->Gphiur_out->total->l_full[l][i];
			double Grut_out = GSF->Grut_out->total->l_full[l][i];

			GSF->beta_out->total->l_full[l][i] = -1.0/(2.0*r0*Omega)*(Gruphi_out - pow(omega_p,2.0)*r0*r0*Gphiur_out + r0*r0*Omega/f0*Grut_out + Omega/f0*hut_out + 1.0/(r0*r0)*huphi_out);

			double hut_in = GSF->hut_in->total->l_full[l][i];
			double huphi_in = GSF->huphi_in->total->l_full[l][i];
			double Gruphi_in = GSF->Gruphi_in->total->l_full[l][i];
			double Gphiur_in = GSF->Gphiur_in->total->l_full[l][i];
			double Grut_in = GSF->Grut_in->total->l_full[l][i];

			GSF->beta_in->total->l_full[l][i] = -1.0/(2.0*r0*Omega)*(Gruphi_in - pow(omega_p,2.0)*r0*r0*Gphiur_in + r0*r0*Omega/f0*Grut_in + Omega/f0*hut_in + 1.0/(r0*r0)*huphi_in);
		}
	}
}

void calculate_hat_beta_from_beta(struct GSF_data *GSF, struct orbital_params *orbit)
{
	int i = 10;
	double r0 = orbit->r_p[i];
	double f0 = 1.0 - 2.0/r0;
	
	double omega_p 	= pow(1.0 - 3.0/r0, 0.5);

	double alpha = -pow(r0*(r0-3.0), -0.5);

	double F_r_out = GSF->Fr_out->total->reg_total[i]/f0;
	double beta_out = GSF->beta_out->total->reg_total[i];

	double beta_hat_out = beta_out - r0*omega_p/2.0 * F_r_out + alpha/(r0*omega_p);

	double F_r_in = GSF->Fr_in->total->reg_total[i]/f0;
	double beta_in = GSF->beta_in->total->reg_total[i];

	double beta_hat_in = beta_in - r0*omega_p/2.0 * F_r_in + alpha/(r0*omega_p);

	if(verbose) printf("hat{beta} (in/out) = %.16e / %.16e. Relative difference = %.3e\n", beta_hat_in, beta_hat_out, 1.0 - fabs(beta_hat_in/beta_hat_out));
	if(!verbose) printf("%.3lf\t%.16e\t%.3e\t%d\t%d\n", r0, beta_hat_out, 1.0 - fabs(beta_hat_in/beta_hat_out), l_max, tail_l_max);
}
















