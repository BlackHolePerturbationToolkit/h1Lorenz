void calculate_tensor_lm_mode_contribution_to_beta(double l_max, double **beta_tensor_full, struct lm_mode_data *lm_mode, struct orbital_params *orbit, int in_or_out)
{
	int l 				= (int)lm_mode->l;
	double m 			= lm_mode->m;
	double omega 	= lm_mode->omega;
	double E 			= orbit.E;
	double L			= orbit.L;

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
		double omega_p = 1/ut;

		struct fields_data fields;
		if(in_or_out == INNER_FIELDS) cset->construct_inner_fields(lm_mode, &fields, i);
		else if(in_or_out == OUTER_FIELDS) cset->construct_outer_fields(lm_mode, &fields, i);
		fields.m 					= m;
		fields.r 					= r0;
		fields.eff				= f0;
		fields.contra_ur 	= orbit_contra_ur(orbit->chi_p[i]); 
		fields.contra_ut	= ut;
		fields.contra_uphi	= uphi;
	
		complex Y 	= normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		complex YV1, YV2, YT1, YT2;
		if(l>0) YV1 = 1.0/(l*(l+1.0)) * normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		if(l>0) YV2 = I*m/(l*(l+1.0)) * normalized_associated_legendre_P(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);
		if(l>1) YT1 = factorial(l-2)/factorial(l+2)*(normalized_associated_legendre_P_theta_2nd_deriv(l,m,M_PI/2.0)* cexp(I*m*orbit->phi_p[i]) + m*m*normalized_associated_legendre_P(l,m, M_PI/2.0)* cexp(I*m*orbit->phi_p[i]));
		if(l>1) YT2 = 2.0*I*m*factorial(l-2)/factorial(l+2)*normalized_associated_legendre_P_theta_deriv(l,m,M_PI/2.0) * cexp(I*m*orbit->phi_p[i]);

		complex hrr 			= 1.0/(2.0*r0*f0*f0)*(fields.h1 - f_p*fields.h6)*Y;
		complex hrphi			= 1.0/(2.0*f0)*(fields.h5*YV2 - fields.h9*YV1);
		complex hrtheta 	= 1.0/(2.0*f0)*(fields.h5*YV1 + fields.h9*YV2);
		complex htphi			= 1.0/2.0*(fields.h4*YV2 - fields.h8*YV1);
		complex htr				=	1.0/(2.0*r0*f0)*fields.h2*Y;
		complex hphiphi		= 1.0/2.0*r0*(fields.h3*Y - fields.h7*YT1 - fields.h10*YT2);

		complex hrphidt		= -I*omega*hrphi;
		complex htrdt			= -I*omega*htr;
		complex htrdphi		= I*m*htr;
		complex hrphidphi	= I*m*hrphi;

		complex httdr			= -htt/(2.0*r0*r0) + 1.0/(2.0*r0)*(fields.h1_rs_deriv*f0 + f0*f0*fields.h6_rs_deriv + df0*fields.h6)*Y;
		complex htphidr		= f0/2.0*(fields.h4_rs_deriv*YV2 - fields.h8_rs_deriv*YV1);
		complex hphiphidr	= 1.0/r0*hphiphi + 1.0/2.0*r0*f0*(fields.h3_rs_deriv*Y - fields.h7_rs_deriv*YT1 - fields.h10_rs_deriv*YT2);

		complex Grut		= (2.0-r0)*E*hrr/pow(r0,3.0) + E*htrdt + (r0-2.0)*L/(2.0*pow(r0,3.0))*hrphidt -0.5*E*httdr + (2.0-r0)*L/(2.0*pow(r0,3.0))*htphidr + (r0-2.0)*L/(2.0*pow(r0,3.0))*htrdphi;

		complex Gphiur 	= -E/pow(r0*(r0-2.0),2.0)*htphi - L/pow(r0,5.0)*hphiphi + E/(2.0*r0*r0-4.0*r0)*hrphidt + E/(2.0*r0*r0-4.0*r0)*htphidr + L/(2.0*pow(r0,4.0)) + E/(4.0*r0-2.0*r0*r0);

		complex Gruphi	= pow(r0-2.0,2.0)*L/pow(r0,3.0)*hrr + 0.5*E*htphidr + (2.0-r0)*L/(2.0*pow(r0,3.0))*hphiphidr + 0.5*E*htrdphi + (r0-2.0)*L/pow(r0,3.0)*hrphidphi;

		complex hut			= ut*htt + uphi*htphi;

		complex huphi		= ut*htphi + uphi*hphiphi;

		complex F_r			= 0.5*(ut*ut*httdr + 2.0*ut*uphi*htphidr + uphi*uphi*hphiphidr);

		double beta			= creal(-1.0/(2.0*r0*omega)*(Gruphi - pow(omega_p,2.0)*r0*r0*Gphiur + r0*r0*omega/f0*Grut + omega/f0*hut +1.0/(r0*r0)*huphi));

		double beta_hat	= creal(beta - r0*omega_p/2.0*F_r);

		if(l<= l_max) beta_tensor_full[l][i] 	+= factor * beta;

	}



}
