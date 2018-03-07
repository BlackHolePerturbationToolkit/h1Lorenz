/* Construct the asymptotically flat monopole */
/* From h1ret-l0m0.nb*/

#include "GSF_ecc.h"

#define Power(x, y) (pow((double)(x), (double)(y)))
#define Sqrt(x)     (sqrt((double)(x)))
#define Log(x)      (log((double)(x)))
#define Pi M_PI

double h1MonoIn(double r0, double r){
	return (-4*Sqrt(Pi)*Power(-2 + r,3)*(2 + r)*(-4 + (-3 + r0)*Log((-2 + r0)/r0)))/(3.*Power(r,3)*Sqrt((-3 + r0)*r0));
}

double h3MonoIn(double r0, double r){
	return (8*Sqrt(Pi)*(-32 + Power(r,3) - (-8 + Power(r,3))*(-3 + r0)*Log((-2 + r0)/r0)))/(3.*Power(r,2)*Sqrt((-3 + r0)*r0));
}

double h6MonoIn(double r0, double r){
	return (4*Sqrt(Pi)*(-64 + 2*Power(r,3) + (16 + Power(r,3))*(-3 + r0)*Log((-2 + r0)/r0)))/(3.*Power(r,2)*Sqrt((-3 + r0)*r0));
}

double h1MonoOut(double r0, double r){
	return (4*Sqrt(Pi)*(-2*(-8 + r0)*(-4 + r0) - 6*Power(r,2)*(-3 + r0) + 4*Power(r,3)*(-3 + r0) + 2*r*(20 + (-8 + r0)*r0) - (-4 + r)*Power(r,3)*(-3 + r0)*Log((-2 + r)/r) + 16*(-1 + r)*(-3 + r0)*Log(r0/(-2 + r))))/(3.*Power(r,3)*Sqrt((-3 + r0)*r0));
}

double h3MonoOut(double r0, double r){
	return (-8*Sqrt(Pi)*(32 + 3*r*(4 + r) - 12*r0 - r*(4 + r)*r0 + Power(r0,2) + Power(r,3)*(-3 + r0)*Log((-2 + r)/r) + 8*(-3 + r0)*Log(r0/(-2 + r))))/(3.*Power(r,2)*Sqrt((-3 + r0)*r0));	
}

double h6MonoOut(double r0, double r){
	return (4*Sqrt(Pi)*(-2*(-8 + r0)*(-4 + r0) + 8*r*(-3 + r0) + 2*Power(r,2)*(-3 + r0) + Power(r,3)*(-3 + r0)*Log((-2 + r)/r) + 16*(-3 + r0)*Log((-2 + r)/r0)))/(3.*Power(r,2)*Sqrt((-3 + r0)*r0));	
}

double dh1MonoIn(double r0, double r){
	return (-4*Sqrt(Pi)*Power(-2 + r,3)*(-4 + (-3 + r0)*Log((-2 + r0)/r0)))/(3.*Power(r,3)*Sqrt((-3 + r0)*r0)) + (4*Sqrt(Pi)*Power(-2 + r,3)*(2 + r)*(-4 + (-3 + r0)*Log((-2 + r0)/r0)))/(Power(r,4)*Sqrt((-3 + r0)*r0)) - 
   (4*Sqrt(Pi)*Power(-2 + r,2)*(2 + r)*(-4 + (-3 + r0)*Log((-2 + r0)/r0)))/(Power(r,3)*Sqrt((-3 + r0)*r0));
}

double dh3MonoIn(double r0, double r){
	return (8*Sqrt(Pi)*(3*Power(r,2) - 3*Power(r,2)*(-3 + r0)*Log((-2 + r0)/r0)))/(3.*Power(r,2)*Sqrt((-3 + r0)*r0)) - (16*Sqrt(Pi)*(-32 + Power(r,3) - (-8 + Power(r,3))*(-3 + r0)*Log((-2 + r0)/r0)))/(3.*Power(r,3)*Sqrt((-3 + r0)*r0));
}

double dh6MonoIn(double r0, double r){
	return (4*Sqrt(Pi)*(6*Power(r,2) + 3*Power(r,2)*(-3 + r0)*Log((-2 + r0)/r0)))/(3.*Power(r,2)*Sqrt((-3 + r0)*r0)) - (8*Sqrt(Pi)*(-64 + 2*Power(r,3) + (16 + Power(r,3))*(-3 + r0)*Log((-2 + r0)/r0)))/(3.*Power(r,3)*Sqrt((-3 + r0)*r0));
}

double dh1MonoOut(double r0, double r){
	return (4*Sqrt(Pi)*((-16*(-1 + r)*(-3 + r0))/(-2 + r) - 12*r*(-3 + r0) + 12*Power(r,2)*(-3 + r0) - ((-((-2 + r)/Power(r,2)) + 1/r)*(-4 + r)*Power(r,4)*(-3 + r0))/(-2 + r) + 2*(20 + (-8 + r0)*r0) - 
        3*(-4 + r)*Power(r,2)*(-3 + r0)*Log((-2 + r)/r) - Power(r,3)*(-3 + r0)*Log((-2 + r)/r) + 16*(-3 + r0)*Log(r0/(-2 + r))))/(3.*Power(r,3)*Sqrt((-3 + r0)*r0)) - 
   (4*Sqrt(Pi)*(-2*(-8 + r0)*(-4 + r0) - 6*Power(r,2)*(-3 + r0) + 4*Power(r,3)*(-3 + r0) + 2*r*(20 + (-8 + r0)*r0) - (-4 + r)*Power(r,3)*(-3 + r0)*Log((-2 + r)/r) + 16*(-1 + r)*(-3 + r0)*Log(r0/(-2 + r))))/
    (Power(r,4)*Sqrt((-3 + r0)*r0));
}

double dh3MonoOut(double r0, double r){
	return (-8*Sqrt(Pi)*(3*r + 3*(4 + r) - (8*(-3 + r0))/(-2 + r) + ((-((-2 + r)/Power(r,2)) + 1/r)*Power(r,4)*(-3 + r0))/(-2 + r) - r*r0 - (4 + r)*r0 + 3*Power(r,2)*(-3 + r0)*Log((-2 + r)/r)))/(3.*Power(r,2)*Sqrt((-3 + r0)*r0)) + 
   (16*Sqrt(Pi)*(32 + 3*r*(4 + r) - 12*r0 - r*(4 + r)*r0 + Power(r0,2) + Power(r,3)*(-3 + r0)*Log((-2 + r)/r) + 8*(-3 + r0)*Log(r0/(-2 + r))))/(3.*Power(r,3)*Sqrt((-3 + r0)*r0));
}

double dh6MonoOut(double r0, double r){
	return (4*Sqrt(Pi)*(8*(-3 + r0) + (16*(-3 + r0))/(-2 + r) + 4*r*(-3 + r0) + ((-((-2 + r)/Power(r,2)) + 1/r)*Power(r,4)*(-3 + r0))/(-2 + r) + 3*Power(r,2)*(-3 + r0)*Log((-2 + r)/r)))/(3.*Power(r,2)*Sqrt((-3 + r0)*r0)) - 
   (8*Sqrt(Pi)*(-2*(-8 + r0)*(-4 + r0) + 8*r*(-3 + r0) + 2*Power(r,2)*(-3 + r0) + Power(r,3)*(-3 + r0)*Log((-2 + r)/r) + 16*(-3 + r0)*Log((-2 + r)/r0)))/(3.*Power(r,3)*Sqrt((-3 + r0)*r0));
}

void construct_static_monopole (struct n_mode_data *n_mode, struct orbital_params *orbit )
{
	double r0 = orbit->p;
	double r;
	
	//Store the fields to the left of the particle
	int i;
	for(i = 0; i <= orbit->r0_grid_index; i++){
		r = orbit->grid[i];		
		n_mode->inhom_data[0][0][i] = h1MonoIn(r0, r);					// i = 1 field
		n_mode->inhom_data[0][1][i] = 0;					
		n_mode->inhom_data[0][2][i] = dh1MonoIn(r0, r);					
		n_mode->inhom_data[0][3][i] = 0;					
		
		n_mode->inhom_data[1][0][i] = h3MonoIn(r0, r);					// i = 3 field
		n_mode->inhom_data[1][1][i] = 0;					
		n_mode->inhom_data[1][2][i] = dh3MonoIn(r0, r);					
		n_mode->inhom_data[1][3][i] = 0;		
		
		n_mode->inhom_data[2][0][i] = h6MonoIn(r0, r);					// i = 6 field
		n_mode->inhom_data[2][1][i] = 0;					
		n_mode->inhom_data[2][2][i] = dh6MonoIn(r0, r);					
		n_mode->inhom_data[2][3][i] = 0;		
	}
	
	//Store the fields to the right of the particle
	for(i = orbit->r0_grid_index+1; i < orbit->gridsize; i++){
		r = orbit->grid[i];		
		
		n_mode->inhom_data[0][0][i] = h1MonoOut(r0, r);					// i = 1 field
		n_mode->inhom_data[0][1][i] = 0;		
		n_mode->inhom_data[0][2][i] = dh1MonoOut(r0, r);		
		n_mode->inhom_data[0][3][i] = 0;		
		
		n_mode->inhom_data[1][0][i] = h3MonoOut(r0, r);					// i = 3 field
		n_mode->inhom_data[1][1][i] = 0;		
		n_mode->inhom_data[1][2][i] = dh3MonoOut(r0, r);		
		n_mode->inhom_data[1][3][i] = 0;		
		
		n_mode->inhom_data[2][0][i] = h6MonoOut(r0, r);					// i = 6 field
		n_mode->inhom_data[2][1][i] = 0;		
		n_mode->inhom_data[2][2][i] = dh6MonoOut(r0, r);		
		n_mode->inhom_data[2][3][i] = 0;		
	}			
	
	// The value of the field at the particle from the right is stored at the end of the inhom_data arrays
	
	i = orbit->gridsize;

	n_mode->inhom_data[0][0][i] = h1MonoOut(r0, r0);					// i = 1 field
	n_mode->inhom_data[0][1][i] = 0;		
	n_mode->inhom_data[0][2][i] = dh1MonoOut(r0, r0);		
	n_mode->inhom_data[0][3][i] = 0;		
	
	n_mode->inhom_data[1][0][i] = h3MonoOut(r0, r0);					// i = 3 field
	n_mode->inhom_data[1][1][i] = 0;		
	n_mode->inhom_data[1][2][i] = dh3MonoOut(r0, r0);		
	n_mode->inhom_data[1][3][i] = 0;		
	
	n_mode->inhom_data[2][0][i] = h6MonoOut(r0, r0);					// i = 6 field
	n_mode->inhom_data[2][1][i] = 0;		
	n_mode->inhom_data[2][2][i] = dh6MonoOut(r0, r0);		
	n_mode->inhom_data[2][3][i] = 0;	

}
