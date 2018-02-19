#include <complex.h>

struct boundary_conditions{
	double omega;
	double l;
	double complex* R;
	double complex* R_rs_deriv;
	double rs_boundary;
};

void monopole_out_bcs(double* a0s, struct boundary_conditions *bcs);
void monopole_in_bcs(double* a0s, struct boundary_conditions *bcs);

void odd_dipole_out_bcs(double* a0s, struct boundary_conditions *bcs);
void odd_dipole_in_bcs(double* a0s, struct boundary_conditions *bcs);

void even_dipole_out_bcs(double* a0s, struct boundary_conditions *bcs);
void even_dipole_in_bcs(double* a0s, struct boundary_conditions *bcs);

void even_out_bcs(double* a0s, struct boundary_conditions *bcs);
void even_in_bcs(double* a0s, struct boundary_conditions *bcs);

void even_static_out_bcs(double* a0s, struct boundary_conditions *bcs);
void even_static_in_bcs(double* a0s, struct boundary_conditions *bcs);

void odd_out_bcs(double* a0s, struct boundary_conditions *bcs);
void odd_in_bcs(double* a0s, struct boundary_conditions *bcs);



