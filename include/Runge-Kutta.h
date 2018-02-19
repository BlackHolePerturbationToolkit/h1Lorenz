#include "GSF_ecc.h"

struct field_data
{
	double complex** R_out;				// The first index is the field index, the second the basis index
	double complex** R_out_rs_deriv;
	double complex** R_in;			
	double complex** R_in_rs_deriv;
};

void integrate_field_equations(struct orbital_params *orbit, struct n_mode_data *n_mode);

void integrate_field_equations_for_scaling(double rs, struct field_data *fields, struct coupled_set *cset, struct orbital_params *orbit, struct n_mode_data *n_mode);
