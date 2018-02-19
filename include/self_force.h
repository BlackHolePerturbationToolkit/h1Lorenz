#include "GSF_ecc.h"

void calculate_lm_mode_contribution_to_Fr(int l_max, struct lm_mode_data *lm_mode, double** Frl_full, struct orbital_params *orbit, int in_or_out);
void calculate_lm_mode_contribution_to_Ft(int l_max, struct lm_mode_data *lm_mode, double** Ftl_full, struct orbital_params *orbit, int in_or_out);
void calculate_lm_mode_contribution_to_huul(double l_max, double **huul, struct lm_mode_data *lm_mode, struct orbital_params *orbit, int in_or_out);

void calculate_tensor_lm_mode_contribution_to_beta(double l_max, double **beta_tensor_full, struct lm_mode_data *lm_mode, struct orbital_params *orbit, int in_or_out);
void calculate_tensor_lm_mode_contribution_to_huul(double l_max, double **huul_full, struct lm_mode_data *lm_mode, struct orbital_params *orbit, int in_or_out);

void calculate_tensor_lm_mode_contribution_to_hR_and_FrR(double l_max, double **huul_R, double **FrlR, struct n_mode_data *n_mode, struct orbital_params *orbit);

void construct_monopole_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
void construct_odd_dipole_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
void construct_even_dipole_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
void construct_odd_modes_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
void construct_even_modes_from_inner_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);

void construct_monopole_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
void construct_odd_dipole_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
void construct_even_dipole_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
void construct_odd_modes_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
void construct_even_modes_from_outer_fields(struct lm_mode_data* lm_mode, struct fields_data* fields, int i);
