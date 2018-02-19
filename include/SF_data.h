struct GSF_component_data_piece{
	double** l_full;								//!< The full field per l mode. The first index denotes the l mode, the second is the chi index.
	double** l_reg;									//!< The regularized field per l mode. The first index denotes the l mode, the second is the chi index.
	double* reg_total;								//!< The total GSF for this component
	double* reg_numerical_sum;						//!< The numerical calculated portion of this component of the GSF
	double* reg_tail;								//!< The tail contribution to this component of the GSF
	double* reg_error;								//!< The estimated tail error to this component of the GSF
};

struct GSF_component_data{
	struct GSF_component_data_piece *total;
	struct GSF_component_data_piece *cons;

	double** l_full_diss;
	double* reg_diss_total; 
};

struct GSF_data{
	struct GSF_component_data *Fr_in;
	struct GSF_component_data *Ft_in;

	struct GSF_component_data *Fr_out;
	struct GSF_component_data *Ft_out;

	struct GSF_component_data *huu_in;
	struct GSF_component_data *huu_out;

	struct GSF_component_data *hut_in;
	struct GSF_component_data *hut_out;

	struct GSF_component_data *huphi_in;
	struct GSF_component_data *huphi_out;

	struct GSF_component_data *Grut_in;
	struct GSF_component_data *Grut_out;

	struct GSF_component_data *Gruphi_in;
	struct GSF_component_data *Gruphi_out;

	struct GSF_component_data *Gphiur_in;
	struct GSF_component_data *Gphiur_out;

	struct GSF_component_data *beta_in;
	struct GSF_component_data *beta_out;

	struct GSF_component_data *beta_tensor_in;
	struct GSF_component_data *beta_tensor_out;
};

void 	setup_GSF_data_structure( struct GSF_data *GSF, int l_max, int num_chi );
void setup_GSF_component_data_structure(struct GSF_component_data *GSF_comp, int l_max, int num_chi);
