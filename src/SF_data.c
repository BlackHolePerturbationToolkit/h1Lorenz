#include <stdlib.h>
#include "SF_data.h"
#include "array_handling.h"

void setup_GSF_component_data_structure(struct GSF_component_data *GSF_comp, int l_max, int num_chi)
{
	GSF_comp->total 					= calloc(1, sizeof(struct GSF_component_data_piece));
	GSF_comp->cons						= calloc(1, sizeof(struct GSF_component_data_piece));

	GSF_comp->total->l_full							= alloc_2D_double_array(l_max+1, 2*num_chi+1);
	GSF_comp->total->l_reg 							= alloc_2D_double_array(l_max+1, 2*num_chi+1);
	GSF_comp->total->reg_total					= calloc(2*num_chi+1, sizeof(double));
	GSF_comp->total->reg_numerical_sum 	= calloc(2*num_chi+1, sizeof(double));
	GSF_comp->total->reg_error					= calloc(2*num_chi+1, sizeof(double));
	GSF_comp->total->reg_tail						= calloc(2*num_chi+1, sizeof(double));

	GSF_comp->cons->l_full							= alloc_2D_double_array(l_max+1, 2*num_chi+1);
	GSF_comp->cons->l_reg 							= alloc_2D_double_array(l_max+1, 2*num_chi+1);
	GSF_comp->cons->reg_total						= calloc(2*num_chi+1, sizeof(double));
	GSF_comp->cons->reg_numerical_sum 	= calloc(2*num_chi+1, sizeof(double));
	GSF_comp->cons->reg_error						= calloc(2*num_chi+1, sizeof(double));
	GSF_comp->cons->reg_tail						= calloc(2*num_chi+1, sizeof(double));

	GSF_comp->l_full_diss								= alloc_2D_double_array(l_max+1, 2*num_chi+1);
	GSF_comp->reg_diss_total						= calloc(2*num_chi+1, sizeof(double));
}


void setup_GSF_data_structure(struct GSF_data *GSF, int l_max, int num_chi)
{
	GSF->Fr_in	 	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Fr_in, l_max, num_chi);

	GSF->Fr_out 	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Fr_out, l_max, num_chi);

	GSF->Ft_in 		= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Ft_in, l_max, num_chi);

	GSF->Ft_out 	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Ft_out, l_max, num_chi);

	GSF->huu_in	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->huu_in, l_max, num_chi);
	
	GSF->huu_out	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->huu_out, l_max, num_chi);

	GSF->hut_in	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->hut_in, l_max, num_chi);
	
	GSF->hut_out	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->hut_out, l_max, num_chi);

	GSF->huphi_in	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->huphi_in, l_max, num_chi);
	
	GSF->huphi_out	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->huphi_out, l_max, num_chi);

	GSF->Grut_in	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Grut_in, l_max, num_chi);
	
	GSF->Grut_out	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Grut_out, l_max, num_chi);

	GSF->Gruphi_in	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Gruphi_in, l_max, num_chi);
	
	GSF->Gruphi_out	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Gruphi_out, l_max, num_chi);

	GSF->Gphiur_in	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Gphiur_in, l_max, num_chi);
	
	GSF->Gphiur_out	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->Gphiur_out, l_max, num_chi);



	GSF->beta_in	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->beta_in, l_max, num_chi);
	
	GSF->beta_out	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->beta_out, l_max, num_chi);

	GSF->beta_tensor_in	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->beta_tensor_in, l_max, num_chi);
	
	GSF->beta_tensor_out	= malloc(sizeof(struct GSF_component_data));
	setup_GSF_component_data_structure(GSF->beta_tensor_out, l_max, num_chi);
}
