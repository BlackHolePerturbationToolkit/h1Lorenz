#include "GSF_ecc.h"

void output_run_info_to_file(double l_max, double run_time, int n_folding, int numprocs, int multi_threading, struct GSF_data *GSF, struct orbital_params *orbit);
void output_reg_data_to_file(int l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_GSF_data(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);


