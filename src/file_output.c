/** \file file_output.c
 *  Contains all the functions that output the computed data to files
 *  For a given orbit the data will be in data/$p/$e/
 */
#include "file_output.h"
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>


void output_Fr_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);
void output_Fr_to_file(struct GSF_data *GSF, struct orbital_params *orbit);

void output_Ft_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_Ft_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_Fr_cons_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);
void output_Fr_cons_to_file(struct GSF_data *GSF, struct orbital_params *orbit);

void output_Ft_cons_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_Ft_cons_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_Fr_diss_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);
void output_Fr_diss_to_file(struct GSF_data *GSF, struct orbital_params *orbit);

void output_Ft_diss_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_Ft_diss_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_huu_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_huu_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_hut_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_hut_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_huphi_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_huphi_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_Grut_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_Grut_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_Gruphi_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_Gruphi_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_Gphiur_to_file(struct GSF_data *GSF, struct orbital_params *orbit);
void output_Gphiur_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_Fl_full_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit);

void output_reg_data_to_file(int l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	output_Fr_to_file(GSF, orbit);
	output_Fr_tests_to_file(l_max, GSF, orbit);

	output_Ft_to_file(GSF, orbit);
	output_Ft_tests_to_file(l_max, GSF, orbit);

	output_Fr_cons_to_file(GSF, orbit);
	output_Fr_cons_tests_to_file(l_max, GSF, orbit);

	output_Ft_cons_to_file(GSF, orbit);
	output_Ft_cons_tests_to_file(l_max, GSF, orbit);

	output_Fr_diss_to_file(GSF, orbit);
	output_Fr_diss_tests_to_file(l_max, GSF, orbit);

	output_Ft_diss_to_file(GSF, orbit);
	output_Ft_diss_tests_to_file(l_max, GSF, orbit);

	output_huu_to_file(GSF, orbit);
	output_huu_tests_to_file(l_max, GSF, orbit);

	output_hut_to_file(GSF, orbit);
	output_hut_tests_to_file(l_max, GSF, orbit);

	output_huphi_to_file(GSF, orbit);
	output_huphi_tests_to_file(l_max, GSF, orbit);

	output_Grut_to_file(GSF, orbit);
	output_Grut_tests_to_file(l_max, GSF, orbit);

	output_Gruphi_to_file(GSF, orbit);
	output_Gruphi_tests_to_file(l_max, GSF, orbit);

	output_Gphiur_to_file(GSF, orbit);
	output_Gphiur_tests_to_file(l_max, GSF, orbit);

}


void output_GSF_data(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	printf("outputing data to files.\n");

	output_Fl_full_to_file(l_max, GSF, orbit);

}



/** Creates a directory structure */
void _mkdir(const char *dir) 
{
	char tmp[256];
	char *p = NULL;
	size_t len;
	snprintf(tmp, sizeof(tmp),"%s",dir);
	len = strlen(tmp);
	if(tmp[len - 1] == '/')
		  tmp[len - 1] = 0;
	for(p = tmp + 1; *p; p++)
		  if(*p == '/') {
		          *p = 0;
		          mkdir(tmp, S_IRWXU);
		          *p = '/';
		  }
	mkdir(tmp, S_IRWXU);
}

void output_run_info_to_file(double l_max, double run_time, int n_folding, int numprocs, int multi_threading, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "run_info");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating run info output file!\n\n");}
	else{
		fprintf(output_file, "Orbital parameters:\n");
			fprintf(output_file, "\tp = %lg\n", orbit->p);
			fprintf(output_file, "\te = %lg\n", orbit->e);
			fprintf(output_file, "\tE = %lg\n", orbit->E);
			fprintf(output_file, "\tL = %lg\n", orbit->L);
			fprintf(output_file, "\tr_max = %lg\n", orbit->r_max);
			fprintf(output_file, "\tr_min = %lg\n", orbit->r_min);
			fprintf(output_file, "\trs_max = %lg\n", orbit->rs_max);
			fprintf(output_file, "\trs_min = %lg\n", orbit->rs_min);
			fprintf(output_file, "\tT_r = %lg\n", orbit->T_r);
			fprintf(output_file, "\tOmega_r = %lg\n", orbit->Omega_r);
			fprintf(output_file, "\tOmega_phi = %lg\n", orbit->Omega_phi);
			fprintf(output_file, "\tOmega_r/Omega_phi = %lg\n", orbit->Omega_r/orbit->Omega_phi);
		fprintf(output_file, "\n\n");
		fprintf(output_file, "l_max = %.0f\n", l_max);
		fprintf(output_file, "run time = %.0fs\n", run_time); 
		fprintf(output_file, "Number of processors used = %d\n", numprocs);
		fprintf(output_file, "Multi-threading the lm mode counter: ");
		if(multi_threading){
			fprintf(output_file, "enabled\n");
		}else{
			fprintf(output_file, "disabled\n");
		}
		fprintf(output_file, "Folding n modes for the m=0 modes: ");
		if(n_folding){
			fprintf(output_file, "enabled\n");
		}else{
			fprintf(output_file, "disabled\n");
		}
		fclose(output_file);
	}

}


void output_Fr_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Frl_in");
	strcat(fileloc2, "Frl_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Fr_in/out output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Fr_in->total->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Fr_out->total->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Frl_in/Frl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Frl_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Fr_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Fr_in->total->l_reg[l][i]/GSF->Fr_out->total->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}


void output_Fl_full_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;	

	char fileloc[200];

	FILE *output_file;
	sprintf(fileloc, "data/p%lf/e%lf/modes/", p, e);
	_mkdir(fileloc);
	int l;
	for(l=0; l <= l_max; l++){
		sprintf(fileloc, "data/p%lf/e%lf/modes/l%d", p, e, l);	
		output_file = fopen(fileloc, "w+");
		if(output_file == NULL) printf("Error creating output file: %s\n", fileloc);
		else{
			fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
			fprintf(output_file, "# format: chi F^r_in F^r_out F^t_in F^t_out huu_in huu_out\n");
			int i;
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.16lf %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", orbit->chi_p[i], GSF->Fr_in->total->l_full[l][i], GSF->Fr_out->total->l_full[l][i], GSF->Ft_in->total->l_full[l][i] ,GSF->Ft_out->total->l_full[l][i], GSF->huu_in->total->l_full[l][i], GSF->huu_out->total->l_full[l][i], GSF->hut_in->total->l_full[l][i], GSF->hut_out->total->l_full[l][i], GSF->huphi_in->total->l_full[l][i], GSF->huphi_out->total->l_full[l][i],  GSF->Grut_in->total->l_full[l][i], GSF->Grut_out->total->l_full[l][i], GSF->Gruphi_in->total->l_full[l][i], GSF->Gruphi_out->total->l_full[l][i], GSF->Gphiur_in->total->l_full[l][i], GSF->Gphiur_out->total->l_full[l][i], GSF->beta_tensor_in->total->l_full[l][i], GSF->beta_tensor_out->total->l_full[l][i]);
			}
			fclose(output_file);
		}

	}


}

void output_Fr_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Fr");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tF^r_in\tError\tF^r(numerical)\tF^r(tail) | F^r_out\tError\tF^r(numerical)\tF^r(tail) | 1-F^r_out/F^r_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->Fr_in->total->reg_total[i], GSF->Fr_in->total->reg_error[i], GSF->Fr_in->total->reg_numerical_sum[i], GSF->Fr_in->total->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->Fr_out->total->reg_total[i], GSF->Fr_out->total->reg_error[i], GSF->Fr_out->total->reg_numerical_sum[i], GSF->Fr_out->total->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->Fr_out->total->reg_total[i]/GSF->Fr_in->total->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_Ft_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Ft");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tF^t_in\tError\tF^t(numerical)\tF^t(tail) | F^t_out\tError\tF^t(numerical)\tF^t(tail) | 1-F^t_out/F^t_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->Ft_in->total->reg_total[i], GSF->Ft_in->total->reg_error[i], GSF->Ft_in->total->reg_numerical_sum[i], GSF->Ft_in->total->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->Ft_out->total->reg_total[i], GSF->Ft_out->total->reg_error[i], GSF->Ft_out->total->reg_numerical_sum[i], GSF->Ft_out->total->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->Ft_out->total->reg_total[i]/GSF->Ft_in->total->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_Ft_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Ftl_in");
	strcat(fileloc2, "Ftl_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Ft_in/out output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Ft_in->total->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Ft_out->total->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Ftl_in/Ftl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Ftl_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Ft_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Ft_in->total->l_reg[l][i]/GSF->Ft_out->total->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////Below is the conservative file output functions////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void output_Fr_cons_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Frl_cons_in");
	strcat(fileloc2, "Frl_cons_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Fr_in/out cons output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Fr_in->cons->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Fr_out->cons->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Frl_in/Frl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Frl_cons_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Fr_in_out_diff cons output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Fr_in->cons->l_reg[l][i]/GSF->Fr_out->cons->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}


void output_Fr_cons_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Fr_cons");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tF^r_in\tError\tF^r(numerical)\tF^r(tail) | F^r_out\tError\tF^r(numerical)\tF^r(tail) | 1-F^r_out/F^r_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->Fr_in->cons->reg_total[i], GSF->Fr_in->cons->reg_error[i], GSF->Fr_in->cons->reg_numerical_sum[i], GSF->Fr_in->cons->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->Fr_out->cons->reg_total[i], GSF->Fr_out->cons->reg_error[i], GSF->Fr_out->cons->reg_numerical_sum[i], GSF->Fr_out->cons->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->Fr_out->cons->reg_total[i]/GSF->Fr_in->cons->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_Ft_cons_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Ft_cons");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tF^t_in\tError\tF^t(numerical)\tF^t(tail) | F^t_out\tError\tF^t(numerical)\tF^t(tail) | 1-F^t_out/F^t_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->Ft_in->cons->reg_total[i], GSF->Ft_in->cons->reg_error[i], GSF->Ft_in->cons->reg_numerical_sum[i], GSF->Ft_in->cons->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->Ft_out->cons->reg_total[i], GSF->Ft_out->cons->reg_error[i], GSF->Ft_out->cons->reg_numerical_sum[i], GSF->Ft_out->cons->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->Ft_out->cons->reg_total[i]/GSF->Ft_in->cons->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_Ft_cons_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Ftl_cons_in");
	strcat(fileloc2, "Ftl_cons_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Ft_in/out cons output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Ft_in->cons->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Ft_out->cons->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Ftl_in/Ftl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Ftl_cons_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Ft_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Ft_in->cons->l_reg[l][i]/GSF->Ft_out->cons->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////Below is the dissipative file output functions/////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void output_Fr_diss_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Fr_diss");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating Fr_diss output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tF^r_in | F^r_out | 1-F^r_out/F^r_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e", orbit->chi_p[i], GSF->Fr_in->reg_diss_total[i]);
			fprintf(output_file, "  |  %.12e", GSF->Fr_out->reg_diss_total[i]);
			fprintf(output_file, " 	|  %.6e\n", 1.0 - GSF->Fr_out->reg_diss_total[i]/GSF->Fr_in->reg_diss_total[i]);
		}
		fclose(output_file);
	}

}

void output_Ft_diss_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Ft_diss");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating Ft_diss output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tF^t_in | F^t_out | 1-F^t_out/F^t_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e", orbit->chi_p[i], GSF->Ft_in->reg_diss_total[i] );
			fprintf(output_file, "  |  %.12e", GSF->Ft_out->reg_diss_total[i]);
			fprintf(output_file, " 	|  %.6e\n", 1.0 - GSF->Ft_out->reg_diss_total[i]/GSF->Ft_in->reg_diss_total[i]);
		}
		fclose(output_file);
	}

}

void output_Fr_diss_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Frl_diss_in");
	strcat(fileloc2, "Frl_diss_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Fr_in/out diss output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Fr_in->l_full_diss[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Fr_out->l_full_diss[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Frl_in/Frl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Frl_diss_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Fr_in_out_diff diss output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Fr_in->l_full_diss[l][i]/GSF->Fr_out->l_full_diss[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}

void output_Ft_diss_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Ftl_diss_in");
	strcat(fileloc2, "Ftl_diss_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Ft_in/out diss output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Ft_in->l_full_diss[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Ft_out->l_full_diss[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Frl_in/Frl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Ftl_diss_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Ft_in_out_diff diss output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Ft_in->l_full_diss[l][i]/GSF->Ft_out->l_full_diss[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// huu /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void output_huu_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "huu");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\thuu_in\tError\thuu_in(numerical)\thuu_in(tail) | huu_out\tError\thuu_out(numerical)\thuu_out(tail) | 1-huu_out/huu_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->huu_in->total->reg_total[i], GSF->huu_in->total->reg_error[i], GSF->huu_in->total->reg_numerical_sum[i], GSF->huu_in->total->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->huu_out->total->reg_total[i], GSF->huu_out->total->reg_error[i], GSF->huu_out->total->reg_numerical_sum[i], GSF->huu_out->total->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->huu_out->total->reg_total[i]/GSF->huu_in->total->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_huu_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "huul_in");
	strcat(fileloc2, "huul_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating huu_in/out cons output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->huu_in->total->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->huu_out->total->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - huul_in/huul_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "huul_cons_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating huu_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->huu_in->total->l_reg[l][i]/GSF->huu_out->total->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// hut /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void output_hut_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "hut");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\thut_in\tError\thut_in(numerical)\thut_in(tail) | hut_out\tError\thut_out(numerical)\thut_out(tail) | 1-hut_out/hut_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->hut_in->total->reg_total[i], GSF->hut_in->total->reg_error[i], GSF->hut_in->total->reg_numerical_sum[i], GSF->hut_in->total->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->hut_out->total->reg_total[i], GSF->hut_out->total->reg_error[i], GSF->hut_out->total->reg_numerical_sum[i], GSF->hut_out->total->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->hut_out->total->reg_total[i]/GSF->hut_in->total->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_hut_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "hutl_in");
	strcat(fileloc2, "hutl_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating hut_in/out cons output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->hut_in->total->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->hut_out->total->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - hutl_in/hutl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "hutl_cons_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating hut_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->hut_in->total->l_reg[l][i]/GSF->hut_out->total->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// huphi /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void output_huphi_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "huphi");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\thuphi_in\tError\thuphi_in(numerical)\thuphi_in(tail) | huphi_out\tError\thuphi_out(numerical)\thuphi_out(tail) | 1-huphi_out/huphi_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->huphi_in->total->reg_total[i], GSF->huphi_in->total->reg_error[i], GSF->huphi_in->total->reg_numerical_sum[i], GSF->huphi_in->total->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->huphi_out->total->reg_total[i], GSF->huphi_out->total->reg_error[i], GSF->huphi_out->total->reg_numerical_sum[i], GSF->huphi_out->total->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->huphi_out->total->reg_total[i]/GSF->huphi_in->total->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_huphi_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "huphil_in");
	strcat(fileloc2, "huphil_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating huphi_in/out cons output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->huphi_in->total->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->huphi_out->total->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - huphil_in/huphil_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "huphil_cons_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating huphi_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->huphi_in->total->l_reg[l][i]/GSF->huphi_out->total->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Grut /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void output_Grut_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Grut");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tGrut_in\tError\tGrut_in(numerical)\tGrut_in(tail) | Grut_out\tError\tGrut_out(numerical)\tGrut_out(tail) | 1-Grut_out/Grut_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->Grut_in->total->reg_total[i], GSF->Grut_in->total->reg_error[i], GSF->Grut_in->total->reg_numerical_sum[i], GSF->Grut_in->total->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->Grut_out->total->reg_total[i], GSF->Grut_out->total->reg_error[i], GSF->Grut_out->total->reg_numerical_sum[i], GSF->Grut_out->total->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->Grut_out->total->reg_total[i]/GSF->Grut_in->total->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_Grut_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Grutl_in");
	strcat(fileloc2, "Grutl_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Grut_in/out cons output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Grut_in->total->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Grut_out->total->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Grutl_in/Grutl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Grutl_cons_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Grut_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Grut_in->total->l_reg[l][i]/GSF->Grut_out->total->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Gruphi /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void output_Gruphi_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Gruphi");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tGruphi_in\tError\tGruphi_in(numerical)\tGruphi_in(tail) | Gruphi_out\tError\tGruphi_out(numerical)\tGruphi_out(tail) | 1-Gruphi_out/Gruphi_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->Gruphi_in->total->reg_total[i], GSF->Gruphi_in->total->reg_error[i], GSF->Gruphi_in->total->reg_numerical_sum[i], GSF->Gruphi_in->total->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->Gruphi_out->total->reg_total[i], GSF->Gruphi_out->total->reg_error[i], GSF->Gruphi_out->total->reg_numerical_sum[i], GSF->Gruphi_out->total->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->Gruphi_out->total->reg_total[i]/GSF->Gruphi_in->total->reg_total[i]);
		}
		fclose(output_file);
	}

}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Gphiur /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void output_Gphiur_to_file(struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	FILE *output_file;
	char fileloc[200];
	sprintf(fileloc, "data/p%lf/e%lf/GSF/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Gphiur");
	output_file = fopen(fileloc, "w+");
	if(output_file == NULL) {printf("Error creating output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n", p, e);
		fprintf(output_file, "# format: chi\tGphiur_in\tError\tGphiur_in(numerical)\tGphiur_in(tail) | Gphiur_out\tError\tGphiur_out(numerical)\tGphiur_out(tail) | 1-Gphiur_out/Gphiur_in\n");
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.6f  |  %.12e\t%.3e\t%.12e\t%.12e", orbit->chi_p[i], GSF->Gphiur_in->total->reg_total[i], GSF->Gphiur_in->total->reg_error[i], GSF->Gphiur_in->total->reg_numerical_sum[i], GSF->Gphiur_in->total->reg_tail[i]);
			fprintf(output_file, "  |  %.12e\t%.3e\t%.12e\t%.12e", GSF->Gphiur_out->total->reg_total[i], GSF->Gphiur_out->total->reg_error[i], GSF->Gphiur_out->total->reg_numerical_sum[i], GSF->Gphiur_out->total->reg_tail[i]);
			fprintf(output_file, " 	|  %.12e\n", 1.0 - GSF->Gphiur_out->total->reg_total[i]/GSF->Gphiur_in->total->reg_total[i]);
		}
		fclose(output_file);
	}

}

void output_Gphiur_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Gphiurl_in");
	strcat(fileloc2, "Gphiurl_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Gphiur_in/out cons output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Gphiur_in->total->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Gphiur_out->total->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Gphiurl_in/Gphiurl_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Gphiurl_cons_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Gphiur_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Gphiur_in->total->l_reg[l][i]/GSF->Gphiur_out->total->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}
void output_Gruphi_tests_to_file(double l_max, struct GSF_data *GSF, struct orbital_params *orbit)
{
	double p = orbit->p;
	double e = orbit->e;

	// Output the Frl_in and Frl_out
	FILE *output_file;
	FILE *output_file2;
	char fileloc[200];
	char fileloc2[200];
	sprintf(fileloc, "data/p%lf/e%lf/test/", p, e);
	sprintf(fileloc2, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc);
	strcat(fileloc, "Gruphil_in");
	strcat(fileloc2, "Gruphil_out");
	output_file = fopen(fileloc, "w+");
	output_file2 = fopen(fileloc2, "w+");
	if(output_file == NULL || output_file2 == NULL) {printf("Error creating Gruphi_in/out cons output file!\n\n");}
	else{
		fprintf(output_file, "# p=%lf, e=%lf.\n#", p, e);
		fprintf(output_file2, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file, "%.12f\t     ", orbit->chi_p[i]);
			fprintf(output_file2, "%.12f\t    ", orbit->chi_p[i]);
		}
		fprintf(output_file, "\n");
		fprintf(output_file2, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file, "%.12e\t", GSF->Gruphi_in->total->l_reg[l][i]);
				fprintf(output_file2, "%.12e\t", GSF->Gruphi_out->total->l_reg[l][i]);
			}
			fprintf(output_file, "\n");
			fprintf(output_file2, "\n");
		}
		fclose(output_file);
		fclose(output_file2);
	}

	//Output 1.0 - Gruphil_in/Gruphil_out
	FILE *output_file3;
	char fileloc3[200];
	sprintf(fileloc3, "data/p%lf/e%lf/test/", p, e);
	_mkdir(fileloc3);
	strcat(fileloc3, "Gruphil_cons_in_out_diff");
	output_file3 = fopen(fileloc3, "w+");
	if(output_file3 == NULL ) {printf("Error creating Gruphi_in_out_diff output file!\n\n");}
	else{
		fprintf(output_file3, "# p=%lf, e=%lf.\n#", p, e);
		int i;
		for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
			fprintf(output_file3, "%.6f\t", orbit->chi_p[i]);
		}
		fprintf(output_file3, "\n");
		int l;
		for(l=0; l <= l_max; l++){
			for(i=0; i < 2*NUM_CHI_VALUES + 1; i++){
				fprintf(output_file3, "%.6e\t", fabs(1.0 - GSF->Gruphi_in->total->l_reg[l][i]/GSF->Gruphi_out->total->l_reg[l][i]));
			}
			fprintf(output_file3, "\n");
		}
		fclose(output_file3);
	}
}
