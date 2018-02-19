#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "GSF_ecc.h"

void read_in_src_dataset(hid_t file_id, double complex **src_ds, int i, int l, int m, int gridsize)
{
    hsize_t     dims[2];
    herr_t      status;

	char dataset_name[30];
	sprintf(dataset_name, "/src i=%d l=%d m=%d", i, l, m);

    /* get the dimensions of the dataset */
    status = H5LTget_dataset_info(file_id, dataset_name, dims, NULL, NULL);
	if(dims[0] != gridsize+1){
		printf("Error: src grid does not match output grid\n\n");
		exit(1);
	}
	
	//printf("dims: %d %d\n", dims[0], dims[1]);
	
	double *tmp = calloc(dims[0]*dims[1], sizeof(double));
	
    /* read in grid */
	 status = H5LTread_dataset_double(file_id, dataset_name, tmp);

	 int k;
	 for(k=0; k < dims[0]; k++){
		 (*src_ds)[k] = tmp[2*k] + I*tmp[2*k+1];
		
	 }
	
	 free(tmp);
}

struct interpolant{
	gsl_spline *spline;
	gsl_interp_accel *acc;
};


double integrand(double r, void *params)
{
	struct interpolant *interp = (struct interpolant *)params;
	
	return gsl_spline_eval(interp->spline, r, interp->acc);
}

double integrate_interpolant(gsl_spline *spline, gsl_interp_accel *acc, double r_min, double r_max)
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  
    double result, error;

	struct interpolant interp = {spline, acc};

    gsl_function F;
    F.function = &integrand;
    F.params = &interp;

	//FIXME put in more sensible error thresholds based upon side of integrand

    gsl_integration_qags (&F, r_min, r_max, 1e-11, 1e-11, 1000, w, &result, &error); 

    gsl_integration_workspace_free (w);
	
	return result;
}



void calculate_h1_R_on_grid(struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	
	 hid_t       file_id;
    /* open file from ex_lite1.c */
    file_id = H5Fopen ("input/src1.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	
	int l = n_mode->l;
	int m = n_mode->m;
	
	int nf = n_mode->cset->num_coupled_fields;
	if((l+m)%2!=0 && m ==0) nf = 1;
	int fields[nf];
	if((l+m) %2 == 0 ){
		if(m!=0){fields[0] = 1; fields[1] = 3; fields[2]= 5; fields[3] = 6; fields[4] = 7;}
		else{fields[0] = 1; fields[1] = 3; fields[2] = 5;}
	}else{
		if(m!=0){fields[0] = 9; fields[1] = 10;}
		else{fields[0] = 8;}
	}
	
	// read in the data from the source file
	int i;
	for(i = 0; i<nf; i++){
		read_in_src_dataset(file_id, &n_mode->src[i], fields[i], l, m, orbit->gridsize);
	}
	//printf("read in test: %.12e %.12e", n_mode->src[0][4080]);
		
    /* close file */
    H5Fclose (file_id);
	
	int lhs = 1;
	
	int j,k;
	for(k = orbit->ra_grid_index; k < orbit->rb_grid_index; k++){
		
		//Form the matrix system

		gsl_matrix_complex * m = gsl_matrix_complex_alloc (2*nf, 2*nf);
		
		int k2 = (orbit->gridsize - orbit->ra_grid_index - 1) - (k - orbit->ra_grid_index);

		for (i = 0; i < nf; i++){
			for (j = 0; j < nf; j++){
				gsl_complex minus_in_field 				= gsl_complex_rect(-n_mode->hom_data_in[i][j][0][k], -n_mode->hom_data_in[i][j][1][k]);
				gsl_complex minus_in_field_rs_deriv 	= gsl_complex_rect(-n_mode->hom_data_in[i][j][2][k], -n_mode->hom_data_in[i][j][3][k]);
				gsl_complex out_field 					= gsl_complex_rect(n_mode->hom_data_out[i][j][0][k2], n_mode->hom_data_out[i][j][1][k2]);
				gsl_complex out_field_rs_deriv 			= gsl_complex_rect(n_mode->hom_data_out[i][j][2][k2], n_mode->hom_data_out[i][j][3][k2]);

				gsl_matrix_complex_set (m, i, j, minus_in_field);
				gsl_matrix_complex_set (m, i, j + nf, out_field);
				gsl_matrix_complex_set (m, i + nf, j, minus_in_field_rs_deriv);
				gsl_matrix_complex_set (m, i + nf, j + nf, out_field_rs_deriv);
			}
		}
		
		// Form the source vector
		gsl_vector_complex * J = gsl_vector_complex_alloc(2*nf);
		for( i = 0; i < nf; i++){
			gsl_complex zero 	= gsl_complex_rect(0,0);
			gsl_complex	source = gsl_complex_rect(creal(n_mode->src[i][k]), cimag(n_mode->src[i][k]));
			if(k == orbit->r0_grid_index && lhs == 0){
				source = gsl_complex_rect(creal(n_mode->src[i][orbit->gridsize]), cimag(n_mode->src[i][orbit->gridsize]));
				//printf("testerismo2: %d %.14e\n", k - orbit->r0_grid_index, n_mode->src[i][orbit->gridsize ]);
			}
			
			gsl_vector_complex_set(J, i, zero);
			gsl_vector_complex_set(J, i+nf, source); 
		}
		
		// Invert the system to get Phi^{-1} * J = C
		gsl_vector_complex *scaling_coeffs = gsl_vector_complex_alloc (2*nf);

		int s;

		gsl_permutation * p = gsl_permutation_alloc (2*nf);
		gsl_linalg_complex_LU_decomp (m, p, &s);
		gsl_linalg_complex_LU_solve (m, p, J, scaling_coeffs);
		
		if(k == orbit->r0_grid_index) {
			for(i=0; i < 2*nf; i++){
	//			printf("result: %.16e\n", GSL_REAL(gsl_vector_complex_get(scaling_coeffs, i)));
			}
		}
		
		for(i=0; i<2*nf; i++){
			if(k <= orbit->r0_grid_index && lhs == 1){
				n_mode->scaling_coeffs_re_lhs[i][k] = GSL_REAL(gsl_vector_complex_get(scaling_coeffs, i));
				n_mode->scaling_coeffs_im_lhs[i][k] = GSL_IMAG(gsl_vector_complex_get(scaling_coeffs, i));
			}else{
				//if(k - orbit->r0_grid_index < 3 && i == 0) printf("testerismo: %d %.14e\n", k - orbit->r0_grid_index, GSL_REAL(gsl_vector_complex_get(scaling_coeffs, i)));
				n_mode->scaling_coeffs_re_rhs[i][k - orbit->r0_grid_index] = GSL_REAL(gsl_vector_complex_get(scaling_coeffs, i));
				n_mode->scaling_coeffs_im_rhs[i][k - orbit->r0_grid_index] = GSL_IMAG(gsl_vector_complex_get(scaling_coeffs, i));
			}
		}
		
		gsl_matrix_complex_free(m);
		gsl_permutation_free(p);
		gsl_vector_complex_free(scaling_coeffs);
	
		if(k == orbit->r0_grid_index && lhs !=0){
			lhs = 0; k--;
		}
		
	}
	
    gsl_interp_accel *acc 	= gsl_interp_accel_alloc ();
    gsl_spline *spline		= gsl_spline_alloc (gsl_interp_cspline, orbit->r0_grid_index + 1);
	
	double r0 = orbit->p;

	double complex *co_total 	= calloc(nf,sizeof(double complex));
	double complex *co_total_dr = calloc(nf,sizeof(double complex));
	for(j = 0; j < nf; j++){
    	gsl_spline_init (spline, orbit->grid, n_mode->scaling_coeffs_re_lhs[nf+ j], orbit->r0_grid_index+1);
		double complex co = integrate_interpolant(spline, acc, r0-2, r0 );
		gsl_spline_init (spline, orbit->grid, n_mode->scaling_coeffs_im_lhs[nf+ j], orbit->r0_grid_index+1);
		co += I * integrate_interpolant(spline, acc, r0-2, r0);
		for(i = 0; i < nf; i++)
		{
			int k0 = (orbit->gridsize - orbit->ra_grid_index - 1) - (orbit->r0_grid_index - orbit->ra_grid_index);
			co_total[i] 	+=  co*(n_mode->hom_data_out[i][j][0][k0] + I * n_mode->hom_data_out[i][j][1][k0]);
			co_total_dr[i] 	+=  co*(n_mode->hom_data_out[i][j][2][k0] + I * n_mode->hom_data_out[i][j][3][k0]);
		}
	}
	
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);

	acc		= gsl_interp_accel_alloc();
	spline 	= gsl_spline_alloc(gsl_interp_cspline, orbit->gridsize - orbit->r0_grid_index);
	
	double *grid_rhs = &(orbit->grid[orbit->r0_grid_index]);
	
	//for(i = 0; i < 10; i++)printf("grid: %.12e %.12e\n", grid_rhs[i], n_mode->scaling_coeffs_re_rhs[0][i]);
		
	double complex *ci_total	= calloc(nf, sizeof(double complex));
	double complex *ci_total_dr = calloc(nf, sizeof(double complex));
	for(j = 0; j < nf; j++){
		gsl_spline_init(spline, grid_rhs, n_mode->scaling_coeffs_re_rhs[j], orbit->gridsize - orbit->r0_grid_index);
		double complex ci = integrate_interpolant(spline, acc, r0, r0+2);
		gsl_spline_init(spline, grid_rhs, n_mode->scaling_coeffs_im_rhs[j], orbit->gridsize - orbit->r0_grid_index);
		ci += I * integrate_interpolant(spline, acc, r0, r0+2);
		for(i = 0; i < nf; i++){
			int k0 = orbit->r0_grid_index;
			ci_total[i]		+= ci*(n_mode->hom_data_in[i][j][0][k0] + I * n_mode->hom_data_in[i][j][1][k0]);
			ci_total_dr[i] 	+= ci*(n_mode->hom_data_in[i][j][2][k0] + I * n_mode->hom_data_in[i][j][3][k0]);
		}
	}
	
	//printf("testization: %.14e\n", n_mode->scaling_coeffs_re_lhs[0][0]);
		
	//double yi = gsl_spline_eval(spline, 10, acc);
	
	//double yi = integrate_interpolant(spline, acc, 8, 10 );
	
	//for(i=0;i<nf;i++) printf("test: %.14e %.14e %.14e\n", creal(co_total[i]), creal(ci_total[i]), creal(co_total[i]+ci_total[i]));
	for(i=0; i<nf; i++){
		//printf("test: %.14e\n", creal(co_total[i] + ci_total[i]));
		n_mode->hR[i][0][orbit->r0_grid_index] = creal(co_total[i] + ci_total[i]);
		n_mode->hR[i][1][orbit->r0_grid_index] = cimag(co_total[i] + ci_total[i]);
		n_mode->hR[i][2][orbit->r0_grid_index] = creal(co_total_dr[i] + ci_total_dr[i]);
		n_mode->hR[i][3][orbit->r0_grid_index] = cimag(co_total_dr[i] + ci_total_dr[i]);
	}
	
	n_mode->cset->construct_gauge_fields_on_grid(n_mode->hR, n_mode, orbit, n_mode->cset->num_coupled_fields);
	for(i=0; i<nf+n_mode->cset->num_gauge_fields; i++){
	//		printf("test: %.14e\n", n_mode->hR[i][0][orbit->r0_grid_index]);
	}
	
	
	gsl_interp_accel_free(acc);
	free(co_total);
	free(ci_total);
	free(co_total_dr);
	free(ci_total_dr);
	
}



































