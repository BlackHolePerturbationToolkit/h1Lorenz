/** \file GSF_ecc.c
 *  Frequency domain Schwarzschild Gravitional Self Force calculator using MPI
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/stat.h>
#include <mpi.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "GSF_ecc.h"
#include "sources.h"
#include "Runge-Kutta.h"
#include "scaling_coefficients.h"
#include "ode_systems.h"
#include "gauge_fields.h"
#include "gauge_fields_on_grid.h"
#include "static_monopole.h"
#include "convergence_tests.h"
#include "file_output.h"
#include "grid.h"
#include "h1R_on_grid.h"

double p 		= 0.0;
double e 		= 0.0;					//!< Although this code has a few routines that work for eccentric orbits, in general it only works for circular ones

double M 		= 1.0;					//!< All black holes have a mass of M = 1
int l_min 		= 0;				
int l_max 		= 0;
int m_min		= 0;
int abs_n_max 	= 20;					//!< The absolute maximum n that will be computed. The code will attempt to measure the convergence but sometimes fields get stuck before they reach the threshold.
int n_min		= 0;			

double MATCHING_INTEGRATION_ACCURACY 		= 1e-12;		//!< Sets the accuracy on the integration that calculates the scaling coefficients (method of EHS)
const double RUNGE_KUTTA_ACCURACY 			= 1e-12; 		//!< Sets the relative and absolute accuracy of the GSL Runge-Kutta algorithm
const double DESIRED_INVERT_ACCURACY 		= 1e-12;		//!< Sets the desired accuacy on the tortoise inverter
const double DESIRED_FUNCTION_ACCURACY 		= 1e-12;		//!< Sets the desired function accuracy for the wrappers around various GSL functions

const double n_convergence_threshold 		= 3e-8;			//!< The desired accuracy when comparing the jump in the fields with the time-domain source

const int convergence_test_type				= INNER_VS_OUTER; //!< Sets the type of convergence test to be used, either INNER_VS_OUTER or OUTER_CONTRIBUTIONS

int multi_threading = 0;									//!< Stores the multi_threading state of the code, which is automatically detected
int n_folding 		= 1;									//!< Set to one to enable n folding for the m = 0 modes
int file_output 	= 0;									//!< Set to one to output data under the data/ folder

int verbose 		= 1;									//!< Set to zero in order reduce the console output of the code

int abs_error_threshold_fix = 1;							//! Set to one in order to implement the proto-type absolute error thresholding fix

int output_hom 		= 0;									//! Set to 1 to also output the basis of homogeneous solutions

static pthread_t thread_id;

struct orbital_params orbit;

void 	rescale_the_field(struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit);
void 	setup_n_mode_data_structure(struct n_mode_data *n_mode, int l, int m, int n, struct orbital_params *orbit, struct coupled_set *cset);
void 	setup_lm_mode_data_structure(struct lm_mode_data *lm_mode, int l, int m);
void 	free_n_mode_data_structure(struct n_mode_data *n_mode, struct coupled_set *cset, struct orbital_params *orbit);
struct 	lm_mode_data* pick_lm_mode(struct lm_mode_data** lm_modes);
void	free_lm_mode_data_structure(struct lm_mode_data *lm_mode);
int 	gamma_to_l(int gamma);
void 	*counter_routine(void * );
void 	init_counter( MPI_Comm comm, MPI_Comm *counter_comm_p ) ;
void 	counter_nxtval( MPI_Comm counter_comm, int *value );
void 	stop_counter( MPI_Comm *counter_comm_p, MPI_Comm *all_but_zero_comm );
void 	_mkdir(const char *dir); 

typedef void (*function_ptr)(double, struct n_mode_data*, struct orbital_params*);

int main(int argc, char *argv[])
{
    char *outdir, *gridfile;
    if( argv[1] != NULL && argv[2] != NULL && argc == 5){
		p = (double)strtod(argv[1], NULL);
		l_max = (double) strtod(argv[2], NULL);
		gridfile = argv[3];
		outdir = argv[4];
    }else{
		printf("Usage: h1Lorenz <p> <lmax> <gridfile> <outdir>\n");
        exit(0);
    }	
	
	int gamma;
	int myid, numprocs;
	MPI_Comm counter_comm;
	int provided;

	/* MPI Initialization */
	MPI_Init_thread( 0, 0, MPI_THREAD_MULTIPLE, &provided );

	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &myid);

	MPI_Comm all_but_zero_comm;
	int color = 0;
	if (myid == 0) color = MPI_UNDEFINED;
	MPI_Comm_split(MPI_COMM_WORLD, color, myid, &all_but_zero_comm);

	set_primary_orbital_params(p, e);
	get_orbit_params(&orbit, NUM_CHI_VALUES);	
	
	// Load the radial grid
	orbit.gridsize = readin_grid(gridfile, &orbit.grid, &orbit.ra_grid_index, &orbit.r0_grid_index, &orbit.rb_grid_index);

	// override the ra and rb grid indexes to get the homogeneous solutions outputted over the entire radial domain
	orbit.ra_grid_index = 0;
	//orbit.rb_grid_index = 3706;		//FIXME: this only works for the r0=10M case
	//orbit.rb_grid_index = 3526;		//FIXME: this only works for the r0=8M case
	orbit.rb_grid_index = orbit.gridsize - 1;
		
	if(provided == MPI_THREAD_MULTIPLE){
		if(myid == 0) printf("\nMPI multithreading available: multithreading the lm mode counter.\n\n");
		multi_threading = 1;
	}else{
		if(numprocs == 1){
			if(myid == 0) printf("\nMPI multithreading not available, but only one core available so multithreading will be enabled\n\n");
			multi_threading = 1;
		}else{
			if(myid == 0) printf("\nMPI multithreading not available: running the counter on a separate core.\n\n");
			multi_threading = 0;
		}
	}
	multi_threading = 0;


	if(myid == 0){
		printf("Number of cores in use: %d\n", numprocs);
		printf("n-folding the m=0 modes: ");
		if(n_folding){
			printf("enabled\n\n");
		}else{
			printf("disabled\n\n");
		}
		print_orbital_params(&orbit);
	}

	struct coupled_set cset_odd;
	struct coupled_set cset_even;
	struct coupled_set cset_odd_dipole;
	struct coupled_set cset_monopole;
	struct coupled_set cset_even_dipole;
	struct coupled_set cset_even_static;
	
	struct coupled_set *csets[5];
	csets[MONOPOLE] 	= &cset_monopole;
	csets[ODD_DIPOLE] 	= &cset_odd_dipole;
	csets[EVEN_DIPOLE] 	= &cset_even_dipole;
	csets[EVEN] 		= &cset_even;
	csets[ODD]			= &cset_odd;


	cset_even.num_coupled_fields 			= 5;
	cset_even.num_gauge_fields 				= 2;
	cset_even.source 						= malloc(cset_even.num_coupled_fields*sizeof(function_ptr));
	cset_even.source[0] 					= R1_source;
	cset_even.source[1]						= R3_source;
	cset_even.source[2]						= R5_source;
	cset_even.source[3]						= R6_source;
	cset_even.source[4] 					= R7_source; 
	cset_even.out_bc_func 					= &even_out_bcs;
	cset_even.in_bc_func					= &even_in_bcs;
	cset_even.ode_system 					= &even_ode_system;
	cset_even.construct_gauge_fields 		= &construct_R2_R4_gauge_fields;
	cset_even.construct_gauge_fields_on_grid = &construct_R2_R4_gauge_fields_on_grid;

	cset_odd.num_coupled_fields 			= 2;
	cset_odd.num_gauge_fields 				= 1;
	cset_odd.source 						= malloc(cset_odd.num_coupled_fields*sizeof(function_ptr));
	cset_odd.source[0] 						= R9_source;
	cset_odd.source[1] 						= R10_source; 
	cset_odd.out_bc_func 					= &odd_out_bcs;
	cset_odd.in_bc_func						= &odd_in_bcs;
	cset_odd.ode_system 					= &odd_ode_system;
	cset_odd.construct_gauge_fields			= &construct_R8_gauge_field;
	cset_odd.construct_gauge_fields_on_grid	= &construct_R8_gauge_field_on_grid;

	cset_odd_dipole.num_coupled_fields 		= 1;
	cset_odd_dipole.num_gauge_fields		= 1;
	cset_odd_dipole.source 					= malloc(cset_odd_dipole.num_coupled_fields*sizeof(function_ptr));
	cset_odd_dipole.source[0] 				= R9_source;
	cset_odd_dipole.out_bc_func 			= &odd_dipole_out_bcs;
	cset_odd_dipole.in_bc_func				= &odd_dipole_in_bcs;
	cset_odd_dipole.ode_system 				= &odd_dipole_ode_system;
	cset_odd_dipole.construct_gauge_fields	= &construct_R8_gauge_field;
	cset_odd_dipole.construct_gauge_fields_on_grid	= &construct_R8_gauge_field_on_grid;

	cset_monopole.num_coupled_fields 		= 3;
	cset_monopole.num_gauge_fields 			= 1;
	cset_monopole.source					= malloc(cset_monopole.num_coupled_fields*sizeof(function_ptr));
	cset_monopole.source[0]					= R1_source;
	cset_monopole.source[1]					= R3_source;
	cset_monopole.source[2]					= R6_source;
	cset_monopole.out_bc_func				= monopole_out_bcs;
	cset_monopole.in_bc_func				= monopole_in_bcs;
	cset_monopole.ode_system				= monopole_ode_system;
	cset_monopole.construct_gauge_fields 	= &construct_R2_R4_gauge_fields;

	cset_even_dipole.num_coupled_fields 	= 4;
	cset_even_dipole.num_gauge_fields		= 2;
	cset_even_dipole.source					= malloc(cset_even_dipole.num_coupled_fields*sizeof(function_ptr));
	cset_even_dipole.source[0]				= R1_source;
	cset_even_dipole.source[1]				= R3_source;
	cset_even_dipole.source[2]				= R5_source;
	cset_even_dipole.source[3]				= R6_source;
	cset_even_dipole.out_bc_func			= even_dipole_out_bcs;
	cset_even_dipole.in_bc_func				= even_dipole_in_bcs;
	cset_even_dipole.ode_system				= even_dipole_ode_system;
	cset_even_dipole.construct_gauge_fields = &construct_R2_R4_gauge_fields;
	cset_even_dipole.construct_gauge_fields_on_grid = &construct_R2_R4_gauge_fields_on_grid;

	cset_even_static.num_coupled_fields		= 3;
	cset_even_static.num_gauge_fields		= 4;
	cset_even_static.source					= malloc(cset_even_static.num_coupled_fields*sizeof(function_ptr));
	cset_even_static.source[0]				= R1_source;
	cset_even_static.source[1]				= R3_source;
	cset_even_static.source[2]				= R5_source;
	cset_even_static.out_bc_func			= even_static_out_bcs;
	cset_even_static.in_bc_func				= even_static_in_bcs;
	cset_even_static.ode_system				= even_static_ode_system;
	cset_even_static.construct_gauge_fields	= &construct_R6_R7_R2_R4_gauge_fields;
	cset_even_static.construct_gauge_fields_on_grid	= &construct_R6_R7_R2_R4_gauge_fields_on_grid;

	init_counter( MPI_COMM_WORLD, &counter_comm );

	int l, m;

	if(multi_threading == 1 || ( multi_threading == 0 && myid != 0)){

		while(1){
			// Get the next gamma 
			counter_nxtval( counter_comm, &gamma );

			 l = gamma_to_l(gamma);
			 m = gamma - l*(l+1)/2;
			 
			 if(m<m_min) continue;

			if(l > l_max + 0*3) break;	// Need only compute three more tensor spherical harmonic l-modes than the required l_max scalar harmonic l-modes

			printf("Calculating tensor mode: l=%d m=%d\n", l, m);		
		
			struct lm_mode_data lm_mode;
			lm_mode.cset = csets[select_correct_coupled_set(l, m)];
			setup_lm_mode_data_structure(&lm_mode, l, m);

			int n = 0;

			int convergence = 0;
			struct coupled_set* cset;

			cset = lm_mode.cset;

			convergence = 0;
			n = n_min;
			while(convergence != CONVERGED){

				struct n_mode_data n_mode;
				setup_n_mode_data_structure(&n_mode, l, m, n, &orbit, cset);

				if(e == 0) {
					if(abs(n) > 0) break;
				}

				if(l == 0 && m == 0 && n == 0) {														// Static monopole requires a different treatment
					construct_static_monopole(&n_mode, &orbit);														
				}else{
				
					if((l+m) % 2 == 0 && m == 0 && n == 0) cset = &cset_even_static;					// Even static modes require a different treatment
					else cset = csets[select_correct_coupled_set(l,m)];	

					integrate_field_equations(&orbit, &n_mode);

					calculate_scaling_coefficients(&n_mode, &orbit);

					rescale_the_field(&lm_mode, &n_mode, &orbit);
				
				}
												
				// Do not perform any convergence testing (this is a relic of the eccentric orbit piece of the code)
				convergence = 1;
				
				if(l >= 0){
					hid_t       file_id; // identifiers 
					hsize_t     dims[2];
					int i,j,k;
					double* data;
					double* data2;
					
					int nf = cset->num_coupled_fields;

					// Create a new directory for the data and open the HDF5 file using default properties. 
					char fileloc[50], filename[50];
					mkdir(outdir, 0700);
					sprintf(filename, "h1-l%dm%d.h5", l, m);
					strcpy(fileloc, outdir);
					strcat(fileloc, "/");
					strcat(fileloc, filename);

					file_id = H5Fcreate(fileloc, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

					// Output the grid
					dims[0] = orbit.gridsize;
					H5LTmake_dataset_double(file_id, "grid", 1, dims, orbit.grid);
					H5LTset_attribute_int(file_id, "grid", "raindex", &orbit.ra_grid_index, 1);
					H5LTset_attribute_int(file_id, "grid", "r0index", &orbit.r0_grid_index, 1);
					H5LTset_attribute_int(file_id, "grid", "rbindex", &orbit.rb_grid_index, 1);
					
					// Output the homogeneous fields
					if(output_hom){
						dims[0] = orbit.rb_grid_index + 1;
						dims[1] = 4*nf*nf;

						data = malloc(sizeof(double) * (4*nf*nf)*(orbit.rb_grid_index + 1));

						for(i = 0; i <= orbit.rb_grid_index; i++){
							for(j = 0; j < nf; j++){
								for(k = 0; k < nf; k++){
									data[i*(4*nf*nf) + (j*4*nf) + k*4 + 0] = n_mode.hom_data_in[j][k][0][i];
									data[i*(4*nf*nf) + (j*4*nf) + k*4 + 1] = n_mode.hom_data_in[j][k][1][i];
									data[i*(4*nf*nf) + (j*4*nf) + k*4 + 2] = n_mode.hom_data_in[j][k][2][i];
									data[i*(4*nf*nf) + (j*4*nf) + k*4 + 3] = n_mode.hom_data_in[j][k][3][i];
								}	
							}
						}

						H5LTmake_dataset_double(file_id, "hom_in", 2, dims, data);
						free(data);

						dims[0] = orbit.gridsize - orbit.ra_grid_index;
						data = malloc(sizeof(double) * (4*nf*nf)*(orbit.gridsize - orbit.ra_grid_index));

						for(i = 0; i < orbit.gridsize - orbit.ra_grid_index; i++){
							for(j = 0; j < nf; j++){
								for(k = 0; k < nf; k++){
									data[i*(4*nf*nf) + (j*4*nf) + k*4 + 0] = n_mode.hom_data_out[j][k][0][i];
									data[i*(4*nf*nf) + (j*4*nf) + k*4 + 1] = n_mode.hom_data_out[j][k][1][i];
									data[i*(4*nf*nf) + (j*4*nf) + k*4 + 2] = n_mode.hom_data_out[j][k][2][i];
									data[i*(4*nf*nf) + (j*4*nf) + k*4 + 3] = n_mode.hom_data_out[j][k][3][i];
								}	
							}
						}

						H5LTmake_dataset_double(file_id, "hom_out", 2, dims, data);
						free(data);
					}


					// Add homogeneous (pure gauge) hbar9 to the odd dipole to make it regular at the horizon (but keep the correct ang. mom.)
					if(l==1 && m ==0){
						for(i = 0; i <= orbit.gridsize; i++){
							double A = -256.0/3.0*sqrt(3.0*M_PI)*orbit.E*orbit.Omega_phi/(orbit.p - 2.0);
							double r = orbit.grid[i];
							if(i==orbit.gridsize) r = orbit.p;
							n_mode.inhom_data[0][0][i] 	+= A/(r*r);
							n_mode.inhom_data[0][2][i] 	+= -2.0*A/(r*r*r);
						}
					}
										
					// Output the inhomogenous data to the left of the particle
					nf = cset->num_coupled_fields + cset->num_gauge_fields;
					dims[0] = orbit.r0_grid_index+1;
					dims[1] = 4*nf;
					data = malloc(sizeof(double) * 4*nf*dims[0]);
					for(i = 0; i < dims[0]; i++){					
						for(j = 0; j < nf; j++){							
							data[i*4*nf + j*4 + 0] = n_mode.inhom_data[j][0][i];
							data[i*4*nf + j*4 + 1] = n_mode.inhom_data[j][1][i];
							data[i*4*nf + j*4 + 2] = n_mode.inhom_data[j][2][i];
							data[i*4*nf + j*4 + 3] = n_mode.inhom_data[j][3][i];					
						}
					}

					H5LTmake_dataset_double(file_id, "inhom_left", 2, dims, data);
					free(data);
			
			
					// Output the inhomogenous data to the right of the particle
					nf = cset->num_coupled_fields + cset->num_gauge_fields;
					dims[0] = orbit.gridsize - orbit.r0_grid_index;
					dims[1] = 4*nf;
					data = malloc(sizeof(double) * 4*nf*dims[0]);
					
					for(i = orbit.r0_grid_index; i < orbit.gridsize + 1; i++){
						for(j = 0; j < nf; j++){
							if(i == orbit.gridsize){ // Use the right-hand derviatives at the particle
								data[0*4*nf + j*4 + 0] = n_mode.inhom_data[j][0][orbit.gridsize];
								data[0*4*nf + j*4 + 1] = n_mode.inhom_data[j][1][orbit.gridsize];
								data[0*4*nf + j*4 + 2] = n_mode.inhom_data[j][2][orbit.gridsize];
								data[0*4*nf + j*4 + 3] = n_mode.inhom_data[j][3][orbit.gridsize];
							}else{
								data[(i-orbit.r0_grid_index)*4*nf + j*4 + 0] = n_mode.inhom_data[j][0][i];
								data[(i-orbit.r0_grid_index)*4*nf + j*4 + 1] = n_mode.inhom_data[j][1][i];
								data[(i-orbit.r0_grid_index)*4*nf + j*4 + 2] = n_mode.inhom_data[j][2][i];
								data[(i-orbit.r0_grid_index)*4*nf + j*4 + 3] = n_mode.inhom_data[j][3][i];	
							}				
						}
					}

					H5LTmake_dataset_double(file_id, "inhom_right", 2, dims, data);
					free(data);
					
					
					//Output the infinity complex amplitudes of the asymptotic metric perturbation
					nf = cset->num_coupled_fields;		// FIXME need to add gauge fields
					dims[0] = 2*nf;						// times 2 for the real and imaginary parts
					data  = calloc(dims[0], sizeof(double));
					data2 = calloc(dims[0], sizeof(double));
					for(i = 0; i < nf; i++){
						for(j = 0; j < nf; j++){
							data[2*i] 	+= creal(n_mode.C_out[j] * n_mode.C_out_hom[i][j]);
							data[2*i+1] += cimag(n_mode.C_out[j] * n_mode.C_out_hom[i][j]);
							
							data2[2*i] 	 += creal(n_mode.C_in[j] * n_mode.C_in_hom[i][j]);
							data2[2*i+1] += cimag(n_mode.C_in[j] * n_mode.C_in_hom[i][j]);
						}						
					}
					
					H5LTmake_dataset_double(file_id, "C_inf", 1, dims, data);
					H5LTmake_dataset_double(file_id, "C_horiz", 1, dims, data2);
					
					free(data);
					free(data2);	
			
					// Close the file.
					H5Fclose(file_id);
					
					

				}

				free_n_mode_data_structure(&n_mode, cset, &orbit);

			}

			free_lm_mode_data_structure(&lm_mode);
		}

		stop_counter( &counter_comm, &all_but_zero_comm );					// Exit the counter thread/core as it is no longer needed

	}

	

	/* MPI Termination */
	MPI_Finalize();

	return 0;
}

/** Rescales the homogeneous fields as perscribed by the method of Extended Homogeneous Solutions
*
* @param lm_mode the lm mode that the rescaled n mode fields will be added to
* @param n_mode the n mode whose fields are to be rescaled
* @param orbit the orbital paramters
*/
void rescale_the_field(struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit)
{
	struct coupled_set *cset = n_mode->cset;

	double l = n_mode->l;
	double m = n_mode->m;
	double n = n_mode->n;

	int i, k, d;								// k is the field index, d the basis index and i the orbit index
	for(i = 0; i < 2*NUM_CHI_VALUES+1; i++){

		int r_chi_index = i;
		if(i > NUM_CHI_VALUES) r_chi_index = 2*NUM_CHI_VALUES - i;

		for(k = 0; k < cset->num_coupled_fields; k++){
			for(d = 0; d < cset->num_coupled_fields; d++){				
				
				//if(i==0) printf("%.12e %.12e\n", creal(n_mode->R_in[k][d][0]), creal(n_mode->R_out[k][d][0]));

				lm_mode->h_in[k][i] 		+= n_fold_function(m, n, n_mode->C_in[d]  * n_mode->R_in[k][d][r_chi_index]*cexp(-I*n_mode->omega * orbit->t_p[i]));
				lm_mode->h_out[k][i] 		+= n_fold_function(m, n, n_mode->C_out[d] * n_mode->R_out[k][d][r_chi_index]*cexp(-I*n_mode->omega * orbit->t_p[i]));

				lm_mode->h_in_rs_deriv[k][i]	+= n_fold_function(m, n, n_mode->C_in[d]  * n_mode->R_in_rs_deriv[k][d][r_chi_index]*cexp(-I*n_mode->omega * orbit->t_p[i]));
				lm_mode->h_out_rs_deriv[k][i] 	+= n_fold_function(m, n, n_mode->C_out[d] * n_mode->R_out_rs_deriv[k][d][r_chi_index]*cexp(-I*n_mode->omega * orbit->t_p[i]));	
			
				lm_mode->h_in_t_deriv[k][i]		+= n_fold_function(m, n, -I*n_mode->omega * n_mode->C_in[d] * n_mode->R_in[k][d][r_chi_index]*cexp(-I*n_mode->omega * orbit->t_p[i]));
				lm_mode->h_out_t_deriv[k][i]	+= n_fold_function(m, n, -I*n_mode->omega * n_mode->C_out[d] * n_mode->R_out[k][d][r_chi_index]*cexp(-I*n_mode->omega * orbit->t_p[i]));
			}
		}
	}
	
	//if(n_mode->l==2) printf("test hom: %e\n", creal(n_mode->R_in[0][0][0]));

	// store the inhomogeneous field everywhere
	int j;
	for(i = 0; i < orbit->gridsize+1; i++){
		//double r = orbit->grid[i];
		for(k = 0; k < cset->num_coupled_fields; k++){
			double complex R = 0;
			double complex dR = 0;
			for(d = 0; d < cset->num_coupled_fields; d++){	
				double**** hom_data; double complex* C;
				if(i <= orbit->r0_grid_index) {hom_data = n_mode->hom_data_in; j = i; C = &n_mode->C_in[d]; }
				else { hom_data = n_mode->hom_data_out; j = orbit->gridsize - 1 - i; C = &n_mode->C_out[d]; }
				
				// Store the data for the right-derivate at the end of the inhom array
				if(i == orbit->gridsize){ hom_data = n_mode->hom_data_out; j = orbit->gridsize - 1 - orbit->r0_grid_index; C = &n_mode->C_out[d];}
				
				R 	+= *C * (hom_data[k][d][0][j] + I*hom_data[k][d][1][j]);
				dR 	+= *C * (hom_data[k][d][2][j] + I*hom_data[k][d][3][j]);
				

	
				n_mode->inhom_data[k][0][i] = creal(R);
				n_mode->inhom_data[k][1][i] = cimag(R);
				n_mode->inhom_data[k][2][i] = creal(dR);
				n_mode->inhom_data[k][3][i] = cimag(dR);
			}
		}
	}

	// Also construct the gauge fields
	cset->construct_gauge_fields(lm_mode, n_mode, orbit, cset->num_coupled_fields);
	if(lm_mode->l > 0){
		cset->construct_gauge_fields_on_grid(n_mode->inhom_data, n_mode, orbit, cset->num_coupled_fields);
	}

}


/** Given gamma= l(l+1)/2 +m return the corresponding l
* 
* @param gamma gamma
* @return spherical harmonic l index
*/
int gamma_to_l(int gamma)
{
	int l = 0;
	while(1){
		if( l*(l+1)/2 > gamma){
			return l-1;
		}
		l++;
	} 
} 

/** Returns gamma = l(l+1)/2 +m. There is a 1-1 mapping between gamma and an lm mode.
*
* @param l the spherical harmonic index l
* @param m the spherical harmonic index m
* @return gamma = l(l+1)/2 + m
*/
int lm_to_gamma(int l, int m)
{
	return l*(l+1)/2 + m;
}	

void *counter_routine( void *p )
{
MPI_Comm *counter_comm_p = (MPI_Comm *)p;

	int incr;
	int ival = lm_to_gamma(l_min, m_min);
	MPI_Status status;

	while (1) {
		MPI_Recv( &incr, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, *counter_comm_p, &status );
		if(status.MPI_TAG == 1) return NULL;
		MPI_Send( &ival, 1, MPI_INT, status.MPI_SOURCE, 0, *counter_comm_p );
		 ival += incr;
	}
}


/* Only process with rank zero creates the thread in this function */
void init_counter( MPI_Comm comm, MPI_Comm *counter_comm_p ) {

	int rank;

	MPI_Comm_dup ( comm, counter_comm_p );
	MPI_Comm_rank( comm, &rank );
	if ( rank == 0 ) { /* only rank 0 */
		if(multi_threading) pthread_create( &thread_id, NULL, counter_routine, counter_comm_p );
		else counter_routine(counter_comm_p);
	}
}
/* Any process can call this to fetch and increment by value */
void counter_nxtval( MPI_Comm counter_comm, int *value )
{
	int incr = 1;

	MPI_Send( &incr, 1, MPI_INT, 0, 0, counter_comm );
	MPI_Recv( value, 1, MPI_INT, 0, 0, counter_comm, MPI_STATUS_IGNORE );
}


/* Stops the mode counter. Every process in counter_comm (including rank 0!) must call stop_counter in order to exit the thread. */
void stop_counter( MPI_Comm *counter_comm_p, MPI_Comm *all_but_zero_comm )
{
	int rank;
	MPI_Comm_rank( *counter_comm_p, &rank );

	if(multi_threading)	MPI_Barrier( *counter_comm_p );
	else MPI_Barrier(*all_but_zero_comm);

	if (rank == 0 || (multi_threading == 0 && rank == 1)) {
		MPI_Send( MPI_BOTTOM, 0, MPI_INT, 0, 1, *counter_comm_p );
		if(multi_threading) pthread_join( thread_id, NULL );
	}
	MPI_Comm_free( counter_comm_p );
}

int select_correct_coupled_set(int l, int m)
{
	if((l+m) % 2 == 0) {							// Even sector
		if(l == 0){
			return MONOPOLE;
		}else if(l == 1){
			return EVEN_DIPOLE;
		}else{
			return EVEN;
		}
	}else{											// Odd sector
		if(l == 1){									
			return ODD_DIPOLE;
		}else{
			return ODD;
		}
	}
}

/** Allocates memory for and sets up the suppied n_mode_data struct.
 *
 * @param l the spherical harmonic l index of the n mode
 * @param m the spherical harminic m index of the n mode
 */
void setup_lm_mode_data_structure(struct lm_mode_data *lm_mode, int l, int m)
{
	struct coupled_set *cset = lm_mode->cset;

	lm_mode->l = l;
	lm_mode->m = m;
	
	lm_mode->h_out 			= alloc_2D_complex_array(cset->num_coupled_fields + cset->num_gauge_fields, 2*NUM_CHI_VALUES+1);
	lm_mode->h_in	 		= alloc_2D_complex_array(cset->num_coupled_fields + cset->num_gauge_fields, 2*NUM_CHI_VALUES+1);

	lm_mode->h_out_rs_deriv = alloc_2D_complex_array(cset->num_coupled_fields + cset->num_gauge_fields, 2*NUM_CHI_VALUES+1);
	lm_mode->h_in_rs_deriv	= alloc_2D_complex_array(cset->num_coupled_fields + cset->num_gauge_fields, 2*NUM_CHI_VALUES+1);

	lm_mode->h_out_t_deriv	= alloc_2D_complex_array(cset->num_coupled_fields + cset->num_gauge_fields, 2*NUM_CHI_VALUES+1);
	lm_mode->h_in_t_deriv	= alloc_2D_complex_array(cset->num_coupled_fields + cset->num_gauge_fields, 2*NUM_CHI_VALUES+1);

}

/** Frees the memory associated with a given lm mode
*
* @param lm_mode the lm_mode to be freed
*/
void free_lm_mode_data_structure(struct lm_mode_data *lm_mode)
{
	struct coupled_set *cset = lm_mode->cset;

	free_alloced_2D_complex_array(lm_mode->h_out, cset->num_coupled_fields + cset->num_gauge_fields);
	free_alloced_2D_complex_array(lm_mode->h_in, cset->num_coupled_fields + cset->num_gauge_fields);

	free_alloced_2D_complex_array(lm_mode->h_out_rs_deriv, cset->num_coupled_fields + cset->num_gauge_fields);
	free_alloced_2D_complex_array(lm_mode->h_in_rs_deriv, cset->num_coupled_fields + cset->num_gauge_fields);
	
	free_alloced_2D_complex_array(lm_mode->h_out_t_deriv, cset->num_coupled_fields + cset->num_gauge_fields);
	free_alloced_2D_complex_array(lm_mode->h_in_t_deriv, cset->num_coupled_fields + cset->num_gauge_fields);
}

/** Allocates memory for and sets up the suppied n_mode_data struct.
 *
 * @param l the spherical harmonic l index of the n mode
 * @param m the spherical harminic m index of the n mode
 * @param n the frequency n index of the n mode
 * @param orbit the orbital parameters
 * @param cset the coupled set corresponding to this lmn mode
 */
void setup_n_mode_data_structure(struct n_mode_data *n_mode, int l, int m, int n, struct orbital_params *orbit, struct coupled_set *cset)
{
	int nf 			= cset->num_coupled_fields;
	
	n_mode->l 		= l;
	n_mode->m 		= m;
	n_mode->n 		= n;
	n_mode->omega 	= m*orbit->Omega_phi + n*orbit->Omega_r;

	n_mode->C_in 			= (double complex*)calloc(cset->num_coupled_fields, sizeof(double complex));
	n_mode->C_out 			= (double complex*)calloc(cset->num_coupled_fields, sizeof(double complex));
	
	n_mode->C_in_hom 		= alloc_2D_complex_array(nf,nf);
	n_mode->C_out_hom		= alloc_2D_complex_array(nf,nf);

	n_mode->R_out 			= alloc_3D_complex_array(cset->num_coupled_fields, cset->num_coupled_fields, NUM_CHI_VALUES+1);
	n_mode->R_out_rs_deriv 	= alloc_3D_complex_array(cset->num_coupled_fields, cset->num_coupled_fields, NUM_CHI_VALUES+1);
	n_mode->R_in			= alloc_3D_complex_array(cset->num_coupled_fields, cset->num_coupled_fields, NUM_CHI_VALUES+1);
	n_mode->R_in_rs_deriv 	= alloc_3D_complex_array(cset->num_coupled_fields, cset->num_coupled_fields, NUM_CHI_VALUES+1);

	n_mode->cset = cset;

	n_mode->hom_data_in		= alloc_4D_double_array(cset->num_coupled_fields, cset->num_coupled_fields, 4, orbit->rb_grid_index + 1);
	n_mode->hom_data_out	= alloc_4D_double_array(cset->num_coupled_fields, cset->num_coupled_fields, 4, orbit->gridsize - orbit->ra_grid_index );

	n_mode->inhom_data		= alloc_3D_double_array(cset->num_coupled_fields + cset->num_gauge_fields, 4, orbit->gridsize + 1);		// +1 as we need to store the left- and right-derivative at the particle (the latter is stored as the last entry in the array)
	n_mode->hR				= alloc_3D_double_array(cset->num_coupled_fields + cset->num_gauge_fields, 4, orbit->gridsize +1 );
	
	n_mode->src 					= alloc_2D_complex_array(nf, orbit->gridsize+1);				// plus one as there are two values at r0 (the latter is stored as the last entry in the array)
	n_mode->scaling_coeffs_re_lhs	= alloc_2D_double_array(2*nf, orbit->r0_grid_index + 1);
	n_mode->scaling_coeffs_im_lhs	= alloc_2D_double_array(2*nf, orbit->r0_grid_index + 1);
	n_mode->scaling_coeffs_re_rhs	= alloc_2D_double_array(2*nf, orbit->gridsize - orbit->r0_grid_index + 1);
	n_mode->scaling_coeffs_im_rhs	= alloc_2D_double_array(2*nf, orbit->gridsize - orbit->r0_grid_index + 1);

}

/** Frees the memory associated with a given n mode
*
* @param n_mode the n_mode to be freed
* @param cset the coupled set associated to this n mode
*/
void free_n_mode_data_structure(struct n_mode_data *n_mode, struct coupled_set *cset, struct orbital_params *orbit)
{
	int nf 			= cset->num_coupled_fields;
	
	free_alloced_3D_complex_array(n_mode->R_out, cset->num_coupled_fields, cset->num_coupled_fields);
	free_alloced_3D_complex_array(n_mode->R_in, cset->num_coupled_fields, cset->num_coupled_fields);
	free_alloced_3D_complex_array(n_mode->R_out_rs_deriv, cset->num_coupled_fields, cset->num_coupled_fields);
	free_alloced_3D_complex_array(n_mode->R_in_rs_deriv, cset->num_coupled_fields, cset->num_coupled_fields);

	free_alloced_4D_double_array(n_mode->hom_data_in, cset->num_coupled_fields, cset->num_coupled_fields, 4);
	free_alloced_4D_double_array(n_mode->hom_data_out, cset->num_coupled_fields, cset->num_coupled_fields, 4);

	free_alloced_3D_double_array(n_mode->inhom_data, cset->num_coupled_fields + cset->num_gauge_fields, 4);
	free_alloced_3D_double_array(n_mode->hR, cset->num_coupled_fields + cset->num_gauge_fields, 4);
	
	free_alloced_2D_complex_array(n_mode->src, nf);
	free_alloced_2D_double_array(n_mode->scaling_coeffs_re_lhs, 2*nf);
	free_alloced_2D_double_array(n_mode->scaling_coeffs_im_lhs, 2*nf);
	free_alloced_2D_double_array(n_mode->scaling_coeffs_re_rhs, 2*nf);
	free_alloced_2D_double_array(n_mode->scaling_coeffs_im_rhs, 2*nf);

	free(n_mode->C_out);
	free(n_mode->C_in);
}


double complex n_fold_function(double m, double n, double complex input)
{
	if(m == 0 && n_folding){
		if(n == 0) return input;
		else return 2.0*creal(input);
	}else{
		return input;
	}
}

void nprintf( const char* format, ... ) {
	if(verbose){
		va_list args;
		va_start( args, format );
		vprintf( format, args );
		va_end( args );
	}
}


