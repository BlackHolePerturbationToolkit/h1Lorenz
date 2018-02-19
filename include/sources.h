#include "GSF_ecc.h"
#include <complex.h>

#define SOURCE_TOTAL 0
#define SOURCE_COEFF 1
#define SOURCE_OSCIL 2

double complex R1_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece);
double complex R3_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece);
double complex R5_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece);
double complex R6_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece);
double complex R7_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece);

double complex R9_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece);
double complex R10_source(double chi, struct n_mode_data *n_mode, struct orbital_params *orbit, int source_piece);
