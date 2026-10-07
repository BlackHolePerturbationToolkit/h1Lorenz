#pragma once
#include <complex.h>
#include "GSF_ecc.h"

/* number of terms of the horizon expansion written to each output file */
#define HORIZON_SERIES_TERMS 12

/* Coefficients c^i_n (out[i*K + n]) of h^i = exp(-i omega r*) sum_n c^i_n (r - 2M)^n for the
 * output fields of one mode (coupled then gauge fields); returns the number of fields, or 0
 * for the even static modes, which h1Lorenz does not compute. */
int horizon_series(struct coupled_set *cset, struct n_mode_data *n_mode, struct orbital_params *orbit, int K, double complex *out);

/* Coefficients b^i_n, c^i_n (out_b/out_c[i*K + n]) of h^i = exp(i omega r*) sum_n (b^i_n + c^i_n log r)/r^n
 * for the output fields of one mode; returns the number of fields, or 0 for the even static modes. */
int infinity_series(struct coupled_set *cset, struct n_mode_data *n_mode, struct orbital_params *orbit, int K,
                    double complex *out_b, double complex *out_c);

