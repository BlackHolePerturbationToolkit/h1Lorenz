#include "GSF_ecc.h"

#define INNER_VS_OUTER 0
#define OUTER_CONTRIBUTIONS 1

#define NOT_CONVERGED 0
#define CONVERGED 1
#define N_MIN_CONVERGED 2

void test_for_convergence(int *converged, struct lm_mode_data *lm_mode, struct n_mode_data *n_mode, struct orbital_params *orbit, int abs_n_max);
