#define NUM_CHI_VALUES 10

#define MONOPOLE 0
#define ODD_DIPOLE 1
#define EVEN_DIPOLE 2
#define EVEN 3
#define ODD 4

#define INNER_FIELDS 0
#define OUTER_FIELDS 1

#ifndef __test
#define __test

#include "boundary_conditions.h"
#include "libSch.h"
#include "math_funcs.h"
#include "array_handling.h"
#include "SF_data.h"



struct coupled_set;

/** Stores the information related to a given lm mode */
struct lm_mode_data;
struct n_mode_data;
struct fields_data;

/** \brief Stores the information related to a given lm mode */
struct lm_mode_data
{
	double l;														//!< The spherical harmonic l index
	double m;														//!< The spherical harmonic m index

	double complex** h_out;												//!< The h_out field. The first index stores denotes the field number and the second the position along the orbit
	double complex** h_in;													//!< The h_in field. The first index stores denotes the field number and the second the position along the orbit
	double complex** h_out_rs_deriv;										//!< The h_out field. The first index stores denotes the field number and the second the position along the orbit
	double complex** h_in_rs_deriv;										//!< The r* derivative of the h_in field. The first index stores denotes the field number and the second the position along the orbit
	double complex** h_out_t_deriv;										//!< The r* derivative of the h_out field. The first index stores denotes the field number and the second the position along the orbit
	double complex** h_in_t_deriv;											//!< The t derivaitive of the h_in field. The first index stores denotes the field number and the second the position along the orbit

	double *** inhom_data;											//!< The inhomogenous field data on the inhom grid, 1 = grid index, 2 = field index, 3 = part index

	struct coupled_set* cset;										//!< The coupled set corresponding to this lm mode
};

/** \brief Stores the information related to a given n mode */
struct n_mode_data
{
	double l;														//!< The spherical harmonic l index
	double m;														//!< The spherical harmonic m index
	double n;														//!< The frequeny n index
	double omega;													//!< The mode frequency

	double complex*** R_out;										//!< The R_out fields. The first index is the field index, the second the basis index, the thirds the chi index.
	double complex*** R_out_rs_deriv;								//!< The r*-derivatives of the R_out fields. The first index is the field index, the second the basis index, the thirds the chi index.
	double complex*** R_in;											//!< The R_in fields. The first index is the field index, the second the basis index, the thirds the chi index.
	double complex*** R_in_rs_deriv;								//!< The r*-derivatives of the R_out fields. The first index is the field index, the second the basis index, the thirds the chi index.

	double complex* C_in;											//!< Array of scaling coefficients for the inner fields. The index is the basis number.
	double complex* C_out;											//!< Array of scaling coefficients for the outer fields. The index is the basis number.

	double complex** C_in_hom;
	double complex** C_out_hom;


	double**** hom_data_in;											//!< Array storing the homogeneous field data,  1 = field index, 2 = basis index, 3 = part (re/im/dre/dim) index, 4 = grid index
	double**** hom_data_out;										//!< Array storing the homogeneous field data,  1 = field index, 2 = basis index, 3 = part (re/im/dre/dim) index, 4 = grid index
	double***  inhom_data;											//!< Array storing the inhomogeneous field data,  1 = field index, 2 = part (re/im/dre/dim) index, 3 = grid index
	
	double***  hR;													//!< Array storing the regularized field data,  1 = field index, 2 = part (re/im/dre/dim) index, 3 = grid index
	
	double complex** src;											//!< Array for storing the source, 1 = field index, 2 = grid index
		
	double** scaling_coeffs_re_lhs;									//!< Array for storing the real part of the scaling coefficients, 1 = coeff index, 2 = grid index
	double** scaling_coeffs_im_lhs;									//!< Array for storing the imaginary part of the scaling coefficients, 1 = coeff index, 2 = grid index
	double** scaling_coeffs_re_rhs;									//!< Array for storing the real part of the scaling coefficients, 1 = coeff index, 2 = grid index
	double** scaling_coeffs_im_rhs;									//!< Array for storing the imaginary part of the scaling coefficients, 1 = coeff index, 2 = grid index

	struct coupled_set* cset;										//!< The coupled set associated with this n mode

};


/** \brief Stores the information related to a given coupled set (coupled odes, their bcs functions etc) */
struct coupled_set
{
	int 	num_coupled_fields;																					//!< The number of coupled fields in this coupled set
	int		num_gauge_fields;																					//!< The number of gauge fields to be algebraically constructed from the integrated fields
	void 	(*out_bc_func)(double[], struct boundary_conditions*);												//!< Pointer to the outer boundary condition function
	void 	(*in_bc_func)(double[], struct boundary_conditions*);												//!< Pointer to the inner boundary condition function
	double complex	(**source)(double, struct n_mode_data*, struct orbital_params*, int source_piece);					//!< Array of pointers to the FD source functions
	int 	(*ode_system) (double, const double[], double[],  void* );											//!< The set of ODEs for this coupled set
	void 	(*construct_gauge_fields)(struct lm_mode_data*, struct n_mode_data*, struct orbital_params*, int);	//!< Pointer to the function that algebraically constructs the gauge fields
	void 	(*construct_inner_fields)(struct lm_mode_data*, struct fields_data*, int);							//!< Pointer to the function that maps the inner fields in the coupled set to the usual numbering
	void 	(*construct_outer_fields)(struct lm_mode_data*, struct fields_data*, int);							//!< Pointer to the function that maps the outer fields in the coupled set to the usual numbering

	void 	(*construct_gauge_fields_on_grid)(double ***data, struct n_mode_data*, struct orbital_params*, int);	//!< Pointer to the function that algebraically constructs the gauge fields on the grid
};



/** \brief Stores the information relating to a given lm mode, in terms of the usual numbering system*/
struct fields_data {
	double m;										//!< The spherical harmonic m index
	double r;										//!< The Schwarzschild radius the value of the fields and their derivatives in this structure are stored at
	double eff;										//!< f(r) for the above radius
	double contra_ur;								//!< The value of u^r (r-component of the contravariant 4-velocity) at the above radius
	double contra_uphi;								//!< The value of u^phi (phi-component of the contravariant 4-velocity) at the above radius
	double contra_ut;								//!< The value of u^t (t-component of the contravariant 4-velocity) at the above radius

	double complex h1;										//!< The value of the h1 field at the specified radius
	double complex h2;										//!< The value of the h2 field at the specified radius
	double complex h3;										//!< The value of the h3 field at the specified radius
	double complex h4;										//!< The value of the h4 field at the specified radius
	double complex h5;										//!< The value of the h5 field at the specified radius
	double complex h6;										//!< The value of the h6 field at the specified radius
	double complex h7;										//!< The value of the h7 field at the specified radius
	double complex h8;										//!< The value of the h8 field at the specified radius
	double complex h9;										//!< The value of the h9 field at the specified radius
	double complex h10;									//!< The value of the h10 field at the specified radius

	double complex h1_t_deriv;								//!< The value of the t-derivative of h1 field at the specified radius
	double complex h1_rs_deriv;							//!< The value of the r*-derivative of h1 field at the specified radius
	double complex h2_t_deriv;								//!< The value of the t-derivative of h2 field at the specified radius
	double complex h2_rs_deriv;							//!< The value of the r*-derivative of h2 field at the specified radius
	double complex h3_t_deriv;								//!< The value of the t-derivative of h3 field at the specified radius
	double complex h3_rs_deriv;							//!< The value of the r*-derivative of h3 field at the specified radius
	double complex h4_t_deriv;								//!< The value of the t-derivative of h4 field at the specified radius
	double complex h4_rs_deriv;							//!< The value of the r*-derivative of h4 field at the specified radius
	double complex h5_t_deriv;								//!< The value of the t-derivative of h5 field at the specified radius
	double complex h5_rs_deriv;							//!< The value of the r*-derivative of h5 field at the specified radius
	double complex h6_t_deriv;								//!< The value of the t-derivative of h6 field at the specified radius
	double complex h6_rs_deriv;							//!< The value of the r*-derivative of h6 field at the specified radius
	double complex h7_t_deriv;								//!< The value of the t-derivative of h7 field at the specified radius
	double complex h7_rs_deriv;							//!< The value of the r*-derivative of h7 field at the specified radius
	double complex h8_t_deriv;								//!< The value of the t-derivative of h8 field at the specified radius
	double complex h8_rs_deriv;							//!< The value of the r*-derivative of h8 field at the specified radius
	double complex h9_t_deriv;								//!< The value of the t-derivative of h9 field at the specified radius
	double complex h9_rs_deriv;							//!< The value of the r*-derivative of h9 field at the specified radius
	double complex h10_t_deriv;							//!< The value of the t-derivative of h10 field at the specified radius
	double complex h10_rs_deriv;							//!< The value of the r*-derivative of h10 field at the specified radius
};



/** Define our own printf so we can direct the output of the code as we wish */
void nprintf( const char* format, ... );

/** \brief Given spherical harmonic indicies l and m this returns the correct coupled set */
int		select_correct_coupled_set(int l, int m);

double complex n_fold_function(double m, double n, double complex input);

#endif

/** \mainpage  Gravitational Self Force in Schwarzschild spacetime via the Frequency Domain
 *      \image html orbit-blkhole.jpg
 *	C software to compute the gravitational Self Force in Schwarzschild spacetime. The code makes use of MPI and can therefore be run on a cluster.<br>
 * By N. Warburton and S. Akcay
 */




