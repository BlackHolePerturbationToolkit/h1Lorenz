#include <complex.h>
extern double M;


double orbit_radial_period();
double orbit_radial_freq();
double orbit_phi_freq();
double orbit_ang_freq();
double orbit_ang_mom();
double orbit_energy();
double orbit_contra_ut(double r);
double orbit_contra_ur(double chi);
double orbit_covariant_ut(double r);
double orbit_covariant_ur(double chi);
double r_max();
double r_min();
double r_of_chi(double chi);
double t_of_chi(double chi);
double phi_of_chi(double chi);
double rs_of_r(double r);
double r_of_rs(double rs);
double f(double r);
double dt_dchi(double chi);

double R8_static_out_solution(int l, double r_p);
double R8_static_out_rs_deriv_solution(int l, double r_p);
double R8_static_in_solution(int l, double r_p);
double R8_static_in_rs_deriv_solution(int l, double r_p);

double locate_effective_horizon();


/** \brief Stores the orbital paramters */
struct orbital_params{
	double e;						//!< The orbital eccentricity
	double p;						//!< The semi-latus rectum
	double T_r;						//!< The radial period 
	double Omega_r;					//!< The radial frequency
	double Omega_phi;				//!< The azimuthal frequency
	double E;						//!< The orbital energy
	double L;						//!< The orbital angular momentum
	double r_max;					//!< The maximum Schwarzschild radius attained by the orbit
	double rs_max;					//!< The maximum tortoise(r*) radius attained by the orbit	
	double r_min;					//!< The minimum Schwarzschild radius attained by the orbit
	double rs_min;					//!< The minimum tortoise(r*) radius attained by the orbit

	double* chi_p;					//!< The array of chi values along the orbit
	double* r_p;					//!< The array of Schwarzschild radius values along the orbit
	double* rs_p;					//!< The array of tortoise(r*) radius values along the orbit
	double* t_p;					//!< The Schwarzschild t coordinate values along the orbit
	double* phi_p;					//!< The Schwarzschild phi coordinate values along the orbit
	
	double* grid;
	int gridsize;
	
	double complex*	hP_gauge;		//!< Array for storing the contributions to the gauge fields from the punctures. Indexed by 2*gamma(l,m)

	int r0_grid_index;
	int ra_grid_index;
	int rb_grid_index;

};


void	set_primary_orbital_params(double p_new, double e_new);
void	get_orbit_params(struct orbital_params *orbit, int num_chi_values);
void 	print_orbital_params(struct orbital_params *orbit);
double 	get_orbital_e();
double 	get_orbital_p();

double circular_orbit_phi_freq();
double circular_orbit_energy();
double circular_orbit_ang_mom();


