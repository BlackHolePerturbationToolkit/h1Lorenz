/*
 * Expansion of the retarded field of one (l, m) mode about the horizon,
 *
 *   h^i(r) = exp(-i omega r*) sum_{n < K} c^i_n (r - 2M)^n,     r* = r + 2M log(r/2M - 1),
 *
 * for every component written to inhom_left/inhom_right (the coupled fields, then the gauge
 * fields, in the same order). The coupled fields follow from the series of the inner boundary
 * conditions of each basis solution (boundary_conditions.c) combined with the amplitudes C_in;
 * the gauge fields h2, h4 (even) and h8 (odd) from the Lorenz gauge conditions as in
 * gauge_fields.c, evaluated on the series. The static monopole and the static odd modes are
 * expanded from their closed forms. The even static modes (l >= 2, m = 0) are not computed here
 * (they are replaced by a separate calculation); horizon_series returns 0 for them.
 */
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include "GSF_ecc.h"

double R8_integrate_static_component(struct orbital_params *orbit, struct n_mode_data *n_mode, int piece);
double calc_alk(int l, int k);

/* --- truncated power series in x = r - 2M (M = 1), K terms ------------------------------ */
static void ser_zero(double complex *a, int K) { int n; for (n = 0; n < K; n++) a[n] = 0; }
static void ser_mul(const double complex *a, const double complex *b, double complex *c, int K)
{
	int n, k;
	double complex *t = malloc(K*sizeof(double complex));
	for (n = 0; n < K; n++) { t[n] = 0; for (k = 0; k <= n; k++) t[n] += a[k]*b[n-k]; }
	memcpy(c, t, K*sizeof(double complex));
	free(t);
}
/* 1/r = 1/(2 + x) */
static void ser_rinv(double complex *a, int K) { int n; for (n = 0; n < K; n++) a[n] = ((n % 2) ? -1.0 : 1.0)/pow(2.0, n + 1); }
/* f = x/(2 + x) */
static void ser_f(double complex *a, int K) { int n; a[0] = 0; for (n = 1; n < K; n++) a[n] = ((n % 2) ? 1.0 : -1.0)/pow(2.0, n); }
/* d/dr (the last coefficient is lost) */
static void ser_deriv(const double complex *a, double complex *d, int K) { int n; for (n = 0; n < K - 1; n++) d[n] = (n + 1)*a[n+1]; d[K-1] = 0; }
static void ser_poly(const double *p, int deg, double complex *a, int K) { int n; ser_zero(a, K); for (n = 0; n <= deg && n < K; n++) a[n] = p[n]; }

/* gauge field from the gauge condition, phase factored out:
 * R = exp(-i w r*) P and R' = exp(-i w r*) (P' - i w P/f), so
 *   (i f/w) (R_a' + g/r) -> exp(-i w r*) [(i f/w) (P_a' + g/r) + P_a] */
static void gauge_field(const double complex *Pa, const double complex *g, double omega, double complex *out, int K)
{
	double complex *f = malloc(K*sizeof(double complex)), *ri = malloc(K*sizeof(double complex)), *t = malloc(K*sizeof(double complex));
	int n;
	ser_f(f, K); ser_rinv(ri, K);
	ser_deriv(Pa, t, K);
	ser_mul(g, ri, out, K);
	for (n = 0; n < K; n++) t[n] += out[n];
	ser_mul(f, t, t, K);
	for (n = 0; n < K; n++) out[n] = I/omega*t[n] + Pa[n];
	free(f); free(ri); free(t);
}

int horizon_series(struct coupled_set *cset, struct n_mode_data *n_mode, struct orbital_params *orbit, int Kout, double complex *out)
{
	const int l = (int) lround(n_mode->l), m = (int) lround(n_mode->m);
	const double omega = n_mode->omega;
	const int K = Kout + 1;                       /* one extra term for the derivatives */
	const int nc = cset->num_coupled_fields, ng = cset->num_gauge_fields, nf = nc + ng;
	double complex *P = calloc((size_t) nf*K, sizeof(double complex));
	int i, d, n, k;

	if (l == 0 && m == 0) {
		/* static monopole (static_monopole.c): h1, h3, h6 inside the orbit, h2 = 0 */
		const double r0 = orbit->p, L = log((r0 - 2.0)/r0), s = sqrt((r0 - 3.0)*r0);
		double complex *ri = malloc(K*sizeof(double complex)), *ri2 = malloc(K*sizeof(double complex)),
		               *ri3 = malloc(K*sizeof(double complex)), *t = malloc(K*sizeof(double complex));
		double p[5];
		ser_rinv(ri, K); ser_mul(ri, ri, ri2, K); ser_mul(ri2, ri, ri3, K);
		/* h1 = -4 sqrt(pi) (r-2)^3 (r+2) (-4 + (r0-3) L) / (3 r^3 s),  (r-2)^3 (r+2) = x^3 (x + 4) */
		p[0] = 0; p[1] = 0; p[2] = 0; p[3] = 4; p[4] = 1;
		ser_poly(p, 4, t, K); ser_mul(t, ri3, t, K);
		for (n = 0; n < K; n++) P[0*K + n] = -4.0*sqrt(M_PI)*(-4.0 + (r0 - 3.0)*L)/(3.0*s)*t[n];
		/* h3 = 8 sqrt(pi) (-32 + r^3 - (r^3 - 8)(r0-3) L) / (3 r^2 s); r^3 = 8 + 12x + 6x^2 + x^3 */
		p[0] = -32 + 8 - 0*(r0 - 3.0)*L; p[1] = 12 - 12*(r0 - 3.0)*L; p[2] = 6 - 6*(r0 - 3.0)*L; p[3] = 1 - (r0 - 3.0)*L;
		p[0] = -24.0;
		ser_poly(p, 3, t, K); ser_mul(t, ri2, t, K);
		for (n = 0; n < K; n++) P[1*K + n] = 8.0*sqrt(M_PI)/(3.0*s)*t[n];
		/* h6 = 4 sqrt(pi) (-64 + 2 r^3 + (16 + r^3)(r0-3) L) / (3 r^2 s) */
		p[0] = -64 + 16 + 24*(r0 - 3.0)*L; p[1] = 24 + 12*(r0 - 3.0)*L; p[2] = 12 + 6*(r0 - 3.0)*L; p[3] = 2 + (r0 - 3.0)*L;
		ser_poly(p, 3, t, K); ser_mul(t, ri2, t, K);
		for (n = 0; n < K; n++) P[2*K + n] = 4.0*sqrt(M_PI)/(3.0*s)*t[n];
		free(ri); free(ri2); free(ri3); free(t);
	} else if (m == 0 && (l % 2) == 0) {
		free(P);
		return 0;                                 /* even static: not computed by h1Lorenz */
	} else if (m == 0) {
		/* static odd: h9 = h10 = 0, h8 = c_in phi_in(r) (gauge_fields.c, libSch.c);
		 * the odd dipole has the pure-gauge h9 = A/r^2 added (GSF_ecc.c) */
		const double c_in = -4.0*16.0*M_PI*orbit->L*normalized_associated_legendre_P_theta_deriv(l, 0, M_PI/2.0)
		                    *R8_integrate_static_component(orbit, n_mode, 0 /* IN_RE */);
		double complex *phi = malloc(K*sizeof(double complex)), *t = malloc(K*sizeof(double complex));
		if (l == 1) {
			double p[3] = {4.0, 4.0, 1.0};      /* r^2 */
			ser_poly(p, 2, phi, K);
		} else {
			/* (x/(1+x)) sum_k a_lk x^k with x = r/2 - 1 = (r - 2)/2, x/(1+x) = (r - 2)/r */
			ser_zero(t, K);
			for (k = 0; k <= l + 1 && k < K; k++) t[k] = calc_alk(l, k)/pow(2.0, k);
			double complex *ri = malloc(K*sizeof(double complex)), *xs = calloc(K, sizeof(double complex));
			ser_rinv(ri, K); xs[1] = 1.0;
			ser_mul(xs, ri, xs, K);
			ser_mul(t, xs, phi, K);
			free(ri); free(xs);
		}
		for (n = 0; n < K; n++) P[(nf - 1)*K + n] = c_in*phi[n];
		if (l == 1) {
			const double A = -256.0/3.0*sqrt(3.0*M_PI)*orbit->E*orbit->Omega_phi/(orbit->p - 2.0);
			double complex *ri = malloc(K*sizeof(double complex));
			ser_rinv(ri, K); ser_mul(ri, ri, t, K);
			for (n = 0; n < K; n++) P[0*K + n] += A*t[n];
			free(ri);
		}
		free(phi); free(t);
	} else {
		/* non-static: coupled fields from the inner boundary-condition series of each basis solution */
		double complex *S = malloc((size_t) nc*K*sizeof(double complex));
		double *a0s = malloc(nc*sizeof(double));
		struct boundary_conditions bcs;
		bcs.R = malloc(nc*sizeof(double complex));
		bcs.R_rs_deriv = malloc(nc*sizeof(double complex));
		bcs.l = n_mode->l;
		bcs.omega = omega;
		for (d = 0; d < nc; d++) {
			for (k = 0; k < nc; k++) a0s[k] = (k == d) ? 1.0 : 0.0;
			for (n = 0; n < nc*K; n++) S[n] = 0;
			bc_record_series(S, K);
			cset->in_bc_func(a0s, &bcs);
			bc_record_series(NULL, 0);
			for (i = 0; i < nc; i++)
				for (n = 0; n < K; n++) P[i*K + n] += n_mode->C_in[d]*S[i*K + n];
		}
		free(S); free(a0s); free(bcs.R); free(bcs.R_rs_deriv);

		double complex *g = malloc(K*sizeof(double complex)), *f = malloc(K*sizeof(double complex)), *t = malloc(K*sizeof(double complex));
		ser_f(f, K);
		if ((l + m) % 2 == 0) {
			/* even: coupled fields (1, 3, 5, 6[, 7]) at 0..nc-1, gauge fields 2, 4 at nc, nc+1 */
			const double complex *P1 = P, *P3 = P + K, *P5 = P + 2*K, *P6 = P + 3*K;
			const double complex *P7 = (nc > 4) ? P + 4*K : NULL;
			/* h2: (i f/w) (R1' - f R3' + (R1 - R5 - f R3 - 2 f R6)/r); with R1 - f R3 as the "a" field:
			 * R1' - f R3' = (R1 - f R3)' + f' R3, f' = (1 - f)/r */
			double complex *Pa = malloc(K*sizeof(double complex)), *ri = malloc(K*sizeof(double complex));
			ser_rinv(ri, K);
			ser_mul(f, P3, t, K);
			for (n = 0; n < K; n++) Pa[n] = P1[n] - t[n];
			/* g = r f' R3 + R1 - R5 - f R3 - 2 f R6 = (1 - f) R3 + R1 - R5 - f R3 - 2 f R6 */
			ser_mul(f, P6, g, K);
			for (n = 0; n < K; n++) g[n] = P3[n] - 2.0*t[n] + P1[n] - P5[n] - 2.0*g[n];
			gauge_field(Pa, g, omega, P + nc*K, K);
			/* h4: (i f/w) (R5' + (2 R5 + l(l+1) R6 - R7)/r) */
			if (ng > 1) {
				for (n = 0; n < K; n++) g[n] = 2.0*P5[n] + l*(l + 1.0)*P6[n] - (P7 ? P7[n] : 0.0);
				gauge_field(P5, g, omega, P + (nc + 1)*K, K);
			}
			free(Pa); free(ri);
		} else {
			/* odd: coupled fields (9, 10), gauge field 8: (i f/w) (R9' + (2 R9 - R10)/r) */
			const double complex *P9 = P, *P10 = (nc > 1) ? P + K : NULL;
			for (n = 0; n < K; n++) g[n] = 2.0*P9[n] - (P10 ? P10[n] : 0.0);
			gauge_field(P9, g, omega, P + nc*K, K);
		}
		free(g); free(f); free(t);
	}

	for (i = 0; i < nf; i++)
		for (n = 0; n < Kout; n++) out[i*Kout + n] = P[i*K + n];
	free(P);
	return nf;
}
