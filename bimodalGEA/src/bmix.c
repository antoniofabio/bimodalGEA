#include <math.h>
#include <gsl/gsl_vector.h>
#include <R.h>
#include "gauss_am.h"

#define THETA_SIZE 4

static gsl_rng* rng;
static gsl_vector* y;
static mcmclib_amh* sampler;
static gsl_vector* x;

#define THETA_UNPACK				\
  double beta = gsl_vector_get(theta, 0);	\
  double mL = gsl_vector_get(theta, 1);		\
  double mH = gsl_vector_get(theta, 2);		\
  double s = gsl_vector_get(theta, 3);

static void theta_print(const gsl_vector* theta) {
  THETA_UNPACK;
  Rprintf("[%.2f; %.2f; %.2f; %.2f]", beta, mL, mH, s);
}

static double Rsquare(const gsl_vector* theta) {
  THETA_UNPACK;
  double d = mH - mL;
  double varBetween = beta * (1.0 - beta) * d * d;
  return varBetween / (varBetween + (s * s));
}

static double phi(double yi, double m, double s2) {
  double d = yi - m;
  double e2 = d * d / s2;
  return exp(- e2) / sqrt(2.0 * M_PI * s2);
}

static double lf(void* ignore, const gsl_vector* theta) {
  ignore = NULL; /*keep the compiler quiet*/
  /* Rprintf("..in function 'lf'..\n"); */
  double out = 0.0;
  THETA_UNPACK;
  double s2 = s * s;
  /*PRIOR*/
  if((beta < 0.01) | (beta > 0.99) |
     (fabs(mL) > 10.0) | (fabs(mH) > 10.0) | (mL > mH) |
     (s < 0.01) | (s > 10.0)) {
    /* Rprintf("we're out of the distribution domain\n"); */
    return(log(0.0));
  }
  out += log(phi(mL, 0, 100.0));
  out += log(phi(mH, 0, 100.0));
  /* out += -11.0 * log(s2); */
  /*LIKELIHOOD*/
  /* Rprintf("computing likelihood\n"); */
  /* theta_print(theta); */
  /* Rprintf("\n"); */
  for(size_t i = 0; i < y->size; i++) {
    double yi = gsl_vector_get(y, i);
    /* Rprintf("y[%zd] = %.3f -> ", i, yi); */
    double li = log(beta * phi(yi, mL, s2) +
		    (1.0 - beta) * phi(yi, mH, s2));
    /* Rprintf("%.2f\n", li); */
    out += li;
  }
  return out;
}

static void Rgsl_handler(const char * reason,
			 const char * file,
			 int line,
			 int gsl_errno) {
  error("gsl error [%d]: %s:%d. %s", gsl_errno, file, line, reason);
}

void bmix_init(double* in_y, int* in_N,
		 double* x0,
		 double* S0,
		 int* t0) {
  gsl_set_error_handler(Rgsl_handler);
  x = gsl_vector_alloc(THETA_SIZE);
  gsl_vector_view x0v = gsl_vector_view_array(x0, THETA_SIZE);
  gsl_vector_memcpy(x, &(x0v.vector));

  gsl_vector_view yv = gsl_vector_view_array(in_y, (size_t) in_N[0]);
  y = gsl_vector_alloc((size_t) in_N[0]);
  gsl_vector_memcpy(y, &(yv.vector));

  rng = gsl_rng_alloc(gsl_rng_default);

  gsl_matrix_view sigma_zero = gsl_matrix_view_array(S0, THETA_SIZE, THETA_SIZE);
  sampler = mcmclib_gauss_am_alloc(rng,
				   lf, NULL,
				   x,
				   &sigma_zero.matrix,
				   (size_t) t0[0]);

  mcmclib_gauss_am_set_sf(sampler, 2.38 * 2.38 / 4.0);
}

void bmix_update(int* N, int* thin, double* out) {
  for(int i = 0; i < N[0]; i++) {
    for(int j = 0; j < thin[0]; j++) {
      mcmclib_amh_update(sampler);
    }
    /* theta_print(x); */
    /* Rprintf("\n"); */
    out[i] = Rsquare(x);
  }
}

void bmix_free() {
  mcmclib_amh_free(sampler);
  gsl_rng_free(rng);
  gsl_vector_free(y);
  gsl_vector_free(x);
}
