/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <assert.h>
#include <gsl/gsl_math.h>
#include "vector_stats.h"
#include "mh.h"

mcmclib_mh* mcmclib_mh_alloc(gsl_rng* r,
			     distrfun_p logdistr, void* logdistr_data,
			     mcmclib_mh_q* q, gsl_vector* x) {
  mcmclib_mh* p = (mcmclib_mh*) malloc(sizeof(mcmclib_mh));
  p->r = r;
  p->logdistr = logdistr;
  p->logdistr_data = logdistr_data;
  p->q = q;
  p->x = x;
  p->x_old = gsl_vector_alloc(x->size);
  p->last_accepted = 0;
  mcmclib_mh_reset(p);
  return p;
}

void mcmclib_mh_reset(mcmclib_mh* p) {
  p->logdistr_old = p->logdistr(p->logdistr_data, p->x);
}

void mcmclib_mh_free(mcmclib_mh* p) {
  gsl_vector_free(p->x_old);
  mcmclib_mh_q_free(p->q);
  free(p);
}

int mcmclib_mh_update(mcmclib_mh* p) {
  gsl_vector_memcpy(p->x_old, p->x);
  mcmclib_mh_q_sample(p->q, p->x);
  if(!mcmclib_vector_finite(p->x))
    GSL_ERROR("sampled a non-finite vector value", GSL_EDOM);
  p->last_accepted = mcmclib_mh_generic_step(p->r, p->x_old, p->x,
					     p->logdistr, p->logdistr_data,
					     &p->logdistr_old,
					     p->q);
  return(p->last_accepted);
}

void mcmclib_mh_update_N(mcmclib_mh* p, size_t N) {
  for(size_t i = 0; i < N; i++) {
    mcmclib_mh_update(p);
  }
}

int mcmclib_mh_generic_step(const gsl_rng* r, gsl_vector* old, gsl_vector* x,
			    distrfun_p logdistr, void* data,
			    double* plogdistr_old, mcmclib_mh_q* q) {
  double logdistr_new;
  double logdistr_old = plogdistr_old[0];
  double mh_offset, mh_ratio;

  if(!gsl_finite(logdistr_old)) {
    GSL_ERROR("non-finite probability density in the starting point", GSL_EDOM);
  }
  mh_offset = mcmclib_mh_q_ratio_offset(q, old, x);
  if(!gsl_finite(mh_offset)) {
    if(mh_offset < 0) {
      gsl_vector_memcpy(x, old);
      return 0;
    } else {
      plogdistr_old[0] = logdistr(data, x);
      return 1;
    }
  }

  logdistr_new = logdistr(data, x);
  if(!gsl_finite(logdistr_new)) {
    gsl_vector_memcpy(x, old);
    return 0;
  }

  mh_ratio = logdistr_new - logdistr_old + mh_offset;
  if((mh_ratio > 0) || (gsl_rng_uniform(r) <= exp(mh_ratio))) {
    plogdistr_old[0] = logdistr_new;
    return 1;
  }

  gsl_vector_memcpy(x, old);
  return 0;
}
