#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

/* expfit.c -- model functions for exponential + background */

struct data {
    size_t n;
    double * y;
    double * sigma;
};

int expb_f (const gsl_vector * x, void *data, gsl_vector * f);
int expb_df (const gsl_vector * x, void *data, gsl_matrix * J);
int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);



