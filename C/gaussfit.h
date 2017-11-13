#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

struct data {
    size_t n;
    double * x;
    double * y;
    double * sigma;
};

int gauss_f (const gsl_vector * x, void *data, gsl_vector * f);
int gauss_df (const gsl_vector * x, void *data, gsl_matrix * J);
int gauss_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);



