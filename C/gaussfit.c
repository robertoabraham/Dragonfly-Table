#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "gaussfit.h"

/* gaussfit.c -- functions for gaussian fits */


/* Residuals from the weighted gaussian model */
int gauss_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    size_t n = ((struct data *)data)->n;
    double *x = ((struct data *)data)->x;
    double *y = ((struct data *)data)->y;
    double *sigma = ((struct data *) data)->sigma;

    double A = gsl_vector_get(p, 0);
    double mu = gsl_vector_get(p, 1);
    double sig = gsl_vector_get(p, 2);

    size_t i;

    for (i = 0; i < n; i++)
    {
        /* Model Yi = A * exp(-lambda * i) + b */
        double Yi = A * exp (-0.5*gsl_pow_2((x[i]-mu)/sig));
        gsl_vector_set(f, i, (Yi - y[i])/sigma[i]);
    }

    return GSL_SUCCESS;
}

/* A function that sets the Jacobian matrix of the system              */
/*                                                                     */
/* The Jacobian matrix is J(i,j) = dfi / dpj, where                    */
/*                                                                     */
/*    fi = (Yi - yi)/sigma[i],                                         */
/*    Yi = A * exp (-0.5*gsl_pow_2((t-mu)/sig));   [this is the model] */
/*    yi are the data                                                  */ 
/*    pj are the parameters (A,mu,sig)                                 */
/*                                                                     */
/* Note that the terms in the Jacobian are easily computed             */
/* using Mathematica. Setting t as the dependent variable,             */
/* the terms are computed as follows:                                  */
/*                                                                     */
/* fi = (A E^(-(t-mu)^2/(2*sig^2))-yi)/sigmai                          */
/* D[fi,A]                                                             */
/* D[fi,mu]                                                            */
/* D[fi,sig]                                                           */
int gauss_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    size_t n = ((struct data *)data)->n;
    double *x = ((struct data *)data)->x;
    double *sigma = ((struct data *) data)->sigma;

    double A = gsl_vector_get(p, 0);
    double mu = gsl_vector_get(p, 1);
    double sig = gsl_vector_get(p, 2);

    size_t i;

    for (i = 0; i < n; i++)
    {
        double t = x[i];
        double s = sigma[i];
        double e = exp(-0.5*gsl_pow_2((t-mu)/sig));
        gsl_matrix_set(J, i, 0, e/s);                                  /* dfi/dA   */
        gsl_matrix_set(J, i, 1, A*e*(t-mu)/gsl_pow_2(sig)/s);          /* dfi/dmu  */
        gsl_matrix_set(J, i, 2, A*e*gsl_pow_2(t-mu)/gsl_pow_3(sig)/s); /* dfi/dsig */
    }
    return GSL_SUCCESS;
}

/* Convenience function to evaluate both the model and its Jacobian */
int gauss_fdf(const gsl_vector *p, void *data, gsl_vector *f, gsl_matrix *J)
{
    gauss_f(p, data, f);
    gauss_df(p, data, J);

    return GSL_SUCCESS;
}



