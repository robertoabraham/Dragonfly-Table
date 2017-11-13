#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

//#include "expfit.h"
#include "gaussfit.h"
#include "table.h"

#define MAXROW 100000
void print_state (size_t iter, gsl_multifit_fdfsolver * s);


char   *help[] = {
"",
"NAME",
"    fithist - fit a gaussian to a histogram",
"",
"SYNOPSIS",
"    % fithist [OPTIONS] xcol ycol < histogram.txt ",
"",
"OPTIONS",
"    -v       Verbose mode", 
"",
"EXAMPLE",
"    Generate an image with gaussian noise (mean=100, stddev=50) and then fit",
"    a normal distribution to its histogram:",
"",
"    % imcalc -c 256 256 'grand 50 * 100 +' rand.fits",
"    % imhist < rand.fits > hist.txt",
"    % fithist PIXVAL COUNT < hist.txt",
"",
"DESCRIPTION",
"",
"    This program fits a gaussian to a histogram provided to the program",
"    via standard input. The histogram must be a SExtractor-format table.",
"    If the histogram is produced by IMHIST then the xcol and ycol strings",
"    are PIXVAL and COUNT.",
"",
"AUTHOR",
"    Roberto Abraham (abraham@astro.utoronto.ca)",
"",
"LAST UPDATE",
"    July 2012",
0};


void print_help()
{
    for (int i = 0; help[i] != 0; i++)
        fprintf(stdout,"%s\n",help[i]);
}


void print_state (size_t iter, gsl_multifit_fdfsolver *s)
{
    printf ("iter: %3u      parameters = % 15.8f % 15.8f % 15.8f "
            "|f(p)| = %g\n",
            (unsigned int)iter,
            gsl_vector_get (s->x, 0), 
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2), 
            gsl_blas_dnrm2 (s->f));
}


int main (int argc, char **argv)
{

    int verbose = 0;
    int quiet = 0;
    char xcolname[64];
    char ycolname[64];
    int narg,c;

    while ((c = getopt (argc, argv, "vqh")) != -1)
        switch (c)
        {
            case 'v':
                verbose = 1;
                break;
            case 'q':
                quiet = 1;
                break;
            case 'h':
                print_help();
                return(0);
                break;
            case '?':
                if (optopt == 'c')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                print_help();
                return 1;
            default:
                abort();
        }
    narg = argc - optind;
    if (narg == 2) 
    { 
        sscanf(argv[optind++],"%s",xcolname);
        sscanf(argv[optind++],"%s",ycolname);
    }
    else
    {
        print_help();
        return(1);
    }

    /* LOAD HISTOGRAM */

    double x[MAXROW];
    double y[MAXROW];
    double sigma[MAXROW];
    int nrow = 0;
    int count = 0;
    int status = 0;

    status = read_xy(xcolname,ycolname,x,y,&nrow);
    if (status)
    {
        fprintf(stderr,"Error reading data table.\n");
        exit(1);
    }
    // Weed out zero entries for now - assuming poisson statistics 
    // they contribute infinite variance
    for (int i=0; i<nrow; i++){
        if (y[i] > 0){
            y[count] = y[i];
            x[count] = x[i];
            sigma[count] = sqrt(y[i]);   /* assume poisson errors */
            count ++;
        }
    }
    nrow = count;
    if (verbose) for(int i=0;i<nrow;i++) printf("%20g %20g\n",x[i],y[i]);


    /* DETERMINE INITIAL GUESSES */

    double norm = 0;
    double mean = 0;
    double rms = 0;
    double npix = 0;
    double mode = 0;
    double maxcount = 0;
    for (int i=0; i<nrow; i++) {
        if (y[i]>maxcount){
            maxcount = y[i];
            mode = x[i];
        }
    }
    norm = maxcount;
    mean = mode;
    rms = sqrt(mode);
    if (!quiet)
        printf("Initial guess: A=%g mu=%g sig=%g\n",norm,mean,rms);

    /* FIT THE DATA */

    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    unsigned int i, iter = 0;
    const size_t p = 3;          /* number of free parameters */
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    struct data d = {nrow, x, y, sigma};
    gsl_multifit_function_fdf f;
    double x_init[3] = {norm, mean, rms};
    gsl_vector_view xx = gsl_vector_view_array (x_init, p);
    const gsl_rng_type * type;
    gsl_rng * r;

    /* Set up GSL */
    gsl_rng_env_setup();
    type = gsl_rng_default;
    r = gsl_rng_alloc (type);

    /* Define the model to fit */
    f.f = &gauss_f;
    f.df = &gauss_df;
    f.fdf = &gauss_fdf;
    f.n = nrow;
    f.p = p;
    f.params = &d;

    /* Define the solver we wish to use */
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(T, nrow, p);
    gsl_multifit_fdfsolver_set(s, &f, &xx.vector);

    if (verbose) print_state (iter, s);

    /* Fit the model to the data */
    do
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate(s);

        if (verbose) printf ("status = %s\n", gsl_strerror (status));
        if (verbose) print_state (iter, s);

        if (status)
            break;

        status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
    while (status == GSL_CONTINUE && iter < 500);

    /* Compute covariance */
    gsl_multifit_covar (s->J, 0.0, covar);

    { 
        double chi = gsl_blas_dnrm2(s->f);
        double dof = nrow - p;
        double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

        if (!quiet) {
            printf ("A    = %.5f +/- %.5f\n", gsl_vector_get(s->x, 0), 
                    c*sqrt(gsl_matrix_get(covar,0,0)));
            printf ("mu   = %.5f +/- %.5f\n", gsl_vector_get(s->x, 1), 
                    c*sqrt(gsl_matrix_get(covar,1,1)));
            printf ("sig  = %.5f +/- %.5f\n", gsl_vector_get(s->x, 2), 
                    c*sqrt(gsl_matrix_get(covar,2,2)));
            printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
        }
        else 
            printf("%.2f %.2f %.2f\n",
                    gsl_vector_get(s->x, 0), 
                    gsl_vector_get(s->x, 1),
                    gsl_vector_get(s->x, 2));

    }

    if (!quiet)
        printf ("status = %s\n", gsl_strerror (status));

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_rng_free (r);
    return 0;
}


