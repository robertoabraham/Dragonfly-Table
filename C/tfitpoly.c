#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <gsl/gsl_multifit.h>
#include "table.h"
  
#define MAXROW 100000

char   *help[] = {
"",
"NAME",
"    tfit - fit a polynomial to table columns",
"",
"SYNOPSIS",
"    % tfit [OPTIONS] xcol ycol [scol] < table.txt ",
"",
"OPTIONS",
"    -v       Verbose mode", 
"    -V       Extra verbose mode (prints input data)", 
"    -n       Order of the polynomial (0=constant, 1=line, 2=parabola)", 
"",
"DESCRIPTION",
"",
"    This program fits a polynomial to named columns in a SExtractor",
"    ASCII table provided via standard input. If two columns are named",
"    the first two columns correspond to X and Y axes, respectively. If",
"    three columns are named the last column is the uncertainty on the Y",
"    axis measurements.",
"",
"AUTHOR",
"    Roberto Abraham (abraham@astro.utoronto.ca)",
"",
"LAST UPDATE",
"    July 31, 2012",
0};


void print_help()
{
    for (int i = 0; help[i] != 0; i++)
        fprintf(stdout,"%s\n",help[i]);
}
  

int main (int argc, char **argv)
{
    double x[MAXROW];
    double y[MAXROW];
    double sigma[MAXROW];
    int nrow = 0;
    int ncol;
    int count = 0;
    int status = 0;
    int i, n;
    double xi, yi, ei, chisq;
    gsl_matrix *X, *cov;
    gsl_vector *yvec, *wvec, *cvec;
    int verbose = 0;
    int extra_verbose = 0;
    char xcolname[64];
    char ycolname[64];
    char scolname[64];
    int has_uncertainties;
    int order = 2;
    int narg,c;

    while ((c = getopt (argc, argv, "vVhn:")) != -1)
        switch (c)
        {
            case 'v':
                verbose = 1;
                break;
            case 'V':
                verbose = 1;
                extra_verbose = 1;
                break;
            case 'n':
                order = atoi(optarg);
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
        has_uncertainties = 0;
    }
    else if (narg == 3) 
    { 
        sscanf(argv[optind++],"%s",xcolname);
        sscanf(argv[optind++],"%s",ycolname);
        sscanf(argv[optind++],"%s",scolname);
        has_uncertainties = 1;
    }
    else
    {
        print_help();
        return(1);
    }


    /* LOAD DATA COLUMNS */
    if (has_uncertainties)
        status = read_xyz(xcolname,ycolname,scolname,x,y,sigma,&nrow);
    else
        status = read_xy(xcolname,ycolname,x,y,&nrow);

    if (status)
    {
        fprintf(stderr,"Error reading data table.\n");
        exit(1);
    }

    /* Define sigma as unity for now */
    if (!has_uncertainties)
        for (int i=0; i<nrow; i++){
            sigma[i] = 1.0;
        }

    if (extra_verbose) {
        printf("# data:\n");
        for(int i=0;i<nrow;i++) printf("%20g %20g %20g\n",x[i],y[i],sigma[i]); 
    }

    /* EXECUTE THE FIT */
    ncol = order + 1;
    X = gsl_matrix_alloc(nrow, ncol);
    yvec = gsl_vector_alloc(nrow);
    wvec = gsl_vector_alloc(nrow);
    cvec = gsl_vector_alloc(ncol);
    cov = gsl_matrix_alloc(ncol, ncol);

    /* Example: when fitting a parabola we want this structure:      */
    /*                                                               */ 
    /* for (i = 0; i < nrow; i++) {                                  */
    /*    gsl_vector_set(yvec, i, y[i]);                             */
    /*    gsl_vector_set(wvec, i, 1.0/(sigma[i]*sigma[i]));          */
    /*    gsl_matrix_set(X, i, 0, 1.0);                              */
    /*    gsl_matrix_set(X, i, 1, x[i]);                             */
    /*    gsl_matrix_set(X, i, 2, x[i]*x[i]);                        */
    /* }                                                             */ 
    /*                                                               */
    /* This is generalized to n'th order below                       */

    for (i = 0; i < nrow; i++) {
        gsl_vector_set(yvec, i, y[i]);
        gsl_vector_set(wvec, i, 1.0/(sigma[i]*sigma[i]));
        for (int j = 0; j <= order; j++) 
            gsl_matrix_set(X, i, j, pow(x[i],j));
    }



    {
        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nrow, ncol);
        gsl_multifit_wlinear(X, wvec, yvec, cvec, cov, &chisq, work);
        gsl_multifit_linear_free(work);
    }

    #define C(i) (gsl_vector_get(cvec,(i)))
    #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

    if (verbose) {
        printf("# best fit: Y = %g", C(0));
        for (i=1;i<=order;i++)
            printf(" + %g X^%d",C(i),i);
        printf("\n");

        printf("# covariance matrix:\n");
        for (int i=0;i<=order;i++){
            for(int j=0; j<=order;j++){
                printf("%+.5e ",COV(i,j));
            }
            printf("\n");
        }

        printf("# chisq = %g\n", chisq);
        printf("# chisq_nu = %g\n", chisq/(nrow - order -1));
    }
    else {
        for (i=0;i<=order;i++)
            printf("%.10g ",C(i));
        printf("\n");

    }

    gsl_matrix_free (X);
    gsl_vector_free (yvec);
    gsl_vector_free (wvec);
    gsl_vector_free (cvec);
    gsl_matrix_free (cov);

    return 0;
}
