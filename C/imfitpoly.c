#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <fitsio.h>
#include "mfits.h"

char   *help[] = {
"",
"NAME",
"    imfitpoly - fit a polynomial to a FITS image",
"",
"SYNOPSIS",
"    imfitpoly [OPTIONS] input.fits",
"",
"OPTIONS",
"    -n            Order of the polynomial (0=constant, 1=linear, 2=quadratic, 3=cubic) [default 1]", 
"    -o file.fits  Output filename [default a.fits]",
"    -s sigma.fits Input error map. This is the sigma_i in $\\Sum(((y-y_i)/sigma_i)^2)$",
"    -h            Print help",
"    -v            Verbose mode", 
"",
"DESCRIPTION",
"    This program fits a polynomial to a FITS image, producing an output FITS image corresponding",
"    to the best-fit model. In verbose mode a listing of the best-fit parameters and also the",
"    goodness-of-fit information (chi squared and reduced chi squared) is output. The default",
"    is to save the output image to a file named 'a.fits', though the output filename can be set",
"    explicitly using the -o option. Image pixels can be weighted by supplying an error map. This",
"    map corresponds to the 1-sigma error on the pixel. This is NOT an inverse variance (weight) map.",
"    In the case where the image has counts in ADU and no systematic sources of error (bad pixels etc)",
"    then the error map should just be the square root of the original image. If you want to",
"    totally de-emphasize a particular pixel give it a huge value on the error map.",
"",
"AUTHOR",
"    Roberto Abraham (abraham@astro.utoronto.ca)",
"",
"LAST UPDATE",
"    Feb 2013",
0};


void print_help()
{
    for (int i = 0; help[i] != 0; i++)
        fprintf(stdout,"%s\n",help[i]);
}
  

int main (int argc, char **argv)
{
    int ndata = 0;
    int npar;
    int count = 0;
    int i, j, n;
    double xi, yi, ei, chisq;
    gsl_matrix *X, *cov;
    gsl_vector *yvec, *sigvec, *cvec;
    int verbose = 0;
    int has_uncertainties;
    int order = 1;
    double *pix,*sigpix;
    double rmode;
    int nx, ny, npts;
    int narg;
    int status = 0;
    int c;
    char imname[1024];
    FILE *outfile, *sigfile;
    char *signame;
    char *outname = "a.fits";
    int use_sigma_map = 0;

    while ((c = getopt (argc, argv, "vn:o:s:h")) != -1)
        switch (c)
        {
            case 'v':
                verbose = 1;
                break;
            case 'n':
                order = atoi(optarg);
                if (order >= 4) {
                    fprintf(stderr,"Order must be less than or equal to 3\n");
                    abort();
                }
                break;
            case 'o':
                outname = optarg;
                break;
           case 's':
                signame = optarg;
                use_sigma_map = 1;
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

    /* Handle non-option arguments */
    narg = argc - optind;
    if (narg != 1)
    {
        print_help();
        return(1);
    }

    /* If no explicit output name has been supplied with -o then use the
     * default one: a.fits (This in intended to be reminiscent of a compiler,
     * by the way, with a.fits playing the role of a.out). If a.fits
     * already exists we nuke it so we can replace it. */
    if ((outfile = fopen(outname, "r")))
    {
        fclose(outfile);
        unlink(outname);
    }

    /* Read in the image */
    sprintf(imname,"%s",argv[optind++]);
    pix = readimage(imname, &nx, &ny, &status);
    if (pix == NULL) 
    {
        printf("Memory allocation error\n");
        return(1);
    }
    ndata = nx*ny;

    /* Weighting is optional */
    if (use_sigma_map) {
        if ((sigfile = fopen(signame, "r"))){
            int wnx, wny;
            sigpix = readimage(signame, &wnx, &wny, &status);
            if (pix == NULL) 
            {
                printf("Memory allocation error (when reading sigma file)\n");
                return(1);
            }
        }
        else {
            fprintf(stderr,"Sigma file not found\n");
            abort();
        }
    }


    /* Now do the heavy lifting! */

    /* Define the number of parameters */
    switch(order)
    {
        case 0: npar = 1;
                break;
        case 1: npar = 4;
                break;
        case 2: npar = 9;
                break;
        case 3: npar = 16;
                break;
    }

    /* Allocate storage */
    X = gsl_matrix_alloc(ndata, npar);
    yvec = gsl_vector_alloc(ndata);
    sigvec = gsl_vector_alloc(ndata);
    cvec = gsl_vector_alloc(npar);
    cov = gsl_matrix_alloc(npar, npar);

    /* Define the model:
     *
     * Substitute so that the following are the dependent variables:
     *
     * Order = 0:
     * Order = 1: {1, x, y, x y}
     * Order = 2: {1, x, x^2, y, x y, x^2 y, y^2, x y^2, x^2 y^2}
     * Order = 3: {1, x, x^2, x^3, y, x y, x^2 y, x^3 y, y^2, x y^2, x^2 y^2, x^3 y^2, y^3, x y^3, x^2 y^3, x^3 y^3}
     *

     * Example: if we were fitting a 1D parabola we would want this structure:      
     *                                                                
     * for (i = 0; i < ndata; i++) {                                 
     *    gsl_vector_set(yvec, i, y[i]);                            
     *    gsl_vector_set(sigvec, i, 1.0/(sigma[i]*sigma[i]));        
     *    gsl_matrix_set(X, i, 0, 1.0);                           
     *    gsl_matrix_set(X, i, 1, x[i]);                         
     *    gsl_matrix_set(X, i, 2, x[i]*x[i]);                   
     * } 
     * */ 

    for (i = 0; i < ndata; i++) {
        double x, y;
        int colnum;

        x = (double) floor(i/ny);
        y = (double) i - x*nx;

        // Load data into 1D vectors
        gsl_vector_set(yvec, i, *(pix+i));              /* Image */
        if (use_sigma_map) {
            gsl_vector_set(sigvec, i, *(sigpix+i));     /* Error map */
        }
        else {                         
            /* Give everything unit weight */
            gsl_vector_set(sigvec, i, 1.0);
        }

        /* Load the model corresponding to this data point. The order of the terms below corresponds
         * to that given in my Mathematica notebook (PolynomialTerms.nb) to which the reader
         * is referred to for reference */

        switch(order) 
        {

            // Order = 1: constant
            // {1}
            case 0:  gsl_matrix_set(X, i, 0, 1.0); 
                     break;

            // Order = 1: plane
            // {1, x, y, x y}
            case 1:  gsl_matrix_set(X, i, 0, 1.0); 
                     gsl_matrix_set(X, i, 1, x);
                     gsl_matrix_set(X, i, 2, y);
                     gsl_matrix_set(X, i, 3, x*y);
                     break;

            // Order = 2: paraboloid 
            // {1, x, x^2, y, x y, x^2 y, y^2, x y^2, x^2 y^2}
            case 2:  gsl_matrix_set(X, i, 0, 1.0); 
                     gsl_matrix_set(X, i, 1, x);
                     gsl_matrix_set(X, i, 2, x*x);
                     gsl_matrix_set(X, i, 3, y);
                     gsl_matrix_set(X, i, 4, x*y);
                     gsl_matrix_set(X, i, 5, x*x*y);
                     gsl_matrix_set(X, i, 6, y*y);
                     gsl_matrix_set(X, i, 7, x*y*y);
                     gsl_matrix_set(X, i, 8, x*x*y*y);
                     break;

            // Order = 3: cubic
            // {1, x, x^2, x^3, y, x y, x^2 y, x^3 y, y^2, x y^2, x^2 y^2, x^3 y^2, y^3, x y^3, x^2 y^3, x^3 y^3}
            case 3:  gsl_matrix_set(X, i, 0, 1.0); 
                     gsl_matrix_set(X, i, 1, x);
                     gsl_matrix_set(X, i, 2, x*x);
                     gsl_matrix_set(X, i, 3, x*x*x);
                     gsl_matrix_set(X, i, 4, y);
                     gsl_matrix_set(X, i, 5, x*y);
                     gsl_matrix_set(X, i, 6, x*x*y);
                     gsl_matrix_set(X, i, 7, x*x*x*y);
                     gsl_matrix_set(X, i, 8, y*y);
                     gsl_matrix_set(X, i, 9, x*y*y);
                     gsl_matrix_set(X, i, 10, x*x*y*y);
                     gsl_matrix_set(X, i, 11, x*x*x*y*y);
                     gsl_matrix_set(X, i, 12, y*y*y);
                     gsl_matrix_set(X, i, 13, x*y*y*y);
                     gsl_matrix_set(X, i, 14, x*x*y*y*y);
                     gsl_matrix_set(X, i, 15, x*x*x*y*y*y);
        }

    }


    // Compute the answer
    {
        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(ndata, npar);
        gsl_multifit_wlinear(X, sigvec, yvec, cvec, cov, &chisq, work);
        gsl_multifit_linear_free(work);
    }


    #define C(i) (gsl_vector_get(cvec,(i)))
    #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

    if (verbose) {
        switch (order) {
            case 0: printf("# {1}\n");
                    break;
            case 1: printf("# {1, x, y, x y}\n");
                    break;
            case 2: printf("# {1, x, x^2, y, x y, x^2 y, y^2, x y^2, x^2 y^2}\n");
                    break;
            case 3: printf("# {1, x, x^2, x^3, y, x y, x^2 y, x^3 y, y^2, x y^2, x^2 y^2, x^3 y^2, y^3, x y^3, x^2 y^3, x^3 y^3}\n");
                    break;
        }

        printf("# best fit parameters:\n");
        for (i=0;i<npar;i++)
            printf("%.10g ",C(i));
        printf("\n");

        printf("# covariance matrix:\n");
        for (int i=0;i<npar;i++){
            for(int j=0; j<npar;j++){
                printf("%+.5e ",COV(i,j));
            }
            printf("\n");
        }

        printf("# chisq = %g\n", chisq);
        printf("# chisq_nu = %g\n", chisq/(ndata - npar -1));
    }
    else {

    }

    // Generate the output image
    {
        double *out;
        double x,y;
        out = (double *) malloc(nx*ny*sizeof(double));
        switch(order) {

            case 0: for (i=0;i<ndata;i++) {
                        x = (double) floor(i/ny);
                        y = (double) i - x*nx;
                        *(out + count) = C(0);
                        count++;
                    };
                    break;

            case 1: for (i=0;i<ndata;i++) {
                        x = (double) floor(i/ny);
                        y = (double) i - x*nx;
                        *(out + count) = C(0) + C(1)*x + C(2)*y + C(3)*x*y;
                        count++;
                    };
                    break;

            case 2: for (i=0;i<ndata;i++) {
                        x = (double) floor(i/ny);
                        y = (double) i - x*nx;
                        *(out + count) = C(0) + C(1)*x + C(2)*x*x +  C(3)*y +  C(4)*x*y +  
                                         C(5)*x*x*y +  C(6)*y*y +  C(7)*x*y*y +  C(8)*x*x*y*y;
                        count++;
                    };
                    break;

            case 3: for (i=0;i<ndata;i++) {
                        x = (double) floor(i/ny);
                        y = (double) i - x*nx;
                        *(out + count) = C(0) + C(1)*x + C(2)*x*x + C(3)*x*x*x + C(4)*y +  C(5)*x*y +  
                                         C(6)*x*x*y  +  C(7)*x*x*x*y +  C(8)*y*y +  C(9)*x*y*y + 
                                         C(10)*x*x*y*y + C(11)*x*x*x*y*y + C(12)*y*y*y + C(13)*x*y*y*y + 
                                         C(14)*x*x*y*y*y + C(15)*x*x*x*y*y*y;
                        count++;
                    };   
        }
        writeimage(outname, out, nx, ny, &status); 
    }

    gsl_matrix_free (X);
    gsl_vector_free (yvec);
    gsl_vector_free (sigvec);
    gsl_vector_free (cvec);
    gsl_matrix_free (cov);

    return 0;
}

