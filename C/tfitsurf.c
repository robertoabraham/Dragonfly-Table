#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include "table.h"

#define MAXROW 100000

char   *help[] = {
"",
"NAME",
"    tfitsurf - fit a polynomial to a set of X, Y, Z data points",
"",
"SYNOPSIS",
"    tfitsurf [OPTIONS] xcol ycol zcol [sigma_col] < table.txt",
"",
"OPTIONS",
"    -n            Order of the polynomial (0=constant, 1=ramp, 2=paraboloid, 3=bicubic) [default 1]", 
"    -h            Print help",
"    -v            Verbose mode", 
"",
"DESCRIPTION",
"    This program fits a polynomial to a set of X,Y,Z data points.",
"",
"AUTHOR",
"    Roberto Abraham (abraham@astro.utoronto.ca)",
"",
"LAST UPDATE",
"    July 2013",
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
    double z[MAXROW];
    double s[MAXROW];
    char xcolname[64];
    char ycolname[64];
    char zcolname[64];
    char scolname[64];

    int nrow = 0;
    int npar;
    int count = 0;
    int i, j, n;
    double xi, yi, ei, chisq;
    gsl_matrix *X, *cov;
    gsl_vector *zvec, *sigvec, *cvec;
    int verbose = 0;
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
    int has_uncertainties = 0;
    int extra_verbose = 0;

    while ((c = getopt (argc, argv, "vn:o:h")) != -1)
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
    if (narg == 3) 
    { 
        sscanf(argv[optind++],"%s",xcolname);
        sscanf(argv[optind++],"%s",ycolname);
        sscanf(argv[optind++],"%s",zcolname);
        has_uncertainties = 0;
    }
    else if (narg == 4) 
    { 
        sscanf(argv[optind++],"%s",xcolname);
        sscanf(argv[optind++],"%s",ycolname);
        sscanf(argv[optind++],"%s",zcolname);
        sscanf(argv[optind++],"%s",scolname);
        has_uncertainties = 1;
    }
    else
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

    /* Load data columns */
    if (has_uncertainties)
        status = read_xyzs(xcolname,ycolname,zcolname,scolname,x,y,z,s,&nrow);
    else
        status = read_xyz(xcolname,ycolname,zcolname,x,y,z,&nrow);

    if (status)
    {
        fprintf(stderr,"Error reading data table.\n");
        exit(1);
    }

    /* Define sigma as unity for now */
    if (!has_uncertainties)
        for (int i=0; i<nrow; i++){
            s[i] = 1.0;
        }

    if (extra_verbose) {
        printf("# data:\n");
        for(int i=0;i<nrow;i++) printf("%20g %20g %20g %20g\n",x[i],y[i],z[i],s[i]); 
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
    X = gsl_matrix_alloc(nrow, npar);
    zvec = gsl_vector_alloc(nrow);
    sigvec = gsl_vector_alloc(nrow);
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
     *    gsl_vector_set(zvec, i, y[i]);                            
     *    gsl_vector_set(sigvec, i, 1.0/(sigma[i]*sigma[i]));        
     *    gsl_matrix_set(X, i, 0, 1.0);                           
     *    gsl_matrix_set(X, i, 1, x[i]);                         
     *    gsl_matrix_set(X, i, 2, x[i]*x[i]);                   
     * } 
     * */ 

    for (i = 0; i < nrow; i++) {
        int colnum;

        // Load data into 1D vectors
        gsl_vector_set(zvec, i, z[i]);            /* Value to fit */
        if (use_sigma_map) {
            gsl_vector_set(sigvec, i, s[i]);      /* Error value */
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
                     gsl_matrix_set(X, i, 1, x[i]);
                     gsl_matrix_set(X, i, 2, y[i]);
                     gsl_matrix_set(X, i, 3, x[i]*y[i]);
                     break;

           // Order = 2: paraboloid 
           // {1, x, x^2, y, x y, x^2 y, y^2, x y^2, x^2 y^2}
            case 2:  gsl_matrix_set(X, i, 0, 1.0); 
                     gsl_matrix_set(X, i, 1, x[i]);
                     gsl_matrix_set(X, i, 2, x[i]*x[i]);
                     gsl_matrix_set(X, i, 3, y[i]);
                     gsl_matrix_set(X, i, 4, x[i]*y[i]);
                     gsl_matrix_set(X, i, 5, x[i]*x[i]*y[i]);
                     gsl_matrix_set(X, i, 6, y[i]*y[i]);
                     gsl_matrix_set(X, i, 7, x[i]*y[i]*y[i]);
                     gsl_matrix_set(X, i, 8, x[i]*x[i]*y[i]*y[i]);
                     break;

            // Order = 3: cubic
            // {1, x, x^2, x^3, y, x y, x^2 y, x^3 y, y^2, x y^2, x^2 y^2, x^3 y^2, y^3, x y^3, x^2 y^3, x^3 y^3}
            case 3:  gsl_matrix_set(X, i, 0, 1.0); 
                     gsl_matrix_set(X, i, 1, x[i]);
                     gsl_matrix_set(X, i, 2, x[i]*x[i]);
                     gsl_matrix_set(X, i, 3, x[i]*x[i]*x[i]);
                     gsl_matrix_set(X, i, 4, y[i]);
                     gsl_matrix_set(X, i, 5, x[i]*y[i]);
                     gsl_matrix_set(X, i, 6, x[i]*x[i]*y[i]);
                     gsl_matrix_set(X, i, 7, x[i]*x[i]*x[i]*y[i]);
                     gsl_matrix_set(X, i, 8, y[i]*y[i]);
                     gsl_matrix_set(X, i, 9, x[i]*y[i]*y[i]);
                     gsl_matrix_set(X, i, 10, x[i]*x[i]*y[i]*y[i]);
                     gsl_matrix_set(X, i, 11, x[i]*x[i]*x[i]*y[i]*y[i]);
                     gsl_matrix_set(X, i, 12, y[i]*y[i]*y[i]);
                     gsl_matrix_set(X, i, 13, x[i]*y[i]*y[i]*y[i]);
                     gsl_matrix_set(X, i, 14, x[i]*x[i]*y[i]*y[i]*y[i]);
                     gsl_matrix_set(X, i, 15, x[i]*x[i]*x[i]*y[i]*y[i]*y[i]);
        }

    }


    // Compute the answer
    {
        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nrow, npar);
        gsl_multifit_wlinear(X, sigvec, zvec, cvec, cov, &chisq, work);
        gsl_multifit_linear_free(work);
    }


    #define C(i) (gsl_vector_get(cvec,(i)))
    #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

    // Output the answer in JSON format. This must pass
    // validation at jsonlint.com

    printf("{\n");
    printf("  \"type\": \"polynomial_surface\",\n");
    printf("  \"order\": %d,\n",order);
    switch (order) {
        case 0: printf("  \"terms\": [\"1\"],\n");
                break;
        case 1: printf("  \"terms\": [\"1\", \"x\", \"y\", \"x*y\"],\n");
                break;
        case 2: printf("  \"terms\": [\"1\", \"x\", \"x**2\", \"y\", \"x*y\", \"(x**2)*y\", \"y**2\", \"x*(y**2)\", \"(x**2)*(y**2)\"],\n");
                break;
        case 3: printf("  \"terms\": [\"1\", \"x\", \"x**2\", \"x**3\", \"y\", \"x*y\", \"(x**2)*y\", \"(x**3)*y\", \"y**2\", \"x*(y**2)\", \"(x**2)*(y**2)\", \"(x**3)*(y**2)\", \"y**3\", \"x*(y**3)\", \"(x**2)*(y**3)\", \"(x**3)*(y**3)\"],\n");
                break;
    }

    printf("  \"coefficients\": [");
    for (i=0;i<npar;i++){
        printf("%.10g",C(i));
        if (i<(npar-1))
            printf(", ");
    }
    printf("],\n");

    printf("  \"equation\": ");
    switch (order) {
        case 0: printf(" \"%.5e\",\n",C(0));
                break;
        case 1: printf(" \"%.5e + %.5e*x + %.5e*y + %.5e*x*y\",\n",
                        C(0),C(1),C(2),C(3));
                break;
        case 2: printf(" \"%.5e + %.5e*x + %.5e*x**2 + %.5e*y + %.5e*x*y + %.5e*(x**2)*y + %.5e*y**2 + %.5e*x*(y**2) + %.5e*(x**2)*(y**2)\",\n",
                        C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8));
                break;
        case 3: printf(" \"%.5e + %.5e*x + %.5e*x**2 + %.5e*x**3 + %.5e*y + %.5e*x*y + %.5e*(x**2)*y + %.5e*(x**3)*y + %.5e*y**2 + %.5e*x*(y**2) + %.5e*(x**2)*(y**2) + %.5e*(x**3)*(y**2) + %.5e*y**3 + %.5e*x*(y**3) + %.5e*(x**2)*(y**3) + %.5e*(x**3)*(y**3)\",\n",
                        C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),C(12),C(13),C(14),C(15));
                break;
    }
 
    printf("  \"covariance matrix\":\n");
    printf("    [\n");
    for (int i=0;i<npar;i++){
        printf("      [");
        for(int j=0; j<npar;j++){
            printf("%.5e",COV(i,j));
            if (j<(npar-1))
                printf(", ");
        }
        printf("]");
        if (i<(npar-1))
            printf(", ");
        printf("\n");
    }
    printf("    ],\n");

    if (has_uncertainties)
        printf("  \"axes\": [\"%s\",\"%s\",\"%s\",\"%s\"],\n",xcolname,ycolname,zcolname,scolname);
    else
        printf("  \"axes\": [\"%s\",\"%s\",\"%s\"],\n",xcolname,ycolname,zcolname);

    printf("  \"chisq\": %g,\n", chisq);
    printf("  \"chisq_nu\": %g,\n", chisq/(nrow - npar -1));
    printf("  \"data_has_uncertainties\": %d,\n", has_uncertainties);
    printf("  \"ndata\": %d,\n",nrow);
    printf("  \"data\":\n");
    printf("    [\n");
    // Triplets
    for (int i=0;i<nrow;i++){
        if (has_uncertainties)
            printf("      [%.5f, %.5f, %.5f, %.5f]",x[i],y[i],z[i],s[i]);
        else 
            printf("      [%.5f, %.5f, %.5f]",x[i],y[i],z[i]);
        if(i<(nrow-1))
            printf(",");
        printf("\n");
    }
    printf("    ],\n");
    // X data points
    printf("  \"xdata\":\n");
    printf("    [\n");
    for (int i=0;i<nrow;i++){
        if (i<(nrow-1))
            printf("      %.5f,\n",x[i]);
        else
            printf("      %.5f\n",x[i]);
    }
    printf("    ],\n");
    // Y data points
    printf("  \"ydata\":\n");
    printf("    [\n");
    for (int i=0;i<nrow;i++){
        if (i<(nrow-1))
            printf("      %.5f,\n",y[i]);
        else
            printf("      %.5f\n",y[i]);
    }
    printf("    ],\n");
    // Z data points
    printf("  \"zdata\":\n");
    printf("    [\n");
     for (int i=0;i<nrow;i++){
        if (i<(nrow-1))
            printf("      %.5f,\n",z[i]);
        else
            printf("      %.5f\n",z[i]);
    }
    printf("    ]\n");

    printf("}\n");
    gsl_matrix_free (X);
    gsl_vector_free (zvec);
    gsl_vector_free (sigvec);
    gsl_vector_free (cvec);
    gsl_matrix_free (cov);

    return 0;
}

