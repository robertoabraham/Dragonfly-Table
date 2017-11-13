#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "table.h"

int is_numeric (const char * s)
{
    char * p;

    if (s == NULL || *s == '\0' || isspace(*s))
      return 0;
    strtod (s, &p);
    return *p == '\0';
}


int read_xyzs(char *xcolname, char *ycolname, char *zcolname, char *scolname, double *x, double *y, double *z, double *s, int *nrow)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t nread;
    int ncol = 0;
    int xcol = -1;
    int ycol = -1;
    int zcol = -1;
    int scol = -1;
    char *dummy1,*dummy2,*keyword,*description;
    const char delims[]=" \n";
    char **tok;

    *nrow = 0;
    while ((nread = getline(&line, &len, stdin)) != -1) {

        if (!strncmp(&line[0],"#",1))
        { 
            // HEADER
       
            // Skip comments
            if (!strncmp(&line[1],"!",1)) continue;
            
            // Determine columns numbers corresponding to desired column names.
            dummy1      = strtok(line,delims);
            dummy2      = strtok(NULL,delims);
            keyword     = strtok(NULL,delims);
            description = strtok(NULL,delims);
            if (!strncmp(xcolname,keyword,64))
                xcol = ncol;
            if (!strncmp(ycolname,keyword,64))
                ycol = ncol;
            if (!strncmp(zcolname,keyword,64))
                zcol = ncol;
            if (!strncmp(scolname,keyword,64))
                scol = ncol;
             ncol++;
        }
        else
        {
            // DATA
            
            if (*nrow==0) {
                if (xcol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",xcolname);
                    return(1);
                }
                if (ycol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",ycolname);
                    return(1);
                }
                if (zcol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",zcolname);
                    return(1);
                }
                if (scol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",scolname);
                    return(1);
                }
                 tok = malloc(ncol*sizeof(char *));
            }

            *tok = strtok(line,delims);
            for (int i=1;i<ncol;i++)
                *(tok + i) = strtok(NULL,delims);

            x[*nrow] = atof(*(tok + xcol));
            y[*nrow] = atof(*(tok + ycol));
            z[*nrow] = atof(*(tok + zcol));
            s[*nrow] = atof(*(tok + scol));

            (*nrow)++; 
        }
    }

    free(line);
    free(tok);

    return (0);
}


int read_xyz(char *xcolname, char *ycolname, char *zcolname, double *x, double *y, double *z, int *nrow)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t nread;
    int ncol = 0;
    int xcol = -1;
    int ycol = -1;
    int zcol = -1;
    char *dummy1,*dummy2,*keyword,*description;
    const char delims[]=" \n";
    char **tok;

    *nrow = 0;
    while ((nread = getline(&line, &len, stdin)) != -1) {

        if (!strncmp(&line[0],"#",1))
        { 
            // HEADER
       
            // Skip comments
            if (!strncmp(&line[1],"!",1)) continue;
            
            // Determine columns numbers corresponding to desired column names.
            dummy1      = strtok(line,delims);
            dummy2      = strtok(NULL,delims);
            keyword     = strtok(NULL,delims);
            description = strtok(NULL,delims);
            if (!strncmp(xcolname,keyword,64))
                xcol = ncol;
            if (!strncmp(ycolname,keyword,64))
                ycol = ncol;
            if (!strncmp(zcolname,keyword,64))
                zcol = ncol;
            ncol++;
        }
        else
        {
            // DATA
            
            if (*nrow==0) {
                if (xcol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",xcolname);
                    return(1);
                }
                if (ycol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",ycolname);
                    return(1);
                }
                if (zcol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",zcolname);
                    return(1);
                }
                tok = malloc(ncol*sizeof(char *));
            }

            *tok = strtok(line,delims);
            for (int i=1;i<ncol;i++)
                *(tok + i) = strtok(NULL,delims);

            x[*nrow] = atof(*(tok + xcol));
            y[*nrow] = atof(*(tok + ycol));
            z[*nrow] = atof(*(tok + zcol));

            (*nrow)++; 
        }
    }

    free(line);
    free(tok);

    return (0);
}



int read_xy(char *xcolname, char *ycolname, double *x, double *y, int *nrow)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t nread;
    int ncol = 0;
    int xcol = -1;
    int ycol = -1;
    char *dummy1,*dummy2,*keyword,*description;
    const char delims[]=" \n";
    char **tok;

    *nrow = 0;
    while ((nread = getline(&line, &len, stdin)) != -1) {

        if (!strncmp(&line[0],"#",1))
        { 
            // HEADER
       
            // Skip comments
            if (!strncmp(&line[1],"!",1)) continue;
            
            // Determine columns numbers corresponding to
            // desired column names.
            dummy1      = strtok(line,delims);
            dummy2      = strtok(NULL,delims);
            keyword     = strtok(NULL,delims);
            description = strtok(NULL,delims);
            if (!strncmp(xcolname,keyword,64))
                xcol = ncol;
            if (!strncmp(ycolname,keyword,64))
                ycol = ncol;
            ncol++;
        }
        else
        {
            // DATA
            
            if (*nrow==0) {
                if (xcol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",xcolname);
                    return(1);
                }
                if (ycol < 0){
                    fprintf(stderr,"Keyword %s not found.\n",ycolname);
                    return(1);
                }
                tok = malloc(ncol*sizeof(char *));
            }

            *tok = strtok(line,delims);
            for (int i=1;i<ncol;i++)
                *(tok + i) = strtok(NULL,delims);

            /*
            if (!is_numeric((*(tok + xcol)))){
                fprintf(stderr,"Non-numeric entry in X column: %s\n",(*(tok + xcol)));
                return(1);
            }
            if (!is_numeric((*(tok + ycol)))){
                fprintf(stderr,"Non-numeric entry in Y column: %s\n",(*(tok + ycol)));
                return(1);
            }
            */

            x[*nrow] = atof(*(tok + xcol));
            y[*nrow] = atof(*(tok + ycol));

            //printf("%g %g\n",x[*nrow],y[*nrow]);
            
            (*nrow)++; 
        }
    }

    free(line);
    free(tok);

    return (0);
}
