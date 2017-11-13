#include <stdlib.h>
#include <stdio.h>

#include "table.h"

#define MAXROW 100000

int main(int argc, char **argv) 
{
    double x[MAXROW];
    double y[MAXROW];
    char xcolname[64];
    char ycolname[64];
    int nrow = 0;
    int status = 0;

    sscanf(argv[1],"%s",xcolname);
    sscanf(argv[2],"%s",ycolname);

    status = read_xy(xcolname,ycolname,x,y,&nrow);
    if (status)
    {
        fprintf(stderr,"Error reading data table.\n");
        exit(1);
    }

    for (int i=0; i<nrow; i++)
        fprintf(stdout,"20%g %20g\n",x[i],y[i]);

    exit(EXIT_SUCCESS);
}
