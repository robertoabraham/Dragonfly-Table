int read_xy(char *xcolname, char *ycolname, double *x, double *y, int *nrow);
int read_xyz(char *xcolname, char *ycolname, char *zcolname, double *x, double *y, double *z, int *nrow);
int read_xyzs(char *xcolname, char *ycolname, char *zcolname, char *scolname, double *x, double *y, double *z, double *s, int *nrow);
int isNumeric (const char * s);
