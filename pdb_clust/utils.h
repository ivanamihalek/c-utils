# ifndef _UTILS_H
# define _UTILS_H
# include <stdio.h>

void usage(char *  use[]);
void error (int errno, char *errstr);
void * emalloc(int	size);
FILE * efopen(char * name, char * mode);
char **chmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
float **fmatrix(long dimension);
void free_chmatrix(char **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

int array_qsort (int * sorted_pos, double * sa, int sequence_length );
# endif
