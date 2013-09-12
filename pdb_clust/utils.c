#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>

#include <sys/types.h>

void error (int errno, char *errstr) {
    fprintf (stderr, errstr);
    exit (errno);
}

void usage(char *  use[])
{

	if (use != NULL) {
		(void) fprintf(stdout, "\n\t%s\n\n", *use);
		while (*++use != NULL)
			(void) fprintf(stdout, "\t%s\n", *use);
		(void) fprintf(stdout, "\n");
	}

	return;
}


void * emalloc(size_t	size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	exit (1);
    }

    return ptr;
}



FILE * efopen(char * name, char * mode)
{

    FILE * fp;


    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,  
	      "efopen: can't open \"%s\" for \"%s\"\n", name, mode);
	exit (1);
    }

    return fp;

}




/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] 
 */
char **chmatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    char **m;
  
    /* allocate pointers to rows */
    m=(char **) malloc((size_t)((nrow+1)*sizeof(char*)));
    if (!m)  {
	fprintf (stderr,"allocation failure 1 in matrix()");
    }
    m += 1;
    m -= nrl;
  
    /* allocate rows and set pointers to them */
    m[nrl]=(char *) calloc( nrow*ncol+1, sizeof(char));
    if (!m[nrl]) {
	fprintf (stderr,"allocation failure 2 in matrix()for char");
    }
    m[nrl] += 1;
    m[nrl] -= ncl;
  
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
    /* return pointer to array of pointers to rows */ 
    return m; 
}



/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] 
 */
int **imatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    int **m;

    /* allocate pointers to rows */
    m=(int **) malloc((size_t)((nrow+1)*sizeof(int*)));
    if (!m) fprintf (stderr,"allocation failure 1 in matrix()");
    m += 1;
    m -= nrl;


    /* allocate rows and set pointers to them */
    m[nrl]=(int *) calloc(nrow*ncol+1,sizeof(int));
    if (!m[nrl]) fprintf (stderr,"allocation failure 2 in matrix() for int");
    m[nrl] += 1;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}


/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
 */
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    double **m;

    /* allocate pointers to rows */
    m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
    if (!m) fprintf (stderr,"allocation failure 1 in matrix()");
    m += 1;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(double *) calloc( nrow*ncol+1,sizeof(double));
    if (!m[nrl]) fprintf (stderr,"allocation failure 2 in matrix()for double");
    m[nrl] += 1;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}



/* allocate a float square matrix 
 */
float **fmatrix(long dimension)
{
    long i;
    float **m;

    /* allocate pointers to rows */
    m=(float **) malloc((size_t)((dimension)*sizeof(float*)));
    if (!m) fprintf (stderr,"allocation failure 1 in matrix()");

    for(i=0;i<dimension;i++)
    {
	m[i] = (float *) calloc(dimension,sizeof(float));
	if (!m[i]) fprintf (stderr,"allocation failure 1 in matrix()");
    }

    /* return pointer to array of pointers to rows */
    return m;
}



/* free a char matrix allocated by chmatrix() 
 */
void free_chmatrix(char **m, long nrl, long nrh, long ncl, long nch)
{
    free( (m[nrl]+ncl-1));
    free( (m+nrl-1));
}



/* free an int matrix allocated by imatrix() 
 */
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
    free( (m[nrl]+ncl-1));
    free( (m+nrl-1));
}



/* free a double matrix allocated by dmatrix() 
 */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free( (m[nrl]+ncl-1));
    free( (m+nrl-1));
}




/* sort array according to the score in the other */
/* I couldn't declare pos_cmp within array_qsort  bcs it   crashed on mac */

double * score_array;

int pos_cmp (const void * a0, const void * b0) {
    
    int * a= (int*) a0;
    int * b= (int*)b0;
    if ( score_array[*a] > score_array[*b]) {
	return 1;
    }
    if ( score_array[*a] < score_array[*b]) {
	return -1;
    }
    return 0;
}


int array_qsort (int * sorted_pos, double * sa, int sequence_length ) {
    /* position comparison function */
    score_array = sa;

    qsort (sorted_pos+1, sequence_length, sizeof(int), pos_cmp);

    return 0;
}
