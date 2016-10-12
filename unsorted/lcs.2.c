# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
double **dmatrix(int rows, int columns);
char **chmatrix(int rows, int columns);
void * emalloc(int	size);
void free_matrix(void **m);

# define FAR_FAR_AWAY 1000000

int main ( int argc, char * argv[]) {

    char A[] = "nematode knowledge";
    char B[] = "empty bottle";
    int m, n;
    int i,j;
    double **F;
 
    m = strlen (A);
    n = strlen (B);
    
    /* allocate F */
    F = dmatrix ( m+1, n+1);

    
    /* fill the table */
    for (i=m-1; i>= 0; i--) {
	for (j=n-1; j>=0; j--) {
	    if ( A[i] == B[j] ) {
		F[i][j] = F[i+1][j+1] + 1;
	    } else if (F[i+1][j] >= F[i][j+1] ) {
		F[i][j] = F[i+1][j];
	    } else {
		F[i][j] = F[i][j+1];
	    }
	   
	}
    }

    /*retrace*/
    i = 0;
    j = 0;
    while ( i<m &&  j<n ) {
	if ( A[i] == B[j] ) {
	    printf ( "%d  %d  %c %c \n", i+1, j+1, A[i], B[j]);
	    i++;
	    j++;
	} else if (F[i+1][j] >= F[i][j+1]) {
	    i++;
	} else {
	    j++;
	}
    }

    /* free */ 
    free_matrix ((void**) F);

    
    return 0;
}
/**********************************************/
double **dmatrix(int rows, int columns){
    double **m;
    int i;
        /* allocate pointers to rows */
    m=(double **) malloc(rows*sizeof(double*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(double *) calloc( rows*columns, sizeof(double));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

char **chmatrix(int rows, int columns){
    char **m;
    int i;
        /* allocate pointers to rows */
    m=(char **) malloc(rows*sizeof(char*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(char *) calloc( rows*columns, sizeof(char));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

void * emalloc(int	size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}


/* free a  matrix  */
void free_matrix(void **m)
{
    free(m[0]);
    free(m);
}


