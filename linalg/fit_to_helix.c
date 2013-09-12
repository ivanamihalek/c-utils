/* find transformation which takes one set of 3D vectors into another*/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define N  3
# define NO_VECTORS  21
# define BUFSIZE    150

void     dgels_ ( char * trans, int * no_rows, int * no_columns,  int * , double ** scratch, int * ,
		  double **A, int *, double * work, int * lwork, int *info);

int main ( int argc, char * argv[]) {

    if ( argc < 2 ) {
	fprintf ( stderr,
		  "Usage:  %s  <coord filename> (+ remember to fix the size of the problem in the src).\n",
		  argv[0]);
	exit (0);
    }

    double x[N][NO_VECTORS],  y[N][NO_VECTORS], T[N];
    double X[N][NO_VECTORS-1],    Y[N][NO_VECTORS-1], A[N][NO_VECTORS-1];
    int component;

    /* input vectors x and y*/
    FILE * fptr = NULL;
    char * filename = argv[1];
    char buffer [BUFSIZE], * retptr;
    fptr = fopen ( filename, "r") ;
    if ( !fptr ) {
	fprintf ( stderr, "could not open %s.\n", filename);
	exit (1);
    }

    int ctr, retval;
    for ( ctr =0; ctr < NO_VECTORS; ctr++ ) {
	memset ( buffer, 0, BUFSIZE);
	retptr = fgets ( buffer, BUFSIZE,  fptr);
	if ( retptr == buffer) {
	    printf ("%s", buffer);
	}
	retval = sscanf (buffer,  "%lf%lf%lf", &x[0][ctr], &x[1][ctr],  &x[2][ctr]);
	if ( retval != N) {
	    fprintf ( stderr, "error reading x%d from  %s  (%d) .\n", ctr, filename, retval);
	    exit (1);
	}
    }
    for ( ctr =0; ctr < NO_VECTORS; ctr++ ) {
	memset ( buffer, 0, BUFSIZE);
	retptr = fgets ( buffer, BUFSIZE,  fptr);
	if ( retptr == buffer) {
	    printf ("%s", buffer);
	}
	retval =  sscanf (buffer,"%lf%lf%lf", &y[0][ctr], &y[1][ctr], &y[2][ctr]);
	if ( retval != N) {
	    fprintf ( stderr, "error reading y%d from  %s.\n", ctr, filename);
	    exit (1);
	}
    }

    /* check: */
    printf (" The vectors read in: \n");
    for ( ctr =0; ctr < NO_VECTORS; ctr++ ) {
	printf ("\t x%1d   %10.4lf  %10.4lf  %10.4lf \n",  ctr, x[0][ctr], x[1][ctr], x[2][ctr]);
    }
    for ( ctr =0; ctr < NO_VECTORS; ctr++ ) {
	printf ("\t y%1d   %10.4lf  %10.4lf  %10.4lf \n",  ctr,  y[0][ctr], y[1][ctr], y[2][ctr]);
    }
     
    /* subtract the last x from the rest, and last y from the rest */
    /* and  write them as matrix - make it transpose to make LAPACK happy*/
    for ( ctr =0; ctr < NO_VECTORS-1; ctr++ ) {
	for ( component=0; component<N; component++) {
	    X[component][ctr] = x[component][ctr]-x[component][NO_VECTORS-1];
	    Y[component][ctr] = y[component][ctr]-y[component][NO_VECTORS-1];
	}
   }
    
    printf (" X-x%d: \n", NO_VECTORS-1);
    for ( component=0; component<N; component++) {
	for ( ctr =0; ctr < NO_VECTORS-1; ctr++ ) {
	    printf ("%8.3lf", X[component][ctr]);
	}
	printf ("\n");
    }
    printf ("\n");

    printf (" Y-y%d: \n", NO_VECTORS-1 );
    for ( component=0; component<N; component++) {
	for ( ctr =0; ctr < NO_VECTORS-1; ctr++ ) {
	    printf ("%8.3lf", Y[component][ctr]);
	}
	printf ("\n");
    }
    printf ("\n");

    
    /* solve the least squares problem */
    char trans= 'N';
    int info, ipiv[N];
    int lwork = 2*NO_VECTORS;
    double work[2*NO_VECTORS];
    double scratch[N][NO_VECTORS-1];
    int  n = N;
    int no_rows = NO_VECTORS-1, no_columns = N;
    
    memcpy (A[0], Y[0], N*(NO_VECTORS-1)*sizeof(double));
    memcpy (scratch[0], X[0], N*(NO_VECTORS-1)*sizeof(double));
    dgels_ ( &trans, &no_rows, &no_columns,  &n, &scratch, &no_rows, &A, &no_rows, work, &lwork, &info);
    printf (" info: %d\n", info);


    printf ("******************************************************\n");
    printf (" solution: \n" );
    for ( ctr =0; ctr < N; ctr++ ) {
	for ( component=0; component<N; component++) {
	    printf ("%10.3lf", A[ctr][component]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    
    /* is the solution orthogonal? */
    int i,j;
    double sum = 0;
    printf (" orthogonal?\n" );
    for ( i =0; i < N; i++ ) {
	for ( j =0; j < N; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<N; component++) {
		sum +=  A[i][component]*A[component][j];
	    }
	    printf ("%10.3lf", sum);
	}
	printf ("\n");
    }
    printf ("\n");
  
    
    /* is the solution ok? */
    printf (" solution ok?\n" );
    for ( i =0; i < N; i++ ) {
	for ( j =0; j < NO_VECTORS-1; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<N; component++) {
		sum +=  A[i][component]*X[component][j];
	    }
	    printf ("%8.3lf", sum);
	}
 	printf ("\n");
   }
    printf ("\n");
 
   

    /* find translation vector */
    printf ("******************************************************\n");
    printf (" translation vector\n" );
    for ( i =0; i < N; i++ ) {
	    sum = 0.0;
	    for ( component=0; component<N; component++) {
		sum +=  A[i][component]*x[component][NO_VECTORS-1];
	    }
	    T[i] = y[i][NO_VECTORS-1] - sum;
	    printf ("%10.3lf", T[i] );
	    
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
   

    /* reconstruct y */
    printf (" reconstruct y \n" );
    for ( j =0; j < NO_VECTORS; j++ ) {
	for ( i =0; i < N; i++ ) {
	    sum = 0.0;
	    for ( component=0; component<N; component++) {
		sum +=  A[i][component]*x[component][j];
	    }
	    printf ("%10.3lf", sum+T[i] );
	}
	printf ("\n");
    }
    printf ("\n");


    printf (" y: \n" );
    for ( ctr =0; ctr < NO_VECTORS; ctr++ ) {
	for ( component=0; component<N; component++) {
	    printf ("%10.3lf", y[component][ctr]);
	}
	printf ("\n");
    }
    printf ("\n");
   
    return 0;
}
