/* find transformation which takes one set of 3D vectors into another*/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define N 3
# define BUFSIZE 150

 void   dgesv_ ( int *,  int *, double ** , int *, int *, double **, int *, int *);
 
int main ( int argc, char * argv[]) {

    if ( argc < 2 ) {
	fprintf ( stderr, "Usage:  %s  <coord filename>.\n", argv[0]);
	exit (0);
    }

    double x[N][N+1],  y[N][N+1], T[N];
    double X[N][N],    Y[N][N], A[N][N];
    int component;

    /* input vectors x and y, there should be 4 of each*/
    FILE * fptr = NULL;
    char * filename = argv[1];
    char buffer [BUFSIZE], * retptr;
    fptr = fopen ( filename, "r") ;
    if ( !fptr ) {
	fprintf ( stderr, "could not open %s.\n", filename);
	exit (1);
    }

    int ctr, retval;
    for ( ctr =0; ctr <= N; ctr++ ) {
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
    for ( ctr =0; ctr <= N; ctr++ ) {
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
    for ( ctr =0; ctr <= N; ctr++ ) {
	printf ("\t x%1d   %10.4lf  %10.4lf  %10.4lf \n",  ctr, x[0][ctr], x[1][ctr], x[2][ctr]);
    }
    for ( ctr =0; ctr <= N; ctr++ ) {
	printf ("\t y%1d   %10.4lf  %10.4lf  %10.4lf \n",  ctr,  y[0][ctr], y[1][ctr], y[2][ctr]);
    }
     
    /* subtract the last x from the rest, and last y from the rest */
    /* and  write them as matrix */
    for ( ctr =0; ctr < N; ctr++ ) {
	for ( component=0; component<N; component++) {
	    X[component][ctr] = x[component][ctr]-x[component][N];
	    Y[component][ctr] = y[component][ctr]-y[component][N];
	}
    }
    
    printf (" X-x3: \n" );
    for ( component=0; component<N; component++) {
	for ( ctr =0; ctr < N; ctr++ ) {
	    printf ("%10.3lf", X[component][ctr]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf (" Y-y3: \n" );
    for ( component=0; component<N; component++) {
	for ( ctr =0; ctr < N; ctr++ ) {
	    printf ("%10.3lf", Y[component][ctr]);
	}
	printf ("\n");
    }
    printf ("\n");
    
   /* solve the system X^TA^T = Y^T for A */
    /* I don't need to transpose  - fortran will kindly do that for me */
    int info, ipiv[N];
    double scratch[N][N];
    int n = N;
    memcpy (A[0], Y[0], N*N*sizeof(double));
    memcpy (scratch[0], X[0], N*N*sizeof(double));
    dgesv_ ( &n,  &n, &scratch, &n, ipiv, &A, &n, &info);
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

    /* is the solution ok? */
    int i,j;
    double sum = 0;
    printf (" solution ok?\n" );
    for ( i =0; i < N; i++ ) {
	for ( j =0; j < N; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<N; component++) {
		sum +=  A[i][component]*X[component][j];
	    }
	    printf ("%10.3lf", sum);
	}
 	printf ("\n");
   }
    printf ("\n");
 
    
    /* is the solution orthogonal? */
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
  

    /* find translation vector */
    printf ("******************************************************\n");
    printf (" translation vector\n" );
    for ( i =0; i < N; i++ ) {
	    sum = 0.0;
	    for ( component=0; component<N; component++) {
		sum +=  A[i][component]*x[component][N];
	    }
	    T[i] = y[i][N] - sum;
	    printf ("%10.3lf", T[i] );
	    
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
   

    /* reconstruct y */
    printf (" reconstruct y \n" );
    for ( j =0; j < N; j++ ) {
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
    for ( ctr =0; ctr < N; ctr++ ) {
	for ( component=0; component<N; component++) {
	    printf ("%10.3lf", y[component][ctr]);
	}
	printf ("\n");
    }
    printf ("\n");
   
    return 0;
}
