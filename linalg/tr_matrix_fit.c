/* find transition matrix by least suqres fit to relative frequency matrix*/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define N  3
# define M  N*N
# define BUFSIZE    150

void     dgelss_ (  int * no_rows, int * no_columns,  int * , double ** scratch, int * ,
		  double **A, int *, double * work, int * lwork, int *info);

int main ( int argc, char * argv[]) {

    if ( argc < 2 ) {
	fprintf ( stderr,
		  "Usage:  %s  <freq filename>.\n",
		  argv[0]);
	exit (0);
    }

    double map_p[N][N],  map_f[N][N];
    double f0[N*N];
    double f[N*(N+1)/2];
    int i,j;
    int ctr, upper_ctr;
    int component;
    double sum = 0;

    /* input target frequencies*/
    f0[N*0+0] = 5; f0[N*0+1] = 3;   f0[N*0+2] = 3;   //f0[N*0+3] = 1;  
    f0[N*1+0] = 4; f0[N*1+1] = 6;   f0[N*1+2] = 3;   //f0[N*1+3] = 2;
    f0[N*2+0] = 2; f0[N*2+1] = 3;   f0[N*2+2] = 4;   //f0[N*2+3] = 1;
    //0[N*3+0] = 1; f0[N*3+1] = 2;   f0[N*3+2] = 1;   f0[N*3+3] = 8;

    /* normalize rowwise: */
    for (i=0; i<N; i++) {
	sum = 0;
	for (j=0; j<N; j++) {
	    sum += f0[N*i+j];
	}
	for (j=0; j<N; j++) {
	    f0[N*i+j]/= sum;
	}
    }
    
    /*print out */
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    printf (" %8.2lf ",  f0[N*i+j]);
	}
	printf ("\n");
    }
    
    /* replace diag elents in f by 1-fii */
    for (i=0; i<N; i++) {
	f0[N*i+i] = 1 - f0[N*i+i];
    }
    
    
    /*out: */
    printf ("\n\n");
    for ( ctr=0; ctr< N*N; ctr++ ) {
	printf (" %8.2lf ",  f0[ctr]);
    }
    printf ("\n\n");



    /*****************************************************************/
    /*****************************************************************/
    /* the map for f (the RHS) */
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    map_f[i][j] = -1;
	}
    }
    ctr = -1;
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    ctr++;
	    map_f [i][j]  = ctr;
	}
    }
    
     /* the map for p (the LHS) */
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    map_p[i][j] = -1;
	}
    }
    ctr = -1;
    for (i=0; i<N; i++) {
	for (j=i+1; j<N; j++) {
	    ctr++;
	    map_p [i][j] = map_p [j][i] = ctr;
	}
    }
    
    
    /* construct matrix A - imposes  sum=1 condition */
    /* one row as a matrix */
    double A[N*N][N*(N-1)/2] = {{0}};
    int k, l;
    int row, column;
    for ( i =0; i < N; i++ ) {
	for ( j = 0; j < N; j++ ) {
	    row = map_f[i][j];
	    if ( row < 0 ) {
		fprintf ( stderr, "map f:  %d  %d \n", i,k);
		exit (1);
	    }
	    if ( i== j ) {
		for ( k =0; k< N; k++ ) {
		    if ( k==i ) continue;
		    column = map_p [i][k];
		    if ( column < 0 ) {
			fprintf ( stderr, "map p:  %d  %d \n", i,k);
			exit (1);
		    }
		    A[row][column] = 1;
		}
	    } else {
		column = map_p [i][j];
		if ( column < 0 ) {
		    fprintf ( stderr, "map p:  %d  %d \n", i,k);
		    exit (1);
		}
		A[row][column] = 1;
	    }
	}
   }

   printf ("\n\n");
   for ( i =0; i < N*N; i++ ) {
	for ( j = 0; j < N*(N-1)/2; j++ ) {
	    printf (" %1d ",  (int) A[i][j]);	    
	}
	printf ("\n");
   }
   printf ("\n\n");
   
   /* solve the least squares problem */
   char trans= 'T';
   int info, ipiv[N];
   int lwork = 2*M;
   double work[2*M];
   double scratch[M];
   int  nrhs = 1;
   int no_rows = M, no_columns = N*(N-1)/2;
    
   memcpy (scratch, f, N*N*sizeof(double));
   dgelss_ ( &no_rows, &no_columns,  &nrhs,  A[0], &no_rows, &scratch, &no_rows, work, &lwork, &info);
   printf (" info: %d\n", info);


   printf ("******************************************************\n");
   printf (" solution: \n" );
   for ( j = 0; j < N*(N-1)/2; j++ ) {
       printf(" %8.2f ", scratch[j]);
   }

   printf ("\n");
   printf ("******************************************************\n");
    
# if 0    
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
# endif
    return 0;
}
