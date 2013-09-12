# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"

void dgeqrf_ (int *M, int *, double *A, int *LDA, double * TAU, double * WORK, int * LWORK, int *INFO );
void dorgrq_ (int *M, int *, int *, double *A,  int *LDA,             double *TAU, double *   WORK,  int *LWORK,  int *INFO);
void dgeqp3_ (int *M, int *,         double *A, int *LDA, int * JPVT, double * TAU, double *WORK, int * LWORK, int *INFO );
    
int main ( int argc, char * argv[]) {
    int n = 3, info= 0;
    double  work[100];
    int lwork = 100, jpvt[3] = {0};
    double **Q, *tau,  R[3][3] = {{0.0}};
    int h, i,j,k,  component;
    double sum, theta;
    double rQ[3][3], rQnew[3][3], H[3][3], v[3];

    tau = emalloc  (3*sizeof(double));
    Q = dmatrix (3, 3);
    
    if ( 0 ) {
	Q[0][0] = Q[1][1] = Q[2][2] = 1.0;
    } else {
	theta = M_PI*0.1;
	Q[0][0] = cos(theta);  Q[0][1] = -sin(theta); /* fortran transposes */
	Q[1][0] = sin(theta);  Q[1][1] =  cos(theta);
	Q[2][2] = 1.0;
    }
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    printf ("%10.3lf",  Q[j][i]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    printf ("initial orthogonality\n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<n; component++) {
		sum +=  Q[i][component]*Q[j][component];
	    }
	    printf ("%10.3lf", sum);
	}
	printf ("\n");
    }
    printf ("\n");

    
    if ( 1) {
	dgeqrf_ ( &n, &n, Q[0], &n, tau, work, &lwork, &info);
	printf ("dgeqrf info: %d\n", info);
    } else {
	dgeqp3_ ( &n, &n,  Q[0], &n, jpvt, tau, work, &lwork, &info);
	printf ("dgeqp3 info: %d\n", info);
    }
 
    
    printf ("Q after dgeqrf\n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    printf ("%10.3lf",  Q[i][j]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    /*extract R*/
    for ( i =0; i < n; i++ ) {
	for ( j =i; j < n; j++ ) {
	    R[i][j] = Q[j][i];
	}
    }
    printf ("R after dgeqrf\n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    printf ("%10.3lf",  R[i][j]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    
    printf ("tau after dgeqrf\n");
    for ( i =0; i < n; i++ ) {
	printf ("%10.3lf",  tau[i]);
    }
    printf ("\n");
    printf ("******************************************************\n");
    
    /* reconstruct Q: */
    //dorgrq_ ( &n, &n, &n,  Q[0], &n, tau, work, &lwork, &info);
    // printf ("dorgrq info: %d\n", info);
    // printf ("******************************************************\n");
    
    memset( rQ[0], 0, n*n*sizeof(double));
    rQ[0][0] = rQ[1][1] = rQ[2][2] = 1.0;
    for ( h =0; h < n; h++ ) {
	/* find vh*/
	for ( i=0; i<h; i++ ) v[i] = 0.0;
	v[h] = 1.0;
	for ( i=h+1; i<n; i++ ) v[i] = Q[h][i];
	
	/* find Hh */
	for ( i =0; i < n; i++ ) {
	    H[i][i] = 1.0 -tau[h]*v[i]*v[i];
	    for ( j =i+1; j < n; j++ ) {
		H[i][j] = H[j][i] = -tau[h]*v[i]*v[j];
	    }
	}
	
	/* multiply rQ by Hi */
	for ( i =0; i < n; i++ ) {
	    for ( j =0; j < n; j++ ) {
		rQnew[i][j] = 0.0;
		for ( k =0; k < n; k++ ) {
		    rQnew[i][j] += rQ[i][k]*H[k][j];
		}
		
	    }
	}
	memcpy ( rQ[0], rQnew[0], n*n*sizeof(double));
    }

    /* to get as close as possible to the original matrix, require that diagonals
       in R be positive (in the limiting case whna the input matrix is already
       orthogonal, R should be I */
    
    for ( i =0; i < n; i++ ) {
	if ( R[i][i] < 0 ) {
	    R[i][i] *= -1;
	    for ( j =0; j < n; j++ ) {
		rQ[j][i] *= -1;
	    }
	    
	} 
    }
  

    
    printf ("Q reconstructed    \n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    printf ("%10.3lf",  rQ[i][j]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    printf ("final orthogonality\n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<n; component++) {
		sum +=  rQ[component][i]*rQ[component][j];
	    }
	    printf ("%10.3lf", sum);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    printf ("QRproduct\n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<n; component++) {
		sum +=  rQ[i][component]*R[component][j];
	    }
	    printf ("%10.3lf", sum);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    
    return 0;
}
