# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>


int main ( int argc, char * argv[]) {
    FILE *fptr;
    double A0[20][20], B[20][20] = {{0.0}};
    double d, d_prev, incr, basic_incr;
    double sum, check;
    int i,j,  accepted, ctr;
    long step, no_steps;
    double distance ( double A0[20][], double B[20][] );
    srand48(9);
    if ( argc < 2 ) {
	fprintf (stderr, "Usage: %s <input matr> [<incr>] [<no_steps>]\n", argv[0]);
	exit (0);
    }
    basic_incr = 0.001;
    if ( argc > 2 ) {
	basic_incr  = atof (argv[2]);
    }
    no_steps = 100*400;
    if ( argc > 3) {
	no_steps  = atoi (argv[3])*400;
    }
    fptr = fopen ( argv[1], "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s\n", argv[1]);
	exit (0);
    }
    for (i=0; i<20; i++ ) {
	for (j=0; j<=i; j++ ) {
	    fscanf (fptr, "%lf", &A0[i][j] );
	    A0[j][i] = A0[i][j];
	}
    }
    fclose (fptr);
    
    for (i=0; i<20; i++ ) {
	B[i][i] = 1.0;;
    }

    d_prev  =  distance (A0, B);
    d = d_prev;
    for (step=1; step <= no_steps; step++ ) {

	/*sweep through all  positions */
	for (i=0; i<20; i++ ) {
	    for (j=0; j<i; j++ ) {
# if 0
		do {
		    i = lrand48() %20;
		} while ( !i);
		do {
		    j  = lrand48() %20;
		} while ( j >= i);
# endif
		incr = (1 - 2*drand48())*basic_incr;
		B[i][j] += incr;
		B[i][i] -= incr;
		B[j][j] -= incr;
		d =  distance (A0, B);
		accepted = ( 	B[i][j] >0 && B[i][i] > 0 && B[j][j]> 0) ;
		accepted =  accepted && ( ( d < d_prev )|| ( drand48() < exp ( -(d-d_prev)/0.0005 ) ) );
		if ( accepted ) {
		    d_prev = d;
		} else {
		    B[i][j] -= incr;
		    B[i][i] += incr;
		    B[j][j] += incr;
		}
	    }
	}
	if ( ! (step%1000) ) printf ( "%5ld  %8.3lf\n", step, d);
	
    }
    //return 0;
    for (i=0; i<20; i++ ) {
	for (j=0; j<i; j++ ) {
	    B[j][i] = B[i][j];
	}
    }
# if 1
    /* initial check: */
    printf ( "\n\n  raw frequency  matrix: ....\n" );
    for (i=0; i<20; i++ ) {
	sum = 0;
	for (j=0; j<20; j++ ) {
	    sum += A0[i][j];
	}
	printf ( "row: %3d  sum: %6.2lf \n", i, sum );
    }

    for (i=0; i<20; i++ ) {
	sum = 0;
	for (j=0; j<20; j++ ) {
	    sum += A0[j][i];
	}
	printf ( "col:  %3d  sum: %6.2lf \n", i, sum );
    }
    printf ( "\n\n  nearest transition matrix: ....\n" );
    /* final check: */
    for (i=0; i<20; i++ ) {
	sum = 0;
	for (j=0; j<20; j++ ) {
	    sum += B[i][j];
	}
	printf ( "row: %3d  sum: %6.2lf \n", i, sum );
    }

    for (i=0; i<20; i++ ) {
	sum = 0;
	for (j=0; j<20; j++ ) {
	    sum += B[j][i];
	}
	printf ( "col:  %3d  sum: %6.2lf \n", i, sum );
    }
# endif
    check = 0;
    ctr = 0;
    for (i=0; i<20; i++ ) {
	for (j=0; j<=i; j++ ) {
	    check += (A0[i][j] -  B[i][j] )*(A0[i][j] -  B[i][j] );
	    ctr++;
	    printf ( "%8.4lf, ", B[i][j] );
	}
	printf ("\n");
    }
    printf ("\n check: %8.3lf \n", sqrt(check/ctr));
    return 0;
}


/********************************************/
double distance ( double A0[20][20], double B[20][20] ) {

    int i, j;
    double dist;
    dist = 0;
    for (i=0; i<20; i++ ) {
	for (j=0; j<=i; j++ ) {
	    dist += fabs ( A0[i][j] - B[i][j]);
	}
    }

    return dist;
}
