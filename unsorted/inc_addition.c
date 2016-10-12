# include <stdio.h>

# define MAXL 15


main () {
    int L, i, j, k, l, sum;
    int Lold, jlo , llo;
    int partial_sum[MAXL][MAXL], partl;
    
    sum = 0;
    partl = 0;
    L = 4;
    for (i = 0; i <L-1; i++ ) {
	for (j = i+1; j <L; j++ ) {
	    for (k = 0; k <L-1; k++ ) {
		for (l = k+1;  l<L; l++ ) {
		    sum += 1;
		}
	    }
	   
	    partial_sum[i][j] = 1;
	}
    }

    printf (" L: %d    sum: %d \n", L, sum) ;

    Lold = L;
    L+=2;
   
    for (k = 0; k <L-1; k++ ) {
	llo = (k<Lold) ? Lold : k+1;
	for (l = llo;  l<L; l++ ) {
	    for (i = 0; i <Lold-1; i++ ) {
		for (j = i+1; j <Lold; j++ ) {
		    sum += partial_sum[i][j] ;
		}
	    }
	}
    }
    
    for (i = 0; i <L-1; i++ ) {
	jlo = (i<Lold) ? Lold : i+1;
	for (j = jlo; j <L; j++ ) {
	    for (k = 0; k <L-1; k++ ) {
		for (l = k+1;  l<L; l++ ) {
		    sum += 1;
		}
	    }
	}
    }
    printf (" L: %d    sum: %d \n", L, sum) ;

   
    
}



# if 0
main () {
    int L, i, j, k, sum;
    int Lold, jlo;
    
    sum = 0;
    L = 4;
    for (i = 0; i <L-1; i++ ) {
	for (j = i+1; j <L; j++ ) {
	    sum += 1;
	}
    }

    printf (" L: %d    sum: %d \n", L, sum) ;

    Lold = L;
    L+=2;

    for (i = 0; i <L-1; i++ ) {
	jlo = ( i<Lold) ? Lold : i+1;
	for (j = jlo; j <L; j++ ) {
	    sum += 1;
	}
    }
   
    printf (" L: %d    sum: %d \n", L, sum) ;

   
    
}
# endif
