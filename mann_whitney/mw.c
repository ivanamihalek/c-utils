#include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>

# define BUFSIZE 150

# define STEPS_PER_POINT 100

/* from Wikipedia: */
/* The Mann-Whitney U test is a non-parametric statistical significance test for assessing whether */
/* the difference in medians between two samples of observations is statistically significant */
/* (whether the distributions of the samples overlap less than would be expected by chance). */
/* The null hypothesis is that medians are equal, that is, that the two samples are drawn from a */
/* single population.   */
/*    1. Choose the sample for which the ranks seem to be smaller (the choice is relevant only */
/*       to ease of computation). Call this "sample 1", and call the other sample "sample 2". */
/*    2. Taking each observation in sample 1, count the number of observations in sample 2 that */
/*       are smaller than it. */
/*    3. The total of these counts is U. */

/*    Here I am calculating the p-value through simulation. */

void nrerror(char error_text[]);


int main ( int argc, char * argv[]) {
    
    char filename[BUFSIZE] = {'\0'};
    char line[BUFSIZE] = {'\0'};
    double *x;
    double pval, aux;
    int ndata, xdata, ydata, * sample, s;
    int i, j, u;
    FILE *fptr;
    double p_value ( double *x, int ndata, int xdata, int u0);
   
   if ( argc < 2) {
	fprintf (stderr, "Usage: %s  <filename>.\n Input format: <score> <sample id>, where sample id = 0 or 1\n", argv[0]);
	exit (1);
    }
    sprintf (filename, "%s", argv[1] );
    
    fptr = fopen ( filename, "r");
    if ( !fptr) {
	fprintf (stderr, "Error opening %s.\n", filename);
	exit (1);
    }
    /*  count data */
    ndata = 0;
    while(fgets(line,BUFSIZE, fptr) )	ndata++;
    
    /* allocate space */
    x = (double*) calloc ( ndata+1, sizeof(double) );
    if (!x) nrerror ("Error allocating.\n");
    sample = (int*) calloc ( ndata+1, sizeof(int) );
    if (!sample) nrerror ("Error allocating.\n");
    
    /* rewind and read in */
    rewind (fptr);
    ndata =  0;
    xdata = ydata = 0;
    while( fgets(line,BUFSIZE, fptr) ) {
	if ( sscanf ( line, "%lf %d", &aux, &s) == 2 ) {
	    x[ndata] = aux;
	    sample[ndata] = s;
	    ndata++;
	    if ( s ) {
		xdata++;
	    } else {
		ydata++;
	    }
	}
    }
    fclose (fptr);

    printf ( "\nread in %d points for sample 1, and %d points for sample 0.\n\n", xdata, ydata);

    /* calculate the U statistics */
    u = 0;
    for ( i=0; i<ndata; i++ ) {
	if ( ! sample[i] ) continue;
	for ( j=0; j<ndata; j++ ) {
	    if (  sample[j] ) continue;
	    u += ( x[j] < x[i] );
	}
    }
    printf ( " u:  %d \n\n", u);

    pval = p_value (x, ndata, xdata, u);
    printf ( " pval:  %8.4le \n\n", pval);
    
    return 0;
} 

/***************************************/

double p_value ( double * x, int ndata, int xdata, int u0) {

    double pval = 0, prob;
    int *sample;
    int init, ctr, i, j, round;
    int u;
    int no_rounds = STEPS_PER_POINT*ndata;
    sample = (int*) calloc ( ndata, sizeof(int) );
    if (!sample) nrerror ("Error allocating.\n");

    if ( ! ndata )  nrerror ("Error: empty sample.\n");

    pval = 0.0;
    prob = (double)xdata/ndata;

    for (round=0; round < no_rounds; round ++ ){
	
	init = lrand48()%ndata;
	ctr = 0;
	memset (sample, 0, ndata*sizeof(int) );
	
	while ( ctr < xdata ) {
	    for (i=0; i<ndata && ctr < xdata; i++){
		j = (init+i)%ndata;
		if ( sample[j]) continue;
		if (drand48() < prob) {
		    sample[j] = 1;
		    ctr ++;
		}
	    }
	}
	u = 0;
	for ( i=0; i<ndata; i++ ) {
	    if ( ! sample[i] ) continue;
	    for ( j=0; j<ndata; j++ ) {
		if (  sample[j] ) continue;
		u += ( x[j] < x[i] );
	    }
	}
	if ( u <= u0) pval += 1.0;
	
    }
    
    return pval/no_rounds;

    
}



/***************************************/
 void nrerror(char error_text[]){
     fprintf ( stderr, "%s\n", error_text);
     exit(1);
}
