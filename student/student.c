#include <math.h>
# include <stdlib.h>
# include <stdio.h>


/***************************************/

void nrerror(char error_text[]);

# define BUFSIZE 150

int main ( int argc, char * argv[]) {
    
    char filename[BUFSIZE] = {'\0'};
    char line[BUFSIZE] = {'\0'};
    double *x;
    double pval, aux, avg[2], avg_sq[2], var[2], avar;
    double estimator[2];
    double df, t;
    int ndata, n[2], * sample, s;
    
    int i, j;
    FILE *fptr;
    double betai(double a, double b, double x);
   
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
    if (!x) {
	fprintf (stderr, "Error opening %s.\n", filename);
	exit (1);
    }
    sample = (int*) calloc ( ndata+1, sizeof(int) );
    if (!sample) {
	fprintf (stderr, "Error opening %s.\n", filename);
	exit (1);
    }
    
    /* rewind and read in */
    rewind (fptr);
    ndata =  0;
    n[0] = n[1] =  0;
    while( fgets(line,BUFSIZE, fptr) ) {
	if ( sscanf ( line, "%lf %d", &aux, &s) == 2 ) {
	    x[ndata] = aux;
	    sample[ndata] = s;
	    ndata++;
	    n[s] ++;
	}
    }
    fclose (fptr);

    printf ( "\nread in %d points for sample 0, and %d points for sample 1.\n\n", n[0], n[1]);

    /* calculate the t statistics */
    avg[0]    = avg[1]    = 0;
    avg_sq[0] = avg_sq[1] = 0;
    
    for ( i=0; i<ndata; i++ ) {
	avg [ sample[i]] +=  x[i];
	avg_sq [ sample[i]] +=  x[i]*x[i];
    }
    for ( s=0; s < 2; s++ ) {
	estimator[s] = 	avg_sq[s] - avg[s]*avg[s]/n[s];
	avg[s] /= n[s];
	avg_sq[s] /= n[s];
	var[s] = avg_sq[s] - avg[s]*avg[s];
    }

    printf ( " %10s    %10s    %10s \n", "  ", "sample 1", "sample 2");
    printf ( " %10s    %10.2le    %10.2le\n",  "avg", avg[0], avg[1]  );
    printf ( " %10s    %10.2le    %10.2le\n",  "avg_sq", avg_sq[0], avg_sq[1]  );
    printf ( " %10s    %10.2le    %10.2le\n",  "stdev", 
	     sqrt (var[0]),   sqrt (var[1]) );
    
    printf ( "\n");
    
    df = n[0] + n[1] - 2;
    
    avar = ( estimator[0] +  estimator[1] )*(1.0/n[0]+1.0/n[1])/df;
    avar = sqrt (avar);
    t = (avg[0] - avg [1])/avar ;
    printf ( " t:  %8.2le   ( avg[0] - avg[1]:  %8.2le   avar:  %8.2le    estimator sum %8.2le  %8.2le %8.2le)  \n\n",
	     t, avg[0]-avg[1], avar,
	     estimator[0] +  estimator[1], (1.0/n[0]+1.0/n[1]), (1.0/n[0]+1.0/n[1])/df );

    pval = betai (0.5*df, 0.5, df/(df + t*t) );
    
    printf ( " pval:  %8.4le \n\n", pval);
    
    return 0;
} 

/************************************************/
/************************************************/

double gammln(double xx)

{

    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
			  24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

#define MAXIT 100
#define EPS 3.0e-15
#define FPMIN 1.0e-30

double betacf(double a, double b, double x)
/* Used by betai: Evaluates continued fraction for incomplete beta function by modified Lentz's method*/

{
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    qab=a+b;  
    qap=a+1.0;   
    qam=a-1.0;
    c=1.0;   
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) {
	m2=2*m;
	aa=m*(b-m)*x/((qam+m2)*(a+m2));
	d=1.0+aa*d;   
	if (fabs(d) < FPMIN) d=FPMIN;
	c=1.0+aa/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	h *= d*c;
	aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
	d=1.0+aa*d;  
	if (fabs(d) < FPMIN) d=FPMIN;
	c=1.0+aa/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	del=d*c;
	h *= del;
	if (fabs(del-1.0) < EPS) break;    
    }
    if (m > MAXIT) {
	fprintf ( stderr, "a or b too big, or MAXIT too small.\n");
	exit (1);
    }
    return h;
}

double betai(double a, double b, double x)
/* Returns the incomplete beta function Ix(a, b). */
{
    double bt;
    if (x < 0.0 || x > 1.0) {
	fprintf (stderr, "Bad x in routine betai");
    }
    if (x == 0.0 || x == 1.0) {
	bt=0.0;
    } else {
	bt = exp(gammln(a+b) - gammln(a) - gammln(b) + a*log(x) + b*log(1.0-x));
    }
    if (x < (a+1.0)/(a+b+2.0))   
	return bt*betacf(a,b,x)/a;
    else   
	return 1.0-bt*betacf(b,a,1.0-x)/b;   
}

