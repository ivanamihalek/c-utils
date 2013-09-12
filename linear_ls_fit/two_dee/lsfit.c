# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define BUFSIZE 150

void nrerror(char error_text[]);

int main ( int argc, char * argv[]) {
    char filename[BUFSIZE] = {'\0'};
    char line[BUFSIZE] = {'\0'};
    double *x, *y;
    double a, b, siga, sigb, r, epsilon, q;
    int ndata;
    FILE *fptr;
    void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
	     double *b, double *siga, double *sigb, double *r,double * epsilon,  double *q);
    
   if ( argc < 2) {
	fprintf (stderr, "Usage: %s  <filename>.\n", argv[0]);
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
    y = (double*) calloc ( ndata+1, sizeof(double) );
    if (!y) nrerror ("Error allocating.\n");
    /* rewind and read in */
    rewind (fptr);
    ndata = 0;
    while( fgets(line,BUFSIZE, fptr) ) {
	ndata++;
	sscanf ( line, "%lf %lf", &x[ndata], &y[ndata]);
    }
    fclose (fptr);
    
    fit(x, y, ndata, NULL, 0, &a, &b,  &siga, &sigb, &r, &epsilon, &q);
    // printf ( "%8.3lf  %8.3lf  %8.3lf  %8.3lf  %8.3lf  %8.3le  \n",
//	     a,  siga, b,sigb, r, epsilon ? log(epsilon):epsilon);
    printf ( "%%%8s  %8s  %8s  %8s  %8s  %8s  \n",
	     "a",  "siga", "b","sigb", "r", "epsilon");
    printf ( " %8.3lf  %8.3lf  %8.3lf  %8.3lf  %8.3lf  %8.3le  \n",
	     a,  siga, b,sigb, r, epsilon);
    return 0;
}


/***********************************************************************/
void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
	 double *b, double *siga, double *sigb, double *r2, double *epsilon, double *q)
{
    double gammq(double a, double x);
    double erfcc(double x);
    double erffc(double x);
    int i;
    double  wt,t,sxoss, syoss, sx=0.0,sy=0.0,st2=0.0,ss,sigdat, chi2;
    double u, su2= 0.0;
    *b=0.0;
    if (mwt) {   // Accumulate sums ...
	ss=0.0;
	for (i=1;i<=ndata;i++) {   // ...with weights
	    wt=1.0/sqrt(sig[i]);
	    ss += wt;
	    sx += x[i]*wt;
	    sy += y[i]*wt;
	}
    } else {
	for (i=1;i<=ndata;i++) {    //...or without weights.
	    sx += x[i];
	    sy += y[i];
	}
	ss=ndata;
    }
    sxoss=sx/ss;
    syoss=sy/ss;
    /* solve for a and b */
    if (mwt) {
	for (i=1;i<=ndata;i++) {
	    t    =  (x[i]-sxoss)/sig[i];
	    st2 +=   t*t;
	    *a  +=   t*y[i]/sig[i];
	}
    } else {
	for (i=1;i<=ndata;i++) {
	    t    =  x[i] - sxoss;
	    u    =  y[i] - syoss;
	    st2 +=  t*t;
	    su2 +=  u*u;
	    *a  +=  t*y[i];
	}
    }
    *a /= st2;    //Solve for a, b, siga, and sig b.
    *b  = (sy-sx*(*a))/ss;
    *sigb = sqrt ( ( 1.0 + sx*sx/(ss*st2) )/ss);
    *siga = sqrt (1.0/st2);

    chi2 = 0.0;    //Calculate chi2.
    *q=1.0;
    *epsilon = 0.0;

    if (mwt) {
	for (i=1;i<=ndata;i++)
	    chi2 += sqrt((y[i]-(*a)-(*b)*x[i])/sig[i]);
    } else {
	for (i=1;i<=ndata;i++) {
	    t  = (y[i]-(*b)-(*a)*x[i]);
	    chi2 += t*t;
	}
	sigdat  = sqrt((chi2)/(ndata-2));    //For unweighted data evaluate typ
	*siga  *= sigdat;    //ical sig using chi2, and ad-
	*sigb  *= sigdat;    //just the standard deviations.
	*r2     = sqrt (1 - chi2/su2);
	*epsilon = erffc( (*r2)*sqrt(ndata/2));
	//*epsilon = 0.0;
    }
    if (ndata>2) *q=gammq(0.5*(ndata-2),0.5*(chi2));   // Equation (15.2.12).
 
}


/*************************************************************/
//Returns the incomplete gamma function P (a, x).
double gammp(double a, double x)
{
    void gcf(double *gammcf, double a, double x, double *gln);
    void gser(double *gamser, double a, double x, double *gln);
    void nrerror(char error_text[]);
    double gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
    if (x < (a+1.0)) {   // Use the series representation.
	gser(&gamser,a,x,&gln);
	return gamser;
    } else {   // Use the continued fraction representation
	gcf(&gammcf,a,x,&gln);
	return 1.0-gammcf;    //and take its complement.
    }
}
/*************************************************************/
//Returns the incomplete gamma function Q(a, x)  1 - P (a, x).
double gammq(double a, double x)
{
    void gcf(double *gammcf, double a, double x, double *gln);
    void gser(double *gamser, double a, double x, double *gln);
    void nrerror(char error_text[]);
    double gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
    if (x < (a+1.0)) {   // Use the series representation
	gser(&gamser,a,x,&gln);
	return 1.0-gamser;   // and take its complement.
    } else {    //Use the continued fraction representation.
	gcf(&gammcf,a,x,&gln);
	return gammcf;
    }
}
/*************************************************************/
//Returns the incomplete gamma function P (a, x) evaluated by its
//series representation as gamser.
//Also returns ln (a) as gln.
#define ITMAX 1000   //Maximum allowed number of iterations.
#define EPS 3.0e-7   //Relative accuracy.
void gser(double *gamser, double a, double x, double *gln)
{
    double gammln(double xx);
    int n;
    double sum,del,ap;
    *gln=gammln(a);
    if (x <= 0.0) {
	if (x < 0.0) nrerror("x less than 0 in routine gser");
	*gamser=0.0;
	return;
    } else {
	ap=a;
	del=sum=1.0/a;
	for (n=1;n<=ITMAX;n++) {
	    ++ap;
	    del *= x/ap;
	    sum += del;
	    if (fabs(del) < fabs(sum)*EPS) {
		*gamser=sum*exp(-x+a*log(x)-(*gln));
		return;
	    }
	}
	nrerror("a too large, ITMAX too small in routine gser");
	return;
    }
}
/*************************************************************/
//Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction 
//representation as gammcf. Also returns ln (a) as gln.
#define FPMIN 1.0e-30   // Number near the smallest representable doubleing-point number.
void gcf(double *gammcf, double a, double x, double *gln)
{
    double gammln(double xx);
    void nrerror(char error_text[]);
    int i;
    double an,b,c,d,del,h;
    *gln=gammln(a);
    b=x+1.0-a;    //Set up for evaluating continued fraction
    c=1.0/FPMIN;   // by modified Lentz's method (§5.2)
    d=1.0/b;    //with b0 = 0.
    h=d;
    for (i=1;i<=ITMAX;i++) {    //Iterate to convergence.
	an = -i*(i-a);
	b += 2.0;
	d=an*d+b;
	if (fabs(d) < FPMIN) d=FPMIN;
	c=b+an/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	del=d*c;
	h *= del;
	if (fabs(del-1.0) < EPS) break;
    }
    if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
    *gammcf=exp(-x+a*log(x)-(*gln))*h;   // Put factors in front.
}
/*************************************************************/
//Returns the value ln[(xx)] for xx > 0.
//Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
//accuracy is good enough.
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
/*************************************************************/
//Returns the complementary error function erfc(x) with fractional error everywhere less than
//1.2 × 10-7.
double erfcc(double x)
{
    double t,z,ans;
    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		    t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
			    t*(-0.82215223+t*0.17087277)))))))));
    // return (x >= 0.0) ? ans : 2.0-ans;
    return ans;
}
//Returns the complementary error function erfc(x).
double erffc(double x)
{
    double gammp(double a, double x);
    double gammq(double a, double x);
    return x < 0.0 ? 1.0+gammp(0.5,x*x) : gammq(0.5,x*x);
}

/***************************************/
 void nrerror(char error_text[]){
     fprintf ( stderr, "%s\n", error_text);
     exit(1);
}
