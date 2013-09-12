# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>

# define NO_POINTS 20

int main ( int argc, char * argv[]) {


    /* test case: vertices of a cube */
    double ** point;
    
    double center_of_mass[3] = {0.0};
    double p[3] = {0.0}; /* the direction vector */
    double I[3][3] = {{0.0}}; /* the "moments of inertia" */
    
    int no_points = NO_POINTS;
    int i, x, y;

    double **dmatrix(int rows, int columns);
    double normalize (double *p);
    /* count the data */
    /* allocate */
    if ( ! (point = dmatrix  ( no_points, 3 ) ) ) {
	fprintf ( stderr, "Alloc error.\n");
	exit (1);
    }
    /* read in the data*/
    double a[3] = {5.1, 6.6, 2.0}, b[3] = {1.0, -1.0, 3.14};
    double t;

    srand48(time(NULL));
    
    t = 0.0;
    for (i=0; i< no_points; i++ ) {
	t += 0.5;
	printf ( " %2d ", i+1);
	for ( x=0; x<3; x++) {
	    point[i][x] = a[x]*t + b[x] + 4*(1-2*drand48());
	    printf ( "   %5.2lf  ", point[i][x]);
	}
	printf ("\n");
    }

    /***************************/
    /* find the center of mass */
    /***************************/
    for (i=0; i< no_points; i++ ) {
	for ( x=0; x<3; x++) {
	    center_of_mass[x] += point[i][x];
	}
    }
    for ( x=0; x<3; x++) {
	center_of_mass[x] /= no_points;
    }
    printf ( "cm:  %5.2lf    %5.2lf    %5.2lf  \n",  
	     center_of_mass[0], center_of_mass[1], center_of_mass[2]);
    
    /***********************************/
    /* move the points to the cm frame */
    /***********************************/
    for (i=0; i< no_points; i++ ) {
	for ( x=0; x<3; x++) {
	    point[i][x] -= center_of_mass[x];
	}
    }

    /**********************************/
    /* find the "moments of inertia"  */
    /**********************************/
    for (i=0; i< no_points; i++ ) {
	for ( x=0; x<3; x++) {  /* modulo = circular permutation */
	    I[x][x] += point[i][(x+1)%3]*point[i][(x+1)%3] +
		point[i][(x+2)%3]*point[i][(x+2)%3];
	    for ( y=x+1; y<3; y++) { /* off diag elements */
		I[x][y] -= point[i][x]*point[i][y];
	    }
	}
    }
    for ( x=0; x<3; x++) { 
	for ( y=x+1; y<3; y++) {
	    I[y][x] =  I[x][y];
	}
    }

    /*****************************************/
    /* diagonalize I[][], pick the direction
       with the smallest moement of inertia,
       and rotate back to the initial frame */
    /*****************************************/
    void dsyev_ ( char * jobz, char * uplo, int* N, double * A, int * leading_dim,
		  double * eigenvalues, double *workspace, int *workspace_size, int * retval);
    char jobz = 'V'; /* find evalues and evectors */
    char uplo = 'L'; /* amtrix is stored as lower (fortran convention) */
    int  N = 3; /* the order of matrix */
    int leading_dim = N;
    int retval;
    double A[N*N];
    double eigenvalues[N];
    double workspace[3*N];
    int workspace_size = 3*N;

    for ( x=0; x<3; x++) {
	for ( y=0; y<3; y++) {
	    A[x*3+y] = I[x][y];
	}
    }
   
    dsyev_ ( &jobz, &uplo, &N, A,  &leading_dim, eigenvalues,
	     workspace, &workspace_size, &retval);

    if ( retval ) {
	fprintf ( stderr, "Dsyev  error: %d.\n", retval);
	exit (1);
    }

    /* the eigenvalues are returned in ascending order, so the first guy i smine: */
    x = 0;
    for ( y=0; y<3; y++) {
	p[y] = A[x*3+y];/*this is  p, the direction vector   */
    }
        
    printf ( "     p:  %5.2lf    %5.2lf    %5.2lf  \n",  
	     p[0], p[1], p[2]);

    /* check: */
    for ( x=0; x<3; x++) { 
	p[x] = a[x];
    }
    normalize (p);
    printf ( " check:  %5.2lf    %5.2lf    %5.2lf  \n",  
	     p[0], p[1], p[2]);

    
   return 0;
}
/********************************************************************************/
/********************************************************************************/

void normalize (double *p) {
    double norm = 0.0;
    int x;
    for ( x=0; x<3; x++) {
	norm += p[x]*p[x];
    }
    norm = sqrt (norm);
    
    
    for ( x=0; x<3; x++) {
	p[x] /= norm;
    }
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



FILE * efopen(char * name, char * mode)
{

    FILE * fp;


    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,  
	      "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }

    return fp;

}
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



/* free a  matrix  */
void free_matrix(void **m)
{
    free(m[0]);
    free(m);
}


