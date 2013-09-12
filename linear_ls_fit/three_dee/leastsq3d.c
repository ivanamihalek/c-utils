# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>


int main ( int argc, char * argv[]) {


    /* test case: vertices of a cube */
    double ** point;
    
    double center_of_mass[3] = {0.0};
    double p[3] = {0.0}; /* the direction vector */
    double I[3][3] = {{0.0}}; /* the "moments of inertia" */
    
    int no_points;
    int i, x, y;

    int  read_points (char * filename, double ***points_ptr, int * no_points_ptr);
    double **dmatrix (int rows, int columns);
    void normalize (double *p);

    if ( argc < 2 ) {
	fprintf (stdout, "Usage: %s <file name>.\n", argv[0]);
	exit (1);
    }
    

    /* read in the data*/
    if ( read_points(argv[1], &point, &no_points ) ) {
	fprintf (stderr, "error reding in data.\n");
	exit (1);
    }
    
    if (0) {
	//printf ( "number of points = %d\n", no_points);
	printf ( " x   y   z\n");
	for (i=0; i< no_points; i++ ) {
	    printf ( " %2d ", i+1);
	    for ( x=0; x<3; x++) {
		printf ( "   %5.2lf  ", point[i][x]);
	    }
	    printf ("\n");
	}
	exit (1);
    }
    
    /***************************/
    /* find the center of mass */
    /***************************/
    for ( x=0; x<3; x++)  center_of_mass[x] = 0.0;
    for (i=0; i< no_points; i++ ) {
	for ( x=0; x<3; x++) {
	    center_of_mass[x] += point[i][x];
	}
    }
    for ( x=0; x<3; x++) {
	center_of_mass[x] /= no_points;
    }
    printf ( "cm:  %7.4lf    %7.4lf    %7.4lf  \n",  
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
    char uplo = 'L'; /* matrix is stored as lower (fortran convention) */
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
        
    printf ( "     p:  %7.4lf    %7.4lf    %7.4lf  \n",  
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

/********************************************************************************/
/********************************************************************************/
#define BUFFLEN 150

int read_points ( char * filename, double ***points_ptr, int * no_points_ptr) {
    FILE * fptr = NULL;
    char line[BUFFLEN];
    double aux1, aux2, aux3;
    double ** points;
    int ptctr;
    
    /* open file */
    fptr = fopen ( filename, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", filename);
	return 1;
    }

    /* count points */
    memset (line,  0, BUFFLEN);
    ptctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( sscanf (line, "%lf %lf %lf ", &aux1, &aux2, &aux3) == 3  ) ptctr++;
    }

    if ( ! (points = dmatrix ( ptctr, 3) ) ){
	fprintf ( stderr, "Allocation error while reading.\n");
	return 1;
    }
    /* read in */
    rewind (fptr);
    ptctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( sscanf (line, "%lf %lf %lf ", &aux1, &aux2, &aux3) == 3  ) {
	    points[ptctr][0] = aux1;
	    points[ptctr][1] = aux2;
	    points[ptctr][2] = aux3;
	    ptctr++;
	}
    }
   
    *no_points_ptr = ptctr;
    *points_ptr  = points;

    fclose (fptr);
    
    return 0;
}
