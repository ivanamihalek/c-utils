/* readin in 2 pdb, presumably roughly aligned structurally, */
/* and a list of residues in the first structure on which the*/
/* improved alignemnt should focus */
/* output the tranformation matrix for the second structure */
/* for the improved alignment */
/* Ivana, Oct, 2010 */
# include "smallig.h"

# define MAX_DIST_TO_CONSIDER 4.0

/* the following functions are defined below main()*/

int smith_waterman (int max_i, int max_j, double **similarity,
		    int *map_i2j, int * map_j2i, double * aln_score);
double  two_point_distance (double point1[3], double point2[3] );

/*************************************************************************/
/*************************************************************************/
int main (int argc, char * argv[]) {


    Protein protein1, protein2;
    Atom* atom[2];
    
    double cm[2][3]= {{0.0}};
    double I[3][3] = {{0.0}};
    double principal_axis[2][3][3], principal_moment[2][3];
    double **x, **y, **R; /* aux vectors fro finding rotations */
    double ** rot_coord;
    double q[4], T[3], rmsd;
    
    char pdbname1 [BUFFLEN] = "\0";
    char pdbname2 [BUFFLEN] = "\0";
    int no_atoms[2] = {0};
    int i,j,a,l;

    double diagonalize (double I[3][3], double evec[3][3], double evals[3]);
    int mappedCoords2rotation (double **x, double **y, int no_vectors,
			       double q[4], double **R, double T[3], double *rmsd);  
    if ( argc < 2 ) {
	printf ( "Usage: %s <pdbname1> <pdbname2>\n", argv[0] );
	exit (1);
    } 
    sprintf ( pdbname1, "%s",  argv[1]);
    sprintf ( pdbname2, "%s",  argv[2]);

    /*********************************************/
    /* read in the two structures                */
    if ( read_pdb ( pdbname1, &protein1, '\0'))  exit (1);
    if ( read_pdb ( pdbname2, &protein2, '\0'))  exit (1);
    if (0) {
	printf ("read in %20s,  %s, length %d, atoms %d\n",
		pdbname1, protein1.sequence[0].res_type, protein1.length, protein1.sequence[0].no_atoms);
	printf ("read in %20s,  %s, length %d, atoms %d\n",
		pdbname2, protein2.sequence[0].res_type, protein2.length, protein2.sequence[0].no_atoms);
    }
    /*********************************************/
    /* calculate moments of inertia for both     */
    atom[0]     = protein1.sequence[0].atom;
    no_atoms[0] = protein1.sequence[0].no_atoms;
    atom[1]     = protein2.sequence[0].atom;
    no_atoms[1] = protein2.sequence[0].no_atoms;
     /* sanity: */
    for (l=0; l<2; l++) {
	if  (! no_atoms[l] ) {
	    fprintf (stderr, "no atoms in %s(?)\n", pdbname1);
	    exit (1);
	}
    }
   
    for (l=0; l<2; l++) {
	/*cm*/
	for (a=0; a<no_atoms[l]; a++) {
	    for (i=0; i<3; i++)  cm[l][i] += atom[l][a].coord[i];
	}
	 for (i=0; i<3; i++)  cm[l][i] /= no_atoms[l];
	/* points in cm */
	for (a=0; a<no_atoms[l]; a++) {
	    for (i=0; i<3; i++)  atom[l][a].coord[i] -= cm[l][i];
	}
	
	/* moment of inertia tensor */
	memset (I[0], 0.0, 3*3*sizeof(double));
	for (a=0; a<no_atoms[l]; a++) {
	    for (i=0; i<3; i++) {
		double d1 =  atom[l][a].coord[(i+1)%3];
		double d2 =  atom[l][a].coord[(i+2)%3];
		I[i][i] += d1*d1 + d2*d2;
		for (j=0; j<3; j++) {
		    I[i][j] -=  atom[l][a].coord[i]*atom[l][a].coord[j];
		}
	    }
	}
	diagonalize (I, principal_axis[l], principal_moment[l]);
    }


    if ( ! (x=dmatrix(3,no_atoms[0]) ) ) exit (1);
    if ( ! (y=dmatrix(3,no_atoms[1]) ) ) exit (1);

    for (i=0; i<3; i++) {
	for (j=0; j<3; j++) {
	    x[j][i] = principal_axis[0][i][j];
	    y[j][i] = principal_axis[1][i][j];
	}
    }

    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */

    mappedCoords2rotation (x, y, 3, q, R, T, &rmsd);
    
    /*******************************************/
    /* output the tfm matrix                   */
    if ( 0 ) {
	for (i=0; i<3; i++) {
	    for (j=0; j<3; j++) {
		printf (" %8.4lf", R[i][j]);
	    }
	    printf (" %8.4lf\n", cm[0][i]);
	}
    }

    /* rotate (the first structure), move back to the old cm,
       and output */
    if ( ! (rot_coord=dmatrix(no_atoms[0],3) ) ) exit (1);
    for (a=0; a<no_atoms[0]; a++) {
       
	for (i=0; i<3; i++) {
	    rot_coord[a][i] = 0.0;
	    for (j=0; j<3; j++) {
		rot_coord[a][i] += R[i][j]*atom[0][a].coord[j] ;
	    }
	}

	/* actually, moe it to the cm of the *other* molecule */
	for (i=0; i<3; i++)  rot_coord[a][i]  += cm[1][i];
	/*       123456789012345678901234567890*/
	printf ("HETATM%5d  %2s  %3s %1c   %4d %8.3f%8.3f%8.3f \n",
		a, atom[0][a].type, protein1.sequence[0].res_type, 'L', 1,
		rot_coord[a][0], rot_coord[a][1], rot_coord[a][2]);
	
    }
    return 0;
}

/*************************************************************************/
/*************************************************************************/
double diagonalize (double I[3][3], double evec[3][3], double evals[3]) {
    
    char jobz = 'V'; /* find evalues and evectors */
    char uplo = 'L'; /* amtrix is stored as lower (fortran convention) */
    int  N = 3; /* the order of matrix */
    int leading_dim = N;
    int i, j, retval;
    double A[N*N];
    double workspace[3*N];
    int workspace_size = 3*N;
    void dsyev_(char * jobz, char * uplo, int* N,
		       double * A, int * leading_dim,
		       double * eigenvalues,
		       double *workspace, int *workspace_size,
		       int * retval);
    for (i=0; i < 3; i++) {
	for (j=0; j < 3; j++) {
	    A[i*3+j] = I[i][j];
	}
    }
   
    dsyev_ ( &jobz, &uplo, &N, A,  &leading_dim, evals,
	     workspace, &workspace_size, &retval);
    if ( retval ) {
	fprintf (stderr, "error in dsyev()\n");
	exit (1);
    }

   for (i=0; i < 3; i++) {
       double norm = 0;
	for (j=0; j < 3; j++) {
	    evec[i][j] = A[i*3+j];
	    norm += evec[i][j]*evec[i][j];
	}
    }
 
    
    return 0;
}

/*********************************************************/
/*********************************************************/

int mappedCoords2rotation (double **x, double **y, int no_vectors,
			   double q[4], double **R, double T[3], double *rmsd) {

    double x_mp[3], y_mp[3];
    int  ctr;
     int  i, j;
 
    double ATA     [4][4] = {{0.0}};
    double prev_ATA[4][4] = {{0.0}};
    double ATA_sum [4][4] = {{0.0}};
    double a[3] = {0.0}, b[3] = {0.0};

    
    int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		       double result[4][4]);
    int construct_ATA (double ATA[4][4], double a[3], double  b[3]);

    /* note how we pass the matrix: pointer to the first element in the block */
    void dsyev_ (char * jobz, char *uplo,  int *n,
		  double *A, int * lda, double * w, double * work, int * lwork, int *info);

    if (!no_vectors) {
	*rmsd = -1;
	return 1;
    }

    memset ( &(q[0]), 0, 4*sizeof(double) );
    /* turn the matching atoms into vectors x and y - use only c-alphas*/
     
    /* check: */
    if (0) {
	printf (" Number of vectors read in: %d. \n", no_vectors);
	for ( ctr =0; ctr < no_vectors; ctr++ ) {
	    printf ("\t x%1d   %10.4lf  %10.4lf  %10.4lf   ",
		    ctr, x[0][ctr], x[1][ctr], x[2][ctr]);
	    printf ("\t y%1d   %10.4lf  %10.4lf  %10.4lf \n",
		    ctr, y[0][ctr], y[1][ctr], y[2][ctr]);
	}
	exit (1);
    }

    /* find the meanpoints: */
    for ( i =0; i < 3; i++ ) {
	x_mp[i] = 0.0;
	y_mp[i] = 0.0;
    }
    for ( ctr =0; ctr < no_vectors; ctr++ ) {
	for ( i =0; i < 3; i++ ) {
	    x_mp[i] += x[i][ctr];
	    y_mp[i] += y[i][ctr];
	}
    }
    for ( i =0; i < 3; i++ ) {
	x_mp[i] /= no_vectors;
	y_mp[i] /= no_vectors;
    }
    /* subtract them from x, y */
    for ( ctr =0; ctr < no_vectors; ctr++ ) {
	for ( i =0; i < 3; i++ ) {
	    x[i][ctr] -= x_mp[i];
	    y[i][ctr] -= y_mp[i];
	}
    }
    /* B = ATA_sum matrix to diagonalize in order to get the quaternion */
    for ( ctr =0; ctr < no_vectors; ctr++ ) {
   	for (i=0; i<3; i++ ) {
	    a[i] = y[i][ctr] + x[i][ctr];
	    b[i] = y[i][ctr] - x[i][ctr];
	}
 	construct_ATA (ATA, a, b);
	add_matrices (prev_ATA, ATA, ATA_sum);
	memcpy (prev_ATA[0], ATA_sum[0], 4*4*sizeof(double));
    }
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA_sum[i][j] /= no_vectors;
	}
    }
    /* diagonalize ATA_sum - the eigenvector corresponsing to the
       smallest lambda is the quaternion we are looking for; the
       eigenvalue is the rmsd*/
    /* use the nomenclature from dsyev*/
    char jobz= 'V'; /*Compute eigenvalues and eigenvectors.*/
    char uplo= 'U'; /* Upper triangle of A (the matrix we are diagonalizing) is stored; */
    int  n = 4;     /* order and the leading dimension of A */
    int  lda = 4;
    double ** A;
    int  info;
    int  lwork = 200;
    double w [4];
    double work[200];
    
    if ( !( A=dmatrix(4,4) ) ) exit (1);
    memcpy (A[0], ATA_sum[0], 4*4*sizeof(double));


   /* note how we pass the matrix: */
    dsyev_ ( &jobz, &uplo,  &n, A[0], &lda, w, work, &lwork, &info);
    if (  ! info) {
	*rmsd = sqrt (w[0]);
	for (i=0; i<4; i++ ) q[i] = A[0][i];
	if (0) {
	    /* w contains the eigenvalues */
	    printf ("\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", w[i]);
	    printf ("\nrmsd: %8.3lf \n", *rmsd);
	    printf ("quat:\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", q[i]);
	    printf ("\n");
	    /* printf (" opt lwork: %d\n", (int) work[0]); */
	}
    } else {
	fprintf (stderr, "Error in dsyev().\n");
	exit (1);
    }
    
    /* construct the rotation matrix R */
    quat_to_R (q,R);
    /* T = y_mp - R x_mp */
    for (i=0; i<3; i++ ) {
	T[i] = y_mp[i];
	for (j=0; j<3; j++ ) {
	    T[i] -= R[i][j]*x_mp[j];
	}
    }
   
    
    free_dmatrix(A);

 
    return 0;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double braket  (double ATA[4][4],double  q[4]){
    int i,j;
    double value = 0.0;

    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    value += q[i]*ATA[i][j]*q[j];
	}
    }
   
    return value;
}
 
/*******************************************************************************/
int construct_ATA (double ATA[4][4], double a[3], double  b[3]){

    int i,j,k;
    double A[4][4] = {{ 0.0, -b[0], -b[1], -b[2]},
		      {b[0],   0.0, -a[2],  a[1]},
		      {b[1],  a[2],   0.0, -a[0]},
		      {b[2], -a[1],  a[0],   0.0}};
    
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA[i][j] = 0.0;
	    for (k=0; k<4; k++ ) {
		ATA[i][j] += A[k][i]*A[k][j];
	    }
	}
    }
     
    return 0;
}

/*******************************************************************************/
int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		   double result[4][4]){
    int i,j;

    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    result[i][j] = matrix1[i][j] + matrix2[i][j];
	}
    }
    return 0;

}

    
/*************************************************************************/
/*************************************************************************/
double  two_point_distance (double point1[3], double point2[3] ) {
    int i;
    double aux;
    double d = 0;
    for (i=0; i<3; i++) {
	aux = point1[i] - point2[i];
	d  += aux*aux;
    }
    return  sqrt (d);
}
/************************************************************/
/************************************************************/
int quat_to_R ( double *quat, double **R ){

    
    double q0 = quat[0];
    double q1 = quat[1];
    double q2 = quat[2];
    double q3 = quat[3];


    R[0][0] = 1 -2*q2*q2 -2*q3*q3;
    R[1][1] = 1 -2*q3*q3 -2*q1*q1 ;
    R[2][2] = 1 -2*q1*q1 -2*q2*q2 ;

    R[0][1] = 2*q1*q2 - 2*q0*q3;
    R[1][0] = 2*q1*q2 + 2*q0*q3;
    
    R[0][2] = 2*q1*q3 + 2*q0*q2;
    R[2][0] = 2*q1*q3 - 2*q0*q2;
    
    R[1][2] = 2*q2*q3 - 2*q0*q1;
    R[2][1] = 2*q2*q3 + 2*q0*q1;

    return 0;

}
/************************************************************/
/************************************************************/
