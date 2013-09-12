# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>


double distance ( double * A, double *B);

double **dmatrix(int rows, int columns);
void * emalloc(int  size);


/************************************************************/
int main (int argc, char * argv[]) {

    double *q;
    double aux, dist_old, dist_new, dist_best;
    double norm,  d_phi = M_PI/180, d_psi;
    double R[3][3], v[3];
    double ** backbone, ** current_bb, ** trial_bb;
    double *target, *translation, angle;
    int bb_length, bb_atom_ctr;
    int done = 0, step = 0, max_step = 500;
    int sign;
    int accepted;
    int i;
    
    /***************************/
    int construct_q (double * init_pt, double * final_pt, double angle, double *q);
    int  rotate_the_rest (double ** current_bb, int bb_length, double * pivot,
			  double R[3][3], double ** trial_bb);
    int quat_to_R  (double quat[4], double R[3][3]);
    int read_pdb   (char * pdbname, double *** backbone_ptr, int * bb_length_ptr);
    int output_pdb (char * pdbname, double ** backbone_ptr);  
    /***************************/
    if ( argc < 5 ) {
	fprintf (stderr, "Usage: %s <pdb_file>  <target coords>.\n", argv[0]);
	exit (1);
    }
    if ( ! (target = emalloc (3*sizeof(double)) ) )  exit (1);
    for (i=0; i<3; i++) {
	target[i] = atof (argv[2+i]);
    }

    
    /****************************/
    /* read in the backbone     */
    /****************************/
    read_pdb ( argv[1],  &backbone, &bb_length);
    /* allocate space for various
       versions of bb we'll keep */
    if ( ! (current_bb = dmatrix (bb_length, 3)))  exit (1);
    if ( ! (  trial_bb = dmatrix (bb_length, 3)))  exit (1);
   
  
    /****************************/
    /* the original dist        */
    /****************************/
    dist_old = distance (backbone[bb_length-1], target);
    //printf ( "original distance btw last atom and the target:   %8.2lf \n",  sqrt(dist_old)); 
  
    /****************************/
    /* init guess               */
    /****************************/
    if ( ! (q = emalloc (4*sizeof(double) )))  exit (1);
    norm  = 0;
    for (i=0; i<3; i++) {
	v[i]  = backbone[bb_length-1][i] + target[i] - 2*backbone[0][i];
	norm += v[i]*v[i];
    }
    norm = sqrt (norm);
    for (i=0; i<3; i++) v[i] /= norm;
    
    q[0] = 0.0;
    for (i=0; i<3; i++) q[i+1] = v[i];
    quat_to_R (q, R);
    
    /* rotate the whole 9 yds */
    memcpy ( current_bb[0], backbone[0], 3*sizeof(double) );
    rotate_the_rest ( backbone+1, bb_length-1, backbone[0], R, current_bb+1);
 
    dist_best = dist_old = distance( current_bb[bb_length-1], target);
 
    //printf ( "first guess distance btw last atom and the target:  %8.2lf \n\n",  dist_best);
    done = ( dist_old < 1.e-1 );

    //exit (0);
     
    /****************************/
    /* MC  loop                 */
    /****************************/
    /* ! NOTE: assumed all N, CA, C present, starting with N */
    srand48(356783);
    if ( ! (translation = emalloc (3*sizeof(double) ) ) )  exit (1);
    //current_bb = bb;
    for (step=1; step<=max_step && !done; step++) {
	for (bb_atom_ctr =1; bb_atom_ctr< bb_length-1; bb_atom_ctr++ ) {
	    /* bb_atom_ctr atom is the pivot (origin) */
	    /* the direction is decided by the previous atom */
	    /* if C is the pivot, that would mean changing the
	       peptide bond angle, so we  skip that */
	    
	    switch ( bb_atom_ctr%3 ) {
	    case 0: 
		continue; /* this is peptide bond */
		break;
	    case 1:
		angle = d_phi;
		break;
	    case 2: 
		angle = d_psi;
		break;
	    }
	    sign = (drand48() < 0.5) ? 1 : -1;

	    construct_q (current_bb[bb_atom_ctr-1], current_bb[bb_atom_ctr],
			 sign*angle, q);
	    quat_to_R (q, R);

	    memcpy (trial_bb[0], current_bb[0], (bb_atom_ctr+1)*3*sizeof(double));
	    
	    rotate_the_rest (current_bb+bb_atom_ctr+1, bb_length-bb_atom_ctr-1,
			     current_bb[bb_atom_ctr], R, trial_bb+bb_atom_ctr+1);
	    
	    dist_new = distance (trial_bb[bb_length-1], target);
	    
	    /* MC acceptability */
	    accepted = ( dist_new < dist_old);
	    /* if accepted, current_bb = new_bb */
	    if ( accepted ) {
		memcpy (current_bb[0], trial_bb[0], bb_length*3*sizeof(double) );
		dist_old = dist_new;
		//printf ( " %4d %4d  %8.2lf  accepted\n", step,bb_atom_ctr, dist_new);
		done = ( dist_new < 1.e-1 );
	    } else {
		//printf ( " %4d %4d  %8.2lf not accepted\n", step,bb_atom_ctr, dist_new);
	    }
	    
	}
	
    }

    
    /* output bb (to be postprocessed by scwrl) */
    output_pdb ( argv[1],  current_bb);
    
    return 0;
}


/**********************************************************/
/**********************************************************/
/**********************************************************/
int  rotate_the_rest (double ** current_bb, int bb_length, double * pivot,
		      double R[3][3], double ** trial_bb ) {
    
    int bb_atom_ctr, i, j;
    
    
    for (bb_atom_ctr =0; bb_atom_ctr < bb_length; bb_atom_ctr++ ) {
	for (i=0; i<3; i++) {
	    trial_bb[bb_atom_ctr][i] = 0;
	    for (j=0; j<3; j++) {
		trial_bb[bb_atom_ctr][i] +=
		    R[i][j]*(current_bb[bb_atom_ctr][j]-pivot[j]);
	    }
	    trial_bb[bb_atom_ctr][i] += pivot[i];
	}
    }

    return 0;
    
}
		       
/**************************/
int construct_q (double * init_pt, double * final_pt, double angle, double *q) {

    double v[3];
    double norm, cosn, sine;
    int i;

    cosn = cos(angle/2);
    sine = sin(angle/2);
    
    norm  = 0;
    for (i=0; i<3; i++) {
	v[i] = final_pt[i] - init_pt[i];
	norm += v[i]*v[i];
    }
    norm = sqrt (norm);
    for (i=0; i<3; i++) v[i] /= norm;
    
    q[0] = cosn;
    for (i=0; i<3; i++) q[i+1] = sine*v[i];

    return 1;
}


/**************************/
int move_point (double * pivot, double *endpoint,
		double d_theta, double d_phi, double *translation ){

    double theta, phi, norm;
    double v[3], new_endpoint[3];
    int i;

    norm = 0;
    for (i=0; i<3; i++){
	v[i] = endpoint[i] - pivot[i];
	norm +=  v[i]*v[i];
    }
    norm = sqrt (norm);

    if ( norm < 0.1 ) return 0; /* pivot & endpoint about the same*/
   
    for (i=0; i<3; i++) v[i] /= norm;

    theta = acos (v[2]);

    if (theta < 1.e-6) {
	theta = 0.0;
	phi   = 0.0;
    } else {
	phi = acos (v[0]/sin(theta) );
	if ( v[1] < 0 )  phi = 2*M_PI - phi;
    }

    theta += d_theta;
    phi   += d_phi;

    v[2] = cos (theta);
    v[0] = sin (theta)*cos(phi);
    v[1] = sin (theta)*sin(phi);

    for (i=0; i<3; i++) {
	v[i] *= norm;
	new_endpoint[i] = pivot[i] + v[i];
	translation[i]  = new_endpoint[i] - endpoint[i];
    }

# if 0
    for (i=0; i<3; i++) {
	printf (" %8.3lf ", translation[i]);
    }
    printf ("\n");
    exit (1);
# endif
    
    return 0;
    
}



/**********************************************************/
/**********************************************************/
/**********************************************************/
# define  BUFFLEN 150
# define         PDB_ATOM_ATOM_NAME         12
# define         PDB_ATOM_ATOM_NAME_LEN     4     
# define        PDB_ATOM_X                  30
# define        PDB_ATOM_X_LEN               8
# define        PDB_ATOM_Y                  38
# define        PDB_ATOM_Y_LEN               8
# define        PDB_ATOM_Z                  46
# define        PDB_ATOM_Z_LEN               8

int output_pdb ( char * pdbname,  double ** backbone) {
    
    FILE    * fptr = NULL;
    char line[BUFFLEN];
    char tmp[BUFFLEN], *auxptr;
    int bb_atom_ctr;
    int ctr, nonblank;
    int isNC, isCA;
    
    /* open file */
    fptr = fopen ( pdbname, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", pdbname);
	return 1;
    }

    /******************************************************/
    /* rewind, and read in the atom coordinates           */
    memset (line,  0, BUFFLEN);
    bb_atom_ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	
	if(  strncmp(line,"ATOM", 4)) continue;
	
	/* read in atom info */
	auxptr = line+ PDB_ATOM_ATOM_NAME;
	memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	/* skip initial blanks*/
	ctr  = 0;
	while ( !(isalpha (*(auxptr + ctr))) &&
		(ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	/* copy alphanum info */
	nonblank = 0;
	while (  isalpha (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
	    tmp[nonblank] =  *(auxptr +ctr);
	    nonblank ++;
	    ctr++;
	}
       
	/* is this a backbone atom?*/
	isNC = (nonblank == 1) && !( strcmp ( tmp, "N") && strcmp ( tmp, "C"));
	isCA = (nonblank == 2) && ! strcmp ( tmp, "CA");
	
	if ( isNC || isCA ) {

	    memset (tmp, 0, PDB_ATOM_X_LEN);
	    sprintf ( tmp, "%8.3lf", backbone[bb_atom_ctr][0]);
	    strncpy ( line+PDB_ATOM_X, tmp, PDB_ATOM_X_LEN);

	    memset (tmp, 0, PDB_ATOM_Y_LEN);
	    sprintf ( tmp, "%8.3lf", backbone[bb_atom_ctr][1]);
	    strncpy ( line+PDB_ATOM_Y, tmp, PDB_ATOM_Y_LEN);

	    memset (tmp, 0, PDB_ATOM_Z_LEN);
	    sprintf ( tmp, "%8.3lf", backbone[bb_atom_ctr][2]);
	    strncpy ( line+PDB_ATOM_Z, tmp, PDB_ATOM_Z_LEN);
	    
	    bb_atom_ctr ++;
	    printf ( "%s", line);
	}
    }

    fclose (fptr);

    return 0;
}



/**********************************************************/

int read_pdb ( char * pdbname,  double *** backbone_ptr, int * bb_length_ptr) {
    
    FILE    * fptr = NULL;
    char line[BUFFLEN];
    char tmp[BUFFLEN], *auxptr;
    double **backbone;
    int bb_atom_ctr;
    int bb_length;
    int ctr, nonblank;
    int isNC, isCA;
    
    /* open file */
    fptr = fopen ( pdbname, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", pdbname);
	return 1;
    }
    /******************************************************/
    /* count bb atoms                                     */
    memset (line,  0, BUFFLEN);
    bb_atom_ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	
	if(  strncmp(line,"ATOM", 4)) continue;
	
	/* read in atom info */
	auxptr = line+ PDB_ATOM_ATOM_NAME;
	memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	/* skip initial blanks*/
	ctr  = 0;
	while ( !(isalpha (*(auxptr + ctr))) &&
		(ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	/* copy alphanum info */
	nonblank = 0;
	while (  isalpha (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
	    tmp[nonblank] =  *(auxptr +ctr);
	    nonblank ++;
	    ctr++;
	}
	
	/* is this a backbone atom?*/
	isNC = (nonblank == 1) && !( strcmp ( tmp, "N") && strcmp ( tmp, "C"));
	isCA = (nonblank == 2) && ! strcmp ( tmp, "CA");
	
	if ( isNC || isCA ) {
	    bb_atom_ctr ++;
	}
    }
    
    bb_length = bb_atom_ctr;
    
    /******************************************************/
    /* allocate backbone space                            */
    if ( ! (backbone=dmatrix (bb_length,3) ) ) return 1;
    

    /******************************************************/
    /* rewind, and read in the atom coordinates           */
    rewind ( fptr);
    memset (line,  0, BUFFLEN);
    bb_atom_ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	
	if(  strncmp(line,"ATOM", 4)) continue;
	
	/* read in atom info */
	auxptr = line+ PDB_ATOM_ATOM_NAME;
	memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	/* skip initial blanks*/
	ctr  = 0;
	while ( !(isalpha (*(auxptr + ctr))) &&
		(ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	/* copy alphanum info */
	nonblank = 0;
	while (  isalpha (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
	    tmp[nonblank] =  *(auxptr +ctr);
	    nonblank ++;
	    ctr++;
	}
       
	/* is this a backbone atom?*/
	isNC = (nonblank == 1) && !( strcmp ( tmp, "N") && strcmp ( tmp, "C"));
	isCA = (nonblank == 2) && ! strcmp ( tmp, "CA");
	
	if ( isNC || isCA ) {
	    
	    strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	    tmp[PDB_ATOM_X_LEN] = '\0';
	    backbone[bb_atom_ctr][0] = atof(tmp);

	    strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	    tmp[PDB_ATOM_Y_LEN] = '\0';
	    backbone[bb_atom_ctr][1] = atof(tmp);

	    strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	    tmp[PDB_ATOM_Z_LEN] = '\0';
	    backbone[bb_atom_ctr][2] = atof(tmp);
	    
	    bb_atom_ctr ++;
	}
    }

    fclose (fptr);
    
    *backbone_ptr = backbone;
    *bb_length_ptr = bb_length;

    return 0;
}

/**********************************************************/
/**********************************************************/
/**********************************************************/

int quat_to_R ( double quat[4], double R[3][3] ){

    
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

/************************************/
double distance ( double * A, double *B) {
    double dist = 0.0, aux;
    int i;
    for (i=0; i<3; i++) {
	aux =  A[i] - B[i];
	dist += aux*aux;
    }
    return sqrt (dist);
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

/************************************/
void * emalloc(int  size) {
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}




/* junkyard */ 
# if 0
    for (bb_atom_ctr =0; bb_atom_ctr< bb_length; bb_atom_ctr++ ) {
	for (i=0; i<3; i++) {
	    printf (" %8.3lf ", backbone[bb_atom_ctr][i]);
	}
	printf ("\n");
    }
    exit (1);	]

    
    for (i=0; i<3; i++) {
	for (j=0; j<3; j++) {
	    printf ( "%8.3lf ", R[i][j]);
	}
	printf ( "%8.3lf \n", origin[i] );
    }
# endif
