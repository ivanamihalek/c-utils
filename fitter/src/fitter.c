/* readin in 2 pdb, presumably roughly aligned structurally, */
/* and a list of residues in the first structure on which the*/
/* improved alignemnt should focus */
/* output the tranformation matrix for the second structure */
/* for the improved alignment */
/* Ivana, Oct, 2010 */
# include "fitter.h"

# define MAX_DIST_TO_CONSIDER 4.0

/* the following functions are defined below main()*/

int smith_waterman (int max_i, int max_j, double **similarity,
		    int *map_i2j, int * map_j2i, double * aln_score);
double  two_point_distance (double point1[3], double point2[3] );

/*************************************************************************/
/*************************************************************************/
int main (int argc, char * argv[]) {


    Protein protein1, protein2;
    double ** similarity;
    double ca1[3], ca2[3];
    double d, d0 = 10.0, aln_score;
    double *q, *T, **R,rmsd;
    char pdbname1  [BUFFLEN] = "\0";
    char pdbname2  [BUFFLEN] = "\0";
    char list_fname [BUFFLEN] = "\0";
    char out_fname [BUFFLEN] = "\0";
    char ** list;
    int *residue_map_i2j, *residue_map_j2i;
    int resctr1, resctr2, list_length;
    int res2_on_the_list;
    int i;


    if ( argc < 5 ) {
	printf ( "Usage:  %s <pdbname1> <pdbname2> <list of residues (filename)> <name_root>.\n\n", argv[0] );
	printf ( "Output: %s will rotate <pdbname1> onto  <pdbname2> \n", argv[0] );
	printf ( "        and output the rotated structure to <name_root>.rot.pdb.\n");
	exit (1);
    } 
    sprintf ( pdbname1,   "%s",  argv[1]);
    sprintf ( pdbname2,   "%s",  argv[2]);
    sprintf ( list_fname, "%s", argv[3]);
    sprintf ( out_fname,  "%s.rot.pdb", argv[4]);

    /*********************************************/
    /* read in the two structures                */
    if ( read_pdb ( pdbname1, &protein1, '\0'))  exit (1);
    if ( read_pdb ( pdbname2, &protein2, '\0'))  exit (1);
    printf ("read in %s, length %d\n", pdbname1, protein1.length);
    printf ("read in %s, length %d\n", pdbname2, protein2.length);
    /* read in the list of residues to be matched in protein 2 */
    read_list ( list_fname, &list, &list_length);
    printf ("list length %d\n", list_length);

    /**********************************************/
    /* allocate space for the "similarity" matrix */
    similarity = dmatrix (protein1.length, protein2.length);
    if ( ! similarity) exit (1);
    /**********************************************/
    /* allocate space for the maps between the    */
    /* sets of Ca's                               */
    if ( ! (residue_map_i2j = emalloc (protein1.length*sizeof(int))) ) exit (1); 
    if ( ! (residue_map_j2i = emalloc (protein2.length*sizeof(int))) ) exit (1);

    /* sanity check: */
    res2_on_the_list = 0;
    for (resctr2=0; resctr2< protein2.length; resctr2++) {
	for ( i=0; i<list_length; i++) {
	    if ( ! strcmp (protein2.sequence[resctr2].pdb_id, list[i]) ) {
		res2_on_the_list = 1;
		break;
	    }
	}
    }
    if ( ! res2_on_the_list ) {
	fprintf (stderr, "No residues from %s found in the list %s.\n",
		 pdbname2, list_fname);
	exit (1);
    }

    
    for (resctr1=0; resctr1<protein1.length; resctr1++) {
	for (resctr2=0; resctr2<protein2.length; resctr2++) {
	    similarity[resctr1][resctr2] = -1;
	}
    }
    /* evaluate the similarity matrix from Ca distances */
    for (resctr2=0; resctr2< protein2.length; resctr2++) {
	
	res2_on_the_list = 0;
	find_Calpha (&protein2, resctr2,  ca2);
	
	for ( i=0; i<list_length; i++) {
	    if ( ! strcmp (protein2.sequence[resctr2].pdb_id, list[i]) ) {
		res2_on_the_list = 1;
		break;
	    }
	}
	if (! res2_on_the_list) continue;
	
	for (resctr1=0; resctr1 < protein1.length; resctr1++) {
	    find_Calpha (&protein1, resctr1,  ca1 );
	    d = two_point_distance (ca1, ca2);

	    
	    if ( d<MAX_DIST_TO_CONSIDER ) { 
		similarity[resctr1][resctr2] = exp (-d/d0);
	    }
	}
    }
    /* find the maps using SW                            */
    smith_waterman (protein1.length, protein2.length, similarity,
		    residue_map_i2j, residue_map_j2i, &aln_score);
    printf ("alignment score: %8.3lf\n", aln_score);

    /* find the quaternion for the tfm  minimizing the rmsd */
    if ( ! (q=emalloc (4*sizeof(double))) ) return 1;
    if ( ! (T=emalloc (3*sizeof(double))) ) return 1;
    if ( ! (R=dmatrix (3,3) ) ) return 1;
    map2rotation (&protein1, &protein2, residue_map_i2j, q, T, &rmsd);
    printf ("rmsd: %8.4lf\n\n", rmsd);
    

    /* translate quaternion into tfm matrix */
    quat_to_R (q, R);

    /*******************************************/
    /* output (the tfm matrix  will bw written on top of the               */
    pdb_output (pdbname1, pdbname2, out_fname, R, T, protein1.sequence, protein1.length);


    free (T);
    free (q);
    free_dmatrix(R);
    return 0;
}


/*************************************************************************/
/*************************************************************************/
int smith_waterman (int max_i, int max_j, double **similarity,
		      int *map_i2j, int * map_j2i, double * aln_score) {

    double **F; /*alignment_scoring table*/
    char   **direction;
    double gap_opening   = -0.2;
    double gap_extension = -0.1;
    double endgap        =  0.0;
    double far_far_away  = -100.0;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    double F_max;
    int use_endgap = 1;
    int F_max_i, F_max_j;
    int i,j;

     /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;


    
    /* fill the table */
    F_max = far_far_away;
    F_max_i = 0;
    F_max_j = 0;
   
    for (i=0; i<= max_i; i++) {
	for (j=0; j<=max_j; j++) {

	    if ( !i && !j ) {
		F[0][0] = 0;
		direction[i][j] = 'd';
		continue;
	    }
	    
	    if ( i && j ){
		
		if ( direction[i-1][j] == 'i' ) {
		    /*  gap extension  */
		    penalty = (use_endgap&&j==max_j) ? endgap : gap_extension;		    
		} else {
		    /*  gap opening  */
		    penalty = (use_endgap&&j==max_j) ? endgap : gap_opening;
		}
		i_sim =  F[i-1][j] + penalty;
		
		if ( direction[i][j-1] =='j' ) {
		    penalty = (use_endgap&&i==max_i) ? endgap : gap_extension;		    
		} else {
		    penalty = (use_endgap&&i==max_i) ? endgap : gap_opening;		    
		}
		j_sim = F[i][j-1] +  penalty;
       	
		
		diag_sim =  F[i-1][j-1] + similarity [i-1][j-1] ;
		
		
	    } else if ( j ) {
		
		if ( use_endgap) {
		    penalty = endgap;
		} else {
		    if ( direction[i][j-1] =='j' ) {
			penalty = gap_extension;
		    } else {
			penalty = gap_opening;
		    }
		}
		j_sim = F[i][j-1] + penalty;
		
		i_sim = diag_sim = far_far_away;

		
	    } else if ( i ) {
		
		if ( use_endgap) {
		    penalty = endgap;
		} else {
		    if ( direction[i-1][j] == 'i' ) {
			penalty =  gap_extension;
		    } else {
		        penalty =  gap_opening;
		    }
		}
		i_sim = F[i-1][j] + penalty;
		
		j_sim = diag_sim = far_far_away;
		

	    } 

	    max_sim = diag_sim;
	    direction[i][j] = 'd';
	    if ( i_sim > max_sim ){
		max_sim = i_sim;
		direction[i][j] = 'i';
	    }
	    if ( j_sim > max_sim ) {
		max_sim = j_sim;
		direction[i][j] = 'j';
	    }

	    if ( max_sim < 0.0 ) max_sim = 0.0;
	    
	    F[i][j] = max_sim;
	    if ( F_max < max_sim ) {
		/* TODO: tie break here */
		F_max = max_sim;
		F_max_i = i;
		F_max_j = j;
		
	    }
	    
	}
    }
    
    /*retrace from the maximum element*/
    i= max_i;
    for( i= max_i; i>F_max_i; i--) map_i2j[i-1] = far_far_away;;
    for( j= max_j; j>F_max_j; j--) map_j2i[j-1] = far_far_away;;
    i = F_max_i;
    j = F_max_j;
    *aln_score = F[i][j]; 
    while ( i>0 ||  j >0 ) {
	//printf (" %4d  %4d  %8.3f  \n", i, j, F[i][j]);
	switch ( direction[i][j] ) {
	case 'd':
	    //printf ( " %4d  %4d \n",  i, j);
	    map_i2j [i-1] = j-1;
	    map_j2i [j-1] = i-1;
	    i--;
	    j--; 
	    break;
	case 'i':
	    //printf ( " %4d  %4d \n",  i, -1);
	    map_i2j [i-1] = far_far_away;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = far_far_away;
	    j--; 
	    break; 
	default: 
	    fprintf (stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
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
int find_Calpha ( Protein *protein, int  resctr, double ca[3] ){

    Residue * res = protein->sequence + resctr;
    int i;
    
    for (i=0; i<res->no_atoms; i++) {
	if ( ! res->atom[i].backbone) continue;
	if ( ! strcmp(res->atom[i].type, "CA")) {
	    ca[0]= res->atom[i].x;
	    ca[1]= res->atom[i].y;
	    ca[2]= res->atom[i].z;
	    break;
	}
    }


    return 0;

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
