# include "fitter.h"

int transform_pdb (double  **tfm_matrix, double * transl_vector,
		   Residue * sequence, int no_res, Residue * new_sequence) ;


/***************************************************************************/
/***************************************************************************/
int pdb_output ( char *pdb_fnm_1, char *pdb_fnm_2, char * out_fnm, double  **tfm_matrix,
		 double * transl_vector, Residue * sequence, int no_res) {

    int i, j;
    int resctr, atomctr, serial;
    FILE * fptr;
    Atom * atomptr,* new_atomptr;
    char atomtype[BUFFLEN] = {'\0'};

    /* rotate the chain: */
    /* alocate           */
    Residue * new_sequence = NULL;
    if ( ! (new_sequence = emalloc ( no_res*sizeof (Residue)))) {
	fprintf (stderr, "Error allocating.\n");
	exit (1);
    }
    /* rotate */
    transform_pdb (tfm_matrix, transl_vector,  sequence, no_res, new_sequence);
    /* open file */
    fptr = fopen ( out_fnm, "w");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", out_fnm);
	return 1;
    }
    /* output the trfm matrix, for the reference */
    fprintf (fptr, "REMARK \n" );
    fprintf (fptr, "REMARK  %s rotated onto %s \n", pdb_fnm_1, pdb_fnm_2);
    fprintf (fptr, "REMARK  tfm matrix: \n" );
    fprintf (fptr, "REMARK \n" );
    for (i=0; i<3; i++) {
	fprintf (fptr, "REMARK ");
	for (j=0; j<3; j++) {
	    fprintf (fptr, "  %8.3lf ", tfm_matrix[i][j]);
	}
	fprintf (fptr, "  %8.3lf \n", transl_vector[i]);	
    }
    fprintf (fptr, "REMARK \n");
    fprintf (fptr, "REMARK \n");

    serial = 0;
    for ( resctr= 0; resctr < no_res; resctr ++ ) {
	for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
	    serial ++;
	    /* 'new' is rotated; the old seqeunce contains the atomtype info */
	    new_atomptr = new_sequence[resctr].atom + atomctr;
	    atomptr     = sequence[resctr].atom + atomctr;
	    memset (atomtype, 0, BUFFLEN);
	    if ( sequence[resctr].is_hetatm ){
		sprintf (atomtype, "%s", "HETATM");
	    } else {
		sprintf (atomtype, "%s", "ATOM");
	    }
	    fprintf (fptr,  "%-6s%5d  %-3s%1s%-3s %1c%4s%1s   %8.3lf%8.3lf%8.3lf\n", 
		     atomtype,  serial ,  atomptr->type,   " ",   sequence[resctr].res_type,
		     sequence[resctr].chain, sequence[resctr].pdb_id,   " " ,
		     new_atomptr->x,  new_atomptr->y,   new_atomptr->z);
	}
    }

    fclose (fptr);
    
    /* free   */
    free(new_sequence);
    return 0;
}


/*******************************************************************************/
int transform_pdb (double  **tfm_matrix, double * transl_vector,
	       Residue * sequence, int no_res, Residue * sequence_new) {
    
    int resctr, atomctr;
    int i,j;
    double  x[3], new_x[3];
    
    for ( resctr= 0; resctr < no_res; resctr ++ ) {
	for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
	    x[0] = sequence[resctr].atom[atomctr].x;
	    x[1] = sequence[resctr].atom[atomctr].y;
	    x[2] = sequence[resctr].atom[atomctr].z;
	    for (i=0; i <3; i++ ) {
		new_x[i] = 0.0;
		for (j=0; j <3; j++ ) {
		    new_x[i] += tfm_matrix[i][j]*x[j];
		}
		new_x[i] += transl_vector[i];
	    }
	    sequence_new[resctr].atom[atomctr].x = new_x[0];
	    sequence_new[resctr].atom[atomctr].y = new_x[1];
	    sequence_new[resctr].atom[atomctr].z = new_x[2];
	}
    }
    return 0;
}
