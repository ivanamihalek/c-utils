# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
#include <ctype.h>


#include "pdb.h"
#include "utils.h"

# define CUTOFF_DIST  4.0
# define CUTOFF_CVG   0.6
# define MAX_LAYER   15

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int layer;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 100 /* so I can handle things like heme */

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+2];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
    int interface;
} Residue;

Residue * sequence1, *sequence2;
int no_res_1, no_res_2;

#define RANKS  "sim.ranks"

int main ( int argc, char * argv[]) {
    

    int read_pdb ( char * pdbname, Residue ** sequence, int *no_res);
    int determine_noc_matrix ( int ** noc_matrix, Residue * sequence, int no_res, double cutoff_dist);
    int  find_distance (int **noc1, int **noc2, int no_res_1, int no_res_2, double *dist_ptr);
 
    char pdbname1 [150] = "\0";
    char pdbname2 [150] = "\0";
    int ** noc1, ** noc2;
    double cutoff_dist = CUTOFF_DIST;
    double dist;

    if ( argc < 3 ) {
	printf ( "Usage: %s <pdbname1> <pdbname2> [<cutoff_dist>] .\n", argv[0] );
	exit (1);
    } 
    sprintf ( pdbname1, "%s", argv[1]);
    sprintf ( pdbname2, "%s", argv[2]);
    /* printf  ("pdbnames: %s%s\n", pdbname1, pdbname2); */
    if ( argc >= 4 ) {
	cutoff_dist = atof ( argv[3]);
    }
   
    
    /* read in the pdbs */
    if ( read_pdb ( pdbname1, &sequence1, &no_res_1)) exit (1);
    if ( read_pdb ( pdbname2, &sequence2, &no_res_2)) exit (1);
    
    if ( no_res_1 != no_res_2 ) {
	fprintf ( stderr, "Error: %s and %s do not have the same no of res ( %d and %d).\n",
		  pdbname1, pdbname2, no_res_1, no_res_2);
	exit (1);
    }

# if 0
    /* test the input*/
     int resctr;
     for ( resctr=0; resctr<no_res_1; resctr++) {
	printf ("resctr: %d     pdbid: %s    restype:  %s    no atoms: %d  \n",
		resctr, sequence1[resctr].pdb_id, sequence1[resctr].res_type,
		sequence1[resctr].no_atoms);
    }
    exit (0);
# endif
    /* assgin noc storage space, and dtermine noc matrix */
    noc1 = imatrix (0, no_res_1 -1, 0, no_res_1 -1);
    noc2 = imatrix (0, no_res_2 -1, 0, no_res_2 -1);
    determine_noc_matrix ( noc1, sequence1, no_res_1, cutoff_dist);
    determine_noc_matrix ( noc2, sequence2, no_res_2, cutoff_dist);
    
    /* calculate the distance*/
    find_distance (noc1, noc2, no_res_1, no_res_2, &dist);
    printf ( "%8.3le\n", dist);
    
    free_imatrix ( noc1, 0, no_res_1 -1, 0, no_res_1 -1);
    free_imatrix ( noc2, 0, no_res_2 -1, 0, no_res_2 -1);
    return 0;
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int    find_distance (int **noc1, int **noc2, int no_res_1, int no_res_2, double *dist_ptr){
    double aux, dist;
    int no_res, resctr1, resctr2, ctr;;
    if ( no_res_1 != no_res_2 ) {
	fprintf ( stderr,
		  "Error in find)dist: noc matrices do not have the same no of res ( %d and %d).\n",
		   no_res_1, no_res_2);
	exit (1);
    }
    no_res = no_res_1;
    dist = 0;
    ctr = 0;
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	for ( resctr2=resctr1+1; resctr2<no_res; resctr2++) {
	    aux = noc1 [resctr2][resctr1] -  noc2 [resctr2][resctr1];
	    dist += aux*aux;
	    ctr++;
	}
    }
    if ( ctr) dist = sqrt (dist/ctr);
    * dist_ptr = dist;
    return 0;
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int determine_noc_matrix ( int ** noc_matrix, Residue * sequence, int no_res, double cutoff_dist){
    
    int resctr1, resctr2;
    int atomctr1, atomctr2;
    int n;
    double dist, aux;
    Atom * atomptr1, * atomptr2;

    /* fast return: */
    if ( !noc_matrix) {
	fprintf (stderr, "In determine_noc_matrix(): noc matrix space must be assigned.\n");
	return 1;
    }
    

    /*diagonal*/
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	noc_matrix [resctr1][resctr1] = 0;
    }

    for ( resctr1=0; resctr1<no_res; resctr1++) {
	
	for ( resctr2=resctr1+1; resctr2<no_res; resctr2++) {

	    n = 0;
	    
	    for ( atomctr1=0; atomctr1 < sequence[resctr1].no_atoms; atomctr1 ++) {
		atomptr1 = sequence[resctr1].atom + atomctr1;

		for ( atomctr2=0; atomctr2 < sequence[resctr2].no_atoms; atomctr2++) {
		    atomptr2 = sequence[resctr2].atom + atomctr2;

		    
		    dist = 0.0;
		    aux = atomptr1->x -  atomptr2->x;
		    dist += aux*aux;
		    aux = atomptr1->y -  atomptr2->y;
		    dist += aux*aux;

		    aux = atomptr1->z -  atomptr2->z;
		    dist += aux*aux;
		    dist = sqrt (dist);
		    if ( dist < cutoff_dist) {
			n++;
		    }
		}
		
	    }
	    noc_matrix [resctr2][resctr1] = noc_matrix [resctr1][resctr2] =  (n>0);
	}
    }
#if 0
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	for ( resctr2=0; resctr2<no_res; resctr2++) {
	    if ( noc_matrix [resctr1][resctr2] )
	    printf ("%4d  %4d   %5d \n", resctr1, resctr2, noc_matrix [resctr1][resctr2]);
	}
    }
    exit (0);
# endif    
    return 0;
}



#define BUFFLEN 250

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int read_pdb ( char * pdbname, Residue ** sequence_ptr, int * no_res_ptr) {

    Residue * sequence;
    FILE * fptr = NULL;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+2]; /* res name: 4 digits + insertion code + \0 */
    char oldrestype [PDB_ATOM_RES_NAME_LEN+2];
    char tmp[PDB_ATOM_X_LEN+1], *auxptr;
    int atomctr, resctr, no_res, ctr, nonblank;

    
    /* open file */
    fptr = fopen ( pdbname, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", pdbname);
	return 1;
    }

    /* count residues */
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
    resctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if( ! strncmp(line,"ATOM", 4) ||  ! strncmp(line,"HETATM", 6)){
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		oldresno[PDB_ATOM_RES_NO_LEN+1] = '\0';
		/* printf ( "New residue number:  %s \n", oldresno); */
		resctr ++;
	    }
	}
    }
    no_res = resctr;
    *no_res_ptr = no_res;
    /* printf ("no residues: %d\n", no_res); */

    /* allocate space */
    sequence = NULL;
    sequence = calloc ( no_res, sizeof (Residue));
    if ( ! sequence ) {
	fprintf ( stderr, "Error allcating sequence space.\n");
	exit (0);
    }
    *sequence_ptr = sequence;

    /* read in the atom */
    rewind ( fptr);
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2); /*  tyring to account for the insertion code */
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+2);
    resctr= -1;
    atomctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if( ! strncmp(line,"ATOM", 4)  ||  ! strncmp(line,"HETATM", 6)){
	    /* adjust the counters */ 
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		strncpy (oldrestype, line+PDB_ATOM_RES_NAME, PDB_ATOM_RES_NAME_LEN);
		oldresno[PDB_ATOM_RES_NO_LEN+1] = '\0';
		oldrestype[PDB_ATOM_RES_NAME_LEN] = '\0';
		resctr ++;
		atomctr = 0;
		
		sequence[resctr].no_atoms = 1;
		strncpy ( sequence[resctr].pdb_id, oldresno, PDB_ATOM_RES_NO_LEN+2);
		sequence[resctr].pdb_id[PDB_ATOM_RES_NO_LEN+1]   = '\0';
		
		strncpy ( sequence[resctr].res_type, oldrestype, PDB_ATOM_RES_NAME_LEN+1);
		sequence[resctr].res_type[PDB_ATOM_RES_NAME_LEN] = '\0';
	   
	    } else {
		atomctr ++;
		sequence[resctr].no_atoms = atomctr + 1;
		if ( atomctr >= MAX_NO_ATOMS ) {
		    fprintf ( stderr, "Error: I thought every aa has < %d atoms.\n",
			      MAX_NO_ATOMS );
		    exit (1);
		}
	    }
	    /* read in atom info */
	    
	    auxptr = line+ PDB_ATOM_ATOM_NAME;
	    memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	    /* skip initial blanks*/
	    ctr  = 0;
	    while ( !(isalpha (*(auxptr + ctr))) &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	    /* copy alphanum info */
	    nonblank = 0;
	    while (  isalpha (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
		tmp[nonblank] =  *(auxptr +ctr);
		nonblank ++;
		ctr++;
	    }

	    strncpy ( sequence[resctr].atom[atomctr].type, tmp, PDB_ATOM_ATOM_NAME_LEN );

	    /* is this a backbone atom?*/
	    sequence[resctr].atom[atomctr].backbone = 0;
	    if ( nonblank == 1) {
		  sequence[resctr].atom[atomctr].backbone =
		      !(  strcmp ( tmp, "N") && strcmp ( tmp, "C") && strcmp ( tmp, "O")  );
	    } else if (  nonblank == 2) {
		  sequence[resctr].atom[atomctr].backbone = ! strcmp ( tmp, "CA" );
	    }
	    /* printf ( " %4d %4d %4s is backbone: %1d \n", resctr, atomctr, */
		     /* sequence[resctr].atom[atomctr].type, sequence[resctr].atom[atomctr].backbone); */
	    strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	    tmp[PDB_ATOM_X_LEN] = '\0';
	    sequence[resctr].atom[atomctr].x=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	    tmp[PDB_ATOM_Y_LEN] = '\0';
	    sequence[resctr].atom[atomctr].y=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	    tmp[PDB_ATOM_Z_LEN] = '\0';
	    sequence[resctr].atom[atomctr].z=atof(tmp);
	   
	}
	
    }

    /* close file */
    fclose (fptr);
    return 0;
}

