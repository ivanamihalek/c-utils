/* average the rho score over the neghborhood of a residue */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
#include <ctype.h>

#include "pdb.h"
#include "utils.h"

# define CUTOFF_DIST  5.0
# define CUTOFF_CVG   0.6
# define MAX_LAYER   15

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int layer;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 40

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+1];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
    double rho, cvg, *radius, *avg_as_fn_of_radius;
    int interface;
} Residue;

Residue * sequence;
int no_res;


#define DEFAULT_RADIUS 5.0

int main ( int argc, char * argv[]) {
    

    int read_pdb ( char * pdbname, Residue ** sequence, int *no_res);
    int read_ranks ( char * pdbname, Residue * sequence, int no_res);
    int find_average (double *radius, int no_of_r);
 
    char pdbname [150] = "\0";
    char ranksfile [150] = "\0";
 
    double *radius;
    int no_of_r, i;
    

    if ( argc < 3 ) {
	printf ( "Usage: %s <pdb_file> <ranks_sorted> [ <r1> ..<rn>].\n", argv[0] );
	exit (1);
    } else {
	sprintf ( pdbname, "%s", argv[1]);
	sprintf ( ranksfile, "%s", argv[2]);
    }

    no_of_r = ( argc > 3) ? argc - 3 : 1;
 
    if ( !(radius = (double *) emalloc ( no_of_r* sizeof(double)) ) )exit(1);
    
    if ( argc > 3) {
	for (i=0; i< no_of_r; i++ ) {
	    radius[i] = atof(argv[3+i]);
	}
    } else {
	radius[0] = DEFAULT_RADIUS;
    }
    
    /* read in the pdbs */
    if ( read_pdb ( pdbname, &sequence, &no_res) )exit (1);
    
    /* read in the ranks */
    if ( read_ranks ( ranksfile, sequence, no_res )) exit (1);


# if 0
    /* test the input*/
     int resctr;
     for ( resctr=0; resctr<no_res; resctr++) {
	printf ("resctr: %d     pdbid: %s    restype:  %s    no atoms: %d    rho: %8.3lf   cvg: %8.3lf\n",
		resctr, sequence[resctr].pdb_id, sequence[resctr].res_type,
		sequence[resctr].no_atoms, 	sequence[resctr].rho, 	sequence[resctr].cvg);
    }
    exit (0);
# endif
    
    /* calculate the score*/
    find_average(radius, no_of_r);
 
    
    return 0;
}



/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int find_average (double *radius, int no_of_r) {


    int resctr1, resctr2;
    int atomctr1, atomctr2;
    int i;
    Atom * atomptr1, * atomptr2;
    double dist, min_dist, aux;
    double ** avg_rho, ** avg_cvg;
    int ** counter;
    
    avg_rho = dmatrix ( 0, no_res, 0, no_of_r);
    avg_cvg = dmatrix ( 0, no_res, 0, no_of_r);
    counter = imatrix ( 0, no_res, 0, no_of_r);
    
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	
	for ( resctr2=resctr1+1; resctr2<no_res; resctr2++) {

	    min_dist = 20000;
	    for (atomctr1=0; atomctr1 < sequence[resctr1].no_atoms; atomctr1++) {
		atomptr1 = sequence[resctr1].atom + atomctr1;
		/* find  dist btw atoms*/
		for (atomctr2=0; atomctr2 < sequence [resctr2].no_atoms; atomctr2++) {
		    atomptr2 = sequence[resctr2].atom + atomctr2;
		    
		    aux = atomptr1->x - atomptr2->x;
		    dist = aux*aux;
		    aux = atomptr1->y - atomptr2->y;
		    dist += aux*aux;
		    aux = atomptr1->z - atomptr2->z;
		    dist += aux*aux;
		    dist = sqrt (dist);
		    if (dist <= min_dist) {
			min_dist = dist;;
		    }
		}
	    }	    

	    for (i=0; i<no_of_r; i++) {
		if ( min_dist < radius[i] ) {
		    avg_rho[resctr1][i] += sequence[resctr2].rho;
		    avg_cvg[resctr1][i] += sequence[resctr2].cvg;
		    counter[resctr1][i] += 1;
		    avg_rho[resctr2][i] += sequence[resctr1].rho;
		    avg_cvg[resctr2][i] += sequence[resctr1].cvg;
		    counter[resctr2][i] += 1;
		}
	    }

	    
	} /* loop over resctr2 */
    }/* loop over resctr1 */

    printf ("%6s%6s%6s%10s%10s", "seqno", " pdbid", "type", "rho ", "cvg ");
    for (i=0; i<no_of_r; i++) {
	printf ("%10s%10s", "avg rho", "avg cvg");
    }
    printf ("\n");
    printf ("%38s", "");
    for (i=0; i<no_of_r; i++) {
	printf ("%5s%5.2lf%5s%5.2lf", "r=", radius[i],  "r=", radius[i]);
    }
    printf ("\n");
    
		
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	 for (i=0; i<no_of_r; i++) {
	     avg_rho[resctr1][i] /= counter[resctr1][i];
	     avg_cvg[resctr1][i] /= counter[resctr1][i];
	 }
	 printf ("%6d%6s%6s%10.3lf%10.3lf",
		resctr1+1, sequence[resctr1].pdb_id, sequence[resctr1].res_type,
		 sequence[resctr1].rho, 	sequence[resctr1].cvg );
	 for (i=0; i<no_of_r; i++) {
	     printf ("%10.3lf%10.3lf",  avg_rho[resctr1][i], avg_cvg[resctr1][i]);
	 }
	 printf ("\n");
	
     }

     free_dmatrix (avg_rho, 0, no_res, 0, no_of_r);
     free_dmatrix (avg_cvg, 0, no_res, 0, no_of_r);
     free_imatrix (counter, 0, no_res, 0, no_of_r);
     return 0;
}




/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

# define BUFFLEN 150

int read_ranks ( char * filename, Residue * sequence, int no_res) {

    FILE * fptr = NULL;
    char line[BUFFLEN], formatstring [BUFFLEN];
    char pdbid[PDB_ATOM_RES_NO_LEN+1];
    char * aux_strptr;
    int retval, resctr, epitope, pss, surface, iface, skipfields, i;
    double rho, cvg;

    /* open file */
    /* sprintf ( ranksfile, "%s.ranks", pdbname); */
    fptr = fopen ( filename, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", filename);
	return 1;
    }
    epitope = 0;
    pss = 0;
    surface = 0;
    iface = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr (line, "residue") ) {
	    /* is there a column indicating epitope? */ 
	    if ( strstr (line, "epitope") ) epitope = 1; 
	    if ( strstr (line, "pss") ) pss = 1; 
	    if ( strstr (line, "surface") ) surface = 1; 
	    if ( strstr (line, "i-face") ) iface = 1; 
	    break;
	}
    }
    skipfields = epitope+surface+pss+iface;
    sprintf ( formatstring,"%s", " %*s  %s  %*s  %*d  %*d   %*s");
    for (i=1; i<=skipfields; i++ ) {
	sprintf ( formatstring,"%s %s ",formatstring, " %*s ");
     }
    sprintf ( formatstring,"%s %s ",formatstring, " %lf  %lf");
	
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	retval = sscanf ( line, formatstring, pdbid, &rho, &cvg);

	if ( retval != 3 ) continue;
	retval = (int) strstr ( pdbid, "-");
	if ( retval) continue;
	for ( resctr=0; resctr<no_res; resctr++) {
	    aux_strptr =  sequence[resctr].pdb_id;
	    while ( isspace ( *aux_strptr ) &&
		    aux_strptr < (sequence[resctr].pdb_id+PDB_ATOM_RES_NO_LEN+1) )
		aux_strptr++;
	    if (  ! (strcmp(  aux_strptr, pdbid)) ) {
		sequence[resctr].rho = rho;
		sequence[resctr].cvg = cvg;
	    }
	}
    }

    /* close file */
    fclose (fptr);
    return 0;
}







/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int read_pdb ( char * pdbname, Residue ** sequence_ptr, int * no_res_ptr) {

    Residue * sequence;
    FILE * fptr = NULL;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+1];
    char oldrestype [PDB_ATOM_RES_NAME_LEN+1];
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
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+1);
    resctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if( ! strncmp(line,"ATOM", 4)){
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN) ) {
		
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN);
		oldresno[PDB_ATOM_RES_NO_LEN] = '\0';
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
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+1);
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+1);
    resctr= -1;
    atomctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if( ! strncmp(line,"ATOM", 4)){
	    /* adjust the counters */ 
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN) ) {
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN);
		strncpy (oldrestype, line+PDB_ATOM_RES_NAME, PDB_ATOM_RES_NAME_LEN);
		oldresno[PDB_ATOM_RES_NO_LEN] = '\0';
		oldrestype[PDB_ATOM_RES_NAME_LEN] = '\0';
		resctr ++;
		atomctr = 0;
		
		sequence[resctr].no_atoms = 1;
		strncpy ( sequence[resctr].pdb_id,oldresno, PDB_ATOM_RES_NO_LEN);
		sequence[resctr].pdb_id[PDB_ATOM_RES_NO_LEN] = '\0';
		strncpy ( sequence[resctr].res_type,oldrestype, PDB_ATOM_RES_NAME_LEN);
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

