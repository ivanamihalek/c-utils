# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include "pdb.h"
# include "geometry.h"
# include "utils.h"
/* read in pdb file and sever chain id's (maybe),
   and output all waters within CUTOFF Angstroms from the chain(s),
   in gro format and some ad hoc protons attached*/

# define CUTOFF 7.0

int main ( int argc, char * argv[]) {

    Atom* solvent;
    int solvent_number;
    int gro_output (Atom*  solvent, int solvent_number);
    int read_pdb (char *argv[], int argc, Atom **solvent_ptr, int * solvent_number_ptr);
    
    if ( argc < 2 ) {
	fprintf (stderr, "Usage: %s <pdb file name> [<chain_id 1> <chain_id 2> ...] \n", argv[0]);
	exit (1);
    }

    if (  read_pdb (argv+1, argc-1, &solvent, &solvent_number) ){
	fprintf (stderr, "Error reading %s.\n", argv[1] );
	exit (1);
    }

    gro_output ( solvent, solvent_number);
    return 0;
}


/****************************************************/
/****************************************************/
# define BUFFLEN 150

int read_pdb (char *argv[], int argc, Atom **solvent_ptr, int * solvent_number_ptr) {

    FILE *fptr = efopen ( argv[0], "r");
    Residue * sequence;
    Atom * solvent, *atomptr;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+2]; /* res name: 4 digits + insertion code + \0 */
    char oldrestype [PDB_ATOM_RES_NAME_LEN+2];
    char tmp[PDB_ATOM_X_LEN+1], *auxptr;
    char chain_id;
    int resctr, solvent_ctr, chain_ok, i, no_res;
    int contact_found;
    int ctr, atomctr, nonblank, near_waters;
    double cutoff_dist_sq = CUTOFF*CUTOFF;
    double dist, aux;
    double x, y, z;
    
    atomctr = 0;

    if ( ! fptr ) return 1;
    
    /* count residues */
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
    resctr = 0;
    solvent_ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if (  ! strncmp(line,"END", 3) ) break;
	if(  strncmp(line,"ATOM", 4) && strncmp(line,"HETATM", 6) ) continue;
	    
	/* solvent? */
	if (  ! strncmp (line+PDB_ATOM_RES_NAME, "HOH",  PDB_ATOM_RES_NAME_LEN)   ) {
	    solvent_ctr ++;
	    continue;
	}
	/* other ligand*/
	if ( ! strncmp(line,"HETATM", 6)) continue;
	    
	/* which chain? */
	chain_id = line[21];
	if ( argc > 1 ) {
	    chain_ok = 0;
	    for (i=1; i<argc && ! chain_ok; i++ ) {
		chain_ok = ( argv[i][0] == chain_id ) ;
	    }
	    if ( ! chain_ok ) continue;
	} 

	
	if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		    
	    strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
	    oldresno[PDB_ATOM_RES_NO_LEN+1] = '\0';
	    /* printf ( "New residue number:  %s \n", oldresno); */
	    resctr ++;
	}
    }
    no_res = resctr;
    //printf ("no residues: %d  no of waters: %d \n", no_res, solvent_ctr);

    /* allocate space */
    sequence = emalloc ( no_res*sizeof (Residue));
    * solvent_ptr = solvent = emalloc ( solvent_ctr*sizeof(Atom) );
    * solvent_number_ptr = solvent_ctr;

    rewind (fptr);
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2); /*  tyring to account for the insertion code */
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+2);
    resctr = -1;
   /* read in residues */
    while(fgets(line, BUFFLEN, fptr)!=NULL){

	
	if (  ! strncmp(line,"END", 3) ) break;
	if(  strncmp(line,"ATOM", 4)  ) continue;
	
	/* which chain? */
	chain_id = line[21];
	if ( argc > 1 ) {
	    chain_ok = 0;
	    for (i=1; i<argc && ! chain_ok; i++ ) {
		chain_ok = ( argv[i][0] == chain_id ) ;
	    }
	    if ( ! chain_ok ) continue;
	}
	    
	if ( line[PDB_ATOM_ATOM_NAME] == 'H' ||  line[PDB_ATOM_ATOM_NAME+1] == 'H') continue;
	/***********/
	/* protein */
	/***********/
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
	    sequence[resctr].res_type_short  = single_letter ( sequence[resctr].res_type );
	    if ( !sequence[resctr].res_type_short ) return 1;
	   
	} else {
	    atomctr ++;
	    sequence[resctr].no_atoms = atomctr + 1;
	    if ( atomctr >= MAX_NO_ATOMS ) {
		fprintf ( stderr, "Error: I thought every aa has < %d atoms.\n",
			  MAX_NO_ATOMS );
		return 1;
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
    

    rewind (fptr);
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2); /*  tyring to account for the insertion code */
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+2);
    /* read in solvent -- if within cutoff from a residue*/
    solvent_ctr = 0;
    near_waters  = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){

	
	if (  ! strncmp(line,"END", 3) ) break;
	if(  strncmp(line,"HETATM", 6) ) continue;
	    
	if (  strncmp (line+PDB_ATOM_RES_NAME, "HOH",  PDB_ATOM_RES_NAME_LEN)   ) continue;
	strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	tmp[PDB_ATOM_X_LEN] = '\0';
	x = solvent[solvent_ctr].x = atof(tmp);
	strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	tmp[PDB_ATOM_Y_LEN] = '\0';
	y = solvent[solvent_ctr].y = atof(tmp);
	strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	tmp[PDB_ATOM_Z_LEN] = '\0';
	z = solvent[solvent_ctr].z = atof(tmp);

	contact_found = 0;
	for ( resctr=0; resctr<no_res && !contact_found; resctr++ ) {
	    for (atomctr=0; atomctr < sequence[resctr].no_atoms && !contact_found; atomctr++) {
		atomptr = sequence[resctr].atom + atomctr;
		aux = atomptr->x - x;
		dist = aux*aux;
		aux = atomptr->y - y;
		dist += aux*aux;
		aux = atomptr->z - z;
		dist += aux*aux;
		if (dist <= cutoff_dist_sq) {
		    solvent[solvent_ctr].near = 1;
		    contact_found = 1;
		    near_waters ++;
		}
	    }
	}
	solvent_ctr++;
	
	
    }

    fclose(fptr);

    //printf ("near waters: %d\n", near_waters);
    
    return 0;
}


/****************************************************/
/****************************************************/
int gro_output (Atom*  solvent, int solvent_number) {

    int ctr, atom_ctr, solvent_ctr;
    double x, y, z;

    ctr = 0;
    atom_ctr = 0;

    for (solvent_ctr =0; solvent_ctr < solvent_number; solvent_ctr++ ) {
	
	if ( ! (solvent+solvent_ctr)->near ) continue;

	ctr ++;
	
	x = (solvent+solvent_ctr)->x;
	y = (solvent+solvent_ctr)->y;
	z = (solvent+solvent_ctr)->z;
	atom_ctr ++;
	printf ("%5d%5s%5s%5d%8.3lf%8.3lf%8.3lf\n",
		ctr, "SOL",  "OW", atom_ctr, x/10, y/10, z/10);

	z  = (solvent+solvent_ctr)->z + 1.0;
	atom_ctr ++;
	printf ( "%5d%5s%5s%5d%8.3lf%8.3lf%8.3lf\n",
	    ctr, "SOL",  "HW1", atom_ctr, x/10, y/10, z/10);

	x  = (solvent+solvent_ctr)->x + 0.948;
	z  = (solvent+solvent_ctr)->z - 0.33;
	atom_ctr ++;
	printf( "%5d%5s%5s%5d%8.3lf%8.3lf%8.3lf\n",
	    ctr, "SOL",  "HW2", atom_ctr, x/10, y/10, z/10);
	
	
    }
    return 0;
}
