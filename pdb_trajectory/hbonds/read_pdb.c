# include "st_hbonds.h"

# define BUFFLEN 150
# define MAX_CHAINS 50

int read_pdb ( char * pdbname, Protein * protein, int * chain_number,
	       Atom ** solvent_ptr, int * solvent_number,
	       Atom ** ligand_ptr, int * ligand_number,
	       char ligand_type [][3], int no_ligand_types,
	       int invert,  int *done) {

    static Residue ** sequence;
    Atom * solvent;
    Atom * ligand;
    static FILE * fptr = NULL;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+2]; /* res name: 4 digits + insertion code + \0 */
    char oldrestype [PDB_ATOM_RES_NAME_LEN+2];
    char tmp[BUFFLEN], *auxptr;
    int atomctr, resctr, no_res[MAX_CHAINS] = {0}, ctr, nonblank;
    int chain_ctr, solvent_ctr, ligand_ctr, i, line_ctr, is_ligand;
    char old_chain_id, chain_id;

    
    /* open file, if not already opened */
    if ( ! fptr ) {
	fptr = fopen ( pdbname, "r");
	if ( !fptr ) {
	    fprintf (stderr, "Cno %s.\n", pdbname);
	    return 1;
	}

	/* count residues */
	memset (line,  0, BUFFLEN);
	memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
	resctr = 0;
	chain_ctr = 0;
	solvent_ctr = 0;
	ligand_ctr = 0;
	old_chain_id= '\0';
	line_ctr = 0;
	while(fgets(line, BUFFLEN, fptr)!=NULL){
	    line_ctr++;
	    if (  ! strncmp(line,"END", 3) ) break;
	    if(  strncmp(line,"ATOM", 4) && strncmp(line,"HETATM", 6) ) continue;
	    
	    /* solvent?  ion? */
	    if (  ! strncmp (line+PDB_ATOM_RES_NAME, "SOL",  PDB_ATOM_RES_NAME_LEN) ||
	        	! strncmp (line+PDB_ATOM_RES_NAME, " Na",  PDB_ATOM_RES_NAME_LEN) ||
		 ! strncmp (line+PDB_ATOM_RES_NAME, " Cl",  PDB_ATOM_RES_NAME_LEN) ) {
		solvent_ctr ++;
		continue;
	    }
	    /* ligand?*/
	    is_ligand = 0;
	    for ( i=0; i< no_ligand_types; i++) {
		if(  !   strncmp (line+PDB_ATOM_RES_NAME, ligand_type[i],  PDB_ATOM_RES_NAME_LEN)   ) {
		    ligand_ctr ++;
		    is_ligand = 1;
		    break;
		}
	    }
	    if (is_ligand) continue;
	    
	    /* which chain? */
	    chain_id = line[21];
	    if (old_chain_id != chain_id ) {
		if ( chain_ctr ) no_res[chain_ctr-1] = resctr;
		resctr = 0;
		chain_ctr ++;
		if ( chain_ctr > MAX_CHAINS) {
		    fprintf ( stderr, "%s  %c  %c \n", line, chain_id, old_chain_id);
		    fprintf ( stderr, "Line %d: More than %d,  protein chains in pdb .\n", line_ctr, MAX_CHAINS);
		    return 1;
		}
		old_chain_id = chain_id;
	    } 

	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		    
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		oldresno[PDB_ATOM_RES_NO_LEN+1] = '\0';
		/* printf ( "New residue number:  %s \n", oldresno); */
		resctr ++;
	    }
	}

	/* printf ("no residues: %d\n", no_res); */

	/* allocate space */
	if ( !chain_ctr ) return 1;
	no_res[chain_ctr-1] = resctr;
	
	
	if ( ! (sequence = emalloc (chain_ctr*sizeof (Residue*))) ) exit (1);
	for (i=0; i<chain_ctr; i++) {
	    if ( ! (sequence[i] = emalloc ( no_res[i]*sizeof (Residue))) ) exit (1);
	}
	/*return values: */
	for (i=0; i<chain_ctr; i++) {
	    (protein+i)->sequence = sequence[i];
	    (protein+i)->length   = no_res[i];
	}
	solvent = emalloc ( solvent_ctr*sizeof(Atom) );
	ligand  = emalloc ( ligand_ctr*sizeof(Atom) );
	* solvent_number = solvent_ctr;
	* ligand_number = ligand_ctr;
	* solvent_ptr = solvent;
	* ligand_ptr = ligand;
	* chain_number = chain_ctr;

	printf ("number of chains : %d   solvent atoms: %d   ligand atoms: %d \n",
		*chain_number, *solvent_number, *ligand_number);

	rewind (fptr);
	
    } else {
	for (i=0; i<*chain_number; i++) {
	    sequence[i] = (protein+i)->sequence;
	    no_res[i]   = (protein+i)->length;
	}
	solvent = * solvent_ptr;
	ligand  = * ligand_ptr;
    }

    /* read in the atoms */
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2); /*  tyring to account for the insertion code */
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+2);
    resctr= -1;
    atomctr = 0;
    chain_ctr = -1;
    solvent_ctr = 0;
    ligand_ctr = 0;
    old_chain_id= '\0';
    line_ctr = 0;
    
    while(fgets(line, BUFFLEN, fptr)!=NULL){

	
	line_ctr++;
	if (  ! strncmp(line,"END", 3) ) break;
	if(  strncmp(line,"ATOM", 4) && strncmp(line,"HETATM", 6) ) continue;
	    
	if ( line[PDB_ATOM_ATOM_NAME] == 'H' ||  line[PDB_ATOM_ATOM_NAME+1] == 'H') continue;
	    
	/* solvent? ion? */
	if (  ! strncmp (line+PDB_ATOM_RES_NAME, "SOL",  PDB_ATOM_RES_NAME_LEN) ||
		! strncmp (line+PDB_ATOM_RES_NAME, " Na",  PDB_ATOM_RES_NAME_LEN) ||
		 ! strncmp (line+PDB_ATOM_RES_NAME, " Cl",  PDB_ATOM_RES_NAME_LEN) ) {

	    
	    
	    strncpy ( tmp, line+PDB_ATOM_ATOM_NO, PDB_ATOM_ATOM_NO_LEN);
	    tmp[PDB_ATOM_ATOM_NO_LEN] = '\0';
	    solvent[solvent_ctr].pdb_id=atoi(tmp);
	    
	    
	    strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	    tmp[PDB_ATOM_X_LEN] = '\0';
	    solvent[solvent_ctr].x=atof(tmp);
	    
	    strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	    tmp[PDB_ATOM_Y_LEN] = '\0';
	    solvent[solvent_ctr].y=atof(tmp);

	    strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	    tmp[PDB_ATOM_Z_LEN] = '\0';
	    solvent[solvent_ctr].z=atof(tmp);
	    
	    solvent_ctr ++;
	    continue;
	}
	/* ligand?*/
	is_ligand = 0;
	for ( i=0; i< no_ligand_types; i++) {
	    if(  ! strncmp (line+PDB_ATOM_RES_NAME, ligand_type[i], PDB_ATOM_RES_NAME_LEN)   ) {
		
		strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
		tmp[PDB_ATOM_X_LEN] = '\0';
		ligand[ligand_ctr].x=atof(tmp);
	    
		strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
		tmp[PDB_ATOM_Y_LEN] = '\0';
		ligand[ligand_ctr].y=atof(tmp);

		strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
		tmp[PDB_ATOM_Z_LEN] = '\0';
		ligand[ligand_ctr].z=atof(tmp);
	    
	    
		ligand_ctr ++;
		is_ligand = 1;
		break;
	    }
	}
	if (is_ligand) continue;
	
	/* which chain? */
	chain_id = line[21];
	if (old_chain_id != chain_id ) {
	    resctr = -1;
	    chain_ctr ++;
	    old_chain_id = chain_id;
	} 

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
		
	    sequence[chain_ctr][resctr].no_atoms = 1;
	    strncpy ( sequence[chain_ctr][resctr].pdb_id, oldresno, PDB_ATOM_RES_NO_LEN+2);
	    sequence[chain_ctr][resctr].pdb_id[PDB_ATOM_RES_NO_LEN+1]   = '\0';
		
	    strncpy ( sequence[chain_ctr][resctr].res_type, oldrestype, PDB_ATOM_RES_NAME_LEN+1);
	    sequence[chain_ctr][resctr].res_type[PDB_ATOM_RES_NAME_LEN] = '\0';
	    sequence[chain_ctr][resctr].res_type_short  = single_letter ( sequence[chain_ctr][resctr].res_type );
	    if ( !sequence[chain_ctr][resctr].res_type_short ) return 1;
	   
	} else {
	    atomctr ++;
	    sequence[chain_ctr][resctr].no_atoms = atomctr + 1;
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

	strncpy ( sequence[chain_ctr][resctr].atom[atomctr].type, tmp, PDB_ATOM_ATOM_NAME_LEN );

	strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	tmp[PDB_ATOM_X_LEN] = '\0';
	sequence[chain_ctr][resctr].atom[atomctr].x=atof(tmp);
	strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	tmp[PDB_ATOM_Y_LEN] = '\0';
	sequence[chain_ctr][resctr].atom[atomctr].y=atof(tmp);
	strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	tmp[PDB_ATOM_Z_LEN] = '\0';
	sequence[chain_ctr][resctr].atom[atomctr].z=atof(tmp);
	   
	
    }

 
    /* clean PDB id tags from spaces
    for (chain_ctr=0; chain_ctr<2; chain_ctr++) {
	for (resctr=0; resctr < no_res[chain_ctr]; resctr ++ ) {
	    string_clean (sequence[chain_ctr][resctr].pdb_id, PDB_ATOM_RES_NO_LEN+1);
	}
    } */

    
    /* close file, if end */
    *done =  feof (fptr);
    if ( *done ) fclose (fptr);


    return 0;
}

