# include "ift.h"

# define BUFFLEN 150



int main ( int argc, char * argv[]) {

    int done, first, frame_ctr, solctr;
    int avg_avg_dist     = 0;
    int avg_avg_dist_sq  = 0;
    int avg_avg_drift    = 0;
    int avg_avg_drift_sq = 0;
    int solvent_number;
    Atom * solvent = NULL, *solvent_init = NULL;
    double * avg_distance_from_init = NULL,  *solvent_neighbors = NULL;
    double * avg_drift = NULL, * old_drift  = NULL;
    
    double z_dist, sigma_dist;
    double z_drift, sigma_drift;
    FILE * fptr;
    int read_pdb ( char * pdbname, Atom ** solvent_ptr,
		   int * solvent_number, int *done);
    int water_motion (Atom * solvent, int solvent_number, Atom * solvent_init,
		  double * avg_distance_from_init, double * solvent_neighbors);
    if ( argc < 2 ) {
	fprintf (stderr, "Usage: %s <pdb file name> \n", argv[0]);
	exit (1);
    }


    done = 0;
    first = 1;
    frame_ctr = 0;
    printf ( "\n");
    while ( !done ) {
	if ( first ) {
	    printf ("frame: %4d\n", frame_ctr+1);
	} else {
	    printf ("frame: %4d\r", frame_ctr+1); fflush(stdout);
	}
	if (  read_pdb (argv[1], &solvent, &solvent_number,  &done) ){
	    fprintf (stderr, "Error reading %s.\n", argv[1] );
	    exit (1);
	}
	if (done) break;
	if ( first ) {


	    if ( ! ( solvent_init = emalloc((solvent_number*sizeof(Atom))) ) ) exit(1);
	    
	    memcpy ( solvent_init, solvent, solvent_number*sizeof(Atom));
	    if ( ! (avg_distance_from_init = emalloc((solvent_number*sizeof(double)) ) ) ) exit(1);
	    if ( ! (solvent_neighbors = emalloc(solvent_number*sizeof(double))) ) exit(1);
	    if ( ! (avg_drift = emalloc(solvent_number*sizeof(double))) ) exit(1);
	    if ( ! (old_drift = emalloc(solvent_number*sizeof(double))) ) exit(1);
	}
	
	first = 0;

	water_motion ( solvent, solvent_number, solvent_init, avg_distance_from_init, solvent_neighbors );

	frame_ctr++;
	
	for ( solctr=0; solctr<solvent_number; solctr++) {
	    avg_drift[solctr] += avg_distance_from_init[solctr]/frame_ctr - old_drift[solctr];
	    old_drift[solctr]  = avg_distance_from_init[solctr]/frame_ctr;
	}
	
	//ne = ( frame_ctr==20);
    }
    printf ("\n");

    avg_avg_dist     = 0;
    avg_avg_dist_sq  = 0;
    avg_avg_drift    = 0;
    avg_avg_drift_sq = 0;
    
    for ( solctr=0; solctr<solvent_number; solctr++) {
	
	avg_distance_from_init[solctr] /= frame_ctr;
	avg_avg_dist += avg_distance_from_init[solctr];
	avg_avg_dist_sq += avg_distance_from_init[solctr]*avg_distance_from_init[solctr];
	
	solvent_neighbors[solctr]      /= frame_ctr;
	
	avg_drift[solctr]             /= frame_ctr;
	avg_avg_drift += avg_drift[solctr] ;
	avg_avg_drift_sq += avg_drift[solctr]*avg_drift[solctr];
    }
    avg_avg_dist /= solvent_number;
    avg_avg_dist_sq /= solvent_number;
    avg_avg_drift /= solvent_number;
    avg_avg_drift_sq /= solvent_number;

    sigma_dist = 0;
    if ( avg_avg_dist_sq >=avg_avg_dist*avg_avg_dist ) {
	if (  avg_avg_dist_sq -  avg_avg_dist*avg_avg_dist > 1.e-6 ) {
	    sigma_dist = sqrt (  avg_avg_dist_sq -  avg_avg_dist*avg_avg_dist );
	}
    } else {
	fprintf ( stderr, "error finding z_score\n");
	exit (1);
    }


    sigma_drift = 0;
    if ( avg_avg_drift_sq >=avg_avg_drift*avg_avg_drift ) {
	if (  avg_avg_drift_sq -  avg_avg_drift*avg_avg_drift > 1.e-6 ) {
	    sigma_drift = sqrt (  avg_avg_drift_sq -  avg_avg_drift*avg_avg_drift );
	}
    } else {
	fprintf ( stderr, "error finding z_score\n");
	exit (1);
    }

    /* solvent */ 
    if ( ! ( fptr = efopen ( "solv.out", "w") ) )  exit (1);
    for ( solctr=0; solctr<solvent_number; solctr++) {
	if ( sigma_dist ) {
	    z_dist =  (avg_distance_from_init[solctr] - avg_avg_dist)/sigma_dist;
	} else {
	    z_dist = 0;
	}
	if ( sigma_drift ) {
	    z_drift =  (avg_drift[solctr] - avg_avg_drift)/sigma_drift;
	} else {
	    z_drift = 0;
	}
	// ( z_dist < -3 || z_drift < -3 ) {
	    fprintf ( fptr, " %5d  %8.3lf   %6.1le    %8.3lf   %6.1le   %8.3lf\n", (solvent+solctr)->pdb_id,
		 avg_distance_from_init[solctr],  z_dist,   avg_drift[solctr],  z_drift, solvent_neighbors[solctr] );
	    //
    }
    fclose (fptr);
    
    return 0;
   
}
/*******************************************************************************/
int water_motion (Atom * solvent, int solvent_number, Atom * solvent_init,
		  double * avg_distance_from_init, double * solvent_neighbors) {

    int solctr;
    Atom * atomptr;
    double x, y, z, aux, dist;
    
     for ( solctr=0; solctr<solvent_number; solctr++) {

	x= (solvent+solctr)->x;
	y= (solvent+solctr)->y;
	z= (solvent+solctr)->z;
	
	/* distance from initial position */

	atomptr = solvent_init + solctr;
	
	aux = x - atomptr->x;
	dist = aux*aux;
	aux = y - atomptr->y;
	dist += aux*aux;
	aux = z - atomptr->z;
	dist += aux*aux;
	avg_distance_from_init[solctr] += sqrt (dist);
	
    }
	    
   return 0;
}


/*******************************************************************************/

int read_pdb ( char * pdbname, 
	       Atom ** solvent_ptr, int * solvent_number,
	       int *done) {

    Atom * solvent;
    static FILE * fptr = NULL;
    char line[BUFFLEN];
    char tmp[BUFFLEN];
    int  solvent_ctr;

    
    /* open file, if not already opened */
    if ( ! fptr ) {
	fptr = fopen ( pdbname, "r");
	if ( !fptr ) {
	    fprintf (stderr, "Cno %s.\n", pdbname);
	    return 1;
	}

	/* count atomss */
	memset (line,  0, BUFFLEN);
	solvent_ctr = 0;
	while(fgets(line, BUFFLEN, fptr)!=NULL){

	    if (  ! strncmp(line,"END", 3) ) break;
	    if(  strncmp(line,"ATOM", 4)  ) continue;
	    
	    /* solvent?  ion? */
	    if (  ! strncmp (line+PDB_ATOM_RES_NAME, "SOL",  PDB_ATOM_RES_NAME_LEN)  ) {
		solvent_ctr ++;
		continue;
	    }
	}

	/* printf ("no residues: %d\n", no_res); */
	solvent = emalloc ( solvent_ctr*sizeof(Atom) );
	* solvent_number = solvent_ctr;

	printf (" solvent atoms: %d  \n",  *solvent_number);

	rewind (fptr);

	* solvent_ptr = solvent;
	
    } else {
	solvent = * solvent_ptr;
    }

    /* read in the atoms */
    memset (line,  0, BUFFLEN);
    solvent_ctr = 0;
    
    while(fgets(line, BUFFLEN, fptr)!=NULL){

	if (  ! strncmp(line,"END", 3) ) break;
	if(  strncmp(line,"ATOM", 4) ) continue;
	    
	if ( line[PDB_ATOM_ATOM_NAME] == 'H' ||  line[PDB_ATOM_ATOM_NAME+1] == 'H') continue;
	    
	/* solvent? ion? */
	if (  ! strncmp (line+PDB_ATOM_RES_NAME, "SOL",  PDB_ATOM_RES_NAME_LEN)  ) {

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

