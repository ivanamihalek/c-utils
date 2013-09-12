# include "postp.h"


int read_score ( char *filename, Protein * protein, int * score2prot, double *score ) {
    
    FILE *fptr = NULL;
    char line[LONGSTRING];
    char pdb_id[PDB_ATOM_RES_NO+2];
    double score_val;
    int res_ctr, ctr;
    
   
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;
    memset ( line, 0, LONGSTRING);
    ctr = 0;
    while(fgets(line,LONGSTRING,fptr)!=NULL){
	sscanf (line, " %s %lf \n", pdb_id,  &score_val);
	for ( res_ctr=0; res_ctr < protein->length; res_ctr++ ) {
	    if ( ! strncmp( protein->sequence[res_ctr].pdb_id, pdb_id, PDB_ATOM_RES_NO_LEN+2) ) {
		score2prot[ctr] = res_ctr;
		score[ctr]      = score_val;
		ctr ++;
	    }
	}
    }
    fclose (fptr);
    return 0;
}
