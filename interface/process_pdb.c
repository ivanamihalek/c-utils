
#include "ifc.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

int process_pdb (Pdb_map * pdb_map, int map_length, int chain_start[], Chain_info ** chain_info_ptr, int * no_of_chains_ptr) {

    char this_residue[10];
    int no_of_chains, c, atom, residue;
    Chain_info * chain_info;
    
    
    /* how many chains ? */
    no_of_chains = 0;
    c = 0;
    while ( chain_start[ c ] >= 0 ) {
	no_of_chains++;
	c++;
    }
# ifdef VERBOSE
    printf (" there are %d chains.\n", no_of_chains);
#endif
    /*allocate space for chain info */
    if ( ! (chain_info = (Chain_info*) calloc ( no_of_chains, sizeof (Chain_info)) ) ) {
	fprintf ( stderr, "Error allocating chain info space.\n" );
	return 1;
    }
    
    /* how many residues in each chain? */
    for (c=0; c<no_of_chains;c++) {
	chain_info[c].start = chain_start[ c ];
	chain_info[c].end =  (c==(no_of_chains-1)) ?
	    map_length-1 :  chain_start[ c+1 ] - 1;
	chain_info[c].id = pdb_map[ chain_start[c] ].chain_id;
	residue = 1;
	atom = chain_info[c].start  ;
	memset ( this_residue, 0, 10);
	strcpy (this_residue, pdb_map[atom].no_res);
	while ( atom <= chain_info[c].end ) {
	    if ( strcmp (pdb_map[atom].no_res, this_residue) ) {
		residue++;
		memset ( this_residue, 0, 10);
		strcpy (this_residue, pdb_map[atom].no_res);

	    }
	    atom++;
	}
	chain_info[c].number_of_residues = residue;
	/* allocate residue position arrays */
	if ( !  (chain_info[c].residue_start = calloc ( residue, sizeof(int))) ) {
	    fprintf ( stderr, "Error allocating chain info space.\n" );
	    return 1;
	}
	if ( !  (chain_info[c].residue_end = calloc ( residue, sizeof(int))) ) {
	    fprintf ( stderr, "Error allocating chain info space.\n" );
	    return 1;
	}

        /* one more run through the whole thing to initialise res positions: */
	residue = 0;
	atom = chain_info[c].start;
	chain_info[c].residue_start[residue] = atom;
	memset ( this_residue, 0, 10);
	strcpy (this_residue, pdb_map[atom].no_res);
	while ( atom <=  chain_info[c].end) {
	    if ( strcmp (pdb_map[atom].no_res, this_residue) ) {
		chain_info[c].residue_end[residue] = atom-1;
		residue++;
		chain_info[c].residue_start[residue] = atom;
		memset ( this_residue, 0, 10);
		strcpy (this_residue, pdb_map[atom].no_res);

	    }
	    atom++;
	}
	chain_info[c].residue_end[residue] = atom-1;;
   }
# if 0
    for (c=0; c<no_of_chains;c++) {
	printf ("Chain %2d, id %c\n", c,  chain_info[c].id );
	printf ("\t               start: %8d\n",  chain_info[c].start);
	printf ("\t                 end: %8d\n",  chain_info[c].end);
 	printf ("\t  number of residues: %8d\n",  chain_info[c].number_of_residues);
	for ( residue=0; residue<chain_info[c].number_of_residues; residue++) {
	    printf ("\t\t\t %4d %8d %8d \n", residue,
		    chain_info[c].residue_start[residue],  chain_info[c].residue_end[residue]);
	    
	} 
	printf ("\n");
    }
# endif    
    

    * no_of_chains_ptr = no_of_chains;
    * chain_info_ptr = chain_info;


   
    return 0;

    
}
