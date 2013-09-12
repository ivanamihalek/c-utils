#include "ifc.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

#define MAX_CHAINS 20

int  read_pdb_file(char filename[], Pdb_map **pdb_map, int *map_length_ptr, int chain_start[]);
int  process_pdb (Pdb_map * pdb_map, int map_length, int chain_start[], Chain_info ** chain_info_ptr, int * no_of_chains_ptr);
int  find_interface (Pdb_map * pdb_map, int map_length, Chain_info *chain, int no_of_chains, double cutoff);
int neighbors (Pdb_map * pdb_map, int map_length, Chain_info * chain_info, int  no_of_chains);

int main ( int argc, char * argv[]) {

   
    Pdb_map * pdb_map;
    Chain_info * chain_info;
    
    int chain_start[MAX_CHAINS], chain_ctr, map_length, no_of_chains;
    double cutoff = 0.0;

    if (argc < 2 ) {
	fprintf (stderr, "Usage: %s <pdb_file_name> [<cutoff_distance>].\n", argv[0]);
	exit(1);
    }
    if ( argc == 3 ) {
	cutoff = atof( argv[2] );
    }
    for ( chain_ctr=0; chain_ctr<MAX_CHAINS; chain_ctr++) {
	chain_start[chain_ctr] = -1;
    }
    if (  read_pdb_file(argv[1], &pdb_map,  &map_length, chain_start) ) {
	fprintf (stderr, "Error reading pdbfile.\n");
	return 1;
    }

    if (  process_pdb  (pdb_map, map_length, chain_start, & chain_info, &no_of_chains )  ) {
	fprintf (stderr, "Error processing pdbfile.\n");
	return 1;
    }

    if ( neighbors(pdb_map, map_length, chain_info, no_of_chains) ) {
	fprintf (stderr, "Error processing neighbor info.\n");
	return 1;
    }			
    if (  find_interface (pdb_map, map_length, chain_info, no_of_chains, cutoff) ) {
	fprintf (stderr, "Error finding interface.\n");
	return 1;
    }			
    
    return 0; 
}
