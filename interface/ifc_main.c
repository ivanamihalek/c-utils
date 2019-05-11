#include "ifc.h"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>


int  read_pdb_file(char filename[], Pdb_map **pdb_map, int *map_length_ptr, int chain_start[]);
int  find_interface (Pdb_map * pdb_map, int map_length, int chain_start[], double cutoff);

int main ( int argc, char * argv[]) {

   
    Pdb_map * pdb_map;
    int chain_start[MAX_CHAINS], chain_ctr, map_length;
    double cutoff = 0.0;

    if (argc < 2 ) {
        fprintf (stderr, "Usage: %s <pdb_file_name> [<cutoff_distance>].\n", argv[0]);
        exit(1);
    }

    if (argc == 3)  cutoff = atof(argv[2]);

    if (read_pdb_file(argv[1], &pdb_map,  &map_length, chain_start) ) {
        fprintf (stderr, "Error reading pdbfile.\n");
        return 1;
    }

    
    printf ("\n chain starts: \n");
    chain_ctr = 0;
    while ( chain_start[chain_ctr]>= 0 ) {
        //printf (" %c: %d ", pdb_map[ chain_start[chain_ctr]].chain_id, chain_start[chain_ctr]);
        //printf (" %d %c\n",  chain_start[chain_ctr], pdb_map[ chain_start[chain_ctr]].chain_id);
        printf (" %d \n",  chain_start[chain_ctr]);
        chain_ctr++;
    }
    printf ("\n\n");
    exit(1);
    
    if (find_interface (pdb_map, map_length, chain_start, cutoff) ) {
        fprintf (stderr, "Error finding interface.\n");
        return 1;
    }			
    
    return 0; 
}
