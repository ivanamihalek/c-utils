
#include "pdb.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


# define CUTOFF_DIST 5.0

int  find_interface (Pdb_map * pdb_map, int map_length, int chain_start[], double cutoff){

    typedef struct {
	char id;
	int start;
	int end;
	int number_of_residues;
	int *residue_start;
	int *residue_end;
    } Chain_info;
    
    int c1, c2; /* chain counters*/
    int no_of_chains;
    int atom1, atom2, done;
    double x, y, z, dist;
    char this_residue[10];
    int residue1, residue2;
    int max_res, * interface_1, *interface_2;
    Chain_info * chain_info;
    char * amino[] = {"GLY", "ALA", "VAL", "LEU","ILE",
                      "MET", "PRO", "TRP",  "PHE", "SER",
                      "CYS", "THR", "ASN",  "GLN", "TYR",
		      "LYS", "ARG",   "HIS",  "ASP", "GLU" };
    int *ifc_pop_1 ,*ifc_pop_2, ifc_size_1, ifc_size_2, a;

    if ( !cutoff ) cutoff = CUTOFF_DIST;
    
 
    /* how many chains ? */
    no_of_chains = 0;
    c1 = 0;
    while ( chain_start[ c1 ] >= 0 ) {
	no_of_chains++;
	c1++;
    }
    printf (" there are %d chains.\n", no_of_chains);

    /*allocate space for chain info */
    if ( ! (chain_info = (Chain_info*) calloc ( no_of_chains, sizeof (Chain_info)) ) ) {
	fprintf ( stderr, "Error allocating chain info space.\n" );
	return 1;
    }
    
    /* how many residues in each chain? */
    for (c1=0; c1<no_of_chains;c1++) {
	chain_info[c1].start = chain_start[ c1 ];
	chain_info[c1].end =  (c1==(no_of_chains-1)) ?
	    map_length-1 :  chain_start[ c1+1 ] - 1;
	chain_info[c1].id = pdb_map[ chain_start[c1] ].chain_id;
	residue1 = 1;
	atom1 = chain_info[c1].start  ;
	memset ( this_residue, 0, 10);
	strcpy (this_residue, pdb_map[atom1].no_res);
	while ( atom1 <= chain_info[c1].end ) {
	    if ( strcmp (pdb_map[atom1].no_res, this_residue) ) {
		residue1++;
		memset ( this_residue, 0, 10);
		strcpy (this_residue, pdb_map[atom1].no_res);

	    }
	    atom1++;
	}
	chain_info[c1].number_of_residues = residue1;
	/* allocate residue position arrays */
	if ( !  (chain_info[c1].residue_start = calloc ( residue1, sizeof(int))) ) {
	    fprintf ( stderr, "Error allocating chain info space.\n" );
	    return 1;
	}
	if ( !  (chain_info[c1].residue_end = calloc ( residue1, sizeof(int))) ) {
	    fprintf ( stderr, "Error allocating chain info space.\n" );
	    return 1;
	}

        /* one more run through the whole thing to initialise res positions: */
	residue1 = 0;
	atom1 = chain_info[c1].start;
	chain_info[c1].residue_start[residue1] = atom1;
	memset ( this_residue, 0, 10);
	strcpy (this_residue, pdb_map[atom1].no_res);
	while ( atom1 <=  chain_info[c1].end) {
	    if ( strcmp (pdb_map[atom1].no_res, this_residue) ) {
		chain_info[c1].residue_end[residue1] = atom1-1;
		residue1++;
		chain_info[c1].residue_start[residue1] = atom1;
		memset ( this_residue, 0, 10);
		strcpy (this_residue, pdb_map[atom1].no_res);

	    }
	    atom1++;
	}
	chain_info[c1].residue_end[residue1] = atom1-1;;
   }

    for (c1=0; c1<no_of_chains;c1++) {
	printf ("Chain %2d, id %c\n", c1,  chain_info[c1].id );
	printf ("\t               start: %8d\n",  chain_info[c1].start);
	printf ("\t                 end: %8d\n",  chain_info[c1].end);
 	printf ("\t  number of residues: %8d\n",  chain_info[c1].number_of_residues);
# if 0
	for ( residue1=0; residue1<chain_info[c1].number_of_residues; residue1++) {
	    printf ("\t\t\t %4d %8d %8d \n", residue1,
		    chain_info[c1].residue_start[residue1],  chain_info[c1].residue_end[residue1]);
	    
	} 
# endif    
	printf ("\n");
    }
    

# if 1



    /* what is max number of residues? */
    max_res = -1;
    for (c1=0; c1<no_of_chains;c1++) {
	if( max_res < chain_info[c1].number_of_residues) {
	    max_res = chain_info[c1].number_of_residues;
	}
    }
    ifc_pop_1 = (int*) calloc ( 20, sizeof(int) ); 
    ifc_pop_2 = (int*) calloc ( 20, sizeof(int) ); 
    interface_1 = (int*) calloc ( max_res, sizeof(int) ); 
    interface_2 = (int*) calloc ( max_res, sizeof(int) ); 
    if ( !interface_1 || !interface_2 ) {
	fprintf ( stderr, "Error allocating interface array  space.\n" );
	return 1;
     }
    
    
    for (c1=0; c1<no_of_chains-1;c1++) {
	for (c2=c1+1; c2< no_of_chains; c2++) {
	    
	    memset (ifc_pop_1, 0, 20*sizeof(int));
	    memset (ifc_pop_2, 0, 20*sizeof(int));
	    memset (interface_1, 0, max_res*sizeof(int));
	    memset (interface_2, 0, max_res*sizeof(int));
	    ifc_size_1 = 0;
	    ifc_size_2 = 0;
	    for ( residue1=0; residue1<chain_info[c1].number_of_residues; residue1++) {
		for ( residue2=0; residue2<chain_info[c2].number_of_residues; residue2++) {
		    if ( interface_1[residue1] && interface_2[residue2] ) {
			continue;
		    }
		    done = 0;
		    for (atom1=chain_info[c1].residue_start[residue1];
			 atom1<=chain_info[c1].residue_end[residue1] && !done ; atom1++) {
			for (atom2=chain_info[c2].residue_start[residue2];
			     atom2<=chain_info[c2].residue_end[residue2] && !done ; atom2++) {
			
			
			    x = pdb_map[atom1].x- pdb_map[atom2].x;
			    y = pdb_map[atom1].y- pdb_map[atom2].y;
			    z = pdb_map[atom1].z- pdb_map[atom2].z;
			    dist = sqrt (x*x+y*y+z*z);
			    if ( dist <= cutoff ) {
				done = 1; /* done with this residue pair */ 
				/* mark residues as belonging to interface */
				interface_1[residue1] = 1;
				interface_2[residue2] = 1;
				ifc_size_1++;
				ifc_size_2++;
				for (a=0; a<20; a++) {
				    if ( ! strcmp (pdb_map[atom1].name_res, amino[a] )) {
					ifc_pop_1[a] ++;
					break;
				    }
				}
				for (a=0; a<20; a++) {
				    if ( ! strcmp (pdb_map[atom2].name_res, amino[a] ) ) {
					ifc_pop_2[a] ++;
					break;
				    }
				}
					    
			    }
			}
		    }
		}

	    }
	    if ( ifc_size_1 && ifc_size_2) {
		for (a=0; a<20; a++) {
		    printf ( " %5.3lf", (double)ifc_pop_1[a]/ifc_size_1);
		}
		printf ( "\n\t");
		for (a=0; a<20; a++) {
		    printf ( " %5.3lf", (double)ifc_pop_2[a]/ifc_size_2);
		}
		printf ( "\n\n");
	    }
	}
    }
# endif
    
    return 0;
}
