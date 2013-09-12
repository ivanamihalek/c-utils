
#include "ifc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


# define CUTOFF_DIST 5.0
int binding_potential (Pdb_map *pdb_map, int atom, int bp []);

int  find_interface (Pdb_map * pdb_map, int map_length, Chain_info *chain_info, int no_of_chains, double cutoff){

    
    int c1, c2; /* chain counters*/
    int atom1, atom2, done, atom;
    double x, y, z, dist;
    int residue1, residue2;
    int max_res;
    char * amino[] = {"GLY", "ALA", "VAL", "LEU","ILE",
                      "MET", "PRO", "TRP",  "PHE", "SER",
                      "CYS", "THR", "ASN",  "GLN", "TYR",
		      "LYS", "ARG",   "HIS",  "ASP", "GLU" };
    int ifc_size[2], h_donor[2], h_acceptor[2], positive[2], negative[2], cysteine[2];
    int * interface_1, *interface_2, i;
    double sign, scale;
    int  bbone_bbone, bbone_other, other;

    if ( !cutoff ) cutoff = CUTOFF_DIST;
    
 

    /* what is max number of residues? */
    max_res = -1;
    for (c1=0; c1<no_of_chains;c1++) {
	if( max_res < chain_info[c1].number_of_residues) {
	    max_res = chain_info[c1].number_of_residues;
	}
    }
    interface_1 = (int*) calloc ( max_res, sizeof(int) ); 
    interface_2 = (int*) calloc ( max_res, sizeof(int) ); 
    if ( !interface_1 || !interface_2 ) {
	fprintf ( stderr, "Error allocating interface array  space.\n" );
	return 1;
     }
    bbone_bbone = bbone_other = other = 0;
    for (c1=0; c1<no_of_chains-1;c1++) {
	for (c2=c1+1; c2< no_of_chains; c2++) {
	    
	    ifc_size[0] = 0;
	    ifc_size[1] = 0;
	    memset (interface_1, 0, max_res*sizeof(int));
	    memset (interface_2, 0, max_res*sizeof(int));
	    h_donor[0] = h_donor[1] = 0;
	    h_acceptor[0] = h_acceptor[1] = 0;
	    positive[0] = negative[1] = 0;
	    positive[1] = negative[0] = 0;
	    cysteine[1] = cysteine[0] = 0;
	    
	    for ( residue1=0; residue1<chain_info[c1].number_of_residues; residue1++) {
		for ( residue2=0; residue2<chain_info[c2].number_of_residues; residue2++) {
		    if ( interface_1[residue1] && interface_2[residue2] ) {
			continue;
		    }
		    done = 0;
		    for (atom1=chain_info[c1].residue_start[residue1];
			 atom1<=chain_info[c1].residue_end[residue1] && !done ; atom1++) {
			//if ( strlen(pdb_map[atom1].name_atom) == 1)  continue;
			//if ( ! strchr(pdb_map[atom1].name_atom,'N')  &&   ! strchr(pdb_map[atom1].name_atom,'O')) continue;
			
			for (atom2=chain_info[c2].residue_start[residue2];
			     atom2<=chain_info[c2].residue_end[residue2] && !done ; atom2++) {
			    //if ( strlen(pdb_map[atom2].name_atom) == 1)  continue;
			    //if ( ! strchr(pdb_map[atom2].name_atom,'N')  &&   ! strchr(pdb_map[atom2].name_atom,'O')) continue;
			    

			    x = pdb_map[atom1].x- pdb_map[atom2].x;
			    y = pdb_map[atom1].y- pdb_map[atom2].y;
			    z = pdb_map[atom1].z- pdb_map[atom2].z;
			    dist = sqrt (x*x+y*y+z*z);
			    if ( dist <= cutoff ) {
# if 0
				printf ("res1: %s   atom1: %3s   nbrs: %3d  dist_to_cm: %8.3f   --- %8.3fA ---   res2: %s   atom2: %3s   nbrs: %3d   dist_to_cm: %8.3f   \n",
					pdb_map[atom1].name_res, pdb_map[atom1].name_atom, pdb_map[atom1].no_of_neighbors,  pdb_map[atom1].dist_to_cm,
					dist,
					pdb_map[atom2].name_res, pdb_map[atom2].name_atom , pdb_map[atom2].no_of_neighbors,  pdb_map[atom2].dist_to_cm);
				memset (bp, 0, NO_BINDING_TYPES*sizeof(int) );
				binding_potential ( pdb_map, atom1, bp );

				printf ("binding potential for atom1: \n");
				for ( b = 0; b <  NO_BINDING_TYPES; b++ ) {
				    if (bp[b]) {
					printf ( "\t type %d x %d \n", b, bp[b] );
				    }
				}
				memset (bp, 0, NO_BINDING_TYPES*sizeof(int) );
				binding_potential ( pdb_map, atom2, bp );
				
				printf ("binding potential for atom2: \n");
				for ( b = 0; b <  NO_BINDING_TYPES; b++ ) {
				    if (bp[b]) {
					printf ( "\t type %d x %d \n", b, bp[b] );
				    }
				}
				printf ("\n");
# endif
				if ( ! interface_1[residue1] ) {
				    ifc_size[0]++;
				}
				if ( ! interface_2[residue2] ) {
				    ifc_size[1]++;
				}
				interface_1[residue1] = 1;
				interface_2[residue2] = 1;
			
					    
			    }/* end  ( dist < cutoff)  case */
			} /* end atom2 loop*/
		    }/* end atom1 loop*/
		} /* end residue2 loop */
	    } /* end residue1 loop */

	    if ( ifc_size[0] > 1 &&  ifc_size[1] > 1 ) {

		int bp[NO_BINDING_TYPES], bp1[NO_BINDING_TYPES], bp2[NO_BINDING_TYPES], b;
		int * bp_larger, *bp_smaller;
		int size_larger, size_smaller;
		double moments1[3] = {0.0}, moments2[3] = {0.0};
		double *mom_smaller, *mom_larger;

		
		/* binding potentials for each chain */
		/* now pretend I do not know who interacts with whom and  find binding potentials */
		memset (bp1, 0, NO_BINDING_TYPES*sizeof(int) );
		for ( residue1=0; residue1<chain_info[c1].number_of_residues; residue1++) {
		    if ( ! interface_1[residue1]) continue;
		    for (atom1=chain_info[c1].residue_start[residue1];
			 atom1<=chain_info[c1].residue_end[residue1] && !done ; atom1++) {
			binding_potential (pdb_map, atom1, bp1);
		    }
		}

		memset (bp2, 0, NO_BINDING_TYPES*sizeof(int) );
		for ( residue2=0; residue2<chain_info[c2].number_of_residues; residue2++) {
		    if ( ! interface_2[residue2]) continue;
		    for (atom2=chain_info[c2].residue_start[residue2];
			 atom2<=chain_info[c2].residue_end[residue2] && !done ; atom2++) {
			binding_potential (pdb_map, atom2, bp2);
		    }
		}
		
		/* moments of inertia */
		moments_of_inertia (pdb_map,chain_info, c1, interface_1, moments1);
		moments_of_inertia (pdb_map,chain_info, c2, interface_2, moments2);
		
		if ( ifc_size[0] >= ifc_size[1] ) {
		    bp_larger = bp1;
		    bp_smaller = bp2;
		    scale = ifc_size[0]/ifc_size[1];
		    size_larger = ifc_size[0];
		    size_smaller = ifc_size[1];
		    mom_larger  = moments1;
		    mom_smaller = moments2;
		} else {
		    bp_larger = bp2;
		    bp_smaller = bp1;
		    scale = ifc_size[1]/ifc_size[1];
		    size_larger = ifc_size[1];
		    size_smaller = ifc_size[0];
		    mom_larger  = moments2;
		    mom_smaller = moments1;
		}
# if 0
		printf (" mi ratios for the chain %d: ", c1);
		for (i=1; i < 3; i++){
		    printf ("%8.4lf ", fabs(moments1[i]/moments1[0]) );
		}
		printf ("\n");
		printf (" mi ratios for the chain %d: ", c2);
		for (i=1; i < 3; i++){
		    printf ("%8.4lf ",  fabs(moments2[i]/moments2[0]) );
		}
		printf ("\n");
		
		printf ( "%2d  %2d  \n",  ifc_size[0], ifc_size[0] );
		for ( b = 1; b <  NO_BINDING_TYPES; b++ ) {
		    if ( b==3) continue;
		    printf ( "    %4d  %4d  %4d  %4d  %8.3lf\n", b+1, bp_larger[b], bp_smaller[NO_BINDING_TYPES-b],
			     bp_larger[b] - bp_smaller[NO_BINDING_TYPES-b], (double)( bp_larger[b] - bp_smaller[NO_BINDING_TYPES-b])/size_larger );
		}
		printf ("\n\n");
# else
		printf ( "%d ",  size_larger);
		for ( b = 1; b <  NO_BINDING_TYPES; b++ ) {
		    if ( b==3) continue;
		    printf ( "    %4d",  bp_larger[b] );
		}
		for (i=1; i < 3; i++){
		    printf ("   %8.2lf ",  100*fabs(mom_larger[i]/mom_larger[0]) );
		}
		printf ( "\n");
		printf ( "%d ",  size_smaller);
		for ( b = 1; b <  NO_BINDING_TYPES; b++ ) {
		    if ( b==3) continue;
		    printf ( "    %4d",  bp_smaller[NO_BINDING_TYPES-b] );
		}
		for (i=1; i < 3; i++){
		    printf ("   %8.2lf ", 100* fabs(mom_smaller[i]/mom_smaller[0]) );
		}
		printf ( "\n\n");
# endif
	    }
	    
	} /* end chain2 loop */
	

    }/* end chain1 loop */

    
    return 0;
}
