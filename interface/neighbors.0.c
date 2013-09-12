
#include "ifc.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* how many nieghbors for each atom within the chain ? how far is the atom from their CM?*/

int neighbors (Pdb_map * pdb_map,int map_length,  Chain_info * chain_info, int  no_of_chains) {

    int c, atom;
    int atom1, atom2;
    int residue1, residue2;
    double x, y, z, dist, cutoff;
    double avg, avg_bb;
    int no_bb;
    cutoff = 3.0;
     for (c=0; c<no_of_chains;c++) {
	 for ( residue1=0; residue1<chain_info[c].number_of_residues-1; residue1++) {
	     
	     for (atom1=chain_info[c].residue_start[residue1];
		  atom1<=chain_info[c].residue_end[residue1]; atom1++) {
		 
		 for ( residue2=residue1+1; residue2<chain_info[c].number_of_residues; residue2++) {
		     
		     for (atom2=chain_info[c].residue_start[residue2];
			  atom2<=chain_info[c].residue_end[residue2] ; atom2++) {

			    x = pdb_map[atom1].x- pdb_map[atom2].x;
			    y = pdb_map[atom1].y- pdb_map[atom2].y;
			    z = pdb_map[atom1].z- pdb_map[atom2].z;
			    dist = sqrt (x*x+y*y+z*z);
			    if ( dist <= cutoff ) {
				pdb_map[atom1].no_of_neighbors++;
				pdb_map[atom2].no_of_neighbors++;
			    }

			 
		     } /* end atom2 loop*/
			
		 }/* end residue2 loop*/
		 
	     }/* end atom1 loop*/
	 }/* end residue1 loop*/
	 
     }/* end chain loop*/


# if 1
     avg = avg_bb = 0;
     no_bb = 0;
     for (atom=0; atom< map_length; atom++ ) {
	 //printf (" %6d   %3d\n", atom, pdb_map[atom].no_of_neighbors);
	 avg += pdb_map[atom].no_of_neighbors;
	 if ( charcount(pdb_map[atom].name_atom,PDB_ATOM_ATOM_NAME_LEN) ==1 ){
	     avg_bb  += pdb_map[atom].no_of_neighbors;
	     no_bb++;
	 }
     }
     printf (" avg no nbrs:  %8.3f\n", avg/map_length);
     printf (" avg no nbrs for backbone:  %8.3f\n", avg_bb/no_bb);

# endif

     return 0;
}
