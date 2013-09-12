
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
    double * cm_x, * cm_y, * cm_z, avg_dist_to_cm;
    double avg_sq, std_dev;

    cm_x = (double*) calloc ( map_length, sizeof(double));
    cm_y= (double*) calloc ( map_length, sizeof(double));
    cm_z = (double*) calloc ( map_length, sizeof(double));

    if ( !cm_x || ! cm_y || !cm_z ) {
	fprintf (stderr, "Error allocating cm arrays.\n");
	exit (1);
    }
    
    cutoff = 2.0;
	     
		     
    for (c=0; c<no_of_chains;c++) {
	for ( residue1=0; residue1<chain_info[c].number_of_residues-1; residue1++) {
	     
	    for (atom1=chain_info[c].residue_start[residue1];
		 atom1<=chain_info[c].residue_end[residue1]; atom1++) {
		 
		for ( residue2=residue1; residue2<chain_info[c].number_of_residues; residue2++) {
		     
		    for (atom2=chain_info[c].residue_start[residue2];
			 atom2<=chain_info[c].residue_end[residue2] ; atom2++) {

			if ( atom2==atom1) continue;
			
			x = pdb_map[atom1].x- pdb_map[atom2].x;
			y = pdb_map[atom1].y- pdb_map[atom2].y;
			z = pdb_map[atom1].z- pdb_map[atom2].z;
			dist = sqrt (x*x+y*y+z*z);
			if ( dist <= cutoff ) {
			    pdb_map[atom1].no_of_neighbors++;
			    pdb_map[atom2].no_of_neighbors++;
			    cm_x [atom1] += pdb_map[atom2].x;
			    cm_y [atom1] += pdb_map[atom2].y;
			    cm_z [atom1] += pdb_map[atom2].z;
			    cm_x [atom2] += pdb_map[atom1].x;
			    cm_y [atom2] += pdb_map[atom1].y;
			    cm_z [atom2] += pdb_map[atom1].z;
			}

		    } /* end atom2 loop*/
			
		}/* end residue2 loop*/
		 
	    }/* end atom1 loop*/
	}/* end residue1 loop*/
	 
    }/* end chain loop*/

			 
    avg_dist_to_cm  = 0.0;
    avg_sq  = 0.0;
    for (atom1=0; atom1 < map_length; atom1++) {
	if ( pdb_map[atom1].no_of_neighbors ) { 
	    cm_x [atom1]/= pdb_map[atom1].no_of_neighbors;
	    cm_y [atom1]/= pdb_map[atom1].no_of_neighbors;
	    cm_z [atom1]/= pdb_map[atom1].no_of_neighbors;
	    x = pdb_map[atom1].x - cm_x [atom1];
	    y = pdb_map[atom1].y - cm_y [atom1];
	    z = pdb_map[atom1].z - cm_z [atom1];
	    dist = sqrt (x*x+y*y+z*z);
	} else {
	    dist = 0;
	}
	pdb_map[atom1].dist_to_cm = dist;
	avg_dist_to_cm  += dist;
	avg_sq += dist*dist;
	//printf ("atom  %d  dist from cm of the neighbors: %8.3lf \n",
	//atom1, dist );
    }






    
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
     avg_dist_to_cm /= map_length;
     avg_sq /= map_length;
     std_dev = sqrt ( avg_sq -  avg_dist_to_cm*avg_dist_to_cm );

     for (atom1=0; atom1 < map_length; atom1++) {
	 pdb_map[atom1].dist_to_cm = 	(pdb_map[atom1].dist_to_cm - avg_dist_to_cm)/std_dev;
     }
    
# if 1
     printf (" avg no nbrs:  %8.3f\n", avg/map_length);
     printf (" avg no nbrs for backbone:  %8.3f\n", avg_bb/no_bb);
     printf (" avg dist_to_cm:  %8.3f\n", avg_dist_to_cm);
     printf (" stdev:  %8.3f\n", std_dev);

# endif

     free (cm_x);
     free (cm_y);
     free (cm_z);
     
     return 0;
}
