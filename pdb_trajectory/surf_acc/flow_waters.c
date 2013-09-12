# include "sac.h" 

int  flow_waters (Protein * protein, Atom *solvent, int solvent_number, int nbr_dist, 
		  Neighbor **old_nbr_list, int * old_list_length,
		  Neighbor **nbr_list, int * curr_list_length, int max_nbr_list_length, double * avg_id_change){

    Atom * atomptr1, * atomptr2;
    Residue * sequence = protein->sequence;
    int no_res = protein->length;
    int atomctr1, atomctr2, resctr;
    int nbr, nbr_found;
    int intersection;
    double aux, dist, min_dist;
    double x, y, z;
    double nbr_dist_sq = nbr_dist*nbr_dist;

    
    for ( resctr=0; resctr<no_res; resctr++) {

        intersection = 0;
	memset (nbr_list[resctr], 0, max_nbr_list_length*sizeof(Neighbor*) );
	curr_list_length[resctr] = 0;
	
	for (atomctr1=0; atomctr1 <solvent_number; atomctr1++) {
	    atomptr1 = solvent + atomctr1;
	    
	    x=atomptr1->x;
	    y=atomptr1->y;
	    z=atomptr1->z;
	    
	    for (atomctr2=0; atomctr2 <  sequence[resctr].no_atoms; atomctr2++) {
		atomptr2 = sequence[resctr].atom + atomctr2;

		if ( atomptr2->backbone) continue;
		
		aux = x - atomptr2->x;
		dist = aux*aux;
		aux = y - atomptr2->y;
		dist += aux*aux;
		aux = z - atomptr2->z;
		dist += aux*aux;
		if ( dist < min_dist) {
		    min_dist = dist;
		}
		if ( dist < nbr_dist_sq ) {
		    
		    
		    /* have we seen this neighbor already? */
		    nbr_found = 0;
		    for (nbr=0; nbr < old_list_length[resctr]; nbr++ ) {
			if ( old_nbr_list[resctr][nbr].id == atomctr1 ) {
			    intersection++;
			}
		    }

		    if ( curr_list_length[resctr] == max_nbr_list_length ) {
			fprintf ( stderr, "Too many neighbors in closest_waters(). Increase MAX_NBRS and recompile.\n" );
			exit (1);
		    }
		    nbr_list[resctr][ curr_list_length[resctr] ].id = atomctr1;
		    nbr_list[resctr][nbr].count = 1;
		    curr_list_length[resctr] ++;
		    
		    break; /* <==== we count each molecule as contacting residue once */
		}
	    }
	}
	//if ( min_dist > 1.e8-1 ) {
	// fprintf ( stderr, "Error in closest_waters(): residue %s.\n",
	//	      sequence[resctr].pdb_id);
	//    exit (1);
	//}
	avg_id_change[resctr] += old_list_length[resctr] + curr_list_length[resctr] - 2*intersection;
	
    }


    return 0;
}
