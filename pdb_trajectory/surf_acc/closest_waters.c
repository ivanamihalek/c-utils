# include "sac.h" 

int  closest_waters (Protein * protein, Atom *solvent, int solvent_number,
		     double *dist_to_water, double nbr_dist, double * avg_nr_of_nbrs,
		     Neighbor **nbr_list, int * curr_list_length, int max_nbr_list_length){

    Atom * atomptr1, * atomptr2;
    Residue * sequence = protein->sequence;
    int no_res = protein->length;
    int atomctr1, atomctr2, resctr;
    int nbr, nbr_found;
    double aux, dist, min_dist;
    double x, y, z;
    double nbr_dist_sq = nbr_dist*nbr_dist;
    
    for ( resctr=0; resctr<no_res; resctr++) {

	min_dist = 1.e8;
	
	for (atomctr1=0; atomctr1 <solvent_number ; atomctr1++) {
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
		    
		    avg_nr_of_nbrs[resctr] ++;
		    
		    /* have we seen this neighbor already? */
		    nbr_found = 0;
		    for (nbr=0; nbr <curr_list_length[resctr]; nbr++ ) {
			if ( nbr_list[resctr][nbr].id == atomptr1->pdb_id) {
			    nbr_found = 1;
			    break;
			}
		    }
		    if (  nbr_found ) {  /*yes*/
			nbr_list[resctr][nbr].count ++;
		    } else {            /* no */
			if ( curr_list_length[resctr] == max_nbr_list_length ) {
			    fprintf ( stderr, "Too many neighbors in closest_waters(). Increase MAX_NBRS and recompile.\n" );
			    exit (1);
			}
			nbr_list[resctr][ curr_list_length[resctr] ].id = atomptr1->pdb_id;
			nbr_list[resctr][nbr].count = 1;
			curr_list_length[resctr] ++;
		    }
		    break; /* <==== we count each molecule as contacting residue once */
		}
	    }
	}
	//if ( min_dist > 1.e8-1 ) {
	// fprintf ( stderr, "Error in closest_waters(): residue %s.\n",
	//	      sequence[resctr].pdb_id);
	//    exit (1);
	//}
	min_dist = sqrt ( min_dist);
	dist_to_water[resctr] += min_dist;
    }


    return 0;
}
