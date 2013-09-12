# include "ift.h" 

int  neighbor_stats (Protein * protein, int *interface, int ** nbr_table,
		     double ** avg_nbr_dist, double * nearest_neighbor,
		      double *avg_density, double cutoff_radius){

    int resctr1, resctr2;
    int atomctr1, atomctr2;
    Atom * atomptr1, * atomptr2;
    double dist,  aux;
    int no_res_1 = protein->length;
    int no_res_2 = (protein+1)->length;
    double cutoff_radius_sq = cutoff_radius*cutoff_radius;
    double cutoff_too_far = cutoff_radius;
    double min_dist, min_dist_btw_atoms;
    int  number_of_neighbors;
    Residue * sequence1 = protein->sequence;
    Residue * sequence2 = (protein+1)->sequence;
    
    for ( resctr1=0; resctr1<no_res_1; resctr1++) {

	if ( ! interface[resctr1] ) continue;
	
	min_dist =  cutoff_too_far*cutoff_too_far;
	number_of_neighbors = 0;
	
	for ( resctr2=0; resctr2<no_res_2; resctr2++) {
	    if (! nbr_table[resctr1][resctr2] ) continue;

	    min_dist_btw_atoms =  cutoff_too_far*cutoff_too_far;
	    for (atomctr1=0; atomctr1 < sequence1[resctr1].no_atoms; atomctr1++) {
		atomptr1 = sequence1[resctr1].atom + atomctr1;
    
		/* find  dist btw atoms*/
		for (atomctr2=0; atomctr2 < sequence2[resctr2].no_atoms; atomctr2++) {
		    atomptr2 = sequence2[resctr2].atom + atomctr2;
		    
		    aux = atomptr1->x - atomptr2->x;
		    dist = aux*aux;
		    aux = atomptr1->y - atomptr2->y;
		    dist += aux*aux;
		    aux = atomptr1->z - atomptr2->z;
		    dist += aux*aux;
		    
		    if (dist <= cutoff_radius_sq) {
			number_of_neighbors ++;
			if (dist < min_dist) {
			    min_dist = dist;
			}
			if (dist < min_dist_btw_atoms) {
			    min_dist_btw_atoms = dist;
			}
		    }
		} /*atomctr2*/
	    }/*atomctr1*/
	    avg_nbr_dist[resctr1][resctr2] += sqrt( min_dist_btw_atoms);
	}
	    
	min_dist = sqrt ( min_dist);
	nearest_neighbor[resctr1] += min_dist;
	avg_density[resctr1]  += number_of_neighbors;
    }

    return 0;
}


/*****************************************************************************/
int make_if_neighbor_table ( Protein * protein, int *interface, double cutoff, int ** nbr_table, int *** contact_type) {

    int resctr1, resctr2;
    int atomctr1, atomctr2;
    Atom * atomptr1, * atomptr2;
    double dist,  aux;
    int no_res_1 = protein->length;
    int no_res_2 = (protein+1)->length;
    double cutoff_sq = cutoff*cutoff;
    double cutoff_too_far = 100.0;
    int too_far;
    Residue * sequence1 = protein->sequence;
    Residue * sequence2 = (protein+1)->sequence;
    

    for ( resctr1=0; resctr1<no_res_1; resctr1++) {
	
	if ( ! interface[resctr1] ) continue;
	
	for ( resctr2=0; resctr2<no_res_2; resctr2++) {

	    too_far =0;
	    nbr_table[resctr1][resctr2] = 0;
	    for (atomctr1=0; atomctr1 < sequence1[resctr1].no_atoms && !too_far; atomctr1++) {
		atomptr1 = sequence1[resctr1].atom + atomctr1;
    
		/* find  dist btw atoms*/
		for (atomctr2=0; atomctr2 < sequence2[resctr2].no_atoms && !too_far; atomctr2++) {
		    atomptr2 = sequence2[resctr2].atom + atomctr2;
		    
		    aux = atomptr1->x - atomptr2->x;
		    if ( aux > cutoff_too_far ) {
			too_far = 1;
			break;
		    }
		    dist = aux*aux;
		    aux = atomptr1->y - atomptr2->y;
		    if ( aux > cutoff_too_far ) {
			too_far = 1;
			break;
		    }
		    dist += aux*aux;
		    aux = atomptr1->z - atomptr2->z;
		    if ( aux > cutoff_too_far ) {
			too_far = 1;
			break;
		    }
		    dist += aux*aux;
		    
		    if (dist <= cutoff_sq ) {
			/* note: this is not assumed symmterica: protein 1 and protein2
			   do not even have to have the same number of residues */
			nbr_table[resctr1][resctr2] = 1;
			if ( atomptr1->backbone ) {
			    if ( atomptr2->backbone ) {
				contact_type[resctr1][resctr2][BB2BB] ++;
			    } else {
				contact_type[resctr1][resctr2][BB2SC] ++;
			    }
			} else {
			    if ( atomptr2->backbone ) {
				contact_type[resctr1][resctr2][SC2BB] ++;
			    } else {
				contact_type[resctr1][resctr2][SC2SC] ++;
			    }
			
			}
		    }
		} /*atomctr2*/
	    }/*atomctr1*/	    
	}
	    
    }

 
    return 0;
}
