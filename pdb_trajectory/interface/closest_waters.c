# include "ift.h" 

int  closest_waters (Protein * protein, int *interface, Atom *solvent, int solvent_number, double *dist_to_water){

    Atom * atomptr1, * atomptr2;
    Residue * sequence = protein->sequence;
    int no_res = protein->length;
    int atomctr1, atomctr2, resctr;
    double aux, dist, min_dist;
    double x, y, z;
    
    for ( resctr=0; resctr<no_res; resctr++) {
        if ( ! interface[resctr] ) continue;

	min_dist = 1.e8;
	
	for (atomctr1=0; atomctr1 < sequence[resctr].no_atoms; atomctr1++) {
	    atomptr1 = sequence[resctr].atom + atomctr1;
	    x=atomptr1->x;
	    y=atomptr1->y;
	    z=atomptr1->z;
	    
	    for (atomctr2=0; atomctr2 < solvent_number; atomctr2++) {
		atomptr2 = solvent + atomctr2;
		aux = x - atomptr2->x;
		dist = aux*aux;
		aux = y - atomptr2->y;
		dist += aux*aux;
		aux = z - atomptr2->z;
		dist += aux*aux;
		if ( dist < min_dist) {
		    min_dist = dist;
		}
	    }
	}
	if ( min_dist > 1.e8-1 ) {
	    fprintf ( stderr, "Error in closest_waters(): residue %s.\n",
		      sequence[resctr].pdb_id);
	    exit (1);
	}
	min_dist = sqrt ( min_dist);
	dist_to_water[resctr] += min_dist;
    }


    return 0;
}
