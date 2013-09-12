# include "ift.h" 

int  if_waters (Protein * protein, int **interface, Atom *solvent, int solvent_number,
		int * if_solvent, double ** min_dist, double cutoff, Atom* solvent_init, double *avg_distance_from_init){

    Atom * atomptr;
    Residue * sequence[2];
    int no_res[2];
    int atomctr, resctr, solctr;
    int  i;
    double aux, dist, md;
    double x, y, z;

    for (i=0; i<2; i++) {
	no_res[i]   = (protein+i)->length;
	sequence[i] = (protein+i)->sequence;
    }
    
    for ( solctr=0; solctr<solvent_number; solctr++) {

	if_solvent[solctr] = 0;
	
	x= (solvent+solctr)->x;
	y= (solvent+solctr)->y;
	z= (solvent+solctr)->z;
	    

	
	for (i=0; i<2; i++) {

	    md = 400;
	    
	    for ( resctr=0; resctr<no_res[i]; resctr++) {
		if ( ! interface[i][resctr] ) continue;
	
		for (atomctr=0; atomctr < sequence[i][resctr].no_atoms; atomctr++) {
		    atomptr = sequence[i][resctr].atom + atomctr;
		    aux = x - atomptr->x;
		    dist = aux*aux;
		    aux = y - atomptr->y;
		    dist += aux*aux;
		    aux = z - atomptr->z;
		    dist += aux*aux;
		    if ( dist < md)  md=dist;
		}
	    }
	    min_dist[i][solctr] += sqrt(md);

	}

	if_solvent[solctr] = 1;


	
	/* distance from initial position */

	atomptr = solvent_init + solctr;
	
	aux = x - atomptr->x;
	dist = aux*aux;
	aux = y - atomptr->y;
	dist += aux*aux;
	aux = z - atomptr->z;
	dist += aux*aux;
	avg_distance_from_init[solctr] += sqrt (dist);
	
    }


    return 0;
}
