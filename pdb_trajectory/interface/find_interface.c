# include "ift.h"

/* find interface on chain 0, as imprinited by chain 1 */

int find_interface (Protein * protein, Atom * ligand, int ligand_number,
		    double cutoff_dist, double cutoff_dist_to_ligand, 
		    int *interface, int if_label) {
    
    int resctr1, resctr2;
    int atomctr1, atomctr2, ligand_ctr;
    Atom * atomptr1, * atomptr2, *ligandptr;
    double dist,  aux;
    int no_res_1 = protein->length;
    int no_res_2 = (protein+1)->length;
    int contact_found;
    double cutoff_dist_sq = cutoff_dist*cutoff_dist;
    double cutoff_dist_to_ligand_sq = cutoff_dist_to_ligand*cutoff_dist_to_ligand;
    Residue * sequence1 = protein->sequence;
    Residue * sequence2 = (protein+1)->sequence;

    if ( if_label ) {
	 sequence1 = (protein+1)->sequence;
	 sequence2 = (protein)->sequence;
	 no_res_1  = (protein+1)->length;
	 no_res_2  = (protein)->length;
    }
    
    for ( resctr1=0; resctr1<no_res_1; resctr1++) {
	
	interface[resctr1] = 0;

	contact_found = 0;
	for ( ligand_ctr=0; ligand_ctr<ligand_number && !contact_found; ligand_ctr++) {
	    ligandptr = ligand+ligand_ctr;
	    for (atomctr1=0; atomctr1 < sequence1[resctr1].no_atoms && !contact_found; atomctr1++) {
		atomptr1 = sequence1[resctr1].atom + atomctr1;
		
		aux = atomptr1->x - ligandptr->x;
		dist = aux*aux;
		if (dist > cutoff_dist_to_ligand_sq) continue;
		
		aux = atomptr1->y - ligandptr->y;
		dist += aux*aux;
		if (dist > cutoff_dist_to_ligand_sq) continue;
		
		aux = atomptr1->z - ligandptr->z;
		dist += aux*aux;
		if (dist > cutoff_dist_to_ligand_sq) continue;
		
		contact_found = 1;
		
	    }
	}
	if (contact_found) break; /* if close to ligand, the atom is not counted as a part of the interface*/
	
	contact_found = 0;
	for ( resctr2=0; resctr2<no_res_2 && !contact_found; resctr2++) {

	    for (atomctr1=0; atomctr1 < sequence1[resctr1].no_atoms && !contact_found; atomctr1++) {
		atomptr1 = sequence1[resctr1].atom + atomctr1;
		//if (atomptr1-> backbone) continue;
    
		/* find  dist btw atoms*/
		for (atomctr2=0; atomctr2 < sequence2[resctr2].no_atoms && !contact_found; atomctr2++) {
		    atomptr2 = sequence2[resctr2].atom + atomctr2;
		    //if (atomptr2-> backbone) continue;
		    
		    aux = atomptr1->x - atomptr2->x;
		    dist = aux*aux;
		    if (dist > cutoff_dist_sq) continue;
		    
		    aux = atomptr1->y - atomptr2->y;
		    dist += aux*aux;
		    if (dist > cutoff_dist_sq) continue;
		    
		    aux = atomptr1->z - atomptr2->z;
		    dist += aux*aux;
		    if (dist > cutoff_dist_sq) continue;
		 
		    interface[resctr1] = 1;
		    contact_found = 1;
		    
		}
	    }	    
	}
	    
    }

     return 0;
}


