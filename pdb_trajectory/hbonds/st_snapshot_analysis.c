# include "st_hbonds.h" 

# define HBOND_DIST 3.9

int  st_snapshot (Protein * protein, Atom  * solvent, int no_solvent_atoms, int bond[][4], int bond_count[], int *  no_bonds_ptr){

    Atom * atomptr1, * atomptr2;
    Residue * sequence = protein->sequence;
    Residue * other_sequence = (protein+1)->sequence;
    int no_res = protein->length;
    int other_no_res = (protein+1)->length;
    int atomctr1, atomctr2, resctr1,  resctr2;
    int no_bonds, bond_ctr, bond_found;
    double aux, dist;
    double x, y, z;
    double hb_dist_sq = HBOND_DIST*HBOND_DIST;

    
    no_bonds = *no_bonds_ptr;
    
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	
	for (atomctr1=0; atomctr1 < sequence[resctr1].no_atoms; atomctr1++) {

	    atomptr1 = sequence[resctr1].atom + atomctr1;
	    
	    if ( ! atomptr1->acceptor && ! atomptr1->donor) continue;

	    
	    x=atomptr1->x;
	    y=atomptr1->y;
	    z=atomptr1->z;

	    /* contacts with dimerizing partner */
	    for ( resctr2=0; resctr2<other_no_res; resctr2++) {
		for (atomctr2=0; atomctr2 < other_sequence[resctr2].no_atoms; atomctr2++) {
		    atomptr2 =  other_sequence[resctr2].atom + atomctr2;
		    if ( ! atomptr2->acceptor && ! atomptr2->donor) continue;
		    
		    aux = x - atomptr2->x;
		    dist = aux*aux;
		    aux = y - atomptr2->y;
		    dist += aux*aux;
		    aux = z - atomptr2->z;
		    dist += aux*aux;
		    if ( dist <= hb_dist_sq) {

			/* have wee seen this bond already?*/
			bond_found = 0;
			for (bond_ctr=0; bond_ctr < no_bonds; bond_ctr++) {
			    if ( bond[bond_ctr][0] == resctr1  &&  bond[bond_ctr][1] == resctr2 &&
				 bond[bond_ctr][2] == atomctr1 &&  bond[bond_ctr][3] == atomctr2     ) {
				bond_found = 1;
				bond_count[bond_ctr] ++;
				break;
			    }
			}
			if ( ! bond_found) {
			    if ( no_bonds == MAX_HBONDS) {
				fprintf (stderr, "Too many H-bonds. Increase MAX_HBONDS and recompile.\n");
				exit (1);
			    }
			    bond[no_bonds][0] = resctr1;
			    bond[no_bonds][1] = resctr2;
			    bond[no_bonds][2] = atomctr1;
			    bond[no_bonds][3] = atomctr2;
			    bond_count[bond_ctr] = 1;
			    no_bonds ++;
			}
		    }
		}
	    }


	    /* contacts with partner through a water intermediate */
	    
	}
	
    }
    
    *no_bonds_ptr = no_bonds;
    
    return 0;
}
