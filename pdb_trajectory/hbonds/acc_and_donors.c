# include "st_hbonds.h" 

int  set_acc_and_donors (Protein * protein,  int  chain_number) {

    int i, resctr, atomctr;
    int no_res;
    Atom * atomptr;
    Residue * sequence;

    for (i=0; i < chain_number; i ++) {
	sequence = (protein+i)->sequence;
	no_res = (protein+i)->length;
	
	for ( resctr=0; resctr<  no_res; resctr++) {
	    for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
		
		atomptr =  sequence[resctr].atom + atomctr;
		atomptr->acceptor = 0;
		atomptr->donor = 0;
		
		if (!  strcmp ( atomptr->type, "O") ) {
		    atomptr->acceptor = 1;
		} else if (!  strcmp ( atomptr->type, "N") ) {
		    atomptr->donor = 1;
		} else  {

		    switch ( sequence[resctr].res_type_short) {
		    case 'D':
			if (! strncmp ( atomptr->type, "OD", 2) ) {
			    atomptr->acceptor = 1;
			}
			break;
		    case 'E':
			if ( !strncmp ( atomptr->type, "OE", 2) ) {
			    atomptr->acceptor = 1;
			}
			break;
		    case 'H':
			if ( !strncmp ( atomptr->type, "ND3", 3) ) {
			    atomptr->acceptor = 1;
			} else if ( !strncmp ( atomptr->type, "NE2", 3) ) {
			    atomptr->donor = 1;
			} else if ( !strncmp ( atomptr->type, "ND1", 3) ) {
			    atomptr->donor = 1;
			}
			break;
		    case 'K':
			if ( !strncmp ( atomptr->type, "NZ", 2) ) {
			    atomptr->donor = 1;
			}
			break;
		    case 'M':
			if ( !strncmp ( atomptr->type, "SD", 2) ) {
			    atomptr->acceptor = 1;
			}
			break;
		    case 'N':
			if ( !strncmp ( atomptr->type, "OD1", 3) ) {
			    atomptr->acceptor = 1;
			} else if ( !strncmp ( atomptr->type, "ND2", 3) ) {
			    atomptr->donor = 1;
			}
			break;
		    case 'Q':
			if ( !strncmp ( atomptr->type, "OE1", 3) ) {
			    atomptr->acceptor = 1;
			} else if ( !strncmp ( atomptr->type, "NE2", 3) ) {
			    atomptr->donor = 1;
			}
			break;
		    case 'R':
			if ( !strncmp ( atomptr->type, "NE", 2) ) {
			    atomptr->donor = 1;
			} else if ( !strncmp ( atomptr->type, "NH1", 3) ) {
			    atomptr->donor = 1;
			}else if ( !strncmp ( atomptr->type, "NH2", 3) ) {
			    atomptr->donor = 1;
			}
			break;
		    case 'S':
			if ( !strncmp ( atomptr->type, "OG", 2) ) {
			    atomptr->acceptor = 1;
			    atomptr->donor = 1;
			}
			break;
		    case 'T':
			if ( !strncmp ( atomptr->type, "OG1", 3) ) {
			    atomptr->acceptor = 1;
			    atomptr->donor = 1;
			}
			break;
		    case 'Y':
			if ( !strncmp ( atomptr->type, "OH", 2) ) {
			    atomptr->acceptor = 1;
			    atomptr->donor = 1;
			}
			break;
		    case 'W':
			if ( !strncmp ( atomptr->type, "NE1", 3) ) {
			    atomptr->donor = 1;
			}
			break;
		    
		    }
		}
	    
	    }
	
	}

    }

    return 0;
}
