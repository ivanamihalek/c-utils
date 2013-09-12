# include "utils.h"
#include "ifc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>




int binding_potential (Pdb_map *pdb_map, int atom, int bp []) {

   
    
    if  ( charcount(pdb_map[atom].name_atom,PDB_ATOM_ATOM_NAME_LEN) ==1) { /*this is a backbone atom*/
//	if (  pdb_map[atom].no_of_neighbors > MAX_BB_NBRS ) {
	if (  pdb_map[atom].dist_to_cm < 0.9 ) {
	} else {
	    if (strchr(pdb_map[atom].name_atom,'N') ) {
		bp [H_DONOR] += 1;
	    } else if (strchr(pdb_map[atom].name_atom,'O') ) {
		bp[H_ACCEPTOR] += 1;
	    } else if (strchr(pdb_map[atom].name_atom,'C') ) {
		bp[H_ACCEPTOR] += 1;
	    } else {
		fprintf ( stderr, "Design error: single-letter atom name not belonging to backbone: %s .\n",
			  pdb_map[atom].name_atom);
		exit (0);
	    }
	}
    } else {
//	if (  pdb_map[atom].no_of_neighbors > MAX_NBRS ) {
	if (  pdb_map[atom].dist_to_cm < 0.9 ) {
	} else {
	    if ( strstr(pdb_map[atom].name_atom,"C") ) {
		//bp [H_ACCEPTOR] += 1;
	    } else if ( !strncmp (pdb_map[atom].name_res, "SER", 3) ) {
		if ( strstr(pdb_map[atom].name_atom,"OG") ) {
		    bp [H_DONOR]    += 1;
		    bp [H_ACCEPTOR] += 1;
		} else if (strstr(pdb_map[atom].name_atom,"CB") ) {
		    bp [H_ACCEPTOR] += 1;
		} else {
		    printf ("atom %d  \"%s\", aa %s.\n",
			    atom, pdb_map[atom].name_atom, pdb_map[atom].name_res);  
		}
		
	    } else if  ( !strncmp (pdb_map[atom].name_res, "THR", 3) ) {
		if ( strstr(pdb_map[atom].name_atom,"OG") ) {
		    bp [H_DONOR]    += 1;
		    bp [H_ACCEPTOR] += 1; /*+= 2; */
		} else if (strstr(pdb_map[atom].name_atom,"C") ) {
		    bp [H_ACCEPTOR] += 1;
		} else {
		    printf ("atom %d  \"%s\", aa %s.\n",
			    atom, pdb_map[atom].name_atom, pdb_map[atom].name_res);  
		}
		
	    } else if  ( ! strncmp (pdb_map[atom].name_res, "TYR", 3) ) {
		if ( strstr(pdb_map[atom].name_atom,"OH") ) {
		    bp [H_DONOR]    += 1;
		    bp [H_ACCEPTOR] += 1; /*+= 2; */
		} else if (strstr(pdb_map[atom].name_atom,"C") ) {
		    bp [H_ACCEPTOR] += 1;
		} else {
		    printf ("atom %d  \"%s\", aa %s.\n",
			    atom, pdb_map[atom].name_atom, pdb_map[atom].name_res);  
		}
		
	    } else if  ( !strncmp (pdb_map[atom].name_res, "ASP", 3) || !strncmp (pdb_map[atom].name_res, "GLU", 3)) {
		if ( strstr(pdb_map[atom].name_atom,"O") ) {
		    bp [H_ACCEPTOR] += 1; /*+= 2; */
		    bp [NEGATIVE]   += 1;
		} else if (strstr(pdb_map[atom].name_atom,"C") ) {
		    bp [H_ACCEPTOR] += 1;
		} else {
		    printf ("atom %d  \"%s\", aa %s.\n",
			    atom, pdb_map[atom].name_atom, pdb_map[atom].name_res);  
		}
		
	    } else if  (  !strncmp (pdb_map[atom].name_res, "ASN", 3) || !strncmp (pdb_map[atom].name_res, "GLN", 3)) {
		if ( strstr(pdb_map[atom].name_atom,"O") ) {
		    bp [H_ACCEPTOR] += 1; /*+= 2; */
		} else if (strstr(pdb_map[atom].name_atom,"N") ) {
		    bp [H_DONOR] += 1; /*+= 2; */
		} else {
		    printf ("atom %d  \"%s\", aa %s.\n",
			    atom, pdb_map[atom].name_atom, pdb_map[atom].name_res);  
		}
		
		
	    } else if  ( !strncmp (pdb_map[atom].name_res, "ARG", 3) ) {
		if ( strstr(pdb_map[atom].name_atom,"NH") ) {
		    bp [H_DONOR]  += 1; /*+= 2; */
		    bp [POSITIVE] += 1;
		} else if (strstr(pdb_map[atom].name_atom,"NE") ) {
		    bp [H_DONOR] += 1;
		} else {
		    printf ("atom %d  \"%s\", aa %s.\n",
			    atom, pdb_map[atom].name_atom, pdb_map[atom].name_res);  
		}
		
	    } else if  ( !strncmp (pdb_map[atom].name_res, "LYS", 3) ) {
		if (strstr(pdb_map[atom].name_atom,"N") ) {
		    bp [H_DONOR]  += 1;/* += 3; */
		    bp [POSITIVE] += 1;
		} else {
		    printf ("atom %d  \"%s\", aa %s.\n",
			    atom, pdb_map[atom].name_atom, pdb_map[atom].name_res);  
		}
		
		
	    } else if  ( strncmp (pdb_map[atom].name_res, "TRP", 3) ) {
		if (strstr(pdb_map[atom].name_atom,"N") ) {
		    bp [H_DONOR]  += 1;
		}
	    } 
	}
	
    }


    return 0;
}
