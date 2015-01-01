# include "pdbclust.h"

# define NONSENSE  1000.0

int  read_residue_selection (char * selected_res_file, Residue * sequence, int no_res, int *selected){
    
    char buf [LONGSTRING];
    char aux [PDB_ATOM_RES_NO_LEN+2];
    char res_type;
    int  found, total, res_ctr;
    FILE * fptr;
   
    fptr   = efopen (selected_res_file, "r" );
    if (! fptr ) {
	fprintf (stderr, "Error opening %s.\n",selected_res_file);
	return 1;
    }
    memset ( selected, 0, no_res);
    total = 0;
    while ( ! feof ( fptr) ) {
 	if ( fgets (buf, LONGSTRING, fptr ) && strlen(buf) > 1 ) {
	    if ( buf[0] == '#' ) continue;
	    if (! (sscanf (buf, "%s %c", aux, &res_type)  ) ) {
		fprintf (stderr, "Error reading %s.\n", selected_res_file);
		fclose (fptr);
		return 1;
	    }
	    found = 0;
	    for ( res_ctr=0; res_ctr < no_res; res_ctr++) {
		if ( ! strcmp (sequence[res_ctr].pdb_id, aux) ) {
		    if ( res_type !=  sequence[res_ctr].res_type_short ) {
			fprintf (stderr, "While reading %s: Type mismatch for res id %s (pdb:%c, here:%c) .\n",
				 selected_res_file, aux, sequence[res_ctr].res_type_short, res_type);
			fclose (fptr);
			return 1;
		    }
		    found = 1;
		    selected[res_ctr]= 1;
		    total++;
		}
	    }
	    if ( ! found ) {
		fprintf (stderr, "While reading %s: Residue with the id %s not found.\n",
			 selected_res_file, aux);
		fclose (fptr);
		return 1;
	    }
	}
    }
    
    if ( !total ) {
	fprintf (stderr, "No positions matching pdb found in %s.\n", selected_res_file);
	return 1;
    }

    fclose (fptr);

    return 0;
}
