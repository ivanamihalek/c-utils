# include "ift.h"

# define MAX_LIGAND_TYPES 10

int main ( int argc, char * argv[]) {
    
    Protein protein[2];
    Atom * solvent = NULL, * solvent_init = NULL;
    Atom * ligand;
    int chain_number, solvent_number, done, first;
    int ctr, ctr2, frame_ctr, resctr, resctr2, i, ok;
    int * interface[2] = {NULL}, *sorted_pos = NULL, pos = 0;
    int * if_solvent = NULL;
    int ** nbr_table = NULL;
    int *** contact_type = NULL;
    int do_nbr_stats;
    int invert, total_contacts;
    int ligand_number, no_ligand_types, solctr;
    char ligand_type[MAX_LIGAND_TYPES][3] = {{'\0'}};
    double cutoff_dist_to_ligand   = 0.0;
    double cutoff_dist   = 7.0;
    double cutoff_radius = 4.0;
    double cutoff_water  = 3.0;
    double if_solvent_cutoff = 3.5;
    double * solvent_min_dist[2] = {NULL};
    double **avg_nbr_dist = NULL;
    double *avg_distance_from_init  = NULL;   
    double * dist_to_water = NULL, *nearest_neighbor = NULL, *avg_density = NULL;
    double * solvent_neighbors = NULL;
   
    int read_pdb (char *pdbname, Protein *protein, int *chain_number,
		   Atom ** solvent, int *solvent_number,   Atom ** ligand, int *ligand_number,
		  char ligand_type [][3], int no_ligand_types, int invert, int *done);
    int find_interface (Protein * protein, Atom * ligand, int ligand_number,
			double cutoff_dist, double cutoff_dist_to_ligand, 
			int *interface, int if_label);
    int make_if_neighbor_table ( Protein * protein, int *interface, double cutoff,
				 int ** nbr_table, int *** contact_type);
    int neighbor_stats (Protein * protein, int *interface, int ** nbr_table,
			double ** avg_nbr_dist, double * nearest_neighbor,
			double *avg_density, double cutoff_radius);
    int water_clusters ( Atom* solvent, int solvent_number, int *if_solvent, double cutoff_dist, double * neighbors);
   
    if ( argc < 2 ) {
	fprintf (stderr, "Usage: %s <pdb file name> [-n] [-i] [-l <lig type 1> <lig type 2> ...] \n", argv[0]);
	exit (1);
    }

    /* options */
    do_nbr_stats = 0;
    invert = 0;
    no_ligand_types =0;
    for (ctr = 2; ctr < argc; ctr++ ) {
	if ( ! strncmp (argv[ctr], "-n", 2) ) {
	    do_nbr_stats = 1;
	} else if (! strncmp (argv[ctr], "-i", 2) ) {
	    invert = 1;
	} else if (! strncmp (argv[ctr], "-l", 2) ) {
	    for (ctr2 = ctr+1; ctr2 < argc; ctr2++ ) {
		if ( argv[ctr2][0] == '-' ) break;
		if ( no_ligand_types >= MAX_LIGAND_TYPES ) {
		    fprintf (stderr, "Too many ligand types. Increase MAX_LIGAND_TYPES and recompile.\n");
		    exit (1);
		}
		sprintf ( ligand_type[no_ligand_types],  "%3s", argv[ctr2]);
		no_ligand_types++;
	    }
	    
	}
    }

    done = 0;
    first = 1;
    frame_ctr = 0;
    //printf ( "\n");
    while ( !done ) {
# if 0
	if ( first ) {
	    printf ("frame: %4d\n", frame_ctr+1);
	} else {
	    printf ("frame: %4d\r", frame_ctr+1);
	    fflush (stdout);
	}
# endif
	if (  read_pdb (argv[1], protein, &chain_number, &solvent, &solvent_number, &ligand, &ligand_number,
			ligand_type, no_ligand_types, invert,  &done) ){
	    fprintf (stderr, "Error reading %s.\n", argv[1] );
	    exit (1);
	}
	if (done) break;
	if ( chain_number != 2 ) {
	    fprintf (stderr, "Error: currently enabled to work with exactly 2 chains (%d chain(s) present in %s)\n",
		     chain_number,argv[1] );
	    exit (1);
	}
	if ( first ) {

	    for (i=0; i<2; i++) {
		int if_ctr, if_size;
		if ( ! ( interface[i] = emalloc((protein+i)->length*sizeof(int)) )) exit(1);
		find_interface ( protein, ligand, ligand_number, cutoff_dist, cutoff_dist_to_ligand, interface[i], i);
		
		/* sanity: how many residues at the interface? */
		if_size = 0;
		for (if_ctr=0; if_ctr<protein->length; if_ctr++) {
		    if_size += interface[i][if_ctr];
		}
		if ( ! if_size ) {
		    fprintf (stderr, "Error: no residues at the interface %d.\n", i+1);
		    exit (1);
		}
		if ( if_size < 5 ) {
		    fprintf ( stderr, "Warning: only %d residues at the interface %d.\n", if_size, i+1);
		}
	    }

	    /* construct neighbor table for the interface */
	    if ( ! ( nbr_table = intmatrix (protein->length, (protein+1)->length) )) exit(1);
	    if ( ! (contact_type = (int ***) emalloc(protein->length*sizeof(int**) )) ) exit(1);
	    
	    for ( ctr = 0; ctr < protein->length; ctr++ ) {
		if ( ! (contact_type[ctr] = (int **) emalloc((protein+1)->length*sizeof(int*) )) ) exit(1);
		for (ctr2=0; ctr2 < (protein+1)->length; ctr2++) {
		    if ( ! (contact_type[ctr][ctr2] = (int *) emalloc(( (SC2SC+1))*sizeof(int) )) ) exit(1);
		}
	    }
	    make_if_neighbor_table (protein, interface[0], cutoff_radius,  nbr_table, contact_type);
	    /* space for different stats for the neighbors */

	    if ( ! ( nearest_neighbor = emalloc(protein->length*sizeof(double)) )) exit(1);
	    if ( ! ( avg_density = emalloc(protein->length*sizeof(double)) )) exit(1);
	    if ( ! ( avg_nbr_dist = dmatrix (protein->length, (protein+1)->length) )) exit(1);
	    
	}
	
	first = 0;
	
	neighbor_stats (protein, interface[0], nbr_table, avg_nbr_dist,  nearest_neighbor,
			    avg_density, cutoff_radius);
	
	frame_ctr++;

	//done = (frame_ctr==5);
    }
    //printf ( "\n");
    if ( done && ! protein->length ) {
	fprintf ( stderr, "Error reading pdb: empty file?\n");
	exit (1);
    }

    /* average  */
    if ( ! ( sorted_pos = emalloc(protein->length*sizeof(int)) )) exit(1);
    for ( resctr=0; resctr<protein->length; resctr++) {
	sorted_pos[resctr] = resctr;
        if ( ! interface[0][resctr] ) continue;
	
	avg_density[resctr]      /= frame_ctr;
	//avg_density[resctr]      /= cutoff_radius*cutoff_radius*cutoff_radius*4.0/3.0*M_PI;
	nearest_neighbor[resctr] /= frame_ctr;
	for ( resctr2=0; resctr2< (protein+1)->length; resctr2++) {
	    if ( ! nbr_table[resctr][resctr2] ) continue;
	    avg_nbr_dist[resctr][resctr2] /= frame_ctr;
	}
    }

    for ( resctr=0; resctr<protein->length; resctr++) {

        if ( ! interface[0][resctr] ) continue;

	pos = resctr;
	
	printf ( "  %5s  %c   %8.2lf   %8.2lf   ", protein->sequence[pos].pdb_id,
		 protein->sequence[pos].res_type_short, 	nearest_neighbor[resctr], avg_density[resctr]);
	
	if (0) {
	    for ( resctr2=0; resctr2< (protein+1)->length; resctr2++) {
	    
		if ( ! nbr_table[pos][resctr2] ) continue;
	    
		total_contacts = 0;
		for (ctr=0; ctr < SC2SC; ctr++ ) {
		    total_contacts += contact_type[pos][resctr2][ctr];
		}
		if ( total_contacts < 2 )  continue;
	    
		printf ( " | %c", (protein+1)->sequence[resctr2].res_type_short);
	    
		if ( contact_type[pos][resctr2][BB2BB] ) {
		    printf ( " bb%d", contact_type[pos][resctr2][BB2BB]);
		}
		if ( contact_type[pos][resctr2][BB2SC] ) {
		    printf (" bs%d", contact_type[pos][resctr2][BB2SC]);
		}
		if ( contact_type[pos][resctr2][SC2BB] ) {
		    printf (" sb%d", contact_type[pos][resctr2][SC2BB]);
		}
		if ( contact_type[pos][resctr2][SC2SC] ) {
		    printf (" ss%d", contact_type[pos][resctr2][SC2SC]);
		}

		printf ( " %6.2lf   ",  avg_nbr_dist[pos][resctr2]);
	    }
		
	    printf ( " |");
	}

	printf ( "\n");

    }


    
    return 0;
}
