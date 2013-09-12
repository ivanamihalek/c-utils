# include "st_hbonds.h"

# define MAX_LIGAND_TYPES 10
# define FILENAME_LENGTH  50

int main ( int argc, char * argv[]) {
    
    Protein protein[2];
    Atom * solvent = NULL;
    Atom * ligand;
    Neighbor *** nbr_list;
    int ** curr_list_length;
    int max_nbr_list_length;
    int chain_number, solvent_number, done, first;
    int ctr, ctr2, frame_ctr, resctr;
    int i, j;
    int do_nbr_stats;
    int invert;
    int ligand_number, no_ligand_types;
    int bond[MAX_HBONDS][4] = {{0}}, no_bonds = 0;
    int bond_count[MAX_HBONDS] = {0};
    int resctr1, resctr2;
    int atomctr1, atomctr2;
    char filename[FILENAME_LENGTH] = {'\0'};
    char ligand_type[MAX_LIGAND_TYPES][3] = {{'\0'}};
    double ** dist_to_water = {NULL};
    double ** avg_nr_of_nbrs = {NULL};
    double nbr_dist = 3.7;
    FILE * fptr;
   
    int read_pdb (char *pdbname, Protein *protein, int *chain_number,
		   Atom ** solvent, int *solvent_number,   Atom ** ligand, int *ligand_number,
		  char ligand_type [][3], int no_ligand_types, int invert, int *done);
    int  set_acc_and_donors (Protein * protein,  int  chain_number);
    int  st_snapshot (Protein * protein, Atom  * solvent, int no_solvent_atoms, int bond[][4], int bond_count[], int *  no_bonds_ptr);
   
    if ( argc < 2 ) {
	fprintf (stderr, "Usage: %s <pdb file name> [-n] [-l <lig type 1> <lig type 2> ...] \n", argv[0]);
	exit (1);
    }

    /* options */
    do_nbr_stats = 0;
    invert = 0;
    no_ligand_types =0;
    for (ctr = 2; ctr < argc; ctr++ ) {
	if ( ! strncmp (argv[ctr], "-n", 2) ) {
	    do_nbr_stats = 1;
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
    printf ( "\n");

    
    while ( !done ) {
	if ( first ) {
	    printf ("frame: %4d\n", frame_ctr+1);
	} else {
	    printf ("frame: %4d\r", frame_ctr+1);
	    fflush (stdout);
	}
	if (  read_pdb (argv[1], protein, &chain_number, &solvent, &solvent_number, &ligand, &ligand_number,
			ligand_type, no_ligand_types, invert,  &done) ){
	    fprintf (stderr, "Error reading %s.\n", argv[1] );
	    exit (1);
	}
	if (done) break;

	if ( first ) {
	    set_acc_and_donors (protein,  chain_number);
	}
	
	first = 0;

	st_snapshot (protein, solvent, solvent_number, bond, bond_count,  &no_bonds );
	frame_ctr++;

/* 	done = (frame_ctr==10); */
    }
    printf ( "\n");

    if ( done && ! (protein->length + (protein+1)->length) ) {
	fprintf ( stderr, "Error reading pdb: empty file?\n");
	exit (1);
    }

    sprintf ( filename, "hbonds.table", i+1);
    if ( ! (fptr = efopen ( filename, "w")) ) exit (1);

    for (ctr=0; ctr< no_bonds; ctr++ ) {
	resctr1 = bond[ctr][0];
	resctr2 = bond[ctr][1]; 
	atomctr1= bond[ctr][2];
	atomctr2= bond[ctr][3];
	fprintf (fptr, "\t  %4d %4s %3s   %4d %4s %3s %7.2lf\n",
		resctr1+1, protein->sequence[resctr1].res_type,
		protein->sequence[resctr1].atom[atomctr1].type,
		resctr2+1, (protein+1)->sequence[resctr2].res_type,
		(protein+1)->sequence[resctr2].atom[atomctr2].type,
		(double)bond_count[ctr]/frame_ctr);
    }
    fclose (fptr);
    


    
    return 0;
}
