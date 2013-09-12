# include "cleanup_waters.h"
/* read in pdb file and sever chain id's (maybe),
   and output all waters within CUTOFF Angstroms from the chain(s),
   in gro format and some ad hoc protons attached*/

# define CUTOFF 2.4


int main ( int argc, char * argv[]) {

    Atom* solvent = NULL;
    Atom * ligand = NULL;
    Protein protein[2];
    int solvent_number, chain_number, ligand_number;
    int gro_output (Atom*  solvent, int solvent_number);
    int read_pdb ( char * pdbname, Protein * protein, int * chain_number,
	       Atom ** solvent_ptr,int * solvent_number,
	       Atom ** ligand_ptr, int * ligand_number) ;
    int    mark_clashing_waters (Protein *protein, int chain_number, Atom *ligand,
				 int ligand_number, Atom * water, int solvent_number,  double min_dist);
    
    if ( argc < 2 ) {
	fprintf (stderr, "%s removes waters clashing with the protein atom coordinates,\n", argv[0]);
	fprintf (stderr, "and outputs water (with some lame hydrogens)  in gro format\n");
	fprintf (stderr, "\nUsage: %s <pdb file> \n\n", argv[0]);
	exit (1);
    }

    /* read protein */
    if (  read_pdb (argv[1], protein, &chain_number, &solvent, &solvent_number,
		    &ligand, & ligand_number) ){
	fprintf (stderr, "Error reading %s.\n", argv[1] );
	exit (1);
    }


    /* mark waters which clash with the protein, the ligand, or among themselves*/
    mark_clashing_waters (protein, chain_number, ligand, ligand_number, solvent, solvent_number, CUTOFF);

    /* output the remaining waters in gro format*/
    gro_output ( solvent, solvent_number);
    return 0;
}



/****************************************************/
/****************************************************/
int    mark_clashing_waters (Protein *protein, int chain_number, Atom *ligand,
			     int ligand_number, Atom * solvent, int solvent_number, double min_dist) {

    Atom * atomptr1, * atomptr2;
    Residue * sequence;
    int i, atomctr1, atomctr2, resctr, no_res;
    double aux, dist, min_dist_sq = min_dist*min_dist;
    double x, y, z;
    
   /* clash with protein ? */
    for (i=0; i<chain_number; i++ ) {
	no_res = (protein+i)->length;
	sequence = (protein+i)->sequence;
	
	for ( resctr=0; resctr<no_res; resctr++) {
	
	    for (atomctr1=0; atomctr1 < sequence[resctr].no_atoms; atomctr1++) {
		atomptr1 = sequence[resctr].atom + atomctr1;
		x=atomptr1->x;
		y=atomptr1->y;
		z=atomptr1->z;
	    
		for (atomctr2=0; atomctr2 < solvent_number; atomctr2++) {
		    atomptr2 = solvent + atomctr2;
		    if ( atomptr2->clash ) continue;
		    aux = x - atomptr2->x;
		    dist = aux*aux;
		    aux = y - atomptr2->y;
		    dist += aux*aux;
		    aux = z - atomptr2->z;
		    dist += aux*aux;
		    atomptr2->clash =   ( dist < min_dist_sq );
		}
	    }
	
	}
    }
    /* clash with ligand ? */
    for (atomctr1=0; atomctr1 < ligand_number; atomctr1++) {
	atomptr1 = ligand + atomctr1;
	x=atomptr1->x;
	y=atomptr1->y;
	z=atomptr1->z;
	for (atomctr2=0; atomctr2 < solvent_number; atomctr2++) {
	    atomptr2 = solvent + atomctr2;
	    if ( atomptr2->clash ) continue;
	    aux = x - atomptr2->x;
	    dist = aux*aux;
	    aux = y - atomptr2->y;
	    dist += aux*aux;
	    aux = z - atomptr2->z;
	    dist += aux*aux;
	    atomptr2->clash =   ( dist < min_dist_sq );
	}
	
    }


    /* clash among themseleves? */
    for (atomctr1=0; atomctr1 < solvent_number; atomctr1++) {
	atomptr1 = solvent + atomctr1;
	x=atomptr1->x;
	y=atomptr1->y;
	z=atomptr1->z;
	for (atomctr2 = atomctr1+1; atomctr2 < solvent_number; atomctr2++) {
	    atomptr2 = solvent + atomctr2;
	    if ( atomptr2->clash ) continue;
	    aux = x - atomptr2->x;
	    dist = aux*aux;
	    aux = y - atomptr2->y;
	    dist += aux*aux;
	    aux = z - atomptr2->z;
	    dist += aux*aux;
	    atomptr2->clash =   ( dist < min_dist_sq );
	}
	
    }
  
    int clashing = 0;;
    for (atomctr1=0; atomctr1 < solvent_number; atomctr1++) {
	atomptr1 = solvent + atomctr1;
	clashing += atomptr1->clash;
    }
    //printf ("marked %d solvent atoms as clashing.\n", clashing);
    
    return 0;
}

/****************************************************/
/****************************************************/
int gro_output (Atom*  solvent, int solvent_number) {

    int ctr, atom_ctr, solvent_ctr;
    double x, y, z;

    ctr = 0;
    atom_ctr = 0;

    for (solvent_ctr =0; solvent_ctr < solvent_number; solvent_ctr++ ) {
	
	if ( (solvent+solvent_ctr)->clash ) continue;

	ctr ++;
	
	x = (solvent+solvent_ctr)->x;
	y = (solvent+solvent_ctr)->y;
	z = (solvent+solvent_ctr)->z;
	atom_ctr ++;
	printf ("%5d%5s%5s%5d%8.3lf%8.3lf%8.3lf\n",
		ctr, "SOL",  "OW", atom_ctr, x/10, y/10, z/10);

	z  = (solvent+solvent_ctr)->z + 1.0;
	atom_ctr ++;
	printf ( "%5d%5s%5s%5d%8.3lf%8.3lf%8.3lf\n",
	    ctr, "SOL",  "HW1", atom_ctr, x/10, y/10, z/10);

	x  = (solvent+solvent_ctr)->x + 0.948;
	z  = (solvent+solvent_ctr)->z - 0.33;
	atom_ctr ++;
	printf( "%5d%5s%5s%5d%8.3lf%8.3lf%8.3lf\n",
	    ctr, "SOL",  "HW2", atom_ctr, x/10, y/10, z/10);
	
	
    }
    return 0;
}
