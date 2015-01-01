# include "pdbclust.h"

# define NONSENSE  1000.0

int determine_dist_matrix ( double ** dist_matrix, Residue * sequence, int no_res){
    
    int resctr1, resctr2;
    int atomctr1, atomctr2;
    double dist, min_dist, aux;
    Atom * atomptr1, * atomptr2;

    /* fast return: */
    if ( !dist_matrix) {
	fprintf (stderr, "In determine_dist_matrix(): distance matrix space must be assigned.\n");
	return 1;
    }
    

    /*diagonal*/
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	dist_matrix [resctr1][resctr1] = 0.0;
    }

    for ( resctr1=0; resctr1<no_res; resctr1++) {
	
	for ( resctr2=resctr1+1; resctr2<no_res; resctr2++) {

	    min_dist = NONSENSE;
	    
	    for ( atomctr1=0; atomctr1 < sequence[resctr1].no_atoms; atomctr1 ++) {
		atomptr1 = sequence[resctr1].atom + atomctr1;

		for ( atomctr2=0; atomctr2 < sequence[resctr2].no_atoms; atomctr2++) {
		    atomptr2 = sequence[resctr2].atom + atomctr2;

		    
		    dist = 0.0;
		    aux = atomptr1->x -  atomptr2->x;
		    dist += aux*aux;
		    aux = atomptr1->y -  atomptr2->y;
		    dist += aux*aux;

		    aux = atomptr1->z -  atomptr2->z;
		    dist += aux*aux;

		    if ( dist < min_dist ) {
			min_dist = dist;
		    }
		}
		
	    }
	    dist_matrix [resctr2][resctr1] = dist_matrix [resctr1][resctr2] =  sqrt (min_dist );
	}
    }
#if 0
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	for ( resctr2=0; resctr2<no_res; resctr2++) {
	    printf ("%8.3lf ", dist_matrix [resctr1][resctr2]);
	}
	printf ("\n");
    }
    exit (0);
# endif    
    return 0;
}
