# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "sqc.h"

int main ( int argc, char * argv[]) {
    
    char distfname[100] = {'\0'};
    double cutoff_dist = 0.0;
    Sequence * sequence = NULL;
    int seq_ctr, no_seqs, seq_len, pos_ctr;
    double ** distmat;
    int ** cluster;
    int seq1, seq2;
    int c;
    FILE * fclust = NULL;
    
    if ( argc < 3 ) {
	fprintf (stderr, "Usage: %s <dist_file> <cutoff_sim> \n", argv[0]);
    }
    sprintf ( distfname, "%s", argv[1]);
    cutoff_dist = 1.0 - atof ( argv[2]);

    /* input seqs */
    if ( read_distf ( distfname, &no_seqs, &seq_len ) ) exit(1);
    printf ("There are %d seqs of length %d in %s.\n",  no_seqs, seq_len, distfname);


    /* allocate space for dist matrix */
    distmat = (double**) dmatrix (0,  no_seqs-1, 0,  no_seqs-1);
    /* calculate dist */
    if ( find_dist_matrix(sequence, no_seqs, seq_len, distmat) ) exit(1);

    /* cluster counting ... */
    {
	int  no_of_clusters, max_size, secnd_max_size;
	int * cluster_count, *mask;
	int ** neighbors;
	    
	cluster_count       =  (int *) emalloc ( (no_seqs+1)*sizeof(int));
	mask                =  (int *) emalloc ( (no_seqs+1)*sizeof(int));
	cluster             =  imatrix (0, no_seqs, 0, no_seqs);
	neighbors           =  imatrix (0, no_seqs-1, 0, no_seqs-1);
	
	for (seq1=0;  seq1 < no_seqs; seq1++ ) {
	    neighbors [seq1][seq1] = 1;
	    for (seq2= seq1+1; seq2 < no_seqs; seq2++ ) {
		neighbors[seq1][seq2] = ( distmat[seq1][seq2] < cutoff_dist);
		neighbors[seq2][seq1] = neighbors[seq1][seq2];
	    }
	}

	for (seq1=0;  seq1 < no_seqs; seq1++ ) {
	    mask[seq1] = 1;
	}
	cluster_counter (no_seqs,  neighbors,  mask, cluster_count, & no_of_clusters,
			 &max_size, &secnd_max_size , cluster);
    }

    /* output */
    fclust = stdout;
    for ( c=0; c <= no_seqs; c++) {
	if ( ! cluster[c][0] ) {
	    continue;
	}
	if ( !c ) {
	    fprintf ( fclust,"\t isolated:\n");
	} else {
	    fprintf ( fclust,"\t cluster size: %3d \n", cluster[c][0]); 
	}
	for ( seq_ctr=1; seq_ctr <=  cluster[c][0]; seq_ctr++) {
	    fprintf ( fclust, "%s  \n", sequence[ cluster[c][seq_ctr] ].name );
	}
	
    }
   
    return 0;

    
    
}
