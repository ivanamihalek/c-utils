# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "sqc.h"

int main ( int argc, char * argv[]) {
    
    char msfname[100] = {'\0'};
    double cutoff_dist = 0.0;
    Sequence * sequence = NULL;
    int no_seqs, seq_len;
    double ** distmat;
    int ** cluster;
    int seq1, seq2;
    int  similarity;
    int  no_of_clusters, max_size, secnd_max_size;
    int * cluster_count, *mask;
    int ** neighbors;
    double avg_no_nbrs ( int no_seqs, int ** neighbors );
  
    if ( argc < 2 ) {
	fprintf (stderr, "Usage: %s <msf_file>  \n", argv[0]);
    }
    sprintf ( msfname, "%s", argv[1]);

    /* input seqs */
    if ( read_msf ( msfname, &sequence, &no_seqs, &seq_len ) ) exit(1);
    printf ("%% There are %d seqs of length %d in %s.\n",  no_seqs, seq_len, msfname);


    /* allocate space for dist matrix */
    distmat = (double**) dmatrix (0,  no_seqs-1, 0,  no_seqs-1);
    /* calculate dist */
    if ( find_dist_matrix(sequence, no_seqs, seq_len, distmat, similarity=0)) exit(1);

    /* cluster counting ... */
    printf ("%%  %10s   %10s  %10s\n", "cutoff", "no clust", "avg_nbrs" );
    for( cutoff_dist= 0.1; cutoff_dist <= 0.75; cutoff_dist += 0.05) {
	    
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
	
	printf ("  %10.2lf  %10d    %10.1lf \n",
		cutoff_dist, no_of_clusters, avg_no_nbrs(no_seqs, neighbors));

	
    }

   
    return 0;

    
    
}
