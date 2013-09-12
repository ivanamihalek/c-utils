#include <stdio.h>
# include "sqc.h"

double avg_no_nbrs ( int no_seqs, int ** neighbors ) {
    int seq1, seq2, no_nbrs;
    double avg = 0.0;
    for ( seq1=0; seq1 < no_seqs; seq1++ ) {
	no_nbrs = 0;
	for ( seq2=0; seq2 < no_seqs; seq2++ ) {
	    if ( seq2== seq1) continue;
	    no_nbrs += neighbors[seq1][seq2];
	}
	avg += no_nbrs;
    }

    return avg/no_seqs;
}
