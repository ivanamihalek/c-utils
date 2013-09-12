# include "postp.h"


int  scoring ( Options *options, Alignment * alignment, double * score) {
    int  entropy ( Alignment * alignment, double *score);
    int  hybrid (Alignment * alignment, double *score);
    int  int_trace (Alignment * alignment, double *score);
    switch  ( options->scoring_method) {
    case ENTROPY:
	entropy ( alignment, score);
	break;
    case IVET:
	int_trace ( alignment, score);
	break;
    case RVET:
	hybrid ( alignment, score);
	break;
    }
    /* sink the gapped positions to the bottom, if requested */
    if ( options->max_gaps ) {
	int ctr;
	double max_score = -1;
	for ( ctr=0; ctr < alignment->length; ctr++ ) {
	    if ( max_score < score[ctr] ) max_score = score[ctr];
	}
	for ( ctr=0; ctr < alignment->length; ctr++ ) {
	    if ((double)alignment->column_gaps[ctr]/alignment->number_of_seqs> options->max_gaps ) {
		score[ctr] = max_score;
	    }
	}
    }

    return 0;
}
