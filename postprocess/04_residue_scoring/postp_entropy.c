# include "postp.h"

int entropy ( Alignment * alignment, double *score){

    int col, seq, ctr;
    int freq[ASCII]; /* ASCII == 128,  ASCII size  */
    double p, entropy;
    
    for (col = 0; col < alignment->length; col++){
	/* find frequencies */
	memset (freq, 0, ASCII*sizeof(int));
	for (seq=0; seq< alignment->number_of_seqs; seq++){
	    freq [ (int) alignment->sequence[seq][col] ] ++;
	}
	/* find entropy */
	entropy = 0.0;
	for ( ctr=0; ctr < 128; ctr++) {
	    if ( freq[ctr] ) {
		p = (double)freq[ctr]/alignment->number_of_seqs;
		entropy -= p*log(p);
	    }
	}
	score[col] = entropy;
    }

    return 0;
}
