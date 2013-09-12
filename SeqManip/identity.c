# include <stdlib.h>
# include <string.h>

# define MAX_SEQS   1000       /* max number of sequences that can be stored */
# define MAX_LENGTH  700       /* max length of each sequence                */

# define SUCCESS        0
# define MEM_ALLOC_ERR -1
# define GENERIC_ERR   -2

int identity (double percent_id, char * sequence, int seq_len, char *ok) {


    static char *storage [MAX_SEQS] = {NULL};
    static int no_of_seqs = 0;
    char *old_seq;
    int s, i;
    int diff, min_diff = (100.0 - percent_id)*0.01*seq_len;
    
    /*compare sequence with the ones already  stored*/
    for ( s=0; s < no_of_seqs; s++) {

	old_seq = storage[s];
	/* count the number of different sites */
	diff = 0;
	for ( i=0; i < seq_len; i++) {
	    if ( old_seq[i] != sequence[i] ) {
		diff++;
		if (diff > min_diff)
		    break;
	    }
	}
	/* check if we have found identity*/
	if ( diff <= min_diff ) {
	    /* the sequence is too similar to something we already have */
	    /* it is not ok to keep; return to caller*/
	    *ok = 'n';
	    return SUCCESS;
	}
	
    }
   
    /* if it's by more  than  perc ent_difference different  from everybody - store */
    /* create storage space*/
    storage [no_of_seqs] = (char *) calloc (MAX_LENGTH, sizeof(char));
    if ( !storage [no_of_seqs] ) {
	return MEM_ALLOC_ERR;
    }
    /* store */
    old_seq = memcpy ( storage [no_of_seqs], sequence, seq_len);
    if ( old_seq != storage [no_of_seqs]) {
	return GENERIC_ERR;
    }
    no_of_seqs++;

    *ok = 'y';
    return SUCCESS;

}


main() {
}
