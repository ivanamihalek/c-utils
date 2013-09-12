# include "postp.h"


int score2rank ( double * score, int * rank_order, int length ) {
    int pos, sort_pos, ro;
    int * sorted_res;
    double old_score;

    if ( ! (sorted_res = (int *) emalloc ( length*sizeof (int) ) ) ) return 1;
   
    /* sort protein residues according to the new array */
    for (pos=0; pos < length; pos++) sorted_res[pos] = pos;
    array_qsort (sorted_res, score, length);

    old_score = -1;
    ro = 0;
    for (sort_pos=0; sort_pos < length; sort_pos++) {
	pos =  sorted_res[sort_pos];
	if ( score [pos] != old_score ) {
	    old_score =  score [pos];
	    ro++;
	}
	rank_order [pos] = ro;
    }

    free (sorted_res);
    return 0;
}
