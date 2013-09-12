# include "postp.h"


/* res_cvg array maps each reaidue into cvg it belongs to */
/* int_cvg is an array of all coverages (expressed as integers)  implied by the score array  */

int coverage ( Protein * protein, int * almt2prot, double * score, int almt_length, int * res_rank, int * int_cvg ) {

    double * protein_score, prev_score; /* "score" refers to alignment positions */
    int * sorted_res;
    int pos, ctr, ctr2;
    int first, cvg_ctr;
    
    /*allocate */
    protein_score = (double *) emalloc ( protein->length*sizeof (double) );
    if (!protein_score ) return 1;
    sorted_res    =    (int *) emalloc ( protein->length*sizeof (int) );
    if (!sorted_res )     return 1;
   
    /* remove gapped positions from the score array */
    pos = 0;
    for (ctr=0; ctr < almt_length; ctr ++ ) {
	if ( almt2prot[ctr] >= 0 ) {
	    protein_score[pos] = score[ctr];
	    pos ++;
	}
    }
    
    /* sort protein residues according to the new array */
    for (pos=0; pos < protein->length; pos++) sorted_res[pos] = pos;
    array_qsort ( sorted_res, protein_score, protein->length);

    /* turn the sorted array to coverage info */
    prev_score = protein_score[ sorted_res[0] ];
    first   = 0;
    cvg_ctr = 0;
    int_cvg[cvg_ctr] = 1;
    for (ctr=1; ctr < protein->length; ctr++) {
	if ( protein_score[ sorted_res[ctr] ] <= prev_score ) {
	    int_cvg[cvg_ctr] ++;
	} else {
	    prev_score  = protein_score[ sorted_res[ctr] ];
	    for (ctr2=first; ctr2 <ctr; ctr2++ ) {
		res_rank[ sorted_res[ctr2] ] = int_cvg[cvg_ctr];
	    }
	    first = ctr;
	    cvg_ctr ++;
	    int_cvg[cvg_ctr] =  int_cvg[cvg_ctr-1] + 1;
	}
    }
    for (ctr2=first; ctr2 <ctr; ctr2++ ) {
	res_rank[ sorted_res[ctr2] ] =  int_cvg[cvg_ctr];
    }
    


    /* free */
    free (protein_score);
    free (sorted_res);

    
    return 0;
}
