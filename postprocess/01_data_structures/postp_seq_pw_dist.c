# include "postp.h"

/* pairwise distances between sequences */

/*********************************************************************************************/
/*********************************************************************************************/

/*********************************************************************************************/
/*********************************************************************************************/

int seq_pw_dist (Alignment * alignment) {
    char *amino_acid_order = "ABCDEFGHIKLMNPQRSTVWXYZ";

    int blosum62[]={
	4,
	-2,  4,
	0, -3,  9,
	-2,  4, -3,  6,
	-1,  1, -4,  2,  5,
	-2, -3, -2, -3, -3,  6,
	0, -1, -3, -1, -2, -3,  6,
	-2,  0, -3, -1,  0, -1, -2,  8,
	-1, -3, -1, -3, -3,  0, -4, -3,  4,
	-1,  0, -3, -1,  1, -3, -2, -1, -3,  5,
	-1, -4, -1, -4, -3,  0, -4, -3,  2, -2,  4,
	-1, -3, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5,
	-2,  3, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6,
	-1, -2, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7,
	-1,  0, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,
	-1, -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5,
	1,  0, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,
	0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,
	0, -3, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4,
	-3, -4, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,
	0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1, -2, -1,
	-2, -3, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2, -1,  7,
	-1,  1, -3,  1,  4, -3, -2,  0, -3,  1, -3, -1,  0, -1,  3,  0,  0, -1, -2, -3, -1, -2,  4};
    
     int i,j,pos, ctr, off_diag;
    int aao_strlen;
    int number_of_similar_sites;
    /* int number_of_common_sites; */
    int similarity[128][128];
    char * seq_i, *seq_j;
    int char_i, char_j;
    double avg;

    if ( ! ( alignment->seq_dist = dmatrix ( alignment->number_of_seqs, alignment->number_of_seqs) ) ) return 1;
    

    /*get the similarity matrix to a usable form */
    aao_strlen = strlen(amino_acid_order);
    
    ctr = 0;
    avg = 0;
    off_diag = 0;
    for(i=0;i<aao_strlen;i++){
	char_i = (int) amino_acid_order [i];
	for (j=0;j<=i;j++){
	    char_j = (int) amino_acid_order [j];
	    similarity[char_i][char_j] = similarity[char_j][char_i] = blosum62[ctr];
	    if ( i != j  && blosum62[ctr] > 0) {
		avg +=  blosum62[ctr];
		off_diag ++;
	    }
	    ctr++;
	}
    }
    
    avg /= off_diag;
    avg = ceil(avg);
    
    int length_i, length_j, shorter;
    
    for (i=0; i<alignment->number_of_seqs; i++) {
	alignment->seq_dist[i][i] = 0.0;
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    alignment->seq_dist[i][j] = alignment->seq_dist[j][i] = 1.0;
	}
    }
    for (i=0; i<alignment->number_of_seqs; i++) {
	seq_i = alignment->sequence[i];
	length_i = alignment->length - alignment->seq_gaps[i];
	if ( ! length_i ) continue;
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    length_j = alignment->length - alignment->seq_gaps[j];
	    if ( ! length_j ) continue;
	    seq_j = alignment->sequence[j];
	    /* number_of_common_sites = 0; */
	    number_of_similar_sites = 0;
	    for (pos=0; pos<alignment->length; pos++) {
		if ( seq_i[pos] == '.' ) continue;
		if ( seq_j[pos] == '.' ) continue;
		/* number_of_common_sites ++; */
		if ( similarity [(int)seq_i[pos]] [(int) seq_j[pos] ] >= avg ) {
		    number_of_similar_sites ++;
		}
	    }
	    shorter = (length_i<length_j) ? length_i : length_j;
	     alignment->seq_dist[i][j] = 1 - (double) number_of_similar_sites/shorter;
	    /*  alignment->seq_dist[i][j] = 1 - (double) number_of_similar_sites/number_of_common_sites; */
	     alignment->seq_dist[j][i] =  alignment->seq_dist[i][j];
	}
    }
# if 0
    for (i=0; i<alignment->number_of_seqs; i++) {
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    printf (" %4d  %4d  %10s %10s  %8.3lf\n",
		    i+1, j+1, alignment->name[i],  alignment->name[j],  alignment->seq_dist[i][j]);
	}
    }
	exit (0);
# endif   
    return 0;
}
