# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "sqc.h"


int find_dist_matrix(Sequence *sequence, int no_seqs, int seq_len, double ** distmat, int sim) {
    
    int seq1, seq2;
    double pairwise_id (char *sequence1, char *sequence2,  int seq_len);
    double pairwise_sim (char *sequence1, char *sequence2,  int seq_len);
    double (*pairwise) (char *sequence1, char *sequence2,  int seq_len);

    if ( sim ) {
	pairwise = pairwise_sim;
    } else {
	pairwise = pairwise_id;
    }

    if ( !distmat) {
	fprintf (stderr, "In find_dist_matrix(): space for the dist matrix must be allocated.\n");
	return 1;
    }
    for (seq1=0;  seq1 < no_seqs; seq1++ ) {
	distmat[seq1][seq1] = 0;
	for (seq2= seq1+1; seq2 < no_seqs; seq2++ ) {
	    distmat[seq1][seq2] = 1 - pairwise( sequence[seq1].position,  sequence[seq2].position, seq_len);
	    distmat[seq2][seq1] = distmat[seq1][seq2];
	} 
    } 
    return 0;
}

double pairwise_sim (char *sequence1, char *sequence2,  int seq_len){
    
    int pos_ctr;
    static char * similarto = NULL;
    int effective_length;
    double similarity;
    void sim_init () {
	int i;
	char aa;
	char *amino_acid_order = "ABCDEFGHIKLMNPQRSTVWXYZ";
	int aao_strlen = strlen(amino_acid_order);

	similarto = (char*) emalloc ( aao_strlen*sizeof(char) );
	for(i=0;i<aao_strlen;i++){
	    aa = amino_acid_order[i];
	    similarto[(int)aa] = aa;
	}
# if 0
	similarto['I'] = 'V';
	similarto['S'] = 'T';
	similarto['D'] = 'E';
	similarto['K'] = 'R';
	similarto['Q'] = 'N';
	similarto['.'] = '.';
# endif
        similarto['A'] = 'V';
        similarto['I'] = 'V';
        similarto['L'] = 'V';
        similarto['M'] = 'V';
        similarto['C'] = 'V';

        similarto['Y'] = 'F';
        similarto['W'] = 'F';
        similarto['H'] = 'F';

        similarto['G'] = 'P';

        similarto['T'] = 'N';
        similarto['S'] = 'N';
        similarto['Q'] = 'N';

        similarto['K'] = 'R';

        similarto['D'] = 'E';


        similarto['.'] = '.';

	
    }

    if ( !similarto) {
	sim_init();
    }

    effective_length = 0;
    similarity = 0.0;
    for (pos_ctr=0; pos_ctr<seq_len; pos_ctr++) {
	if ( sequence1[pos_ctr] != '.' && sequence2[pos_ctr] != '.') {
	    effective_length++;
	    if ( similarto[(int)sequence1[pos_ctr]] == similarto[(int)sequence2[pos_ctr]] ) {
		similarity ++;
	    }
	}
    }
    similarity  /= effective_length;
    
    return similarity;
}

double pairwise_id (char *sequence1, char *sequence2,  int seq_len){
    
    int pos_ctr;
    int effective_length;
    double similarity;

    effective_length = 0;
    similarity = 0.0;
    for (pos_ctr=0; pos_ctr<seq_len; pos_ctr++) {
	if ( sequence1[pos_ctr] != '.' && sequence2[pos_ctr] != '.') {
	    effective_length++;
	    if ( sequence1[pos_ctr]  ==  sequence2[pos_ctr]  ) {
		similarity ++;
	    }
	}
    }
    similarity  /= effective_length;
    
    return similarity;
}

