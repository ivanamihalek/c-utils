# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "sqc.h"

int main ( int argc, char * argv[]) {
    
    char msfname[100] = {'\0'};
    char outname[20] = {'\0'};
    double cutoff_dist = 0.0;
    Sequence * sequence = NULL;
    int no_seqs, seq_len;
    double ** distmat;
    int seq1, seq2;
    FILE * fnames = NULL;
    
    if ( argc < 3 ) {
	fprintf (stderr, "Usage: %s <msf_file> <cutoff_sim> [<outname>] \n", argv[0]);
    }
    sprintf ( msfname, "%s", argv[1]);
    cutoff_dist = 1.0 - atof ( argv[2]);

    if ( argc >= 4 ) {
	sprintf (outname, "%s", argv[3]);
    } else {
	sprintf (outname, "%s", "outcont");
    }

    /* input seqs */
    if ( read_msf ( msfname, &sequence, &no_seqs, &seq_len ) ) exit(1);
    printf ("There are %d seqs of length %d in %s.\n",  no_seqs, seq_len, msfname);


    /* allocate space for dist matrix */
    distmat = (double**) dmatrix (0,  no_seqs-1, 0,  no_seqs-1);
    /* calculate dist */
    if ( find_dist_matrix(sequence, no_seqs, seq_len, distmat) ) exit(1);


    /* output */
    sprintf ( msfname, "%s.names", outname);
    if ( ! (fnames = fopen ( msfname,"w") ) ) exit (1);
    for (seq1= 0; seq1 < no_seqs; seq1++ ) {
	fprintf ( fnames, " %4d  %15s\n", seq1+1, sequence[seq1].name);
    }
    fclose (fnames);
    
    sprintf ( msfname, "%s.contacts", outname);
    if ( ! (fnames = fopen ( msfname,"w") ) ) exit (1);
    for (seq1= 0; seq1 < no_seqs; seq1++ ) {
	for (seq2= seq1+1; seq2 < no_seqs; seq2++ ) {
	    //if ( distmat[seq1][seq2] < cutoff_dist) {
		fprintf ( fnames, " %4d  %15s      %4d  %15s    %8.3lf \n",
			  seq1+1, sequence[seq1].name, seq2+1, sequence[seq2].name, distmat[seq1][seq2]);
		//}
	}
    }
    fclose (fnames);
    
    return 0;

    
    
}
