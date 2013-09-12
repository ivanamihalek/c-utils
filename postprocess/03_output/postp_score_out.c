# include "postp.h"
int output_score ( char * base_filename, Protein * protein, Alignment * alignment,
		   double * score, int * almt2prot, int *res_rank) {

    int almt_pos;
    FILE * fptr;
    char filename[BUFFLEN], pdbid[PDB_ATOM_RES_NO_LEN];
    double cvg;
    char aa;

    sprintf (filename, "%s.score", base_filename);

    fptr = efopen (filename, "w");
    if ( !fptr) return 1;
    
    fprintf (fptr, "%%%6s%6s%6s%8s%8s%8s\n", "almt", "pdb ", "aa ", "score", "cvg ", "gaps");
    for (almt_pos = 0; almt_pos < alignment->length; almt_pos++) {
	if ( almt2prot && almt2prot[almt_pos] >= 0 ) {
	    sprintf (pdbid, "%s", protein->sequence[ almt2prot[almt_pos] ].pdb_id );
	    aa = protein->sequence[ almt2prot[almt_pos] ].res_type_short;
	    cvg = (double)res_rank [ almt2prot[almt_pos] ]/protein->length;
	    
	} else {
	    sprintf (pdbid, "%6s", "-");
	    aa = '.';
	    cvg = 1;
	}
	fprintf (fptr, "%6d%6s%6c%8.2lf%8.3lf%8.3lf\n",
		 almt_pos+1, pdbid, aa, score[almt_pos], cvg,
		 (double)alignment->column_gaps[almt_pos]/alignment->number_of_seqs);
    }

    fclose (fptr);

    return 0;
}
