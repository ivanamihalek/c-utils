# include "postp.h"
int  struct_almt_mapping (Protein * protein, Alignment * alignment, char *query_name,  int * prot2almt, int * almt2prot){
    int prot_pos, almt_pos, query;
    char * query_seq;
    Residue * prot_seq;
    /* locate query in the alignment */
    query=0;
    while (query < alignment->number_of_seqs && strcmp (alignment->name[query], query_name)) query++;
    if ( query >= alignment->number_of_seqs ) {
	fprintf (stderr, "Query %s not found in the alignment.\n", query_name);
	return 1;
    }
    query_seq = alignment->sequence[query];

    /*compare */
    prot_pos = 0;
    prot_seq = protein->sequence;
    for (almt_pos=0; almt_pos < alignment->length; almt_pos++ ) {
	if ( query_seq [almt_pos] == '.' ) {
	    almt2prot [almt_pos] = -1;
	} else {
	    if ( prot_seq[prot_pos].res_type_short ==  query_seq [almt_pos] ) {
		prot2almt[prot_pos] = almt_pos;
		almt2prot[almt_pos] = prot_pos;
	    } else {
		fprintf (stderr, "Structure/alignment mismatch,\n");
		fprintf (stderr, "\t structure: pdbid %s  value %c \n", prot_seq[prot_pos].pdb_id,  prot_seq[prot_pos].res_type_short);
		fprintf (stderr, "\t alignment: pos %d  value %c \n",   almt_pos,  query_seq[almt_pos]);
		return 1;
	    }
	    prot_pos++;
	}
    }
    return 0;
}
