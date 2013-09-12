# include "postp.h"

int read_clustalw ( char * cwname, Alignment * alignment){
    
    FILE * fptr = NULL;
    char line[BUFFLEN];
    int  number_of_seqs, almt_length, ctr;
    int * seq_pos, pos;
    char ** sequence;
    char ** name;
    char curr_name[BUFFLEN];
     
    /* open file */
    fptr = efopen ( cwname, "r");
    if ( !fptr ) return 1;
    
    /* find the alignment length info */
    almt_length = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr(line, "MSF:" ) ){
	    sscanf (line, "%*s %d", &almt_length);
	    break;
	}
    }
    if ( almt_length ) {
	/* printf ( "Alignment length in %s is %d.\n", cwname, almt_length); */
    } else {
	fprintf ( stderr, "Alignment length info not found in %s. Is the format gcg?\n", cwname);
	return 1;
    }

    /* determine the number of sequences */
    number_of_seqs = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( ! strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) number_of_seqs++;
    }
    if ( number_of_seqs ) {
	/* printf ( "Number of sequences in %s is %d.\n", cwname, number_of_seqs); */
    } else {
	fprintf ( stderr, "No sequences found in %s. Is the format gcg?\n", cwname);
	return 1;
    } 
    
    /* allocate */
    sequence = chmatrix (number_of_seqs, almt_length);
    if ( !sequence ) return 1;
    name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
    if ( !name ) return 1;
    seq_pos = (int *) calloc ( number_of_seqs, sizeof(int));
    if ( !seq_pos ) return 1;
    
    /* read in */
    rewind(fptr);
    ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if (!  strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    sscanf (line, "%*s %s", name[ctr]);
	    ctr ++;
	}
    }
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( isspace (line[0] ) ) continue;
	sscanf (line, "%s", curr_name);
	ctr = 0;
	while (  ctr <number_of_seqs &&  strcmp (name[ctr], curr_name) ) ctr++;
	if ( ctr >= number_of_seqs ) {
	    fprintf ( stderr, "The name %s not found in the header of %s.\n", curr_name, cwname);
	    return 1;
	}
	pos = 0;
	while ( ! isspace(line[pos]) ) pos++;
	while  (line[pos] != '\n' && pos < BUFFLEN) {
	    if ( !  isspace(line[pos] ) ){
		if ((line[pos]>=97)&&(line[pos]<=122)) {line[pos] -= 32;} /* --> turn to uppercase */
		if ( line[pos]==126)                   {line[pos]  = 46;} /* turn tweedle to dot */
		sequence [ctr] [ seq_pos[ctr] ] = line[pos];
		seq_pos[ctr]++;
	    }
	    pos ++;
	}
    }
    fclose(fptr);

    /* sanity check */
    for (ctr=0; ctr < number_of_seqs; ctr++ ) {
	if ( seq_pos[ctr] != almt_length ) {
	    fprintf (stderr, "Sequence %s is shorter (%d position) than the alignment.\n", name[ctr],  seq_pos[ctr]);
	    return 1;
	}
    }

    /* return values */
    alignment->number_of_seqs = number_of_seqs;
    alignment->length         = almt_length;
    alignment->sequence       = sequence;
    alignment->name           = name;

    /* free */
    free (seq_pos);

    count_gaps (alignment);
    
    return 0;
}
