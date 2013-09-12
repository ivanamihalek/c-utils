# include "postp.h"

int count_gaps (Alignment * alignment) {

    int s, c;
    alignment->seq_gaps    = (int *) emalloc (alignment->number_of_seqs*sizeof(int));
    if (!alignment->seq_gaps) return 1;
    alignment->column_gaps = (int *) emalloc (alignment->length*sizeof(int));
    if (!alignment->column_gaps) return 1;
    alignment->total_gaps = 0;
    for ( s=0; s<alignment->number_of_seqs; s++ ) {
	for ( c=0; c<alignment->length; c++) {
	    if ( alignment->sequence[s][c] == '.' ) {
		alignment->column_gaps[c] ++;
		alignment->seq_gaps[s] ++;
		alignment->total_gaps  ++;
	    }
	}
    }
    return 0;
}

double  avg_dist_to_special ( Options * options, Alignment * alignment) {

    int i, special = -1;
    double avg;

    if ( ! options->special || options->special[0] == '\0' ) return -1.0;
    
    /* find the special */
    for ( i=0; i < alignment->number_of_seqs; i++ ) {
	if ( strcmp ( alignment->name[i], options->special) ) continue;
	special = i;
	break;
    }
    if ( special  == -1) return -1.0;

    /* determine the avg dist to it */
    avg = 0;
    for ( i=0; i < alignment->number_of_seqs; i++ ) {
	if ( i== special ) continue;
	avg += alignment->seq_dist [i][special];
    }
    avg /= (alignment->number_of_seqs - 1);
    return avg;
}

/*****************************************************************************/
int reconstruct_alignment (char * namesfile, Alignment * big_almtptr,  Alignment * almtptr) {

     FILE * fptr = NULL;
     char line[BUFFLEN];
     char token[MAX_TOK][MEDSTRING] = {{'\0'}};
     int  number_of_seqs, almt_length, ctr, ctr2;
     int retval, max_token;
     int found;
     char comment_char, *name;
     
     /* open file */
     fptr = efopen ( namesfile, "r");
     if ( !fptr ) return 1;

     /* the length is the length of the big (overall) alignment */ 
     almt_length =  almtptr->length = big_almtptr->length;
     
     /* determine the number of sequences */
     number_of_seqs = 0;
     while(fgets(line, BUFFLEN, fptr)!=NULL){
	 retval = tokenize ( token, &max_token, line, comment_char= '!' );
	 switch ( retval ) {
	 case  TOK_TOOMNY:
	     fprintf (stderr, "Too many strings on a signle line in %s.\n", namesfile);
	     return 1;
	     break;
	 case TOK_TOOLONG:
	     fprintf (stderr, "Sequence name too long in %s.\n", namesfile);
	     return 1;
	     break;
	 }
	 if ( max_token < 0 ) continue;
	 number_of_seqs ++;
    }

     if ( number_of_seqs ) {
	/*  printf ( "Number of names in %s is %d.\n", namesfile, number_of_seqs); */
     } else {
	 fprintf ( stderr, "No name found in %s.\n", namesfile);
	 return 1;
     }
     almtptr->number_of_seqs = number_of_seqs;
     
     /* allocate */
     almtptr->name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
     if ( !almtptr->name ) return 1;
     almtptr->sequence = (char**) emalloc ( number_of_seqs*sizeof (char *) );
     if ( !almtptr->sequence ) return 1;
   
     /* read the names in */
     rewind(fptr);
     ctr = 0;
     while(fgets(line, BUFFLEN, fptr)!=NULL){
	 retval = tokenize ( token, &max_token, line, comment_char= '!' );
	 if ( max_token < 0 ) continue;
	 sprintf ( almtptr->name[ctr], "%s", token[0]);
	 ctr++;
     }
     fclose (fptr);

     /* direct the sequence ptrs to the corresponding places in the big alignment */
     for (ctr=0; ctr < number_of_seqs; ctr++) {
	 name = almtptr->name[ctr];
	 found = 0;
	 for (ctr2=0; ctr2<big_almtptr->number_of_seqs; ctr2++ ) {
	     if ( ! strcmp(name, big_almtptr->name[ctr2] ) ) {
		 almtptr->sequence[ctr] = big_almtptr->sequence[ctr2];
		 found = 1;
		 break;
	     }
	 }
	 if ( !found ) {
	     fprintf ( stderr, "Name %s  from %s not found in the overall alignment.\n",
		       name, namesfile);
	     return 1;
	 }
     }

     count_gaps (almtptr);
    
     
     return 0;
}
