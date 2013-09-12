# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "sqc.h"



int read_msf ( char msfname[], Sequence ** seqptr, int * noseqptr, int * seqlenptr) {

    FILE * fptr = NULL;
    char line [BUFLEN] = {'\0'};
    int seq_ctr, no_seqs;
    int seq_len;
    char *token, name[NAMELEN], *pos;
    Sequence * sequence;
    int pos_ctr;

    
    fptr = fopen (msfname, "r");
    if ( !fptr) {
	fprintf (stderr, "Cannot open %s.\n", msfname);
	return 1;
    }

    seq_ctr = 0;
    seq_len = 0;
    while(fgets(line,BUFLEN,fptr)!=NULL){
	if (  strstr (line,"Type")  ) {
	    strtok (line, " ");
	    *seqlenptr = seq_len = atoi ( strtok (NULL, " ") );
	}
	if (  strstr (line,"Name")  ) {
	    seq_ctr ++;
	}
    }
    *noseqptr = no_seqs = seq_ctr;

    /* sanity */
    if ( !seq_len ) {
	    fprintf (stderr, "Zero sequence length in %s?!\n", msfname);
	    return 1;
    }
    if ( !no_seqs ) {
	    fprintf (stderr, "No sequences in %s?!\n", msfname);
	    return 1;
    }

    /* allocate space */
    *seqptr = sequence = (Sequence*) calloc ( no_seqs, sizeof(Sequence) );
    if ( !sequence) {
	fprintf (stderr, "error allocating sequence space.\n");
	return 1;
    }
    for (seq_ctr=0; seq_ctr < no_seqs; seq_ctr++) {
	sequence[seq_ctr].position = (char *) calloc (seq_len, sizeof(char));
	if ( !sequence[seq_ctr].position ) {
	    fprintf (stderr, "error allocating sequence space.\n");
	    return 1;
	}
    }
    

    /* read in the names */
    rewind (fptr);
    seq_ctr = -1;
    while(fgets(line,BUFLEN,fptr)!=NULL){
	if (  strstr (line,"Name") ) {
	    seq_ctr++;
	    strtok (line, " ");
	    sprintf ( sequence[seq_ctr].name, "%s",  strtok (NULL, " ")  );
	} 
    }

    /* look for lines containing info */
    rewind(fptr);
    while(fgets(line,BUFLEN,fptr)!=NULL){
	
	sprintf ( name, "%s",  strtok (line, " "));
	
	for (seq_ctr=0; seq_ctr < no_seqs; seq_ctr++) {
	
	    if (  ! strcmp (name, sequence[seq_ctr].name) ) {
		
		/* write to the seqence called  <name> */
		pos = sequence[seq_ctr].position;
		while ( *pos) pos++;
		if (pos >= sequence[seq_ctr].position+seq_len) {
		    fprintf (stderr,"While reading %s: %s is longer  than declared (%d should be max).\n",
			     msfname, sequence[seq_ctr].name, seq_len);
		    return 1;
		}
		while (  (token=strtok (NULL, " \n")) ) {
		    strncpy ( pos, token, strlen(token));
		    pos += strlen(token);
		    *pos = '\0';
		}
		break;
	    }
	}
    }


    /*sanity*/
    for (seq_ctr=0; seq_ctr < no_seqs; seq_ctr++) {
	if ( strlen (sequence[seq_ctr].position) != seq_len ) {
	    fprintf (stderr, "While reading %s: %s is shorter (%d positions) than declared (%d).\n",
		     msfname, sequence[seq_ctr].name, strlen (sequence[seq_ctr].position), seq_len);
	    return 1;
	}
    }
   
    /* gap notation: */
    for (seq_ctr=0; seq_ctr < no_seqs; seq_ctr++) {
	for (pos_ctr=0; pos_ctr<seq_len; pos_ctr++) {
	    if ( sequence[seq_ctr].position[pos_ctr] == '=') {
		sequence[seq_ctr].position[pos_ctr] = '.' ;
	    }
	}
	
    }
    
    fclose (fptr);
    return 0;
}
