# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "pdbclust.h"

/****************************************************************************/
int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
	     char * fmt, char * warnstr) {

    fprintf ( fptr, "Error on line %3d:     %s", line_ctr, line);
    fprintf ( fptr, fmt, warnstr);
    return 0;
}
/****************************************************************************/
int read_selection (Residue *sequence, int no_res, char * filename, int * selection) {
    FILE * fptr, *log = stdout;
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    char comment_char;
    int  line_ctr = 0, retval;
    int  max_token;
    /***************/

    int i;
    char res_short;
    
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;
    memset ( line, 0, LONGSTRING);
    while(fgets(line,LONGSTRING,fptr)!=NULL){
 	line_ctr++;
 	/* tokenize */
	retval = tokenize ( token, &max_token, line, comment_char= '!' );
	switch ( retval ) {
	case  TOK_TOOMNY:
	    errmsg ( log, line_ctr, line, "\t\t %s\n", "Too many tokens.");
	    fclose (log);
	    break;
	case TOK_TOOLONG:
	    errmsg ( log, line_ctr, line, "\t\t %s\n", "Token too long.");
	    fclose (log);
	    break;
	}
	if ( max_token < 0 ) continue;
	res_short = token[1][0];
	for (i=0; i<no_res; i++) {
	    if (! strcmp(sequence[i].pdb_id, token[0]) ){
		/* check */
		if ( res_short != sequence[i].res_type_short) {
		    fprintf (stderr, "type mismatch: %s  %c  %c \n",
			    sequence[i].pdb_id, res_short, sequence[i].res_type_short);
		    exit(1);
		}
		selection[i] = 1;
		break;
	    }
	}
    }

    return 0;
}
