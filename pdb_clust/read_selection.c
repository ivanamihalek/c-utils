# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "pdbclust.h"
/*
  input format
F  3
K  48
R  50
D  59
C  62
G  71
A  78
Y  79
E  82
G  91
A  97
C  100
G  102
L  109
etc 
 */


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
	res_short = token[0][0];
	for (i=0; i<no_res; i++) {
	    if (! strcmp(sequence[i].pdb_id, token[1]) ){
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
    int no_selected = 0;
    for (i=0; i<no_res; i++) no_selected += selection[i];
    printf ("selected: %d\n", no_selected);
    
    return 0;
}
