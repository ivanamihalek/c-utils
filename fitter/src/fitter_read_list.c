# include "fitter.h"

int read_list ( char * filename, char *** list, int *length_ptr) {

    FILE * fptr, *log = stdout; 
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    char comment_char;
    char **my_list;
    int  max_token;
    int line_ctr, retval, list_length;
    int ctr;
    int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
		 char * fmt, char * warnstr) ;

    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;
    
    /* count the number of residues in the list */
    line_ctr = 0;
    memset ( line, 0, LONGSTRING);
    list_length = 0;
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

	list_length += max_token+1;
    }

    if ( ! (my_list = chmatrix (list_length, MEDSTRING))) return 1;

    rewind (fptr );

    line_ctr = 0;
    memset ( line, 0, LONGSTRING);
    list_length = 0;
    while(fgets(line,LONGSTRING,fptr)!=NULL){
	line_ctr++;
	/* tokenize */
	retval = tokenize ( token, &max_token, line, comment_char= '!' );
	
	if ( max_token < 0 ) continue;

	for (ctr=0;  ctr<= max_token; ctr++) {
	    sprintf ( my_list[list_length], "%s", token[ctr] );
	    list_length ++;
	}
    }
   
    *list = my_list;
    *length_ptr = list_length;
   
    fclose (fptr);
    return 0;
    
    
}

/****************************************************************************/
int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
	     char * fmt, char * warnstr) {

    fprintf ( fptr, "\tError on line %3d:     %s", line_ctr, line);
    fprintf ( fptr, fmt, warnstr);
    return 0;
}
