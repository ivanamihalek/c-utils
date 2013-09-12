# include "postp.h"

int read_cmd_file (char *filename, Options * options) {
    FILE * fptr, *log = stdout; 
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    char comment_char;
    int  max_token, current_char = 0;
    int namectr, line_ctr, retval;
    int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING], char * fmt, char * warnstr) ;
   
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;
    
    memset (options, 0, sizeof(Options));

    /* count the number of name_files */ 
    
    line_ctr = 0;
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
	
	
	/* turn  first token to lowercase */
	current_char = 0;
	while ( token[0][current_char] &&  (current_char < SHORTSTRING) ){
	    if ( token[0][current_char] < 97)  token[0][current_char] += 32;
	    current_char++;
	}
	if (  ! strncmp (token[0], "name", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a file name (full path).\n", token[0]);
		return 1;
	    }
	    options->no_of_alignments ++;
	}	
    }
    /* allocate space for the msf file names */
    if ( !( options->namefile = emalloc ( options->no_of_alignments*sizeof(char*) ) ) ) return 1;
    for ( namectr = 0; namectr < options->no_of_alignments; namectr++ ) {
	 if ( !( options->namefile[namectr] = emalloc (BUFFLEN*sizeof(char) ) ) ) return 1;
    }

    rewind (fptr);
    namectr = 0;
    line_ctr = 0;
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
	
	
	/* turn  first token to lowercase */
	current_char = 0;
	while ( token[0][current_char] &&  (current_char < SHORTSTRING) ){
	    if ( token[0][current_char] < 97)  token[0][current_char] += 32;
	    current_char++;
	}

	/* check first token for meaning (first token should be a keyword)*/
	if (  ! strncmp (token[0], "align", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a file name (full path).\n", token[0]);
		return 1;
	    }
	    sprintf ( options->almtname, "%s", token[1]);
	} else if  (  ! strncmp (token[0], "chai", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a chain identifier.\n", token[0]);
		return 1;
	    }
	    options->chain = token[1][0];
	    
	} else if  (  ! strncmp (token[0], "insc", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a file name (full path).\n", token[0]);
		return 1;
	    }
	    sprintf ( options->scorename, "%s", token[1]);
	    
	} else if (  ! strncmp (token[0], "meth", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a name.\n", token[0]);
		return 1;
	    }
	    /* turn  second token to lowercase */
	    current_char = 0;
	    while ( token[1][current_char] &&  (current_char < SHORTSTRING) ){
		if ( token[1][current_char] < 97)  token[1][current_char] += 32;
		current_char++;
	    }
	    if ( ! strncmp ( token[1], "entr", 4)  )  {
		options->scoring_method = ENTROPY;
	    } else if ( ! strncmp ( token[1], "ivet", 4)  )  {
		options->scoring_method = IVET;
	    } else if ( ! strncmp ( token[1], "rvet", 4)  )  {
		options->scoring_method = RVET;
	    }
	    
	} else if (  ! strncmp (token[0], "name", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a file name (full path).\n", token[0]);
		return 1;
	    }
	    sprintf ( options->namefile[namectr], "%s", token[1] );
	    namectr++;
	
	} else if (  ! strncmp (token[0], "outn", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a name.\n", token[0]);
		return 1;
	    }
	    sprintf ( options->outname, "%s", token[1]);
	    
	} else if  (  ! strncmp (token[0], "pdbf", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a file name (full path).\n", token[0]);
		return 1;
	    }
	    sprintf ( options->pdbname, "%s", token[1]);
	    
	} else if (  ! strncmp (token[0], "query", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\t\t Keyord %s should be followed by a name.\n", token[0]);
		return 1;
	    }
	    sprintf ( options->query, "%s", token[1]);
	    
	} else if (  ! strncmp (token[0], "sink", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
		        "\t\t Keyord %s should be followed  by a percentage of gaps which \"sinks\" a position.\n",
			 token[0]);
		return 1;
	    }
	    options->max_gaps = atof (token[1]);
	    
	} else if (  ! strncmp (token[0], "spec", 4)  ) {
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
		        "\t\t Keyord %s should be followed  by the name of the \"special\" sequence\n",
			 token[0]);
		return 1;
	    }
	    sprintf ( options->special, "%s", token[1]);
	    
	} else {
	    errmsg ( log, line_ctr, line, "\t\t Keyword %s not recognized.\n", token[0]);
	    return 1;
		    
	}
	memset (line, 0, LONGSTRING);
    }
    fclose (fptr);
    
    /* checking and default setting: */
     if ( !options->query[0] && options->pdbname[0] ) {
	fprintf ( stderr, "Pdb sequence should be part of the alignment, and the query name provided to do the mapping.\n");
	return 1;
    }
     if ( !options->scorename[0] && !options->almtname ) {
	fprintf ( stderr, "Neither score nor the alignment are provided.(What am I supposed to do?).\n");
	return 1;
    }
   if ( !options->outname[0]) sprintf (options->outname, "%s", "postp");
   if ( !options->scoring_method ) options->scoring_method = RVET;


  return 0;

}

/***************************************************************************************************/
int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING], char * fmt, char * warnstr) {

    fprintf ( fptr, "\tError on line %3d:     %s", line_ctr, line);
    fprintf ( fptr, fmt, warnstr);
    return 0;
}
