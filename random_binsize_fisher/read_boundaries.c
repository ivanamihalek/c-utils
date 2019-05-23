
# include <stdio.h>
# include <stdlib.h>

#include "utils.h"
/****************************************************************************/
int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
             char * fmt, char * warnstr) {

        fprintf ( fptr, "Error on line %3d:     %s", line_ctr, line);
        fprintf ( fptr, fmt, warnstr);
        return 0;
}
/****************************************************************************/
int read_boundaries (char * filename, int ** bds_ptr, int * bds_size ) {
        FILE * fptr, *log = stdout;
        char line[LONGSTRING];
        char token[MAX_TOK][MEDSTRING] = {{'\0'}};
        char comment_char;
        int line_ctr = 0, retval;
        int max_token;
        /***************/

        fptr   = efopen ( filename, "r" );
        if (!fptr ) return 1;
        memset ( line, 0, LONGSTRING);
        *bds_size = 0;
        /* get the size and check the contents */
        while(fgets(line,LONGSTRING,fptr)!=NULL) {
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
                (*bds_size) ++;
        }
        int* boundaries = (int*)ecalloc(*bds_size, sizeof(int));
        *bds_ptr = boundaries;
        rewind ( fptr);
        while(fgets(line,LONGSTRING,fptr)!=NULL) {
                tokenize ( token, &max_token, line, comment_char= '!' );
                if ( max_token < 0 ) continue;
                *boundaries = atoi(token[0]);
                boundaries++;
        }
        fclose(fptr);

        return 0;
}
