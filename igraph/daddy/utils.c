#include "igraph_daddy.h"

// hardcoded
int curr_input_buf_size = 0;


// generic helper functions
FILE * efopen(char * name, char * mode) {

    FILE * fp;
    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,    "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }
    return fp;
}

void * emalloc(int size) {
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %d bytes", size);
	exit(1);
    }
   return ptr;
}

void increase_buf_size (char ** buffer, int *curr_input_buf_size_ptr, char ** current_write_pos) {
    char  *tmp_buffer_ptr;
    int curr_input_buf_size = *curr_input_buf_size_ptr;
    tmp_buffer_ptr          = emalloc(curr_input_buf_size + BUF_BLOCK);
    memcpy (tmp_buffer_ptr, *buffer, curr_input_buf_size);
    *current_write_pos =  tmp_buffer_ptr + curr_input_buf_size;
    free (*buffer);
    *buffer = tmp_buffer_ptr;
    curr_input_buf_size += BUF_BLOCK;
    *curr_input_buf_size_ptr = curr_input_buf_size;
    return;
}

/***************************************************************************/
int contains_nondigit (char * string) {
    int pos;
    for (pos=0; pos<strlen(string); pos ++ ) {
	if (!isdigit(string[pos])) return 1;
    }
    return 0;
}
/***************************************************************************/
int tokenize ( char token[MAX_TOK][TOKENLENGTH], int * max_token,
	       char * line , char comment_char) {
    /* assumes the tokens to be no bigger than TOKENLENGTH */ 
    
    char * chrptr, *last; 
    int current_token, current_char = 0;
    int reading;
   
    memset (token[0], 0, MAX_TOK*TOKENLENGTH*sizeof(char)); 
    chrptr = line;
    last   = chrptr + strlen (line);
    current_token = -1;
    current_char  =  0;
    reading = 0;
    while ( chrptr <= last) {
	if ( *chrptr == comment_char ) break;
	if ( *chrptr == '\n' ) break;
	if ( *chrptr && ! isspace(*chrptr) ) {
	    if ( ! reading ) {
		reading = 1;
		current_char = 0;
		current_token++;
		if ( current_token >= MAX_TOK ) {
		    return 1; /* defined in possum_utils.h */
		}
	    }
	    if ( current_char >= TOKENLENGTH ) {
		return 1;
	    }
	    token[current_token][current_char] = *chrptr;
	    current_char++;
	} else {
	    if ( reading ) {
		reading = 0;
	    }
	}
	chrptr++;
    }
    *max_token = current_token;

    return 0;
    
}
