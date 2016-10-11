#include "igraph_daddy.h"

// hardcoded
int BUF_BLOCK = 4;
int curr_input_buf_size = 0;


// generic helper functions
void * emalloc(int size) {
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %d bytes", size);
	exit(1);
    }
   return ptr;
}
void increase_buf_size (char ** buffer, char ** current_write_pos) {
    char  *tmp_buffer_ptr;
    tmp_buffer_ptr = emalloc(curr_input_buf_size + BUF_BLOCK);
    memcpy (tmp_buffer_ptr, *buffer, curr_input_buf_size);
    *current_write_pos =  tmp_buffer_ptr + curr_input_buf_size;
    free (*buffer);
    *buffer = tmp_buffer_ptr;
    curr_input_buf_size += BUF_BLOCK;
    return;
}
