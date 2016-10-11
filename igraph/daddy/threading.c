#include "igraph_daddy.h"

////////////////////////////////////
int  create_socket_connection() {
    struct sockaddr_un addr;
    int fd = socket(AF_UNIX, SOCK_STREAM, 0);
    if (fd==-1) {
 	perror("socket error");
	exit(-1);
    }
    memset(&addr, 0, sizeof(addr));
    addr.sun_family = AF_UNIX;
    *addr.sun_path = '\0';
    strncpy(addr.sun_path+1, SOCKET_PATH+1, sizeof(addr.sun_path)-2);
    if (bind(fd, (struct sockaddr*)&addr, sizeof(addr)) == -1) {
	perror("bind error");
	exit(-1);
    }
    if (listen(fd, 5) == -1) {
	perror("listen error");
	exit(-1);
    }
    return fd;
}
////////////////////////////////////
void *thread_handler(void *inptr) {
    printf ("Hello from thread_handler\n");
    thread_handler_arg *thread_handler_arg_ptr = (thread_handler_arg *) inptr;
    // read from connection
    char *input_buffer  = emalloc (BUF_BLOCK);  curr_input_buf_size=BUF_BLOCK;
    char *current_write_pos = input_buffer;
    char *output_buffer;
    // First, read the length of the text message from the socket.
    // If read returns zero, the client closed the connection.
    int client_socket_id = thread_handler_arg_ptr->client_connection_id;
    int retval;
    int done = 0;
    int panic_ctr = 0;
    while (!done) {
	if ( recv(client_socket_id, current_write_pos, BUF_BLOCK, 0) <  0) {
	    perror ("error reading");
	    goto cleanup_and_exit;
	}
	panic_ctr ++;
	done = (strchr (input_buffer, '\n') != NULL) || (panic_ctr > 10); 
	if(!done) increase_buf_size (&input_buffer, &current_write_pos);
    }
    // parse the input
    igraph_args ig_args;
    memset (ig_args, 0, sizeof(igraph_args));
    if (input_parser (input_buffer, &ig_args)) goto cleanup_and_exit;
    // solve
    // write to connection
    output_buffer = emalloc(curr_input_buf_size+BUF_BLOCK);
    sprintf(output_buffer, "%s pong: %d ", input_buffer, getpid());
    int size  =  strlen(output_buffer)+1;
    
    sigignore(SIGPIPE); // we don't want to die here if the client quit
    printf ("going to send\n");
    if (send (client_socket_id, output_buffer, size, 0) < 0 ) {
 	perror ("write error");
	goto cleanup_and_exit;
    }
    printf ("back from send\n");
    printf ("send retval: %d\n", retval); 

cleanup_and_exit:
    if (ig_args.node_list) free(ig_args.node_list);
    free (input_buffer);
    free (output_buffer);
    free (thread_handler_arg_ptr);
    return NULL;
}

