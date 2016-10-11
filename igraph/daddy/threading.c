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
    char *output_buffer = NULL;
    // First, read the length of the text message from the socket.
    // If read returns zero, the client closed the connection.
    int client_socket_id = thread_handler_arg_ptr->client_connection_id;
    int done = 0;
    int panic_ctr = 0;
    while (!done) {
	int retv = recv(client_socket_id, current_write_pos, BUF_BLOCK, 0);
	if ( retv  <  0) {
	    perror ("error reading");
	    goto cleanup_and_exit;
	}
	if (++panic_ctr > 10) goto cleanup_and_exit;
	done = (retv<BUF_BLOCK); 
    }
    // parse the input
    igraph_arg ig_args;
    memset (&ig_args, 0, sizeof(igraph_arg));
    if (input_parser (input_buffer, &ig_args)) goto cleanup_and_exit;
    ig_args.graph_ptr = thread_handler_arg_ptr->graph_ptr;
    // solve - if an errr exists and can be classified, it will be in the return buffer
    if (solver(&ig_args, &output_buffer)) goto cleanup_and_exit;
    if (!output_buffer) goto cleanup_and_exit;
    // write to connection
    int size  =  strlen(output_buffer)+1;
    // MSG_NOSIGNAL flag Requests not to send SIGPIPE on errors on stream oriented sockets
    // when the other end breaks the connection. 
    printf ("going to send, size= %d \n", size);
    if (send (client_socket_id, output_buffer, size, MSG_NOSIGNAL) < 0 ) {
 	perror ("write error");
	goto cleanup_and_exit;
    }

cleanup_and_exit:
    if (ig_args.node_list) free(ig_args.node_list);
    if (input_buffer) free (input_buffer);
    if (output_buffer) free (output_buffer);
    free (thread_handler_arg_ptr);
    return NULL;
}

