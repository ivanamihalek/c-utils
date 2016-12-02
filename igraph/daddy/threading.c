#include "igraph_daddy.h"
#include <sys/unistd.h>
#include <sys/fcntl.h>

////////////////////////////////////
int  create_socket_connection(const char * socket_path) {
    struct sockaddr_un addr;
    int fd = socket(PF_UNIX, SOCK_STREAM, 0);
    if (fd==-1) {
 	perror("socket error");
	exit(-1);
    }
    unlink(socket_path);
    memset(&addr, 0, sizeof(addr));
    addr.sun_family = PF_UNIX;
    //? what is the size of addr.sun_path?
    //The C library function char *strncpy(char *dest, const char *src, size_t n)
    // copies up to n characters from the string pointed to, by src to dest.
    //In a case where the length of src is less than that of n, the remainder of dest will be padded with null bytes.
    strncpy(addr.sun_path, socket_path, sizeof(addr.sun_path)-1);
	
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
// undo max numbre of nbrs, max time
// I am not getting in here
// redo valgrind
// freeinf of attribute space
//
////////////////////////////////////
void send_timeout_message (void * thread_args_ptr_cast_to_void) {
    // I am getting here at the end of the thread, even if it wasn't canceled
    // I'll have to figure it out some other time - for now
    // I'm just relying on  the job done flag being set
    printf ("in send timeout message\n");
    thread_args * thread_args_ptr = (thread_args *)thread_args_ptr_cast_to_void;
    printf ("in send timeout message:  %d \n", thread_args_ptr->job_done);
    if (thread_args_ptr->job_done ) {
    } else {
	int client_socket_id = thread_args_ptr->client_connection_id;
	char * output_buffer = emalloc(100);
	if (output_buffer) {
	    printf ("sending \n");
	    printf ("sending %s\n", output_buffer);
	    sprintf (output_buffer, "iGraph operation timed out");
	    int size  =  strlen(output_buffer)+1;
	    send (client_socket_id, output_buffer, size, MSG_NOSIGNAL);
	    free (output_buffer);
	}
	close(client_socket_id);
    }
}


/////////////////
void *thread_handler(void *inptr) {

    thread_args *thread_args_ptr = (thread_args *) inptr;
    int oldtype; // a return value we don't really need here
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &oldtype);
    printf("in thread\n");
     
    // I am counting on thread being detached to free all resources,
    // so the only thing we do on canclation is sent the info to the client
    // faterwards, the thread is terminate  - it exits, see
    //http://man7.org/linux/man-pages/man3/pthread_cancel.3.html
       
    // read from connection
    char *input_buffer  = emalloc (BUF_BLOCK);
    int curr_input_buf_size=BUF_BLOCK;
    char *current_write_pos = input_buffer;
    char *output_buffer = NULL;
    // First, read the length of the text message from the socket.
    // If read returns zero, the client closed the connection.
    int client_socket_id = thread_args_ptr->client_connection_id;
    int done = 0;
    int panic_ctr = 0;
    while (!done) {
	int retv = recv(client_socket_id, current_write_pos, BUF_BLOCK, 0);
	if ( retv  <  0) {
	    perror ("error reading");
	    goto cleanup_and_exit;
	}
	if (++panic_ctr > 10) {
	    char * errmsg =  "Error: Input too long.\n";
	    send (client_socket_id, errmsg, strlen(errmsg), MSG_NOSIGNAL);
	    goto cleanup_and_exit;
	}
	done = (retv<BUF_BLOCK);
	if (!done) increase_buf_size (&input_buffer, &curr_input_buf_size, &current_write_pos);
    }
    // parse the input
    igraph_arg ig_args;
    memset (&ig_args, 0, sizeof(igraph_arg));
    if (input_parser (input_buffer, &ig_args)) goto cleanup_and_exit;
    ig_args.graph_ptr = thread_args_ptr->graph_ptr;
    
    // solve - if an errr exists and can be classified, it will be in the return buffer
    pthread_cleanup_push (send_timeout_message, (void*) thread_args_ptr);
    
    if (solver(&ig_args, &output_buffer) || !output_buffer) {
    	output_buffer = emalloc(100);
    	sprintf(output_buffer, "Unspecified error runnig iGraph");
    }
    thread_args_ptr->job_done = 1;
    // cleanup_pop removes the routine at the top of
    //   the stack of clean-up handlers, and optionally executes it if execute
    //   is nonzero
    pthread_cleanup_pop(1);
     
    // write to connection
    int size  =  strlen(output_buffer)+1;
    // MSG_NOSIGNAL flag Requests not to send SIGPIPE on errors on stream oriented sockets
    // when the other end breaks the connection. 
    if (send (client_socket_id, output_buffer, size, MSG_NOSIGNAL) < 0 ) {
 	perror ("write error");
	goto cleanup_and_exit;
    }

cleanup_and_exit:
    close(client_socket_id);
    if (ig_args.node_list) free(ig_args.node_list);
    if (input_buffer) free (input_buffer);
    if (output_buffer) free (output_buffer);
    thread_args_ptr -> job_done = 1;
    return NULL;
}

