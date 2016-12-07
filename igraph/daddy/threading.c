#include "igraph_daddy.h"
#include <sys/unistd.h>
#include <sys/fcntl.h>

extern  pool threadpool;
extern  igraph_t graph;


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
    thread_args * thread_args_ptr = (thread_args *)thread_args_ptr_cast_to_void;
# ifdef VERBOSE
    printf ("in send timeout message\n");
    printf (" *** job done:  %d ,  time:  %d\n", thread_args_ptr->job_done, (int)(clock()-thread_args_ptr->start_time));
# endif
    if (thread_args_ptr->job_done ) {
    } else {
	int client_socket_id = thread_args_ptr->client_connection_id;
	char * output_buffer = emalloc(100);
    	sprintf(output_buffer,"iGraph operation timed out");
	int size  =  strlen(output_buffer)+1;
# ifdef VERBOSE
        printf  ("sending %s, size %d\n\n", output_buffer, size);  fflush(stdout);
# endif
	send (client_socket_id, output_buffer, size, MSG_NOSIGNAL);
	close(client_socket_id);
	free(output_buffer);
    }
}


/////////////////
void *thread_payload_function(void *inptr) {

    thread_args *main_thread_args_ptr = (thread_args *) inptr;
    int oldtype; // a return value we don't really need here
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &oldtype);
     
    // I am counting on thread being detached to free all resources,
    // so the only thing we do on canclation is send the info to the client
    // faterwards, the thread is terminate  - it exits, see
    // http://man7.org/linux/man-pages/man3/pthread_cancel.3.html
    
    pthread_rwlock_wrlock (&threadpool.args_lock);
    thread_args * thread_args_ptr = emalloc (sizeof(thread_args));
# ifdef VERBOSE
    printf ("thread associated with  socket id %d  allocated  %d bytes to  %p\n",
    	    main_thread_args_ptr->client_connection_id, (int)sizeof(thread_args), thread_args_ptr);
# endif
    memcpy (thread_args_ptr, main_thread_args_ptr, sizeof(thread_args));
    main_thread_args_ptr -> copied_to = (void*) thread_args_ptr;
    thread_args_ptr -> pthread_id = pthread_self();
    pthread_rwlock_unlock (&threadpool.args_lock);
     
    // read from connection
    char *output_buffer = NULL;
    int client_connection_id = thread_args_ptr->client_connection_id;
    char *input_buffer  = emalloc (BUF_BLOCK);
    int curr_input_buf_size = BUF_BLOCK;
    char *current_write_pos = input_buffer;
    // First, read the length of the text message from the socket.
    // If read returns zero, the client closed the connection.
    int done = 0;
    int panic_ctr = 0;
    while (!done) {
    	int retv = recv(client_connection_id, current_write_pos, BUF_BLOCK, 0);
    	if ( retv  <  0) {
    	    perror ("error reading");
    	    goto cleanup_and_exit;
    	}
    	if (++panic_ctr > 10) {
    	    char * errmsg =  "Error: Input too long.\n";
    	    send (client_connection_id, errmsg, strlen(errmsg), MSG_NOSIGNAL);
    	    goto cleanup_and_exit;
    	}
    	done = (retv<BUF_BLOCK);
    	if (!done) increase_buf_size (&input_buffer, &curr_input_buf_size, &current_write_pos);
    }
    // parse the input
    igraph_arg ig_args;
    memset (&ig_args, 0, sizeof(igraph_arg));
    if (input_parser (input_buffer, &ig_args)) goto cleanup_and_exit;
    ig_args.graph_ptr = &graph;

    // solve - if an errr exists and can be classified, it will be in the return buffer
    pthread_cleanup_push (send_timeout_message, (void*) thread_args_ptr);    
    if (solver(&ig_args, &output_buffer) || !output_buffer) {
    	output_buffer = emalloc(100);
    	sprintf(output_buffer, "Unspecified error runnig iGraph");
    }
    thread_args_ptr->job_done = 1;
    // cleanup_pop removes the routine at the top of
    // the stack of clean-up handlers, and optionally executes it if execute
    //  is nonzero
    pthread_cleanup_pop(1);
     
    // write to connection
    int size  =  strlen(output_buffer)+1;
    // MSG_NOSIGNAL flag Requests not to send SIGPIPE on errors on stream oriented sockets
    // when the other end breaks the connection. 
    if (send (client_connection_id, output_buffer, size, MSG_NOSIGNAL) < 0 ) {
 	perror ("write error");
	goto cleanup_and_exit;
    }

cleanup_and_exit:
    close(client_connection_id);
    //if (ig_args.node_list) free(ig_args.node_list);
    //if (input_buffer) free (input_buffer);
    if (output_buffer) free (output_buffer);
    thread_args_ptr -> job_done = 1;
    return NULL;
}

