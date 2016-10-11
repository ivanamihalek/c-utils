#include "igraph_daddy.h"

// hardcoded
char *SOCKET_PATH = "\0igraph_daddy";

//////////////////////////////////////////////////
int main ( int argc, char * argv[]) {
    // arguments: path to graph
    // the graph is just a list of from-to identifier pairs, one per line
    if (argc != 2) {
	fprintf (stderr, "Usage: %s  <path to graph file>\n",  argv[0]);
	exit(1);
    }
    // connect to the rest of the world
    int socket_id = create_socket_connection();
    
    // read in and reconstruct the graph
    igraph_t graph;
    construct_graph(argv[1], &graph);
    
    // main loop acepting the connection and solving the task    
    while (1) {
	int client_connection_id = accept(socket_id, NULL, NULL);
	if (client_connection_id < 0)  continue; // accept error
	pthread_t solver_thread;
	// thread_handler arg structure containts both the client_connection_id, 
	// and the pointer to the graph structure
	// we want to allocate it anew each time, and will let
	// the receiving thread free the mem - not ideal, but not sure
	// if it is worth further fussing
	thread_handler_arg * thread_handler_arg_ptr = emalloc (sizeof(thread_handler_arg));
	thread_handler_arg_ptr -> client_connection_id = client_connection_id;
	thread_handler_arg_ptr -> graph_ptr = &graph;
	/* int pthread_create(pthread_t *thread, const pthread_attr_t *attr,
                          void *(*start_routine) (void *), void *arg);
	  The pthread_create() function starts a new thread in the calling
          process.  The new thread starts execution by invoking
          start_routine(); arg is passed as the sole argument of start_routine().
	*/
	int pthred_retval = pthread_create (&solver_thread, NULL,
					    thread_handler, thread_handler_arg_ptr);
	if (pthred_retval<0) {
	    perror("could not create thread");
	    return 1;
	}
    }
    return 0;
}
