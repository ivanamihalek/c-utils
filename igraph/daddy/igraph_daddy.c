#include "igraph_daddy.h"

// hardcoded
char *SOCKET_PATH = "/tmp/igraph_daddy";
int MAX_THREADS   = 10;
    
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

    printf ("daddy's listening ... \n");
    
    // main loop acepting the connection and solving the task
    thread_handler_arg **thread_pool = emalloc(MAX_THREADS*sizeof(thread_handler_arg));
    int number_of_threads = 0;
    while (1) {
	int client_connection_id = accept(socket_id, NULL, NULL);
	if (client_connection_id < 0)  continue; // accept error
        printf (" opening connection %d \n", client_connection_id);

	// if the number of assigned threads is akready maxed out, wait for the queue to drain
	int done = 0;
	while (!done) {
	    if (number_of_threads < MAX_THREADS) {
		assign_thread (thread_handler_arg **thread_pool, int * number_of_threads);
		done = 1;
	    } else {
		thread_pool_cleanup(thread_pool, &number_of_threads);
	    }
	    // else check the timee: if timeout, done
	    //timeout = ;
	}
	
    }
    return 0;
}
//rename thread_handler_arg to thread_args; re-checkk with valgrind when done

int thread_timeout(thread_handler_arg * thread_handler_arg_ptr) {
	// if time(now) - thread_handler_arg-> time(started) > TIME_MAX return true
	// otherwise return false
	
}
    

int thread_pool_cleanup (thread_handler_arg **thread_pool, int * number_of_threads) {
    int t;
    int new_number_of_threads = number_of_threads;
    for (t=0; t< number_of_threads; t++) {
	int done = thread_pool[t] -> job_done;
	if ( !done   && ! thread_timeout(thread_pool[t]) continue;
	     // free the mem allocated at the pointer
	     free(thread_pool[t]);
	     // rearrange the pool array
	     memmove(thread_pool+t, thread_pool+t+1, (MAX_THREADS-t)*sizeof(thread_handler_arg) );
	     new_number_of_threads --;
	     }
	return 0;
    }



    int assign_thread (thread_handler_arg **thread_pool, int * number_of_threads) {

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
	thread_pool[number_of_threads] = thread_handler_arg_ptr;
	number_of_threads ++;
	return 0;
    
    }
