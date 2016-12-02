#include "igraph_daddy.h"

const char * SOCKET_PATH =  "/tmp/igraph_daddy";

pthread_mutex_t lock;

// send an apoptotic signal to a thread that has been running for too long
int thread_timeout(thread_args * thread_args_ptr) {
    if ( clock()-thread_args_ptr->start_time > TIME_MAX) {
	printf ("cancelling a thread\n");
	pthread_cancel(thread_args_ptr->pthread_id);
	return 1;
    }
    return 0;
}

// get rid of threads stuck in a too long calculation
void *killer_thread_handler(void  * threads_cast_to_void) {
    printf ("starting the garbage collector\n");
    while(1) {
	sleep(5);
	threadpool * threads = (threadpool *) threads_cast_to_void;
	pthread_mutex_lock(&lock);
	printf ("the beast is awake\n");
	int t;
	int new_number_of_threads = threads->number;
	for (t=MAX_THREADS-1; t>=0; t--) {
	    if (threads->args[t] == NULL) continue;
	    if (!thread_timeout(threads->args[t])) continue;
	    memmove(threads->args+t, threads->args+t+1, (MAX_THREADS-t-1)*sizeof(thread_args) );
	    threads->args[MAX_THREADS-1] = NULL;
	    new_number_of_threads --;
	}
	threads->number = new_number_of_threads;
	pthread_mutex_unlock(&lock);
    }
    return NULL;
}

// inidcate that more threads are available for work
int thread_pool_cleanup (threadpool * threads) {
    pthread_mutex_lock(&lock);
    printf ("the beast is awake\n");
    int t;
    int new_number_of_threads = threads->number;
    for (t=MAX_THREADS-1; t>=0; t--) {
	if (threads->args[t] == NULL) continue;
	int done = threads->args[t] -> job_done;
	if ( !done ) continue;
	memmove(threads->args+t, threads->args+t+1, (MAX_THREADS-t-1)*sizeof(thread_args) );
	threads->args[MAX_THREADS-1] = NULL;
	new_number_of_threads --;
    }
    threads->number = new_number_of_threads;
    pthread_mutex_unlock(&lock);
    return 0;
}

//
int assign_thread (threadpool * threads, int client_connection_id, igraph_t * graph_ptr) {

    pthread_t solver_thread;
    // thread_handler arg structure containts both the client_connection_id, 
    // and the pointer to the graph structure
    // we want to allocate it anew each time, and will let
    // the receiving thread free the mem - not ideal, but not sure
    // if it is worth further fussing
    thread_args * thread_args_ptr = emalloc (sizeof(thread_args));
    thread_args_ptr -> client_connection_id = client_connection_id;
    thread_args_ptr -> graph_ptr =  graph_ptr;
    thread_args_ptr -> start_time = clock();
    /* int pthread_create(pthread_t *thread, const pthread_attr_t *attr,
       void *(*start_routine) (void *), void *arg);
       The pthread_create() function starts a new thread in the calling
       process.  The new thread starts execution by invoking
       start_routine(); arg is passed as the sole argument of start_routine().
    */
    //pthread_attr_t attr; // thread attribute
    //pthread_attr_init(&attr); // use in place of hte seond argument of pthread_create
    // do not detach later - create as detachable - then I have to worry about freeing the attr structure
    //pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    int pthred_retval = pthread_create (&solver_thread, NULL, thread_handler, thread_args_ptr);
    if (pthred_retval<0) {
	free(thread_args_ptr);
	perror("could not create thread");
	return 1;
    }
    pthread_detach(solver_thread);
    thread_args_ptr -> pthread_id = solver_thread;
    
    pthread_mutex_lock(&lock);
    threads->args[ threads->number ] = thread_args_ptr;
    (threads->number) ++;
    pthread_mutex_unlock(&lock);
    
    return 0;
    
}

//////////////////////////////////////////////////
int threadpool_init(threadpool * threads) {
    threads->args = emalloc(MAX_THREADS*sizeof(thread_args));
    if (threads->args) {
	threads->number = 0;
	return 0;
    }
    fprintf (stderr, "Failed to allocate threadpool space\n");
    return 1;
}

int threadpool_shutdown(threadpool * threads) {
    if (threads && threads->args) free(threads->args);
    return 0;
}

//////////////////////////////////////////////////
int main ( int argc, char * argv[]) {
    // arguments: path to graph
    // the graph is just a list of from-to identifier pairs, one per line
    if (argc != 2) {
	fprintf (stderr, "Usage: %s  <path to graph file>\n",  argv[0]);
	exit(1);
    }
    // connect to the rest of the world
    int socket_id = create_socket_connection(SOCKET_PATH);
    
    // read in and reconstruct the graph
    igraph_t graph;
    construct_graph(argv[1], &graph);
    printf ("daddy's listening ... \n");

    //  initialize the locking device - to be used between the main and the killer thread
    if (pthread_mutex_init(&lock, NULL) != 0) {
        printf("\n mutex init failed\n");
        return 1;
    }
    
    // we maintain the pool max number of threads, so as not to swamp the system
    threadpool threads;
    threadpool_init(&threads);
    
    // I think this is suposed to be bad programming, but then iGraph
    // should not block, and yet it does
    // therefore, I have a killer thread, with the only
    // purpose of killing iGraph solver when it gets stuck chasing its tail
    pthread_t killer_thread_id;
    int pthred_retval = pthread_create (&killer_thread_id, NULL, (void*) killer_thread_handler, &threads);
    if (pthred_retval<0) {
	threadpool_shutdown(&threads);
	perror("could not create thread");
	return 1;
    }    
    
    // main loop acepting the connection and solving the task
    int number_of_threads = 0;
    while (1) {
	printf ("number of threads: %d\n", number_of_threads);
	//note: accept blocks, it won't return untill somebody calls
	int client_connection_id = accept(socket_id, NULL, NULL);
	printf ("connection: %d\n", client_connection_id);
	if (client_connection_id <  0) continue; // something wnet wrong with accepting the connection
	// wait for an available thread if needed
	if (threads.number) thread_pool_cleanup(&threads);	
	assign_thread (&threads, client_connection_id, &graph);
    }

    threadpool_shutdown(&threads);
    return 0;
}
