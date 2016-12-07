#include "igraph_daddy.h"


igraph_t graph;
fifo connection_queue;
pool threadpool;


//
int assign_thread (int client_connection_id) {

    if (threadpool.population >= MAX_THREADS) return 1;
    pthread_t solver_thread;
    // thread_handler arg structure containts both the client_connection_id, 
    // and the pointer to the graph structure
    // we want to allocate it anew each time, and will let
    // the receiving thread free the mem - not ideal, but not sure
    // if it is worth further fussing
    thread_args * thread_args_ptr = emalloc (sizeof(thread_args));
# ifdef VERBOSE
    printf ("main thread allocated  %ld bytes to  %p\n",  sizeof(thread_args), (void *)thread_args_ptr);
# endif
    thread_args_ptr -> client_connection_id = client_connection_id;
    thread_args_ptr -> start_time = clock();
    /* int pthread_create(pthread_t *thread, const pthread_attr_t *attr,
       void *(*start_routine) (void *), void *arg);
       The pthread_create() function starts a new thread in the calling
       process.  The new thread starts execution by invoking
       start_routine(); arg is passed as the sole argument of start_routine().
    */
    pthread_attr_t attr; // thread attribute
    pthread_attr_init(&attr); // use in place of hte seond argument of pthread_create
    // do not detach later - create as detachable - do  I have to worry about freeing the attr structure
    // i.e will the thing copy its attributes?
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    int pthred_retval = pthread_create (&solver_thread, &attr, thread_payload_function, thread_args_ptr);
    if (pthred_retval<0) {
    	free(thread_args_ptr);
    	perror("could not create thread");
    	return 1;
    }
    pthread_attr_destroy(&attr);
    
    int args_copied = 0;
    while (! args_copied) {
    	pthread_rwlock_rdlock(&(threadpool.args_lock));
    	args_copied = (thread_args_ptr->copied_to!=NULL);
    	pthread_rwlock_unlock(&(threadpool.args_lock));
    }
# ifdef VERBOSE
     printf ("main thread args copied  to  %p\n",   (void *)thread_args_ptr->copied_to);
# endif
    
    pthread_mutex_lock(&(threadpool.lock));
    threadpool.args[ threadpool.population ] =  (thread_args*) thread_args_ptr->copied_to;
    (threadpool.population) ++;
    pthread_mutex_unlock(&(threadpool.lock));
# ifdef VERBOSE
    // the thread should have copied its arguments
    printf ("number of threads: %d \n", threadpool.population); fflush(stdout);
    printf ("main thread freeing  %p\n", thread_args_ptr); fflush(stdout);
# endif
    free(thread_args_ptr);
    
    return 0;
    
}

//////////////////////////////////////////////////
int  dequeue () {
    int conn_id;
    pthread_mutex_lock (&(connection_queue.lock));
    if (connection_queue.population ==0 ) {
	conn_id =  -1;
    } else {
	conn_id  = connection_queue.position[0];
	connection_queue.population--;
	memmove (connection_queue.position, connection_queue.position+1,
		 connection_queue.population*sizeof(int));
    }
    pthread_mutex_unlock (&(connection_queue.lock));

    return conn_id;
}


//////////////////////////////////////////////////
// send an apoptotic signal to a thread that has been running for too long
int thread_timeout(thread_args * thread_args_ptr) {
    if ( clock()-thread_args_ptr->start_time > TIME_MAX) {
# ifdef VERBOSE
	printf ("cancelling  thread %d:  %d  limit  %d\n",
		thread_args_ptr->client_connection_id, (int) (clock()-thread_args_ptr->start_time), TIME_MAX);
# endif
	pthread_cancel(thread_args_ptr->pthread_id);
	return 1;
    }
    return 0;
}

//////////////////////////////////////////////////
int thread_pool_cleanup () {
# ifdef VERBOSE    
    printf ("in thread_pool_cleanup: pop size = %d\n", threadpool.population);
# endif
    if (threadpool.population==0) return 0;
    pthread_mutex_lock(&(threadpool.lock));
    int t;
    int new_number_of_threads = threadpool.population;
    for (t=MAX_THREADS-1; t>=0; t--) {
	if (threadpool.args[t] == NULL) continue;
	int done = threadpool.args[t]->job_done;
	if ( !done  && !thread_timeout(threadpool.args[t])) continue;
	memmove (threadpool.args+t, threadpool.args+t+1, (MAX_THREADS-t-1)*sizeof(thread_args*) );
	threadpool.args[MAX_THREADS-1] = NULL;
	new_number_of_threads --;
    }
    threadpool.population = new_number_of_threads;
    pthread_mutex_unlock(&(threadpool.lock));

    return 0;
}
//////////////////////////////////////////////////
void flush_queue () {
    int conn_id;
    while ( (conn_id=dequeue()) >= 0) {
	// give chance to threadpool to free some threads
	while (assign_thread (conn_id)) {
	    sleep(1);
	    thread_pool_cleanup ();
	}
    }
}
    


//////////////////////////////////////////////////
// systemctl stop is sending SIGTERM
// SIGKILL apparently cannot be  intercepted
int daddy_shutdown = 0;
void signal_handler(int signo) {
    daddy_shutdown = 1;
}

//////////////////////////////////////////////////
void  send_and_close (int client_socket_id,  char * message) {
    char * output_buffer = emalloc(100);
    strncpy (output_buffer, message, 99);
    int size  =  strlen(output_buffer)+1;
    if (send (client_socket_id, output_buffer, size, MSG_NOSIGNAL) < 0 ) {
 	perror ("write error");
    }
    free(output_buffer);
    return;
}

//////////////////////////////////////////////////
void queue (int conn_id) {
    pthread_mutex_lock (&(connection_queue.lock));
    if (connection_queue.population < CONNECTION_QUEUE_SIZE) {
	connection_queue.position[connection_queue.population] = conn_id;
	connection_queue.population ++;
    } else {
	send_and_close (conn_id, "Connection dropped");
    }
    pthread_mutex_unlock (&(connection_queue.lock));
}

//////////////////////////////////////////////////
void * fifo_controller () {
  
    connection_queue.position   = emalloc(CONNECTION_QUEUE_SIZE*sizeof(int));
    connection_queue.population = 0;
    pthread_mutex_init  (&(connection_queue.lock), NULL);
    threadpool.args       = emalloc(MAX_THREADS*sizeof(thread_args*));
    threadpool.population = 0;
    pthread_mutex_init(&(threadpool.lock), NULL);
    pthread_rwlock_init (&(threadpool.args_lock), NULL);
    
    while (!daddy_shutdown) {
	sleep(20);
	thread_pool_cleanup();
	flush_queue();
    }
    
    pthread_rwlock_destroy (&(threadpool.args_lock));
    pthread_mutex_destroy  (&(threadpool.lock));
    pthread_mutex_destroy  (&(connection_queue.lock));
    free (connection_queue.position);
    free (threadpool.args);
    
    return NULL;
}

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

//////////////////////////////////////////////////
int main ( int argc, char * argv[]) {
    // arguments: path to graph
    // the graph is just a list of from-to identifier pairs, one per line
    if (argc != 3) {
	fprintf (stderr, "Usage: %s  <socket path>  <path to graph file>\n",  argv[0]);
	exit(1);
    }
    // connect to the rest of the world
    int socket_id = create_socket_connection(argv[1]);    
    // read in and reconstruct the graph
    construct_graph(argv[2], &graph); // graph globally declared
    // handle termination gracefully
    signal(SIGTERM, signal_handler);
    printf ("daddy's listening ... \n");

    // we maintain the pool max number of threads, so as not to swamp the system
    // I think this is suposed to be bad programming, but then iGraph
    // should not block, and yet it does
    // therefore, I have a killer thread, with the only
    // purpose of killing iGraph solver when it gets stuck chasing its tail
    pthread_t fifo_controller_thread_id;
    int pthred_retval = pthread_create (&fifo_controller_thread_id, NULL, (void*) fifo_controller, NULL);
    if (pthred_retval<0) {
	perror("could not create fifo controller thread");
	return 1;
    }    
    
    // main loop acepting the connection and solving the task
    while (!daddy_shutdown) {
 	int client_connection_id = accept(socket_id, NULL, NULL);
	if (client_connection_id <  0) continue; // something went wrong with accepting the connection
 	// post  connection id 
	queue (client_connection_id); // if the queue is full, th connection will be dropped
	thread_pool_cleanup();
	flush_queue();
   }
    // signal to fifo controller
    // wait for the queue to drain and fifo_thread controller to join
    daddy_shutdown = 1;
    pthread_join(fifo_controller_thread_id, NULL);
    return 0;
}
