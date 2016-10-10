#include<stdio.h>
#include<string.h>    //strlen
#include<stdlib.h>    //strlen
#include<sys/socket.h>
#include<sys/un.h>
#include<arpa/inet.h> //inet_addr
#include<unistd.h>    //write
#include<pthread.h>   //for threading , link with lpthread
#include <sys/signal.h>

// hardcoded
char *SOCKET_PATH = "\0igraph_daddy";
int BUF_SIZE = 4;

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
    tmp_buffer_ptr = emalloc(sizeof(*buffer) + BUF_SIZE);
    memcpy (tmp_buffer_ptr, *buffer, sizeof(*buffer));
    *current_write_pos =  tmp_buffer_ptr +  sizeof(*buffer);
    free (*buffer);
    *buffer = tmp_buffer_ptr;
    return;
}

// arguments parser

// connection creator
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
// graph reader

//the thread function
typedef struct _solver_arg {
    // connection id
    int client_connection_id;
    // graph pointer
    void * graph_ptr;
    //graph * graph_ptr;
}  solver_arg;

////////////////////////////////////

void *solver(void *inptr) {
    printf ("Hello from solver\n");
    solver_arg *solver_arg_ptr = (solver_arg *) inptr;
    // read from connection
    char *input_buffer  = emalloc (BUF_SIZE); // will initialize to 0
    char *current_write_pos = input_buffer;
    char *output_buffer = emalloc (2*BUF_SIZE);
    // First, read the length of the text message from the socket.
    // If read returns zero, the client closed the connection.
    int client_socket_id = solver_arg_ptr->client_connection_id;
    int retval;
    printf ("  4  starting the read\n");  fflush(stdout);
    int done = 0;
    int panic_ctr = 0;
    while (!done) {
	if ( retval = recv(client_socket_id, current_write_pos, BUF_SIZE, 0) <  0) {
	    perror ("error reading");
	    return NULL;
	}
	panic_ctr ++;
	done = (strchr (input_buffer, '\n') != NULL) || (panic_ctr > 10); 
	if(!done) increase_buf_size (&input_buffer, &current_write_pos);
	
    }
    printf ("*** input: %s\n", input_buffer);
    // parse the input
    // solve
    // write to connection
    sprintf(output_buffer, "%s  pong: %d ", input_buffer, getpid());
    int size  =  strlen(output_buffer)+1;
    sigignore(SIGPIPE); // we don't want to die here if the client quit
    printf ("going send\n");
    if ( retval = send (client_socket_id, output_buffer, size, 0) < 0 ) {
 	perror ("write error");
	return NULL;
    }
    printf ("back from send\n");
    printf ("send retval: %d\n", retval); 
    
    free (input_buffer);
    free (output_buffer);
    free (solver_arg_ptr);
    return NULL;
}

//////////////////////////////////////////////////
int main ( int argc, char * argv[]) {
    // arguments: path to graph

    // create socket connection
    int socket_id = create_socket_connection();
    // read in and reconstruct the graph
    void * graph_ptr;
    // main lopp acepting the connection and solving the task
    
    while (1) {
	int client_connection_id = accept(socket_id, NULL, NULL);
	if (client_connection_id < 0)  continue; // accept error
	pthread_t solver_thread;
	// solver arg structure containts both the client_connection_id, 
	// and the pointer to the graph structure
	// we want to allocate it anew each time, and will let
	// the receiving thread free the mem - not ideal, but not sure
	// if it si worth futher fussing
	solver_arg * solver_arg_ptr = emalloc (sizeof(solver_arg));
	solver_arg_ptr -> client_connection_id = client_connection_id;
	solver_arg_ptr -> graph_ptr = graph_ptr;
	/* int pthread_create(pthread_t *thread, const pthread_attr_t *attr,
                          void *(*start_routine) (void *), void *arg);
	  The pthread_create() function starts a new thread in the calling
          process.  The new thread starts execution by invoking
          start_routine(); arg is passed as the sole argument of start_routine().
	*/
	printf ("... creating thread ...\n");
	int pthred_retval = pthread_create (&solver_thread, NULL, solver, solver_arg_ptr);
	if (pthred_retval<0) {
	    perror("could not create thread");
	    return 1;
	}
    }
    return 0;
}
