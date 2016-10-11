#ifndef _igraph_daddy_h
#define _igraph_daddy_h 

#include<stdio.h>
#include<string.h>    //strlen
#include<stdlib.h>    //strlen
#include<sys/socket.h>
#include<sys/un.h>
#include<arpa/inet.h> //inet_addr
#include<unistd.h>    //write
#include<pthread.h>   //for threading , link with lpthread
#include<sys/signal.h>
#include<igraph.h>


// hardcoded
extern char *SOCKET_PATH;
extern int BUF_BLOCK;
extern int curr_input_buf_size;

//thread handler arguments - must be wrapped in a single structure
typedef struct {
    // connection id
    int client_connection_id;
    // graph pointer
    void * graph_ptr;
    //graph * graph_ptr;
}  thread_handler_arg;

typedef enum {NEIGHBORS, PATH} igraph_method;

typedef struct {
    igraph_method ig_method;
    int order;
    int * node_list;
    int node_list_length; 
} igraph_args;


// utility functions
FILE * efopen(char * name, char * mode);
void * emalloc(int size);
void increase_buf_size (char ** buffer, char ** current_write_pos);

// threading related functions
int  create_socket_connection();
void *thread_handler(void *inptr);

// igraph-related functions
void construct_graph(char *filename, igraph_t *graph_ptr);
int input_parser (char *input_buffer, igraph_arg *ig_args);
int solver( igraph_arg *ig_args, char * output_buffer);
#endif
