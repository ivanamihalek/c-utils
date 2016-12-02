#ifndef _igraph_daddy_h
#define _igraph_daddy_h 

#include<stdio.h>
#include <ctype.h>  // isspace
#include<string.h>  //strlen, memcpy
#include<stdlib.h>  
#include<sys/socket.h>
#include<sys/un.h>
#include<arpa/inet.h> //inet_addr
#include<unistd.h>    //write
#include<pthread.h>   //for threading , link with lpthread
#include<sys/signal.h>
#include<igraph.h>
#include <time.h>
// input tokenizer array
#define MAX_TOK 1010
#define TOKENLENGTH 30
#define BUF_BLOCK 1024

# define MAX_THREADS   10
//clock_t TIME_MAX = 60*CLOCKS_PER_SEC; //I have to make sure that the requests are small enough
# define  TIME_MAX  100 //I have to make sure that the requests are small enough



//thread handler arguments - must be wrapped in a single structure
typedef struct {
    // ids
    pthread_t pthread_id;
    int client_connection_id;
    // graph pointer
    void * graph_ptr;
    // flag to signal job competion
    int job_done;
    // the job will be timed out if needed:
    clock_t start_time;
} thread_args;

typedef struct {
    thread_args ** args;
    int number;
} threadpool;




typedef enum {UNK, NEIGHBORS, PATH} igraph_method;

typedef struct {
    igraph_method ig_method;
    int order;
    long int * node_list;
    int node_list_length;
    void * graph_ptr;
    char * error_msg;
} igraph_arg;

// utility functions
FILE * efopen(char * name, char * mode);
void * emalloc(int size);
void increase_buf_size (char ** buffer, int*curr_input_buf_size_ptr,  char ** current_write_pos);
int contains_nondigit (char * string);
int tokenize (char token[MAX_TOK][TOKENLENGTH], int * max_token, char * line , char comment_char);

// threading related functions
int  create_socket_connection(const char * socket_path);
void *thread_handler(void *inptr);

// igraph-related functions
void construct_graph (char *filename, igraph_t *graph_ptr);
int input_parser (char *input_buffer, igraph_arg *ig_args);
int solver (igraph_arg *ig_args, char ** output_buffer_ptr);

#endif
