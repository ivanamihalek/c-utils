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
#include <sys/signal.h>

// hardcoded
extern char *SOCKET_PATH;
extern int BUF_BLOCK;
extern int curr_input_buf_size;

//thread handler arguments - must be wrapped in a single structure
typedef struct _thread_handler_arg {
    // connection id
    int client_connection_id;
    // graph pointer
    void * graph_ptr;
    //graph * graph_ptr;
}  thread_handler_arg;

// utility functions
void * emalloc(int size);
void increase_buf_size (char ** buffer, char ** current_write_pos);

// threading related functions
int  create_socket_connection();
void *thread_handler(void *inptr);

#endif
