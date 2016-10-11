#include "igraph_daddy.h"

#define MAX_ORDER 5 // no more than 5 hops away
void construct_graph(char *filename, igraph_t *graph_ptr){

    FILE * infile  =  efopen(filename, "r") ;
    if (!infile) exit(1); // efopen outputs the error msg
    igraph_read_graph_edgelist (graph_ptr, infile,  0, IGRAPH_UNDIRECTED);
    fclose(infile);

    return;
}

int input_parser (char *input_buffer, igraph_arg *ig_args) {

    char token[MAX_TOK][TOKENLENGTH];
    int max_token;
    
    tokenize(token, &max_token, input_buffer, '#');
    if (strcmp(token[0], "neighbors")==0) {
	ig_args->ig_method =  NEIGHBORS;
    } else  if (strcmp(token[0], "path")==0) {
	ig_args->ig_method =  PATH;	
    } else {
	return 1
    }
    // max token is the last index - the token list is  max_token + 1
    // the first two tokens are method name and 'order' - the number of hops
    ig_args->node_list_length = max_token + 1 - 2 ;
    if (ig_args->ig_method ==  NEIGHBORS && ig_args->node_list_length < 1) {
	return 1;
    }
    if (ig_args->ig_method ==  PATH && ig_args->node_list_length < 2) {
	return 1;
    }
    
    ig_args->order = atoi(token[1]);
    if ( ig_args->order < 0  ||  ig_args->order>MAX_ORDER)  return 1;

    ig_args->node_list = emalloc(ig_args->node_list_length * sizeof(int));

    int i;
    for(i=0; i<ig_args->node_list_length; i++) {
	ig_args->node_list[i] = atoi(token[2+i]);	
    }
    
    return 0;
}
////////////////////////////////////////////////////////

int solver( igraph_arg *ig_args, char * output_buffer) {

    printf ("in solver: \n");
    if (ig_args->ig_method == NEIGHBORS) {
	 printf ("\t method: neighors \n");
    } else if (ig_args->ig_method == PATH) {
	 printf ("\t method: path \n");
    } else { // we shouldn't be here ...
	return 1;
    }
    
    return 0;
}
