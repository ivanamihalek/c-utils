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
    
    if (tokenize(token, &max_token, input_buffer, '#')) return 1;
    if (strcmp(token[0], "neighbors")==0) {
	ig_args->ig_method = NEIGHBORS;
    } else  if (strcmp(token[0], "path")==0) {
	ig_args->ig_method = PATH;	
    } else {
	return 1;
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

    ig_args->node_list = emalloc(ig_args->node_list_length * sizeof(long int));

    int i;
    for(i=0; i<ig_args->node_list_length; i++) {
	ig_args->node_list[i] = atol(token[2+i]);	
    }
    
    return 0;
}
////////////////////////////////////////////////////////
int node_list2vs_vector (igraph_vs_t *vs, long int * node_list, int node_list_length) {
    long int i, n = node_list_length ;
    vs->type=IGRAPH_VS_VECTOR;
    vs->data.vecptr=igraph_Calloc(1, igraph_vector_t);
    if (vs->data.vecptr==0) {
	IGRAPH_ERROR("Cannot create vertex selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)vs->data.vecptr);
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)vs->data.vecptr,n);
    
    for (i=0; i<n; i++) {
	VECTOR(*vs->data.vecptr)[i] = node_list[i];
    }
    IGRAPH_FINALLY_CLEAN(2);
    return 0;  
}
///
void print_vector (char ** outbuf_ptr, igraph_vector_ptr_t *result_ptr) {
    igraph_vector_t *v;    
    long int j, i, l;
    int n_res = igraph_vector_ptr_size(result_ptr);
    int alloc_size = 0;
    for (j=0; j<n_res; j++) {
        v = VECTOR(*result_ptr)[j];
	l = igraph_vector_size(v);
	alloc_size += l*TOKENLENGTH;
    }
    
    printf ("allocating %d  \n", alloc_size);
    char * outbuf = emalloc(alloc_size); // can it ever have more than this many digits?
    for (j=0; j<n_res; j++) {
        v = VECTOR(*result_ptr)[j];
	l = igraph_vector_size(v);
	for (i=0; i<l; i++)
	    sprintf(outbuf + strlen(outbuf), "  %li", (long int) VECTOR(*v)[i]);
	printf("\n");
    }
    *outbuf_ptr = outbuf;
}
////
int neighbors (igraph_arg *ig_args, char ** output_buffer_ptr) {
    
    igraph_t * graph_ptr   = ig_args->graph_ptr;
    igraph_integer_t order = ig_args->order;
    igraph_vs_t  vids;
    node_list2vs_vector (&vids, ig_args->node_list, ig_args->node_list_length);
    igraph_vector_ptr_t result;
    // the second argument here is the number of "seed"  vertices
    igraph_vector_ptr_init (&result, ig_args->node_list_length);
    igraph_neighborhood (graph_ptr, &result, vids,  order, IGRAPH_ALL);
    /// output
    print_vector (output_buffer_ptr, &result);    
    
    igraph_vector_ptr_destroy(&result);

    return 0;
}
///
int solver( igraph_arg *ig_args, char ** output_buffer_ptr) {

    printf ("in solver: \n");
    if (ig_args->ig_method == NEIGHBORS) {
	 printf ("\t method: neighbors \n");
	 return neighbors(ig_args,output_buffer_ptr);
    } else if (ig_args->ig_method == PATH) {
	 printf ("\t method: path \n");
    } else { // we shouldn't be here ...
	return 1;
    }
    return 0;
}
//////
