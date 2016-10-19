#include <igraph.h>


/************************************************/
FILE * efopen(char * name, char * mode) {

    FILE * fp;
    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,    "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }
    return fp;
}
/************************************************/
void print_vector(igraph_vector_t *v) {
    long int i, l=igraph_vector_size(v);
    for (i=0; i<l; i++)  printf("  %li", (long int) VECTOR(*v)[i]);
    printf("\n");
}


/************************************************/
int main(int argc, char*argv[]) {

    // the graph is just a list og from to identifier pairs, one per line
    if (argc != 2) {
	fprintf (stderr, "Usage: %s  <path to graph>\n", argv[0]);
	exit(1);
    }
    /* turn on attribute handling */
    igraph_i_set_attribute_table(&igraph_cattribute_table);

 
    igraph_t graph;
    printf ("number of vertices:  %d\n", igraph_vcount(&graph));
    /*The number of vertices in the graph. If smaller than the largest integer in
      the file it will be ignored. It is thus safe to supply zero here. */
    igraph_integer_t zero = 0;

    FILE * infile  =  efopen(argv[1], "r") ;
    
    igraph_bool_t directed  = IGRAPH_UNDIRECTED;
    igraph_read_graph_edgelist (&graph, infile,  zero, directed);
    printf ("number of vertices:  %d\n", igraph_vcount(&graph));
    
    /*int igraph_decompose(const igraph_t *graph, igraph_vector_ptr_t *components, 
		     igraph_connectedness_t mode,
		     long int maxcompno, long int minelements);
    */
    igraph_vector_ptr_t components;
    igraph_vector_ptr_init (&components, 0);

    igraph_decompose (&graph, &components, IGRAPH_WEAK, -1, 2);
    
    int number_of_components  =  igraph_vector_ptr_size(&components);
    printf ("number of components:  %d\n", number_of_components);
        
    int i;
    for (i=0; i<igraph_vector_ptr_size(&components); i++) { 
	igraph_t *comp=VECTOR(components)[i];
    }
    /*Don't forget to call igraph_destroy() and free() on the elements of
      this pointer vector to free unneeded memory. Alternatively, you can simply
      call igraph_decompose_destroy() that does this for you. 
    */
    
    igraph_vector_ptr_destroy(&components);
    igraph_destroy(&graph);

    if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;
  
    return 0;
}
