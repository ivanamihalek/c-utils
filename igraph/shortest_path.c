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
    if (argc != 4) {
	fprintf (stderr, "Usage: %s  <path to graph>  <from id>  <to id>\n", argv[0]);
	exit(1);
    }
    
    // 57572 DOCK6
    // 7049 TGFBR3
    // 2263 FGFBP3
    igraph_vs_t  vertex_to; 
    igraph_integer_t vertex_from;
    vertex_from = atol(argv[2]);
    // in principle, unumber of "to" vertices is arbitrary,
    // but that is a piece of code I will write some other time
    //igraph_vs_vector_small(&vertex_to, 2263, 7049,  -1);
    igraph_vs_vector_small(&vertex_to, atol(argv[3]),  -1);
   
    igraph_t graph;
    /*The number of vertices in the graph. If smaller than the largest integer in
      the file it will be ignored. It is thus safe to supply zero here. */
    igraph_integer_t zero = 0;

    FILE * infile  =  efopen(argv[1], "r") ;
    
    igraph_bool_t directed  = IGRAPH_UNDIRECTED;
    igraph_read_graph_edgelist (&graph, infile,  zero, directed);

    /* shortest path calculation: */
    /* http://igraph.org/c/doc/igraph-Structural.html#igraph_get_all_shortest_paths
       int igraph_get_all_shortest_paths(const igraph_t *graph,
				  igraph_vector_ptr_t *res, 
				  igraph_vector_t *nrgeo,
				  igraph_integer_t from, const igraph_vs_t to,
				  igraph_neimode_t mode);
    */
    igraph_vector_ptr_t result;
    igraph_vector_t nrgeo;
    igraph_integer_t i;
    // the second argument here is the number of "to"  vertices
    igraph_vector_ptr_init (&result, 1);
    igraph_vector_init(&nrgeo, 0);
    igraph_get_all_shortest_paths (&graph, &result, &nrgeo,
			    vertex_from,  vertex_to,
			   IGRAPH_ALL);
    for (i=0; i<igraph_vector_ptr_size(&result); i++) {
	print_vector(VECTOR(result)[i]);
    }

    igraph_vector_ptr_destroy(&result);
    igraph_vs_destroy(&vertex_to);
    igraph_destroy(&graph);

    if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;
  
    return 0;
}
