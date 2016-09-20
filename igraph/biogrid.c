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
int main(void) {
    
    igraph_t graph;
    /*The number of vertices in the graph. If smaller than the largest integer in
      the file it will be ignored. It is thus safe to supply zero here. */
    igraph_integer_t zero = 0;

    FILE * infile  =  efopen("/mnt/databases/biogrid/biogrid_homo.uniq.txt", "r") ;
    
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
    igraph_integer_t vertex_from;
    igraph_vs_t  vertex_to; 
    igraph_integer_t i;
 
    // 57572 DOCK6
    // 7049 TGFBR3
    // 2263 FGFBP3
    vertex_from = 57572; 
    igraph_vs_vector_small(&vertex_to, 2263, 7049,  -1);

    igraph_vector_ptr_init (&result, 2);
    igraph_vector_init(&nrgeo, 0);
    igraph_get_all_shortest_paths (&graph, &result, &nrgeo,
			    vertex_from,  vertex_to,
			   IGRAPH_ALL);
    printf ("number of geodesics 0  %8.3lf\n", VECTOR(nrgeo)[0]);
    for (i=0; i<igraph_vector_ptr_size(&result); i++) {
	printf ("**************%7d   %3d ********************\n", vertex_from, i);
	print_vector(VECTOR(result)[i]);
    }

    igraph_vector_ptr_destroy(&result);
    igraph_vs_destroy(&vertex_to);
    igraph_destroy(&graph);

    if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;
  
    return 0;
}
