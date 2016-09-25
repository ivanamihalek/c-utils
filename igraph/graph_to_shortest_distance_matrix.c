#include <igraph.h>
/************************************************/
FILE* efopen(char * name, char * mode);
int print_matrix(FILE* fptr, const igraph_matrix_t *m);

/************************************************/
int main(int argc, char*argv[]) {

    // the graph is just a list og from to identifier pairs, one per line
    if (argc != 3) {
	fprintf (stderr, "Usage: %s  <path to graph>  <path to matrix (output)>\n", argv[0]);
	exit(1);
    }
      
    igraph_t graph;
    /*The number of vertices in the graph. If smaller than the largest integer in
      the file it will be ignored. It is thus safe to supply zero here. */
    igraph_integer_t zero = 0;
    FILE * infile  =  efopen(argv[1], "r") ;
    FILE * outfile =  efopen(argv[2], "w") ;
    
    igraph_read_graph_edgelist (&graph, infile,  zero, IGRAPH_UNDIRECTED);
    fclose(infile);
    // vertex selector - initialization
    //igraph_vs_t vs; //A vertex selector object.
    //igraph_vit_t vit;//an uninitialized vertex iterator object. 
    //igraph_vit_create (&graph,  vs, &vit);
    /*int igraph_independent_vertex_sets(const igraph_t *graph,
				   igraph_vector_ptr_t *res,
				   igraph_integer_t min_size,
				   igraph_integer_t max_size);
      min_size:
             Integer giving the minimum size of the sets to be returned.
	     If negative or zero, no lower bound will be used.

      max_size:
	     Integer giving the maximum size of the sets to be returned.
	     If negative or zero, no upper bound will be used.
    */
    igraph_vector_ptr_t  result;
    long int n;
    igraph_vector_ptr_init (&result, 0);
    igraph_independent_vertex_sets(&graph, &result, 0, 0);
    n = igraph_vector_ptr_size(&result);
    printf("%ld independent sets found\n", (long)n);
  
#if 0
    igraph_matrix_t result;
    igraph_matrix_init(&result, zero, zero);

    igraph_shortest_paths (&graph, &result, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
    printf ("output ...\n");
    print_matrix(outfile, &result);

    igraph_matrix_destroy(&result);
#endif
    
    fclose(outfile);
    igraph_destroy(&graph);
    if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;
  
    return 0;
}

/************************************************/
FILE* efopen(char * name, char * mode) {

    FILE * fp;
    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,    "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }
    return fp;
}
/************************************************/
int print_matrix(FILE* fptr, const igraph_matrix_t *m) {
    long int nrow=igraph_matrix_nrow(m);
    long int ncol=igraph_matrix_ncol(m);
    long int i, j;
    igraph_real_t val;

    for (i=0; i<nrow; i++) {
	fprintf (fptr, "%li:", i);
	for (j=0; j<ncol; j++) {
	    val = MATRIX(*m, i, j);
	    if (igraph_is_inf(val)) {
		if (val < 0) {
		    fprintf (fptr, "-inf");
		} else {
		    fprintf (fptr, " inf");
		}
	    } else {
		fprintf (fptr, " %3.0f", val);
	    }
	}
	fprintf (fptr, "\n");
    }
    return 0;
}

