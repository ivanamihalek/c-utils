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
 int vector_ptr2vs_vector (igraph_vs_t *vs, igraph_vector_ptr_t result){
    long int i, j, k, size,  n = igraph_vector_ptr_size(&result) ;
    igraph_vector_t *v;
	
    vs->type=IGRAPH_VS_VECTOR;
    vs->data.vecptr=igraph_Calloc(1, igraph_vector_t);
    if (vs->data.vecptr==0) {
	IGRAPH_ERROR("Cannot create vertex selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)vs->data.vecptr);
    k = 0;
    size = 2; // for the first and the last vertex
    for (i=0; i<n; i++) {
	v = VECTOR(result)[i];
	long int j, m=igraph_vector_size(v);
	for (j=1; j<m-1; j++) size++;
    }
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)vs->data.vecptr,size);

    k = 0;
    v = VECTOR(result)[0];
    // the frist and the lst vertex are the same for all paths,
    // that's why we store them only once
    VECTOR(*vs->data.vecptr)[k] = (long int) VECTOR(*v)[0];  k++;
    // I should also be removing duplicates here ...
    for (i=0; i<n; i++) {
	v = VECTOR(result)[i];
	long int j, m=igraph_vector_size(v);
	for (j=1; j<m-1; j++) {
	    VECTOR(*vs->data.vecptr)[k] = (long int)VECTOR(*v)[j]; k++;
	}
	if (i==n-1){
	    j = m-1;
	    VECTOR(*vs->data.vecptr)[k] = (long int)VECTOR(*v)[j]; k++;
	}
    }


    IGRAPH_FINALLY_CLEAN(2);
    return 0;  
}

/************************************************/
int main(int argc, char*argv[]) {

    // the graph is just a list og from to identifier pairs, one per line
    if (argc != 5) {
	fprintf (stderr, "Usage: %s  <path to graph>   <from id>  <to id>  <order> \n", argv[0]);
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
    
    igraph_integer_t order = atoi(argv[4]);
   
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
    igraph_integer_t number_of_paths = igraph_vector_ptr_size(&result);
    if (number_of_paths==0) {
	printf ("no paths found\n");
	exit(0);
    }
    // I'm not sure what's with all the idioutc formtats these
    // guys are defining; at any rate, turn vector_ptr to vs_vector
    // because that is what we need to find neighbors
    igraph_vs_t  vids; 
    vector_ptr2vs_vector (&vids, result);
    // clean up the old junk
    igraph_vs_destroy(&vertex_to);
    igraph_vector_ptr_destroy(&result);
    // create new result vector
    igraph_vector_ptr_init (&result, number_of_paths);
    // finally find all nbrs on the path
    igraph_neighborhood (&graph, &result, vids,  order, IGRAPH_ALL);
    for (i=0; i<igraph_vector_ptr_size(&result); i++) {
	print_vector(VECTOR(result)[i]);
    }

    igraph_vector_ptr_destroy(&result);
    igraph_destroy(&graph);

    if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;
  
    return 0;
}
