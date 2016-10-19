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
long int find_vector_by_name(igraph_t * g){
    /*** do we have the name attribute ? ***/
    // could not get this to work
    /* printf ("  have names?  %d\n", 
	    (int) igraph_cattribute_has_attr (g, IGRAPH_ATTRIBUTE_VERTEX, "name"));
    */

   

    
    return 0;
}


/************************************************/
void print_vector(igraph_vector_t *v) {
    long int i, l=igraph_vector_size(v);
    for (i=0; i<l; i++)  sprintf("  %li", (long int) VECTOR(*v)[i]);
    printf("\n");
}


/************************************************/
int args2vs_vector (igraph_vs_t *vs, char *argv[], int arg_from, int arg_to) {
    long int i, n = arg_to-arg_from+1 ;
    vs->type=IGRAPH_VS_VECTOR;
    vs->data.vecptr=igraph_Calloc(1, igraph_vector_t);
    if (vs->data.vecptr==0) {
	IGRAPH_ERROR("Cannot create vertex selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)vs->data.vecptr);
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)vs->data.vecptr,n);
    
    for (i=0; i<n; i++) {
	VECTOR(*vs->data.vecptr)[i]=atol(argv[arg_from+i]);
    }
    IGRAPH_FINALLY_CLEAN(2);
    return 0;  
}    

/************************************************/
int main(int argc, char*argv[]) {

    // the graph is just a list of from-to identifier pairs, one per line
    if (argc < 4) {
	fprintf (stderr, "Usage: %s  <path to graph file> <order> <id1>  [<id2> <id3>  ... ]\n",
		 argv[0]);
	exit(1);
    }
    FILE * infile  =  efopen(argv[1], "r") ;
    igraph_t graph;
    int retval  =  igraph_read_graph_ncol (&graph, infile, NULL, 1,
					   IGRAPH_ADD_WEIGHTS_NO,  IGRAPH_UNDIRECTED);
    fclose(infile);
    if (retval == IGRAPH_PARSEERROR) {
	fprintf (stderr, "Error reading %s\n", argv[1]);
    }
    printf ("number of vertices:  %d\n",  igraph_vcount(&graph));
    /* setup iterator over all vertices: selector */
    igraph_vs_t all_vertices_selector;
    igraph_vit_t vertex_iterator;
    igraph_vs_all  (&all_vertices_selector);
    igraph_integer_t size;
    igraph_vs_size (&graph, &all_vertices_selector, &size);
    printf ("%li ", (long int) size);
  
    //igraph_vit_create (&graph, vertices, &vertex_iterator);
    
    //find_vector_by_name (&graph);
    exit(1);

    int order = atoi(argv[2]);
    // 57572 DOCK6
    // 7049 TGFBR3
    // 2263 FGFBP3
    igraph_vs_t  vids; 
    args2vs_vector (&vids, argv, 3, argc-1);
    
     fclose(infile);
    //printf ("done reading in the graph\n"); fflush(stdout);
    /* neighborhood calculation: */
    /* int igraph_neighborhood(const igraph_t *graph, igraph_vector_ptr_t *res,
			igraph_vs_t vids, igraph_integer_t order,
			igraph_neimode_t mode);
      vids: The vertices for which the calculation is performed. 
    */
    igraph_vector_ptr_t result;
    // the second argument here is the number of "seed"  vertices
    igraph_vector_ptr_init (&result, argc-3);
    igraph_neighborhood (&graph, &result, vids,  order, IGRAPH_ALL);
    int i;
    for (i=0; i<igraph_vector_ptr_size(&result); i++) {
	print_vector(VECTOR(result)[i]);
    }
    
    igraph_vit_destroy (&vertex_iterator);
    igraph_vector_ptr_destroy(&result);
    igraph_destroy(&graph);

    if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;
  
    return 0;
}
