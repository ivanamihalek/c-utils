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


int print_attributes(const igraph_t *g) {

    igraph_vector_t gtypes, vtypes, etypes;
    igraph_strvector_t gnames, vnames, enames;
    long int i;

    igraph_vector_t vec;
    igraph_strvector_t svec;
    long int j;

    igraph_vector_init(&gtypes, 0);
    igraph_vector_init(&vtypes, 0);
    igraph_vector_init(&etypes, 0);
    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);

    igraph_cattribute_list(g, &gnames, &gtypes, &vnames, &vtypes, 
			   &enames, &etypes);

    /* vertex attributes */
    for (i=0; i<igraph_vcount(g); i++) {
	long int j;
	printf("Vertex %li: ", i);
	for (j=0; j<igraph_strvector_size(&vnames); j++) {
	    printf("%s=", STR(vnames, j));
	    if (VECTOR(vtypes)[j]==IGRAPH_ATTRIBUTE_NUMERIC) {
		igraph_real_printf(VAN(g, STR(vnames,j), i));
		putchar(' ');
	    } else {
		printf("\"%s\" ", VAS(g, STR(vnames,j), i));
	    }
	}
	printf("\n");
    }



    /* Check vector-based query functions 
       igraph_vector_init(&vec, 0);
       igraph_strvector_init(&svec, 0);
  
       for (j=0; j<igraph_strvector_size(&vnames); j++) {
       if (VECTOR(vtypes)[j]==IGRAPH_ATTRIBUTE_NUMERIC) {
       igraph_cattribute_VANV(g, STR(vnames, j), igraph_vss_all(), &vec);
       for (i=0; i<igraph_vcount(g); i++) {
       igraph_real_t num=VAN(g, STR(vnames, j), i);
       if (num != VECTOR(vec)[i] &&
       (!isnan(num) || !isnan(VECTOR(vec)[i]))) {
       exit(51);
       }
       }
       } else {
       igraph_cattribute_VASV(g, STR(vnames, j), igraph_vss_all(), &svec);
       for (i=0; i<igraph_vcount(g); i++) {
       const char *str=VAS(g, STR(vnames, j), i);
       if (strcmp(str,STR(svec, i))) {
       exit(52);
       }
       }
       }
       }
    */
 

    igraph_strvector_destroy(&svec);
    igraph_vector_destroy(&vec);

    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    igraph_vector_destroy(&etypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&gtypes);

    return 0;
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
  
    print_attributes (&graph);
    
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
    
    igraph_vector_ptr_destroy(&result);
    igraph_destroy(&graph);

    if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;
  
    return 0;
}
