
#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>

void check_attr(igraph_t *graph, int offset) {

    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "name")) {
	printf("No graph attribute `name`\n");
	exit(offset + 2);
    }
    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "type")) {
	printf("No graph attribute `type`\n");
	exit(offset + 3);
    }
    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "p")) {
	printf("No graph attribute `p`\n");
	exit(offset + 4);
    }
    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX,
				    "name")) {
	printf("No vertex attribute `id`\n");
	exit(offset + 5);
    }
    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE,
				    "weight")) {
	printf("No edge attribute `weight'\n");
	exit(offset + 6);
    }
}


/******************************************************/
int main(int argc, char*argv[]) {

    // the graph is just a list of from-to identifier pairs, one per line
    if (argc != 2) {
	fprintf (stderr, "Usage: %s  <path to graph file>\n", argv[0]);
	exit(1);
    }
    FILE * infile  =  efopen(argv[1], "r") ;
    igraph_t graph;
    int retval  =  igraph_read_graph_ncol (&graph, infile, NULL, 1,
					   IGRAPH_ADD_WEIGHTS_NO,  IGRAPH_UNDIRECTED);
    fclose(infile);
    
    igraph_error_handler_t* oldhandler;
    int result;

    igraph_i_set_attribute_table(&igraph_cattribute_table);
   
    check_attr(&graph, 10);


    igraph_destroy(&graph);

    return 0;
}
