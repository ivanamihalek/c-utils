# include "postp.h"


int build_tree (Alignment * alignment, Node ** root_ptr, Node** leaf_ptr) {
    double **seq_dist = alignment->seq_dist, distance;
    int retval;
    int node_ctr, no_seqs = alignment->number_of_seqs, tree_size;
    int closest1, closest2,* current_list;
    Node * node;
    int closest_in_curr_list (double **seq_dist, Node * node, int tree_size, int *  current_list,
			      int *closest1, int *closest2, double * dist_ptr);

    /*****************************************/
    /*****************************************/
    /*****************************************/
    /*          build tree                   */
    /*****************************************/
    /*****************************************/
    /* allocate space */
    tree_size = 2*no_seqs -1;
    node = (Node *) emalloc ( tree_size*sizeof(Node) );
    /* initialize leaves to point to sequences from the alignment*/
    for ( node_ctr=0; node_ctr < no_seqs; node_ctr++ ) {
	node[node_ctr].id   = node_ctr;
	node[node_ctr].type = LEAF;
	node[node_ctr].seq  = alignment->sequence[node_ctr];
        node[node_ctr].name = alignment->name[node_ctr];
	node[node_ctr].number_of_leaves = 1;
    }
    /* initialize current list of nodes whose distance needs to be compared */
    current_list = (int*) emalloc ( tree_size*sizeof(int));
    for ( node_ctr=0; node_ctr < no_seqs; node_ctr++ ) {
	current_list[node_ctr] = 1; 
    }

    /* find children for each of the remaining nodes */
    for ( node_ctr=no_seqs; node_ctr < tree_size; node_ctr++ ) {
	retval = closest_in_curr_list (seq_dist, node, tree_size, current_list, &closest1, &closest2, &distance);
	if (retval) return retval;
	/* fill in the new  node fields */ 
	node[node_ctr].left  = &node[closest1];
	node[node_ctr].right = &node[closest2];
	node[node_ctr].type  = INNER;
	node[node_ctr].id    = tree_size - node_ctr; /* this will server as the "rank" */
	node[node_ctr].number_of_leaves = node[closest1].number_of_leaves +  node[closest2].number_of_leaves;
	node[closest1].parent = & node[node_ctr];
	node[closest2].parent = & node[node_ctr];
	node[closest1].dist_to_parent = distance/2;
	node[closest2].dist_to_parent = distance/2;
	/* remove the two from the current list, and replace them with the parent */ 
	current_list[closest1] = 0;
	current_list[closest2] = 0;
	current_list[node_ctr] = 1;
    }
   
    
    free (seq_dist);
    free (current_list);
    (node + tree_size -1)->type = ROOT;
    *root_ptr = node + tree_size -1;
    *leaf_ptr = node;
# if 0
    print_tree (stdout, *root_ptr);
    exit(0);
# endif
   
    return 0;
}
/********************************************************************************/
/********************************************************************************/
int closest_in_curr_list (double **seq_dist, Node * node, int tree_size, int *  current_list,
			  int * closest1, int *closest2, double * dist_ptr){
    int ctr1, ctr2;
    double distance, min_distance;
    double  node_distance ( double **seq_dist, Node* node1, Node* node2 );

    min_distance = 1000;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
        if ( ! current_list[ctr1]) continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    distance = node_distance ( seq_dist,  node+ctr1, node+ctr2 );
	    /* we need the average */ 
	    distance /= (node+ctr1)->number_of_leaves*(node+ctr2)->number_of_leaves;
	    if ( distance > 1 ) {
		printf ( "%d   %d  %d   %d  %8.2lf\n", ctr1, ctr2,  (node+ctr1)->number_of_leaves,
			 (node+ctr2)->number_of_leaves, distance);
		exit (0);
	    }
	    if ( distance < min_distance) {
		min_distance = distance;
		*closest1 = ctr1;
		*closest2 = ctr2;
	    }
	}
    }
    *dist_ptr = min_distance;
    return 0;
}
/********************************************************************************/
/********************************************************************************/
double  node_distance ( double **seq_dist, Node* node1, Node* node2 ){

    
    if (  (node1->type == LEAF) &&  (node2->type == LEAF) ) {
	return seq_dist[node1->id][node2->id];
    } else {
	double distance = 0.0;
	if ( node1->type != LEAF ) { /* left side dpeth first */
	    distance += node_distance( seq_dist, node1->left,  node2);
	    distance += node_distance( seq_dist, node1->right, node2);
	} else 	if ( node2->type != LEAF ) {
	    distance += node_distance( seq_dist, node1,  node2->left);
	    distance += node_distance( seq_dist, node1,  node2->right);
	}
	return distance;
    }
}
/********************************************************************************/
/********************************************************************************/
int print_tree (FILE * fptr,Node * node ) {
    
    if ( node->type != LEAF ) {
	fprintf ( fptr, "(");
	/* down the left branch*/ 
	print_tree ( fptr, node->left);
	fprintf ( fptr,","); 
	/* down the right  branch*/ 
	print_tree ( fptr, node->right);
	/* handle unrooted tree: */
	if (node->type != ROOT && node->parent->type==LEAF ) {
	    fprintf ( fptr,","); 
	    print_tree ( fptr, node->parent);
	}
	
	fprintf ( fptr, ")");
	fprintf ( fptr,"%d", node->id);
	if ( node->dist_to_parent > 0.0001 ) {
	    fprintf ( fptr,":%f ", node->dist_to_parent );
	} else {
	    fprintf ( fptr,":0.0001");
	}
	
    } else {
        fprintf ( fptr,"%s", node->name);
	if ( node->dist_to_parent > 0.0001 ) {
	    fprintf ( fptr,":%f ", node->dist_to_parent);
	} else {
	    fprintf ( fptr,":0.0001");
	}
    }

    return 0;
}
/********************************************************/
void print_debug_tree (FILE * fptr, Node *node) {

  
    if (node->type != LEAF ) {

	fprintf ( fptr, "(");
	/* down the left branch*/ 
	print_debug_tree ( fptr, node->left);
	fprintf ( fptr,","); 
	/* down the right  branch*/ 
	print_debug_tree ( fptr, node->right);
	/* handle unrooted tree: */
	if (node->type != ROOT && node->parent->type==LEAF ) {
	    fprintf ( fptr,","); 
	    print_debug_tree ( fptr, node->parent);
	}
	fprintf ( fptr, ")");
	
	fprintf ( fptr,"%d", node->id);
	if ( node->dist_to_parent > 0.0001 ) {
	    fprintf ( fptr,":%f \n", node->dist_to_parent );
	} else {
	    fprintf ( fptr,":0.0001\n");
	}
	
    } else {
	
        fprintf ( fptr,"%s [%d] ", node->name, node->id);
	if ( node->dist_to_parent > 0.0001 ) {
	    fprintf ( fptr,":%f \n", node->dist_to_parent);
	} else {
	    fprintf ( fptr,":0.0001\n");
	}
	
    }
    
    return ;
}
