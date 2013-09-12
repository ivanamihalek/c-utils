# include "postp.h"

int hybrid (Alignment * alignment, double *score){
    
    int retval, pos, tree_size, no_seqs = alignment->number_of_seqs;
    int rank, group_tag, g;
    double subsum, rho;
    Node *root, *leaf,  ***group_root;
    int build_tree (Alignment * alignment, Node ** root, Node**leaf);
    int entropy_recursive (Node *node, int pos) ;
    int find_group_roots (Node *node, Node*** group_root,  int rank, int * group_tag );
    /* build the seq similarity tree */
    retval  = build_tree(alignment, &root, &leaf);
    if (retval) return retval;

    /* find groups */
    group_root = node_matrix (no_seqs, no_seqs);
    for ( rank = 1; rank < no_seqs; rank++ ) {
	group_tag = -1;
	find_group_roots (root,  group_root,  rank, & group_tag );
    }
    /* calculate the score */
    tree_size = 2*no_seqs -1;
    for ( pos=0; pos < alignment->length; pos++) {
	entropy_recursive(root, pos);
	rho = 1.0;
	for ( rank = 1; rank < no_seqs; rank++ ) {
	    subsum = 0.0;
	    for ( g=0; g<rank; g++) {
		subsum  +=  group_root[rank][g]->entropy;
	    }
	    rho += subsum/rank;
	}
 	score[pos] =  rho;
    }

    free_matrix ((void**)group_root);
    return 0;
    
}

/******************************************************************************/

int entropy_recursive (Node *node, int pos) {

    if ( node->type == LEAF ) {
	memset (node->population, 0, ASCII*sizeof(int));
	node->population[(int) node->seq[pos] ] = 1;
	node->entropy = 0.0;
    } else {
	int ctr;
	int no_of_leaves = node->number_of_leaves;
	double fr;
	entropy_recursive (node->left, pos);
	entropy_recursive (node->right, pos);
	/* for each aa type: population = population left + pop right */
	node->entropy = 0.0;	
	for (ctr=0; ctr < ASCII; ctr++ ) {
	    node->population[ctr] = node->left->population[ctr] + node->right->population[ctr];
	    /* frequency is pop/number of seqs */
	    if ( node->population[ctr] ) {
		fr =  (double)node->population[ctr]/no_of_leaves;
		node->entropy -= fr*log(fr);
	    }
	}
   }

    return 0;
}
/******************************************************************************/

int find_group_roots (Node *node, Node*** group_root,  int rank, int * group_tag ) {
    if ( node->type==LEAF ) {
    } else {
	
	/* determining the roots for all the groups at the rank "rank" */
	/* remember:  at rank N, N is the smallest group root possible */
	if ( node->id == rank ) {
	    ++(*group_tag);
	    group_root[rank][*group_tag] = node;
	} else if( node->id < rank ){
	    
	    Node *left, *right;
	    left  = node->left;
	    right = node->right;

	    if ( left->id >= rank  || left->type==LEAF) {
		++(*group_tag);
		group_root[rank][*group_tag] = left;
	    } else {
		find_group_roots ( left, group_root, rank, group_tag);
	    }
	    if ( right->id >= rank  || right->type==LEAF) {
		++(*group_tag);
		group_root[rank][*group_tag] = right;
	    } else {
		find_group_roots ( right, group_root, rank, group_tag);
	    }
	}
	
    }
    return 0;
}
