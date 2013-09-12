# include "postp.h"


int int_trace (Alignment * alignment, double *score){
    
    int retval, pos, tree_size, no_seqs = alignment->number_of_seqs;
    int node_ctr;
    Node *root, *leaf, *node_ptr;
    int  build_tree (Alignment * alignment, Node ** root, Node**leaf);
    /* build the seq similarity tree */
    retval  = build_tree(alignment, &root, &leaf);
    if (retval) return retval;
    
    /* calculate the score */
    tree_size = 2*no_seqs -1;
    for ( pos=0; pos < alignment->length; pos++) {
	score[pos] = 1;
	/* for a leaf, the consensus is aa at the position pos */
	for ( node_ctr=0; node_ctr <no_seqs ; node_ctr++ ) {
	    node_ptr = leaf+node_ctr;
	    node_ptr->consensus = node_ptr->seq[pos];
	}
	/* for an inner node, the consensus is aa if the same in both children */
	for ( node_ctr= no_seqs; node_ctr < tree_size; node_ctr++ ) {
	    node_ptr = leaf+node_ctr;
	    if ( node_ptr->left->consensus ==  node_ptr->right->consensus ) {
		node_ptr->consensus = node_ptr->left->consensus;
	    } else {
		score[pos] = node_ptr->id+1;
		break;
	    }
	}
    }
    return 0;
}
