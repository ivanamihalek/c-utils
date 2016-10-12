# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define PANIC( errmsg )                      \
   fprintf (stderr,"%s\n" errmsg);            \
   exit(1)

void *  emalloc ( size_t size) {
    void * retval;
    if (  (retval=  malloc(size))  ) {
	memset ( retval, 0, size);
	return retval;
    } else {
	PANIC ("Memory allocation failure.\n");
	return NULL;
    }
}

# define NODE_PRINT( fptr, node )                      \
    fprintf ( fptr, "id:%3d     type:%3d      value:%10x    parent:%10x    left:%10x   right:%10x \n",\
	    node ->id,  node ->type, node, node ->parent,  node ->left,  node ->right);\



# define ROOT  1
# define INNER 2
# define LEAF  4

typedef struct Node {

    struct Node *left, *right, *parent;
    int id;
    int type;

} Node;


int  tree_init (Node *root, Node *left, Node *right);
int  tree_sanity_check ( Node * node);
int  tree_print ( Node * node);
int  tree_print_nhx ( FILE * fptr, Node * node);
Node*  node_insert (Node *old, Node *new);

int main ( int argc, char * argv[]) {

    Node root = {0};
    Node left = {0}, right = {0};
    Node new  = {0};
    Node* curr_root_ptr = & root;
    Node* retnode;
    int retval;


    root.id = 1;
    left.id = 2;
    right.id = 3;
    
    tree_init (&root, &left, &right);

# if 1
    new. id = 4;
    retnode = node_insert (&root, &new);
    if ( retnode ) {
	curr_root_ptr = retnode;
    } else {
	curr_root_ptr = &root;
    }
# endif
    
    if ( retval = tree_sanity_check ( curr_root_ptr ) ) {
	printf ("Tree trouble: %d.\n", retval);
	return 1;
    } else {
	printf ("Tree passes sanity check.\n");
	tree_print (curr_root_ptr );
	tree_print_nhx (stdout, curr_root_ptr );
    }
    return 0;
}


int  tree_init (Node *root, Node *left, Node *right) {

    printf ("Initalizing tree ...\n");
    
    root->type = ROOT;
    root->left      = left;
    root->right     = right;
    root->parent    = NULL;

    left->type = LEAF;
    left->left      = NULL;
    left->right     = NULL;
    left->parent    = root;
    
    right->type = LEAF;
    right->left      = NULL;
    right->right     = NULL;
    right->parent    = root;
    
    printf ("                 ... OK\n");

    return 0;
    
}


int tree_print ( Node * node) {
    /* preorder print */
    
    if ( !node) return;
    
    NODE_PRINT(stdout,  node );      

    tree_print ( node->left);
    tree_print ( node->right);
    
    return 0;
}

int  tree_print_nhx (FILE* fptr,  Node * node){
    /*postorder*/
    if ( node->type == LEAF ) {
	fprintf ( fptr,  "%d", node->id );
    } else {
	fprintf ( fptr, "(" );
	if ( ! node->left ) {
	    NODE_PRINT (stderr, node);
	    PANIC ("Error in tree_print_nhx.\n");
	}
	tree_print_nhx ( fptr, node->left );
	fprintf ( fptr, "," );
	if ( ! node->right ) {
	    NODE_PRINT (stderr, node );
	    PANIC ("Error in tree_print_nhx.\n");
	}
	tree_print_nhx ( fptr, node->right );
	fprintf (fptr,  ")" );
    }
    if ( node->type == ROOT ) {
	fprintf (fptr, "\n");
    }
    return 0;
    
}


int tree_sanity_check ( Node * node) {
    /* preorder  */

    int retval, ok;
    
    if ( !node) return 0;

    switch (node->type) {
    case ROOT:
	ok = !(node->parent) &&  (node->left) && (node->right);
	break;
    case INNER:
	ok = (node->parent) &&  (node->left) && (node->right);
	break;
    case LEAF:
	ok = (node->parent) &&  !(node->left) && !(node->right);
    }
    if (!ok ) {
	printf ("Sanity check failure:\n");
	NODE_PRINT( stderr, node );      
	return 1;
	    
    }
    

    retval = tree_sanity_check ( node->left);
    if ( retval ) {
	return retval;
    }
    
    retval = tree_sanity_check ( node->right);
    if ( retval ) {
	return retval;
    }

    return 0;
}


Node *  node_insert (Node *old, Node *new){

    /* create new inner node above old */
    Node * new_inner   ;
    new_inner = (Node*) emalloc ( sizeof (Node) );
    printf ("\nMaybe I don't want node_insert to allocate - that might become expensive \n");
    printf ("when exploring new trees. Rather, I'll want to reuse the existing pool of nodes.\n");
    printf ("Also, I will need node deletion.\n");
    printf ("Then,  I can have some generic UPGMA and NJ, for the heck of it.\n");
    printf ("\n");
    
    
    /* attach old as its left child */
    new_inner->left = old;

    /* attach new as its right child */
    new_inner->right = new;

    /* adjust the types and pointers */
    new->parent = new_inner;
    new->type = LEAF;
    
    if ( old->type == ROOT ) {
	new_inner->type = ROOT;
	old->type = INNER;

	new_inner->parent = NULL;
	old->parent = new_inner;
	return new_inner;
	
    } else {
	new_inner->type = INNER;
	new_inner->parent = old->parent;
	/* parent of the old node now must point to the new_inner node as its  child*/
	if ( old->parent->left == old ) {
	    old->parent->left  = new_inner;
	} else {
	    old->parent->right = new_inner;
	}
	old->parent = new_inner;
	return NULL;
    }


}


