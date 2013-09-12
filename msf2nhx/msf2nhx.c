# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>
# include <math.h>

/*************************************************************/
/*************************************************************/
/*************************************************************/
/* type declarations: */        
#define BUFFLEN 150


# define UPGMA 1
# define NEIGHBOR_JOINING 2
# define ALMT_NAME_LENGTH 30

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int * seq_gaps;
    int * column_gaps;
    double **seq_dist;
    char * query_seq;
    char * query_seq_bkp;
    char query_name[ALMT_NAME_LENGTH];
    int query_gaps;
}  Alignment;
# define NAME_LENGTH 30
# define ASCII 128 /* number of characters */

# define ROOT  1
# define INNER 2
# define LEAF  4
# define CENTRAL 8 /*for handling unrooted trees*/  

# define UP    1
# define DOWN  2
# define LEFT  4
# define RIGHT 8

typedef struct Node {

    struct Node *left, *right, *parent;
    int   id;
    int   type;
    int  number_of_leaves; /* in the corresponding subtree */
    int marked;
    int  bin; /*similarity bin the node belongs to*/
    double dist_to_parent, dist_to_left, dist_to_right;
    double avg_sim;
    char * seq;
    char * name;
    char consensus;
    int population[ASCII]; /* population of amino acid types */
    double entropy;
    double entropy_complement;

} Node;

typedef struct Tree{

    Node *root;
    Node *leaf; /* this is actually the beginning of the node storage array */
    Node  ***group_root,  ***group_root_sim; 
    int size; /* leaves + inner nodes */
    int no_of_leaves;

} Tree;

Node ***node_matrix(int rows, int columns);
/*************************************************************/
/*************************************************************/
/*************************************************************/
/* function declarations: */        
int read_clustalw ( char * filename, Alignment * alignment);
int build_tree (Alignment * alignment, Tree * tree, int method);
int print_tree (FILE * fptr,Node * node );
int process_almt (Alignment * alignment);

/*************************************************************/
/*************************************************************/
/*************************************************************/
int main ( int argc, char * argv[]) {

    int method;
    Alignment alignment;
    Tree tree;
    int retval;
    if ( argc < 2 ) {
	fprintf (stderr, "Usage: %s <msf file> [nj]\n", argv[0]);
	exit (0);
    }
    method = UPGMA;
    if ( argc > 2  && !strncmp(argv[2], "nj", 2) ) {
	method = NEIGHBOR_JOINING;
    }
    
    /* read in the alignment */
    /* find  number of gaps, seq distances;  find query in the alignment*/
    retval = read_clustalw (argv[1], &alignment);
    if (retval) exit(retval);
    /* find seq distances*/
    retval =  process_almt ( &alignment);
    if (retval) exit(retval);
    /* build the seq similarity tree */
    memset ( &tree, 0, sizeof(Tree) );  
    retval  = build_tree (&alignment, &tree, method);
    if (retval) return retval;
   
    print_tree (stdout, tree.root);

    return 0;
 
}



/*************************************************************/
/*************************************************************/
/*************************************************************/
/*************************************************************/
/*************************************************************/
/*************************************************************/
/***********************************************************/
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






/*************************************************************/
/***********************************************************/

char **chmatrix(int rows, int columns){
    char **m;
    int i;
        /* allocate pointers to rows */
    m=(char **) malloc(rows*sizeof(char*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(char *) calloc( rows*columns, sizeof(char));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

double **dmatrix(int rows, int columns){
    double **m;
    int i;
        /* allocate pointers to rows */
    m=(double **) malloc(rows*sizeof(double*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(double *) calloc( rows*columns, sizeof(double));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

void * emalloc(int	size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}


/****************************************************************/

int read_clustalw ( char * filename, Alignment * alignment){
    
    FILE * fptr = NULL;
    char line[BUFFLEN];
    int  number_of_seqs, almt_length, ctr;
    int * seq_pos, pos;
    int * pos_ctr;
    char * seq_ptr;
    char ** sequence;
    char ** name;
    char curr_name[BUFFLEN];
     
    /* open file */
    fptr = fopen ( filename, "r");
    if ( !fptr ) {
	fprintf ( stderr, "Could not open %s.\n", filename);
	return 1;
    }

    memset (alignment, 0, sizeof(Alignment) );
    
     
    /* find the alignment length info */
    almt_length = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr(line, "MSF:" ) ){
	    sscanf (line, "%*s %d", &almt_length);
	    break;
	}
    }
    if ( almt_length ) {
	/* printf ( "Alignment length in %s is %d.\n", cwname, almt_length); */
    } else {
	fprintf ( stderr,
		  "Alignment length info not found in %s. Is the format gcg?\n",
		  filename );
	return 1;
    }

    /* determine the number of sequences */
    number_of_seqs = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( ! strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    number_of_seqs++;
	}
    }
    if ( number_of_seqs ) {
    } else {
	fprintf ( stderr, "No sequences found in %s. Is the format gcg?\n",
		  filename);
	return 1;
    }
    
    
    /* allocate */
    sequence = chmatrix (number_of_seqs, almt_length);
    if ( !sequence ) return 1;
    
    alignment->seq_dist = dmatrix ( number_of_seqs, number_of_seqs);
    if ( !alignment->seq_dist ) return 1;
    
    name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
    if ( !name ) return 1;

    seq_pos = (int *) emalloc ( number_of_seqs*sizeof(int));
    if ( !seq_pos ) return 1;
    
    /* read in */
    rewind(fptr);
    ctr = 0;

    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if (!  strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    sscanf (line, "%*s %s", curr_name);
	    if ( strcmp (curr_name, alignment->query_name) ) {
		sprintf (name[ctr], "%s", curr_name);
		ctr ++;
	    }
	}
    }
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( isspace (line[0] ) ) continue;
	sscanf (line, "%s", curr_name);
	
	ctr = 0;
	while (  ctr <number_of_seqs &&  strcmp (name[ctr], curr_name) ) ctr++;
	if ( ctr >= number_of_seqs ) {
	    fprintf ( stderr, "The name %s not found in the header of %s.\n",
		      curr_name,  filename);
	    return 1;
	}
	seq_ptr = sequence [ctr];
	pos_ctr = seq_pos + ctr;
	
	pos = 0;
	while ( ! isspace(line[pos]) ) pos++;
	while  (line[pos] != '\n' && pos < BUFFLEN) {
	    if ( !  isspace(line[pos] ) ){
                /* --> turn to uppercase */
		if ((line[pos]>=97)&&(line[pos]<=122)) {line[pos] -= 32;}
                /* turn tweedle to dot */
		if ( line[pos]==126)                   {line[pos]  = 46;} 
		seq_ptr [ *pos_ctr ] = line[pos];
		(*pos_ctr)++;
	    }
	    pos ++;
	}
    }
    fclose(fptr);

    /* sanity check */
    for (ctr=0; ctr < number_of_seqs; ctr++ ) {
	if ( seq_pos[ctr] != almt_length ) {
	    fprintf (stderr,
		     "Sequence %s is shorter (%d position) than the alignment.\n",
		     name[ctr],  seq_pos[ctr]);
	    return 1;
	}
    }

    /* return values */
    alignment->number_of_seqs = number_of_seqs;
    alignment->length         = almt_length;
    alignment->sequence       = sequence;
    alignment->name           = name;

    /* free */
    free (seq_pos);

    return 0;
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
int upgma_and_nj (Alignment * alignment, Tree * tree, int method);
/********************************************************************************/

int build_tree (Alignment * alignment, Tree * tree, int method) {
    

    int retval;
    int upgma_and_nj (Alignment * alignment, Tree * tree, int  method); 
     
    switch (method) {
    case UPGMA:
    case NEIGHBOR_JOINING:
	retval =  upgma_and_nj (alignment, tree, method); /* should fill the tree->size value*/
	if ( retval ) return retval;
	break;
    default:
	fprintf ( stderr,"Unrecognized tree building method.\n");
	return 1;
    }
    

    return 0;
}


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

/**********************************************************************************************/
int place_root ( Node *leaf, int no_of_nodes, Node * root) {
    
    int edge_ctr = 0, longest_ctr = 0, i;
    int find_path ( Node ** path, int * edge_ctr, Node ** longest_path,
		    int * longest_ctr,Node * node, Node* from);
    Node * node;

    Node** path, ** longest_path;
    if ( ! (path = emalloc (no_of_nodes*sizeof(Node*) ) ) ) exit (1);
    if ( ! (longest_path = emalloc (no_of_nodes*sizeof(Node*) ) ) ) exit (1);
    /* print all possible paths: */
    /* find any leaf to start fom: */
    longest_ctr = 0;
    edge_ctr = 0;
    
    for (i = 0; i < no_of_nodes; i++ ) {
	node = leaf + i;
    }
   
    for (i = 0; i < no_of_nodes; i++ ) {
	if ( (leaf + i) -> type != LEAF  ) continue;
	find_path (path, &edge_ctr, longest_path,
		   &longest_ctr, leaf + i, NULL);   }
    
    /* place the root above the middle node on that path */
    {
	int insert_node (Node * new_node,
			 Node * old_node_1, Node * old_node_2);
	int reorient_nodes ( Node * node, Node * from );
	root->type = ROOT;
	insert_node (root, longest_path[longest_ctr/2-1],
		     longest_path[longest_ctr/2]);
 	reorient_nodes ( root, NULL);
    }
    return 0;
}

/**********************************************************************************************/
int insert_node ( Node * new_node, Node * old_node_1, Node * old_node_2) {

    new_node->left = old_node_1;
    new_node->right= old_node_2;

    if  ( old_node_1 ==  old_node_2->left ) {
	old_node_2->left = new_node;
    } else if  ( old_node_1 ==  old_node_2->right ){
	old_node_2->right = new_node;
    } else if  ( old_node_1 ==  old_node_2->parent ){
	old_node_2->parent = new_node;
	old_node_2->dist_to_parent = old_node_1->dist_to_parent =  old_node_2->dist_to_parent/2;
	
    } else {
  	fprintf ( stderr, "Error inserting node.\n");
	exit (1);
    }

    if  ( old_node_2 ==  old_node_1->left ) {
	old_node_1->left = new_node;
    } else if  ( old_node_2 ==  old_node_1->right ){
	old_node_1->right = new_node;
    } else if  ( old_node_2 ==  old_node_1->parent ){
	old_node_1->parent = new_node;
 	old_node_2->dist_to_parent = old_node_1->dist_to_parent =  old_node_1->dist_to_parent/2;
   } else {
  	fprintf ( stderr, "Error inserting node.\n");
	exit (1);
    }

    new_node->dist_to_left  =  old_node_1 -> dist_to_parent;
    new_node->dist_to_right =  old_node_2 -> dist_to_parent;
    return 0;
}
/**********************************************************************************************/
int reorient_nodes ( Node * node, Node * from ) {

    Node *left, *right;
    double dist_to_left, dist_to_right, dist_to_parent;
    /* decide on "left" and "right" */
    if (  from ==  node->parent  ) {
	
	left = node->left;
	right = node->right;
	
	dist_to_left   = node->dist_to_left;
	dist_to_right  = node->dist_to_right;
	dist_to_parent =  node->dist_to_parent;
	
    } else if (from  ==  node->left  ) {
	
	right = node->right;
	left  = node->parent;
	
	dist_to_right  = node->dist_to_right;
	dist_to_left   = node->dist_to_parent;
	dist_to_parent = node->dist_to_left;
	
    } else if ( from == node->right  ) {
	
	left   = node->left;
	right  = node->parent;
	
	dist_to_left   = node->dist_to_left;
	dist_to_right  = node->dist_to_parent;
	dist_to_parent = node->dist_to_right;

	
	
    } else {
	fprintf ( stderr, "Error traversing tree in reorient_nodes.\n"); 
	exit (1);
    }
    node->parent = from; 
    node->left   = left; 
    node->right  = right; 
    node->dist_to_left   = dist_to_left;
    node->dist_to_right  = dist_to_right;
    node->dist_to_parent = dist_to_parent;
 
    if ( node->type != LEAF ) { 
	reorient_nodes (left,  node); 
	reorient_nodes (right, node);
    }
    return 0;
}
/**********************************************************************************************/
int find_path ( Node ** path, int * edge_ctr, Node ** longest_path,
		int * longest_ctr,Node * node, Node* from) {

    if ( ! node )  {
	fprintf ( stderr, "Error traversing tree in find_path (incoming node = 0x0).\n");
	exit (1);
    }
    
    /* add yourself to the path */
    path[*edge_ctr] = node;
    /* increase the counter     */
    (*edge_ctr) ++;
    if (  node->type != LEAF ) {
	Node *left, *right;
	/* decide on "left" and "right" */
	if (  node->parent == from ) {
	    left = node->left;
	    right = node->right;
	} else if ( node->left == from ) {
	    right = node->right;
	    left  = node->parent;
	} else if ( node->right == from ) {
	    right = node->parent;
	    left  = node->left;
	} else {
	    fprintf ( stderr, "Error traversing tree in find_path.\n");
	    exit (1);
	}
	/* go down the left */
	find_path (path, edge_ctr, longest_path, longest_ctr,  left, node);
	/* go down the   right */
	find_path (path, edge_ctr, longest_path, longest_ctr,  right, node);
    } else {
	if ( *edge_ctr == 1 ) { /* we are at the beginning */
	    find_path (path, edge_ctr, longest_path, longest_ctr,node->parent, node);
	} else {  /* we are at the end */
	    if ( *edge_ctr > *longest_ctr ) {
		*longest_ctr = *edge_ctr;
		memcpy ( longest_path, path, (*longest_ctr)*sizeof(Node*));
	    }
	}
    }
    /* get yourself off the path */
    path[*edge_ctr] = NULL;
    /* decrease the counter */
    (*edge_ctr) --;
  
    return 0;
}

double * score_array;

int pos_cmp (const void * a0, const void * b0) {
    
    int * a= (int*) a0;
    int * b= (int*)b0;
    if ( score_array[*a] > score_array[*b]) {
	return 1;
    }
    if ( score_array[*a] < score_array[*b]) {
	return -1;
    }
    return 0;
}

/********************************************************************************************/

int array_qsort (int * sorted_pos, double * sa, int sequence_length ) {
    /* position comparison function */
    score_array = sa;

    qsort (sorted_pos, sequence_length, sizeof(int), pos_cmp);

    return 0;
}


/**********************************************************************************************/
int rank_inner_nodes (Node *root, int  no_of_nodes) {

    int i, inner_node_ctr;
    int  * node_label_sorted;
    double *distance;
    Node ** label2node;
    int find_leaves_dist ( Node *node, int *inner_node_ctr, int * node_label_sorted, Node ** label2node, double * distance );
    
    /* allocate */
    if ( ! ( node_label_sorted = emalloc (no_of_nodes*sizeof(int)))) return 1;
    if ( ! ( distance  = emalloc (no_of_nodes*sizeof(double)))) return 1;
    if ( ! ( label2node  = emalloc (no_of_nodes*sizeof(Node*)))) return 1;
    
    /*find distance to the root from each inner node */
    inner_node_ctr = 0;
    find_leaves_dist ( root, &inner_node_ctr, node_label_sorted, label2node, distance );
   

    /* sort inner nodes by that distance */
    array_qsort ( node_label_sorted, distance, inner_node_ctr);
    /* assign them their sorted number */
    for (i = 0; i < inner_node_ctr; i++ ) {
	label2node[node_label_sorted[i]]->id = i+1;
    }
   
    free (node_label_sorted);
    free (label2node);
    free (distance);

    return 0;
}
/**********************************************************************************************/
int find_leaves_dist ( Node *node, int *inner_node_ctr, int * node_label_sorted,
		       Node ** label2node, double * distance ){
    double distance_to_parent (Node * node);
 
    if ( node -> type == LEAF ) {
    } else {
	node_label_sorted[*inner_node_ctr] = *inner_node_ctr;
	label2node[*inner_node_ctr] = node;
	distance[*inner_node_ctr] = distance_to_parent (node);
 	(*inner_node_ctr) ++;
	find_leaves_dist ( node->left, inner_node_ctr, node_label_sorted, label2node, distance );
	find_leaves_dist ( node->right, inner_node_ctr, node_label_sorted, label2node, distance );
    }

    return 0;
}



/**********************************************************************************************/
double distance_to_parent (Node * node){

    Node * current = node ;
    double distance= 0.0;
    
    while ( current->parent ) {
	distance += current -> dist_to_parent;
	current = current->parent;
    };

    return distance;
    
}

/********************************************************************************/
/********************************************************************************/
int  set_precedence_table ( Node* node, int number_of_nodes, int **is_precedent){

    int node_id;
    
    if (node->left->type != LEAF )  {
	is_precedent[node->id][node->left->id] = 1;
	set_precedence_table ( node->left, number_of_nodes, is_precedent);
	for ( node_id = node->left->id + 1; node_id <= number_of_nodes; node_id++ ) {
	    is_precedent[node->id][node_id] |= is_precedent[node->left->id][node_id];
	}
    }
    if (node->right->type != LEAF ) {
	is_precedent[node->id][node->right->id] = 1;
	set_precedence_table ( node->right, number_of_nodes, is_precedent);
	for ( node_id = node->right->id + 1; node_id <= number_of_nodes; node_id++ ) {
	    is_precedent[node->id][node_id] |= is_precedent[node->right->id][node_id];
	}
    }
     
    return 0;
}

/**************************************************************/


int upgma_and_nj (Alignment * alignment, Tree * tree, int method) {
    double **seq_dist, distance1, distance2;
    int retval;
    int node_ctr, no_seqs = alignment->number_of_seqs, tree_size;
    int closest1, closest2,* current_list;
    int upper, last_ctr;
    Node * node;
    int ( *closest_in_curr_list ) (double **seq_dist, Node * node, int tree_size, int new_node, 
				   int *  current_list, int * closest1, int *closest2,
				   double * dist_1_ptr, double * dist_2_ptr );
    int closest_in_curr_list_nj (double **seq_dist, Node * node, int tree_size, int new_node, 
			         int *  current_list, int * closest1, int *closest2,
				 double * dist_1_ptr, double * dist_2_ptr );
    /* I need dummy here to make upgma and nj defs look formally the same */ 
    int closest_in_curr_list_upgma (double **seq_dist, Node * node, int tree_size, int dummy, 
			     int *  current_list, int * closest1, int *closest2,
				    double * dist_1_ptr, double * dist_2_ptr );

    switch (method) {
    case UPGMA:
	closest_in_curr_list = closest_in_curr_list_upgma;
	break;
    case NEIGHBOR_JOINING:
	closest_in_curr_list = closest_in_curr_list_nj;
	break;
    default:
	fprintf ( stderr,"Unrecognized tree building method.\n");
	return 1;
    }

    /*****************************************/
    /*****************************************/
    /* find pairwise distances btw seqs      */
    /*****************************************/
    /*****************************************/
    seq_dist = alignment->seq_dist;
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
    upper = (method == UPGMA) ? tree_size : tree_size - 1;
    for ( node_ctr=no_seqs; node_ctr < upper; node_ctr++ ) {
	retval = closest_in_curr_list (seq_dist, node, tree_size, node_ctr,
				       current_list, &closest1, &closest2, &distance1, &distance2);
	if (retval) return retval;
	/* hack, so I can order the nodes in the tree: */
	if ( distance1 <= 0 ) distance1 = 0.001;
	if ( distance2 <= 0 ) distance2 = 0.001;
	
	/* fill in the new  node fields */ 
	node[node_ctr].left  = &node[closest1];
	node[node_ctr].dist_to_left   = distance1;
	node[node_ctr].right = &node[closest2];
	node[node_ctr].dist_to_right  = distance2;
	node[node_ctr].type  = INNER;
	node[node_ctr].id    = tree_size - node_ctr; /* this will server as the "rank" */
	node[node_ctr].number_of_leaves = node[closest1].number_of_leaves +
	    node[closest2].number_of_leaves;
	node[closest1].parent = & node[node_ctr];
	node[closest2].parent = & node[node_ctr];
	node[closest1].dist_to_parent = distance1;
	node[closest2].dist_to_parent = distance2;
	/* remove the two from the current list, and replace them with the parent */ 
	current_list[closest1] = 0;
	current_list[closest2] = 0;
	current_list[node_ctr] = 1;
    }

    last_ctr = node_ctr - 1;
    tree->leaf = node;
    (node + tree_size -1)->type = ROOT;
    tree->root = node + tree_size - 1;

    
    if(method == NEIGHBOR_JOINING) {
	
	int place_root ( Node *leaf, int no_of_nodes, Node * root) ;
	int rank_inner_nodes (Node *root, int  no_of_nodes);
	/* make the last two nodes  parents of each other */
	/* place the root */
	
	for ( node_ctr=0; node_ctr < upper; node_ctr++ ) {
	    if ( current_list[node_ctr] ) {
		node[node_ctr].parent =  & node[last_ctr];
		node[last_ctr].parent =  & node[node_ctr];
		break;
	    }
	}
	place_root ( tree->leaf,  tree_size-1, tree->root);
	rank_inner_nodes ( tree->root,  tree_size);
    }

    tree->size = tree_size;
    tree->no_of_leaves = no_seqs; 
    free (current_list);
   
    return 0;
}
/********************************************************************************/
/********************************************************************************/
int closest_in_curr_list_upgma (double **seq_dist, Node * node, int tree_size, int dummy, 
				int *  current_list, int * closest1,
				int *closest2, double * dist_1_ptr, double * dist_2_ptr ){
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
	 
	    if ( distance < min_distance) {
		min_distance = distance;
		*closest1 = ctr1;
		*closest2 = ctr2;
	    }
	}
    }
    *dist_1_ptr = min_distance/2;
    *dist_2_ptr = min_distance/2;
    return 0;
}
/********************************************************************************/
/********************************************************************************/
int closest_in_curr_list_nj (double **seq_dist, Node * node, int tree_size, int new_node, 
			     int *  current_list, int * closest1, int *closest2,
			     double * dist_1_ptr, double * dist_2_ptr ){
    /* actually, for nj, they are not the closest, but rather the pair which,
       if joined next, gives the minimum total sum of branch lengths */
    /* new_node  is where I'll put the parent for the closest pair I find */
    
    int ctr1, ctr2, curr_list_length;
    int no_leaves = (tree_size+1)/2;
    double sum, min_sum;
    double avg1, avg2, avg1_min, avg2_min;
    static double ** dist_table = NULL;
    double  nj_sum_of_branch_lengths (double ** dist_table,  int * current_list,
				   int  node_ctr1, int node_ctr2,
				   int tree_size, double *avg1_ptr, double *avg2_ptr);

    if ( ! dist_table ) { /* initialize distance table */
	/* allocate */
	dist_table = dmatrix (tree_size, tree_size);
	/* copy distances btw the leaves */
	for (ctr1=0; ctr1 < no_leaves; ctr1++ ) {
	    for (ctr2=ctr1+1; ctr2 < no_leaves; ctr2++ ) {
		dist_table[ctr1][ctr2] = dist_table[ctr2][ctr1] = seq_dist[ctr1][ctr2];
	    }
	}
    }
    
    curr_list_length = 0;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	curr_list_length += current_list[ctr1];
    }

    avg1_min = 0;
    avg2_min = 0;
    min_sum = 10000;
    
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    sum = nj_sum_of_branch_lengths ( dist_table, current_list,
					     ctr1, ctr2, tree_size, &avg1, &avg2);
	    if ( sum < min_sum) {
		min_sum   = sum;
		*closest1 = ctr1;
		*closest2 = ctr2;
		avg1_min  = avg1; 
		avg2_min  = avg2; 
	    }
	}
    }
    
    
    *dist_1_ptr = ( dist_table[*closest1][*closest2]+ avg1_min - avg2_min)/2;
    *dist_2_ptr = ( dist_table[*closest1][*closest2]+ avg2_min - avg1_min)/2;

    dist_table[new_node][*closest1] =  dist_table[*closest1][new_node] = *dist_1_ptr;
    dist_table[new_node][*closest2] =  dist_table[*closest2][new_node] = *dist_2_ptr;

    /* update the distance table */
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	if ( ctr1 == new_node )   continue;
	if ( ctr1 == *closest1 ||  ctr1 == *closest2) continue;
	dist_table[new_node][ctr1] = dist_table[ctr1][new_node] =
	    ( dist_table[*closest1][ctr1] + dist_table[*closest2][ctr1]
	      -  dist_table[*closest2][*closest1] )* 0.5;
    }
	
         
    return 0;
}


/***************************************************************************************/
double  nj_sum_of_branch_lengths (double ** dist_table,  int * current_list,
				  int  node_ctr1, int node_ctr2, int tree_size,
				  double *avg1_ptr, double *avg2_ptr) {
    double sum, term;
    double avg[2];
    int current_pair[2] = {node_ctr1,node_ctr2};
    int ctr1, ctr2, n;
    
    sum =  dist_table[node_ctr1][node_ctr2]/2;


    for (ctr1 = 0; ctr1 < 2; ctr1++ ) {
	n =0;
	avg[ctr1] = 0;
	for (ctr2= 0; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    if ( ctr2 == node_ctr1 || ctr2  == node_ctr2 )
		continue;
	    avg[ctr1] += dist_table[ current_pair[ctr1]][ ctr2];
	    n ++;
	}
	if ( !n ) {
	    printf ( "** %d %d \n", node_ctr1, node_ctr2);
	    printf ( "available:  \n");
	    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
		if ( ! current_list[ctr1]) continue;
		printf ( "%4d", ctr1);
	    }
	    printf ( "\n");
	    exit (1);
	}
    
 	    
	avg[ctr1] /= n;
    }

    sum += (avg[0] + avg[1])/2;
    
    *avg1_ptr = avg[0];
    *avg2_ptr = avg[1];
    
    term = 0;
    n = 0;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	if ( ctr1 == node_ctr1 ||  ctr2  == node_ctr2 )
	    continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    if ( ctr2 == node_ctr1 ||  ctr2 == node_ctr2 )
		continue;
	    term += dist_table[ctr1][ctr2];
	    n ++;
	}
    }
    sum += term/n;
		
	    
    return sum;
}
/*********************************************************************************************/
/*********************************************************************************************/

int seq_pw_dist (Alignment * alignment) {
    char *amino_acid_order = "ABCDEFGHIKLMNPQRSTVWXYZ";

    int blosum62[]={
	4,
	-2,  4,
	0, -3,  9,
	-2,  4, -3,  6,
	-1,  1, -4,  2,  5,
	-2, -3, -2, -3, -3,  6,
	0, -1, -3, -1, -2, -3,  6,
	-2,  0, -3, -1,  0, -1, -2,  8,
	-1, -3, -1, -3, -3,  0, -4, -3,  4,
	-1,  0, -3, -1,  1, -3, -2, -1, -3,  5,
	-1, -4, -1, -4, -3,  0, -4, -3,  2, -2,  4,
	-1, -3, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5,
	-2,  3, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6,
	-1, -2, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7,
	-1,  0, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,
	-1, -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5,
	1,  0, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,
	0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,
	0, -3, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4,
	-3, -4, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,
	0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1, -2, -1,
	-2, -3, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2, -1,  7,
	-1,  1, -3,  1,  4, -3, -2,  0, -3,  1, -3, -1,  0, -1,  3,  0,  0, -1, -2, -3, -1, -2,  4};
    
     int i,j,pos, ctr, off_diag;
    int aao_strlen;
    int number_of_similar_sites;
    /* int number_of_common_sites; */
    int similarity[128][128];
    char * seq_i, *seq_j;
    int char_i, char_j;
    double avg;


    /*get the similarity matrix to a usable form */
    aao_strlen = strlen(amino_acid_order);
    
    ctr = 0;
    avg = 0;
    off_diag = 0;
    for(i=0;i<aao_strlen;i++){
	char_i = (int) amino_acid_order [i];
	for (j=0;j<=i;j++){
	    char_j = (int) amino_acid_order [j];
	    similarity[char_i][char_j] = similarity[char_j][char_i] = blosum62[ctr];
	    if ( i != j  && blosum62[ctr] > 0) {
		avg +=  blosum62[ctr];
		off_diag ++;
	    }
	    ctr++;
	}
    }
    
    avg /= off_diag;
    avg = ceil(avg);
    
    int length_i, length_j, shorter;
    
    for (i=0; i<alignment->number_of_seqs; i++) {
	seq_i = alignment->sequence[i];
	length_i = alignment->length - alignment->seq_gaps[i];
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    length_j = alignment->length - alignment->seq_gaps[j];
	    seq_j = alignment->sequence[j];
	    /* number_of_common_sites = 0; */
	    number_of_similar_sites = 0;
	    for (pos=0; pos<alignment->length; pos++) {
		if ( seq_i[pos] == '.' ) continue;
		if ( seq_j[pos] == '.' ) continue;
		/* number_of_common_sites ++; */
		if ( similarity [(int)seq_i[pos]] [(int) seq_j[pos] ] >= avg ) {
		    number_of_similar_sites ++;
		}
	    }
	    shorter = (length_i<length_j) ? length_i : length_j;
	    alignment->seq_dist[i][j] = 1 - (double) number_of_similar_sites/shorter;
	    /* alignment->seq_dist[i][j] = 1 - (double) number_of_similar_sites/number_of_common_sites; */
	    alignment->seq_dist[j][i] = alignment->seq_dist[i][j];
	}
    }
# if 0
    for (i=0; i<alignment->number_of_seqs; i++) {
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    printf (" %4d  %4d  %10s %10s  %8.3lf\n",
		    i+1, j+1, alignment->name[i],  alignment->name[j], alignment->seq_dist[i][j]);
	}
    }
	exit (0);
# endif   
    return 0;
}
/************************************************************/

int count_gaps (Alignment * alignment) {

    int s, c;
    alignment->seq_gaps    = (int *) emalloc (alignment->number_of_seqs*sizeof(int));
    if (!alignment->seq_gaps) return 1;
    alignment->column_gaps = (int *) emalloc (alignment->length*sizeof(int));
    if (!alignment->column_gaps) return 1;
    for ( s=0; s<alignment->number_of_seqs; s++ ) {
	for ( c=0; c<alignment->length; c++) {
	    if ( alignment->sequence[s][c] == '.' ) {
		alignment->column_gaps[c] ++;
		alignment->seq_gaps[s] ++;
	    }
	}
    }
    return 0;
}
/************************************************************/
/************************************************************/

int process_almt ( Alignment *alignment) {
    
    int retval;
    
    /* gaps */
     int count_gaps (Alignment * alignment);

    /* gaps */
    count_gaps (alignment);

    /* seq dist */
    alignment->seq_dist = NULL;

    alignment->seq_dist =
	dmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    
    retval   = seq_pw_dist (alignment);
    return retval;
}
/************************************************************/
