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
#define NAME_LENGTH 50


# define UPGMA 1
# define NEIGHBOR_JOINING 2
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
int build_tree (int number_of_names, char ** name_list,double ** dist_table,
		Tree * tree, int method);
int print_tree (FILE * fptr, Node * node);
int read_in_dist ( char * filename, int * number_of_names_ptr,
		   char *** name_list_ptr, double *** dist_table_ptr);

/*************************************************************/
/*************************************************************/
/*************************************************************/
int main ( int argc, char * argv[]) {

    int method;
    int number_of_names;
    char **name_list;
    double **dist_table;
    Tree tree;
    int retval;

    printf ("this is unfinished and untested.\n");
    exit (1);
    
    if ( argc < 2 ) {
	fprintf (stderr, "Usage: %s <distance file> [nj]\n", argv[0]);
	exit (0);
    }
    method = UPGMA;
    if ( argc > 2  && !strncmp(argv[2], "nj", 2) ) {
	method = NEIGHBOR_JOINING;
    }
    
    /* read in the distances*/
    retval =  read_in_dist ( argv[1], &number_of_names, &name_list, &dist_table);
    if (retval) exit(retval);
    /* build the seq similarity tree */
    memset ( &tree, 0, sizeof(Tree) );  
    retval  = build_tree (number_of_names, name_list, dist_table, &tree, method);
    if (retval) return retval;
   
    print_tree (stdout, tree.root);

    return 0;
 
}



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

/*************************************************************/
int read_in_dist ( char * filename, int * number_of_names_ptr,
		   char *** name_list_ptr, double *** dist_table_ptr) {

    int number_of_names, declared_number_of_names;
    int aux, aux1, name_index;
    char   ** name_list;
    double **dist_table;
    double aux_double;
    FILE * fptr = NULL;
    char line[BUFFLEN];
    char name[NAME_LENGTH];

    /* open file */
    fptr = fopen ( filename, "r");
    if ( !fptr ) {
	fprintf ( stderr, "Could not open %s.\n", filename);
	return 1;
    }

    declared_number_of_names = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr(line, "end" ) ) break;
	if ( !declared_number_of_names) {
	    if ( sscanf (line , "%d ", &aux) ) {
		declared_number_of_names = aux;
	    }
	}
    }
    if ( !declared_number_of_names) {
	fprintf (stderr, "List length declaration not found.\n");
	fclose (fptr);
	return 1;
    }
    name_list = chmatrix (declared_number_of_names, NAME_LENGTH);
    if ( !name_list ) return 1;
    dist_table = dmatrix (declared_number_of_names, declared_number_of_names);
    if ( !dist_table ) return 1;
    
    number_of_names = 0;
    memset (name, 0, NAME_LENGTH);
    
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr(line, "end" ) ) break;
	if ( sscanf (line , "%d %s", &aux, name ) == 2) {
	    name_index = aux;
	    if ( name_index > declared_number_of_names  ) {
		fprintf (stderr, "index for name %s larger than declared: %d.",
			 name, declared_number_of_names);
		return 1;
	    }
	    if ( name_list[name_index-1][0] && ! strcmp (name, name_list[name_index-1]) ) {
		fprintf (stderr, "Two names found with the same index %d: %s  %s\n",
			 name_index, name_list[name_index-1], name);
		return 1;
	    }
	    sprintf ( name_list[name_index-1], "%s", name);
	    memset (name, 0, NAME_LENGTH);
	}
    }
    
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr(line, "end" ) ) break;
 	if ( sscanf (line , "%d %d %lf", &aux, &aux1, &aux_double) == 3) {
	    if (  aux > declared_number_of_names  ) {
		fprintf (stderr, "index %d larger than declared list length: %d.",
			 aux, declared_number_of_names);
		return 1;
	    }
	    if (  aux1 > declared_number_of_names  ) {
		fprintf (stderr, "index %d larger than declared list length: %d.",
			 aux1, declared_number_of_names);
		return 1;
	    }
	    aux--; aux1--;
	    dist_table[aux][aux1] =  dist_table[aux1][aux] = aux_double;	    
	}
   }
	
    fclose (fptr);
    * number_of_names_ptr = number_of_names; 
    * name_list_ptr  = name_list;
    * dist_table_ptr = dist_table;
    
    return 0;
}




/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/********************************************************************************/

int build_tree (int number_of_names, char ** name_list,double ** dist_table,
		Tree * tree, int method) {
    

    int retval;
    int upgma_and_nj (int number_of_names, char ** name_list,
		      double ** dist_table, Tree * tree, int  method); 
     
    switch (method) {
    case UPGMA:
    case NEIGHBOR_JOINING:
	/* this function should fill the tree->size value*/
	retval =  upgma_and_nj (number_of_names, name_list,dist_table, tree, method);
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


int upgma_and_nj (int number_of_leaves, char ** name_list,
		  double ** dist_table, Tree * tree, int method) {
    
    double  distance1, distance2;
    int retval;
    int node_ctr, tree_size;
    int closest1, closest2,* current_list;
    int upper, last_ctr;
    Node * node;
    int ( *closest_in_curr_list ) (double **dist_table, Node * node, int tree_size, int new_node, 
				   int *  current_list, int * closest1, int *closest2,
				   double * dist_1_ptr, double * dist_2_ptr );
    int closest_in_curr_list_nj (double **dist_table, Node * node, int tree_size, int new_node, 
			         int *  current_list, int * closest1, int *closest2,
				 double * dist_1_ptr, double * dist_2_ptr );
    /* I need dummy here to make upgma and nj defs look formally the same */ 
    int closest_in_curr_list_upgma (double **dist_table, Node * node, int tree_size, int dummy, 
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
    /*          build tree                   */
    /*****************************************/
    /*****************************************/
    /* allocate space */
    tree_size = 2*number_of_leaves -1;
    node = (Node *) emalloc ( tree_size*sizeof(Node) );
    /* initialize leaves to point to sequences from the alignment*/
    for ( node_ctr=0; node_ctr < number_of_leaves; node_ctr++ ) {
	node[node_ctr].id   = node_ctr;
	node[node_ctr].type = LEAF;
	node[node_ctr].name = name_list[node_ctr];
	node[node_ctr].number_of_leaves = 1;
    }
    /* initialize current list of nodes whose distance needs to be compared */
    current_list = (int*) emalloc ( tree_size*sizeof(int));
    for ( node_ctr=0; node_ctr < number_of_leaves; node_ctr++ ) {
	current_list[node_ctr] = 1; 
    }

    /* find children for each of the remaining nodes */
    upper = (method == UPGMA) ? tree_size : tree_size - 1;
    for ( node_ctr=number_of_leaves; node_ctr < upper; node_ctr++ ) {
	retval = closest_in_curr_list (dist_table, node, tree_size, node_ctr,
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
    tree->no_of_leaves = number_of_leaves; 
    free (current_list);
   
    return 0;
}
/********************************************************************************/
/********************************************************************************/
int closest_in_curr_list_upgma (double **dist_table, Node * node, int tree_size, int dummy, 
				int *  current_list, int * closest1,
				int *closest2, double * dist_1_ptr, double * dist_2_ptr ){
    int ctr1, ctr2;
    double distance, min_distance;
    double  node_distance ( double **dist_table, Node* node1, Node* node2 );

    min_distance = 1000;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
        if ( ! current_list[ctr1]) continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    distance = node_distance ( dist_table,  node+ctr1, node+ctr2 );
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
int closest_in_curr_list_nj (double **leaf_dist_table, Node * node,
			     int tree_size, int new_node, 
			     int *  current_list, int * closest1, int *closest2,
			     double * dist_1_ptr, double * dist_2_ptr ){
    /* actually, for nj, they are not the closest, but rather the pair which,
       if joined next, gives the minimum total sum of branch lengths */
    /* new_node  is where I'll put the parent for the closest pair I find */
    
    int ctr1, ctr2, curr_list_length;
    int no_leaves = (tree_size+1)/2;
    double sum, min_sum;
    double avg1, avg2, avg1_min, avg2_min;
    static double ** node_dist_table = NULL;
    double  nj_sum_of_branch_lengths (double ** dist_table,  int * current_list,
				   int  node_ctr1, int node_ctr2,
				   int tree_size, double *avg1_ptr, double *avg2_ptr);

    if ( ! node_dist_table ) { /* initialize distance table */
	/* allocate */
	node_dist_table = dmatrix (tree_size, tree_size);
	/* copy distances btw the leaves */
	for (ctr1=0; ctr1 < no_leaves; ctr1++ ) {
	    for (ctr2=ctr1+1; ctr2 < no_leaves; ctr2++ ) {
		node_dist_table[ctr1][ctr2] = node_dist_table[ctr2][ctr1]
		    = leaf_dist_table[ctr1][ctr2];
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
	    sum = nj_sum_of_branch_lengths ( node_dist_table, current_list,
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
    
    
    *dist_1_ptr = ( node_dist_table[*closest1][*closest2]+ avg1_min - avg2_min)/2;
    *dist_2_ptr = ( node_dist_table[*closest1][*closest2]+ avg2_min - avg1_min)/2;

    node_dist_table[new_node][*closest1] =  node_dist_table[*closest1][new_node] = *dist_1_ptr;
    node_dist_table[new_node][*closest2] =  node_dist_table[*closest2][new_node] = *dist_2_ptr;

    /* update the distance table */
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	if ( ctr1 == new_node )   continue;
	if ( ctr1 == *closest1 ||  ctr1 == *closest2) continue;
	node_dist_table[new_node][ctr1] = node_dist_table[ctr1][new_node] =
	    ( node_dist_table[*closest1][ctr1] + node_dist_table[*closest2][ctr1]
	      -  node_dist_table[*closest2][*closest1] )* 0.5;
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
