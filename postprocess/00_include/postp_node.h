# define NODE_PRINT( fptr, node )                      \
    fprintf ( fptr, "id:%3d     type:%3d      value:%10p    parent:%10p    left:%10p  right:%10p \n",\
	    node ->id,  node ->type, node, node ->parent,  node ->left,  node ->right);\

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
    double dist_to_parent;
    char * seq;
    char * name;
    char consensus;
    int population[ASCII]; /* population of amino acid types */
    double entropy;

} Node;

# define NODE_PRINT( fptr, node )                      \
    fprintf ( fptr, "id:%3d     type:%3d      value:%10p    parent:%10p    left:%10p  right:%10p \n",\
	    node ->id,  node ->type, node, node ->parent,  node ->left,  node ->right);\

Node ***node_matrix(int rows, int columns);
