# include "postp.h"

Node  ***node_matrix(int rows, int columns){
    Node ***m;
    int i;
        /* allocate pointers to rows */
    m=(Node ***) malloc(rows*sizeof(Node**));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(Node **) calloc( rows*columns, sizeof(Node*));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}
