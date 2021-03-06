# include "st_hbonds.h"
#include <errno.h>
#include <fcntl.h>

#include <sys/types.h>


void * emalloc(int	size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}



FILE * efopen(char * name, char * mode)
{

    FILE * fp;


    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,  
	      "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }

    return fp;

}




/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] 
 */
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

int **intmatrix(int rows, int columns){
    int **m;
    int i;
        /* allocate pointers to rows */
    m=(int **) malloc(rows*sizeof(int*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in intmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(int *) calloc( rows*columns, sizeof(int));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in intmatrix().\n");
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
	fprintf (stderr,"row allocation failure  in dmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(double *) calloc( rows*columns, sizeof(double));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in dmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}



/* free a  matrix  */
void free_matrix(void **m)
{
    free(m[0]);
    free(m);
}




/* sort array according to the score in the other */
/* I couldn't declare pos_cmp within array_qsort  bcs it   crashed on mac */

double * score_array;

int pos_cmp (const void * a0, const void * b0) {
    
    int * a= (int*) a0;
    int * b= (int*)b0;
    if ( score_array[*a] < score_array[*b]) {
	return 1;
    }
    if ( score_array[*a] > score_array[*b]) {
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

/**********************************************************/
/* get rid of spaces in a string */
int  string_clean ( char* string, int length) {
    int ctr;
    for (ctr = 0; ctr < length; ctr ++) {
	if ( isspace (string[ctr]) ) string[ctr] = '\0';
    }
    ctr=0;
    while ( !string[ctr]) ctr++;
    if ( ctr ) {
	memmove (string, string+ctr, length-ctr);
	memset ( string+length-1-ctr, 0, ctr);
    }

    return 0;
}

/**********************************************************/
char single_letter ( char code[]){

    switch ( code[0] ) {
    case 'A':
	switch ( code [1]) {
	case 'L':
	    return 'A';
	    break;
	case 'R':
	    return 'R';
	    break;
	case 'S':
	    switch ( code[2] ) {
	    case 'N':
		return 'N';
		break;
	    case 'P':
		return  'D';
		break;
	    }
	    break;
	}
	break;
    case 'C':
	return 'C'; 
	break;
    case 'G':
	/* the second letter is always L */ 
	switch ( code[2] ) {
	case 'U':
	    return 'E';
	    break;
	case 'N':
	    return  'Q';
	    break;
	case 'Y':
	    return 'G';
	    break;
	}
	break;
    case 'H':
	return  'H';
	break;
    case 'I':
	return  'I';
	break;
    case 'L':
	switch ( code [1]) {
	case 'E':
	    return 'L';
	    break;
	case 'Y':
	    return 'K';
	    break;
	}
	break;
    case 'M':
	return 'M';
	break;
    case 'P':
	switch ( code [1]) {
	case 'H':
	    return 'F';
	    break;
	case 'R':
	    return 'P';
	    break;
	}
	break;
    case 'S':
	return 'S';
	break;
    case 'T':
	switch ( code [1]) {
	case 'H':
	    return 'T';
	    break;
	case 'R':
	    return 'W';
	    break;
	case 'Y':
	    return 'Y';
	    break;
	}
	break;
    case 'V':
	return 'V';
	break;
	
    }


    fprintf (stdout, "Unrecognized amino acid code: %s.\n", code);
    return 0;
}

