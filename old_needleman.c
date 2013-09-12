# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

#define FAR_FAR_AWAY  -10

int needleman_wunsch (int max_i, int max_j, double **similarity,
		      int *map_i2j, int * map_j2i);
char   **chmatrix(int rows, int columns);
double **dmatrix(int rows, int columns);
void *   emalloc(int	size);
void     free_cmatrix(char **m);
void     free_dmatrix(double **m);

int main ( int argc, char * argv[]) {

    // READ IN TWO DNA STRINGS (we expect shortish one here)
    // 
    double ** similarity;
    int * map_i2j, * map_j2i;
    int len1 ;
    int len2;
    int i, j;
    int total_aligned;
    int total_matching;
    if ( argc < 3 ) {
	fprintf (stderr, "Usage: %s <string 1>  <string 2>\n", argv[0]);
	exit (-1);
    }
    
    len1 = strlen ( argv[1] );
    len2 = strlen ( argv[2] );


    if ( !(similarity = dmatrix (len1, len2) )) return -1;
    if ( !(map_i2j    = emalloc (len1*sizeof(int)))) return -1;
    if ( !(map_j2i    = emalloc (len2*sizeof(int)))) return -1;
    

    for (i=0; i < len1; i++ ) {
	for (j=0; j < len2; j++ ) {
	    similarity[i][j] = (argv[1][i] == argv[2][j] )? 1 : -1;
	}
    }
    
    needleman_wunsch (len1, len2, similarity, map_i2j, map_j2i);

    total_aligned = 0;
    total_matching = 0;
    i = 0; j = 0;
    while (i< len1 ||  j < len2 ) {
	while (i< len1 &&  map_i2j[i] == FAR_FAR_AWAY) {
	    printf ( " %c   %c\n", argv[1][i], '.');
	    i ++;
	}
	while (j < len2 && map_j2i[j] == FAR_FAR_AWAY) {
	    printf ( " %c   %c\n", '.', argv[2][j] );
	    j ++;
	}
	if (i< len1 &&  j < len2 )  {
	    if ( map_i2j[i] != j || map_j2i[j] != i ) {
		fprintf (stderr, "alignment error (?)\n");
		exit (1);
	    }
	    printf ( " %c   %c\n", argv[1][i], argv[2][j] );
	    i ++;
	    j ++;
	    total_aligned ++;
	    if ( argv[1][i] ==  argv[2][j]) total_matching++;
	    
	}
	
    }
    printf ( "%d %d \n", len1, len2);
    printf ( "%d %d \n", total_aligned, total_matching);
    
    free_dmatrix (similarity);
    free ( map_i2j);
    free ( map_j2i);
	 
    
    return 0;
}

/////////////////////////////////////////////////

int needleman_wunsch (int max_i, int max_j, double **similarity,
		      int *map_i2j, int * map_j2i) {

    double **F; /*alignment_scoring table*/
    char ** direction;
    double gap_opening = -1;
    double gap_extension = -1;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    int i,j;

     /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;

    /* fill the table */
    for (i=0; i<= max_i; i++) {
	for (j=0; j<=max_j; j++) {

	    if ( !i && !j ) {
		F[0][0] = 0;
		direction[i][j] = 'd';
		continue;
	    }
	    
	    if ( i && j ){
		if ( direction[i-1][j] == 'i' ) {
		    /*  gap extension  */
		    i_sim = F[i-1][j] +  gap_extension;
		} else {
		    /*  gap opening  */
		    i_sim = F[i-1][j] +  gap_opening;
		}
		if ( direction[i][j-1] =='j' ) {
		    j_sim = F[i][j-1] +  gap_extension;
		} else {
		    j_sim = F[i][j-1] +  gap_opening;
		}
		diag_sim =  F[i-1][j-1] + similarity [i-1][j-1] ;
		
	    } else if ( j ) {
		if ( direction[i][j-1] =='j' ) {
		    j_sim = F[i][j-1] +  gap_extension;
		} else {
		    j_sim = F[i][j-1] +  gap_opening;
		}
	
		i_sim = diag_sim =  FAR_FAR_AWAY;
	    } else if ( i ) {
		if ( direction[i-1][j] == 'i' ) {
		    /*  gap extension  */
		    i_sim = F[i-1][j] +  gap_extension;
		} else {
		    /*  gap opening  */
		    i_sim = F[i-1][j] +  gap_opening;
		}
		
		j_sim = diag_sim =  FAR_FAR_AWAY;
	    } 

	    max_sim = diag_sim;
	    direction[i][j] = 'd';
	    if ( i_sim > max_sim ){
		max_sim = i_sim;
		direction[i][j] = 'i';
	    }
	    if ( j_sim > max_sim ) {
		max_sim = j_sim;
		direction[i][j] = 'j';
	    }

	    F[i][j] = max_sim;
	    /* printf ("**  %4d  %4d  %8.3f    %c \n", i, j, */
            /*       F[i][j], direction[i][j]); */
	    
	}
    }
    
    /*retrace*/
    i = max_i;
    j = max_j;
    while ( i>0 ||  j >0 ) {
	//printf (" %4d  %4d  %8.3f  \n", i, j, F[i][j]);
	switch ( direction[i][j] ) {
	case 'd':
	    //printf ( " %4d  %4d \n",  i, j);
	    map_i2j [i-1] = j-1;
	    map_j2i [j-1] = i-1;
	    i--;
	    j--; 
	    break;
	case 'i':
	    //printf ( " %4d  %4d \n",  i, -1);
	    map_i2j [i-1] = FAR_FAR_AWAY;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = FAR_FAR_AWAY;
	    j--; 
	    break; 
	default: 
	    fprintf ( stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
    return 0; 
   
    
}
/////////////////////////////////////////////////////////////
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
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(int *) calloc( rows*columns, sizeof(int));
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

void * emalloc(int  size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}





/* free a  matrix  */
void free_cmatrix(char **m)
{
    free(m[0]);
    free(m);
}
void free_imatrix(int **m)
{
    free(m[0]);
    free(m);
}
void free_dmatrix(double **m)
{
    free(m[0]);
    free(m);
}
