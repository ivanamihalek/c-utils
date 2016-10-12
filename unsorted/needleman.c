/** SSE alignment **/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>


#define  GAP_OPENING    -0.6
#define  GAP_EXTENSION  -0.3
#define  X_SIM           0.0
#define  ENDGAP          0.0
# define BUFFLEN  150
# define NAME_LENGTH 50


#define  FAR_FAR_AWAY  -10000

double  needleman_wunsch (int max_i, int max_j, double **similarity,
		      int *map_i2j, int * map_j2i, int use_endgap);
char   **chmatrix(int rows, int columns);
double **dmatrix(int rows, int columns);
void *   emalloc(int	size);
void     free_cmatrix(char **m);
void     free_dmatrix(double **m);
int      read_seq    (char *filename, char ** seq, int *len, char * name);
int      write_afa   (char *filename, char ** seq, int length, char name[2][NAME_LENGTH] );

int main ( int argc, char * argv[]) {

    double ** similarity;
    int * map_i2j, * map_j2i;
    int len1 ;
    int len2;
    int i, j;
 
    double sim_matrix[4][4] = /* the assumed order: HSCX */
	{ {  1.5, -1, -1, X_SIM},
	  { -1,  1, -1, X_SIM},
	  { -1, -1,  1, X_SIM},
	  { X_SIM, X_SIM, X_SIM,  1}};
    
    int sse_index[200];
    int sse_index_i, sse_index_j;
    int use_endgap = 0;
    int almt_pos = 0;
    int gap_open, aln_length;
    char * seq[2] = {NULL};
    char ** aligned_seq;
    char   name[2][NAME_LENGTH] = {{'\0'}};
    char filename[NAME_LENGTH]  = {'\0'};
    double score;
   
    if ( argc < 3 ) {
	fprintf (stderr, "Usage: %s <sse file 1> <sse file 2> [endgap] \n", argv[0]);
	exit (-1);
    }
    

    if ( argc == 4 && ! (strncmp (argv[3], "end", 3 ) ) ) {
	use_endgap = 1;
    }

    if ( read_seq (argv[1], seq, &len1, name[0]) )  return -1;
    if ( read_seq (argv[2], seq+1, &len2, name[1]) )  return -1;

    if ( !(aligned_seq=chmatrix(2,len1+len2 )))  return -1;
    
 
    if ( !(similarity = dmatrix (len1, len2) )) return -1;
    if ( !(map_i2j    = emalloc (len1*sizeof(int)))) return -1;
    if ( !(map_j2i    = emalloc (len2*sizeof(int)))) return -1;
    

    for (i=0; i < 200; i++ ) sse_index[i] = -1;
    sse_index['H'] = sse_index['h'] = 0;
    sse_index['S'] = sse_index['s'] = 1;
    sse_index['C'] = sse_index['c'] = 2;
    sse_index['X'] = sse_index['x'] = 3;
    
    for (i=0; i < len1; i++ ) {
	sse_index_i = sse_index[ (int) seq[0][i]  ];
	if ( sse_index_i < 0 ) {
	    fprintf (stderr, "Unrecognized nucleotide type: %c in %s @ pos %d\n",
		     seq[0][i], name[0], i);
	    exit (1);
	}
	for (j=0; j < len2; j++ ) {
	    sse_index_j = sse_index[ (int)  seq[1][j]  ];
	    if ( sse_index_j < 0 ) {
		fprintf (stderr, "Unrecognized nucleotide type: %c in %s @ pos %d\n",
			 seq[1][j], name[1], j);
		exit (1);
	    }
	    similarity[i][j] =  sim_matrix[sse_index_i][sse_index_j];
	}
    }
  
    score = needleman_wunsch (len1, len2, similarity, map_i2j, map_j2i, use_endgap);


    
    i = 0; j = 0;
    almt_pos = 0;
    while (i< len1 ||  j < len2 ) {
	gap_open = 0;
	while (i< len1 &&  map_i2j[i] == FAR_FAR_AWAY) {
	    aligned_seq[0][almt_pos] = seq[0][i];
	    aligned_seq[1][almt_pos] = '.';
	    almt_pos ++;
	    i ++;
	}
	gap_open = 0;
	while (j < len2 && map_j2i[j] == FAR_FAR_AWAY) {
	    aligned_seq[0][almt_pos] = '.';
	    aligned_seq[1][almt_pos] = seq[1][j];
	    almt_pos ++;
	    j ++;
	}
	if (i< len1 &&  j < len2 )  {
	    if ( map_i2j[i] != j || map_j2i[j] != i ) {
		fprintf (stderr, "alignment error (?)\n");
		exit (1);
	    }
	    aligned_seq[0][almt_pos] = seq[0][i];
	    aligned_seq[1][almt_pos] = seq[1][j];
	    almt_pos ++;
	    i ++;
	    j ++;
	}
	
    }
     
    aligned_seq[0][almt_pos] = '\0';
    aligned_seq[1][almt_pos] = '\0';
    aln_length = almt_pos;

    sprintf   ( filename, "%s_%s.sse_afa", name[0], name[1]);
    write_afa ( filename, aligned_seq, aln_length, name);
    
    printf ( "score  %8.2lf\n", score);
    
    free_dmatrix (similarity);
    free ( map_i2j);
    free ( map_j2i);
	 
    
    return 0;
}

/////////////////////////////////////////////////
int write_afa (char *filename, char ** seq, int length, char name[2][NAME_LENGTH] ) {

    FILE * fptr = NULL;
    int seq_ctr, pos_ctr, c;
    FILE * efopen(char * name, char * mode);
    
    /* open file */
    fptr = efopen ( filename, "w");
    if ( !fptr ) return 1;

    
    for (seq_ctr=0; seq_ctr < 2; seq_ctr++ ) {

	fprintf (fptr,  ">%s\n", name[seq_ctr]);
	
	c = 0; /* counter for newline insertion */ 
	for (pos_ctr=0; pos_ctr < length; pos_ctr++) {
	    c++;  if ( c && !(c%50) ) fprintf (fptr, "\n");
	    fprintf ( fptr, "%c", seq[seq_ctr][pos_ctr] );
	}
	fprintf (fptr, "\n");
    }
 
    fclose (fptr);
    return 0;
}
/////////////////////////////////////////////////
int read_seq (char *filename, char ** seq_ptr, int *len, char * name ) {
    
    FILE * fptr = NULL;
    char * seq;
    char line[BUFFLEN];
    int length, ctr;
    int seq_pos, begin_name;

    char **chmatrix(int rows, int columns);
    FILE * efopen(char * name, char * mode);
    void * emalloc(int	size);

    /* open file */
    fptr = efopen ( filename, "r");
    if ( !fptr ) return 1;
    
    /* find the alignment length */
    length = 0;

    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0] == '>' ) {
	} else {
	    ctr =0;
	    while  ( ctr < BUFFLEN && line[ctr] != '\n' )  {
		if ( !isspace (line[ctr]) ) length ++;
		ctr++;
	    }
	}
    }

    /* allocate */
    seq = emalloc (length*sizeof(char));
    if ( !seq ) return 1;
    
    /* read in */
    rewind(fptr);
    seq_pos = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0] == '>' ){
	    seq_pos = 0;
	    /* chomp*/
	    ctr = 1;
	    while ( ctr < BUFFLEN && line[ctr] != '\n' ) ctr++;
	    if ( ctr < BUFFLEN ) line[ctr] = '\0';
	    
	    ctr = 1;
	    while ( isspace (line[ctr]) ) ctr++;
	    begin_name = ctr;
	    while ( !isspace (line[ctr]) ) ctr++;
	    line[ctr] = '\0';
	    /* make sure the name is not too long */
	    if ( strlen ( &line[ctr]) > NAME_LENGTH ) {
		line[NAME_LENGTH+ctr-1] = '\0';
	    }
	    sprintf ( name, "%s", &line[begin_name]);
	    
	} else  {
	    
	    ctr =0;
	    while  ( ctr < BUFFLEN && line[ctr] != '\n' )  {
		if ( !isspace (line[ctr]) ) {
		    if (line[ctr] == '-' ) {
			seq[seq_pos] = '.';
		    } else {
			seq[seq_pos] = line[ctr];
		    }
		    seq_pos ++;
		}
		ctr ++;
	    }
	}
    }
    
  
    printf ( "%s length = %d\n", name, length);
    
    *len = length;
    *seq_ptr = seq;
    fclose (fptr);
    
    return 0;
    
}


/////////////////////////////////////////////////
double needleman_wunsch (int max_i, int max_j, double **similarity,
		      int *map_i2j, int * map_j2i, int use_endgap) {

    double **F; /*alignment_scoring table*/
    char ** direction;
    double gap_opening = GAP_OPENING;
    double gap_extension = GAP_EXTENSION;
    double endgap = ENDGAP;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    double score = 0.0;
    int i,j;

     /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;

    /* fill the table */
    for (i=0; i<= max_i; i++) {
	for (j=0; j<=max_j; j++) {

	    if ( !i && !j ) { /* upper left corner */
		F[0][0] = 0;
		direction[i][j] = 'd';
		continue;
	    }
	    
	    if ( i && j ){ 
		if ( direction[i-1][j] == 'i' ) {
		    /*  gap extension  */
		    penalty = (use_endgap&&j==max_j) ? endgap : gap_extension;		    
		} else {
		    /*  gap opening  */
		    penalty = (use_endgap&&j==max_j) ? endgap : gap_opening;
		}
		i_sim =  F[i-1][j] + penalty;

		
		if ( direction[i][j-1] =='j' ) {
		    penalty = (use_endgap&&i==max_i) ? endgap : gap_extension;		    
		} else {
		    penalty = (use_endgap&&i==max_i) ? endgap : gap_opening;		    
		}
		j_sim = F[i][j-1] +  penalty;
		
		
		diag_sim =  F[i-1][j-1] + similarity [i-1][j-1] ;
		
	    } else if ( j ) {
		
		if ( use_endgap) {
		    penalty = endgap;
		} else {
		    if ( direction[i][j-1] =='j' ) {
			penalty =  gap_extension;
		    } else {
			penalty =  gap_opening;
		    }
		}
		j_sim = F[i][j-1] + penalty;
		
		i_sim = diag_sim =  FAR_FAR_AWAY;

	    } else if ( i ) {
		if ( use_endgap) {
		    penalty = endgap;
		} else {
		    if ( direction[i-1][j] == 'i' ) {
			penalty =  gap_extension;
		    } else {
		        penalty =  gap_opening;
		    }
		}
		i_sim = F[i-1][j] + penalty;
		
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
	    
	}
    }

    score = F[max_i][max_j];

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
    
    return score; 
   
    
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

//////////////////////////////////////////////////////
FILE * efopen(char * name, char * mode) {

    FILE * fp;


    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,  
	      "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }

    return fp;

}
