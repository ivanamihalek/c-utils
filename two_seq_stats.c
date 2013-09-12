/************************************************/
/*      Ivana, BII, Singapore, 2007             */
/************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>

# define DEFAULT_MAX_ID 0.98
# define DEFAULT_MIN_ID 0.30
# define DEFAULT_MIN_QRY_ID 0.30
# define BUFFLEN  150
# define ALMT_NAME_LENGTH 50

# define WINDOW_LENGTH 20

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int * seq_gaps;
    int * column_gaps;
    int * marked;
}  Alignment;

/************************************************/
int main ( int argc, char * argv[]) {

    char infile[BUFFLEN]= {'\0'};
    Alignment alignment;
    int read_afa (char *filenname, Alignment * alignment);
    int process_almt (Alignment * alignment);
   
    if ( argc < 2 ) {
	printf ("Usage: %s  <afa file> \n", argv[0]);
	exit (1);
    }
    sprintf (infile, "%s", argv[1]);

    if ( read_afa( infile, &alignment) ) exit (1);
    process_almt (&alignment);
    
     
    return 0;
}



/************************************************/
int process_almt (Alignment * alignment){

    int almt_size;
    int conserved_columns;
    int length[2], longer_length;
    double  pct_conserved_columns;
    /***************************/
    int count_gaps (Alignment * alignment);
    int count_conserved_columns (Alignment * alignment);
    
    almt_size = alignment->length*alignment->number_of_seqs;
    count_gaps (alignment);
    
    length[0] = alignment->length - alignment->seq_gaps[0];
    length[1] = alignment->length - alignment->seq_gaps[1];

    if ( length[0] > length[1] ) {
	longer_length = length[0];
    } else {
	longer_length = length[1];
    }
    
    conserved_columns = count_conserved_columns (alignment);
    pct_conserved_columns =
	(double)conserved_columns/longer_length;

    

    printf (" pct_cons_columns    %8.2lf\n", pct_conserved_columns);
     
    
    return 0; 
    
} 

# define ASCII 128 

/*************************************************/
int count_similar_columns (Alignment * alignment){

    char val;
    int s, c;
    int similar, similar_count = 0;
    static int similarity_set = 0;
    static char similar_to[ASCII] = {'\0'};;

    if ( !similarity_set) {
	similarity_set = 1;
	int i;
	for(i=0;i<ASCII;i++){
	    similar_to[i] = i;
	}
	similar_to['I'] = 'V';
	similar_to['L'] = 'V';
	similar_to['S'] = 'T';
	similar_to['D'] = 'E';
	similar_to['K'] = 'R';
	similar_to['Q'] = 'N';
	similar_to['.'] = '.';

	similar_to['A'] = 'V';
	similar_to['M'] = 'V';
	similar_to['G'] = 'V';
	similar_to['F'] = 'Y';
	similar_to['H'] = 'R';	
    }

    for ( c=0; c<alignment->length; c++) {
	similar = 1;
	val = alignment->sequence[0][c];
	if( val == '.' ) continue;
	for ( s=1; s<alignment->number_of_seqs; s++ ) {
	    if  ( similar_to[(int)alignment->sequence[s][c]] !=
		  similar_to[(int)val] ) {
		similar = 0;
		break;
	    }
	}
	similar_count += similar;
	    
    }

    return similar_count;
}

/*************************************************/
int count_conserved_columns (Alignment * alignment){

    char val;
    int s, c;
    int conserved, conserved_count = 0;

    for ( c=0; c<alignment->length; c++) {
	conserved = 1;
	val = alignment->sequence[0][c];
	if( val == '.' ) continue;
	for ( s=1; s<alignment->number_of_seqs; s++ ) {
	    if  ( alignment->sequence[s][c] != val ) {
		conserved = 0;
		break;
	    }
	}
	conserved_count += conserved;
	    
    }

    return conserved_count;
}


/************************************************/
int count_gaps (Alignment * alignment) {

    int s, c, total_gaps = 0;
    void * emalloc(int	size);
    alignment->seq_gaps    = (int *) emalloc (alignment->number_of_seqs*sizeof(int));
    if (!alignment->seq_gaps) return 1;
    alignment->column_gaps = (int *) emalloc (alignment->length*sizeof(int));
    if (!alignment->column_gaps) return 1;
    for ( s=0; s<alignment->number_of_seqs; s++ ) {
	for ( c=0; c<alignment->length; c++) {
	    if ( alignment->sequence[s][c] == '.' ) {
		total_gaps ++;
		alignment->column_gaps[c] ++;
		alignment->seq_gaps[s] ++;
	    }
	}
    }
    return total_gaps;
}

/************************************************/
int read_afa (char *filename, Alignment * alignment) {

    FILE * fptr = NULL;
    char line[BUFFLEN];
    int number_of_seqs, almt_length, ctr;
    int reading, seq_ctr;
    int seq_pos, begin_name;
    int * marked = NULL;
    char ** sequence;
    char ** name;
    char **chmatrix(int rows, int columns);
    FILE * efopen(char * name, char * mode);
    void * emalloc(int	size);
    /* open file */
    fptr = efopen ( filename, "r");
    if ( !fptr ) return 1;
    
    /* find the alignment length  and number of seqs info */
    almt_length = 0;
    number_of_seqs = 0;
    reading = 0;
    ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0] == '>' ){
	    number_of_seqs++;
	    reading = (! almt_length);
	} else if ( reading ) {
	    ctr =0;
	    while  ( ctr < BUFFLEN && line[ctr] != '\n' )  {
		if ( !isspace (line[ctr]) ) almt_length ++;
		ctr ++;
	    }
	}
    }
   
     /* allocate */
    sequence = chmatrix (number_of_seqs, almt_length);
    if ( !sequence ) return 1;
    name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
    if ( !name ) return 1;
    marked = (int * )emalloc ( number_of_seqs*sizeof(int));
    if ( !marked ) return 1;
    
    
    /* read in */
    rewind(fptr);
    seq_ctr = -1;
    seq_pos = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0] == '>' ){
	    seq_pos = 0;
	    seq_ctr++;
	    
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
	    if ( strlen ( &line[ctr]) > ALMT_NAME_LENGTH ) {
		line[ALMT_NAME_LENGTH+ctr-1] = '\0';
	    }
	    sprintf ( name[seq_ctr], "%s", &line[begin_name]);
	} else  {
	    ctr =0;
	    while  ( ctr < BUFFLEN && line[ctr] != '\n' )  {
		if ( !isspace (line[ctr]) ) {
		    if ( seq_pos >= almt_length ) {
			fprintf (stderr, "Error: sequence %s longer than the first sequence.\n",
				 name[seq_ctr]);
			exit (1);
		    }
		    if (line[ctr] == '-' ) {
			sequence[seq_ctr][seq_pos] = '.';
		    } else {
			sequence[seq_ctr][seq_pos] = line[ctr];
		    }
		    seq_pos ++;
		}
		ctr ++;
	    }
	}
    }

    
   
    /* return values */
    alignment->number_of_seqs = number_of_seqs;
    alignment->length         = almt_length;
    alignment->sequence       = sequence;
    alignment->name           = name;
    alignment->marked           = marked;

    
    fclose (fptr);
    
    return 0;
}


/************************************************/
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

/************************************************/
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



/************************************************/
void * emalloc(int	size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}


