# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>

# define DEFAULT_MAX_ID 0.99
# define BUFFLEN  150
# define ALMT_NAME_LENGTH 50

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

    char infile[BUFFLEN]= {'\0'}, outfile[BUFFLEN]= {'\0'};
    char query[ALMT_NAME_LENGTH];
    double max_id;
    Alignment alignment;
    int read_afa (char *filenname, Alignment * alignment);
    int compare_seqs (Alignment * alignment, char *query, double max_id);
    int write_afa (char *filename, Alignment * alignment);
   
    if ( argc < 3 ) {
	printf ("Usage: %s  <infile>  <outfile> [<query> [<max_id>]].\n", argv[0]);
	exit (1);
    }
    sprintf (infile, "%s", argv[1]);
    sprintf (outfile, "%s", argv[2]);
    max_id = DEFAULT_MAX_ID;
    if (argc > 3 ) sprintf (query, "%s", argv[3]);
    if (argc > 4 ) max_id = atof(argv[4]);


    if ( read_afa( infile, &alignment) ) exit (1);
    compare_seqs (&alignment, query, max_id);
    if ( write_afa ( outfile, &alignment) ) exit (1);

    
     
    return 0;
}



/************************************************/
int compare_seqs (Alignment * alignment, char *query, double max_id){
    
    int seq_ctr_1,  seq_ctr_2, pos_ctr;
    int common_length, identical_length;
    int is_query;
    char * seq1,  *seq2;
    double pid;
    int count_gaps (Alignment * alignment);

    if ( count_gaps(alignment ) ) exit(1);

    is_query = -1;
    for (seq_ctr_1=0; seq_ctr_1 < alignment->number_of_seqs; seq_ctr_1++ ) {
	if ( ! strncmp (alignment->name[seq_ctr_1], query, strlen(query) ) &&
	     ! strncmp (alignment->name[seq_ctr_1], query, strlen(alignment->name[seq_ctr_1]) ) ) {
	    is_query = seq_ctr_1;
	    break;
	}
    }
   
    for (seq_ctr_1=0; seq_ctr_1 < alignment->number_of_seqs; seq_ctr_1++ ) {
	if ( alignment->marked[seq_ctr_1] ) continue;
	seq1 = alignment->sequence[seq_ctr_1];
	for (seq_ctr_2=seq_ctr_1+1; seq_ctr_2 < alignment->number_of_seqs; seq_ctr_2++ ) {
	    if ( alignment->marked[seq_ctr_2] ) continue;
	    seq2 = alignment->sequence[seq_ctr_2];
	    common_length = identical_length = 0;
	    for (pos_ctr=0; pos_ctr < alignment->length; pos_ctr++) {
		if ( seq1[pos_ctr] == '.' || seq2[pos_ctr] == '.' ) continue;
		common_length ++;
		if ( seq1[pos_ctr] ==  seq2[pos_ctr] ) identical_length++;
	    }
	    pid = common_length ? (double)identical_length/common_length:0.0;
	    if ( pid > max_id ) {
		printf ( " %s  %s  %8.2lf \n", alignment->name[seq_ctr_1],
			 alignment->name[seq_ctr_2], pid );
		if (query &&  seq_ctr_1 == is_query) {
		    alignment->marked[seq_ctr_2] = 1;
		} else if ( query &&  seq_ctr_2 == is_query) {
		    alignment->marked[seq_ctr_1] = 1;
		    break; /* no longer interested in comparinf to seq1 */
		} else if ( alignment->seq_gaps[seq_ctr_1] >= alignment->seq_gaps[seq_ctr_2] ) {
		    alignment->marked[seq_ctr_1] = 1;
		    break; /* no longer interested in comparinf to seq1 */
		} else {
		    alignment->marked[seq_ctr_2] = 1;
		}
	    }
	}
    }

    return 0;
}

/************************************************/
int count_gaps (Alignment * alignment) {

    int s, c;
    void * emalloc(int	size);
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

/************************************************/
int write_afa (char *filename, Alignment * alignment) {

    FILE * fptr = NULL;
    int seq_ctr, pos_ctr;
    FILE * efopen(char * name, char * mode);
    /* open file */
    fptr = efopen ( filename, "w");
    if ( !fptr ) return 1;

    for (seq_ctr=0; seq_ctr < alignment->number_of_seqs; seq_ctr++ ) {
	if ( alignment->marked[seq_ctr] ) continue;
	fprintf (fptr,  "> %s\n", alignment->name[seq_ctr]);
	for (pos_ctr=0; pos_ctr < alignment->length; pos_ctr++) {
	    if ( pos_ctr && !(pos_ctr%50) ) fprintf (fptr, "\n");
	    fprintf ( fptr, "%c", alignment->sequence[seq_ctr][pos_ctr] );
	}
	fprintf (fptr, "\n");
    }
 
    fclose (fptr);
    return 0;
}
/************************************************/
int read_afa (char *filename, Alignment * alignment) {

    FILE * fptr = NULL;
    char line[BUFFLEN];
    int number_of_seqs, almt_length, ctr;
    int reading, seq_ctr;
    int seq_pos;
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
    printf ( "number of seqs = %d\n",   number_of_seqs);
    printf ( "almt_length = %d\n", almt_length);
    
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
	    /* make sure the name is not too long */
	    if ( strlen ( &line[ctr]) > ALMT_NAME_LENGTH ) {
		line[ALMT_NAME_LENGTH+ctr-1] = '\0';
	    }
	    sprintf ( name[seq_ctr], "%s", &line[ctr]);
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


