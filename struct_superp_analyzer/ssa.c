# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include "pdb.h"
# include "utils.h"

# define TOK_TOOMNY  1 /* tokenazier error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define  LONGSTRING  250
# define  MEDSTRING  100
# define  SHORTSTRING  25
# define ALMT_NAME_LENGTH 30
# define BUFFLEN  250

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int * seq_gaps;
    int * column_gaps;
}  Alignment;

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 100 /* so I can handle things like heme */

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+2];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
    int interface;
} Residue;

typedef struct {
    int length;
    Residue * sequence;
} Protein;

int tokenize ( char token[MAX_TOK][MEDSTRING], int * max_token, char * line , char comment_char);
FILE * efopen(char * name, char * mode);
int read_clustalw ( char * cwname, Alignment * alignment);

int main ( int argc, char * argv[]) {

    Alignment alignment;
    Protein * protein;
    
    
    if ( argc <2 ) {
	fprintf ( stderr, "Usage: %s  <msf file>\n", argv[0]);
	exit (1);
    }

    /* read in the alignment */
    read_clustalw ( argv[0], &alignment);

    /* assign space for as amny structures as there are names in the alignment */
    protein = emalloc ( alignment.number_of_seqs*sizeof(Protein) );
    /*  for each name in the alignment*/
    /* read in the pdb file */
    /* map the pdb to the alignment */

    /*for each position in the alignment, find rmsd for c-alphas and gap entropy */

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

/*****************************************/
int read_clustalw ( char * cwname, Alignment * alignment){
    
    FILE * fptr = NULL;
    char line[BUFFLEN];
    int  number_of_seqs, almt_length, ctr;
    int * seq_pos, pos;
    char ** sequence;
    char ** name;
    char curr_name[BUFFLEN];
     
    /* open file */
    fptr = efopen ( cwname, "r");
    if ( !fptr ) return 1;
    
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
	fprintf ( stderr, "Alignment length info not found in %s. Is the format gcg?\n", cwname);
	return 1;
    }

    /* determine the number of sequences */
    number_of_seqs = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( ! strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) number_of_seqs++;
    }
    if ( number_of_seqs ) {
	/* printf ( "Number of sequences in %s is %d.\n", cwname, number_of_seqs); */
    } else {
	fprintf ( stderr, "No sequences found in %s. Is the format gcg?\n", cwname);
	return 1;
    } 
    
    /* allocate */
    sequence = chmatrix (number_of_seqs, almt_length);
    if ( !sequence ) return 1;
    name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
    if ( !name ) return 1;
    seq_pos = (int *) calloc ( number_of_seqs, sizeof(int));
    if ( !seq_pos ) return 1;
    
    /* read in */
    rewind(fptr);
    ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if (!  strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    sscanf (line, "%*s %s", name[ctr]);
	    ctr ++;
	}
    }
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( isspace (line[0] ) ) continue;
	sscanf (line, "%s", curr_name);
	ctr = 0;
	while (  ctr <number_of_seqs &&  strcmp (name[ctr], curr_name) ) ctr++;
	if ( ctr >= number_of_seqs ) {
	    fprintf ( stderr, "The name %s not found in the header of %s.\n", curr_name, cwname);
	    return 1;
	}
	pos = 0;
	while ( ! isspace(line[pos]) ) pos++;
	while  (line[pos] != '\n' && pos < BUFFLEN) {
	    if ( !  isspace(line[pos] ) ){
		if ((line[pos]>=97)&&(line[pos]<=122)) {line[pos] -= 32;} /* --> turn to uppercase */
		if ( line[pos]==126)                   {line[pos]  = 46;} /* turn tweedle to dot */
		sequence [ctr] [ seq_pos[ctr] ] = line[pos];
		seq_pos[ctr]++;
	    }
	    pos ++;
	}
    }
    fclose(fptr);

    /* sanity check */
    for (ctr=0; ctr < number_of_seqs; ctr++ ) {
	if ( seq_pos[ctr] != almt_length ) {
	    fprintf (stderr, "Sequence %s is shorter (%d position) than the alignment.\n", name[ctr],  seq_pos[ctr]);
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

    { 
	int count_gaps (Alignment * alignment);
	count_gaps (alignment);
    }
    return 0;
}


