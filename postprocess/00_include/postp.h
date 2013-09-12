# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include "postp_pdb.h"
# include "postp_alignment.h"
# include "postp_geometry.h"
# include "postp_utils.h"
# include "postp_node.h"

# define VERSION "1.0"

/*****************************/
/*   log   modes             */
/*****************************/
# define INTRO 1
# define WARN  2
# define STATUS  3
# define NOTE  4

# define BUFFLEN  250
/******************************/
/*   tokenizer                */
/******************************/

# define TOK_TOOMNY  1 /* tokenazier error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define LONGSTRING  250
# define MEDSTRING  100
# define SHORTSTRING  25

int tokenize ( char  token[MAX_TOK][], int * max_token, char * line, char comment_char);
/******************************/
/*   geometry                 */
/******************************/
/* distance to be called a nighbor in the protein structure */
# define CUTOFF_DIST 4.0

/******************************/
/*   user options:            */
/******************************/
/* please update echo options int the logger */
typedef struct {
    char pdbname   [BUFFLEN];
    char almtname  [BUFFLEN];
    char **namefile;
    char query     [BUFFLEN];
    char special   [BUFFLEN];
    char outname   [BUFFLEN];
    char scorename [BUFFLEN]; /* file with the score already provided */
    int  scoring_method;
    int  no_of_alignments;
    char chain;
    double max_gaps;
} Options; 

/* scoring methods */
# define ENTROPY 1
# define IVET 2
# define RVET 3

/******************************/
/* function declarations :    */
/******************************/
double area_over_coverage (int * int_cvg, double * value, int no_of_res );
double avg_dist_to_special ( Options * options, Alignment * alignment);
int  clustering ( Protein *protein,  int* res_rank, int * int_cvg, double *clustering_score);
int  coverage ( Protein * protein, int * almt2prot, double * score, int almt_length, int * res_rank, int * int_cvg );
int  determine_adj_matrix ( int ** adj_matrix, Residue * sequence, int no_res, double cutoff_dist);
int  logger (Options * options, int mode, char *msg );
int  output_clustering (char * base_filename,  int * int_cvg, double *clustering_score, int length);
int  output_score ( char * base_filename, Protein * protein, Alignment * alignment,
		    double * score, int * almt2prot, int *res_rank);
int  print_tree (FILE * fptr,Node * node );
int  read_pdb      ( char * pdbname, Protein   * protein,  char chain);
int  read_clustalw ( char * cwname, Alignment * alignment);
int  read_cmd_file (char *filename, Options * options);
int  read_score ( char *filename, Protein * protein, int * score2prot, double *score );
int  reconstruct_alignment (char * namesfile, Alignment * big_almtptr, Alignment * almtptr);
int  scoring ( Options *optins, Alignment * alignment, double * score);
int  score2rank ( double * score, int * rank_order, int length ) ;
int  seq_pw_dist (Alignment * alignment) ;
char single_letter ( char code[]);
double  spearman ( int * rank1, int * rank2, int length );
int  struct_almt_mapping (Protein * protein, Alignment * alignment, char *query_name,  int * prot2almt, int * almt2prot);

# endif
