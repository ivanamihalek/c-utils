# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>
# include <assert.h>
# include "smallig_pdb.h"


# define  BUFFLEN 150

# define TOK_TOOMNY  1 /* tokenizer error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25

# include "smallig_utils.h"

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double coord[3];
    int layer;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 400

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+1];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
} Residue;

typedef struct {
    int length;
    Residue * sequence;
} Protein;


/***********************************/


int almt_out (FILE *fptr, Protein * protein1, Protein *protein2,
	      int * residue_map_i2j, int * residue_map_j2i);
int find_Calpha ( Protein *protein, int  resctr, double ca[3] );
int map2rotation (Protein *protein1, Protein *protein2, int *map_i2j,
		   double *q, double *T, double *rmsd);
int quat_to_R ( double *q, double **R );
int read_list ( char * filename, char *** list, int * length) ;
int read_pdb ( char * pdbname, Protein * protein, char chain);


# endif
