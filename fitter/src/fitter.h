# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>
# include <assert.h>
# include "fitter_pdb.h"


# define  BUFFLEN 150

# define TOK_TOOMNY  1 /* tokenazier error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25

# include "fitter_utils.h"

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int layer;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 100

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+1];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    char chain;
    int no_atoms;
    int is_hetatm;
    Atom  atom[MAX_NO_ATOMS];
} Residue;

typedef struct {
    int length;
    Residue * sequence;
} Protein;


/***********************************/


int almt_out (FILE *fptr, Protein * protein1, Protein *protein2,
	      int * residue_map_i2j, int * residue_map_j2i);
int find_Calpha  (Protein *protein, int  resctr, double ca[3] );
int map2rotation (Protein *protein1, Protein *protein2, int *map_i2j,
		   double *q, double *T, double *rmsd);
int quat_to_R ( double *q, double **R );
int read_list ( char * filename, char *** list, int * length) ;
int read_pdb ( char * pdbname, Protein * protein, char chain);
int pdb_output (char *pdb_fnm_1, char *pdb_fnm_2, char * out_fnm,
		 double  **tfm_matrix, double * transl_vector, Residue * sequence, int no_res);

# endif
