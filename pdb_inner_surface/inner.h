# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>
# include <assert.h>

#include "inner_pdb.h"


# define  BUFFLEN 150

# define TOK_TOOMNY  1 /* tokenizer error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25

# include "inner_utils.h"

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double coord[3];
    double rho, theta, z;
    int layer;
    int backbone;
    void * parent_residue;
} Atom;

# define  MAX_NO_ATOMS 400

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+1];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    char chain;
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
} Residue;

typedef struct {
    int length;
    Residue * residue;
    Atom** calpha;
} Protein;


/***********************************/
int quat_to_R ( double *q, double **R );
int read_list ( char * filename, char *** list, int * length) ;
int read_pdb ( char * pdbname, Protein * protein, char chain);
int find_principal_axes (Protein *protein, double principal_moment[3], double principal_axis[3][3]);
int max_index (double vec[3]);
int min_index (double vec[3]);
int rotate_structure(Protein *protein, double new_z_direction[3]);

int coords2cylindrical(Protein * protein);
int bin_atoms(Protein *protein, double z_step, int  number_of_theta_bins, Atom **** bin_ptr, int **bin_size_ptr, int *number_of_z_bins_ptr);
int is_pointing_inward(Residue * res);
int distribution_of_rho (Protein* protein, double *avg_ptr, double *stdev_ptr);

# endif









