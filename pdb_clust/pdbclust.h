# ifndef UMBRELLA_H
# define UMBRELLA_H

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "pdb.h"
# include "utils.h"

# define PANIC(msg, str) {			\
    fprintf (stderr, "%s %s.\n",msg, str);	\
    exit(1);					\
}

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
} Atom;
# define  MAX_NO_ATOMS 100
typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+1];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
} Residue;

Residue * sequence1, *sequence2;
int no_res_1, no_res_2;

typedef struct {
    char aa;
    double distance;
} Point;


extern char *amino_acid_order;
# define AA_ORDER_LENGTH 21


void cluster_counter (int  no_of_things,  int *neighbors[], int * mask,
		      int cluster_count_per_size[], int * no_of_clusters,
		      int * max_size, int * secnd_max_size , int * clusters[]);
int cluster_score (int no_of_res, int *seq, int ** adj_matrix,double *score);
int determine_dist_matrix ( double ** dist_matrix, Residue * sequence, int no_res);
int read_pdb (char * pdbname, char *chain_id_ptr, Residue ** sequence_ptr, int * no_res_ptr);

int read_residue_selection (char * selected_res_file, Residue * sequence, int no_res, int *selected);
int std_dev_over_S (int L, int M, int ** adj_matrix, double *avg, double * std_dev, int first);

# endif
