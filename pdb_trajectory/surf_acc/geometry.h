#include "sac.h"

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    int pdb_id;
    double x,y,z;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 100 /* so I can handle things like heme */

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+2];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    int  no_atoms;
    Atom  atom[MAX_NO_ATOMS];
    int interface;
} Residue;

typedef struct {
    int length;
    Residue * sequence;
} Protein;
