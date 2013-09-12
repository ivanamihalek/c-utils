#include "pdb.h"

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int backbone;
    int clash;
    int pdb_id;
} Atom;

# define  MAX_NO_ATOMS 1000 /* so I can handle things like heme */

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
