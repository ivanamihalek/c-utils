
#ifndef PDB_H
# define	PDB_H

# include	<stdio.h>

#define        PDB_ATOM_TEMP_LEN           6  
    

/* 73 - 76        LString(4)      segID         Segment identifier, left-justified. */
#define        PDB_ATOM_SEGID             72
#define        PDB_ATOM_SEGID_LEN          4
    

/* 77 - 78        LString(2)      element       Element symbol, right-justified. */
#define        PDB_ATOM_ELMT              76
#define        PDB_ATOM_ELMT_LEN           2
     

/* 79 - 80        LString(2)      charge        Charge on the atom. */
#define        PDB_ATOM_CHARGE            78
#define        PDB_ATOM__CHARGE_LEN        2  
    

/* field start and length in atom record: */
    
/*  1 -  6        Record name     "ATOM  " */
# define        PDB_ATOM_REC_NAME          0
# define        PDB_ATOM_REC_NAME_LEN      6
/*     7 - 11        Integer         serial        Atom serial number. */
#define         PDB_ATOM_ATOM_NO           6       
#define         PDB_ATOM_ATOM_NO_LEN       5
  
/* 13 - 16        Atom            name          Atom name. */
#define         PDB_ATOM_ATOM_NAME         12
#define         PDB_ATOM_ATOM_NAME_LEN     4     
    

/* 17             Character       altLoc        Alternate location indicator. */
#define         PDB_ATOM_ALTLOC            16  
#define         PDB_ATOM_ALTLOC_LEN        1  

/* 18 - 20        Residue name    resName       Residue name. */
# define        PDB_ATOM_RES_NAME          17
# define        PDB_ATOM_RES_NAME_LEN      3

/* 22             Character       chainID       Chain identifier. */
#define        PDB_ATOM_CHAINID            21
#define        PDB_ATOM_CHAINID_LEN         1

/* 23 - 26        Integer         resSeq        Residue sequence number */
#define        PDB_ATOM_RES_NO             22
#define        PDB_ATOM_RES_NO_LEN         4

/* 27             AChar           iCode         Code for insertion of residues. */
#define        PDB_ATOM_INS_CODE           26
#define        PDB_ATOM_INS_CODE_LEN       1
    
/* 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in */
/*                                              Angstroms. */
#define        PDB_ATOM_X                  30
#define        PDB_ATOM_X_LEN               8
    
/* 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in */
/*                                              Angstroms. */
#define        PDB_ATOM_Y                  38
#define        PDB_ATOM_Y_LEN               8
    

/* 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in */
/*                                              Angstroms. */
#define        PDB_ATOM_Z                  46
#define        PDB_ATOM_Z_LEN               8
    


/* 55 - 60        Real(6.2)       occupancy     Occupancy. */
#define        PDB_ATOM_OCC                54
#define        PDB_ATOM_OCC_LEN             6
    

/* 61 - 66        Real(6.2)       tempFactor    Temperature factor. */
#define        PDB_ATOM_TEMP              60



/* binding capabilities */
# define NO_BINDING_TYPES 6
# define NONE       0
# define VDW        2
# define H_DONOR    1
# define H_ACCEPTOR 5
# define POSITIVE   2
# define NEGATIVE   4

# define MAX_BB_NBRS 5
# define MAX_NBRS    1


typedef struct {
    int no_atom;
    char name_atom[5];
    char name_res[4];
    char chain_id;
    char  no_res[6];
    double x;
    double y;
    double z;
    int no_of_neighbors;
    double dist_to_cm;
} Pdb_map;


typedef struct {
    char id;
    int start;
    int end;
    int number_of_residues;
    int *residue_start;
    int *residue_end;
} Chain_info;

#define MAX_CHAINS 30

# endif
