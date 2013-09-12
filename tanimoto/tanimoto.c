# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
#include <ctype.h>

#include "pdb.h"

double covalent_radius[127] = {0.};

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int layer;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 100 /* so I can handle things like heme */
# define  DELTA  0.05  /*grid distance in angstroms */

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+2];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
    int interface;
} Residue;

Residue * sequence1, *sequence2;
int no_res_1, no_res_2;


int main ( int argc, char * argv[]) {
    

    int read_pdb ( char * pdbname, Residue ** sequence, int *no_res);
    int set_covalent();
    int tanimoto ();
    
    char pdbname1 [150] = "\0";
    char pdbname2 [150] = "\0";

    if ( argc < 3 ) {
	printf ( "Usage: %s <pdbname1> <pdbname2> .\n", argv[0] );
	exit (1);
    } 
    sprintf ( pdbname1, "%s", argv[1]);
    sprintf ( pdbname2, "%s", argv[2]);
    /* printf  ("pdbnames: %s%s\n", pdbname1, pdbname2); */
   
    
    /* read in the pdbs */
    if ( read_pdb ( pdbname1, &sequence1, &no_res_1)) exit (1);
    if ( read_pdb ( pdbname2, &sequence2, &no_res_2)) exit (1);
    
    set_covalent();

# if 0
    /* test the input*/
     int resctr;
     for ( resctr=0; resctr<no_res_2; resctr++) {
	printf ("resctr: %d     pdbid: %s    restype:  %s    no atoms: %d  \n",
		resctr, sequence2[resctr].pdb_id, sequence2[resctr].res_type,
		sequence2[resctr].no_atoms);
    }
    exit (0);
# endif
    
    /* calculate the score*/
    tanimoto ();
    return 0;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int tanimoto () {

    double overlap (Residue * sequence1, Residue *sequence2);
    double ovlp12, self1, self2, tanimoto;
    
    self1    = overlap(sequence1, sequence1);
    self2    = overlap(sequence2, sequence2);
    ovlp12   = overlap(sequence1, sequence2);
    tanimoto = ovlp12/ (self1+self1-ovlp12);

    printf ( "   self1:  %8.3lf \n", self1);
    printf ( "   self2:  %8.3lf \n", self2);
    printf ( "  ovlp12:  %8.3lf \n", ovlp12);
    printf ( "tanimoto:  %8.3lf \n", tanimoto);


    return 0;
    
}
/*******************************************************************************/

double overlap (Residue * sequence1, Residue *sequence2) {
    /* evalute the expression for the volume overlap of two densites
       which are gausians cenetered at the atom centers and with their
       Gaussian const equal to the covalent radius */ 
    int resctr1, resctr2;
    int atomctr1, atomctr2;
    Atom * atomptr1, * atomptr2;
    double x1, y1, z1 ,a1, a2, aux, dist, sum, asq, term;

    sum = 0.0;
    for ( resctr1=0; resctr1<no_res_1; resctr1++) {
	for (atomctr1=0; atomctr1 < sequence1[resctr1].no_atoms; atomctr1++) {
	    atomptr1 = sequence1[resctr1].atom + atomctr1;
	    
	    x1 = atomptr1->x;
	    y1 = atomptr1->y;
	    z1 = atomptr1->z;
	    a1 = covalent_radius[ (int)atomptr1->type[0] ];
	    if ( !a1 ) {
		/* if in doubt, use C */
		a1 = covalent_radius['C' ];

	    }
	     
	    for ( resctr2=0; resctr2<no_res_2; resctr2++) {
		for (atomctr2=0; atomctr2 < sequence2[resctr2].no_atoms; atomctr2++) { 
		    		    
		    atomptr2 = sequence2[resctr2].atom + atomctr2; 
		    a2 = covalent_radius[ (int)atomptr2->type[0] ]; 
		    if ( !a2) {
			/* if in doubt, use C */
			a2 = covalent_radius['C' ];
		
		    }
	     
		    
		    aux =  (x1 - atomptr2->x); 
		    dist = aux*aux;
		    aux =  (y1 - atomptr2->y); 
		    dist += aux*aux; 
		    aux =  (z1 - atomptr2->z); 
		    dist += aux*aux;
		    asq = (a1*a1 + a2*a2);
		    term = pow ( a1*a2/asq, 3.0)*exp ( -dist/asq);
		    /* printf ( "  %8.3lf   %8.3lf   %8.3lf  %8.3lf   %8.3lf  %8.3lf  %8.3lf \n",  */
			    /*  a1, a2, asq, dist, pow ( a1*a2/asq, 3.0), exp ( -dist/asq), term); */
		    sum += term; 
		}
	    }
	}
    }


    return sum;
    
}



#define BUFFLEN 250

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int read_pdb ( char * pdbname, Residue ** sequence_ptr, int * no_res_ptr) {

    Residue * sequence;
    FILE * fptr = NULL;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+2]; /* res name: 4 digits + insertion code + \0 */
    char oldrestype [PDB_ATOM_RES_NAME_LEN+2];
    char tmp[PDB_ATOM_X_LEN+1], *auxptr;
    int atomctr, resctr, no_res, ctr, nonblank;

    
    /* open file */
    fptr = fopen ( pdbname, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", pdbname);
	return 1;
    }

    /* count residues */
    resctr = 1;
    
    no_res = resctr;
    *no_res_ptr = no_res;
    /* printf ("no residues: %d\n", no_res); */

    /* allocate space */
    sequence = NULL;
    sequence = calloc ( no_res, sizeof (Residue));
    if ( ! sequence ) {
	fprintf ( stderr, "Error allocating sequence space.\n");
	exit (0);
    }
    *sequence_ptr = sequence;

    /* read in the atom */
    rewind ( fptr);
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2); /*  tyring to account for the insertion code */
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+2);
    resctr= 0;
    atomctr = -1;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if( ! strncmp(line,"ATOM", 4)  ||  ! strncmp(line,"HETATM", 6)){
	    atomctr ++;
	    if ( atomctr >= MAX_NO_ATOMS ) {
		fprintf ( stderr, "Error: max atoms in tanimoto: %d \n",
			  MAX_NO_ATOMS );
		exit (1);
	    }
	    sequence[resctr].no_atoms = atomctr + 1;
	    
	    /* read in atom info */
	    auxptr = line+ PDB_ATOM_ATOM_NAME;
	    memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	    /* skip initial blanks*/
	    ctr  = 0;
	    while ( !(isalpha (*(auxptr + ctr))) &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	    /* copy alphanum info */
	    nonblank = 0;
	    while (  isalpha (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
		tmp[nonblank] =  *(auxptr +ctr);
		nonblank ++;
		ctr++;
	    }
	    
	    if ( tmp[0] == 'C' && tmp[1] ) {
		tmp[0] = tmp[1];
		tmp[1] = 0;
	    } else if ( ! strncmp  (tmp, "FE", 2) ) {
		tmp[0] = 'E'; tmp[1] = 0;
	    } else if ( ! strncmp  (tmp, "SI", 2) ) {
		tmp[0] = 's'; tmp[1] = 0;
	    } else if ( ! strncmp  (tmp, "AS", 2) ) {
		tmp[0] = 'a'; tmp[1] = 0;
	    } else if ( ! strncmp  (tmp, "SE", 2) ) {
		tmp[0] = 'e'; tmp[1] = 0;
	    } else if ( ! strncmp  (tmp, "BR", 2) ) {
		tmp[0] = 'r'; tmp[1] = 0;
	    } else if ( ! strncmp  (tmp, "TE", 2) ) {
		tmp[0] = 't'; tmp[1] = 0;
	    }
	    strncpy ( sequence[resctr].atom[atomctr].type, tmp, PDB_ATOM_ATOM_NAME_LEN );

	    /* is this a backbone atom?*/
	    sequence[resctr].atom[atomctr].backbone = 0;
	    if ( nonblank == 1) {
		  sequence[resctr].atom[atomctr].backbone =
		      !(  strcmp ( tmp, "N") && strcmp ( tmp, "C") && strcmp ( tmp, "O")  );
	    } else if (  nonblank == 2) {
		  sequence[resctr].atom[atomctr].backbone = ! strcmp ( tmp, "CA" );
	    }
	    /* printf ( " %4d %4d %4s is backbone: %1d \n", resctr, atomctr, */
		     /* sequence[resctr].atom[atomctr].type, sequence[resctr].atom[atomctr].backbone); */
	    strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	    tmp[PDB_ATOM_X_LEN] = '\0';
	    sequence[resctr].atom[atomctr].x=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	    tmp[PDB_ATOM_Y_LEN] = '\0';
	    sequence[resctr].atom[atomctr].y=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	    tmp[PDB_ATOM_Z_LEN] = '\0';
	    sequence[resctr].atom[atomctr].z=atof(tmp);
	   
	}
	
    }

    /* close file */
    fclose (fptr);
    return 0;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int set_covalent (){

    covalent_radius['H'] = 	0.320;
    covalent_radius['C'] = 	0.720;
    covalent_radius['N'] = 	0.680;
    covalent_radius['O'] = 	0.860;
    covalent_radius['P'] = 	1.036;
    covalent_radius['S'] = 	1.020;
    covalent_radius['F'] = 	0.640;
    covalent_radius['A'] = 	0.992; /*Ca*/
    covalent_radius['E'] = 	1.420; /* Fe*/
    covalent_radius['Z'] = 	1.448; /*Zn*/
    covalent_radius['D'] = 	1.688; /*Cd*/
    covalent_radius['I'] = 	1.400; 
    
    covalent_radius['B'] = 	0.830; /*B, boron*/
    covalent_radius['L'] = 	0.680; /*Cl, chlorine*/
    covalent_radius['s'] = 	0.860; /* Si */
    covalent_radius['a'] = 	1.036; /* As, arsenic */
    covalent_radius['e'] = 	1.036; /* Se, selenium */
    covalent_radius['r'] = 	1.036; /* Br, bromine*/
    covalent_radius['t'] = 	1.036; /* Te*/
   
       
/* Element  	Covalent Radius  	van der Waals Radius  	United Atom Radius (includes hydrogen) */
/* H  	0.320 Angstroms  	1.100 Angstroms  	(not applicable) */
/* C 	0.720 	1.548 	1.872 */
/* N 	0.680 	1.400 	1.507 */
/* O 	0.680 	1.348 	1.400 */
/* P 	1.036 	1.880 	(not applicable) */
/* S 	1.020 	1 */
/* Ca  	0.992  	1.948  	(not applicable) */
/* Fe 	1.420 	1.948 	(not applicable) */
/* Zn 	1.448 	1.148 	(not applicable) */
/* Cd 	1.688 	1.748 	(not applicable) */
/* I 	1.400 	1.748 	(not applicable) */

    return 0;

}
