/* this subroutine read the pdb file */
/* for ATOM line format in a pdb file, see
   http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html */
#include "ifc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define LONGSTRING 300

int  read_pdb_file(char filename[], Pdb_map ** pdb_map_ptr, int *map_length_ptr, int chain_start[])

{
  int atom_ctr, i, j;
  int map_length;
  char line[LONGSTRING];
  char tmp[LONGSTRING];
  char old_chain_id;
  int chain_ctr;
  FILE * fptr;
  Pdb_map *pdb_map;

# ifdef VERBOSE
  printf("reading the pdb file ...\n");
# endif
  
  if((fptr=fopen(filename,"r"))==NULL) {
      fprintf ( stderr, "Error opening %s.\n", filename);
      return 1;
  }

  /* init */
  for (chain_ctr=0; chain_ctr<MAX_CHAINS; chain_ctr++) chain_start[chain_ctr] = -1;

  /* count the number of atoms */
  atom_ctr=0;
  while(fgets(line,LONGSTRING,fptr)!=NULL){
       // the line must start with either ATOM or HETATM
       if(strncmp(line,"ATOM", 4) &&  strncmp(line,"HETATM", 6)) continue;
       /* if it's a hydrogen - skip */
       if (line[PDB_ATOM_ATOM_NAME] == 'H' ||  line[PDB_ATOM_ATOM_NAME+1] == 'H')  continue;

 	   atom_ctr++;
  }

  map_length=atom_ctr;/**** --> map_length == number of atoms ****/
  if ( map_length < 1){
     fprintf(stderr,"\n ERROR : no atom found in %s file\n\n",filename);
     return 1;
  } else {
# ifdef VERBOSE
      printf ("\t %d atoms\n", map_length);
# endif
  }

  /* allocate the space for the pdb_map */
  pdb_map = (Pdb_map *) calloc((map_length+1),sizeof(Pdb_map));

  /* rewind the pdb file */
  rewind(fptr);
  
  /* now reread the pdb file and store all information */
  atom_ctr=0;
  old_chain_id = '*';
  chain_ctr = 0;
  memset (line, LONGSTRING, 0);
  while(fgets(line,LONGSTRING,fptr)!=NULL){
       // the line must start with either ATOM or HETATM
       if (strncmp(line,"ATOM", 4) &&  strncmp(line,"HETATM", 6)) continue;
       /* if it's a hydrogen - skip */
       if (line[PDB_ATOM_ATOM_NAME] == 'H' ||  line[PDB_ATOM_ATOM_NAME+1] == 'H')  continue;

       strncpy ( tmp, line+PDB_ATOM_ATOM_NO, PDB_ATOM_ATOM_NO_LEN);
	   tmp[PDB_ATOM_ATOM_NO_LEN] = '\0';
       pdb_map[atom_ctr].no_atom=atoi(tmp);

	   if ( ! isspace (line[PDB_ATOM_ALTLOC]) ) {
	        fprintf (stderr, "Warning: atom %d has alternative location in %s.\n",
		    atoi(tmp), filename);
	   }

       strncpy (tmp, line+PDB_ATOM_ATOM_NAME,  PDB_ATOM_ATOM_NAME_LEN);
	   tmp[PDB_ATOM_ATOM_NAME_LEN] = '\0';
	   j= 0;
	   for (i=0; i<PDB_ATOM_ATOM_NAME_LEN; i++) {
	       if (isalnum(tmp[i])) {
		   pdb_map[atom_ctr].name_atom[j] = tmp[i];
		   j++;
	       }
	   }

       strncpy ( tmp, line+PDB_ATOM_RES_NAME, PDB_ATOM_RES_NAME_LEN);
	   tmp[PDB_ATOM_RES_NAME_LEN] = '\0';
	   
       strcpy(pdb_map[atom_ctr].name_res,tmp);
       pdb_map[atom_ctr].chain_id = line[ PDB_ATOM_CHAINID];

	   if ( line[ PDB_ATOM_CHAINID] != old_chain_id ) {
# ifdef VERBOSE
	       printf ("new chain found : %c \n", line[ PDB_ATOM_CHAINID]);
# endif
	       old_chain_id =  line[ PDB_ATOM_CHAINID];
	       chain_start [chain_ctr] =  atom_ctr;
	       chain_ctr++;
	       if (chain_ctr>=MAX_CHAINS) {
	            fprintf(stderr, "Max number of chains exceeded. Increase MAX_CHAINS and recompile.\n");
	            exit(1);
	       }
	   }
	   
	   /* crappy thing about pdb: "residue number" may  
	      carry an insertion code, which is a letter */
	   memcpy( pdb_map[atom_ctr].no_res, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN );
	   pdb_map[atom_ctr].no_res[PDB_ATOM_RES_NO_LEN] = line[PDB_ATOM_INS_CODE];
	   /* some sanity checking: */
	   /* if the name of the residue changes, but the residue number does not, yell foul */
	   /* (it was seen in some pdb's ... ) */ 
	   if (atom_ctr && strncmp(pdb_map[atom_ctr-1].name_res, pdb_map[atom_ctr].name_res,PDB_ATOM_RES_NAME_LEN) &&
		    !strncmp(pdb_map[atom_ctr-1].no_res,  pdb_map[atom_ctr].no_res, PDB_ATOM_RES_NO_LEN+PDB_ATOM_INS_CODE) ) {
               fprintf(stderr,"\nInconsistency in  %s file: \n",filename);
               fprintf(stderr,"\tatoms %d and %d carry the same residue number (%s),\n",
               pdb_map[atom_ctr-1].no_atom,  pdb_map[atom_ctr].no_atom, pdb_map[atom_ctr-1].no_res);
               fprintf(stderr,"\tbut not the same type (%s and %s, respectively).\n",
                  pdb_map[atom_ctr-1].name_res, pdb_map[atom_ctr].name_res );
               return 1;
	   }
	   

	   strncpy ( tmp, line+PDB_ATOM_X , PDB_ATOM_X_LEN);
	   tmp[8] = '\0';
           pdb_map[atom_ctr].x=atof(tmp);
	   strncpy ( tmp, line+PDB_ATOM_Y , PDB_ATOM_Y_LEN );
	   tmp[8] = '\0';
           pdb_map[atom_ctr].y=atof(tmp);
	   strncpy ( tmp, line+PDB_ATOM_Z , PDB_ATOM_Z_LEN);
	   tmp[8] = '\0';
           pdb_map[atom_ctr].z=atof(tmp);   
   
	   atom_ctr++;
   }

   *map_length_ptr = map_length;
   *pdb_map_ptr = pdb_map;
   fclose (fptr);
   return 0;

} 


