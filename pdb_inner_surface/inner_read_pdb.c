#include "inner.h"

# define BUFFLEN 150

/************************************************************/
/************************************************************/
int find_calpha ( Protein *protein){

    protein->calpha = emalloc (protein->length*sizeof(Atom*));
    int resctr;
    for (resctr=0; resctr<protein->length; resctr++) {
        Residue * res = protein->residue + resctr;
        protein->calpha[resctr] = NULL;
        int i;
        for (i=0; i< res->no_atoms; i++) {
			if ( ! strcmp(res->atom[i].type, "CA")) {
				protein->calpha[resctr] = &(res->atom[i]);
				break;
			}
        }
        if (!protein->calpha[resctr]) {
        	fprintf (stderr, "C-alpha not found for residue %s%s\n",res->res_type, res->pdb_id);
        	exit(1);
        }
    }
    return 0;

}
/************************************************************/
int read_pdb(char * pdbname, Protein * protein, char chain_id) {

	Residue * sequence;
	FILE * fptr = NULL;
	char line[BUFFLEN];
	char oldresno[PDB_ATOM_RES_NO_LEN + 1];
	char oldrestype[PDB_ATOM_RES_NAME_LEN + 1];
	char tmp[PDB_ATOM_X_LEN + 1];
	int atomctr, resctr, no_res, skip;

	/* open file */
	fptr = fopen(pdbname, "r");
	if (!fptr) {
		fprintf(stdout, "Cno %s.\n", pdbname);
		return 1;
	}
	/* check for NMR and theory files */
	while (fgets(line, BUFFLEN, fptr) != NULL) {
		if (!strncmp(line, "EXPDTA", 6)) {
			if (strstr(line, "NMR")) {
				fprintf(stdout, "%s is an NMR file.\n", pdbname);
				return 1;
			} else if (strstr(line, "THEO")) {
				fprintf(stdout, "%s is a theoretical  file.\n", pdbname);
				return 1;
			}
		}
	}

	/* count residues */
	memset(line, 0, BUFFLEN);
	memset(oldresno, 0, PDB_ATOM_RES_NO_LEN + 1);
	resctr = 0;
	rewind(fptr);
	while (fgets(line, BUFFLEN, fptr) != NULL) {
		if (!strncmp(line, "ATOM", 4) || !strncmp(line, "HETATM", 6)) {
			/* if chain specified, check if it is the right one: */
			if (chain_id && chain_id != line[PDB_ATOM_CHAINID])
				continue;

			if (strncmp(line + PDB_ATOM_RES_NO, oldresno,
					PDB_ATOM_RES_NO_LEN)) {

				strncpy(oldresno, line + PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN);
				oldresno[PDB_ATOM_RES_NO_LEN] = '\0';
				/* printf ( "New residue number:  %s \n", oldresno); */
				resctr++;
			}

		}
	}

	no_res = resctr;
	protein->length = no_res;
	printf("reading chain %c   no residues: %d\n", chain_id, no_res);

	/* allocate space */
	sequence = NULL;
	sequence = calloc(no_res, sizeof(Residue));
	if (!sequence) {
		fprintf( stderr, "Error allocating sequence space.\n");
		exit(1);
	}
    protein->residue= sequence;

	/* read in the atoms */
	rewind(fptr);
	memset(line, 0, BUFFLEN);
	memset(oldresno, 0, PDB_ATOM_RES_NO_LEN + 1);
	memset(oldrestype, 0, PDB_ATOM_RES_NAME_LEN + 1);
	resctr = -1;
	atomctr = 0;
	skip = 0;
	while (fgets(line, BUFFLEN, fptr) != NULL) {
		if (skip) { /*skipping the alt locations; see below */
			skip = 0;
			continue;
		}
		if (!strncmp(line, "ATOM", 4) || !strncmp(line, "HETATM", 6)) {

			/* if chain specified, check if it is the right one: */
			if (chain_id && chain_id != line[PDB_ATOM_CHAINID])
				continue;
			/* skip hydrohens */
			if (line[PDB_ATOM_ATOM_NAME] == 'H'
					|| line[PDB_ATOM_ATOM_NAME + 1] == 'H') {
				continue;
			}
			/* adjust the counters */
			if (strncmp(line + PDB_ATOM_RES_NO, oldresno,
					PDB_ATOM_RES_NO_LEN)) {
				strncpy(oldresno, line + PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN);
				strncpy(oldrestype, line + PDB_ATOM_RES_NAME,
						PDB_ATOM_RES_NAME_LEN);
				oldresno[PDB_ATOM_RES_NO_LEN] = '\0';
				oldrestype[PDB_ATOM_RES_NAME_LEN] = '\0';
				resctr++;
				atomctr = 0;

				sequence[resctr].no_atoms = 1;
				strncpy(sequence[resctr].pdb_id, oldresno, PDB_ATOM_RES_NO_LEN);
				sequence[resctr].pdb_id[PDB_ATOM_RES_NO_LEN] = '\0';
				strncpy(sequence[resctr].res_type, oldrestype,
						PDB_ATOM_RES_NAME_LEN);
				sequence[resctr].res_type[PDB_ATOM_RES_NAME_LEN] = '\0';

			} else {
				atomctr++;
				sequence[resctr].no_atoms = atomctr + 1;
				if (atomctr >= MAX_NO_ATOMS) {
					fprintf( stdout,
							"Error: I thought every aa has < %d atoms.\n",
							MAX_NO_ATOMS); /*stdout so perl can capture it easily */
					exit(1);
				}
			}
			sequence[resctr].chain =  line[PDB_ATOM_CHAINID];
			/* read in atom info */
			strncpy(sequence[resctr].atom[atomctr].type,
					line + PDB_ATOM_ATOM_NAME,
					PDB_ATOM_ATOM_NAME_LEN);
			sequence[resctr].atom[atomctr].type[ PDB_ATOM_ATOM_NAME_LEN] = '\0';
			strncpy(tmp, line + PDB_ATOM_X, PDB_ATOM_X_LEN);
			tmp[PDB_ATOM_X_LEN] = '\0';
			sequence[resctr].atom[atomctr].coord[0] = atof(tmp);
			strncpy(tmp, line + PDB_ATOM_Y, PDB_ATOM_Y_LEN);
			tmp[PDB_ATOM_Y_LEN] = '\0';
			sequence[resctr].atom[atomctr].coord[1] = atof(tmp);
			strncpy(tmp, line + PDB_ATOM_Z, PDB_ATOM_Z_LEN);
			tmp[PDB_ATOM_Z_LEN] = '\0';
			sequence[resctr].atom[atomctr].coord[2] = atof(tmp);

			sequence[resctr].atom[atomctr].parent_residue = sequence + resctr;

			if (!isspace(line[PDB_ATOM_ALTLOC])) {
				/* atom   has alternative location  */
				/* skip the next line */
				skip = 1;
			}

		}

	}
	/* close file */
	fclose(fptr);
	/* clean PDB id tags from spaces */
	for (resctr = 0; resctr < no_res; resctr++) {
		Residue * res = sequence + resctr;
		if (string_clean(res->pdb_id, PDB_ATOM_RES_NO_LEN + 1)) {
			fprintf(stderr,
					"Error in read_pdb(): empty id string for residue with sequential no %d.\n",
					resctr);
			return 1;
		}
		if (string_clean(res->res_type, PDB_ATOM_RES_NAME_LEN)) {
			fprintf(stderr,
					"Error in read_pdb(): empty id string for residue with sequential no %d.\n",
					resctr);
			return 1;
		}

		res->res_type_short = single_letter(res->res_type);
        int i;
        for (i=0; i< res->no_atoms; i++) {
        	string_clean(res->atom[i].type, PDB_ATOM_ATOM_NAME_LEN +1);
        }
	}


	find_calpha(protein);

	return 0;
}


