
#include "ifc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

int moments_of_inertia ( Pdb_map * pdb_map, Chain_info * chain_info, int chain, int interface[], double moments[]) {

    int c = chain;
    int atom, residue, i, j;
    double x[3] = {0.0};
    double dist_sq;
    double I[3][3] = {{0.0}};
    
    
    for ( residue=0; residue<chain_info[c].number_of_residues-1; residue++) {
	
	if ( ! interface[residue] )  continue;
	for (atom =chain_info[c].residue_start[residue];
	     atom <=chain_info[c].residue_end[residue ]; atom ++) {

		 x[0] = pdb_map[atom ].x;
		 x[1] = pdb_map[atom ].y;
		 x[2] = pdb_map[atom ].z;
		 
		 dist_sq = 0.0;
		 for ( i=0; i<3; i++) {
		     dist_sq =  x[i]* x[i];
		 }
		 for ( i=0; i<3; i++) {
		     for ( j=0; j<3; j++) {
			 if ( i==j ) {
			     I[i][j] += dist_sq;
			 }
			 I[i][j] -= x[i]*x[j]; 
		     }
		 }
	}
    }

    {
# define SIZE 3
# define  WORK_SIZE 1000
	double VL[SIZE][SIZE], VR[SIZE][SIZE];
	double WR[SIZE], WI[SIZE];
	double WORK [WORK_SIZE];
	int LWORK;
	int INFO, N, LDA, LDVL, LDVR;
	char JOBVL, JOBVR;
	
	JOBVL = JOBVR = 'V';
	N = LDA = LDVL = LDVR = SIZE;
	LWORK = WORK_SIZE;
	dgeev_ ( &JOBVL, &JOBVR, &N, I, &LDA, WR, WI, VL, &LDVL, VR, &LDVR,
		 WORK, &LWORK, &INFO );
	for (i=0; i < SIZE; i++){
	    moments[i] = WR[i];
	}
    }

    return 0;
}
