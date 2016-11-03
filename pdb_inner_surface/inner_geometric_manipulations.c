#include "inner.h"

/*************************************************************************/
double diagonalize(double I[3][3], double evec[3][3], double evals[3]) {

	char jobz = 'V'; /* find evalues and evectors */
	char uplo = 'L'; /* amtrix is stored as lower (fortran convention) */
	int N = 3; /* the order of matrix */
	int leading_dim = N;
	int i, j, retval;
	double A[N * N];
	double workspace[3 * N];
	int workspace_size = 3 * N;
	void dsyev_(char * jobz, char * uplo, int* N, double * A, int * leading_dim,
			double * eigenvalues, double *workspace, int *workspace_size,
			int * retval);
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			A[i * 3 + j] = I[i][j];
		}
	}
	dsyev_(&jobz, &uplo, &N, A, &leading_dim, evals, workspace, &workspace_size,
			&retval);
	if (retval) {
		fprintf(stderr, "error in dsyev()\n");
		exit(1);
	}

	for (i = 0; i < 3; i++) {
		double norm = 0;
		for (j = 0; j < 3; j++) {
			evec[i][j] = A[i * 3 + j];
			norm += evec[i][j] * evec[i][j];
		}
	}
	return 0;
}

/*************************************************************************/
int find_principal_axes (Protein *protein, double principal_moment[3], double principal_axis[3][3]) {
	/*********************************************/
	/* calculate moments of inertia              */
	Atom **atom = protein->calpha;
	int no_atoms = protein->length;
	double cm[3] =  { 0.0 };
	double I[3][3] = { { 0.0 } };
	int a, i, j;

    /*cm*/
	for (a = 0; a < no_atoms; a++) {
		for (i = 0; i < 3; i++) {
			cm[i] += atom[a]->coord[i];
		}

	}
	for (i = 0; i < 3; i++) cm[i] /= no_atoms;

	/* points in cm */
	for (a = 0; a < no_atoms; a++) {
		for (i = 0; i < 3; i++) {
			atom[a]->coord[i] -= cm[i];
		}
	}

	/* moment of inertia tensor */
	memset(I[0], 0.0, 3 * 3 * sizeof(double));
	for (a = 0; a < no_atoms; a++) {
		for (i = 0; i < 3; i++) {
			double d1 = atom[a]->coord[(i + 1) % 3];
			double d2 = atom[a]->coord[(i + 2) % 3];
			I[i][i] += d1 * d1 + d2 * d2;
			for (j = i+1; j < 3; j++) {
				I[i][j] -= atom[a]->coord[i] * atom[a]->coord[j];
			}
		}
	}
	for (i = 0; i < 3; i++) {
		for (j = i+1; j < 3; j++) {
			I[j][i] = I[i][j];
		}
	}


	diagonalize(I, principal_axis, principal_moment);
	if (0) for (i=0; i<3; i++) {
		printf ("%5.2lf  %5.2lf   ", cm[i], principal_moment[i]);
		for (j=0; j<3; j++) printf ("%8.5lf  ", principal_axis[i][j]);
		printf ("\n");
	}

	return 0;
}

/************************************************************/
int quat_to_R(double quat[4], double **R) {

	double q0 = quat[0];
	double q1 = quat[1];
	double q2 = quat[2];
	double q3 = quat[3];

	R[0][0] = 1 - 2 * q2 * q2 - 2 * q3 * q3;
	R[1][1] = 1 - 2 * q3 * q3 - 2 * q1 * q1;
	R[2][2] = 1 - 2 * q1 * q1 - 2 * q2 * q2;

	R[0][1] = 2 * q1 * q2 - 2 * q0 * q3;
	R[1][0] = 2 * q1 * q2 + 2 * q0 * q3;

	R[0][2] = 2 * q1 * q3 + 2 * q0 * q2;
	R[2][0] = 2 * q1 * q3 - 2 * q0 * q2;

	R[1][2] = 2 * q2 * q3 - 2 * q0 * q1;
	R[2][1] = 2 * q2 * q3 + 2 * q0 * q1;

	return 0;

}

/*************************************************************************/
int rotate_structure(Protein *protein, double new_z_direction[3]) {

	int i;

	int already_z = 1;
	for (i=0; i<2; i++) {
		if (fabs(new_z_direction[i])>1.e-6) already_z=0;
	}
	if (fabs(new_z_direction[2]-1)>1.e-6) already_z=0;

	if (!already_z) {
		// my current case is already oriented
		fprintf (stderr, "Rotation not implemented -please do so.\n");
		exit(1);
	}

# if 0
	double **R;
	double ** rot_coord;
	double q[4], T[3], rmsd;
	if (!(R = dmatrix(3, 3))) return 1; /* compiler is bugging me otherwise */
	quat_to_R(q, R);

	/*******************************************/
	/* output the tfm matrix                   */
	if (1) {
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				printf(" %8.4lf", R[i][j]);
			}
		}
	}

	/* rotate (the first structure), move back to the old cm,
	 and output */
	for (a = 0; a < no_atoms[0]; a++) {

		for (i = 0; i < 3; i++) {
			rot_coord[a][i] = 0.0;
			for (j = 0; j < 3; j++) {
				rot_coord[a][i] += R[i][j] * atom[0][a].coord[j];
			}
		}

	}
	free_dmatrix(R);
#endif
	return 0;
}

/*************************************************************************/
int coords2cylindrical(Protein *protein) {

	int resctr;
    for (resctr=0; resctr<protein->length; resctr++) {
        Residue * res = protein->residue + resctr;
        int i;
        for (i=0; i< res->no_atoms; i++) {
            double x = res->atom[i].coord[0];
            double y = res->atom[i].coord[1], z = res->atom[i].coord[2];
			res->atom[i].rho = sqrt(x*x+y*y);
			if (res->atom[i].rho >1.e-6) {
				res->atom[i].theta = acos(x/res->atom[i].rho);
				if (y<0) res->atom[i].theta =  2*M_PI - res->atom[i].theta;

			} else	 {
				res->atom[i].theta = 0; // it shouldn't matter in this case
			}
			res->atom[i].z = z;
        }
    }
	return 0;
}

/*************************************************************************/
int is_pointing_inward(Residue * res) {

	//finc ca and cb
	int i;
	double * calpha = NULL;
	double *  cbeta = NULL;

    for (i=0; i< res->no_atoms; i++) {
		if ( ! strcmp(res->atom[i].type, "CA")) {
			calpha = res->atom[i].coord;
		} else if  (! strcmp(res->atom[i].type, "CB")) {
			cbeta = res->atom[i].coord;
		}
    }
    if ( !calpha || !cbeta) return 0;


    double diff[3], dotprod;
    // we are interested in x and y only - the whole story is taking place in cylindrical geometry
    for (i=0; i<2; i++) diff[i]= cbeta[i] - calpha[i];
    dotprod = 0.0; for (i=0; i<2; i++) dotprod += diff[i]* (- calpha[i]);
    if (dotprod > 0)  return 1;
	return 0; // here this means no, and not SUCCESS
}

int distribution_of_rho (Protein* protein, double *avg_ptr, double *stdev_ptr) {
	double avg = 0;
	double avg_sq = 0;
    int resctr, ct= 0;
    for (resctr=0; resctr<protein->length; resctr++) {
        Residue * res = protein->residue + resctr;
        protein->calpha[resctr] = NULL;
        int i;
        for (i=0; i< res->no_atoms; i++) {
        	double rho = res->atom[i].rho;
        	avg += rho;
        	avg_sq += rho*rho;
        	ct += 1;
        }
     }
     if (!ct) {
    	 fprintf (stderr, "no atomes ?\n");
    	 exit(1);
     }

     avg /= ct;
     avg_sq /= ct;
     if (avg_sq < avg*avg) {
    	 fprintf (stderr, "huh?\n");
    	 exit(1);
     }
     *avg_ptr = avg;
     *stdev_ptr = sqrt(avg_sq - avg*avg);

    return 0;
}

