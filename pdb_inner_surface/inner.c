/* readin in 2 pdb, presumably roughly aligned structurally, */
/* and a list of residues in the first structure on which the*/
/* improved alignemnt should focus */
/* output the tranformation matrix for the second structure */
/* for the improved alignment */
/* Ivana, Oct, 2010 */
#include "inner.h"

# define MAX_DIST_TO_CONSIDER 4.0

/* the following functions are defined below main()*/


/*************************************************************************/
/*************************************************************************/
int main(int argc, char * argv[]) {

	Protein protein;
	double principal_moment[3], principal_axis[3][3];
	char pdbname[BUFFLEN] = "\0";

	if (argc < 2) {
		printf("Usage: %s <pdbname> [min/max]\n", argv[0]);
		exit(1);
	}
	sprintf(pdbname, "%s", argv[1]);
	int use_max =1; // for now using the axis with max moment of inertia is the default
	if (argc>2 && !strcmp("min", argv[2])) use_max = 0;

	/*********************************************/
	/* read in the structure                     */
	if (read_pdb(pdbname, &protein, '\0'))
		exit(1);

	if (1) {
		printf("read in %20s, number of residues %d \n", pdbname, protein.length);
		/* sanity: */
		if (!protein.residue[0].no_atoms) {
			fprintf(stderr, "no atoms in %s(?)\n", pdbname);
			exit(1);
		}
	}
	
	/*********************************************/
	/* rotate so the desired  principal axis is z */
	// no:t the function will move center of mass to origin
	find_principal_axes (&protein, principal_moment,  principal_axis);
	int moment_index = use_max ? max_index(principal_moment) :  min_index(principal_moment);
	rotate_structure (&protein, principal_axis[moment_index]);

	/* coordinates to cylindrical */
	coords2cylindrical(&protein);
	/* bin atoms by z and theta */
	double del_z = 1.0;
	int number_of_theta_bins = 36, number_of_z_bins;
	Atom *** bin; // we are binning atom pointers
	bin_atoms(&protein, del_z, number_of_theta_bins,  &bin, &number_of_z_bins);
	/* find min rho atom for each bin */
	/* find set of residues that belong to min rho */
	/* half the bin size and repeat */
	/* if the set of residues remained the same, output the set of residues */

	return 0;
}
