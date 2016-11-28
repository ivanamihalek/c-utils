/* residues lining the inside of a molecule of more or less cylindrical geometry */

#include "inner.h"

# define MAX_DIST_TO_CONSIDER 4.0



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

	if (0) { // sounds like a good idea, but its not working
		double del_z = 1.0;
		int number_of_theta_bins = 72, number_of_z_bins;
		Atom *** bin; // we are binning atom pointers
		int *bin_size; // number of atoms in each bin
		bin_atoms(&protein, del_z, number_of_theta_bins,  &bin,  &bin_size, &number_of_z_bins);
		/* find min rho atom for each bin */
		/* find set of residues that belong to min rho */
		/* half the bin size and repeat */
		/* if the set of residues remained the same, output the set of residues */
		int i, tot_bins = number_of_theta_bins * number_of_z_bins;
		for (i = 0; i < tot_bins; i++) {
			int a;
			double min_rho = 10000;
			Atom * min_atom_ptr = NULL;
			for (a = 0; a < bin_size[i]; a++) {
				double rho = bin[i][a]->rho;
				if (rho < min_rho) {
					min_rho = rho;
					min_atom_ptr = bin[i][a];
				}
			}
			if (!min_atom_ptr) continue;
			Residue *res = min_atom_ptr->parent_residue;
			// check direction of the sidechain - keep only if facing inwards
			// if not pointing inward residue continue
		    if (is_pointing_inward(res) )
		    	printf(" %c %s %5.1lf \n", res->chain, res->pdb_id, min_rho);
		}
	}

	// other approach - distribution in rho, take anything that is on standard deviation or more
	// closer to the center than the average
	double avg_rho, stdev_rho;
	distribution_of_rho(&protein, &avg_rho, &stdev_rho);
	printf("  avg rho %5.1lf  stdev %5.1lf \n",  avg_rho, stdev_rho);
	int resctr;
    for (resctr=0; resctr<protein.length; resctr++) {
        Residue * res = protein.residue + resctr;
        int inner = 0, i;
        double min_rho = 0;
        for (i=0; i< res->no_atoms; i++) {
			//if (  (res->atom[i].rho - avg_rho)/stdev_rho > 1) {
			if (  (res->atom[i].rho - avg_rho)/stdev_rho < -1) {
				min_rho = res->atom[i].rho;
				inner=1;
				break;
			}
        }
	    if (inner )
	    	printf(" %c %s %5.1lf \n", res->chain, res->pdb_id, min_rho);

     }

	return 0;
}
