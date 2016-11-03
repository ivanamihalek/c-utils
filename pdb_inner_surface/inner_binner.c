#include "inner.h"

#define MAX_ATOM_BINS  100
#define SMALL_BUF 50
int bin_atoms(Protein *protein, double z_step, int  number_of_theta_bins, Atom **** bin_ptr, int *number_of_z_bins_ptr){

	if (z_step<1.e-4) {
		fprintf (stderr, "z_tep too small. Convergence problem?\n");
		exit(1);
	}

	Atom ***bin;
	int *bin_size;
	double min_z =  10000;
	double max_z = -10000;
    int i;

	int resctr;
    for (resctr=0; resctr<protein->length; resctr++) {
        Residue * res = protein->residue + resctr;
        for (i=0; i< res->no_atoms; i++) {
        	double z = res->atom->z;
			if (z < min_z) min_z = z;
			if (z > max_z) max_z = z;
        }
    }
    min_z -= z_step/2;
    max_z += z_step/2;
    char smallbuf[SMALL_BUF] = {'\0'};
    sprintf(smallbuf, "%5.0lf", (max_z-min_z)/z_step);
    int number_of_z_bins = atoi(smallbuf);
    int tot_bins = number_of_z_bins*number_of_theta_bins;
    bin = emalloc(tot_bins*sizeof(Atom**));
    bin_size = emalloc(tot_bins*sizeof(int));
    for (i=0; i<tot_bins; i++) {
    	bin[i] = emalloc(MAX_ATOM_BINS*sizeof(Atom*));
    }

    /* find my bin */
    for (resctr=0; resctr<protein->length; resctr++) {
        Residue * res = protein->residue + resctr;
        for (i=0; i< res->no_atoms; i++) {
           	double z = res->atom->z;
           	double theta = res->atom->theta;
        	int z_index = (z - min_z)/(max_z-min_z)*number_of_z_bins;
        	int theta_index = theta/(2*M_PI)*number_of_theta_bins;
        	int onedim_index = z_index*number_of_theta_bins + theta_index;
        	if (onedim_index >= tot_bins) {
        		fprintf (stderr, "Underdim bin array:  %d  %d \n", onedim_index, tot_bins);
        		exit(1);

        	}
         	int atom_ctr = bin_size[onedim_index];
        	bin_size[onedim_index] ++;
        	if (atom_ctr >= MAX_ATOM_BINS) {
        		fprintf (stderr, "Underdim MAX_ATOM_BINS\n");
        		exit(1);
        	}
        	bin[onedim_index][atom_ctr] = res->atom;
        }
    }

    for (i=0; i<tot_bins; i++) {
    	int a;
    	double min_rho = 10000;
    	Atom * min_atom_ptr = NULL;
    	for (a=0; a<bin_size[i]; a++) {
    		double rho = bin[i][a]->rho;
    		if (rho < min_rho) {
    			min_rho = rho;
    			min_atom_ptr = bin[i][a];
    		}
    	}
    	if (!min_atom_ptr) continue;
    	Residue *res = min_atom_ptr->parent_residue;
    	printf ( " %c %s \n", res->chain, res->pdb_id);

    }

	*bin_ptr = bin;
	*number_of_z_bins_ptr = number_of_z_bins;
	return 0;
}
