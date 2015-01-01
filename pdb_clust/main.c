# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "pdbclust.h"

int main ( int argc, char * argv[]) {
    
    Residue * sequence;
    char pdbname[150]        = {'\0'};
    char selected_res_file[150]   = {'\0'};
    Boolean cluster_counting = FALSE, cluster_scoring = FALSE;
    Boolean ok_so_far;
    double cutoff_dist       = 5.0;
    int no_res, res_ctr;
    int a, c, res1, res2; // assorted counters
    int * selected           = NULL;
    double **distmat         = NULL;
    int ** adj_matrix        = NULL;
    int ** cluster           = NULL;
    char chain_id            = '\0';
    FILE * fclust            = NULL; // currently, we just set this to stdout somewhere below

    if ( argc < 5 ) {
	fprintf (stderr, "Usage: %s <pdbfile> <chain> <cutoff_dist> <selected residues (file)>  [-c] [-s] \n", argv[0]);
	fprintf (stderr, "To ignore the chain label, provide it as \"-\".\n");
	exit (1);
    }
    
    sprintf ( pdbname, "%s", argv[1]);
    chain_id =  argv[2][0]=='-' ? '\0' :  argv[2][0];
    cutoff_dist = atof(argv[3]);
    sprintf (selected_res_file, "%s", argv[4]);
    for (a=5; a < argc; a++) {
	if (argv[a][0] != '-') continue;
	switch (argv[a][1]) {
	case 'c':
	    cluster_counting = TRUE;
	    break;
	case 's':
	    cluster_scoring = TRUE;
	}
    }
    
    /* input the structure */ 
    if ( read_pdb (pdbname, &chain_id, &sequence, &no_res) ) exit (1);
    printf ("there are %d residues in %s, chain %c.\n",  no_res, pdbname, chain_id);
    
    /* input residue selection */
    selected = (int*) emalloc(no_res*sizeof(int));
    ok_so_far = (selected != NULL);
    if (ok_so_far) {
	ok_so_far = read_residue_selection (selected_res_file, sequence, no_res, selected)? FALSE:TRUE ;
    }

    if (ok_so_far) {
	/* allocate space for dist matrix */
	distmat = (double**) dmatrix (0,  no_res-1, 0,  no_res-1);
	ok_so_far = (distmat != NULL);
    }

    if (ok_so_far) {
        /* calculate dist */
	ok_so_far = determine_dist_matrix(distmat,  sequence, no_res) ? FALSE:TRUE;
    }

    /* adjacency matrix, needed for eithre cluster counting or cluster scoring */
    if (ok_so_far && (cluster_counting || cluster_scoring) ) {
	adj_matrix  =  imatrix (0, no_res-1, 0, no_res-1);
	ok_so_far   =  (adj_matrix != NULL);
	if ( ok_so_far) {
	    for (res1=0;  res1 < no_res; res1++ ) {
		adj_matrix [res1][res1] = 1;
		for (res2= res1+1; res2 < no_res; res2++ ) {
		    adj_matrix[res1][res2] = ( distmat[res1][res2] < cutoff_dist);
		    adj_matrix[res2][res1] = adj_matrix[res1][res2];
		}
	    }

	}
    }
    
    /* cluster counting ... */
    if (ok_so_far && cluster_counting) {
	int  no_of_clusters, max_size, secnd_max_size;
	int * cluster_count;
	    
	cluster_count       =  (int *) emalloc ( (no_res+1)*sizeof(int));
	cluster             =  imatrix (0, no_res, 0, no_res);
	
	cluster_counter (no_res, adj_matrix, selected, cluster_count, & no_of_clusters,
			 &max_size, &secnd_max_size, cluster);

	/* output */
	fclust = stdout;
	for ( c=0; c <= no_res; c++) {
	    if ( ! cluster[c][0] ) {
		continue;
	    }
	    if ( !c ) {
		fprintf ( fclust,"\t isolated:\n");
	    } else {
		fprintf ( fclust,"\t cluster size: %3d \n", cluster[c][0]); 
	    }
	    for ( res_ctr=1; res_ctr <=  cluster[c][0]; res_ctr++) {
		fprintf ( fclust, "%s  \n", sequence[ cluster[c][res_ctr] ].pdb_id );
	    }
	
	}
    }

    /* cluster scoring ... */
    if (ok_so_far && cluster_scoring) {
	double score, avg, std_dev, z;
	int number_of_selected_residues = 0;
	cluster_score (no_res, selected, adj_matrix, &score);
	/* find avg and stddev in the set of random picks */
	for (res1=0; res1< no_res; res1++) number_of_selected_residues += selected[res1];
	std_dev_over_S (no_res, number_of_selected_residues, adj_matrix, &avg, &std_dev, TRUE);
	
	/* evaluate and store the z-score */
	z = (std_dev>1.e-5) ? (score - avg)/std_dev : 0.0;

	fprintf ( stdout,"\n z_score (for the overall clustering pattern: %8.2f \n", z); 

    }

    if (sequence) free (sequence);
    if (selected) free (selected);
    if (distmat)  free_dmatrix(distmat, 0,  no_res-1, 0,  no_res-1);
    
    return 0;

    
    
}
