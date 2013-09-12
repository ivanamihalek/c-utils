# include "postp.h"

int main ( int argc, char * argv[]) {

    Options options;
    Protein protein;
    Alignment big_alignment;
    Alignment * alignment = NULL;
    int retval;
    int almtctr1, almtctr2;
    double **score = NULL, corr,  pctg_gaps;
    double **clustering_score = NULL;
    double *area, *distance;
    int **rank_order= NULL,**res_rank=NULL,**int_cvg=NULL ;
    int ** correlated = NULL, **almt2prot = NULL, **prot2almt = NULL;
    /* command file is required */
    if ( argc < 2 ) {
	fprintf ( stderr, "Usage: %s <command file>.\n", argv[0]);
	exit (0);
    }
    retval = read_cmd_file ( argv[1], &options);
    if (retval) exit(retval);
    retval = logger (&options, INTRO, "");
    if (retval) exit(retval);
   
    
    /*******************************************/
    /*                                         */
    /*  PDB input                              */
    /*                                         */
    /*******************************************/
    if ( ! options.pdbname[0]) {
	fprintf (stderr, "%s cannot work without structure (cmd file was %s).\n",
		 argv[0], argv[1]);
	exit (1);
	
    } else {

	/* warn if no chain given */
	if ( !options.chain) {
	    retval = logger (&options, WARN, "No chain specified. Using the first one.");
	    if ( retval) exit (1);
	}
	if (retval) exit(retval);
	/* read in the structure */
	retval = read_pdb (options.pdbname, &protein, options.chain);
	if (retval) exit(retval);

   }

    
    /*******************************************/
    /*                                         */
    /*  alignment scoring                      */
    /*                                         */
    /*******************************************/
    if ( ! ( alignment = emalloc ( options.no_of_alignments*sizeof(Alignment)) )) return 1;
    if ( ! ( score = emalloc ( options.no_of_alignments*sizeof(double*)) )) return 1;
    if ( ! ( rank_order = emalloc ( options.no_of_alignments*sizeof(int*)) )) return 1;
    if ( ! ( clustering_score = emalloc ( options.no_of_alignments*sizeof(double*)) )) return 1;
    if ( ! ( res_rank = emalloc ( options.no_of_alignments*sizeof(int*)) )) return 1;
    if ( ! ( int_cvg = emalloc ( options.no_of_alignments*sizeof(int*)) )) return 1;
    if ( ! ( almt2prot = emalloc ( options.no_of_alignments*sizeof(int*)) )) return 1;
    if ( ! ( prot2almt = emalloc ( options.no_of_alignments*sizeof(int*)) )) return 1;
    if ( ! ( area = emalloc ( options.no_of_alignments*sizeof(double)) )) return 1;
    if ( ! ( distance = emalloc ( options.no_of_alignments*sizeof(double)) )) return 1;

    printf ( "\t%8s   %20s  %8s  %8s  %8s  \n", "almt#", "name        ",  "<dist to qry>", "%gaps", "area");
    
    /* read in the overall alignment */
    retval = read_clustalw (options.almtname, &big_alignment);
    if (retval) exit(retval);
    
    for ( almtctr1 = 0; almtctr1 < options.no_of_alignments; almtctr1++) {

	/* read in the alignment */
	retval = reconstruct_alignment (options.namefile[almtctr1], &big_alignment, alignment + almtctr1);
	if (retval) exit(retval);

	/* pairwise distances btw the seqs */
	retval   = seq_pw_dist (alignment+almtctr1);
	if ( retval) return retval;
	/* average dist to the query in this alignment: */ 
	distance[almtctr1] = avg_dist_to_special (&options, alignment + almtctr1);
	/* percentage of gaps in the alignment: */
	pctg_gaps = (double) alignment->total_gaps/ ( (alignment+almtctr1)->length*(alignment+almtctr1)->number_of_seqs);
	/* make the residue scoring array */
	score[almtctr1] = emalloc ( alignment[almtctr1].length*sizeof(double));
	/* fill in the score array */ 
	scoring (&options,  alignment+almtctr1, score[almtctr1]);
	
	/* translate the scoring into rank order */
	rank_order[almtctr1] = emalloc ( alignment[almtctr1].length*sizeof(int));
	score2rank (score[almtctr1], rank_order[almtctr1], alignment[almtctr1].length);
	
	/* mapping between the protein and the alignment almtctr1 */
	if ( ! (almt2prot[almtctr1] = (int *) emalloc (alignment[almtctr1].length*sizeof(int))) )exit (1);
	if ( ! (prot2almt[almtctr1] = (int *) emalloc (protein.length*sizeof(int))) )exit (1);
	retval    = struct_almt_mapping (&protein, alignment+almtctr1, options.query,  prot2almt[almtctr1], almt2prot[almtctr1]);
	if (retval) exit(retval);
	
	/* find coverage info implied by the scoring array */
	if ( ! (res_rank[almtctr1] = (int*) emalloc (protein.length*sizeof(int))) ) exit (1);
	if ( ! (int_cvg[almtctr1] =  (int*) emalloc (protein.length*sizeof(int))) ) exit (1);
	coverage ( &protein, almt2prot[almtctr1], score[almtctr1], alignment[almtctr1].length,
		   res_rank[almtctr1], int_cvg[almtctr1] );
	/*clustering score*/
	clustering_score[almtctr1]  =  (double*) emalloc (protein.length*sizeof(double));
	if (!clustering_score[almtctr1]) exit(retval);
	clustering ( &protein,  res_rank[almtctr1], int_cvg[almtctr1], clustering_score[almtctr1]);
	/* cumulative clustering score*/
	area[almtctr1]  = area_over_coverage (int_cvg[almtctr1], clustering_score[almtctr1], protein.length);
					     
	printf ( "\t   %4d   %20s   %8.3lf     %8.3lf  %8.3lf \n",
		 almtctr1, options.namefile[almtctr1], distance[almtctr1], pctg_gaps, area[almtctr1]);
    }
    

    /* find the table of correlations */
    if ( ! (correlated = intmatrix ( options.no_of_alignments, options.no_of_alignments) ) ) return 1;
    for ( almtctr1 = 0; almtctr1 < options.no_of_alignments -1; almtctr1++) {
	correlated[almtctr1][almtctr1] = 1;
	for ( almtctr2 = almtctr1+1; almtctr2 < options.no_of_alignments; almtctr2++) {
	    if ( alignment[almtctr1].length != alignment[almtctr2].length  ) {
		fprintf ( stderr, "Error alignments in the files %s and %s ",
			  options.namefile[almtctr1], options.namefile[almtctr2]);
		fprintf ( stderr, "seem to be of unequal length: %d and %d.\n",
			  alignment[almtctr1].length ,  alignment[almtctr2].length);
		return 1;
	    }
	    corr = spearman ( rank_order[almtctr1], rank_order[almtctr2], alignment[almtctr1].length );
	    printf ( " %3d  %3d  %8.4lf\n", almtctr1, almtctr2, corr);
	    correlated[almtctr1][almtctr2] = ( corr > 0.9 );
	}
    }

    
    /* find corelated clusters (of sequence selections)*/
    {
	int  *cluster_count_per_size;
	int  no_of_clusters;
	int  max_size, secnd_max_size , ** cluster;
	int size = options.no_of_alignments;
	int i,j;
	double dist, ar, max_area, dist_at_max_area;
	double min_dist_at_max_area, min_dist, max_area_at_min_dist;
	int almt_no, min_dist_almt;
	int cluster_counter (int  no_of_things,  int *neighbors[],
			      int cluster_count_per_size[], int * no_of_clusters,
			      int * max_size, int * secnd_max_size , int * cluster[]);
	
	
	if ( ! ( cluster_count_per_size = emalloc (size*sizeof(int)))) return 1; 
	if ( ! (cluster = intmatrix ( size+1, size+1) ) ) return 1;
	retval = cluster_counter (size,  correlated,  cluster_count_per_size,  &no_of_clusters,
			 & max_size,  &secnd_max_size , cluster);
	if ( retval ) return 1;

	printf ( "number of clusters: %d \n", no_of_clusters);
	for (i=0; i<=size; i++ ) {
	    if ( ! cluster[i][0] ) continue;
	    if ( !i ) {
		printf ( "\t isolated:\n");
	    } else {
		printf ("\t cluster size: %3d \n", cluster[i][0]); 
	    }
	    for (j=1; j <= cluster[i][0]; j++ ) {
		printf ( "%3d ", cluster[i][j] );
	    }
	    printf ( "\n");
	}

	
	/* which cluster is the closest to the singled out sequence ("special") */
	min_dist_at_max_area = dist_at_max_area = 10;
	max_area_at_min_dist = min_dist = -10;
	min_dist_almt = -1;
	for (i=0; i<=size; i++ ) {
	    if ( ! cluster[i][0] ) continue;
	    
	    max_area = -100;
	    almt_no =  dist_at_max_area = -1;
	    
	    for (j=1; j <= cluster[i][0]; j++ ) {
		dist = distance[cluster[i][j]] ;
		ar =  area[cluster[i][j]] ;
		if ( max_area < ar ) {
		    max_area = ar;
		    dist_at_max_area = dist;
		    almt_no = cluster[i][j];
		}
	    }
	    if ( almt_no < 0 ) {
		fprintf ( stderr, "Error selecting the alignment (1)\n");
		exit (1);
	    }
	    
	    if ( min_dist_at_max_area > dist_at_max_area ) {
		min_dist = dist_at_max_area;
		max_area_at_min_dist = max_area;
		min_dist_almt = almt_no;
	    }
	}
	if ( min_dist_almt < 0 ) {
	    fprintf ( stderr, "Error selecting the alignment (2)\n");
	    exit (1);
	}
	
	printf ( "choosing alignment %d %s (distance: %5.3f  area: %6.3f)\n",
		min_dist_almt, options.namefile[min_dist_almt],  min_dist, max_area_at_min_dist);
	
	
	free (cluster_count_per_size);
	free_matrix ( (void **) cluster);
    }
    free (score);

    logger ( &options, NOTE, "");
    return 0;
}
