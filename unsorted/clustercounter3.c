# include <stdio.h>
# include <stdlib.h>
# include <limits.h>

# define PANIC( errmsg )			\
   fprintf (stderr,"%s\n" errmsg);		\
   exit(1)

void cluster_counter (int  no_of_res,  int *neighbors[], int * mask,
		     int cluster_count_per_size[],
		     int * no_of_clusters, int * max_size, int * secnd_max_size){
	

    /* arrays */ 
    static int*  * flag_ptr;          /* array of pointers to cluster flags*/
    static int   * flags;             /* array of available flags */
    static int   * bin  ;             /* for counting different clusters */
    /* counters, booleans etc */
    int flag_ctr, this_residue, other_residue, repeat_ctr, coverage_ctr, i; 
    int neighbor_found, new_flag, cluster_count, isolated_count;
    int this_value, other_value, max_cluster, second_max;

    /* for allocation purposes */
    static int first = 1;

    if ( first ) { /* do allocation */
	first = 0;
	/* flag ptrs   */
	flag_ptr     = calloc (no_of_res, sizeof(int*));
	/* flags        */
	flags        = calloc (no_of_res, sizeof(int));
	/* bins         */
	bin          = calloc (no_of_res, sizeof(int));
    }
    /* check if all alive: */ 
    if ( !( flag_ptr && flags && bin) ) {
	PANIC ("Error allocating memory in ClusterCounter."); 
    }

	
    /* set all the flags to 0 */
    memset (flags, 0, no_of_res*sizeof(int));
    /* the number of times new flag is assigned:*/
    new_flag = 0;
    /* set all the flag ptrs  to NULL */
    memset (flag_ptr, 0, no_of_res*sizeof(int*));
    /* color by cluster */ 
     
    for (this_residue=0; this_residue < no_of_res; this_residue++) {
	if (  mask [this_residue] ) {
	    for (other_residue=this_residue+1; other_residue < no_of_res; other_residue++) {
		if (  mask [other_residue] && neighbors[this_residue][other_residue]) {
		    if (flag_ptr[this_residue]){
			if (flag_ptr[other_residue]){ /*if both ptrs assigned*/
				/*************************************************/
				/* look at the flag values they are assigned to: */
			    if ( *flag_ptr[this_residue]  !=  *flag_ptr[other_residue] ) { 
				/* i.e. do something only if they differ*/
				this_value   = *flag_ptr[this_residue];
				other_value  = *flag_ptr[other_residue];
				for ( flag_ctr=0; flag_ctr < new_flag; flag_ctr++ ) {
				    if ( flags[flag_ctr] == other_value) {
					flags[flag_ctr] = this_value;
				    }
				}
				    
			    }
			} else {                       /* one not assigned*/ 
				/*************************************************/
			    flag_ptr[other_residue] = flag_ptr[this_residue];
			}
		    } else {
			if (flag_ptr[other_residue]){ /* one not assigned*/
				/*************************************************/
			    flag_ptr[this_residue]  = flag_ptr[other_residue];
			} else {                      /* both null*/
				/*************************************************/
			    /*  create new flag*/
			    flags[new_flag] = new_flag;
				/*  make both ptrs point there*/
			    flag_ptr[this_residue]  = flag_ptr[other_residue] = &flags[new_flag];
			    new_flag++;
			}
		    
		    }

		}
	    }
	}
    }

    /*count the clusters*/
    memset (bin, 0, no_of_res*sizeof(int));
    cluster_count = 0;
    isolated_count = 0;
    for (this_residue=0; this_residue < no_of_res; this_residue++) {
	if (  mask [this_residue] ) {
	    if ( !flag_ptr[this_residue] ) {
		isolated_count++; 
	    } else {
		if ( ! bin[*flag_ptr[this_residue]] ){
		    cluster_count ++;
		}
		bin[*flag_ptr[this_residue]] ++;
	    }
	}
    }
    /* find max cluster */
    if (isolated_count == 0 ) {
	second_max = max_cluster = 0;
    } else {
	second_max = max_cluster = 1;
    }
    memset ( cluster_count_per_size, 0, no_of_res * sizeof(int));
	    
    for ( flag_ctr=0; flag_ctr < new_flag; flag_ctr++ ) {
	cluster_count_per_size[ bin[flag_ctr] ] ++;
	if ( bin[flag_ctr] >= max_cluster ) {
	    second_max = max_cluster;
	    max_cluster = bin[flag_ctr];
	}
    }
    cluster_count_per_size[1] = isolated_count;
    
    /* save the count and the max cluster */
    * no_of_clusters = cluster_count+isolated_count;
    * max_size = max_cluster;
    * secnd_max_size = second_max;
    return;
}

