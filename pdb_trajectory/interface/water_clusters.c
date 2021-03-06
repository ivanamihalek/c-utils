# include "ift.h"

int water_clusters ( Atom* solvent, int solvent_number, int *if_solvent, double cutoff_dist, double * neighbors) {

    int  c1, c2;
    double cutoff_dist_sq = cutoff_dist*cutoff_dist;
    double dist, aux;
    Atom * atomptr1, * atomptr2;


    /* neighboring waters */
    for ( c1=0; c1 < solvent_number; c1++) {
	if ( ! if_solvent[c1] ) continue;
	atomptr1 = solvent + c1;
	
	for ( c2=c1+1; c2 <  solvent_number; c2++) {
	    if ( ! if_solvent[c2] ) continue;
	    atomptr2 = solvent + c2;
		    
	    aux = atomptr1->x - atomptr2->x;
	    dist = aux*aux;
	    if (dist >  cutoff_dist_sq) continue;
	    aux = atomptr1->y - atomptr2->y;
	    dist += aux*aux;
	    if (dist >  cutoff_dist_sq) continue;
	    aux = atomptr1->z - atomptr2->z;
	    dist += aux*aux;
	    if (dist >  cutoff_dist_sq) continue;
	    
	    neighbors[c1] ++;
	    neighbors[c2] ++;
	    
	    
	}
    }

    
		 
    return 0;
}

/***************************************************************************/
/***************************************************************************/
/* I can't do that - too many waters */
void cluster_counter (int  no_of_things,  int *neighbors[], 
		     int cluster_count_per_size[], int * no_of_clusters,
		      int * max_size, int * secnd_max_size , int * clusters[]){
	

    /* arrays */ 
    static int*  * flag_ptr;          /* array of pointers to cluster flags*/
    static int   * flags;             /* array of available flags */
    static int   * bin  ;             /* for counting different clusters */
    /* counters, booleans etc */
    int flag_ctr, this_thing, other_thing; 
    int new_flag, cluster_count, isolated_count;
    int this_value, other_value, max_cluster, second_max;
    int color;

    /* for allocation purposes */
    static int first = 1;

    if ( first ) { /* do allocation */
	first = 0;
	/* flag ptrs   */
	flag_ptr     = calloc (no_of_things, sizeof(int*));
	/* flags        */
	flags        = calloc (no_of_things, sizeof(int));
	/* bins         */
	bin          = calloc (no_of_things, sizeof(int));
    }
    /* check if all alive: */ 
    if ( !( flag_ptr && flags && bin) ) {
	fprintf ( stderr, "Error allocating memory in ClusterCounter.");
	exit (1);
    }

	
    /* set all the flags to 0 */
    memset (flags, 0, no_of_things*sizeof(int));
    /* the number of times new flag is assigned:*/
    new_flag = 0;
    /* set all the flag ptrs  to NULL */
    memset (flag_ptr, 0, no_of_things*sizeof(int*));
    /* color by cluster */ 
     
    for (this_thing=0; this_thing < no_of_things; this_thing++) {
	for (other_thing=this_thing+1; other_thing < no_of_things; other_thing++) {
	    if (neighbors[this_thing][other_thing]) {
		if (flag_ptr[this_thing]){
		    if (flag_ptr[other_thing]){ /*if both ptrs assigned*/
			/*************************************************/
			/* look at the flag values they are assigned to: */
			if ( *flag_ptr[this_thing]  !=  *flag_ptr[other_thing] ) { 
			    /* i.e. do something only if they differ*/
			    this_value   = *flag_ptr[this_thing];
			    other_value  = *flag_ptr[other_thing];
			    for ( flag_ctr=0; flag_ctr < new_flag; flag_ctr++ ) {
				if ( flags[flag_ctr] == other_value) {
				    flags[flag_ctr] = this_value;
				}
			    }
				    
			}
		    } else {                       /* one not assigned*/ 
			/*************************************************/
			flag_ptr[other_thing] = flag_ptr[this_thing];
		    }
		} else {
		    if (flag_ptr[other_thing]){ /* one not assigned*/
			/*************************************************/
			flag_ptr[this_thing]  = flag_ptr[other_thing];
		    } else {                      /* both null*/
			/*************************************************/
			/*  create new flag*/
			flags[new_flag] = new_flag;
			/*  make both ptrs point there*/
			flag_ptr[this_thing]  = flag_ptr[other_thing] = &flags[new_flag];
			new_flag++;
		    }
		    
		}

	    }
	}
    }

    /*count the clusters*/
    memset (bin, 0, no_of_things*sizeof(int));
    memset (clusters[0], 0, (no_of_things+1)*(no_of_things+1)*sizeof(int));
    cluster_count = 0;
    isolated_count = 0;
    for (this_thing=0; this_thing < no_of_things; this_thing++) {
	if ( !flag_ptr[this_thing] ) {
	    isolated_count++;
	    clusters [0][0]++;
	    clusters [0][ clusters [0][0] ] = this_thing;
	} else {
	    color = *flag_ptr[this_thing];
	    if ( ! bin[color] ){
		cluster_count ++;
	    }
	    bin[color] ++;
	    color += 1;
	    clusters [color][0]++;
	    clusters [color][ clusters [color][0] ] = this_thing;
	}
	
    }



    /* find max cluster */
    if (isolated_count == 0 ) {
	second_max = max_cluster = 0;
    } else {
	second_max = max_cluster = 1;
    }
    memset ( cluster_count_per_size, 0, no_of_things * sizeof(int));
	    
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

