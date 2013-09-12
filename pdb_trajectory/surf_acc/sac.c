# include "sac.h"

# define MAX_LIGAND_TYPES 10
# define FILENAME_LENGTH  50

int main ( int argc, char * argv[]) {
    
    Protein protein[2];
    Atom * solvent = NULL;
    Atom * ligand;
    Neighbor *** nbr_list;
    Neighbor *** old_nbr_list;
    Neighbor ** tmp_nbr_list;
    int ** curr_list_length;
    int ** old_list_length;
    int max_nbr_list_length;
    int chain_number, solvent_number, done, first;
    int ctr, ctr2, frame_ctr, resctr;
    int i, j;
    int do_nbr_stats;
    int invert;
    int ligand_number, no_ligand_types;
    int special[MAX_SPECIAL_RESIDUES] = {0}, no_special_residues;
    int special_chain[MAX_SPECIAL_RESIDUES] = {0};
    int chain, residue, adjusted_list_length;
    char filename[FILENAME_LENGTH] = {'\0'};
    char ligand_type[MAX_LIGAND_TYPES][3] = {{'\0'}};
    double ** dist_to_water = {NULL};
    double ** avg_nr_of_nbrs = {NULL};
    double ** avg_id_change = {NULL};
    double score;
    double nbr_dist = 3.7;
    FILE * fptr;
    FILE * fptr_ids;
   
    int read_pdb (char *pdbname, Protein *protein, int *chain_number,
		  Atom ** solvent, int *solvent_number,   Atom ** ligand, int *ligand_number,
		  char ligand_type [][3], int no_ligand_types, int invert, int *done);
    int closest_waters (Protein * protein, Atom *solvent, int solvent_number,
			 double *dist_to_water, double nbr_dist, double * avg_nr_of_nbrs,
			 Neighbor **nbr_list, int * curr_list_length, int max_nbr_list_length);
    int flow_waters (Protein * protein, Atom *solvent, int solvent_number, int nbr_dist, 
		      Neighbor **old_nbr_list, int * old_list_length,
		      Neighbor **nbr_list, int * curr_list_length, int max_nbr_list_length,  double * avg_id_change);
 
    if ( argc < 3 ) {
	fprintf (stderr, "Usage: %s <pdb file name> <output name> [-n] [-l <lig type 1> <lig type 2> ...] [-s <chain res 1>   <chain res 2> ...]\n", argv[0]);
	exit (1);
    }

    /* options */
    do_nbr_stats = 0;
    invert = 0;
    no_ligand_types =0;
    no_special_residues =0;
    for (ctr = 3; ctr < argc; ctr++ ) {
	if ( ! strncmp (argv[ctr], "-n", 2) ) {
	    do_nbr_stats = 1;
	} else if (! strncmp (argv[ctr], "-l", 2) ) {
	    for (ctr2 = ctr+1; ctr2 < argc; ctr2++ ) {
		if ( argv[ctr2][0] == '-' ) break;
		if ( no_ligand_types >= MAX_LIGAND_TYPES ) {
		    fprintf (stderr, "Too many ligand types. Increase MAX_LIGAND_TYPES and recompile.\n");
		    exit (1);
		}
		sprintf ( ligand_type[no_ligand_types],  "%3s", argv[ctr2]);
		no_ligand_types++;
	    }
	    
	} else if (! strncmp (argv[ctr], "-s", 2) ) {
	    for (ctr2 = ctr+1; ctr2 < argc; ctr2+=2 ) {
		if ( argv[ctr2][0] == '-' ) break;
		if ( no_special_residues >= MAX_SPECIAL_RESIDUES) {
		    fprintf (stderr, "Too many special residues. Increase MAX_SPECIAL_RESIDUES and recompile.\n");
		    exit (1);
		}
		if ( ctr2+1 >= argc ) {
		    fprintf (stderr, "Chain w/o residues ?\n");
		    exit (1);
		}
		special_chain[no_special_residues] =  atoi(argv[ctr2]);
		special [no_special_residues] = atoi(argv[ctr2+1]);
		no_special_residues++;
	    }
	}
    }

    done = 0;
    first = 1;
    frame_ctr = 0;
    printf ( "\n");
    while ( !done ) {
	if ( first ) {
	    printf ("frame: %4d\n", frame_ctr+1);
	} else {
	    printf ("frame: %4d\r", frame_ctr+1);
	    fflush (stdout);
	}
	if (  read_pdb (argv[1], protein, &chain_number, &solvent, &solvent_number, &ligand, &ligand_number,
			ligand_type, no_ligand_types, invert,  &done) ){
	    fprintf (stderr, "Error reading %s.\n", argv[1] );
	    exit (1);
	}
	if (done) break;
	if ( first ) {
	    max_nbr_list_length = MAX_FRAMES*MAX_NBRS;
	    if ( ! ( dist_to_water = emalloc( chain_number*sizeof(double *)) )) exit(1);
	    if ( ! ( avg_nr_of_nbrs = emalloc( chain_number*sizeof(double *)) )) exit(1);
	    if ( ! ( avg_id_change = emalloc( chain_number*sizeof(double *)) )) exit(1);
	    if ( ! ( nbr_list  = emalloc(chain_number*sizeof(Neighbor **)) )) exit(1);
	    if ( ! ( curr_list_length  = emalloc(chain_number*sizeof(int *)) )) exit(1);
	    if ( ! ( old_nbr_list  = emalloc(chain_number*sizeof(Neighbor **)) )) exit(1);
	    if ( ! ( old_list_length  = emalloc(chain_number*sizeof(int *)) )) exit(1);
	    for (i=0; i< chain_number; i++ ) {
		/* space for dist to water */
		if ( ! ( dist_to_water[i]  = emalloc((protein+i)->length*sizeof(double)) )) exit(1);
		if ( ! ( avg_nr_of_nbrs[i] = emalloc((protein+i)->length*sizeof(double)) )) exit(1);
		if ( ! ( avg_id_change[i] = emalloc((protein+i)->length*sizeof(double)) )) exit(1);
		/* space for neighbor list */
		if ( ! ( nbr_list[i] = emalloc((protein+i)->length*sizeof(Neighbor*)) )) exit(1);
		if ( ! ( old_nbr_list[i] = emalloc((protein+i)->length*sizeof(Neighbor*)) )) exit(1);
		for ( j=0; j< (protein+i)->length; j++ ) {
		    if ( ! (nbr_list[i][j] =   emalloc(max_nbr_list_length*sizeof(Neighbor) )) )  exit(1);
		    if ( ! (old_nbr_list[i][j] =   emalloc(max_nbr_list_length*sizeof(Neighbor) )) )  exit(1);
		}
		if ( ! ( curr_list_length[i] = emalloc((protein+i)->length*sizeof(int)) )) exit(1);		
		if ( ! ( old_list_length[i]  = emalloc((protein+i)->length*sizeof(int)) )) exit(1);		
	    }
	}
	

	for (i=0; i< chain_number; i++ ) {
	    closest_waters (protein+i, solvent, solvent_number, dist_to_water[i],  nbr_dist, avg_nr_of_nbrs[i],
	    		 nbr_list[i], curr_list_length[i],  max_nbr_list_length);

# if 0
	    flow_waters (protein+i, solvent, solvent_number, nbr_dist, 
			 old_nbr_list[i], old_list_length[i],
			 nbr_list[i], curr_list_length[i], max_nbr_list_length,  avg_id_change[i]);
	    tmp_nbr_list = old_nbr_list[i];
	    old_nbr_list[i] = nbr_list[i];
	    nbr_list[i] = tmp_nbr_list ;

	    old_list_length[i] = curr_list_length[i];
	    if (  first ) memset ( avg_id_change[i], 0, (protein+i)->length*sizeof(double) );
# endif

	}
	
	
	first = 0;
	if (no_special_residues)  {
	    printf (" %4d", frame_ctr + 1);
	    for (i=0; i<no_special_residues; i++ ) {
		chain   = special_chain[i] - 1;
		residue = special[i] - 1;
		printf ("\t %2d  %4d", special_chain[i],  special[i]);
		for (j=0; j <  curr_list_length[chain][residue]; j++) {
		    printf (" %4d",  nbr_list[chain][residue][j].id);
		}
		printf ("\n");
	    }
	    printf ("\n\n");
	}



	
	frame_ctr++;

	//done = (frame_ctr==5);
    }
    printf ( "\n");
    if ( done && ! (protein->length + (protein+1)->length) ) {
	fprintf ( stderr, "Error reading pdb: empty file?\n");
	exit (1);
    }


    /****** get rid of this for heterodimers: ******/
    //chain_number = 1;
    
    
    /* average  */
   for (i=0; i< chain_number; i++ ) {
	for ( resctr=0; resctr<(protein+i)->length; resctr++) {
	    dist_to_water[i][resctr]  /= frame_ctr;
	    avg_nr_of_nbrs[i][resctr] /= frame_ctr;
	    if ( frame_ctr > 1) avg_id_change[i][resctr] /=  (frame_ctr-1);
	}
   }
     
   for (i=0; i< chain_number; i++ ) {
       sprintf ( filename, "%s_acc.%d.table", argv[2],  i+1);
       //sprintf ( filename, "water_flow.%d.table", i+1);
	if ( ! (fptr = efopen ( filename, "w")) ) exit (1);
	
       sprintf ( filename, "%s_ids.%d.table",argv[2],  i+1);
	if ( ! (fptr_ids = efopen ( filename, "w")) ) exit (1);
	for ( resctr=0; resctr< (protein+i)->length; resctr++) {

	    
	    //fprintf ( fptr,"  %5s %5d  %c   %8.3lf ", (protein+i)->sequence[resctr].pdb_id, resctr+1,
	    
	    //(protein+i)->sequence[resctr].res_type_short, avg_id_change[i][resctr] );
	    adjusted_list_length  =  curr_list_length[i][resctr];
	    score = 0.0;
	    for (j=0; j <  curr_list_length[i][resctr]; j++) {
	    	if(  (double) nbr_list[i][resctr][j].count/frame_ctr > 0.8 ) {
		   adjusted_list_length --;
	    	}
		score += 1.0/ nbr_list[i][resctr][j].count;
	    } 
	    score /= frame_ctr;
	    
	    fprintf ( fptr,"  %5s %5d  %c  %8.3lf   %8.3lf  %8.3lf  %8.3lf   %4d  ",
		      (protein+i)->sequence[resctr].pdb_id, resctr+1,
	    	     (protein+i)->sequence[resctr].res_type_short, dist_to_water[i][resctr],  avg_nr_of_nbrs[i][resctr],
	    	     score, (double)adjusted_list_length/frame_ctr,  curr_list_length[i][resctr]);
	    
	    for (j=0; j <  curr_list_length[i][resctr]; j++) {
	    	if(  (double) nbr_list[i][resctr][j].count/frame_ctr > 0.3 ) {
	    	    fprintf ( fptr, " %3d", nbr_list[i][resctr][j].count);
	    	}
	    } 
	    fprintf ( fptr, "\n");

	    /* water ids - for illustration */ 
	    for (j=0; j <  curr_list_length[i][resctr]; j++) {
	    	if(  (double) nbr_list[i][resctr][j].count/frame_ctr <  0.3 ) {
	    	    fprintf ( fptr_ids, " %6d  %s\n", nbr_list[i][resctr][j].id, (protein+i)->sequence[resctr].pdb_id);
	    	}
	    } 
	     
	}  
	fclose (fptr);
	fclose (fptr_ids);
    } 




   
    
    return 0; 
} 
