# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"


# define BUFLEN  300
# define MAX_NAME_LEN 100

int  read_distf (char distfname[100], char *** names_ptr, double *** distmat_ptr,int * no_names_ptr );
void cluster_counter (int  no_of_things,  int *neighbors[], int * mask,
		     int cluster_count_per_size[], int * no_of_clusters,
		      int * max_size, int * secnd_max_size , int * clusters[]);

int main ( int argc, char * argv[]) {
    
    char distfname[100] = {'\0'};
    double cutoff_dist = 0.0;
    char  ** name = NULL;
    int number_of_names;
    double ** distmat;
    int ** cluster;
    int ctr, ctr1, ctr2;
    int c;
    FILE * fclust = NULL;
    
    if ( argc < 3 ) {
	fprintf (stderr, "Usage: %s <dist_file> <cutoff_sim> \n", argv[0]);
	exit (1);
    }
    sprintf ( distfname, "%s", argv[1]);
    cutoff_dist = atof ( argv[2]);

    /* input seqs */
    /* first line must be the number of different names in the file */
    if ( read_distf ( distfname, &name, &distmat, &number_of_names ) ) {
	fprintf ( stderr, "Error reading distance matrix.\n");
	exit(1);
    }
    printf ("There are %d names in %s.\n",  number_of_names, distfname);


    /* cluster counting ... */
    {
	int  no_of_clusters, max_size, secnd_max_size;
	int * cluster_count, *mask;
	int ** neighbors;
	    
	cluster_count       =  (int *) emalloc ( (number_of_names+1)*sizeof(int));
	mask                =  (int *) emalloc ( (number_of_names+1)*sizeof(int));
	cluster             =  intmatrix ( number_of_names+1,  number_of_names+1);
	neighbors           =  intmatrix ( number_of_names, number_of_names);
	
	for (ctr1=0;  ctr1 < number_of_names; ctr1++ ) {
	    neighbors [ctr1][ctr1] = 1;
	    for (ctr2= ctr1+1; ctr2 < number_of_names; ctr2++ ) {
		neighbors[ctr1][ctr2] = ( distmat[ctr1][ctr2] < cutoff_dist);
		neighbors[ctr2][ctr1] = neighbors[ctr1][ctr2];
	    }
	}

	for (ctr1=0;  ctr1 < number_of_names; ctr1++ ) {
	    mask[ctr1] = 1;
	}
	cluster_counter (number_of_names,  neighbors,  mask, cluster_count, & no_of_clusters,
			 &max_size, &secnd_max_size , cluster);
    }

    /* output */
    fclust = stdout;
    for ( c=0; c <= number_of_names; c++) {
	if ( ! cluster[c][0] ) {
	    continue;
	}
	if ( !c ) {
	    fprintf ( fclust,"\t isolated:\n");
	} else {
	    fprintf ( fclust,"\t cluster size: %3d \n", cluster[c][0]); 
	}
	for ( ctr=1; ctr <=  cluster[c][0]; ctr++) {
	    fprintf ( fclust, "%s  \n", name [ cluster[c][ctr] ]  );
	}
	
    }

    
    return 0;

    
    
}



int  read_distf (char  distfname[100], char *** names_ptr, double *** distmat_ptr, int * no_names_ptr ) {

    char ** name;
    FILE * fptr = NULL;
    char line [BUFLEN] = {'\0'};
    int no_names, pos1, pos2;
    char name1[MAX_NAME_LEN], name2[MAX_NAME_LEN];
    double **distmat, dist;
    int  find_pos (char ** name, char * name1, int  no_names,  int * pos_ptr);

    
    fptr = fopen (distfname, "r");
    if ( !fptr) {
	fprintf (stderr, "Cannot open %s.\n", distfname);
	return 1;
    }
    if (fgets(line,BUFLEN,fptr) != NULL){
	sscanf (line, "%d", &no_names);
    } else {
	return 0;
    }
    name = chmatrix (no_names, MAX_NAME_LEN);
    distmat = dmatrix ( no_names, no_names);

    
    while(fgets(line,BUFLEN,fptr)!=NULL){
	sscanf (line, "%s  %s  %lf", name1, name2, &dist);

	/* process first name */
	if ( find_pos (name, name1, no_names, &pos1) ) return 1;
	
	/* process second  name */ 
	if ( find_pos (name, name2, no_names, &pos2) ) return 1;

	distmat[pos1][pos2] = distmat[pos2][pos1] = dist;

	
    }

    
    *names_ptr    = name;
    *no_names_ptr = no_names;
    *distmat_ptr = distmat;
    return 0;
}



int  find_pos (char ** name, char * name1, int  no_names,  int * pos_ptr) {
    
    int pos;
    int ctr;
    ctr = 0;
    pos = -1;
    while ( ctr < no_names && name[ctr][0] ) {
	if ( ! strcmp ( name[ctr], name1) ) {
	    pos = ctr;
	    break;
	}
	ctr++;
    }
    if ( pos < 0 ) {
	if ( ctr >= no_names ) {
	    fprintf ( stderr,
		      "Error: number of names larger than declared (%s).\n", name1);
	    return 1;
	}
	sprintf ( name[ctr], "%s", name1 );
	pos = ctr;
    }

    *pos_ptr = pos;
    
    return 0;
}
