# include "postp.h"

int output_clustering (char * base_filename,  int * int_cvg, double *clustering_score, int length){

    int ctr;
    FILE * fptr;
    char filename[BUFFLEN];
    
    sprintf (filename, "%s.clustering", base_filename);

    fptr = efopen (filename, "w");
    if ( !fptr) return 1;
    
    fprintf ( fptr, "%%%8s%8s%12s\n", "cum res", "cvg ", "z  ");
  
    for (ctr=0; ctr < length &&int_cvg[ctr] ; ctr++ ) {
	fprintf ( fptr, "%8d%8.3lf%12.3le\n",
		  int_cvg[ctr], (double)int_cvg[ctr]/length, clustering_score[ctr]);
	
    }
    fclose (fptr);
    
    return 0;
}
