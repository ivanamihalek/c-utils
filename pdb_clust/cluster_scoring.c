# include "pdbclust.h"
/************************************************************************************/
int cluster_score (int no_of_res, int *seq, int ** adj_matrix,double *score){
    int i,j, dist;
    int  size;
    double sum;

    size = no_of_res;
    sum = 0.0;
    for (i=0; i<size-1 ; i++) {
	if ( seq[i] ) {
	    for (j=i+1; j<size; j++) {
		if (seq[j] && adj_matrix[i][j]) {
		    dist = j-i;
		    sum += dist;
		}
	    }
	}
    }
    *score = sum;
    return 0;    
}

/************************************************************************************/
int std_dev_over_S (int L, int M, int ** adj_matrix, double *avg, double * std_dev, int first) {
    
    int i,j,k,l,n;
    double  std_dev_thry, avg_thry;
    double ratio[3];
    double aux;
    static double subsum[3], bare_avg_thry;

    if (first) {
	bare_avg_thry = 0;
	subsum[0] = subsum[1] = subsum[2] = 0.0;
	for (i=0; i<L-1; i++) {
	    for (j=i+1; j<L; j++) {
			
		if ( adj_matrix[i][j]) {
			    
		    bare_avg_thry +=  (j-i);
			    
		    for (k=0; k<L-1; k++) {
			for (l=k+1; l<L; l++) {
				    
			    if  ( adj_matrix[k][l]) {
				n = (k==i) + (l==j) + (k==j) + (l==i);
				subsum [2-n] += (j-i)*(l-k);
			    }
			}
		    }
		}
	    }
	}
        bare_avg_thry /= (L*(L-1));
    }
    
    ratio[0] = (double) M*(M-1)/(L*(L-1));
    ratio[1] = ratio[0]*(M-2)/(L-2);
    ratio[2] = ratio[1]*(M-3)/(L-3);
    std_dev_thry = 0.0;
    for (i =0; i<3; i++) {
	std_dev_thry += subsum[i]*ratio[i];
    }
   
    avg_thry      = bare_avg_thry*M*(M-1);

    aux = (std_dev_thry-avg_thry*avg_thry)/std_dev_thry;
    if ( aux  < -1.0e-4){
	fprintf (stderr,"Unspecified error in std dev (S) calculation.\n");
	fprintf (stderr, " %8.3e  %8.3e  %8.3e \n", std_dev_thry, avg_thry*avg_thry, aux);
	return 1;
    } else if ( aux < 0 ) {
	*std_dev = 0;
    } else {
	std_dev_thry -= avg_thry*avg_thry;
	*std_dev = sqrt(std_dev_thry);
    }
    *avg = avg_thry;

    return 0;
    
}

