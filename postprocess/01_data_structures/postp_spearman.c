# include "postp.h"

double  spearman ( int * rank1, int * rank2, int length ) {
    double corr, avg1, avg2, sum_sq1, sum_sq2;
    double r, s;
    int i;
   
    avg1 = avg2 = sum_sq1 =  sum_sq2 = 0;
    for (i=0; i<  length; i++) {
	avg1 +=  rank1[i];
	avg2 +=  rank2[i];
    }
    avg1 /= length;
    avg2 /= length;
     
    corr = 0.0;
    for (i=0; i< length; i++) {
	r =  (rank1[i]-avg1);
	s = (rank2[i]-avg2);
	corr += r*s;
	sum_sq1 += r*r;
	sum_sq2 += s*s;
    }
    if (  sum_sq1*sum_sq2 < 1.e-3) {
	corr = 0.0;
    } else {
	corr /= sqrt ( sum_sq1*sum_sq2);
    }

    return corr;
}
