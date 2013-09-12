# include "postp.h"

double area_over_coverage (int * int_cvg, double * value, int no_of_res ) {
     
     double area, prev_value;
     int cvg_ctr, bin_ctr, empty_bins;
     

     cvg_ctr = 0;
     prev_value = 0.0;
     empty_bins = 0;
     area = 0.0;
     for (bin_ctr = 0; bin_ctr <no_of_res; bin_ctr++ ) {
	 if ( bin_ctr == int_cvg[cvg_ctr] ) {
	     area += value[cvg_ctr];
	     if ( empty_bins) {
		 if ( prev_value < value[cvg_ctr] ) {
		     area += empty_bins*prev_value;
		 } else {
		     area += empty_bins*value[cvg_ctr];
		 }
		 empty_bins = 0;
	     }
	     prev_value =  value[cvg_ctr];
	     cvg_ctr ++;
	 } else {
	     empty_bins ++;
	 }
     }

     
     
     area /= no_of_res;
     return area;
		
}
