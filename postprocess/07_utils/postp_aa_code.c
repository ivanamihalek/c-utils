# include "postp.h"
char single_letter ( char code[]){
    char      *aa [20] = {"GLY","ALA","VAL", "LEU","ILE","MET","PRO","TRP",
		       "PHE", "TYR", "SER", "CYS","THR", "ASN", "GLN","HIS",
		       "ASP","GLU","LYS","ARG"};
    char      *a1 = "GAVLIMPWFYSCTNQHDEKR";
    int ctr;
    for (ctr=0; ctr < 20; ctr++ ) {
	if ( ! strncmp ( aa[ctr], code, 3) ) {
	    return a1[ctr];
	}
    }
    return 'X';

}

