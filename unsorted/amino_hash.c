# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define PRIM 127

int hashfun ( char * aa) {
    return  (aa[0]*7+ aa[1]*3+aa[2]-748)/2;
}

int main ( int argc, char * argv[]) {

    char *aa [20] = {"GLY","ALA","VAL", "LEU","ILE","MET","PRO","TRP",
		       "PHE", "TYR", "SER", "CYS","THR", "ASN", "GLN","HIS",
		       "ASP","GLU","LYS","ARG"};
    char *a1 = "GAVLIMPWFYSCTNQHDEKR";
    char translation [100];
    int i, val;

# if 1
    for (i=0; i<20; i++) {
	val = hashfun (aa[i]);
	translation[val] = a1[i];
    }
# endif
    
    printf ("\n");
    for (i=0; i<20; i++) {
	val = hashfun (aa[i]);
	printf ("\t%2d  %3s   %5d %1c\n", i, aa[i], val,translation[val]  );
    }
    printf ("\n");
    
    return 0;
}
