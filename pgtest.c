# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <unistd.h>
# include <cpgplot.h>

# define N 5

int main ( int argc, char * argv[]) {
    float xs[N] = {1.0, 2.0, 3.0, 4.0, 5.0};
    float ys[N];
    int ctr;
    
    printf ("\nThe  env variable PGPLOT_DIR must be set. \n\n");
    cpgbeg (0,"?",1,1);
    cpgenv (0.0, 10.0, 0.0, 20.0, 0,1);
    cpglab ("x", "y", "title");
    for (ctr=0; ctr <N; ctr++) {
	ys[ctr] = xs[ctr]*xs[ctr]/2;
    }
    cpgpt ( N, xs, ys, 9);
    
    cpgline ( N, xs, ys);
    
    cpgsci (2);
    cpgpt1 (xs[N-1], ys[N-1], 9);
    
    cpgend ();
    
    return 0;
}
