# include <math.h>
# include <stdio.h>

int main() {

    int x, y;
    for (x=0; x<=10; x++) {
	for (y=0; y<=10; y++) {
	    printf ("%3d %3d %15.4lf\n", x, y, sin(M_PI/10*x)*sin(M_PI/10*y) );
	}
    }
    return 0;
}
