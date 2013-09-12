# include <math.h>
# include <stdio.h>

int main() {

    int x, y;
    for (x=100; x>=0; x--) {
	    printf ("%3d %3d %15.4lf\n", 100-x, 1, (double)x );
    }
    return 0;
}
