# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>


int main ( int argc, char * argv[]) {
    printf ("\n\n");
    printf ( "\t  %10s  %4d \n", "char", sizeof(char) );
    printf ( "\t  %10s  %4d \n", "int", sizeof(int) );
    printf ( "\t  %10s  %4d \n", "long", sizeof(long long) );
    printf ( "\t  %10s  %4d \n", "double", sizeof(double) );
    printf ( "\t  %10s  %4d \n", "int*", sizeof(int*) );
    printf ("\n\n");
    return 0;
}
