

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>


int main ( int argc, char * argv[]) {


    FILE * fptr, *fptrout = stdout;
    int ctr, aux1, aux2, aux3;
    char buf[150];

    if  ( argc <2 ) exit(1);
    fptr = fopen ( argv[1], "r" ) ;
    if ( !fptr) exit(1);

    ctr = 0;
    while ( !feof (fptr) && ctr<4) {
	fgets (buf, 150, fptr);
	//fprintf (fptrout, "%s", buf);
	ctr++;
    }
    fprintf(fptrout, "P6\n");
    fprintf(fptrout,"# created by Ivana using a C program\n");
    fprintf(fptrout,"1240  383\n");
    fprintf(fptrout,"255\n");
    while ( !feof (fptr) ) {
	fscanf (fptr, "%d%d%d", &aux1, &aux2, &aux3);
	aux1 = 1;
	aux2 =  100;
	aux3 = 200;
	//printf ( "%d  %d  %d \n", aux1, aux2, aux3);
#if 1
	fputc ( (char) aux1, fptrout);
	fputc ( (char) aux2, fptrout);
	fputc ( (char) aux3, fptrout);
# else 
	fputc ( 0, fptrout);
	fputc ( 0, fptrout);
	fputc ( 0, fptrout);
# endif
    }
    return 0;
}
