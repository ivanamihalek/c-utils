/* turn text (data) file into ppm file */
/* note that the matrix values are assumed POSITIVE */ 
#include "cpgplot.h"

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <limits.h>
# include <math.h>
# define MATRIX_FILE "matrix.dat" /* input data are here*/
# define PPM_FILE "matrix.ppm"     /* output file name */
# define MATRIX_FORMAT "%d %d %lf" /* format in the input file */ 

# define COLOR_RANGE 100
# define WIDTH  750
# define HEIGHT 750


# define PANIC(msg, str) {			\
    fprintf (stderr, "%s %s.\n",msg, str);	\
    exit(1);					\
}
/*  255 255 255 is white; 0 0 0 is black; 255 0 0 is red; 0 255 0 is green; 0 0 255 is blue.  */
void palett (int* TYPE, float* CONTRA, float*  BRIGHT);

main (int argc, char * argv[]) {

    FILE * matrfile, *ppmfile;

    int x,y, x_offset, y_offset, x_old, y_old, x_size, y_size;
    float ** matrix, tfm[6];
    int width = WIDTH,height = HEIGHT, red[COLOR_RANGE+1] = {0};
    int green[COLOR_RANGE+1]= {0}, blue[COLOR_RANGE+1]= {0};
    int raster[HEIGHT][WIDTH];
    int  first, value, ctr, linectr;
    double xstep, ystep, dummy, max, min;
    char infilename [5] = {'\0'};
    int type;
    float contra, bright;

    /* *****************************************************/
    /* ******** INPUT;  INITIALIZATION         *************/
    /* *****************************************************/
    /*read in the datafile name */
    if ( argc>1) {
	sprintf (infilename,"%s", argv[1]);
    } else {
	sprintf (infilename,"%s", MATRIX_FILE);
    }
    /* see if the datafile there; open */
    if ( ! (matrfile= fopen ( infilename, "r" )  ) ) {
	PANIC ("Cno ", infilename);
    }

    /* max value: */ 
    if ( argc >2) { /* read in the max value */
	max = atof (argv[2]);
    } else {        /* max will be determined from the input data*/ 
	max = -1; 
    }

    /* min value: */ 
    if ( argc >3) { /* read in the max value */
	min = atof (argv[3]);
    } else {        /* max will be determined from the input data*/ 
	min = 101; 
    }

    
    /* find out the matrix size*/
    first = 1;
    x_size = y_size =0;
    x_offset = y_offset = INT_MAX;
    linectr = 0;
    while ( ! feof (matrfile)) {
	
	if ( fscanf (matrfile, MATRIX_FORMAT, &x, &y, &dummy) == 3 ) {
	    linectr ++;
	    if ( x > x_size ){
		x_size = x;
	    }
	    if ( x < x_offset ){
		x_offset = x; 
	    }
	    if ( y > y_size ){
		y_size = y;
	    }
	    if ( y <  y_offset ){
		y_offset = y;
	    }
	}
	
    }
    x_size -= (x_offset-1);
    y_size -= (y_offset-1);

    /* allocate matrix space */
    if ( ! (matrix = malloc ( x_size*sizeof(float*) ) ) ){
	PANIC ("mem alloc failure","");
    }
    for (x = 0; x < x_size; x++)  {
	if ( ! (matrix[x] = calloc ( y_size, sizeof(float) ) ) ){
	    PANIC ("mem alloc failure","");
	}
    }
    rewind (matrfile);

    /* read in the data */
    while ( ! feof (matrfile)) {
	if ( fscanf (matrfile, MATRIX_FORMAT, &x, &y, &dummy) == 3 ) {
	    x -= x_offset;
	    y -= y_offset;
	    /* distance plain */
	    matrix[x][y] = dummy;
	}
    }
    fclose (matrfile);
    
    /* determine max, if needed: */
    if ( max < 0 ) {
	for (x = 0;x < x_size;x++)  {
	    for (y = 0;y < y_size;y++)  {
		if ( matrix [x][y] > max ) {
		    max = matrix[x][y];
		}
	    }
	}
    }
    if ( min > 100 ) {
	for (x = 0;x < x_size;x++)  {
	    for (y = 0;y < y_size;y++)  {
		if ( matrix [x][y] <  min ) {
		    min = matrix[x][y];
		}
	    }
	}
    }
    printf ("max is %lf \n", max);

    

    /*actual /plotting: */
   if(cpgbeg(0, "/xw", 1, 1) != 1)
	exit(EXIT_FAILURE);
    cpgask(1);

    tfm[0] = 0; tfm[1] = 1; tfm[3] = 0; 
    tfm[3] = 0; tfm[4] = 0; tfm[5] = 1;

    cpgpage();
    cpgvstd();
    cpgwnad (-1.0, 1.0, -1.0, 1.0);
    cpgsci(1);
    type = 2; contra = 0.5; bright = 1.0;
    palett_ (&type, &contra, &bright );
    cpgimag ( matrix[0], x_size, y_size, 1, x_size, 1, y_size, min, max, tfm );
    cpgend();


    return 0;  

 }



