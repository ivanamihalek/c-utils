/* turn text (data) file into ppm file */
/* note that the matrix values are assumed POSITIVE */ 
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

main (int argc, char * argv[]) {

    FILE * matrfile, *ppmfile;

    int x,y, x_offset, y_offset, x_old, y_old, x_size, y_size;
    double ** matrix;
    int width = WIDTH,height = HEIGHT,  red_space[2*COLOR_RANGE+1] = {0};
    int green_space[2*COLOR_RANGE+1]= {0}, blue_space[2*COLOR_RANGE+1]= {0};
    int raster[HEIGHT][WIDTH];
    int *red, *blue, *green;
    int  first, value, ctr;
    double xstep, ystep, dummy, max, min, large;
    char infilename [5] = {'\0'};


    red =  &red_space[COLOR_RANGE];
    blue =  &blue_space[COLOR_RANGE];
    green = & green_space[COLOR_RANGE];
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
	max = 0; 
    }
    min=0;

    
    /* find out the matrix size*/
    first = 1;
    x_size = y_size =0;
    x_offset = y_offset = INT_MAX;
    while ( ! feof (matrfile)) {
	
	if ( fscanf (matrfile, MATRIX_FORMAT, &x, &y, &dummy) == 3 ) {
/*	    printf ( " %d  %d  %lf\n", x, y, dummy);	*/

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
/*		printf ( "\told: %d  new: %d \n", y_offset, y);  */
		y_offset = y;
	    }
	}
	
    }
    x_size -= (x_offset-1);
    y_size -= (y_offset-1);

    /* allocate matrix space */
    if ( ! (matrix = malloc ( x_size*sizeof(double*) ) ) ){
	PANIC ("mem alloc failure","");
    }
    for (x = 0; x < x_size; x++)  {
	if ( ! (matrix[x] = calloc ( y_size, sizeof(double) ) ) ){
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
    
    /* determine max: */
	for (x = 0;x < x_size;x++)  {
	    for (y = 0;y < y_size;y++)  {
		if ( matrix [x][y] > max ) {
		    max = matrix[x][y];
		}
		if ( matrix [x][y] < min ) {
		    min = matrix[x][y];
		}
	    }
	}
    printf ("max is %lf \n", max);
    printf ("min is %lf \n", min);



    /* *****************************************************/
    /* ******* COLOR  CALCULATION                    ******/
    /* *****************************************************/

    /* spread the matrix on the raster array */
    xstep = (double)HEIGHT/x_size;
    ystep = (double)WIDTH/y_size;
   
    large=max;
    if (-min>max)
	    large=-min;

    for (x = 0;x < height;x++)  {
	for (y = 0;y < width;y++)  {
	    raster [x][y] = matrix[(int)(x/xstep)][(int) (y/ystep)]
		/large*COLOR_RANGE;
	}
    }
    
    /* set the pallette: */
    
    for ( ctr=0; ctr <= COLOR_RANGE; ctr++ ) {
	red[ctr]=255;
	green[ctr]=blue[ctr] = 255-255*ctr/COLOR_RANGE;
    }

     for ( ctr= -1; ctr >=  -COLOR_RANGE; ctr-- ) {
	blue[ctr]=255;
	green[ctr]=red[ctr] = 255+255*ctr/COLOR_RANGE;
/*printf ("ctr red green blue = %d %d %d %d\n", ctr, red[ctr], green[ctr], blue[ctr]); */
    }

    /* *****************************************************/
    /* ******* OUTPUT                                 ******/
    /* *****************************************************/
    /* spread the matrix on the raster array */
   
    if ( ! (ppmfile= fopen ( PPM_FILE, "w" )  ) ) {
	PANIC ("Cno ", PPM_FILE);
    }
    fprintf(ppmfile, "P6\n");
    fprintf(ppmfile,"# created by Ivana using a C program\n");
    fprintf(ppmfile,"%d %d\n", width+2, height+2);
    fprintf(ppmfile,"255\n");

   
    
    /*frame:*/
    x = 0; {
        for ( y=0; y < width+2; y++)  {
	    fputc(0,ppmfile);
	    fputc(0,ppmfile);
	    fputc(0,ppmfile);
       }
    }
    for (x = 0; x < height;x++)  {
	/*frame:*/
	fputc(0,ppmfile);
	fputc(0,ppmfile);
	fputc(0,ppmfile);
	/*pic:*/
        for ( y=0; y < width; y++)  {
	    value = raster[x][y];
	    fputc((char)red[value],ppmfile);
	    fputc((char)green[value],ppmfile);
	    fputc((char)blue[value],ppmfile);
       }
	/*frame:*/
	fputc(0,ppmfile);
	fputc(0,ppmfile);
	fputc(0,ppmfile);
    }
    x = height; {
    /*frame:*/
        for ( y=0; y < width+2; y++)  {
	    fputc(0,ppmfile);
	    fputc(0,ppmfile);
	    fputc(0,ppmfile);
       }
    }

    fclose (ppmfile);
    
    return 0;
}
