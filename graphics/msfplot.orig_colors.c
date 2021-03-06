/* turn text (data) file into ppm file */
/* note that the matrix values are assumed POSITIVE */ 
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <limits.h>
# include <math.h>
# include <ctype.h>
# define MATRIX_FILE "msf.msf" /* input data are here*/
# define PPM_FILE "msf.ppm"     /* output file name */

# define MATRIX_FORMAT "%*d %*lf %le %le %le" /* format in the input file */ 

# define COLOR_RANGE 127 /* ascii table size */
# define WIDTH  1000
# define HEIGHT 750


# define BUFSIZE 150


# define PANIC(msg, str) {			\
    fprintf (stderr, "%s %s.\n",msg, str);	\
    exit(1);					\
}
/*  255 255 255 is white; 0 0 0 is black; 255 0 0 is red; 0 255 0 is green; 0 0 255 is blue.  */
#define MAXSTRING  90
#define LONGSTRING 300
typedef char    bigstring[LONGSTRING];
typedef char    string[MAXSTRING];
typedef struct  {
    string        name;
    int           length;
    char          *position;
    int           marked; /* can be marked for any purpose */
} Sequence;

int read_msf (char *filename, Sequence ** alignment_ptr, int * no_seqs_ptr) ;

int main (int argc, char * argv[]) {

    FILE  *ppmfile;
    Sequence * sequence;
    int no_seqs;
    int x,y,  x_size, y_size;
    int width = WIDTH,height = HEIGHT, red[COLOR_RANGE+1] = {0};
    int green[COLOR_RANGE+1]= {0}, blue[COLOR_RANGE+1]= {0};
    int raster[HEIGHT][WIDTH];
    int value, ctr;
    double xstep, ystep;
    char infilename [5] = {'\0'};
 
    /* *****************************************************/
    /* ******** INPUT;  INITIALIZATION         *************/
    /* *****************************************************/
    /*read in the datafile name */
    if ( argc>1) {
	sprintf (infilename,"%s", argv[1]);
    } else {
	sprintf (infilename,"%s", MATRIX_FILE);
    }
    /* read in the msf */
    if (  read_msf (infilename, &sequence,  &no_seqs)   ) {
	PANIC ("Error processing ", infilename);
    }


    
    /* the matrix size*/
    x_size = no_seqs; /* each column 10 pixels wide */
    y_size = sequence[0].length;
    printf ( "%d  %d \n", x_size, y_size );




    /* *****************************************************/
    /* ******* COLOR  CALCULATION                    ******/
    /* *****************************************************/
    
    

    /* spread the alignment on the raster array */
    xstep = (double)HEIGHT/x_size;
    ystep = (double)WIDTH/y_size;
   
    for (x = 0;x < height;x++)  {
	for (y = 0;y < width;y++)  {
	    raster [x][y] = sequence[(int) (x/xstep)].position[(int)(y/ystep)];
	}
    }
    
    /* set the pallette: */
    for ( ctr=0; ctr <= COLOR_RANGE; ctr++ ) {
	green[ctr] = blue [ctr] = red [ctr] = 255; /* white, unless explicitly mad otherwise */
    }
    /* black to greyish-green  F, Y, W*/
    red['W'] =  0;  green['W'] =  0;  blue['W'] =  0;
    red['Y'] = 85;  green['Y'] = 85;  blue['Y'] = 85;
    red['F'] = 191;  green['F'] = 191;  blue['F'] = 191;

    /* dark green - C */ 
    red['C'] = 0;  green['C'] = 140;  blue['C'] = 0;
    
    /* dirty yellow - M */
    red['M'] = 200;  green['M'] = 200; blue['M'] =  0;
    
    /* yellow V, I, L */
    red['V'] = 255;  green['V'] = 255;  blue['V'] =   0;
    red['I'] = 255;  green['I'] = 230;  blue['I'] =   0;
    red['L'] = 255;  green['L'] = 175;  blue['L'] =   0;

    /* orange A */ 
    red['A'] = 255;  green['A'] = 170;   blue['A'] =  0;

    /* brigt red P */ 
    red['P'] = 255;  green['P'] = 0;   blue['P'] = 0;

    /* dull red G */ 
    red['G'] = 180;  green['G'] = 0;   blue['G'] = 0;

    /* purple: S, T */ 
    red['S'] = 255;  green['S'] = 0;   blue['S'] = 255;
    red['T'] = 212;  green['T'] = 0;   blue['T'] = 255;

    /* dark blue N,Q */
    red['N'] = 125;  green['N'] = 0;   blue['N'] = 255;
    red['Q'] = 50;  green['Q'] = 0;   blue['Q'] = 255;

    /* light blue D, E */ 
    red['D'] = 0;  green['D'] = 125;   blue['D'] = 255;
    red['E'] = 0;  green['E'] = 255;   blue['E'] = 255;

    /* blue-green K,R,H */ 
    red['K'] = 155;  green['K'] = 255;   blue['K'] = 155; 
    red['R'] =  100;  green['R'] = 255;   blue['R'] = 100; 
    red['H'] =   0;  green['H'] = 255;   blue['H'] =   0; 


    
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

    fprintf(ppmfile,"\n");
    fclose (ppmfile);
    
    return 0;
}



int read_msf (char *filename, Sequence ** alignment_ptr, int * no_seqs_ptr) {


    int       conversion_num;
    int       no_of_seqs;
    int       sequence_length, max_sequence_length;
    int       i, name_found, no_of_residues;
    char      input_line[LONGSTRING], *strptr;
    string    name, keyword1, keyword2;
    FILE *    fptr;
    Sequence *sequence;
    
    
    if( ! (fptr = fopen(filename,"r"))  ) {
	fprintf ( stderr, "Error opening %s.\n", filename);
	return 1;
    }

    no_of_seqs              = max_sequence_length = 0;
  
   
    while( fgets (input_line, LONGSTRING, fptr)  ) {

	memset ( keyword1, 0 , MAXSTRING);
	memset ( keyword2, 0 , MAXSTRING);
	memset ( name, 0 , MAXSTRING);
	
	conversion_num = sscanf(input_line,"%s%s%s%d", keyword1, name, keyword2,  &sequence_length);

	/* is the first string in the line Name: */
	
	if( conversion_num == 4 && !strcmp(keyword1,"Name:")  ) {
	    
	    if ( !strcmp(keyword2,"Len:") ) {
		
		++no_of_seqs;
		if( sequence_length > max_sequence_length) {
		    max_sequence_length = sequence_length;
		}
		
	    } else {
	    
		fprintf (stderr, "Msffile header error: expecting format \"Name: <name>  Len: <length> [anything]\".\n");
		return 1;
	    
	    } 
	} 
    }

    printf (" %d sequences of length %d \n" , no_of_seqs, max_sequence_length );

    if ( ! (sequence = (Sequence *) calloc ( no_of_seqs, sizeof (Sequence)) ) ) {
	fprintf (stderr, "Error allocating alignment space.\n");
	return 1;
    }
    for ( i=0; i< no_of_seqs; i++ ) {
	if ( ! (  sequence[i].position = (char*)  calloc ( max_sequence_length, sizeof (char) ) ) ) {
	    fprintf (stderr, "Error allocating alignment space.\n");
	    return 1;
	}
	sequence[i].length = max_sequence_length;
    }
    
    rewind (fptr);

    /* read in the names: */
    i = 0;
    while (fgets (input_line, LONGSTRING, fptr)  ) {

	conversion_num = sscanf(input_line,"%s%s%s%d", keyword1, name, keyword2,  &sequence_length);
	if( conversion_num == 4 && !strcmp(keyword1,"Name:")  ) {
	    
	    if ( !strcmp(keyword2,"Len:") ) {
	       
		sprintf ( sequence[i].name, "%s", name);
		sequence[i].length = 0;
		i ++;
		
	    } else {
	    
		fprintf (stderr, "Msffile header error: expecting format \"Name: <name>  Len: <length> [anything]\".\n");
		return 1;
	    
	    } 
	} 
	if  ( !strncmp(keyword1,"//",2) ) { /* delimits the header from the rest of the file */
		break;
	}
    }
    
    /* read in the seqeunces: */ 
    while (fgets (input_line, LONGSTRING, fptr)  ) {
	
    	if((conversion_num=sscanf(input_line,"%s",name))==1) {
	    
	    name_found = 0;
	    for ( i=0; i< no_of_seqs; i++ ) {
		if( ! (strcmp(name, sequence[i].name)) ) {
		    
		    name_found = 1;
		    strptr = strstr(input_line, name) + strlen(name);
		    while ( isspace (*strptr) ) strptr++;
		    /* keep reading until the end of the line */
		    no_of_residues = sequence[i].length;
		    while(*strptr!='\n') {
	    
			/* if it is a nonwhite space, add it to the sequence */
			if(! isspace(*strptr) ) {

			    no_of_residues++;
			    if( no_of_residues > max_sequence_length) {
				fprintf (stderr,"Sequence %d (%s) has too many residues. ",
					 i,sequence[i].name);
				fprintf (stderr," (%d residues, seq length %d, value: '%c') \n",
					 no_of_residues, max_sequence_length, *strptr);
				return 1;
			    }
			    if ((*strptr>=97)&&(*strptr<=122)) {*strptr -= 32;} /* --> turn to uppercase */
			    if (*strptr==126) {*strptr = 46;} /* turn tweedle to dot */

			   
			   
			    sequence[i].position[no_of_residues-1] = *strptr;
			    
			}
	    
			/* increment strptr */
			++ strptr;
		    }
		    sequence[i].length = no_of_residues; /* no_of_res read in so far */
		    break;
		}
		
	    } /* end for i loop */
	    if ( !name_found) {
		fprintf ( stderr, "the sequence name \"%s\" was not seen  in the name list.\n", name);
		return 1;
	    }
	}

    }

    fclose(fptr);

    

    *alignment_ptr  = sequence;
    *no_seqs_ptr = no_of_seqs;
    return 0;



}
