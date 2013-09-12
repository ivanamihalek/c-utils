# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>


# define XDIM  10
# define YDIM  80

int main ( int argc, char * argv[]) {
    unsigned long mask = 0, anded, carry, old_carry;
    unsigned long twoDbitmap[YDIM][XDIM] = {{0}};
    unsigned long other_twoDbitmap[YDIM][XDIM] = {{0}};
    unsigned long * array;
    int i, j, array_i, array_j, total, value;
    int x_shift, y_shift;
    time_t time_now, time_cum;
    int write_bitmap (unsigned long map[YDIM][XDIM], int i, int j, int value );
    printf ("unsigned long is %d bytes long.\n\n",
	    sizeof(mask));
    /* physical map has dimensions YDIMxXDIM*16 */
    /* writing at physical position i, j */
    /* array coordinate j/16, integer arithmetic */
    /* I want to write 01 for +1 or 10 for -1 */
    /* using 2 bits at a time */
    /* bit coordinate: (j%16)*2 [+1 for negative] */

    /* fill two maps with random numbers,
       and see how long it takes to compare the two */
    for (i=0; i<YDIM; i++) {
	for (j=0; j<XDIM*16; j++) {
	    value = ( drand48() < 0.5 ) ? 1: -1;
	    write_bitmap ( twoDbitmap, i, j, value);
	    value = ( drand48() < 0.5 ) ? 1: -1;
	    write_bitmap ( other_twoDbitmap, i, j, value);
	}
    }

    time(&time_cum);
    /* compare */
    total = 0;

    for (x_shift = 0; x_shift<XDIM*16; x_shift++ ) {

	if ( x_shift) {
	    carry = 0;
	    old_carry = 0;
	    for (i=0; i<YDIM; i++) {
		array = other_twoDbitmap [i];
		for (array_j=0; array_j<XDIM; array_j++) {
		    
		    carry = array[array_j] & 3;
		    array[array_j] >>= 2;
		    array[array_j] |= old_carry;
		    old_carry  = carry << 30;
		    
		}
	    }
	    array[0] |= old_carry;
	}
	
	for (y_shift = 0; y_shift<YDIM; y_shift++ ) {
	    
	    for (i=0; i<YDIM; i++) {
		array_i = (i+y_shift)%YDIM;
	    
		for (array_j=0; array_j<XDIM; array_j++) {
		    anded = twoDbitmap [i][array_j] &
			other_twoDbitmap [array_i][array_j];
		    if (anded) {
			mask = 1;
			do {
			    if ( anded & mask ) total++;
			    mask <<=1;
			} while (mask);
		    }
		}
	    }
	}
    }
    time (&time_now);
    printf ("Time:  %5d (s) \n",  (int)(time_now -time_cum) );
    
    return 0;
}

int write_bitmap (unsigned long map[YDIM][XDIM], int i, int j, int value ) {

    int array_j, word_j;
    unsigned long mask;

    if ( value ) {
	array_j = j/16;
	word_j = (j%16)*2;
	    
	mask = 1 << word_j;
	if ( value < 0 ) {
	    mask <<= 1;
	}
	map [i][array_j] |= mask;

    }

    return 0;
    
}
