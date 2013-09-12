# include "pdbclust.h"

int restype2number( char code[]){

    switch ( code[0] ) {
    case 'A':
	switch ( code [1]) {
	case 'L':
	    return 0;
	    break;
	case 'R':
	    return 15;
	    break;
	case 'S':
	    switch ( code[2] ) {
	    case 'N':
		return 12;
		break;
	    case 'P':
		return  3;
		break;
	    }
	    break;
	}
	break;
    case 'C':
	return 2;
	break;
    case 'G':
	/* the second letter is always L */ 
	switch ( code[2] ) {
	case 'U':
	    return 4;
	    break;
	case 'N':
	    return  14;
	    break;
	case 'Y':
	    return 6;
	    break;
	}
	break;
    case 'H':
	return  7;
	break;
    case 'I':
	return  8;
	break;
    case 'L':
	switch ( code [1]) {
	case 'E':
	    return 10;
	    break;
	case 'Y':
	    return 9;
	    break;
	}
	break;
    case 'M':
	return 11;
	break;
    case 'P':
	switch ( code [1]) {
	case 'H':
	    return 5;
	    break;
	case 'R':
	    return 13;
	    break;
	}
	break;
    case 'S':
	return 16;
	break;
    case 'T':
	switch ( code [1]) {
	case 'H':
	    return 17;
	    break;
	case 'R':
	    return 19;
	    break;
	case 'Y':
	    return 1;
	    break;
	}
	break;
    case 'V':
	return 18;
	break;
	
    }


    fprintf (stdout, "Unrecognized amino acid code: %s.\n", code);
    return 20;
}

