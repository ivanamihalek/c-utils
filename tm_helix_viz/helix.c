/* a simple way to produce right handed protein-like alpha helix
   in  z-direction - not sure if it is the cheapest way to do things if
   the direction is arbitrary */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>


# define BUFSIZE  150
# define NO_STEPS 20
# define STEP_SIZE 3.8  /*typcal distance btw c-alpha*/
# define CLIMB_THETA  acos(1.5/STEP_SIZE)   /* helix climb per C-alpha in units of STEP_SIZE*/
# define TURN (2*M_PI/3.6)  /* helix turn per c-alpha */

int helix ( FILE *fptr, int no_steps, double theta,  double* x_ptr, double* y_ptr, double* z_ptr);
int loop ( FILE *fptr, int l, int direction, double* x_ptr, double* y_ptr, double* z_ptr);
int pdbout ( FILE *, double, double, double);

int main ( int argc, char * argv[]) {
    int step, no_steps= NO_STEPS, l;
    int ctr, direction;
    char buf[BUFSIZE];
    double x, y, z;
    double d;
    double theta;
    FILE *fptr = stdout, *infptr;

    int chain_length, no_helices, * helix_start , *helix_end;
    int helix_length, shortest_helix;

    if ( argc!=2 ) {
	printf ("Usage: %s <input_file>.\n", argv[0]);
	exit (1);
    }
    
    infptr = fopen ( argv[1], "r" );
    if ( ! infptr) {
	fprintf ( stderr,"Cno argv[1].\n");
	exit (1);
    }

    fgets ( buf, BUFSIZE, infptr);    
    sscanf ( buf," %d",  & chain_length);
    
    fgets ( buf, BUFSIZE, infptr);    
    sscanf ( buf," %d",  & no_helices);
    
    helix_start = (int *) calloc ( no_helices, sizeof(int));
    helix_end = (int *) calloc ( no_helices, sizeof(int));
    if ( ! helix_start || ! helix_end ) {
	fprintf (stderr, "memalloc error.\n");
	exit (1);
    }
     
    for ( ctr=0; ctr < no_helices; ctr++ ) {
	fgets ( buf, BUFSIZE, infptr);    
	sscanf ( buf," %d%d", &helix_start[ctr],  &helix_end[ctr]); 
    }
    
    /* find the shortest helix */
    shortest_helix =  chain_length;
    for ( ctr=0; ctr < no_helices; ctr++ ) {
	helix_length =  helix_end[ctr] -  helix_start[ctr] + 1;
	if ( shortest_helix > 	helix_length ) {
	    shortest_helix =  helix_length;
	}
    }
    /* membrane thickness  - shortest helix*/
    d = shortest_helix;

    /* do something about initial loop, if exists */
    if ( helix_start[0] > 1 ) {
	l =  helix_start[0] -1 ;
	if ( l < 3 ) {
	    fprintf (stderr, "Fix me: cannot handle init loop shorter than 3.\n");
	    exit (1); 
	}
	direction = -1; 
	x = -2; 
	y = z = 0;
	loop (fptr, l, direction,  &x, &y, &z);
    }

    direction = 1;
    x = y = z = 0;
    for ( ctr=0; ctr < no_helices; ctr++ ) {

	
	/****  helix *************/
	/* number of steps through the membrane */
	no_steps =  helix_end[ctr] -  helix_start[ctr] + 1;
	/* tilt angle: */
	theta =  acos (d/no_steps);
	if ( direction < 0 ) theta = M_PI - theta;
	helix (fptr, no_steps, theta,  &x, &y, &z);

	
	/****** loop *****/
	/* loop length */
	if ( ctr < no_helices-1 ) {
	    l = helix_start[ctr+1] - helix_end[ctr] - 1;
	} else {
	    l = chain_length -  helix_end[ctr];
	}
	loop (fptr, l, direction,  &x, &y, &z);
	
	
	/* toggle climb direction */
	direction = - direction;
   }
    return 0;
}

/**************************************************************************************/
/**************************************************************************************/
int loop ( FILE *fptr, int l, int direction, double* x_ptr, double* y_ptr, double* z_ptr) {

    double x, y, z;
    int ctr;
    
    x = *x_ptr;
    y = *y_ptr;
    z = *z_ptr;

    
    for (ctr=1; ctr<=(l-1)/2; ctr++) {
	y += 1*direction;
	pdbout (fptr,  x, y, z );
    }
    if ( l%2 )  { /* the loop length is odd*/
	x+=1.0;
	pdbout (fptr,  x, y, z );
	x+=1.0;
    } else { /* the loop length is even make a small "dome" at the top*/
	x+=0.5;
	y+=sqrt(3.0)/2*direction;
	pdbout (fptr,  x, y, z );
	x+=1.0;
	pdbout (fptr,  x, y, z );
	x+=0.5;
	y-=sqrt(3.0)/2*direction;
    }
    for (ctr=(l-1)/2; ctr >=1; ctr--) {
	pdbout (fptr,  x, y, z );
	y -= 1*direction;
    }

    
    // if ( cytoplasmic ) y =-y;
    //pdbout (fptr, step, x, y, z );
    *x_ptr = x;
    *y_ptr = y;
    *z_ptr = z;
    
    return 0;
}
/**************************************************************************************/
/**************************************************************************************/
int helix ( FILE *fptr, int no_steps, double theta, double* x_ptr, double* y_ptr, double* z_ptr) {
    double x, y, z;
    double x_rot, y_rot, z_rot;
    double x_transl = 0, y_transl = 0, z_transl = 0l;
    double phi;
    double cos_theta, sin_theta;
    double projection, climb;
    int step;
    
    /* helix constants */
    climb = cos (CLIMB_THETA);
    projection = sqrt ( (1-climb*climb)/(2*(1-cos(TURN))));

    /* starting position */
    x = 0;
    y = 0;
    z = 0;
    

    phi = 0;
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for (step=1; step<=no_steps; step++) {
	phi += TURN;
	y   += climb;
	z    = projection*cos(phi);
	x    = projection*sin(phi);

	
	//printf ( " %8.3f   %8.3f   %8.3f \n", x, y, z);

	/* rotate by an angle theta from the z direction (in x-y plane) */
	z_rot = z;
	y_rot = y*cos_theta - x* sin_theta;
	x_rot = y*sin_theta + x* cos_theta;

	if ( step==1 ) {
	    x_transl = *x_ptr - x_rot;
	    y_transl = *y_ptr - y_rot;
	    z_transl = *z_ptr - z_rot;
	}
	/* translate to the tip of the last helix */
	x_rot += x_transl;
	y_rot += y_transl; 
	z_rot += z_transl;
	
	
	pdbout (fptr, x_rot, y_rot, z_rot );
    }
   
    *x_ptr = x_rot;
    *y_ptr = y_rot;
    *z_ptr = z_rot;
    return 0;
}

/**************************************************************************************/
/**************************************************************************************/
int pdbout ( FILE *fptr,  double x, double y, double z) {
    int static ctr = 0;
    ctr++;
    if ( !fptr) {
	fprintf ( stderr, "Error pdbout needs valid file ptr.\n");
	return 1;
    }
    x*=STEP_SIZE;
    y*=STEP_SIZE;
    z*=STEP_SIZE;
    fprintf ( fptr, "ATOM      2  CA  ALA  %4d    %8.3lf%8.3lf%8.3lf\n", ctr, x, y, z);
    return 0;
}
