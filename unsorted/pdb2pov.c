/****************************************************************************
 * pdb2pov.c
 * Author Doug Wagner
 * Copyright 1989, Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Permission use, copy, and modify this software and its documentation
 * without fee for personal use or use within your organization is hereby
 * granted, provided that the above copyright notice is preserved in all
 * copies and that that copyright and this permission notice appear in
 * supporting documentation.  Permission to redistribute this software to
 * other organizations or individuals is not granted;  that must be
 * negotiated with the PSC.  Neither the PSC nor Carnegie Mellon
 * University make any representations about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *
 * Modified by Andrew Haveland-Robinson to work for PoV
 *****************************************************************************/

#include <stdio.h>
#include <math.h>

/* Define file that holds size, color and shading info. */
#ifndef ELEMENT_TABLE
#define ELEMENT_TABLE "element.tab"
#endif

#define n_elements 120
#define mxcoord 100.0
#define maxparticles 5000

/*
Object boundaries - these will get set to proper values as atoms are parsed.
*/
static float xmax= -mxcoord, xmin= mxcoord;
static float ymax= -mxcoord, ymin= mxcoord;
static float zmax= -mxcoord, zmin= mxcoord;

/* Count of atoms */
static int natoms= 0;

/*
struct to hold atom specifications.
*/
struct atom {
	int 	serial_num;	/* serial number */
	char 	name[5];	/* name */
	char 	ali;		/* alternate location indicator */
	char 	res[4];		/* residue name */
	char 	chain_id;	/* chain identifier */
	int	seq_num;	/* residue sequence number */
	char	res_i_code;	/* residue insertion code */
	float	x,y,z;		/* x, y, z coordinate */
	float	occupancy;	/* occupancy */
	float	temp_fac;	/* temperature factor */
	int	footnote;	/* footnote number */
	int	atomic_num;	/* atomic number */
	float	r,g,b;		/* red, green, blue intensities */
	float	radius;		/* atomic radius */
	char	shade[10];	/* shader */
};

/*
Table to facilitate translation from chemical element names to
atomic numbers.
*/
struct { char name[2]; int num; } name_list[]= {
	{" H",1},{" C",6},{" N",7},{" O",8},{" S",16},{"CU",29},
      	{"HE",2},{"LI",3},{"BE",4},{" B",5},{" F",9},{"NE",10},
	{"NA",11},{"MG",12},{"AL",13},{"SI",14},{" P",15},{"CL",17},
	{"AR",18},{" K",19},{"CA",20},{"SC",21},{"TI",22},{" V",23},
	{"CR",24},{"MN",25},{"FE",26},{"CO",27},{"NI",28},{"ZN",30},
	{"GA",31},{"CE",32},{"AS",33},{"SE",34},{"BR",35},{"KR",36},
	{"RB",37},{"SR",38},{" Y",39},{"ZR",40},{"NB",41},{"MO",42},
	{"TC",43},{"RU",44},{"RH",45},{"PD",46},{"AG",47},{"CD",48},
	{"IN",49},{"SN",50},{"SB",51},{"TE",52},{" I",53},{"XE",54},
	{"CS",55},{"BA",56},{"LA",57},{"CE",58},{"PR",59},{"ND",60},
	{"PM",61},{"SM",62},{"EU",63},{"GD",64},{"TB",65},{"DY",66},
	{"HO",67},{"ER",68},{"TM",69},{"YB",70},{"LU",71},{"HF",72},
	{"TA",73},{" W",74},{"RE",75},{"OS",76},{"IR",77},{"PT",78},
	{"AU",79},{"HG",80},{"TL",81},{"PB",82},{"BI",83},{"PO",84},
	{"AT",85},{"RN",86},{"FR",87},{"RA",88},{"AC",89},{"TH",90},
	{"PA",91},{" U",92},{"NP",93},{"PU",94},{"AM",95},{"CM",96},
	{"BK",97},{"CF",98},{"ES",99},{"FM",100},{"MD",101},{"NO",102},
	{"LR",103},{"RF",104},{"HA",105},
	{"\0\0",0},
};

/*
Table to hold the table of element attributes by atomic number.
*/
struct { float red,green,blue,radius;
	 char shade[10];
	 int is_defined; } el_table[n_elements];

/*
Input buffer definition follows.
*/
#define maxrec 256
static char inbuf[maxrec];

main(argc,argv)
int argc;
char **argv;
{
	int ierr= 0;
	FILE *pdbfile, *output, *fsetup();
	char dummy[128], scn[128];

	if (argc != 2 && argc != 3)
		{
		fprintf(stderr,"%s: usage: %s pdbfilename [elementfile]\n",
			argv[0],argv[0]);
		exit(0);
		};

	if (argc > 3)
		load_el_table(argv[2]);
	else
		load_el_table(ELEMENT_TABLE);

	/* Get the name of the output file and open it */
	strcpy(scn,argv[1]);
	change_end(scn);
	output= fsetup(scn,"w",&ierr);
	if (ierr) exit(2);

	/* Write the heading of the program */
	preamble(output);

	/* Open PDBfile and scan for elements, close PDBfile */
	pdbfile= fsetup(argv[1],"r",&ierr);
	if (ierr) exit(2);
	el_scan(pdbfile,output);
	fclose(pdbfile);

	/* Open PDBfile again, and finish writing the program */
	pdbfile= fsetup(argv[1],"r",&ierr);
	if (ierr) exit(2);
	intro(output);
	make_molecule(pdbfile,output);
	fclose(pdbfile);
	finale(output);
	fclose(output);
}

make_molecule(pdbfile,output)
FILE *pdbfile, *output;
{
	while (!feof(pdbfile)) do_rec(pdbfile,output);
	fprintf(stdout,"%d atoms.\n",natoms);
}

el_scan(pdbfile,output)
FILE *pdbfile, *output;
{
	while (!feof(pdbfile)) do_el(pdbfile,output);
}

load_el_table(fname)
char fname[];
{
	FILE *table, *fsetup();
	int ierr= 0,i;
	float r,g,b,radius;
	char string[10];

	table= fsetup(fname,"r",&ierr);
	if (ierr) exit(2);

	if ( !fgets(inbuf,maxrec,table) )
		{
		fprintf(stderr,"Error reading element table input record.\n");
		exit(2);
		};

	while (!feof(table) && fgets(inbuf,maxrec,table))
		{
		sscanf(inbuf,"%d %f %f %f %f %s",&i,&r,&g,&b,&radius,string);
		el_table[i-1].red= r;
		el_table[i-1].green= g;
		el_table[i-1].blue= b;
		el_table[i-1].radius= radius;
		trim_quotes(string);
		strcpy(el_table[i-1].shade,string);
		el_table[i-1].is_defined= 0;
		};

	if (ferror(table))
		fprintf(stderr,
			"warning: error reading element table data.\n");

	fclose(table);
}

do_rec(pdbfile,output)
FILE *pdbfile, *output;
{
	struct atom thisatom;
	char typestring[7];
	float get_float();
	int get_int(), i;
	char name[2];
	float x, y, z;

	fgets(inbuf,maxrec,pdbfile);
	if ( ferror(pdbfile) )
		{
		fprintf(stderr,"Error reading PDB input record.\n");
		exit(2);
		};
	/* Make sure it's an 'ATOM' entry */
	strncpy(typestring,&inbuf[0],6);
	typestring[7]= '\0';
	if ( strncmp(typestring,"ATOM  ",6) ) return(1);
// XXX HETATM 	if ( strncmp(typestring,"ATOM  ",6) ) return(1);

	/* Increment the atom counter */
	natoms++;

	/* Extract and store the data in the ATOM record */
	thisatom.serial_num= get_int(inbuf+6,5);
	strncpy(thisatom.name,&inbuf[12],4); thisatom.name[4]= '\0';
	thisatom.ali= inbuf[16];
	strncpy(thisatom.res,&inbuf[17],3); thisatom.res[3]= '\0';
	thisatom.chain_id= inbuf[21];
	thisatom.seq_num= get_int(inbuf+22,4);
	thisatom.res_i_code= inbuf[26];
	thisatom.x= get_float(inbuf+30,8);
	thisatom.y= get_float(inbuf+38,8);
	thisatom.z= get_float(inbuf+47,8);
	thisatom.occupancy= get_float(inbuf+54,6);
	thisatom.temp_fac= get_float(inbuf+60,6);
	thisatom.footnote= get_int(inbuf+67,3);

	/* Look up the element to see if it is is_defined */
	if (!el_lookup(&thisatom,0)) return(1);

	/* Push out the object boundaries if necessary */
	if ( ((float) fabs((double)thisatom.x )) > mxcoord )
		{
		fprintf(stderr,
		  "atom %d (%s) has x coordinate %f: out of bounds\n",
		  thisatom.serial_num,thisatom.name,thisatom.x);
		return(1);
		};
	if ( ((float) fabs((double)thisatom.y)) > mxcoord)
		{
		fprintf(stderr,
		  "atom %d (%s) has y coordinate %f: out of bounds\n",
		  thisatom.serial_num,thisatom.name,thisatom.y);
		return(1);
		};
	if ( ((float) fabs((double)thisatom.z)) > mxcoord)
		{
		fprintf(stderr,
		  "atom %d (%s) has z coordinate %f: out of bounds\n",
		  thisatom.serial_num,thisatom.name,thisatom.z);
		return(1);
		};
	if (thisatom.x > xmax) xmax= thisatom.x;
	if (thisatom.x < xmin) xmin= thisatom.x;
	if (thisatom.y > ymax) ymax= thisatom.y;
	if (thisatom.y < ymin) ymin= thisatom.y;
	if (thisatom.z > zmax) zmax= thisatom.z;
	if (thisatom.z < zmin) zmin= thisatom.z;

	/* Add the atom to the model */
	x = thisatom.x;
	y = thisatom.y;
	z = thisatom.z;
	strncpy(name,thisatom.name,2);
	if (name[0]==' ')
		fprintf (output,"object { ATOM_%c translate <%g %g %g>}\n",name[1],x,y,z);
	else
		fprintf (output,"object { ATOM_%c%c translate <%g %g %g>}\n",name[0],name[1],x,y,z);
}

do_el(pdbfile,output)
FILE *pdbfile, *output;
{
	char typestring[7];
	struct atom thisatom;

	fgets(inbuf,maxrec,pdbfile);
	if ( ferror(pdbfile) )
		{
		fprintf(stderr,"Error reading PDB input record.\n");
		exit(2);
		};
	/* Make sure it's an 'ATOM' entry */
	strncpy(typestring,&inbuf[0],6);
	typestring[7]= '\0';
	if ( strncmp(typestring,"ATOM  ",6) ) return(1);

	/* Lookup the atom and extract info from el_table */
	strncpy(thisatom.name,&inbuf[12],4); thisatom.name[4]= '\0';
	if (!el_lookup(&thisatom,1)) return(1);

	/* Define a function for that atom. */
	def_atom(&thisatom,output);
}

el_lookup(thisatom,flag)
struct atom *thisatom;
int flag;
{
	int i= 0,atomic_num;

 	while ( strncmp(thisatom->name,name_list[i].name,2)
		&& name_list[i].num ) i++;

	if (name_list[i].num) /* The named element was found */
		{
		if (!flag) return(1);
		atomic_num= name_list[i].num;
		if (el_table[atomic_num-1].is_defined) return(0);
		el_table[atomic_num-1].is_defined = 1;
		thisatom->r= el_table[atomic_num-1].red;
		thisatom->g= el_table[atomic_num-1].green;
		thisatom->b= el_table[atomic_num-1].blue;
		thisatom->radius= el_table[atomic_num-1].radius;
		strcpy(thisatom->shade,el_table[atomic_num-1].shade);
		}
	else /* Named element not found */
		{
		if (flag)
		 	fprintf(stderr,"Couldn't find element entry for %s\n",
				thisatom->name);
		return(0);
		}

	return(1);
}

preamble(output)
FILE *output;
{
fprintf(output,"/* Converted from PDB */\n\
#include \"colors.inc\"\n\
#include \"textures.inc\"\n\n");
}

def_atom(thisatom,output)
struct atom *thisatom;
FILE *output;
{
	char name[2];
	float r,g,b,d;

	strncpy(name,thisatom->name,2);
	r = thisatom->r;
	g = thisatom->g;
	b = thisatom->b;
	d = thisatom->radius;

	/* Write the definition of the atom to the SCN file */
	if (name[0]==' ')
		fprintf(output,"#declare ATOM_%c =\nobject {\n",name[1]);
	else
		fprintf(output,"#declare ATOM_%c%c =\nobject {\n ",
		        name[0],name[1]);
	fprintf(output,"   sphere {<0.0 0.0 0.0> %g\n",d);
	fprintf(output,"          texture { %s\n",thisatom->shade);
	fprintf(output,"                  colour red %g green %g blue %g\n",r,g,b);
	fprintf(output,"                  }\n");
	fprintf(output,"          }\n");
	fprintf(output,"}\n\n");
}

intro(output)
FILE *output;
{
	fprintf(output,"// Atoms\n");
}

finale(output)
FILE *output;

/* This routine calculates the camera position information, and writes
   the camera command.
*/

{
#define tan_fov_angle 0.46630766 /* tangent of roughly half the fovea angle */
	float cam_x, cam_y, cam_z, look_x, look_y, look_z, hither, yon;
	float light_x, light_y, light_z;

	look_x= (xmax+xmin)/2.0;
	look_y= (ymax+ymin)/2.0;
	look_z= (zmax+zmin)/2.0;

	cam_x= look_x;
	cam_y= ymax + (5/tan_fov_angle)*
		( (xmax-xmin > zmax-zmin) ? (xmax-xmin) : (zmax-zmin) );
/* was 0.5/ */

	cam_z= look_z;


	hither= 0.25*( ymax - cam_y );
	yon= 2.0*( ymin - cam_y );

	light_x= -cam_x;
	light_y= (cam_y+look_y)/2.0;
	light_z= cam_z;

	fprintf(output,"camera {\n");
fprintf(output," location <0 0 0>\n");
fprintf(output," translate <%g %g %g>\n",cam_x, cam_y, cam_z);
fprintf(output," direction <0 10 0>\n");
fprintf(output," look_at <%g %g %g>\n",look_x, look_y, look_z);
fprintf(output," up <0 0 1>\n");
fprintf(output," sky <0 0 1>\n");
fprintf(output," right <1.33 0 0>\n");
fprintf(output," // hither %g\n",hither);
fprintf(output," // yon %g\n",yon);
fprintf(output,"}\n\n");
fprintf(output,"object {light_source { <%g %g %g> colour White}}\n",light_x, light_y, light_z);
fprintf(output,"object {light_source { <%g %g %g> colour White}}\n",-light_x, light_y, light_z);
}

/*
camera {
 location <0 0 0>
 translate <0 160 18>
 direction <0 1.6 0>
 right <1.33 0 0>
 look_at <0 0 18>
 // hither -4.67502
 // yon -131.26
}

camera {
 location <0 0 0>
 translate <18.575 0 60.7401>
 direction <0 0 1.6>
 look_at <0 0 18.575>
 right <1.33 0 0>
 // hither -4.67502
 // yon -131.26
}


*/




trim_quotes(string)
char *string;
{
	char *next;

	next= string+1;
	while ( (*next) && (*next != '\'') ) *string++ = *next++;
	*string= '\0';
}

change_end(s)
char s[];
{
	int i;

	i=strlen(s)-1;
	while(s[i] != '.') i--;
	s[++i]='p';
	s[++i]='o';
	s[++i]='v';
	s[++i]='\0';
}

float get_float(string,length)
char *string;
int length;
{
	char buffer[maxrec];
	float val;

	strncpy(buffer,string,length);
	buffer[length]= '\0';
	sscanf(buffer,"%f",&val);

	return(val);
}

int get_int(string,length)
char *string;
int length;
{
	char buffer[maxrec];
	int val;

	strncpy(buffer,string,length);
	buffer[length]= '\0';
	sscanf(buffer,"%d",&val);

	return(val);
}

FILE *fsetup(fname,fmode,ierr)
char fname[],fmode[];
int *ierr;
/*    fsetup safely opens a file and returns its pointer. */
{
     FILE *fp;
     if (fmode[0]=='w')
	    {
            if ( strcmp(fname,"-") == 0 ) fp= stdout;
            else fp= fopen(fname,fmode);
	    if ( fp == NULL )
	          { fprintf(stderr,"Error opening %s\n",fname); *ierr= 1;};
	    }
      else
            {
            if ( (fp= fopen(fname,fmode)) == NULL )
	          { fprintf(stderr,"Error opening %s\n",fname); *ierr= 1;};
            };
	return(fp);
}
