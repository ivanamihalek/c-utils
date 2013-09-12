/*
 * $Id: template.c,v 1.4 2001/07/23 15:28:29 lindahl Exp $
 */

/* definition of t_trxframe:
   in /usr/local/gromacs/include/gromacs/types/trx.h */

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
# include <string.h>


int main(int argc,char *argv[])
{
  static char *desc[] = {
    "this is a small test program meant to serve as a template ",
    "when writing your own analysis tools. The advantage of ",
    "using gromacs for this is that you have access to all ",
    "information in the topology, and your program will be ",
    "able to handle all types of coordinates and trajectory ",
    "files supported by gromacs. Go ahead and try it! ",
    "This test version just writes the coordinates of an ",
    "arbitrary atom to standard out for each frame. You can ",
    "select which atom you want to examine with the -n argument."
  };
  
  static int n=1;

  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
  t_pargs pa[] = {
    { "-n", FALSE, etINT, {&n},
      "Plot data for atom number n (starting on 1)"
    }
  };
  
  t_topology top;
  char       title[STRLEN];
  t_trxframe fr;
  rvec       *xtop;
  matrix     box;
  int        status;
  int        flags = TRX_READ_X;
  int grpctr, subctr;

  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  };
  
#define NFILE asize(fnm)

  //CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
   * calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  /* We don't need any topology information to write the coordinates,
   * but to show how it works we start by writing the name and
   * charge of the selected atom. It returns a boolean telling us
   * whether the topology was found and could be read
   */
  
  memset ( &top, 0, sizeof(t_topology));
  for (grpctr=0;  grpctr< egcNR; grpctr++ ) {
      top.atoms.grps[grpctr].nr = 0;
  }
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,box,TRUE);
  sfree(xtop);

  n=1; /* Our enumeration started on 1, but C starts from 0 */
  /* check that this atom exists */
  
  /* The first time we read data is a little special */
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);

  {
      int n;
      t_atoms *atoms = fr.atoms;
      for ( n=0; n<atoms->nr; n++) {
	  printf ( "\t %6d  %5s   %5s  %5s \n", n, *(atoms->atomname[n]), *(atoms->atomtype[n]),  *(atoms->resname[n]));
      }
  }
   
  /* This is the main loop over frames */
  do {
    /* coordinates are available in the vector fr.x
     * you can find this and all other structures in
     * the types directory under the gromacs include dir.
     * Note how flags determines wheter to read x/v/f!
     */
      // for ( n=0; n<fr.natoms; n++) {
//	  printf("Coordinates at t=%8.3f : %8.5f %8.5f %8.5f\n",fr.time,fr.x[n][XX],fr.x[n][YY],fr.x[n][ZZ]);
      //  }
      //printf ( "\n time=%8.3f   number of groups: %d\n", fr.time, fr.ngrpname);
      //exit (0);
  } while(read_next_frame(status,&fr));

  
  return 0;
}

