# ifndef _IFT_H
# define _IFT_H


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include "pdb.h"
# include "geometry.h"
# include "utils.h"


/* contact types */

# define BB2BB 0
# define BB2SC 1
# define SC2BB 2
# define SC2SC 3

# define MAX_NBRS 30
# define MAX_FRAMES 150
# define MAX_HBONDS 400

typedef struct {
    int id;
    int count;
} Neighbor;

# endif
