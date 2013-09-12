# include "postp.h"
# include <time.h>
   /* echo options to the log file*/
int logger (Options * options, int mode, char *msg )  {
    
    char filename[LONGSTRING];
    FILE * fptr;
    sprintf ( filename, "%s.log", options->outname);

    switch ( mode) {
    case INTRO:
    {
	fptr = efopen ( filename, "w");
	time_t    time_now;
	time(&time_now);     /* get time in seconds */
	fprintf (fptr, "\n");
	fprintf (fptr, "  POSTP version %s.\n", VERSION);
	fprintf (fptr, "  by Ivana Mihalek, 2005\n");
	fprintf (fptr, "  File produced on %s", asctime(localtime(&time_now)));
	fprintf (fptr, "\n");
	fprintf (fptr, "  Options:\n");
	if (options->pdbname[0])   fprintf (fptr, "  \t pdb file:               %s\n", options->pdbname);
	if (options->chain)        fprintf (fptr, "  \t pdb chain:              %c\n", options->chain);
	if (options->almtname[0])  fprintf (fptr, "  \t alignment  file:        %s\n", options->almtname);
	if (options->query[0])     fprintf (fptr, "  \t query  name:            %s\n", options->query);
	if (options->scorename[0]) fprintf (fptr, "  \t score file:             %s\n", options->scorename);
	if (options->max_gaps)     fprintf (fptr, "  \t max fraction of gaps:   %4.2lf\n", options->max_gaps);
	switch ( options->scoring_method) {
	case ENTROPY:
	    fprintf ( fptr, "  \t residue scoring method:  entropy\n");
	    break;
	case IVET:
	    fprintf ( fptr, "  \t residue scoring method:  integer valued trace\n");
	    break;
	case RVET:
	    fprintf ( fptr, "  \t residue scoring method:  real valued trace\n");
	    break;
	}
	fclose (fptr);
    }
    break;
    case WARN:
    {
	fptr = efopen ( filename, "a");
	fprintf (fptr, "  Warning:\n");
	fprintf (fptr, "  \t %s\n", msg);
	fclose (fptr);
    }
    break;
    case STATUS:
    {
	fptr = efopen ( filename, "a");
	fprintf (fptr, "  Status:\n");
	fprintf (fptr, "  \t %s\n", msg);
	fclose (fptr);
    }
    break;
    case NOTE:
    {
	fptr = efopen ( filename, "a");
	fprintf (fptr, "  Notes:\n");
	fprintf (fptr, "  \t Remember that the results will depend on the input sequence selection,\n");
	fprintf (fptr, "  \t as well as the alignment quality and proportion of gaps in the alignment.\n");
	fclose (fptr);
    }
   
    }


    return 0;
}
 
