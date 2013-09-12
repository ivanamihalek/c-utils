/* compile with gcc -o one one.c -lmpich */

# include<stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <unistd.h>  /* sleep */
# include <time.h>  
# include <sys/types.h>
# include <sys/wait.h>


#include <mpi.h>



# define MSG_LEN 20
char * pdb_name[] = {"1a59", "1ad3A", "1ai2", "1aln", "1an9A", "1btoA", "1cg0A", "1digB", "1dqrA", "1dqxA", "1e2dA", "1e5qH", "1e7yA", "1e98A", "1ek4A", "1h16A", "1h7tB", "1hkvA", "1j79A", "1kc3A", "1kerB", "1l5wA", "1l9wA", "1lbxA", "1lxyA", "1m0sA", "1m4nA", "1m7pA", "1m9nA", "1nywB", "1qinA", "1tb5A", "1tc2A", "1umpC", "1w1uA", "1y6rB", "2bifA", "2bwoD", "2dorA", "6gstA"};
int nr_of_names = 40;

int main ( int argc, char * argv[]) {

    int myrank, size;
    int rank;
    char  name[MPI_MAX_PROCESSOR_NAME];
    char msg[MSG_LEN]; 
    int  length;
    MPI_Status status;
    
    /* Initialize MPI */
    MPI_Init (&argc, &argv);

    /* Find out the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* Find out my identity in the default communicator */
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Get_processor_name(name, &length);

    /******************************************************************/
    /*         MASTER                                                 */
    /******************************************************************/
    if ( myrank == 0) { /*master */
	int have_message = 0;
	int done = 0, ctr = 0;
	
	/* keep sending jobs, untill done */
	while ( ! done ) {
	    for ( rank=1; rank < size && ! done; rank++) {
		
		MPI_Iprobe(rank, 99, MPI_COMM_WORLD, &have_message, &status);
		
		if ( have_message ) {
		    MPI_Recv(msg, MSG_LEN, MPI_CHAR, rank, 99, MPI_COMM_WORLD, &status);

		    strcpy(msg, pdb_name[ctr]); 
		    MPI_Send(msg, strlen(msg), MPI_CHAR, rank, 99, MPI_COMM_WORLD);
		    ctr ++;
		    done = ( ctr== nr_of_names);
		    printf("rank: %d   name: %s  collected %d msgs\n", myrank, name, ctr);
		    fflush(stdout);
		    
		}
	    }
	}
        /* when  done, inform all the processes about it */
	for (rank=1; rank < size; rank++) {
	    /*send message*/
	    strcpy (msg,"done"); 
	    MPI_Send(msg, strlen(msg), MPI_CHAR, rank, 99, MPI_COMM_WORLD); 
	}
	
    /******************************************************************/
    /*         SLAVE                                                  */
    /******************************************************************/
    } else {
	int done = 0;
	srand48 ( time(NULL)*myrank );
	/* say hi to master */ 
	memset   (msg, 0, MSG_LEN);
	strcpy   (msg, "I'm up"); 
	MPI_Send (msg, strlen(msg), MPI_CHAR, 0, 99, MPI_COMM_WORLD); 
	/* loop until  done */
	while ( ! done ) {
	    /* receive a mesage */
	    memset (msg, 0, MSG_LEN);
	    MPI_Recv(msg, MSG_LEN, MPI_CHAR, 0, 99, MPI_COMM_WORLD, &status); 
	    printf("rank: %d   name: %s   received: %s\n", myrank, name, msg);
	    fflush(stdout);
	    if ( !strcmp(msg, "done") ){
		done = 1;
	    } else {
		/* do the job */
		
		pid_t pid;
		int retval;
		char cmd[250];
		char local_name[50];

		memset  (local_name, 0, 50);
		sprintf (local_name, "%s", msg);
		memset  (msg, 0, MSG_LEN);
		
		switch(pid=fork()) {
		case -1: /* failure */
		    strcpy   (msg, "fail fork");
		    break;

		case 0:  /* child process */
		    memset (cmd, 0, 250);
		    sprintf ( cmd,
			      "/home/imihalek/projects/rate4site/src_18_10_04/rate4site -s %s.hssp.pruned.98_15.ph -a %s  -o %s.r4s ",
			      local_name, local_name,  local_name);
		    retval = system (cmd);
		    exit(retval);

		default: /*parent */
		    wait(&retval);
		    if ( retval ) {
			strcpy   (msg, "fail r4s");
		    } else {
			strcpy   (msg, "success");
		    }
		    
		}

		
		/* send the message back */
		printf("rank: %d   name: %s   %s\n", myrank, name, msg);
		fflush(stdout);
		MPI_Send (msg, strlen(msg), MPI_CHAR, 0, 99, MPI_COMM_WORLD); 
	    }
	}
    }

    /* Shut down MPI */
    MPI_Finalize();   
    return 0;

}
