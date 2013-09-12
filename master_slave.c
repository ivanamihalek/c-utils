/* compile with gcc -o one one.c -lmpich */

# include<stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <unistd.h>  /* sleep */
# include <time.h>  

#include <mpi.h>



# define MSG_LEN 20

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
	int done = 0, ctr;
	for (rank=1; rank < size; rank++) {
	    /*send message*/
	    strcpy(msg,"Hello there"); 
	    MPI_Send(msg, strlen(msg), MPI_CHAR, rank, 99, MPI_COMM_WORLD); 
	}
	/* if done, send it another job */
	ctr = 0;
	while ( ! done ) {
	    for ( rank=1; rank < size && ! done; rank++) {
		
		MPI_Iprobe(rank, 99, MPI_COMM_WORLD, &have_message, &status);
		
		if ( have_message ) {
		    MPI_Recv(msg, MSG_LEN, MPI_CHAR, rank, 99, MPI_COMM_WORLD, &status); 
		    strcpy(msg,"Hello again"); 
		    MPI_Send(msg, strlen(msg), MPI_CHAR, rank, 99, MPI_COMM_WORLD);
		    ctr ++;
		    done = ( ctr==10 );
		    printf("rank: %d   name: %s  collected %d msgs\n", myrank, name, ctr);
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
	int sleep_length;
	srand48 ( time(NULL)*myrank );
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
		sleep_length =  (int)(drand48()*10) + 1;
		printf("rank: %d   name: %s  sleeping for %ds\n", myrank, name, sleep_length);
		sleep ( sleep_length );
		/* send the message back */
		printf("rank: %d   name: %s   done\n", myrank, name, msg);
		memset   (msg, 0, MSG_LEN);
		strcpy   (msg, "done"); 
		MPI_Send (msg, strlen(msg), MPI_CHAR, 0, 99, MPI_COMM_WORLD); 
	    }
	}
    }

    /* Shut down MPI */
    MPI_Finalize();   
    return 0;

}
