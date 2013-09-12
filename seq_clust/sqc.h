# ifndef SQC_H
# define SQC_H

#define BUFLEN  150
#define NAMELEN 50
typedef struct {
    char name[NAMELEN];
    char *position;
} Sequence;


int read_msf ( char msfname[], Sequence ** seqptr, int * noseqptr, int * seqlenptr);
int find_dist_matrix(Sequence *sequence, int no_seqs, int seq_len, double ** simmat, int sim);
void cluster_counter (int  no_of_things,  int *neighbors[], int * mask,
		     int cluster_count_per_size[], int * no_of_clusters,
		     int * max_size, int * secnd_max_size , int * clusters[]);
	



#endif
