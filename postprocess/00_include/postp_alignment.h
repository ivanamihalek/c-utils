# define ALMT_NAME_LENGTH 30

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int * seq_gaps;
    int * column_gaps;
    int total_gaps;
    double ** seq_dist;
}  Alignment;

int count_gaps (Alignment * alignment);
