# ifndef _UTILS_H
# define _UTILS_H

int      array_qsort (int * sorted_pos, double * sa, int sequence_length );
char   **chmatrix(int rows, int columns);
double **dmatrix(int rows, int columns);
void *   emalloc(int	size);
FILE *   efopen(char * name, char * mode);
void     free_matrix(void **m);
int    **intmatrix(int rows, int columns);
int      string_clean ( char* string, int length);
char single_letter ( char code[]);

# endif
