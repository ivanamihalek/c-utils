# ifndef _UTILS_H
# define _UTILS_H
# include <stdio.h>

int      array_qsort (int * sorted_pos, double * sa, int sequence_length );
char   **chmatrix(int rows, int columns);
double **dmatrix(int rows, int columns);
void *   emalloc(int	size);
FILE *   efopen(char * name, char * mode);
void     free_matrix(void **m);
void     free_cmatrix(char **m);
void     free_imatrix(int **m);
void     free_dmatrix(double **m);
void     free_d3matrix(double ***m);
void     free_strmatrix(char ***m);
int    **intmatrix(int rows, int columns);
char     single_letter ( char code[]);
int      string_clean ( char* string, int length);
int      tokenize ( char token[MAX_TOK][MEDSTRING],
		    int * max_token, char * line , char comment_char);
# endif
