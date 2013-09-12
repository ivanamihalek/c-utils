# include <ctype.h>

int charcount (char * word, int len) {
    int i,ct;
    ct = 0;
    for (i=0; i<len; i++) {
	ct += isalnum ( word[i]) ? 1: 0;
    }
    return ct;
}
