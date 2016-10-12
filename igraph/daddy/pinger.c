#include <sys/socket.h>
#include <sys/un.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

char *socket_path = "\0igraph_daddy";
int MAX_BUF=1024;

void * emalloc(int size) {
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %d bytes", size);
	exit(1);
    }
   return ptr;
}

/////////////////////////////////////
int main(int argc, char *argv[]) {

    if (argc < 2) {
	fprintf (stderr, "Usage: %s \"igraph command line\".\n", argv[0]);
	exit(1);
    }
    
#ifdef LONG_INPUT
    #ifdef TOO
    MAX_BUF *= 50;
    # else
    MAX_BUF *= 5;
    # endif
#endif
    
    struct sockaddr_un addr;
    char * buf, * aux_buf;
    int fd,rc;

    buf = emalloc(MAX_BUF);
    
    if ( (fd = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
	perror("socket error");
	exit(-1);
    }

    memset(&addr, 0, sizeof(addr));
    addr.sun_family = AF_UNIX;
    if (*socket_path == '\0') {
	*addr.sun_path = '\0';
	strncpy(addr.sun_path+1, socket_path+1, sizeof(addr.sun_path)-2);
    } else {
	strncpy(addr.sun_path, socket_path, sizeof(addr.sun_path)-1);
    }

    if (connect(fd, (struct sockaddr*)&addr, sizeof(addr)) == -1) {
	perror("connect error");
	exit(-1);
    }

#ifdef NWLN
    sprintf (buf, "%s\n", argv[1]);
#else
    sprintf (buf, "%s", argv[1]);
#endif
    
#ifdef LONG_INPUT
    while (strlen(buf)+strlen(" 12345 6789 ")<MAX_BUF) sprintf(buf+strlen(buf),  " 12345 6789 ");
#endif

    int size = strlen(buf)+1;
    int retval;

    //send some data
# ifndef LONG_INPUT
    printf (" ... sending: %s\n", buf);
# else
    aux_buf = emalloc(100);
    memcpy (aux_buf, buf, 30);
    printf (" ... sending big input (%d bytes)  %30s ...\n", size, aux_buf);
    free(aux_buf);
# endif
    retval = send(fd, buf, size, 0);
    //}
    printf ("done sending\n");
    if (retval < 0) {
 	perror ("Send failed");
	exit(-1);
    }

    // wait for the answer
#ifndef SKIP_RECV
    memset (buf, 0, MAX_BUF);
    if ( recv (fd, buf,  MAX_BUF, 0) < 0) {
 	perror ("error reading");
	exit(-1);
    }
    #ifdef LONG_INPUT
    sprintf (buf+100, " ...");
    buf[strlen(buf)] = '\0';
    #endif    
    printf ("%s  \n", buf);
#endif    
    close(fd);

    free(buf);
    return 0;
}
