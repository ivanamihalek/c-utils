#include <sys/socket.h>
#include <sys/un.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

//char *socket_path = "./socket";
//char *socket_path = "\0hidden";
char *socket_path = "\0igraph_daddy";
int MAX_BUF=150;
int main(int argc, char *argv[]) {
    
    struct sockaddr_un addr;
    char buf[150];
    int fd,rc;

    if (argc > 1) socket_path=argv[1];

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

    int  size;
    int retval;
    memset (buf, 0, MAX_BUF);
    sprintf (buf, "ping: %d\n",  getpid());
    size = strlen(buf)+1;
    //send some data
    //while (
    printf (" ... sending\n");
    retval = send(fd, buf, size, 0);
    //}
    printf ("done sending\n");
    if (retval < 0) {
 	perror ("Send failed");
	exit(-1);
    }

    // wait for the answer
    memset (buf, 0, MAX_BUF);
    if ( recv (fd, buf,  MAX_BUF, 0) < 0) {
 	perror ("error reading");
	exit(-1);
    }
    printf (" %s  \n", buf);
    
    close(fd);
    return 0;
}
