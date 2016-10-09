#include <stdio.h>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <stdlib.h>

//char *socket_path = "./socket";
/*Linux has a special feature: if the pathname for a UNIX domain socket begins with a null byte \0, its name is not mapped into the filesystem. Thus it won’t collide with other names in the filesystem. Also, when a server closes its UNIX domain listening socket in the abstract namespace, its file is deleted; with regular UNIX domain sockets, the file persists after the server closes it.

 */

char *socket_path = "\0hidden";
int MAX_BUF=150;

int main(int argc, char *argv[]) {
    struct sockaddr_un addr;
    char buf[MAX_BUF];
    int fd,cl,rc;

    if (argc > 1) socket_path=argv[1];

    if ( (fd = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
	perror("socket error");
	exit(-1);
    }

    memset(&addr, 0, sizeof(addr));
    addr.sun_family = AF_UNIX;
    //http://troydhanson.github.io/network/Unix_domain_sockets.html
    
    if (*socket_path == '\0') {
	*addr.sun_path = '\0';
	strncpy(addr.sun_path+1, socket_path+1, sizeof(addr.sun_path)-2);
    } else {
	strncpy(addr.sun_path, socket_path, sizeof(addr.sun_path)-1);
	unlink(socket_path);
    }

    if (bind(fd, (struct sockaddr*)&addr, sizeof(addr)) == -1) {
	perror("bind error");
	exit(-1);
    }

    if (listen(fd, 5) == -1) {
	perror("listen error");
	exit(-1);
    }

    while (1) {
	if ( (cl = accept(fd, NULL, NULL)) == -1) {
	    perror("accept error");
	    continue;
	}
	memset (buf, 0, MAX_BUF);
	while ( (rc=read(cl,buf,sizeof(buf))) > 0) {
	    printf("%s\n",  buf); fflush (stdout);
	}
        printf("==================\n");  fflush (stdout);
	if (rc == -1) {
	    perror("read");
	    exit(-1);
	}
	else if (rc == 0) {
	    close(cl);
	}
    }


    return 0;
}

