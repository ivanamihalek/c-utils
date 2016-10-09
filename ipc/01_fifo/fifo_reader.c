// you can probe it by running it, in bg perhaps,
// and then pinging it by echo "something" > /tmp/myfifo

#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h> // access declared here

#define MAX_BUF 1024

int main()
{
    int fd;
    char * myfifo = "/tmp/myfifo";
    char buf[MAX_BUF];

    /* check file exists and readable */
    if( access( myfifo, R_OK ) == -1 ) {
	fprintf (stderr, "%s does not exist or is not readable.\n", myfifo);
	return 1;
    }     

    /* open, read, and display the message from the FIFO */
     while (1) {
        fd = open(myfifo, O_RDONLY);
	memset (buf, 0, MAX_BUF);
        read (fd, buf, MAX_BUF);
        printf("Received: %s", buf);// it looks like \n travels along
        close(fd);
    }
    return 0;
}
