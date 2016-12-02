#include <pthread.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#define handle_error_en(en, msg)					\
    do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)

static int done = 0;
static int cnt = 0;

void
cleanup_handler(void *arg)
{
    printf("Called clean-up handler\n");
    cnt = 0;
}

void *
thread_start(void *arg)
{
    time_t start, curr;
    int oldtype; // a return value we don't really need here
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &oldtype);

    printf("New thread started\n");

    pthread_cleanup_push(cleanup_handler, NULL);

    //pthread_testcancel();
    sleep(60);
    

    pthread_cleanup_pop(1);
    return NULL;
}

int
main(int argc, char *argv[])
{
    pthread_t thr;
    int s;
    void *res;

    s = pthread_create(&thr, NULL, thread_start, NULL);
    pthread_detach(thr);
    if (s != 0)
	handle_error_en(s, "pthread_create");
   

    sleep(2);           /* Allow new thread to run a while */

    printf("Canceling thread\n");
    if (pthread_cancel(thr)) handle_error_en(s, "pthread_cancel");
    
    sleep(2);           /* make sure the thread is not killed by the parent process exit */
    printf ("Main exiting.\n");
    return 0;
}

