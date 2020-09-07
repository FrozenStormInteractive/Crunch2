#include "crn_core.h"

#include <signal.h>
#include <sys/time.h>

#include "crn_darwin_pthreads.h"

/* A structure of type timeoutDetails is passed to the thread used to implement the timeout. */
typedef struct
{
    /* Specifies the delay, relative to now */
    struct timespec delay;
    /* The thread doing the sem_wait call */
    pthread_t callingThread;
    /* Address of a flag set to indicate that the timeout was triggered. */
    volatile short* timedOutShort;
} timeoutDetails;

/*
 * A structure of type cleanupDetails is passed to the thread cleanup
 * routine which is called at the end of the routine or if the thread calling
 * it is cancelled.
 */
typedef struct
{
    /* Address of the variable that holds the Id of the timeout thread. */
    pthread_t* threadIdAddr;
    /* Address of the old signal action handler. */
    struct sigaction* sigHandlerAddr;
    /* Address of a flag set to indicate that the timeout was triggered. */
    volatile short* timedOutShort;
} cleanupDetails;

static void* timeoutThreadMain(void* passedPtr);
static void timeoutThreadCleanup(void* passedPtr);
static int triggerSignal(int Signal, pthread_t Thread);
static void ignoreSignal(int Signal);

/*
 * sem_timedwait() implementation for Apple.
 */
int sem_timedwait(sem_t* sem, const struct timespec* abs_timeout)
{
    /* Code returned by this routine 0 or -1 */
    int result = 0;

    /*  "Under no circumstances shall the function fail if the semaphore
     *  can be locked immediately". So we try to get it quickly to see if we
     *  can avoid all the timeout overheads.
     */

    if (sem_trywait(sem) == 0)
    {
        result = 0;
    }
    else
    {
        /*  No, we've got to do it with a sem_wait() call and a thread to run
         *  the timeout. First, work out the time from now to the specified
         *  timeout, which we will pass to the timeout thread in a way that can
         *  be used to pass to nanosleep(). So we need this in seconds and
         *  nanoseconds. Along the way, we check for an invalid passed time,
         *  and for one that's already expired.
         */
        if ((abs_timeout->tv_nsec < 0) || (abs_timeout->tv_nsec > 1000000000))
        {
            /* Passed time is invalid */
            result = -1;
            errno = EINVAL;
        }
        else
        {

            struct timeval currentTime;                              /* Time now */
            long secsToWait, nsecsToWait;            /* Seconds and nsec to delay */
            gettimeofday(&currentTime, nullptr);
            secsToWait = abs_timeout->tv_sec - currentTime.tv_sec;
            nsecsToWait = (abs_timeout->tv_nsec - (currentTime.tv_usec * 1000));
            while (nsecsToWait < 0)
            {
                nsecsToWait += 1000000000;
                secsToWait--;
            }
            if ((secsToWait < 0) || ((secsToWait == 0) && (nsecsToWait < 0)))
            {

                /*  Time has passed. Report an immediate timeout. */

                result = -1;
                errno = ETIMEDOUT;

            }
            else
            {

                /*  We're going to have to do a sem_wait() with a timeout thread.
                 *  The thread will wait the specified time, then will issue a
                 *  SIGUSR2 signal that will interrupt the sem_wait() call.
                 *  We pass the thread the id of the current thread, the delay,
                 *  and the address of a flag to set on a timeout, so we can
                 *  distinguish an interrupt caused by the timeout thread from
                 *  one caused by some other signal.
                 */

                volatile short timedOut;                /* Flag to set on timeout */
                timeoutDetails details;     /* All the stuff the thread must know */
                struct sigaction oldSignalAction;       /* Current signal setting */
                pthread_t timeoutThread;                  /* Id of timeout thread */
                cleanupDetails cleaningDetails; /* What the cleanup routine needs */
                int oldCancelState;                /* Previous cancellation state */
                int ignoreCancelState;               /* Used in call, but ignored */
                int createStatus;              /* Status of pthread_create() call */

                /*  If the current thread is cancelled (and CML does do this)
                 *  we don't want to leave our timer thread running - if we've
                 *  started the thread we want to make sure we join it in order
                 *  to release its resources. So we set a cleanup handler to
                 *  do this. We pass it the address of the structure that will
                 *  hold all it needs to know. While we set all this up,
                 *  we prevent ourselves being cancelled, so all this data is
                 *  coherent.
                 */

                pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldCancelState);
                timeoutThread = (pthread_t)0;
                cleaningDetails.timedOutShort = &timedOut;
                cleaningDetails.threadIdAddr = &timeoutThread;
                cleaningDetails.sigHandlerAddr = &oldSignalAction;
                pthread_cleanup_push(timeoutThreadCleanup, &cleaningDetails);

                /*  Set up the details for the thread. Clear the timeout flag,
                 *  record the current SIGUSR2 action settings so we can restore
                 *  them later.
                 */

                details.delay.tv_sec = secsToWait;
                details.delay.tv_nsec = nsecsToWait;
                details.callingThread = pthread_self();
                details.timedOutShort = &timedOut;
                timedOut = CRNLIB_FALSE;
                sigaction(SIGUSR2, nullptr, &oldSignalAction);

                /*  Start up the timeout thread. Once we've done that, we can
                 *  restore the previous cancellation state.
                 */

                createStatus = pthread_create(&timeoutThread, nullptr, timeoutThreadMain, (void*)&details);
                pthread_setcancelstate(oldCancelState, &ignoreCancelState);

                if (createStatus < 0)
                {

                    /* Failed to create thread. errno will already be set properly */

                    result = -1;

                }
                else
                {
                    if (sem_wait(sem) == 0)
                    {
                        /*  Got the semaphore OK. We return zero, and all's well. */
                        result = 0;
                    }
                    else
                    {

                        /*  If we got a -1 error from sem_wait(), it may be because
                         *  it was interrupted by a timeout, or failed for some
                         *  other reason. We check for the expected timeout
                         *  condition, which is an 'interrupted' status and the
                         *  timeout flag set by the timeout thread. We report that as
                         *  a timeout error. Anything else is some other error and
                         *  errno is already set properly.
                         */

                        result = -1;
                        if (errno == EINTR)
                        {
                            if (timedOut) errno = ETIMEDOUT;
                        }
                    }

                }

                /*  The cleanup routine - timeoutThreadCleanup() - packages up
                 *  any tidying up that is needed, including joining with the
                 *  timer thread. This will be called if the current thread is
                 *  cancelled, but we need it to happen anyway, so we set the
                 *  execute flag true here as we remove it from the list of
                 *  cleanup routines to be called. So normally, this line amounts
                 *  to calling timeoutThreadCleanup().
                 */

                pthread_cleanup_pop(CRNLIB_TRUE);
            }
        }
    }

    return (result);
}

void timeoutThreadCleanup(void* passedPtr)
{
    /*  Get what we need from the structure we've been passed. */

    cleanupDetails* detailsPtr = (cleanupDetails*)passedPtr;
    short timedOut = *(detailsPtr->timedOutShort);
    pthread_t timeoutThread = *(detailsPtr->threadIdAddr);

    /*  If we created the thread, stop it - doesn't matter if it's no longer
     *  running, pthread_cancel can handle that. We make sure we wait for it
     *  to complete, because it is this pthread_join() call that releases any
     *  memory the thread may have allocated. Note that cancelling a thread is
     *  generally not a good idea, because of the difficulty of cleaning up
     *  after it, but this is a very simple thread that does nothing but call
     *  nanosleep(), and that we can cancel quite happily.
     */

    if (!timedOut)
    {
        pthread_cancel(timeoutThread);
    }
    pthread_join(timeoutThread, nullptr);

    /*  The code originally restored the old action handler, which generally
     *  was the default handler that caused the task to exit. Just occasionally,
     *  there seem to be cases where the signal is still queued and ready to
     *  trigger even though the thread that presumably sent it off just before
     *  it was cancelled has finished. I had thought that once we'd joined
     *  that thread, we could be sure of not seeing the signal, but that seems
     *  not to be the case, and so restoring a handler that will allow the task
     *  to crash is not a good idea, and so the line below has been commented
     *  out.
     *
     *  sigaction (SIGUSR2,detailsPtr->sigHandlerAddr,nullptr);
     */
}

static void* timeoutThreadMain(void* passedPtr)
{
    void* Return = (void*)0;

    /*  We grab all the data held in the calling thread right now. In some
     *  cases, we find that the calling thread has vanished and released
     *  its memory, including the details structure, by the time the timeout
     *  expires, and then we get an access violation when we try to set the
     *  'timed out' flag.
     */

    timeoutDetails details = *((timeoutDetails*)passedPtr);
    struct timespec requestedDelay = details.delay;

    /*  We do a nanosleep() for the specified delay, and then trigger a
     *  timeout. Note that we allow for the case where the nanosleep() is
     *  interrupted, and restart it for the remaining time. If the
     *  thread that is doing the sem_wait() call gets the semaphore, it
     *  will cancel this thread, which is fine as we aren't doing anything
     *  other than a sleep and a signal.
     */

    for (;;)
    {
        struct timespec remainingDelay;
        if (nanosleep(&requestedDelay, &remainingDelay) == 0)
        {
            break;
        }
        else if (errno == EINTR)
        {
            requestedDelay = remainingDelay;
        }
        else
        {
            Return = (void*)errno;
            break;
        }
    }

    /*  We've completed the delay without being cancelled, so we now trigger
     *  the timeout by sending a signal to the calling thread. And that's it,
     *  although we set the timeout flag first to indicate that it was us
     *  that interrupted the sem_wait() call. One precaution: before we
     *  try to set the timed-out flag, make sure the calling thread still
     *  exists - this may not be the case if things are closing down a bit
     *  messily. We check this quickly using a zero test signal.
     */

    if (pthread_kill(details.callingThread, 0) == 0)
    {
        *(details.timedOutShort) = CRNLIB_TRUE;
        if (triggerSignal(SIGUSR2, details.callingThread) < 0)
        {
            Return = (void*)errno;
        }
    }

    return Return;
}

static int triggerSignal(int Signal, pthread_t Thread)
{
    int Result = 0;
    struct sigaction SignalDetails;
    SignalDetails.sa_handler = ignoreSignal;
    SignalDetails.sa_flags = 0;
    (void)sigemptyset(&SignalDetails.sa_mask);
    if ((Result = sigaction(Signal, &SignalDetails, nullptr)) == 0)
    {
        Result = pthread_kill(Thread, Signal);
    }
    return Result;
}

static void ignoreSignal(int Signal)
{
    Signal = 0;
}
