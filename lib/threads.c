/*********************************************************************
Functions to facilitate using threads.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <time.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/threads.h>

#include <nproc.h>         /* from Gnulib, in Gnuastro's source */





/*****************************************************************/
/*********    Implementation of pthread_barrier    ***************/
/*****************************************************************/
/* Re-implementation of the example code given in:
http://blog.albertarmea.com/post/47089939939/using-pthread-barrier-on-mac-os-x
 */
#if GAL_CONFIG_HAVE_PTHREAD_BARRIER == 0

/* Initialize the barrier structure. A barrier is a high-level way to wait
   until several threads have finished. */
int
pthread_barrier_init(pthread_barrier_t *b, pthread_barrierattr_t *attr,
                     unsigned int limit)
{
  int err;

  /* Sanity check: */
  if(limit==0)
    {
      errno = EINVAL;
      error(EXIT_FAILURE, errno, "in pthread_barrier_init, limit is zero");
    }

  /* Initialize the mutex: */
  err=pthread_mutex_init(&b->mutex, 0);
  if(err)
    error(EXIT_FAILURE, err, "inializing mutex in pthread_barrier_init");

  /* Initialize the condition variable: */
  err=pthread_cond_init(&b->cond, 0);
  if(err)
    {
      pthread_mutex_destroy(&b->mutex);
      error(EXIT_FAILURE, err, "inializing cond in pthread_barrier_init");
    }

  /* set the values: */
  b->limit=limit;
  b->condfinished=b->count=0;

  return 0;
}





/* Suspend the calling thread (tell it to wait), until the limiting number
   of barriers is reached by the other threads. When the number isn't
   reached yet (process goes into the `else'), then we use the
   `pthread_cond_wait' function, which will wait until a broadcast is
   announced by another thread that succeeds the `if'. After the thread no
   longer needs the condition variable, we increment `b->condfinished' so
   `pthread_barrier_destroy' can know if it should wait (sleep) or
   continue.*/
int
pthread_barrier_wait(pthread_barrier_t *b)
{
  pthread_mutex_lock(&b->mutex);
  ++b->count;
  if(b->count >= b->limit)
    {
      pthread_cond_broadcast(&b->cond);
      ++b->condfinished;
      pthread_mutex_unlock(&b->mutex);
      return 1;
    }
  else
    {
      /* Initially b->count is smaller than b->limit, otherwise control
         would never have reached here. However, when it get to
         `pthread_cond_wait', this thread goes into a suspended period and
         is only awaken when a broad-cast is made. But we only want this
         thread to finish when b->count >= b->limit, so we have this while
         loop. */
      while(b->count < b->limit)
        pthread_cond_wait(&b->cond, &b->mutex);
      ++b->condfinished;
      pthread_mutex_unlock(&b->mutex);
      return 0;
    }
}





/* Wait until all the threads that were blocked behind this barrier have
   finished (don't need the mutex and condition variable anymore. Then
   destroy the two. After this function, you can safely re-use the
   pthread_barrier_t structure. */
int
pthread_barrier_destroy(pthread_barrier_t *b)
{
  struct timespec request, remaining;

  /* Wait until no more threads need the barrier. */
  request.tv_sec=0;
  request.tv_nsec=GAL_THREADS_BARRIER_DESTROY_NANOSECS;
  while( b->condfinished < b->limit )
    nanosleep(&request, &remaining);

  /* Destroy the condition variable and mutex */
  pthread_cond_destroy(&b->cond);
  pthread_mutex_destroy(&b->mutex);
  return 0;
}

#endif  /* GAL_CONFIG_HAVE_PTHREAD_BARRIER == 0 */




















/*******************************************************************/
/************              Thread utilities           **************/
/*******************************************************************/
size_t
gal_threads_number()
{
  return num_processors(NPROC_CURRENT);
}





/* We have `numactions` jobs and we want their indexs to be divided
   between `numthreads` CPU threads. This function will give each index to
   a thread such that the maximum difference between the number of
   images for each thread is 1. The results will be saved in a 2D
   array of `outthrdcols` columns and each row will finish with a
   (size_t) -1, which is larger than any possible index!. */
void
gal_threads_dist_in_threads(size_t numactions, size_t numthreads,
                            size_t **outthrds, size_t *outthrdcols)
{
  size_t *sp, *fp;
  size_t i, *thrds, thrdcols;
  *outthrdcols = thrdcols = numactions/numthreads+2;

  errno=0;
  thrds=*outthrds=malloc(numthreads*thrdcols*sizeof *thrds);
  if(thrds==NULL)
    error(EXIT_FAILURE, errno, "allocating thrds in prepindexsinthreads");

  /* Initialize all the elements to NONINDEX. */
  fp=(sp=thrds)+numthreads*thrdcols;
  do *sp=GAL_THREADS_NON_THRD_INDEX; while(++sp<fp);

  /* Distribute the labels in the threads.  */
  for(i=0;i<numactions;++i)
    thrds[ (i%numthreads)*thrdcols+(i/numthreads) ] = i;

  /* In case you want to see the result:
  for(i=0;i<numthreads;++i)
    {
      size_t j;
      printf("\n\n############################\n");
      printf("THREAD %zu: \n", i);
      for(j=0;thrds[i*thrdcols+j]!=GAL_THREADS_NON_THRD_INDEX;j++)
        printf("%zu, ", thrds[i*thrdcols+j]);
      printf("\b\b.\n");
    }
  exit(0);
  */
}





void
gal_threads_attr_barrier_init(pthread_attr_t *attr, pthread_barrier_t *b,
                              size_t limit)
{
  int err;

  err=pthread_attr_init(attr);
  if(err) error(EXIT_FAILURE, 0, "thread attr not initialized");
  err=pthread_attr_setdetachstate(attr, PTHREAD_CREATE_DETACHED);
  if(err) error(EXIT_FAILURE, 0, "thread attr not detached");
  err=pthread_barrier_init(b, NULL, limit);
  if(err) error(EXIT_FAILURE, 0, "thread barrier not initialized");
}
