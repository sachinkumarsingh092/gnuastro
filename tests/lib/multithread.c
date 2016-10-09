/*********************************************************************
A test program to multithreaded building using Gnuastro's helpers.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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
#include <stdio.h>
#include <stdlib.h>

#include "gnuastro/threads.h"


/* The params structure will keep the input information for each thread. As
   you see below, we will actually be defining an array of these
   structures, one for each thread. The reason for this is that the
   function that spins off threads only passes one argument, so as that
   argument, we will be passing the pointer to this structure. You can
   easily add new elements to this structure to use in the threads
   function. */
struct params
{
  size_t            id; /* Id of this thread.                            */
  size_t       *indexs; /* Indexes of actions to be done in this thread. */
  pthread_barrier_t *b; /* Pointer the barrier for all threads.          */
};





/* Worker function on threads */
void *
worker(void *inparam)
{
  /* The first thing to do is to say what the input pointer actually is. */
  struct params *prm=(struct params*)inparam;

  /* Now you can go onto do defining the function like any other
     function: first you define the variables and so on... */
  size_t i;

  /* Go over the jobs indexed for this thread: */
  for(i=0; prm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    {
      /* The indexes of the actions will make it possible to point to
         whatever data structure or input you want. So in this test, we
         will just print the thread ID and action id. */
      printf("thread %zu: %zu\n", prm->id, prm->indexs[i]);
    }

  /* Wait until all other threads finish. When there was only one thread,
     we explicitly set the pointer to the barrier structure to NULL, so
     only wait when a barrier is actually defined.*/
  if(prm->b)
    pthread_barrier_wait(prm->b);

  /* Return the NULL pointer. */
  return NULL;
}




/* This is the thread spinner function. */
int
main(void)
{
  int err;
  pthread_t t;          /* All thread ids saved in this, not used. */
  struct params *prm;
  size_t numbarriers;
  pthread_attr_t attr;
  pthread_barrier_t b;
  size_t i, *indexs, thrdcols;
  size_t numthreads=8, numactions=1000;

  /* Allocate the array of `param' structures. Note that in most cases, the
     number of threads will not be a constant like this simple case, it
     will be a variable passed to the thread-spinner. So we are using
     dynamic allocation for more general use as a tutorial. */
  prm=malloc(numthreads*sizeof *prm);
  if(prm==NULL)
    {
      fprintf(stderr, "%zu bytes could not be allocated for prm.",
              numthreads*sizeof *prm);
      exit(EXIT_FAILURE);
    }

  /* Distribute the actions into the threads: */
  gal_threads_dist_in_threads(numactions, numthreads, &indexs, &thrdcols);

  /* Do the job: when only one thread is necessary, there is no need to
     spin off one thread, just call the function directly (spinning off
     threads is expensive). This is for the generic thread spinner
     function, not this simple function where `numthreads' is a
     constant. */
  if(numthreads==1)
    {
      prm[0].id=0;
      prm[0].b=NULL;
      prm[0].indexs=indexs;
      worker(&prm[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
         (that spinns off the nt threads) is also a thread, so the
         number the barriers should be one more than the number of
         threads spinned off. */
      numbarriers = (numactions<numthreads ? numactions : numthreads) + 1;
      gal_threads_attr_barrier_init(&attr, &b, numbarriers);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
        if(indexs[i*thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            prm[i].id=i;
            prm[i].b=&b;
            prm[i].indexs=&indexs[i*thrdcols];
            err=pthread_create(&t, &attr, worker, &prm[i]);
            if(err)
              {
                fprintf(stderr, "can't create thread %zu", i);
                exit(EXIT_FAILURE);
              }
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }

  /* Clean up. */
  free(prm);
  free(indexs);

  /* Return. */
  return EXIT_SUCCESS;
}
