/********************************************************************
mkprof (MakeProfiles) - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include "checkset.h"
#include "astrthreads.h"
#include "fitsarrayvv.h"

#include "main.h"
#include "mkprof.h"




/*

Basic Strategy
==============

N: Number of threads.


BUILDER THREADS
---------------

N-1 threads build profiles and fill an internal queue identical to the
global one but only for new models they have built since the last time
they added their constructs to the queue.

Upon building each profile, using pthread_mutex_trylock, they check to
see if they can lock the mutex for the global queue. If so, they add
their queue to the global one and free the mutex. Then they reset
their internal queue. If they can't get the lock, they continue adding
to their internal queue until they get the lock.


ADDER THREAD
------------

One thread is in charge of adding the separate profiles to the output
array using the global queue. Using a conditional, it waits until the
queue has some elements to start building. It starts building from the
end of the queue (note that the threads add elements from the start of
the queue). Every element it builts, it frees the array within it and
the element its self, then goes onto the next queue. If it reaches the
start of the queue and all the threads are finished, then the job is
done.


 */



















/**************************************************************/
/************            The builders             *************/
/**************************************************************/
/* Several threads (one minus the total number of threads) will call
   this function. It will build the profiles that are listed in its
   indexs parameter. */
void *
build(void *inparam)
{
  struct mkonthread *mkp=(struct mkonthread *)inparam;
  struct mkprofparams *p=mkp->p;

  size_t i;

  /* Make each profile that was specified for this thread. */
  for(i=0;mkp->indexs[i]!=NONTHRDINDEX;++i)
    {

    }

  /* Wait until all other threads finish. */
  if(p->cp.numthreads>1)
    pthread_barrier_wait(mkp->b);

  return NULL;
}




















/**************************************************************/
/************              The writer             *************/
/**************************************************************/
void
write(struct mkprofparams *p)
{

}




















/**************************************************************/
/************           Outside function          *************/
/**************************************************************/
void
mkprof(struct mkprofparams *p)
{
  int err;
  pthread_t t;		 /* Thread id not used, all are saved here. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct mkonthread *mkp;
  size_t nt=p->cp.numthreads;
  size_t i, *indexs, thrdcols;


  /* Allocate the arrays to keep the thread and parameters for each
     thread. Note that we only want nt-1 threads to do the
     building. */
  errno=0;
  mkp=malloc((nt-1)*sizeof *mkp);
  if(mkp==NULL)
    error(EXIT_FAILURE, errno,
	  "%lu bytes in mkprof (mkprof.c) for mkp", (nt-1)*sizeof *mkp);


  /* Distribute the different profiles for different threads. Note
     that one thread is left out for writing, while nt-1 are left
     for building. */
  distinthreads(p->cs0, nt-1, &indexs, &thrdcols);


  /* Build the profiles: */
  if(nt==1)
    {
      mkp[0].p=p;
      mkp[0].indexs=indexs;
      build(&mkp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that we want to build on nt-1
	 threads and write the arrays on one thread (this one that
	 spins off the nt-1 builder threads). So in total, the barrier
	 should stop nt threads. */
      attrbarrierinit(&attr, &b, nt);

      /* Initialize the condition variable and mutex. */
      err=pthread_mutex_init(&p->qlock, NULL);
      if(err) error(EXIT_FAILURE, 0, "Mutex not initialized.");
      err=pthread_cond_init(&p->qready, NULL);
      if(err) error(EXIT_FAILURE, 0, "Condition variable not initialized.");

      /* Spin off the threads: */
      for(i=0;i<nt-1;++i)
	if(indexs[i*thrdcols]!=NONTHRDINDEX)
	  {
	    mkp[i].p=p;
	    mkp[i].b=&b;
	    mkp[i].indexs=&indexs[i*thrdcols];
	    err=pthread_create(&t, &attr, build, &mkp[i]);
	    if(err)
	      error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
	  }
    }


  /* Write the created arrays into the image. */
  write(p);


  /* If numthreads>1, then wait for all the jobs to finish and destroy
     the attribute and barrier. */
  if(nt>1)
    {
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
      pthread_cond_destroy(&p->qready);
      pthread_mutex_destroy(&p->qlock);
    }
}
