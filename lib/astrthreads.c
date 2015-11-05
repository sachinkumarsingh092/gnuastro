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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include "astrthreads.h"







/*****************************************************************/
/*********    Implementation of pthread_barrier    ***************/
/*****************************************************************/
/* Re-implementation of the example code given in:
http://blog.albertarmea.com/post/47089939939/using-pthread-barrier-on-mac-os-x
 */
#ifndef HAVE_PTHREAD_BARRIER
int
pthread_barrier_init(pthread_barrier_t *b, pthread_barrierattr_t *attr,
                     unsigned int count)
{
  int err, junk=*attr;

  /* Sanity check: */
  junk=junk+1;               /* So there is no unused variable warning. */
  if(count==0)
    {
      errno = EINVAL;
      error(EXIT_FAILURE, errno, "in pthread_barrier_init, count is zero");
    }

  /* Initialize the mutex: */
  err=pthread_mutex_init(&b->mutex, 0);
  if(err)
    error(EXIT_FAILURE, err, "Inializing mutex in pthread_barrier_init");

  /* Initialize the condition variable: */
  err=pthread_cond_init(&b->cond, 0);
  if(err)
    {
      pthread_mutex_destroy(&b->mutex);
      error(EXIT_FAILURE, err, "Inializing cond in pthread_barrier_init");
    }

  /* set the values: */
  b->tripCount=count;
  b->count=0;

  return 0;
}





int
pthread_barrier_destroy(pthread_barrier_t *b)
{
  pthread_cond_destroy(&b->cond);
  pthread_mutex_destroy(&b->mutex);
  return 0;
}





int
pthread_barrier_wait(pthread_barrier_t *b)
{
  pthread_mutex_lock(&b->mutex);
  ++(b->count);
  if(b->count >= b->tripCount)
    {
      b->count = 0;
      pthread_cond_broadcast(&b->cond);
      pthread_mutex_unlock(&b->mutex);
      return 1;
    }
  else
    {
      pthread_cond_wait(&b->cond, &(b->mutex));
      pthread_mutex_unlock(&b->mutex);
      return 0;
    }
}
#endif




















/*******************************************************************/
/************     Distribute job indexs in threads    **************/
/*******************************************************************/
/* We have `nindexs` jobs and we want their indexs to be divided
   between `nthrds` CPU threads. This function will give each index to
   a thread such that the maximum difference between the number of
   images for each thread is 1. The results will be saved in a 2D
   array of `outthrdcols` columns and each row will finish with a
   (size_t) -1, which is larger than any possible index!. */
void
gal_threads_dist_in_threads(size_t nindexs, size_t nthrds, size_t **outthrds,
	      size_t *outthrdcols)
{
  size_t *sp, *fp;
  size_t i, *thrds, thrdcols;
  *outthrdcols = thrdcols = nindexs/nthrds+2;

  errno=0;
  thrds=*outthrds=malloc(nthrds*thrdcols*sizeof *thrds);
  if(thrds==NULL)
    error(EXIT_FAILURE, errno, "Allocating thrds in prepindexsinthreads");

  /* Initialize all the elements to NONINDEX. */
  fp=(sp=thrds)+nthrds*thrdcols;
  do *sp=NONTHRDINDEX; while(++sp<fp);

  /* Distribute the labels in the threads.  */
  for(i=0;i<nindexs;++i)
    thrds[ (i%nthrds)*thrdcols+(i/nthrds) ] = i;

  /* In case you want to see the result:
  for(i=0;i<nthrds;++i)
    {
      size_t j;
      printf("\n\n############################\n");
      printf("THREAD %lu: \n", i);
      for(j=0;thrds[i*thrdcols+j]!=NONTHRDINDEX;j++)
	printf("%lu, ", thrds[i*thrdcols+j]);
      printf("\b\b.\n");
    }
  exit(0);
  */
}





void
gal_threads_attr_barrier_init(pthread_attr_t *attr, pthread_barrier_t *b,
		size_t numthreads)
{
  int err;

  err=pthread_attr_init(attr);
  if(err) error(EXIT_FAILURE, 0, "Thread attr not initialized.");
  err=pthread_attr_setdetachstate(attr, PTHREAD_CREATE_DETACHED);
  if(err) error(EXIT_FAILURE, 0, "Thread attr not detached.");
  err=pthread_barrier_init(b, NULL, numthreads);
  if(err) error(EXIT_FAILURE, 0, "Thread barrier not initialized.");
}
