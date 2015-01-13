/*********************************************************************
Functions to facilitate using threads.
This is part of GNU Astronomy Utilities (AstrUtils) package.

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
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <pthread.h>

#include "astrthreads.h"










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
distinthreads(size_t nindexs, size_t nthrds, size_t **outthrds,
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
attrbarrierinit(pthread_attr_t *attr, pthread_barrier_t *b,
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
