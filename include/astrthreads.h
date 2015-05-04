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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef ASTRTHREADS_H
#define ASTRTHREADS_H





#include <pthread.h>





#define NONTHRDINDEX (size_t)(-1)










/*****************************************************************/
/*********    Implementation of pthread_barrier    ***************/
/*****************************************************************/
#ifndef HAVE_PTHREAD_BARRIER

typedef int pthread_barrierattr_t;

typedef struct
{
  pthread_mutex_t mutex;
  pthread_cond_t   cond;
  size_t          count;
  size_t      tripCount;
} pthread_barrier_t;

int
pthread_barrier_init(pthread_barrier_t *b, pthread_barrierattr_t *attr,
                     unsigned int count);
int
pthread_barrier_destroy(pthread_barrier_t *b);

int
pthread_barrier_wait(pthread_barrier_t *b);

#endif





















/*****************************************************************/
/****************      gnuastro functions       ******************/
/*****************************************************************/
void
distinthreads(size_t nindexs, size_t nthrds, size_t **outthrds,
	      size_t *outthrdcols);

void
attrbarrierinit(pthread_attr_t *attr, pthread_barrier_t *b,
		size_t numthreads);

#endif
