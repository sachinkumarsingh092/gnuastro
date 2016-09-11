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
#ifndef __GAL_THREADS_H__
#define __GAL_THREADS_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <pthread.h>



/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */



/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */



/* Constant to use for non-existant index */
#define GAL_THREADS_NON_THRD_INDEX (size_t)(-1)



/*****************************************************************/
/*********    Implementation of pthread_barrier    ***************/
/*****************************************************************/
#ifndef HAVE_PTHREAD_BARRIER

/* Integer number of nano-seconds that `pthread_barrier_destroy' should
   wait for a check to see if all barriers have been reached. */
#define GAL_THREADS_BARRIER_DESTROY_NANOSECS 1000

typedef int pthread_barrierattr_t;

typedef struct
{
  pthread_mutex_t mutex;
  pthread_cond_t   cond;
  size_t          count;
  size_t          limit;
  size_t   condfinished;
} pthread_barrier_t;

int
pthread_barrier_init(pthread_barrier_t *b, pthread_barrierattr_t *attr,
                     unsigned int limit);

int
pthread_barrier_wait(pthread_barrier_t *b);

int
pthread_barrier_destroy(pthread_barrier_t *b);

#endif



/*****************************************************************/
/****************      gnuastro functions       ******************/
/*****************************************************************/
void
gal_threads_dist_in_threads(size_t nindexs, size_t nthrds, size_t **outthrds,
                            size_t *outthrdcols);

void
gal_threads_attr_barrier_init(pthread_attr_t *attr, pthread_barrier_t *b,
                              size_t numthreads);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_THREADS_H__ */
