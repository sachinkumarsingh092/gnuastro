/*********************************************************************
Functions to facilitate using threads.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <gnuastro/blank.h>

/* When we are within Gnuastro's building process, 'IN_GNUASTRO_BUILD' is
   defined. In the build process, installation information (in particular
   'GAL_CONFIG_HAVE_PTHREAD_BARRIER' that we need below) is kept in
   'config.h'. When building a user's programs, this information is kept in
   'gnuastro/config.h'. Note that all '.c' files must start with the
   inclusion of 'config.h' and that 'gnuastro/config.h' is only created at
   installation time (not present during the building of Gnuastro).*/
#ifndef IN_GNUASTRO_BUILD
#include <gnuastro/config.h>
#endif

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





/*****************************************************************/
/*********    Implementation of pthread_barrier    ***************/
/*****************************************************************/
#if GAL_CONFIG_HAVE_PTHREAD_BARRIER == 0

/* Integer number of nano-seconds that 'pthread_barrier_destroy' should
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

#endif  /* GAL_CONFIG_HAVE_PTHREAD_BARRIER == 0 */





/*******************************************************************/
/************              Thread utilities           **************/
/*******************************************************************/
size_t
gal_threads_number();

void
gal_threads_dist_in_threads(size_t numactions, size_t numthreads,
                            size_t **outthrds, size_t *outthrdcols);

void
gal_threads_attr_barrier_init(pthread_attr_t *attr, pthread_barrier_t *b,
                              size_t limit);




/*******************************************************************/
/************     Run a function on multiple threads  **************/
/*******************************************************************/
struct gal_threads_params
{
  size_t            id; /* Id of this thread.                            */
  void         *params; /* Input structure for higher-level settings.    */
  size_t       *indexs; /* Indexes of actions to be done in this thread. */
  pthread_barrier_t *b; /* Pointer the barrier for all threads.          */
};

void
gal_threads_spin_off(void *(*worker)(void *), void *caller_params,
                     size_t numactions, size_t numthreads);


__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_THREADS_H__ */
