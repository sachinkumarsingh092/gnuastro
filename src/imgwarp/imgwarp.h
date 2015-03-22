/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (gnuastro) package.

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
#ifndef IMGTRANSFORM_H
#define IMGTRANSFORM_H

#include "astrthreads.h"

#define MAXPOLYGONCORNERS 50

/* A function-like macro to multiply a 2 element vector with a
   transformation matrix and put the result in the 2 element output
   array. It is assumed that the input is flat. */
#define mappoint(V, T, O)                                       \
  {                                                             \
    (O)[0]=( ( (T)[0]*(V)[0] + (T)[1]*(V)[1] + (T)[2] )         \
             / ( (T)[6]*(V)[0] + (T)[7]*(V)[1] + (T)[8] ) );    \
    (O)[1]=( ( (T)[3]*(V)[0] + (T)[4]*(V)[1] + (T)[5] )         \
             / ( (T)[6]*(V)[0] + (T)[7]*(V)[1] + (T)[8] ) );    \
  }

struct iwpparams
{
  /* General input parameters: */
  struct imgwarpparams *p;

  /* Thread parameters. */
  size_t          *indexs;    /* Indexs to be used in this thread.     */
  pthread_barrier_t    *b;    /* Barrier to keep threads waiting.      */
};

void
imgwarp(struct imgwarpparams *p);

#endif
