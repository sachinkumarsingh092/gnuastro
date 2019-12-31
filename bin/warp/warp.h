/*********************************************************************
Warp - Warp images using projective mapping.
Warp is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef IMGTRANSFORM_H
#define IMGTRANSFORM_H

#include <gnuastro/threads.h>


/* Limits to account for floating point errors: */
#define ABSOLUTEFLTERROR 1e-10
#define RELATIVEFLTERROR 1e-6


/* Internal structure. */
struct iwpparams
{
  /* General input parameters: */
  struct warpparams *p;

  /* Thread parameters. */
  size_t          *indexs;    /* Indexs to be used in this thread.     */
  pthread_barrier_t    *b;    /* Barrier to keep threads waiting.      */
};


/* Extenal functions. */
void
warp(struct warpparams *p);

#endif
