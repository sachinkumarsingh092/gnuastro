/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef CONVOLVE_H
#define CONVOLVE_H

#include <gnuastro/threads.h>
#include <gsl/gsl_fft_complex.h>

struct fftonthreadparams
{
  /* Operating info: */
  size_t                id; /* The number of this thread.               */
  struct convolveparams *p; /* Pointer to main program structure.       */
  int   forward1backwardn1; /* Operate on one or two images.            */
  size_t            stride; /* 1D FFT on rows or columns?               */

  /* Pointers to GSL FFT structures: */
  gsl_fft_complex_wavetable *ps0wave;
  gsl_fft_complex_wavetable *ps1wave;
  gsl_fft_complex_workspace *ps0work;
  gsl_fft_complex_workspace *ps1work;

  /* Thread parameters. */
  size_t          *indexs;  /* Indexs to be used in this thread.        */
  pthread_barrier_t    *b;  /* Barrier to keep threads waiting.         */
};


void
convolve(struct convolveparams *p);


#endif
