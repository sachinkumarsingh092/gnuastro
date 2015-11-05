/*********************************************************************
SpatialConvolve - Convolve an image in the spatial domain.
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
#ifndef SPATIALCONVOLVE_H
#define SPATIALCONVOLVE_H





/* For pthread_barrier_t: */
#include "astrthreads.h"





/* Main structure: */
struct sconvparams
{
  /* General input parameters: */
  float           *input;     /* Input image array.                    */
  float          *kernel;     /* Kernel array.                         */
  float             *out;     /* Output image.                         */
  size_t             is0;     /* Image size along first C axis.        */
  size_t             is1;     /* Image size along second C axis.       */
  size_t             ks0;     /* Kernel size along first C axis.       */
  size_t             ks1;     /* Kernel size along second C axis.      */
  int     edgecorrection;     /* Correct the edges of the image.       */
  long       fpixel_i[2];     /* First pixel in input image.           */
  long       lpixel_i[2];     /* Last pixel in input image.            */
  long       fpixel_o[2];     /* First pixel in kernel.                */
  long       lpixel_o[2];     /* Last pixel in kernel.                 */

  /* Thread parameters. */
  size_t      numthreads;     /* Number of threads.                    */
  size_t         *indexs;     /* Indexs to be used in this thread.     */
  pthread_barrier_t   *b;     /* Barrier to keep threads waiting.      */
};





/* Functions: */
void
gal_spatialconvolve_pparams(float *input, size_t is0, size_t is1, float *kernel,
                            size_t ks0, size_t ks1, size_t nt,
                            int edgecorrection, float *out, size_t *indexs,
                            struct sconvparams *scp);

void *
gal_spatialconvolve_thread(void *inparam);

void
gal_spatialconvolve_convolve(float *input, size_t is0, size_t is1,
                             float *kernel, size_t ks0, size_t ks1,
                             size_t nt, int edgecorrection, float **out);

#endif
