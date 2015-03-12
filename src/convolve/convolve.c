/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (gnuastro) package.

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
#include <config.h>

#include <stdio.h>
#include <stdlib.h>

#include "astrthreads.h"
#include "fitsarrayvv.h"
#include "spatialconvolve.h"

#include "main.h"

void
frequencyconvolve()
{

}




















void
convolve(struct convolveparams *p)
{
  float *convolved;

  /* Do the convolution. */

  if(p->spatial)
    spatialconvolve(p->input, p->is0, p->is1, p->kernel, p->ks0,
                    p->ks1, p->cp.numthreads, p->edgecorrection,
                    &convolved);
  else
    frequencyconvolve();

  if(p->spatial)  /* TEMPORARY: until the frequency domain is completed,
                     this condition is temporary. */
    arraytofitsimg(p->cp.output, "Convolved", FLOAT_IMG, convolved,
                   p->is0, p->is1, p->wcs, SPACK_STRING);

  /* Free the output array: */
  free(convolved);
}
