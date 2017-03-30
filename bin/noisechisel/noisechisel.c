/*********************************************************************
NoiseChisel - Detect and segment signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/convolve.h>

#include <timing.h>

#include "main.h"

#include "detection.h"




static void
noisechisel_convolve(struct noisechiselparams *p)
{
  struct timeval t1;
  struct gal_tile_two_layer_params *tl=&p->cp.tl;

  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  p->conv=gal_convolve_spatial(tl->tiles, p->kernel, p->cp.numthreads,
                               1, tl->workoverch);

  if(!p->cp.quiet) gal_timing_report(&t1, "Convolved with kernel.", 1);
  if(p->detectionname)
    {
      gal_fits_img_write(p->input, p->detectionname, NULL, PROGRAM_STRING);
      gal_fits_img_write(p->conv, p->detectionname, NULL, PROGRAM_STRING);
    }
}





void
noisechisel(struct noisechiselparams *p)
{
  /* Convolve the image. */
  noisechisel_convolve(p);

  /* Do the initial detection: */
  detection_initial(p);
}
