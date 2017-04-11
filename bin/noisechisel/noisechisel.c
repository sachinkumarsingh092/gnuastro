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

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "ui.h"
#include "sky.h"
#include "detection.h"
#include "threshold.h"
#include "segmentation.h"










/***********************************************************************/
/*************  Wrapper functions (for clean high-level) ***************/
/***********************************************************************/
static void
noisechisel_convolve(struct noisechiselparams *p)
{
  struct timeval t1;
  struct gal_tile_two_layer_params *tl=&p->cp.tl;

  /* Do the convolution. */
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





static void
noisechisel_find_sky_subtract(struct noisechiselparams *p)
{
  /* Free the initial Sky and its STD estimates. */
  gal_data_free(p->sky);
  gal_data_free(p->std);

  /* Find the final Sky value. */
  sky_and_std(p, p->skyname);

  /* Abort if the user only wanted to see until this point.*/
  if(p->skyname && !p->continueaftercheck)
    ui_abort_after_check(p, p->skyname, NULL,
                         "derivation of final Sky (and its STD) value");

  /* Subtract the Sky from the Input and Convolved (necessary for
     segmentation) images. */
  sky_subtract(p);
}





/* If convolution was not done while respecting channel edges (when there
   is more than one channel, pixels outside the edge weren't used in the
   convolution), then correct it. */
static void
noisechisel_convolve_correct_ch_edges(struct noisechiselparams *p)
{
  struct gal_tile_two_layer_params *tl=&p->cp.tl;

  /* Correct the convolved image if necessary. */
  if( tl->totchannels>1 && tl->workoverch==0 )
    {
      /* Do the correction. */
      gal_convolve_spatial_correct_ch_edge(tl->tiles, p->kernel,
                                           p->cp.numthreads, 1, p->conv);

      /* Inform the user. */
      if(!p->cp.quiet)
        gal_timing_report(NULL, "Corrected convolution of touching channel "
                          "edges", 1);
    }
}




















/***********************************************************************/
/*************             High level function           ***************/
/***********************************************************************/
void
noisechisel(struct noisechiselparams *p)
{
  /* Convolve the image. */
  noisechisel_convolve(p);

  /* Do the initial detection. */
  detection_initial(p);

  /* Remove false detections. */
  detection(p);

  /* Find the Sky value and subtract it from the input and convolved
     images. */
  noisechisel_find_sky_subtract(p);

  /* Correct the convolved image channel edges if necessary. */
  noisechisel_convolve_correct_ch_edges(p);

  /* Do the segmentation. */
  segmentation(p);
}
