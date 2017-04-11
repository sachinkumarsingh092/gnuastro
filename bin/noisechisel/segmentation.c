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
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/blank.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "ui.h"
#include "clumps.h"
#include "segmentation.h"


void
segmentation(struct noisechiselparams *p)
{
  uint32_t *l, *lf;
  float *f, snthresh;

  /* To keep the user up to date. */
  if(!p->cp.quiet)
    gal_timing_report(NULL, "Starting over-segmentation (finding clumps).",
                      1);


  /* If a check segmentation image was requested, then put in the
     inputs. */
  if(p->segmentationname)
    {
      gal_fits_img_write(p->input, p->segmentationname, NULL, PROGRAM_STRING);
      gal_fits_img_write(p->conv, p->segmentationname, NULL, PROGRAM_STRING);
      gal_fits_img_write(p->olabel, p->segmentationname, NULL,
                         PROGRAM_STRING);
    }


  /* Allocate the clump labels image. */
  p->clabel=gal_data_alloc(NULL, p->olabel->type, p->olabel->ndim,
                           p->olabel->dsize, p->olabel->wcs, 1,
                           p->cp.minmapsize, NULL, NULL, NULL);


  /* Set any possibly existing NaN values to blank. */
  f=p->input->array; lf=(l=p->clabel->array)+p->clabel->size;
  do if(isnan(*f++)) *l=GAL_BLANK_UINT32; while(++l<lf);


  /* Find the clumps over the un-detected regions of the input. */
  snthresh=clumps_on_undetected_sn(p);


  /* If the user wanted to check the segmentation and hasn't called
     `continueaftercheck', then stop NoiseChisel. */
  if(p->segmentationname && !p->continueaftercheck)
    ui_abort_after_check(p, p->segmentationname, NULL,
                         "showing all segmentation steps");
}
