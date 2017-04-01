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
#include <gnuastro/binary.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "ui.h"
#include "threshold.h"





/****************************************************************
 ************           Initial detection            ************
 ****************************************************************/
void
detection_initial(struct noisechiselparams *p)
{
  char *msg;
  uint8_t *b, *bf;
  struct timeval t0, t1;


  /* Get the starting time. */
  if(!p->cp.quiet)
    {
      gal_timing_report(NULL, "Starting to find initial detections.", 1);
      gettimeofday(&t0, NULL);
    }


  /* Find and apply the threshold on the input. */
  threshold_quantile_find_apply(p);
  if(p->detectionname)
    {
      p->binary->name="THRESHOLDED";
      gal_fits_img_write(p->binary, p->detectionname, NULL, PROGRAM_STRING);
      p->binary->name=NULL;
    }


  /* Erode the image. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  gal_binary_erode(p->binary, p->erode, p->erodengb, 1);
  if(!p->cp.quiet)
    {
      asprintf(&msg, "Eroded %zu time%s (%zu-connectivity).", p->erode,
               p->erode>1?"s":"", p->erodengb);
      gal_timing_report(&t1, msg, 2);
    }
  if(p->detectionname)
    {
      p->binary->name="ERODED";
      gal_fits_img_write(p->binary, p->detectionname, NULL, PROGRAM_STRING);
      p->binary->name=NULL;
    }


  /* Correct the no-erode values. */
  bf=(b=p->binary->array)+p->binary->size; do *b = *b>0; while(++b<bf);


  /* Do the opening. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  gal_binary_open(p->binary, p->opening, p->openingngb, 1);
  if(!p->cp.quiet)
    {
      asprintf(&msg, "Opened (depth: %zu, %s connectivity).",
              p->opening, p->openingngb==4 ? "4" : "8");
      gal_timing_report(&t1, msg, 2);
    }
  if(p->detectionname)
    {
      p->binary->name="OPENED";
      gal_fits_img_write(p->binary, p->detectionname, NULL, PROGRAM_STRING);
      p->binary->name=NULL;
    }


  /* Label the connected components. */
  p->numobjects=gal_binary_connected_components(p->binary, &p->olabel, 1);
  if(p->detectionname)
    {
      p->olabel->name="LABELED";
      gal_fits_img_write(p->olabel, p->detectionname, NULL, PROGRAM_STRING);
      p->olabel->name=NULL;
    }


  /* Report the ending of initial detection. */
  if(!p->cp.quiet)
    {
      asprintf(&msg, "%zu initial detections found.", p->numobjects);
      gal_timing_report(&t0, msg, 1);
      free(msg);
    }
}




















/****************************************************************
 ************     Pseudo detection S/N threshold     ************
 ****************************************************************/
/* We have the thresholded image (with blank values for regions that should
   not be used). Find the pseudo-detections in those regions. */
void
detection_find_pseudos(struct noisechiselparams *p, gal_data_t workbin)
{

}




















/****************************************************************
 ************        Removing false detections       ************
 ****************************************************************/
/*  */
static void
detection_keep_sky_or_det(struct noisechiselparams *p, uint8_t *w, int s0d1)
{
  uint32_t *l=p->olabel->array;
  uint8_t *b=p->binary->array, *bf=b+p->binary->size;

  if(s0d1)
    /* Set all sky regions (label equal to zero) to blank. */
    do *w++ = *l++ ? *b : GAL_BLANK_UINT8; while(++b<bf);
  else
    /* Set all detected pixels (label larger than zero) to blank. */
    do *w++ = *l++ ? GAL_BLANK_UINT8 : *b; while(++b<bf);
}





/* The initial detection has been done, now we want to remove false
   detections. */
void
detection_remove_false(struct noisechiselparams *p)
{
  char *msg;
  struct timeval t1;
  gal_data_t *workbin;

  /* Find the Sky and its Standard Deviation from the initial detectios. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  threshold_sky_and_std(p);
  if(!p->cp.quiet)
    gal_timing_report(&t1, "Initial (crude) Sky and its STD found.", 2);


  /* Apply the sky threshold. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  threshold_apply(p, p->sky->array, p->std->array, THRESHOLD_SKY_STD);
  if(!p->cp.quiet)
    {
      asprintf(&msg, "Pseudo-detection thresh (%.3f sigma) applied.",
               p->dthresh);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* Allocate the space for the pseudo-detection threshold. */
  workbin=gal_data_alloc(NULL, GAL_TYPE_UINT8, p->input->ndim,
                         p->input->dsize, p->input->wcs, 0, p->cp.minmapsize,
                         NULL, NULL, NULL);


  /* Set all the initial detected pixels to blank values. */
  detection_keep_sky_or_det(p, workbin->array, 0);
  if(p->detectionname)
    {
      workbin->name="DTHRESH-ON-SKY";
      gal_fits_img_write(workbin, p->detectionname, NULL, PROGRAM_STRING);
      workbin->name=NULL;
    }


  /* Clean up. */
  gal_data_free(workbin);


  /* If the user wanted to check the threshold and hasn't called
     `continueaftercheck', then stop NoiseChisel. */
  if(p->detectionname && !p->continueaftercheck)
    ui_abort_after_check(p, p->detectionname, "showing all detection steps");
}
