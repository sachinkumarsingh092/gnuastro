/*********************************************************************
NoiseChisel - Detect and segment signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2018, Free Software Foundation, Inc.

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
#include <gnuastro/blank.h>
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

  /* Convovle with sharper kernel. */
  if(p->conv==NULL)
    {
      /* Make the convolved image. */
      if(!p->cp.quiet) gettimeofday(&t1, NULL);
      p->conv = gal_convolve_spatial(tl->tiles, p->kernel, p->cp.numthreads,
                                     1, tl->workoverch);

      /* Report and write check images if necessary. */
      if(!p->cp.quiet)
        {
          if(p->widekernel)
            gal_timing_report(&t1, "Convolved with sharper kernel.", 1);
          else
            gal_timing_report(&t1, "Convolved with given kernel.", 1);
        }
    }

  /* Set a fixed name for the convolved image (since it will be used in
     many check images). */
  if(p->conv->name) free(p->conv->name);
  gal_checkset_allocate_copy( ( p->widekernel
                                ? "CONVOLVED-SHARPER"
                                : "CONVOLVED" ), &p->conv->name);

  /* Save the convolution step if necessary. */
  if(p->detectionname)
    {
      gal_fits_img_write(p->input, p->detectionname, NULL, PROGRAM_NAME);
      gal_fits_img_write(p->conv, p->detectionname, NULL, PROGRAM_NAME);
    }

  /* Convolve with wider kernel (if requested). */
  if(p->widekernel)
    {
      if(!p->cp.quiet) gettimeofday(&t1, NULL);
      p->wconv=gal_convolve_spatial(tl->tiles, p->widekernel,
                                    p->cp.numthreads, 1, tl->workoverch);
      gal_checkset_allocate_copy("CONVOLVED-WIDER", &p->wconv->name);

      if(!p->cp.quiet)
        gal_timing_report(&t1, "Convolved with wider kernel.", 1);
    }
}





static void
noisechisel_find_sky_subtract(struct noisechiselparams *p)
{
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
/*************                   Output                  ***************/
/***********************************************************************/
/* The input image has been sky subtracted for further processing. So we'll
   need to copy the input image directly into the output. */
static void
noisechisel_output_copy_input(struct noisechiselparams *p)
{
  int status=0;
  fitsfile *in, *out;
  char card[FLEN_CARD];


  /* Create/open the output file. */
  out=gal_fits_open_to_write(p->cp.output);


  /* Open the input FITS file in the proper extension. */
  in=gal_fits_hdu_open(p->inputname, p->cp.hdu, READWRITE);


  /* Copy the input HDU into the output. */
  if( fits_copy_hdu(in, out, 0, &status) )
    gal_fits_io_error(status, "copying input hdu into first output hdu");


  /* If an extension name exists in the input HDU, then don't touch it. If
     the input doesn't have any, then make an EXTNAME keyword for it. Note
     that `fits_read_card' will return a non-zero if it doesn't find the
     keyword. */
  if( fits_read_card(out, "EXTNAME", card, &status) )
    {
      status=0;
      fits_write_key(out, TSTRING, "EXTNAME", "INPUT", "", &status);
    }


  /* Close the two files. */
  fits_close_file(in, &status);
  fits_close_file(out, &status);
}





/* Write the output file. */
static void
noisechisel_output(struct noisechiselparams *p)
{
  gal_fits_list_key_t *keys=NULL;

  /* Copy the input image into the first extension. */
  noisechisel_output_copy_input(p);


  /* Write the object labels and useful information into it's header. */
  if(p->onlydetection==0)
    gal_fits_key_list_add(&keys, GAL_TYPE_STRING, "WCLUMPS", 0, "yes", 0,
                          "Generate catalog with clumps?", 0, "bool");
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "NUMLABS", 0,
                        &p->numobjects, 0, "Total number of labels "
                        "(inclusive)", 0, "counter");
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "DETSN", 0, &p->detsnthresh,
                        0, "Minimum S/N of true pseudo-detections", 0,
                        "ratio");
  p->olabel->name = p->onlydetection ? "DETECTIONS" : "OBJECTS";
  gal_fits_img_write(p->olabel, p->cp.output, keys, PROGRAM_NAME);
  p->olabel->name=NULL;
  keys=NULL;


  /* Write the clumps labels and useful information into it's header. Note
     that to make the clumps image more easily viewable, we will set all
     sky pixels to blank. Only clump pixels that have an overlapping object
     pixel will be use anyway, so the sky pixels are irrelevant. */
  if(p->onlydetection==0)
    {
      p->clabel->name="CLUMPS";
      gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "NUMLABS", 0,
                            &p->numclumps, 0, "Total number of clumps", 0,
                            "counter");
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "CLUMPSN", 0,
                            &p->clumpsnthresh, 0, "Minimum S/N of true clumps",
                            0, "ratio");
      gal_fits_img_write(p->clabel, p->cp.output, keys, PROGRAM_NAME);
      p->clabel->name=NULL;
      keys=NULL;
    }


  /* Write the Sky image into the output */
  if(p->sky->name) free(p->sky->name);
  p->sky->name="SKY";
  gal_tile_full_values_write(p->sky, &p->cp.tl, 1, p->cp.output,
                             NULL, PROGRAM_NAME);
  p->sky->name=NULL;


  /* Write the Sky standard deviation into the output. */
  p->std->name="SKY_STD";
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "MAXSTD", 0, &p->maxstd, 0,
                        "Maximum raw tile standard deviation", 0,
                        p->input->unit);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "MINSTD", 0, &p->minstd, 0,
                        "Minimum raw tile standard deviation", 0,
                        p->input->unit);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "MEDSTD", 0, &p->medstd, 0,
                        "Median raw tile standard deviation", 0,
                        p->input->unit);
  gal_tile_full_values_write(p->std, &p->cp.tl, 1, p->cp.output, keys,
                             PROGRAM_NAME);
  p->std->name=NULL;
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

  /* If the user only wanted detection, ignore the segmentation steps. */
  if(p->onlydetection==0)
    {
      /* Correct the convolved image channel edges if necessary. */
      noisechisel_convolve_correct_ch_edges(p);

      /* Do the segmentation. */
      segmentation(p);
    }
  else
    p->numobjects=p->numdetections;

  /* Write the output. */
  noisechisel_output(p);
}
