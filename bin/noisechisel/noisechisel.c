/*********************************************************************
NoiseChisel - Detect signal in a noisy dataset.
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
#include <errno.h>
#include <error.h>
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
      /* Do the convolution if a kernel was requested. */
      if(p->kernel)
        {
          /* Make the convolved image. */
          if(!p->cp.quiet) gettimeofday(&t1, NULL);
          p->conv = gal_convolve_spatial(tl->tiles, p->kernel,
                                         p->cp.numthreads, 1, tl->workoverch);

          /* Report and write check images if necessary. */
          if(!p->cp.quiet)
            {
              if(p->widekernel)
                gal_timing_report(&t1, "Convolved with sharper kernel.", 1);
              else
                gal_timing_report(&t1, "Convolved with given kernel.", 1);
            }
        }
      else
        p->conv=p->input;
    }

  /* Set a fixed name for the convolved image (since it will be used in
     many check images). */
  if(p->conv!=p->input)
    {
      if(p->conv->name) free(p->conv->name);
      gal_checkset_allocate_copy( ( p->widekernel
                                    ? "CONVOLVED-SHARPER"
                                    : "CONVOLVED" ), &p->conv->name);
    }

  /* Save the convolution step if necessary. */
  if(p->detectionname)
    {
      gal_fits_img_write(p->input, p->detectionname, NULL, PROGRAM_NAME);
      if(p->input!=p->conv)
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




















/***********************************************************************/
/*************                   Output                  ***************/
/***********************************************************************/
static void
noisechisel_params_in_keywords(struct noisechiselparams *p)
{
  gal_fits_list_key_t *keys=NULL;
  struct gal_tile_two_layer_params *tl=&p->cp.tl, *ltl=&p->ltl;

  /* Define the Keywords to write, in the same order as `main.h'. */
  gal_fits_key_write_filename("input", p->inputname, &keys);
  gal_fits_key_list_add(&keys, GAL_TYPE_STRING, "hdu", 0, p->cp.hdu,
                        0, "Extension name or number of input data.", 0,
                        NULL);
  if(p->kernelname)
    {
      gal_fits_key_write_filename("kernel", p->kernelname, &keys);
      if(p->khdu)
        gal_fits_key_list_add(&keys, GAL_TYPE_STRING, "khdu", 0, p->khdu,
                              0, "HDU/extension of kernel.", 0, NULL);
    }
  if(p->convolvedname)
    {
      gal_fits_key_write_filename("convolved", p->convolvedname, &keys);
      if(p->chdu)
        gal_fits_key_list_add(&keys, GAL_TYPE_STRING, "chdu", 0, p->chdu,
                              0, "HDU of convolved input.", 0, NULL);
    }
  if(p->widekernelname)
    {
      gal_fits_key_write_filename("widekernel", p->widekernelname, &keys);
      if(p->whdu)
        gal_fits_key_list_add(&keys, GAL_TYPE_STRING, "whdu", 0, p->whdu,
                              0, "HDU of wide kernel.", 0, NULL);
    }
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "tilesize_d1", 0,
                        &tl->tilesize[0], 0,
                        "Regular tile size on dim.1 (FITS order).", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "tilesize_d2", 0,
                        &tl->tilesize[1], 0,
                        "Regular tile size on dim.2 (FITS order).", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "largetilesize_d1", 0,
                        &ltl->tilesize[0], 0,
                        "Regular large tile size on dim.1 (FITS order).",
                        0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "largetilesize_d2", 0,
                        &ltl->tilesize[1], 0,
                        "Regular large tile size on dim.2 (FITS order).",
                        0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "numchannels_d1", 0,
                        &tl->numchannels[0], 0,
                        "No. of channels in dim.1 (FITS order).", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "numchannels_d2", 0,
                        &tl->numchannels[1], 0,
                        "No. of channels in dim.2 (FITS order).", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "remainderfrac", 0,
                        &tl->remainderfrac, 0,
                        "Fraction of remainder to split last tile.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_UINT8, "workoverch", 0,
                        &tl->workoverch, 0,
                        "Work (not tile) over channel edges.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_UINT8, "interponlyblank", 0,
                        &p->cp.interponlyblank, 0,
                        "Only interpolate over the blank tiles.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "interpnumngb", 0,
                        &p->cp.interpnumngb, 0,
                        "No. of neighbors to use for interpolation.",
                        0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "mirrordist", 0,
                        &p->mirrordist, 0,
                        "Max. dist. (error multip.) to find mode.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "modmedqdiff", 0,
                        &p->modmedqdiff, 0,
                        "Max. mode and median quant diff. per tile.",
                        0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "qthresh", 0,
                        &p->qthresh, 0,
                        "Quantile threshold on convolved image.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "qthreshtilequant", 0,
                        &p->qthreshtilequant, 0,
                        "Remove tiles at higher quantiles.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "smoothwidth", 0,
                        &p->smoothwidth, 0,
                        "Flat kernel width to smooth interpolated.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "erode", 0,
                        &p->erode, 0,
                        "Number of erosions after thresholding.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "erodengb", 0,
                        &p->erodengb, 0,
                        "4 or 8 connectivity in erosion.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "noerodequant", 0,
                        &p->noerodequant, 0,
                        "Quantile for no erosion.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "opening", 0,
                        &p->opening, 0,
                        "Depth of opening after erosion.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "openingngb", 0,
                        &p->openingngb, 0,
                        "4 or 8 connectivity in opening.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT64, "sigmaclipmultip", 0,
                        &p->sigmaclip[0], 0,
                        "Multiple of sigma for sigma-clipping.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT64, "sigmaclipend", 0,
                        &p->sigmaclip[1], 0,
                        "Termination criteria for sigma-clipping", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "minskyfrac", 0,
                        &p->minskyfrac, 0,
                        "Min. fraction of undetected area in tile.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "dthresh", 0,
                        &p->dthresh, 0,
                        "Sigma threshold for Pseudo-detections.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "snminarea", 0,
                        &p->snminarea, 0,
                        "Min. pseudo-detection area for S/N dist.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "minnumfalse", 0,
                        &p->minnumfalse, 0,
                        "Minimum number for S/N estimation.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "snquant", 0,
                        &p->snquant, 0,
                        "Quantile in pseudo-det. to define true.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "detgrowquant", 0,
                        &p->detgrowquant, 0,
                        "Minimum quant. to expand true detections.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "detgrowmaxholesize", 0,
                        &p->detgrowmaxholesize, 0,
                        "Max. area of holes after growth to fill.", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_UINT8, "cleangrowndet", 0,
                        &p->cleangrowndet, 0,
                        "Remove small S/N grown detections.", 0, NULL);

  /* Reverse the list and write them. */
  gal_fits_key_list_reverse(&keys);
  gal_fits_key_write_version(&keys, "NoiseChisel input parameters",
                             p->cp.output, "0");

}





/* Write the output file. */
static void
noisechisel_output(struct noisechiselparams *p)
{
  gal_fits_list_key_t *keys=NULL;


  /* Put a copy of the input into the output (when necessary). */
  if(p->rawoutput==0)
    {
      /* Subtract the Sky value. */
      sky_subtract(p);

      /* Correct the name of the input and write it out. */
      if(p->input->name) free(p->input->name);
      p->input->name="INPUT-NO-SKY";
      gal_fits_img_write(p->input, p->cp.output, NULL, PROGRAM_NAME);
      p->input->name=NULL;
    }


  /* Write the object labels and useful information into it's header. */
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "DETSN", 0, &p->detsnthresh,
                        0, "Minimum S/N of true pseudo-detections", 0,
                        "ratio");
  if(p->label)
    {
      gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "NUMLABS", 0,
                            &p->numdetections, 0, "Total number of labels "
                            "(inclusive)", 0, "counter");
      p->olabel->name = "DETECTIONS";
      gal_fits_img_write(p->olabel, p->cp.output, keys, PROGRAM_NAME);
      p->olabel->name=NULL;
    }
  else
    {
      p->binary->name = "DETECTIONS";
      gal_fits_img_write(p->binary, p->cp.output, keys, PROGRAM_NAME);
      p->binary->name=NULL;
    }
  keys=NULL;


  /* Write the Sky image into the output */
  if(p->sky->name) free(p->sky->name);
  p->sky->name="SKY";
  gal_tile_full_values_write(p->sky, &p->cp.tl, !p->ignoreblankinsky,
                             p->cp.output, NULL, PROGRAM_NAME);
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
  gal_tile_full_values_write(p->std, &p->cp.tl, !p->ignoreblankinsky,
                             p->cp.output, keys, PROGRAM_NAME);
  p->std->name=NULL;

  /* Write NoiseChisel's parameters as keywords into the first extension of
     the output. */
  noisechisel_params_in_keywords(p);

  /* Let the user know that the output is written. */
  if(!p->cp.quiet)
    printf("  - Output written to `%s'.\n", p->cp.output);
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

  /* If we have any detections, find the Sky value and subtract it from the
     input and convolved images. */
  if(p->numdetections)
    {
      /* Find the final Sky and Sky STD values. */
      sky_and_std(p, p->skyname);

      /* Abort if the user only wanted to see until this point.*/
      if(p->skyname && !p->continueaftercheck)
        ui_abort_after_check(p, p->skyname, NULL,
                             "derivation of final Sky (and its STD) value");

      /* Write the output. */
      noisechisel_output(p);
    }
  else
    {
      if(p->cp.quiet)
        error(0, 0, "no output file created: no detections could found "
              "in `%s' with given parameters", p->inputname);
      else
        gal_timing_report(NULL, "NO OUTPUT FILE CREATED (try with "
                          "`--checkdetection' to see why)", 1);
    }
}
