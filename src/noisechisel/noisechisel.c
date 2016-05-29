/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/timing.h>
#include <gnuastro/arraymanip.h>
#include <gnuastro/fitsarrayvv.h>

#include "main.h"

#include "sky.h"
#include "label.h"
#include "binary.h"
#include "thresh.h"
#include "detection.h"
#include "noisechisel.h"
#include "segmentation.h"










/*********************************************************************/
/******************         NoiseChisel        ***********************/
/*********************************************************************/
void
makeoutput(struct noisechiselparams *p)
{
  float *fpt;
  double *dpt;
  long num[0];
  float *sky=NULL, *std=NULL;
  struct gal_fitsarray_header_ll *keys=NULL;
  size_t s0=p->smp.s0, s1=p->smp.s1;


  /* First put a copy of the input image. */
  gal_fitsarray_array_to_fits_img(p->cp.output, "Input", FLOAT_IMG, p->img,
                                  s0, s1, 0, p->wcs, NULL, SPACK_STRING);


  /* The object labels image with a keyword mentioning the number of
     objects. */
  num[0]=p->numobjects-1;
  gal_fitsarray_add_to_fits_header_ll(&keys, TLONG, "NOBJS", 0, num, 0,
                                      "Number of objects in the image.", 0,
                                      NULL);
  dpt=&p->detsn;
  gal_fitsarray_add_to_fits_header_ll(&keys, TDOUBLE, "DETSN", 0, dpt, 0,
                                      "Signal to noise of true "
                                      "pseudo-detections.", 0, NULL);
  gal_fitsarray_array_to_fits_img(p->cp.output, "Objects", LONG_IMG, p->olab,
                                  s0, s1, p->anyblank, p->wcs, keys,
                                  SPACK_STRING);
  keys=NULL;     /* keys was freed after writing. */


  /* The clump labels, mentioning the total number of clumps. Note
     that even if no segmentation is to be done, we need this
     extension filled so the sky and its standard deviation can remain
     on the 3rd and 4th extensions. */
  num[0] = p->detectonly ? 0 : p->numclumps-1;
  gal_fitsarray_add_to_fits_header_ll(&keys, TLONG, "NCLUMPS", 0, num, 0,
                                      "Number of clumps in the image.", 0,
                                      NULL);
  dpt=&p->clumpsn;
  gal_fitsarray_add_to_fits_header_ll(&keys, TDOUBLE, "CLUMPSN", 0, dpt, 0,
                                      "Signal to noise of true clumps.", 0,
                                      NULL);
  gal_fitsarray_array_to_fits_img(p->cp.output, "Clumps", LONG_IMG, p->clab, s0,
                                  s1, p->anyblank, p->wcs, keys, SPACK_STRING);
  keys=NULL;

  /* The sky image: */
  gal_mesh_check_garray(&p->smp, &sky, &std);
  gal_fitsarray_array_to_fits_img(p->cp.output, "Sky", FLOAT_IMG, sky,
                                  s0, s1, 0, p->wcs, NULL, SPACK_STRING);


  /* The sky standard deviation image. Note that since this is a
     linked list, the median will be written to the FITS header
     first. Also note that these values were found before
     interpolating or smoothing the image, so we have put a "raw" in
     the comments of each variable.*/
  fpt=&p->maxstd;
  gal_fitsarray_add_to_fits_header_ll(&keys, TFLOAT, "MAXSTD", 0, fpt, 0,
                                      "Maximum raw mesh sky standard "
                                      "deviation.", 0, NULL);
  fpt=&p->minstd;
  gal_fitsarray_add_to_fits_header_ll(&keys, TFLOAT, "MINSTD", 0, fpt, 0,
                                      "Minimum raw mesh sky standard "
                                      "deviation.", 0, NULL);
  fpt=&p->medstd;
  gal_fitsarray_add_to_fits_header_ll(&keys, TFLOAT, "MEDSTD", 0, fpt, 0,
                                      "Median raw mesh standard "
                                      "deviation.", 0, NULL);
  gal_fitsarray_array_to_fits_img(p->cp.output, "Standard deviation",
                                  FLOAT_IMG, std, s0, s1, 0, p->wcs,
                                  keys, SPACK_STRING);
  keys=NULL;


  /* Clean up: */
  free(sky);
  free(std);
}





void
noisechisel(struct noisechiselparams *p)
{
  struct gal_mesh_params *smp=&p->smp, *lmp=&p->lmp;

  float *imgcopy;
  struct timeval t1;
  int verb=p->cp.verb;
  size_t i, s0=smp->s0, s1=smp->s1;
  char report[GAL_TIMING_VERB_MSG_LENGTH_V], *oreport;


  /* Convolve the image: */
  if(verb) gettimeofday(&t1, NULL);
  gal_mesh_spatial_convolve_on_mesh(smp, &p->conv);
  if(p->detectionname)
    {
      gal_fitsarray_array_to_fits_img(p->detectionname, "Input", FLOAT_IMG,
                                      smp->img, s0, s1, p->anyblank, p->wcs,
                                      NULL, SPACK_STRING);
      gal_fitsarray_array_to_fits_img(p->detectionname, "Convolved", FLOAT_IMG,
                                      p->conv, s0, s1, p->anyblank, p->wcs,
                                      NULL, SPACK_STRING);
    }
  if(verb) gal_timing_report(&t1, "Convolved with kernel.", 1);



  /* Do the initial detection: */
  if(verb)
    {
      gal_timing_report(NULL, "Starting to find initial detections.", 1);
      gettimeofday(&t1, NULL);
    }
  initialdetection(p);
  if(verb)
    {
      sprintf(report, "%lu initial detections found.", p->numobjects-1);
      gal_timing_report(&t1, report, 1);
    }



  /* Remove the false detections */
  if(verb)
    {
      gal_timing_report(NULL, "Starting to find and remove false detections.",
                        1);
      gettimeofday(&t1, NULL);
    }
  onlytruedetections(p);
  if(verb)
    {
      sprintf(report, "%lu true detections identified.", p->numobjects-1);
      gal_timing_report(&t1, report, 1);
    }



  /* Dilate the byt array and find the new number of detections. Note
     that the connectivity has to be 8 connected. This is because we
     check the eight neighbors of every river pixel within each
     detection during the segmentation. Therefore if 4 connectivity is
     used here, two detection might be four connected and so their
     labels will be mixed during the checking of neighbors in the
     segmentation. */
  if(verb) gettimeofday(&t1, NULL);
  if(p->dilate)
    {
      for(i=0;i<p->dilate;++i)
        dilate0_erode1_8con(p->byt, s0, s1, 0);
      p->numobjects=BF_concmp(p->byt, p->olab, s0, s1, p->anyblank, 8);
      if(verb)
        {
          sprintf(report, "%lu detections after %lu dilation%s",
                  p->numobjects-1, p->dilate, p->dilate>1 ? "s." : ".");
          gal_timing_report(&t1, report, 1);
        }
    }
  if(p->detectionname)
    gal_fitsarray_array_to_fits_img(p->detectionname, "Dilated", LONG_IMG,
                                    p->olab, s0, s1, p->anyblank, p->wcs,
                                    NULL, SPACK_STRING);
  if(p->maskdetname)
    {
      gal_fitsarray_array_to_fits_img(p->maskdetname, "Input", FLOAT_IMG,
                                      p->img, s0, s1, p->anyblank, p->wcs,
                                      NULL, SPACK_STRING);
      gal_arraymanip_float_copy(p->img, s0*s1, &imgcopy);
      maskbackorforeground(imgcopy, s0*s1, p->byt, 0);
      gal_fitsarray_array_to_fits_img(p->maskdetname, "Undetected masked",
                                      FLOAT_IMG, imgcopy, s0, s1, p->anyblank,
                                      p->wcs, NULL, SPACK_STRING);
      free(imgcopy);
      gal_arraymanip_float_copy(p->img, s0*s1, &imgcopy);
      maskbackorforeground(imgcopy, s0*s1, p->byt, 1);
      gal_fitsarray_array_to_fits_img(p->maskdetname, "Detected masked",
                                      FLOAT_IMG, imgcopy, s0, s1, p->anyblank,
                                      p->wcs, NULL, SPACK_STRING);
      free(imgcopy);
    }

  /* Correct convolution if it was done on edges independently: */
  if(smp->nch>1 && smp->fullconvolution==0 && p->detectonly==0)
    {
      if(verb) gettimeofday(&t1, NULL);
      gal_mesh_change_to_full_convolution(smp, p->conv);
      if(verb)
        gal_timing_report(&t1, "Convolved image internals corrected.", 1);
    }

  /* Find and subtract the final sky value for the input and convolved
     images.

     VERY IMPORTANT:
     ###############

     The sky value for the input image should be found after the
     convolved image. This is because findavestdongrid, will also set
     the p->cpscorr value and we want that from the input image, not
     the convolved image.  */
  if(verb) gettimeofday(&t1, NULL);
  if(p->detectonly==0) findsubtractskyconv(p);
  findavestdongrid(p, p->skyname);
  if(p->detectonly==0) subtractskyimg(p);
  if(verb)
    gal_timing_report(&t1, "Final sky and its STD found and subtracted.", 1);




  /* Segment the detections if segmentation is to be done. */
  if(p->detectonly)
    clabwithnoseg(p->olab, p->clab, s0*s1, p->anyblank);
  else
    {
      if(verb)
        {
          gal_timing_report(NULL, "Starting to find clumps and objects.", 1);
          gettimeofday(&t1, NULL);
        }
      segmentation(p);
      if(verb)
        {
          sprintf(report, "%lu object%s""containing %lu clump%sfound.",
                  p->numobjects-1, p->numobjects==2 ? " " : "s ",
                  p->numclumps-1,  p->numclumps ==2 ? " " : "s ");
          gal_timing_report(&t1, report, 1);
        }
    }



  /* Make the output: */
  if(verb) gettimeofday(&t1, NULL);
  makeoutput(p);
  if(verb)
    {
      errno=0; oreport=malloc(strlen(p->cp.output)+100*sizeof *oreport);
      if(oreport==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for oreport in noisechisel "
              "(noisechisel.c)", strlen(p->cp.output)+100*sizeof *oreport);
      sprintf(oreport, "Output written to %s.", p->cp.output);
      gal_timing_report(&t1, oreport, 1);
      free(oreport);
    }

  /* Clean up: */
  free(p->conv);           /* Allocated in gal_mesh_spatial_convolve_on_mesh. */
  gal_mesh_free_mesh(smp); /* Allocated here.                                 */
  gal_mesh_free_mesh(lmp); /* Allocated here.                                 */
}
