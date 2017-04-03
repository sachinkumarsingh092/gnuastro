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
#include <gnuastro/binary.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

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


  /* Label the connected components. */
  p->numobjects=gal_binary_connected_components(p->binary, &p->olabel, 1);
  if(p->detectionname)
    {
      p->olabel->name="OPENED-LABELED";
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
 ************            Pseudo detections           ************
 ****************************************************************/
/* Set all the pixels we don't need to Nan. */
static void
detection_pseudo_sky_or_det(struct noisechiselparams *p, uint8_t *w, int s0d1)
{
  uint32_t *l=p->olabel->array;
  uint8_t *b=p->binary->array, *bf=b+p->binary->size;

  if(s0d1)
    /* Set all sky regions (label equal to zero) to zero. */
    do *w++ = *l++ ? *b : 0; while(++b<bf);
  else
    /* Set all detected pixels (label larger than zero) to blank. */
    do *w++ = *l++ ? GAL_BLANK_UINT8 : *b; while(++b<bf);
}





/* We have the thresholded image (with blank values for regions that should
   not be used). Find the pseudo-detections in those regions. */
static size_t
detection_pseudo_find(struct noisechiselparams *p, gal_data_t *workbin,
                      gal_data_t *worklab, int s0d1)
{
  uint8_t *b, *bf;


  /* Set all the initial detected pixels to blank values. */
  detection_pseudo_sky_or_det(p, workbin->array, s0d1);
  if(p->detectionname)
    {
      workbin->name="DTHRESH-ON-SKY";
      gal_fits_img_write(workbin, p->detectionname, NULL, PROGRAM_STRING);
      workbin->name=NULL;
    }


  /* Fill the four-connected bounded holes. */
  gal_binary_fill_holes(workbin);
  if(p->detectionname)
    {
      workbin->name="HOLES-FILLED";
      gal_fits_img_write(workbin, p->detectionname, NULL, PROGRAM_STRING);
      workbin->name=NULL;
    }


  /* Open the image. */
  gal_binary_open(workbin, 1, 4, 1);
  if(p->detectionname)
    {
      workbin->name="OPENED";
      gal_fits_img_write(workbin, p->detectionname, NULL, PROGRAM_STRING);
      workbin->name=NULL;
    }


  /* Label all regions, but first, deal with the blank pixels in the
     `workbin' dataset when working on the Sky. Recall that the blank
     pixels are the other domain (detections if working on Sky and sky if
     working on detections). On the Sky image, blank should be set to 1
     (because we want the detected objects to have the same labels as the
     pseudo-detections that cover them). This will allow us to later remove
     these pseudo-detections. */
  if(s0d1==0)
    {
      bf=(b=workbin->array)+workbin->size;
      do if(*b==GAL_BLANK_UINT8) *b = !s0d1; while(++b<bf);
    }
  return gal_binary_connected_components(workbin, &worklab, 1);
}





#define PSN_EXTNAME "PSEUDOS-FOR-SN"
static gal_data_t *
detection_pseudo_sn(struct noisechiselparams *p, gal_data_t *worklab,
                         size_t num, int s0d1)
{
  float *snarr;
  uint8_t *flag;
  size_t tablen=num+1;
  gal_data_t *sn, *snind;
  uint32_t *plabend, *indarr;
  double ave, err, *xy, *brightness;
  struct gal_linkedlist_stll *comments=NULL;
  size_t ind, ndim=p->input->ndim, xyncols=1+ndim;
  size_t i, *area, counter=0, *dsize=p->input->dsize;
  size_t *coord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);
  float *img=p->input->array, *f=p->input->array, *ff=f+p->input->size;
  uint32_t *plab = worklab->array, *dlab = s0d1 ? NULL : p->olabel->array;


  /* Sanity check. */
  if(p->input->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "the input dataset to `detection_pseudo_sn' "
          "must be float32 type, it is %s",
          gal_type_to_string(p->input->type, 1));
  if(!isnan(GAL_BLANK_FLOAT32))
    error(EXIT_FAILURE, 0, "currently `detection_pseudo_sn' only "
          "recognizes a NaN value for blank floating point data types, the "
          "blank value is defined to be %f", GAL_BLANK_FLOAT32);
  if(ndim!=2)
    error(EXIT_FAILURE, 0, "currently `detection_pseudo_sn' only "
          "works on 2D datasets, your input is %zu dimensions", ndim);


  /* Allocate all the necessary arrays, note that since we want to put each
     object's information into the same index, the number of allocated
     spaces has to be `tablen=num+1'. */
  area       = gal_data_calloc_array(GAL_TYPE_SIZE_T,  tablen          );
  brightness = gal_data_calloc_array(GAL_TYPE_FLOAT64, tablen          );
  xy         = gal_data_calloc_array(GAL_TYPE_FLOAT64, xyncols*tablen  );
  flag       = s0d1==0 ? gal_data_calloc_array(GAL_TYPE_UINT8, tablen) : NULL;
  sn         = gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &tablen, NULL, 1,
                              p->cp.minmapsize, "SIGNAL-TO-NOISE", "ratio",
                              NULL);
  snind      = ( p->checkdetsn==0 ? NULL
                 : gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &tablen, NULL, 1,
                                  p->cp.minmapsize, "LABEL", "counter",
                                  NULL) );


  /* Go over all the pixels and get the necessary information. */
  do
    {
      /* All this work os only necessary when we are actually on a
         pseudo-detection label? */
      if(*plab)
        {
          /* For Sky pseudo-detections we'll start to see if it has already
             been determined that the object lies over a detected object or
             not. If it does, then just ignore it. */
          if(s0d1==0)
            {
              if( flag[*plab] ) { ++plab; ++dlab; continue; }
              else if(*dlab)    /* We are on a detection. */
                { flag[*plab]=1; area[*plab]=0; ++plab; ++dlab; continue; }
            }

          /* If we are on a blank pixel, ignore this pixel. */
          if( isnan(*f) ) { ++plab; if(s0d1==0) ++dlab; continue; }

          /* Save all the necessary values. */
          ++area[*plab];
          brightness[*plab] += *f;
          if( *f > 0.0f )  /* For calculatiing the approximate center, */
            {              /* necessary for calculating Sky and STD.   */
              xy[*plab*xyncols  ] += *f;
              xy[*plab*xyncols+1] += (double)((f-img)/dsize[1]) * *f;
              xy[*plab*xyncols+2] += (double)((f-img)%dsize[1]) * *f;
            }
        }

      /* Increment the other two labels. */
      ++plab;
      if(s0d1==0) ++dlab;
    }
  while(++f<ff);

  /* A small sanity check.
  {
    size_t i;
    for(i=1;i<num+1;++i)
      printf("%zu (%u): %-5zu %-13.3f %-13.3f %-13.3f %-13.3f\n", i, flag[i],
             area[i], brightness[i], xy[i*xyncols], xy[i*xyncols+1],
             xy[i*xyncols+2]);
  }
  */


  /* If the user wants to see the steps (on the background), remove all the
     pseudo-detections that will not be used in the final quantile
     calcluation. */
  if(p->detectionname)
    {
      plabend = (plab=worklab->array) + worklab->size;
      do
        if( *plab && ( area[*plab]<p->detsnminarea || brightness[*plab]<0) )
          *plab=0;
      while(++plab<plabend);
      worklab->name=PSN_EXTNAME;
      gal_fits_img_write(worklab, p->detectionname, NULL,
                         PROGRAM_STRING);
      worklab->name=NULL;
    }


  /* calculate the signal to noise for successful detections: */
  snarr=sn->array;
  if(snind) indarr=snind->array;
  if(s0d1) { snarr[0]=NAN; if(snind) indarr[0]=GAL_BLANK_UINT32; }
  for(i=1;i<num+1;++i)
    {
      ave=brightness[i]/area[i];
      if( area[i]>p->detsnminarea && ave>0.0f && xy[i*xyncols]>0.0f )
        {
          /* Get the flux weighted center coordinates. */
          coord[0]=GAL_DIMENSION_FLT_TO_INT( xy[i*xyncols+1]/xy[i*xyncols] );
          coord[1]=GAL_DIMENSION_FLT_TO_INT( xy[i*xyncols+2]/xy[i*xyncols] );

          /* Calculate the Sky and Sky standard deviation on this tile. */
          ave -= ((float *)(p->sky->array))[
                         gal_tile_full_id_from_coord(&p->cp.tl, coord) ];
          err  = ((float *)(p->std->array))[
                         gal_tile_full_id_from_coord(&p->cp.tl, coord) ];

          /* If the image was already sky subtracted, the second power of
             the error needs to be doubled. */
          err *= p->skysubtracted ? err : 2.0f*err;

          /* Correct the index in the sn to store the Signal to noise
             ratio. When we are dealing with the noise, we only want the
             non-zero signal to noise values, so we will just use a
             counter. But for initial detections, it is very important that
             their Signal to noise ratio be placed in the same index as
             their label. */
          ind = s0d1 ? i : counter++;
          if(snind) indarr[ind]=i;
          snarr[ind] = ( sqrt( (float)(area[i])/p->cpscorr )
                         * ave / sqrt(ave+err) );
        }
      else
        /* In detection pseudo-detections, order matters, so we will set
           all non-usable values to blank. */
        if(s0d1)
          {
            snarr[i]=NAN;
            if(snind) indarr[i]=GAL_BLANK_UINT32;;
          }
    }


  /* If we are in Sky mode, the sizes have to be corrected */
  if(s0d1==0)
    {
      sn->dsize[0]=sn->size=counter;
      if(snind) snind->dsize[0]=snind->size=counter;
    }


  /* If the user wanted a list of S/N values for all pseudo-detections,
     save it. */
  if(snind)
    {
      /* Make the comments, then write the table. */
      gal_linkedlist_add_to_stll(&comments, "See also: `"PSN_EXTNAME"' HDU "
                                 "of output with `--checkdetection'", 1);
      gal_linkedlist_add_to_stll(&comments, s0d1 ? "Pseudo-detection S/N "
                                 "over initially detected regions."
                                 : "Pseudo-detection S/N over initially "
                                 "undetected regions.", 1);
      threshold_write_sn_table(p, sn, snind, ( s0d1 ? p->detsn_d_name
                                               : p->detsn_s_name ), comments);
      gal_linkedlist_free_stll(comments, 1);
    }


  /* Clean up and return. */
  free(xy);
  free(area);
  free(brightness);
  if(flag) free(flag);
  if(snind) gal_data_free(snind);
  return sn;
}





static void
detection_pseudo_remove_low_sn(struct noisechiselparams *p,
                               gal_data_t *workbin, gal_data_t *worklab,
                               gal_data_t *sn, float snthresh)
{
  size_t i;
  float *snarr=sn->array;
  uint8_t *b=workbin->array;
  uint32_t *l=worklab->array, *lf=l+worklab->size;
  uint8_t *keep=gal_data_calloc_array(GAL_TYPE_UINT8, sn->size);

  /* Specify the new labels for those that must be kept/changed. Note that
     when an object didn't have an S/N, its S/N was given a value of NaN
     (which will fail on any condition), so it acts as if it had an S/N
     lower than the required value. */
  for(i=0;i<sn->size;++i)
    if( snarr[i] > snthresh ) keep[i]=1;


  /* Go over the pseudo-detection labels and only keep those that must be
     kept (using the new labels). */
  do *b++ = keep[ *l++ ] > 0; while(l<lf);


  /* If the user wanted to see the steps. */
  if(p->detectionname)
    {
      workbin->name="TRUE-PSEUDO-DETECTIONS";
      gal_fits_img_write(workbin, p->detectionname, NULL,
                         PROGRAM_STRING);
      workbin->name=NULL;
    }

  /* Clean up: */
  free(keep);
}






/* Find and do the necessary work on pseudo-detections. */
static gal_data_t *
detection_pseudo_real(struct noisechiselparams *p)
{
  char *msg;
  float snthresh;
  size_t numpseudo;
  struct timeval t1;
  gal_data_t *sn, *quant, *workbin, *worklab;

  /* Allocate the space for the working datasets. */
  worklab=gal_data_copy(p->olabel);
  workbin=gal_data_alloc(NULL, GAL_TYPE_UINT8, p->input->ndim,
                         p->input->dsize, p->input->wcs, 0, p->cp.minmapsize,
                         NULL, NULL, NULL);

  /* Over the Sky: find the pseudo-detections and make the S/N table. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  numpseudo=detection_pseudo_find(p, workbin, worklab, 0);
  sn=detection_pseudo_sn(p, worklab, numpseudo, 0);


  /* Get the S/N quantile and report it if we are in non-quiet mode. */
  quant=gal_statistics_quantile(sn, p->detquant, 1);
  snthresh = *((float *)(quant->array));
  if(!p->cp.quiet)
    {
      asprintf(&msg, "Pseudo-det S/N: %.2f (%.3f quant of %zu).",
               snthresh, p->detquant, sn->size);
      gal_timing_report(&t1, msg, 2);
    }
  gal_data_free(sn);
  gal_data_free(quant);


  /* Over the detections: find pseudo-detections and make S/N table. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  numpseudo=detection_pseudo_find(p, workbin, worklab, 1);
  sn=detection_pseudo_sn(p, worklab, numpseudo, 1);


  /* Remove the pseudo detections with a low S/N. */
  detection_pseudo_remove_low_sn(p, workbin, worklab, sn, snthresh);


  /* Clean up and return. */
  gal_data_free(sn);
  gal_data_free(worklab);
  return workbin;
}



















/****************************************************************
 ************        Removing false detections       ************
 ****************************************************************/
static size_t
detection_remove_false_initial(struct noisechiselparams *p,
                               gal_data_t *workbin)
{
  size_t i;
  uint8_t *b=workbin->array;
  uint32_t *l=p->olabel->array, *lf=l+p->olabel->size, curlab=1;
  uint32_t *newlabels=gal_data_calloc_array(GAL_TYPE_UINT32, p->numobjects+1);

  /* Find the new labels for all the existing labels. Recall that
     `newlabels' was initialized to zero, so any label that is not given a
     new label here will be automatically removed. After the first pixel of
     a label overlaps with dbyt[i], we don't need to check the rest of that
     object's pixels. At this point, tokeep is only binary: 0 or 1.

     Note that the zeroth element of tokeep can also be non zero, this is
     because the holes of the labeled regions might be filled during
     filling the holes, but have not been filled in the original labeled
     array. They are not important so you can just ignore them. */
  do
    {
      if(*l)
        {
          newlabels[ *l ] =
            newlabels[ *l ]     /* Have we already checked this label?    */
            ? 1                 /* Yes we have. Just set it to 1.         */
            : *b;               /* No we haven't, check pseudo-detection. */
        }
      ++b;
    }
  while(++l<lf);
  newlabels[0]=0;


  /* Now that we know which labels to keep, set the new labels for the
     detections that must be kept. */
  for(i=0;i<p->numobjects;++i) if(newlabels[i]) newlabels[i] = curlab++;


  /* Replace the byt and olab values with their proper values. If the
     user doesn't want to dilate, then change the labels in `lab'
     too. Otherwise, you don't have to worry about the label
     array. After dilation a new labeling will be done and the whole
     labeled array will be re-filled.*/
  b=workbin->array; l=p->olabel->array;
  if(p->dilate)
    do
      {
        if(*l!=GAL_BLANK_UINT32)
          *b = newlabels[ *l ] > 0;
        ++b;
      }
    while(++l<lf);
  else
    do
      if(*l!=GAL_BLANK_UINT32)
        *l = newlabels[ *l ];
    while(++l<lf);


  /* Clean up and return. */
  free(newlabels);
  return curlab-1;
}





/* The initial detection has been done, now we want to remove false
   detections. */
void
detection(struct noisechiselparams *p)
{
  char *msg;
  gal_data_t *workbin;
  struct timeval t0, t1;
  size_t num_true_initial;

  /* Report for the user. */
  if(!p->cp.quiet)
    {
      gal_timing_report(NULL, "Starting to find/remove false detections.", 1);
      gettimeofday(&t0, NULL);
    }


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


  /* Find the real pseudo-detections. */
  workbin=detection_pseudo_real(p);


  /* Only keep the initial detections that overlap with the real
     pseudo-detections. */
  num_true_initial=detection_remove_false_initial(p, workbin);
  if(!p->cp.quiet)
    {
      asprintf(&msg, "%zu false initial detections removed.",
               p->numobjects - num_true_initial);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }



  /* If the user asked for dilation, then apply it. */
  if(p->dilate)
    {
      gal_binary_dilate(workbin, p->dilate, 8, 1);
      p->numobjects = gal_binary_connected_components(workbin, &p->olabel, 8);
    }
  else p->numobjects=num_true_initial;
  if(!p->cp.quiet)
    {
      asprintf(&msg, "%zu detections after %zu dilation%s",
              p->numobjects, p->dilate, p->dilate>1 ? "s." : ".");
      gal_timing_report(&t0, msg, 1);
    }
  if(p->detectionname)
    {
      p->olabel->name="TRUE-INITIAL-DETECTIONS";
      gal_fits_img_write(p->olabel, p->detectionname, NULL,
                         PROGRAM_STRING);
      p->olabel->name=NULL;
    }



  /*  */


  /* If the user wanted to check the threshold and hasn't called
     `continueaftercheck', then stop NoiseChisel. */
  if(p->detectionname && !p->continueaftercheck)
    ui_abort_after_check(p, p->detectionname, "showing all detection steps");
}
