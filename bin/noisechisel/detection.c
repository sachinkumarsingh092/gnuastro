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
#include <string.h>

#include <gnuastro/fits.h>
#include <gnuastro/binary.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "ui.h"
#include "sky.h"
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
      free(msg);
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
      free(msg);
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





/* Copy the space of this tile into the full/large array. */
static void
detection_write_in_large(gal_data_t *tile, gal_data_t *copy)
{
  uint8_t *c=copy->array;
  GAL_TILE_PARSE_OPERATE({*i=*c++;}, tile, NULL, 0, 0);
}





/* Fill the holes and open on multiple threads to find the
   pseudo-detections. Ideally both these should be done immediately after
   each other on each large tile, but when the user wants to check the
   steps, we need to break out of the threads at each step. */
struct fho_params
{
  int                    step;
  uint8_t          *copyspace;
  gal_data_t         *workbin;
  gal_data_t         *worklab;
  struct noisechiselparams *p;
};

static void *
detection_fill_holes_open(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct fho_params *fho_prm=(struct fho_params *)(tprm->params);
  struct noisechiselparams *p=fho_prm->p;

  void *tarray;
  gal_data_t *tile, *copy, *tblock;
  size_t i, dsize[]={1,1,1,1,1,1,1,1,1,1}; /* For upto 10-Dimensions! */


  /* A temporary data structure to wrap around the copy space. Note that
     the initially allocated space for this tile is only 1 pixel! */
  copy=gal_data_alloc(NULL, GAL_TYPE_UINT8, p->input->ndim, dsize,
                      NULL, 0, -1, NULL, NULL, NULL);
  free(copy->array);
  copy->array=&fho_prm->copyspace[p->maxltcontig*tprm->id];


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    {
      /* For easy reading. */
      tile=&p->ltl.tiles[tprm->indexs[i]];

      /* Change the tile pointers (temporarily). */
      tarray=tile->array;
      tblock=tile->block;
      tile->array=gal_tile_block_relative_to_other(tile, fho_prm->workbin);
      tile->block=fho_prm->workbin;

      /* Copy the tile into the contiguous patch of memory to work on, but
         first reset the size element so `gal_data_copy_to_allocated' knows
         there is enough space. */
      copy->size=p->maxltcontig;
      gal_data_copy_to_allocated(tile, copy);

      /* Fill the holes in this tile. */
      gal_binary_fill_holes(copy);
      if(fho_prm->step==1)
        {
          detection_write_in_large(tile, copy);
          tile->array=tarray;
          tile->block=tblock;
          continue;
        }

      /* Open all the regions. */
      gal_binary_open(copy, 1, 4, 1);

      /* Write the copied region back into the large input and AFTERWARDS,
         correct the tile's pointers, the pointers must not be corrected
         before writing the copy. */
      detection_write_in_large(tile, copy);
      tile->array=tarray;
      tile->block=tblock;
    }


  /* Clean up. */
  copy->array=NULL;
  gal_data_free(copy);


  /* Wait until all the threads finish and return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* We have the thresholded image (with blank values for regions that should
   not be used). Find the pseudo-detections in those regions. */
static size_t
detection_pseudo_find(struct noisechiselparams *p, gal_data_t *workbin,
                      gal_data_t *worklab, int s0d1)
{
  uint8_t *b, *bf;
  gal_data_t *bin;
  struct fho_params fho_prm={0, NULL, workbin, worklab, p};


  /* Set all the initial detected pixels to blank values. */
  detection_pseudo_sky_or_det(p, workbin->array, s0d1);
  if(p->detectionname)
    {
      workbin->name="DTHRESH-ON-SKY";
      gal_fits_img_write(workbin, p->detectionname, NULL, PROGRAM_STRING);
      workbin->name=NULL;
    }


  /* Allocate the space necessary to work on each tile (to avoid having to
     allocate it it separately for each tile and within each
     thread. `maxltcontig' is the maximum contiguous patch of memory needed
     to store all tiles. Finally, since we are working on a `uint8_t' type,
     the size of each element is only 1 byte. */
  fho_prm.copyspace=gal_data_malloc_array(GAL_TYPE_UINT8,
                                          p->cp.numthreads*p->maxltcontig);


  /* Fill the holes and open on each large tile. When no check image is
     requested, the two steps can be done independently on each tile, but
     when a check image is requested, we need to break out of the thread
     spinning function to save the full image then continue it. */
  if( p->detectionname )
    {
      /* Necessary initializations. */
      bin=gal_data_copy(workbin); /*  - Temporary array for demonstration.  */
      fho_prm.workbin=bin;        /*  - To pass onto the thread.            */
      fho_prm.step=1;             /*  - So we can break out of the threads. */

      /* Do each step. */
      while(fho_prm.step<3)
        {
          /* Put a copy of `workbin' into `bin' for every step (only
             necessary for the second step and after). For the first time
             it was already copied.*/
          if(fho_prm.step>1)
            memcpy(bin->array, workbin->array, workbin->size);

          /* Do the respective step. */
          gal_threads_spin_off(detection_fill_holes_open, &fho_prm,
                               p->ltl.tottiles, p->cp.numthreads);

          /* Set the extension name based on the step. */
          switch(fho_prm.step)
            {
            case 1:
              bin->name="HOLES-FILLED";
              break;
            case 2:
              bin->name="OPENED";
              break;
            default:
              error(EXIT_FAILURE, 0, "a bug! the value %d is not recognized "
                    "in `detection_pseudo_find'. Please contact us at %s so "
                    "we can address the issue", fho_prm.step,
                    PACKAGE_BUGREPORT);
            }

          /* Write the temporary array into the check image. */
          gal_fits_img_write(bin, p->detectionname, NULL, PROGRAM_STRING);

          /* Increment the step counter. */
          ++fho_prm.step;
        }

      /* Clean up: the array in `bin' should just be replaced with that in
         `workbin' because it is used in later steps. */
      free(workbin->array);
      workbin->array=bin->array;
      bin->name=bin->array=NULL;
      gal_data_free(bin);
    }
  else
    gal_threads_spin_off(detection_fill_holes_open, &fho_prm,
                         p->ltl.tottiles, p->cp.numthreads);


  /* Clean up. */
  free(fho_prm.copyspace);


  /* Label all regions, but first, deal with the blank pixels in the
     `workbin' dataset when working on the Sky. Recall that in this case,
     the blank pixels are the detections. On the Sky image, blank should be
     set to 1 (because we want the detected objects to have the same labels
     as the pseudo-detections that cover them). This will allow us to later
     remove these pseudo-detections. */
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
  uint32_t *plabend, *indarr=NULL;
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


  /* Calculate the signal to noise for successful detections: */
  snarr=sn->array;
  if(snind) indarr=snind->array;
  if(s0d1) { snarr[0]=NAN; if(snind) indarr[0]=GAL_BLANK_UINT32; }
  for(i=1;i<tablen;++i)
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

      /* Abort NoiseChisel if the user asked for it. */
      if(s0d1 && !p->continueaftercheck)
        ui_abort_after_check(p, p->detsn_s_name, p->detsn_d_name,
                             "pseudo-detection S/N values in a table");
    }


  /* Clean up and return. */
  free(xy);
  free(area);
  free(coord);
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
      free(msg);
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
  else                               /* We need the binary array even when */
    do                               /* there is no dilation: the binary   */
      {                              /* array is used for estimating the   */
        if(*l!=GAL_BLANK_UINT32)     /* Sky and its STD. */
          *b = ( *l = newlabels[ *l ] ) > 0;
        ++b;
      }
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
  sky_and_std(p, p->detskyname);
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
      free(msg);
    }
  if(p->detectionname)
    {
      p->olabel->name="TRUE-INITIAL-DETECTIONS";
      gal_fits_img_write(p->olabel, p->detectionname, NULL,
                         PROGRAM_STRING);
      p->olabel->name=NULL;
    }


  /* p->binary was used to keep the initial pseudo-detection threshold. But
     we don't need it any more, so we'll just free it and put the `workbin'
     array in its place. Note that `workbin' has a map of all the detected
     objects, which is still necessary during NoiseChisel. */
  gal_data_free(p->binary);
  p->binary=workbin;


  /* If the user wanted to check the threshold and hasn't called
     `continueaftercheck', then stop NoiseChisel. */
  if(p->detectionname && !p->continueaftercheck)
    ui_abort_after_check(p, p->detectionname, NULL,
                         "showing all detection steps");
}
