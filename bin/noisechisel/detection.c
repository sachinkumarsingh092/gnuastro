/*********************************************************************
NoiseChisel - Detect signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <gnuastro/label.h>
#include <gnuastro/blank.h>
#include <gnuastro/binary.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

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
  float *f;
  char *msg;
  uint8_t *b, *bf;
  int resetblank=0;
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
      gal_fits_img_write(p->binary, p->detectionname, NULL, PROGRAM_NAME);
      p->binary->name=NULL;
    }


  /* Remove any blank elements from the binary image if requested. */
  if(p->blankasforeground==0 && gal_blank_present(p->binary,0))
    {
      resetblank=1;
      bf=(b=p->binary->array)+p->binary->size;
      do *b = *b==GAL_BLANK_UINT8 ? 0 : *b; while(++b<bf);
    }


  /* Erode the image. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  gal_binary_erode(p->binary, p->erode, p->erodengb==4 ? 1 : 2, 1);
  if(!p->cp.quiet)
    {
      if( asprintf(&msg, "Eroded %zu time%s (%zu-connected).", p->erode,
                   p->erode!=1?"s":"", p->erodengb)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }
  if(p->detectionname)
    {
      p->binary->name="ERODED";
      gal_fits_img_write(p->binary, p->detectionname, NULL, PROGRAM_NAME);
      p->binary->name=NULL;
    }


  /* Correct the no-erode values. */
  bf=(b=p->binary->array)+p->binary->size;
  do *b = *b==THRESHOLD_NO_ERODE_VALUE ? 1 : *b; while(++b<bf);


  /* Do the opening. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  gal_binary_open(p->binary, p->opening, p->openingngb==4 ? 1 : 2, 1);
  if(!p->cp.quiet)
    {
      if( asprintf(&msg, "Opened (depth: %zu, %zu-connected).",
                   p->opening, p->openingngb)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* Reset the blank values (if requested). */
  if(resetblank)
    {
      f=p->input->array;
      bf=(b=p->binary->array)+p->binary->size;
      do *b = isnan(*f++) ? GAL_BLANK_UINT8 : *b; while(++b<bf);
    }


  /* Label the connected components. */
  p->numinitialdets=gal_binary_connected_components(p->binary, &p->olabel,
                                                    p->binary->ndim);
  if(p->detectionname)
    {
      p->olabel->name="OPENED_AND_LABELED";
      gal_fits_img_write(p->olabel, p->detectionname, NULL, PROGRAM_NAME);
      p->olabel->name=NULL;
    }


  /* Report the ending of initial detection. */
  if(!p->cp.quiet)
    {
      if( asprintf(&msg, "%zu initial detections found.",
                   p->numinitialdets)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
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
  int32_t *l=p->olabel->array;
  uint8_t *b=p->binary->array, *bf=b+p->binary->size;

  if(s0d1)
    /* Set all sky regions (label equal to zero) to zero, since a blank
       pixel is also non-zero, we don't need to check for blanks in this
       case. */
    do *w++ = *l++ ? *b : 0; while(++b<bf);
  else
    /* Set all detected pixels to 1. */
    do
      {
        *w++ = *l ? *l==GAL_BLANK_INT32 ? GAL_BLANK_UINT8 : 1 : *b;
        ++l;
      }
    while(++b<bf);
}





/* Copy the space of this tile into the full/large array. */
static void
detection_write_in_large(gal_data_t *tile, gal_data_t *copy)
{
  uint8_t *c=copy->array;
  GAL_TILE_PARSE_OPERATE(tile, NULL, 0, 0, {*i=*c++;});
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
  uint8_t *b, *bf;
  gal_data_t *tile, *copy, *tblock;
  size_t i, dsize[]={1,1,1,1,1,1,1,1,1,1}; /* For upto 10-Dimensions! */


  /* A temporary data structure to wrap around the copy space. Note that
     the initially allocated space for this tile is only 1 pixel! */
  copy=gal_data_alloc(NULL, GAL_TYPE_UINT8, p->input->ndim, dsize,
                      NULL, 0, -1, 1, NULL, NULL, NULL);
  free(copy->array);
  copy->array=&fho_prm->copyspace[p->maxltcontig*tprm->id];


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. */
      tile=&p->ltl.tiles[tprm->indexs[i]];

      /* Change the tile pointers (temporarily). */
      tarray=tile->array;
      tblock=tile->block;
      tile->array=gal_tile_block_relative_to_other(tile, fho_prm->workbin);
      tile->block=fho_prm->workbin;

      /* Copy the tile into the contiguous patch of memory to work on, but
         first reset the size element so 'gal_data_copy_to_allocated' knows
         there is enough space. */
      copy->flag=0;
      copy->size=p->maxltcontig;
      gal_data_copy_to_allocated(tile, copy);

      /* Take blank values to the background (set them to zero) if
         necsesary. */
      if( p->blankasforeground==0
          && gal_blank_present(p->input,0)
          && gal_blank_present(copy, 1) )
        {
          bf=(b=copy->array)+copy->size;
          do *b = *b==GAL_BLANK_UINT8 ? 0 : *b; while(++b<bf);
        }

      /* Fill the holes in this tile: holes with maximal connectivity means
         that they are most strongly bounded. */
      gal_binary_holes_fill(copy, p->holengb==4 ? 1 : 2, -1);
      if(fho_prm->step==1)
        {
          detection_write_in_large(tile, copy);
          tile->array=tarray;
          tile->block=tblock;
          continue;
        }

      /* Open all the regions. */
      gal_binary_open(copy, p->dopening, p->dopeningngb==4 ? 1 : 2, 1);

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
  float *f;
  uint8_t *b, *bf;
  gal_data_t *bin;
  int con=p->pseudoconcomp==4 ? 1 : 2;
  struct fho_params fho_prm={0, NULL, workbin, worklab, p};


  /* Set all the initial detected pixels to blank values. */
  detection_pseudo_sky_or_det(p, workbin->array, s0d1);
  if(p->detectionname)
    {
      workbin->name = s0d1 ? "DTHRESH-ON-DET" : "DTHRESH-ON-SKY";
      gal_fits_img_write(workbin, p->detectionname, NULL, PROGRAM_NAME);
      workbin->name=NULL;
    }


  /* Allocate the space necessary to work on each tile (to avoid having to
     allocate it it separately for each tile and within each
     thread. 'maxltcontig' is the maximum contiguous patch of memory needed
     to store all tiles. Finally, since we are working on a 'uint8_t' type,
     the size of each element is only 1 byte. */
  fho_prm.copyspace=gal_pointer_allocate(GAL_TYPE_UINT8,
                                         p->cp.numthreads*p->maxltcontig, 0,
                                         __func__, "fho_prm.copyspace");


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
          /* Put a copy of 'workbin' into 'bin' for every step (only
             necessary for the second step and after). For the first time
             it was already copied.*/
          if(fho_prm.step>1)
            memcpy(bin->array, workbin->array, workbin->size);

          /* Do the respective step. */
          gal_threads_spin_off(detection_fill_holes_open, &fho_prm,
                               p->ltl.tottiles, p->cp.numthreads);

          /* Reset the blank values (if they were changed). */
          if( p->blankasforeground==0 && gal_blank_present(p->input,0) )
            {
              f=p->input->array;
              bf=(b=workbin->array)+workbin->size;
              do *b = isnan(*f++) ? GAL_BLANK_UINT8 : *b; while(++b<bf);
            }

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
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so "
                    "we can address the issue. the value %d is not "
                    "recognized.", __func__, PACKAGE_BUGREPORT, fho_prm.step);
            }

          /* Write the temporary array into the check image. */
          gal_fits_img_write(bin, p->detectionname, NULL, PROGRAM_NAME);

          /* Increment the step counter. */
          ++fho_prm.step;
        }

      /* Clean up: the array in 'bin' should just be replaced with that in
         'workbin' because it is used in later steps. */
      if(workbin->mmapname)
        {
          /* Delete the memory mapped file and set the filename of 'bin'
             for 'workbin'. */
          remove(workbin->mmapname);
          free(workbin->mmapname);
          workbin->mmapname=bin->mmapname;
          bin->mmapname=NULL;
        }
      else
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
     'workbin' dataset when working on the Sky. Recall that in this case,
     the blank pixels are the detections. On the Sky image, blank should be
     set to 1 (because we want the detected objects to have the same labels
     as the pseudo-detections that cover them). This will allow us to later
     remove these pseudo-detections.
  if(s0d1==0)
    {
      bf=(b=workbin->array)+workbin->size;
      do if(*b==GAL_BLANK_UINT8) *b = !s0d1; while(++b<bf);
    }
  */
  return gal_binary_connected_components(workbin, &worklab, con);
}





/* Write the S/N tables to a file. */
static void
detection_sn_write_to_file(struct noisechiselparams *p, gal_data_t *sn,
                           gal_data_t *snind, int s0d1D2)
{
  char *str, *extname=NULL;
  gal_list_str_t *comments=NULL;

  /* Comment for extension on further explanation. */
  if( asprintf(&str, "See also: '%s' HDU of output with "
               "'--checkdetection'", ( s0d1D2<2
                                       ? "PSEUDOS-FOR-SN": "DILATED" ))<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_str_add(&comments, str, 0);


  /* Description comment. */
  str = ( s0d1D2
          ? ( s0d1D2==2
              ? "S/N of grown detections."
              : "Pseudo-detection S/N over initial detections." )
          : "Pseudo-detection S/N over initial undetections.");
  gal_list_str_add(&comments, str, 1);


  /* Set the file name and write the table. */
  str = ( s0d1D2
          ? ( s0d1D2==2 ? p->detsn_D_name : p->detsn_d_name )
          : p->detsn_s_name );
  if( p->cp.tableformat!=GAL_TABLE_FORMAT_TXT )
    extname = ( s0d1D2
                ? ( s0d1D2==2 ? "GROWN_DETECTION_SN" : "DET_PSEUDODET_SN" )
                : "SKY_PSEUDODET_SN" );
  threshold_write_sn_table(p, sn, snind, str, comments, extname);
  gal_list_str_free(comments, 1);
}






static gal_data_t *
detection_sn(struct noisechiselparams *p, gal_data_t *worklab, size_t num,
             int s0d1D2, char *extname)
{
  float *snarr;
  uint8_t *flag;
  size_t tablen=num+1;
  gal_data_t *sn, *snind;
  int32_t *plabend, *indarr=NULL;
  double ave, err, *xy, *brightness;
  float s, ss, *f, *ff, *fs, *sky=p->sky->array;
  size_t ind, ndim=p->input->ndim, xyncols=1+ndim;
  size_t i, *area, counter=0, *dsize=p->input->dsize;
  int32_t *plab = worklab->array, *dlab = s0d1D2 ? NULL : p->olabel->array;
  size_t *coord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                      "coord");

  /* Sanity check. */
  if(p->input->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: the input dataset must be float32 type, "
          "it is %s", __func__, gal_type_name(p->input->type, 1));
  if(!isnan(GAL_BLANK_FLOAT32))
    error(EXIT_FAILURE, 0, "%s: only a NaN value is recognized for blank "
          "floating point data types, the blank value is defined to be %f",
          __func__, GAL_BLANK_FLOAT32);
  if(ndim!=2)
    error(EXIT_FAILURE, 0, "%s: only 2D datasets are acceptable, your input "
          "is %zu dimensions", __func__, ndim);


  /* Allocate all the necessary arrays, note that since we want to put each
     object's information into the same index, the number of allocated
     spaces has to be 'tablen=num+1'. */
  area       = gal_pointer_allocate(GAL_TYPE_SIZE_T,  tablen, 1, __func__,
                                    "area");
  brightness = gal_pointer_allocate(GAL_TYPE_FLOAT64, tablen, 1, __func__,
                                     "brightness");
  xy         = gal_pointer_allocate(GAL_TYPE_FLOAT64, xyncols*tablen, 1,
                                     __func__, "xy");
  flag       = ( s0d1D2==0
                 ? gal_pointer_allocate(GAL_TYPE_UINT8, tablen, 1, __func__,
                                         "flag")
                 : NULL );
  sn         = gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &tablen, NULL, 1,
                              p->cp.minmapsize, p->cp.quietmmap,
                              "SIGNAL-TO-NOISE", "ratio", NULL);
  snind      = ( p->checksn==0 ? NULL
                 : gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &tablen, NULL, 1,
                                  p->cp.minmapsize, p->cp.quietmmap, "LABEL",
                                  "counter", NULL) );

  /* Go over all the pixels and get the necessary information. */
  fs = f = p->input->array;
  ff = f + p->input->size;
  do
    {
      /* All this work is only necessary when we are actually on a
         pseudo-detection label: it is non-zero and not blank. */
      if(*plab && ( (p->input->flag | GAL_DATA_FLAG_HASBLANK)
                    && *plab!=GAL_BLANK_INT32 ) )
        {
          /* For Sky pseudo-detections we'll start to see if it has already
             been determined that the object lies over a detected object or
             not. If it does, then just ignore it. */
          if(s0d1D2==0)
            {
              if( flag[*plab] ) { ++plab; ++dlab; continue; }
              else if(*dlab)    /* We are on a detection. */
                { flag[*plab]=1; area[*plab]=0; ++plab; ++dlab; continue; }
            }

          /* If we are on a blank pixel, ignore this pixel. */
          if( isnan(*f) ) { ++plab; if(s0d1D2==0) ++dlab; continue; }

          /* Save all the necessary values. */
          ++area[*plab];
          gal_dimension_index_to_coord(f-fs, ndim, dsize, coord);
          s  = sky[ gal_tile_full_id_from_coord(&p->cp.tl, coord) ];
          ss = *f-s;
          brightness[*plab] += ss;
          if( ss > 0.0f )  /* For calculatiing the approximate center, */
            {              /* necessary for calculating Sky and STD.   */
              xy[*plab*xyncols  ] += ss;
              xy[*plab*xyncols+1] += (double)coord[0] * ss;
              xy[*plab*xyncols+2] += (double)coord[1] * ss;
            }
        }

      /* Increment the other two labels. */
      ++plab;
      if(s0d1D2==0) ++dlab;
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


  /* If the user wants to see the steps (on the background) and we are
     working on pseudo-detections, remove those that will not be used in
     the final quantile calculation. */
  if(p->detectionname)
    {
      if(s0d1D2<2)
        {
          plabend = (plab=worklab->array) + worklab->size;
          do
            if( *plab!=GAL_BLANK_INT32
                && ( area[*plab]<p->snminarea || brightness[*plab]<0) )
              *plab=0;
          while(++plab<plabend);
        }
      worklab->name=extname;
      gal_fits_img_write(worklab, p->detectionname, NULL, PROGRAM_NAME);
      worklab->name=NULL;
    }


  /* Calculate the signal to noise for successful detections: */
  snarr=sn->array;
  if(snind) indarr=snind->array;
  if(s0d1D2) { snarr[0]=NAN; if(snind) indarr[0]=GAL_BLANK_INT32; }
  for(i=1;i<tablen;++i)
    {
      ave=brightness[i]/area[i];
      if( area[i]>p->snminarea && ave>0.0f && xy[i*xyncols]>0.0f )
        {
          /* Get the flux weighted center coordinates. */
          coord[0]=GAL_DIMENSION_FLT_TO_INT( xy[i*xyncols+1]/xy[i*xyncols] );
          coord[1]=GAL_DIMENSION_FLT_TO_INT( xy[i*xyncols+2]/xy[i*xyncols] );

          /* Get the Sky standard deviation on this tile. */
          err  = ((float *)(p->std->array))[
                         gal_tile_full_id_from_coord(&p->cp.tl, coord) ];

          /* Correct the index in the sn to store the Signal to noise
             ratio. When we are dealing with the noise, we only want the
             non-zero signal to noise values, so we will just use a
             counter. But for initial detections, it is very important that
             their Signal to noise ratio be placed in the same index as
             their label. */
          ind = s0d1D2 ? i : counter++;
          if(snind) indarr[ind]=i;
          snarr[ind] = ( sqrt( (float)(area[i])/p->cpscorr )
                         * ave / sqrt( ave + err*err ) );
        }
      else
        /* In detection pseudo-detections, order matters, so we will set
           all non-usable values to blank. */
        if(s0d1D2)
          {
            snarr[i]=NAN;
            if(snind) indarr[i]=GAL_BLANK_INT32;;
          }
    }


  /* A small sanity check. */
  if( s0d1D2==0 && counter==0 )
    error(EXIT_FAILURE, 0, "no sky pseudo-detections.");


  /* If we are in Sky mode, the sizes have to be corrected */
  if(s0d1D2==0)
    {
      sn->dsize[0]=sn->size=counter;
      if(snind) snind->dsize[0]=snind->size=counter;
    }


  /* If the user wanted a list of S/N values for all pseudo-detections,
     save it. */
  if(snind) detection_sn_write_to_file(p, sn, snind, s0d1D2);


  /* Clean up and return. */
  free(xy);
  free(area);
  free(coord);
  free(brightness);
  if(flag) free(flag);
  if(snind) gal_data_free(snind);
  return sn;
}





/* ONLY FOR PSEUDO DETECTIONS: remove pseudo-detections that have a small
   S/N from the binary image (the labeled image will be left untouched). */
static void
detection_pseudo_remove_low_sn(struct noisechiselparams *p,
                               gal_data_t *workbin, gal_data_t *worklab,
                               gal_data_t *sn)
{
  size_t i;
  float *snarr=sn->array;
  uint8_t *b=workbin->array;
  int32_t *l=worklab->array, *lf=l+worklab->size;
  uint8_t *keep=gal_pointer_allocate(GAL_TYPE_UINT8, sn->size, 1, __func__,
                                     "keep");

  /* Specify the new labels for those that must be kept/changed. Note that
     when an object didn't have an S/N, its S/N was given a value of NaN
     (which will fail on any condition), so it acts as if it had an S/N
     lower than the required value. */
  for(i=0;i<sn->size;++i)
    if( snarr[i] > p->detsnthresh ) keep[i]=1;


  /* Go over the pseudo-detection labels and only keep those that must be
     kept (using the new labels) in the binary array. */
  if( p->input->flag & GAL_DATA_FLAG_HASBLANK )
    do
      *b++ = *l == GAL_BLANK_INT32 ? GAL_BLANK_UINT8 : keep[ *l ] > 0;
    while(++l<lf);
  else
    do *b++ = keep[ *l ] > 0; while(++l<lf);


  /* If the user wanted to see the steps. */
  if(p->detectionname)
    {
      workbin->name="TRUE-PSEUDOS";
      gal_fits_img_write(workbin, p->detectionname, NULL,
                         PROGRAM_NAME);
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
  size_t numpseudo;
  struct timeval t1;
  gal_data_t *sn, *quant, *workbin, *worklab;

  /* Allocate the space for the working datasets. */
  worklab=gal_data_copy(p->olabel);
  workbin=gal_data_alloc(NULL, GAL_TYPE_UINT8, p->input->ndim,
                         p->input->dsize, p->input->wcs, 0, p->cp.minmapsize,
                         p->cp.quietmmap, NULL, NULL, NULL);
  workbin->flag=p->input->flag;


  /* If a manual S/N is given, then don't bother with using
     pseudo-detections. */
  if( isnan(p->snthresh) )
    {
      /* Over the Sky: find the pseudo-detections and make the S/N
         table. */
      if(!p->cp.quiet) gettimeofday(&t1, NULL);
      numpseudo=detection_pseudo_find(p, workbin, worklab, 0);
      sn=detection_sn(p, worklab, numpseudo, 0, "PSEUDOS-FOR-SN");


      /* A small sanity check */
      if( sn->size < p->minnumfalse)
        error(EXIT_FAILURE, 0, "only %zu pseudo-detections could be found "
              "over the sky region to estimate an S/N. This is less than "
              "%zu (value to '--minnumfalse' option). Please adjust "
              "parameters like '--dthresh', '--snminarea', or make sure "
              "that there actually is sufficient sky area after initial "
              "detection. You can use '--checkdetection' to see every step "
              "until this point", sn->size, p->minnumfalse);


      /* Get the S/N quantile and report it if we are in non-quiet mode. */
      quant=gal_statistics_quantile(sn, p->snquant, 1);
      p->detsnthresh = *((float *)(quant->array));
      if(!p->cp.quiet)
        {
          if( asprintf(&msg, "Pseudo-det S/N: %.3f (%g quant of %zu).",
                       p->detsnthresh, p->snquant, sn->size)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_timing_report(&t1, msg, 2);
          free(msg);
        }
      gal_data_free(sn);
      gal_data_free(quant);
    }
  else
    p->detsnthresh=p->snthresh;


  /* Over the detections: find pseudo-detections and make S/N table. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  numpseudo=detection_pseudo_find(p, workbin, worklab, 1);
  sn=detection_sn(p, worklab, numpseudo, 1, "PSEUDOS-FOR-SN");


  /* Remove the pseudo detections with a low S/N. */
  detection_pseudo_remove_low_sn(p, workbin, worklab, sn);


  /* Clean up and return. */
  gal_data_free(sn);
  gal_data_free(worklab);
  return workbin;
}





/* This is for the final detections (grown) detections. */
static size_t
detection_final_remove_small_sn(struct noisechiselparams *p,
                                gal_data_t *workbin, size_t num)
{
  size_t i;
  int8_t *b;
  float *snarr;
  gal_data_t *sn;
  int32_t *l, *lf, curlab=1;
  int32_t *newlabs=gal_pointer_allocate(GAL_TYPE_INT32, num+1, 1, __func__,
                                        "newlabs");

  /* Get the Signal to noise ratio of all detections. */
  sn=detection_sn(p, p->olabel, num, 2, "DILATED");

  /* Only keep the objects with an S/N above the pseudo-detection limit. */
  snarr=sn->array;
  for(i=1;i<num+1;++i)
    newlabs[i] = snarr[i] > p->detsnthresh ? curlab++ : 0;

  /* Go over the labeled image and correct the labels. */
  b=workbin->array;
  lf=(l=p->olabel->array)+p->olabel->size;
  if( p->input->flag & GAL_DATA_FLAG_HASBLANK )
    {
      do
        {
          if( *l != GAL_BLANK_INT32 )
            *b = (*l=newlabs[ *l ]) > 0;
          ++b;
        }
      while(++l<lf);
    }
  else do *b++ = (*l=newlabs[ *l ]) > 0; while(++l<lf);

  /* Clean up and return. */
  free(newlabs);
  gal_data_free(sn);
  return curlab-1;
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
  float *e_th, *arr=p->conv->array;
  int32_t *l=p->olabel->array, *lf=l+p->olabel->size, curlab=1;
  int32_t *newlabels=gal_pointer_allocate(GAL_TYPE_UINT32,
                                          p->numinitialdets+1, 1, __func__,
                                          "newlabels");

  /* Find the new labels for all the existing labels. Recall that
     'newlabels' was initialized to zero, so any label that is not given a
     new label here will be automatically removed. After the first pixel of
     a label overlaps with dbyt[i], we don't need to check the rest of that
     object's pixels. At this point, tokeep is only binary: 0 or 1.

     Note that the zeroth element of 'tokeep' can also be non zero, this is
     because the holes of the labeled regions might be filled during
     filling the holes, but have not been filled in the original labeled
     array. They are not important so you can just ignore them. */
  do
    {
      if( *l && *l!=GAL_BLANK_INT32 )
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
  for(i=0;i<p->numinitialdets;++i) if(newlabels[i]) newlabels[i] = curlab++;


  /* Replace the byt and olab values with their proper values. Note that we
     need the binary array when there is no growth also. The binary array
     is later used in estimating the sky and its STD.

     See if growth is necessary or not. Note that even if the user has
     asked for growth, if the edges of the objects in the image are sharp
     enough, no growth will be necessary (and thus the labeled image won't
     be re-written during growth). So it is necessary to check for growth
     here and later do it in 'detection_quantile_expand'. */
  p->numexpand=0;
  b=workbin->array;
  l=p->olabel->array;
  if(p->detgrowquant==1.0)
    do
      {                                       /* No growth will happen.  */
        if(*l!=GAL_BLANK_INT32)               /* So set both the binary  */
          *b = ( *l = newlabels[ *l ] ) > 0;  /* AND label images.       */
        ++b;
      }
    while(++l<lf);
  else
    {
      /* Expand the threshold values (from one value for each tile) to the
         whole dataset. Since we know the tiles cover the whole image, we
         don't neeed to initialize or check for blanks. */
      p->exp_thresh_full=gal_tile_block_write_const_value(p->expand_thresh,
                                                          p->cp.tl.tiles, 0,
                                                          0);

      /* Remove the false detections and count how many pixels need to
         grow. */
      e_th=p->exp_thresh_full->array;
      do                                    /* Growth is necessary later.  */
        {                                   /* So there is no need to set  */
          if(*l==GAL_BLANK_INT32)           /* the labels image, but we    */
            *b=GAL_BLANK_UINT8;             /* have to count the number of */
          else                              /* pixels to (possibly) grow.  */
            {
              *b = newlabels[ *l ] > 0;
              if( *b==0 && *arr>*e_th )
                ++p->numexpand;
            }
          ++b; ++arr; ++e_th;
        }
      while(++l<lf);


      /* If there aren't any pixels to later expand, then reset the labels
         (remove false detections in the labeled image). */
      if(p->numexpand==0)
        {
          l=p->olabel->array;
          do if(*l!=GAL_BLANK_INT32) *l = newlabels[ *l ]; while(++l<lf);
        }
    }


  /* Clean up and return. */
  free(newlabels);
  return curlab-1;
}





/* Expand the initial detections based on the quantile threshold and then
   label the connected regions. If expansion is not possible, then return
   the 'GAL_BLANK_SIZET'.*/
static size_t
detection_quantile_expand(struct noisechiselparams *p, gal_data_t *workbin)
{
  int32_t *o, *of;
  size_t *d, numexpanded=0;
  gal_data_t *diffuseindexs;
  float *i, *e_th, *arr=p->conv->array;
  uint8_t *b=workbin->array, *bf=b+workbin->size;

  /* Only continue if there actually are any pixels to expand (this can
     happen!). */
  if(p->numexpand)
    {
      /* Allocate the space necessary to keep the index of all the pixels
         that must be expanded and re-initialize the necessary pointers. */
      diffuseindexs=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &p->numexpand,
                                   NULL, 0, p->cp.minmapsize, p->cp.quietmmap,
                                   NULL, NULL, NULL);

      /* Fill in the diffuse indexs and initialize the objects dataset. */
      b    = workbin->array;
      arr  = p->conv->array;
      d    = diffuseindexs->array;
      e_th = p->exp_thresh_full->array;
      of   = (o=p->olabel->array) + p->olabel->size;
      do
        {
          /* If the binary value is 1, then we want an initial label of 1
             (the object is already detected). If it isn't, then we only
             want it if it is above the threshold. */
          *o = *b==1 ? 1 : ( *arr>*e_th ? GAL_LABEL_INIT : 0);
          if(*b==0 && *arr>*e_th)
            *d++ = o - (int32_t *)(p->olabel->array);

          /* Increment the points and go onto the next pixel. */
          ++b;
          ++arr;
          ++e_th;
        }
      while(++o<of);

      /* Expand the detections. Note that because we are only concerned
         with those regions that are touching a detected region, it is
         irrelevant to sort the dataset. */
      gal_label_grow_indexs(p->olabel, diffuseindexs, 0, p->olabel->ndim);

      /* Only keep the 1 valued pixels in the binary array and fill its
         holes. */
      o=p->olabel->array;
      bf=(b=workbin->array)+workbin->size;
      do *b = (*o++ == 1); while(++b<bf);
      workbin=gal_binary_dilate(workbin, 1, 1, 1);
      gal_binary_holes_fill(workbin, 1, p->detgrowmaxholesize);

      /* Get the labeled image. */
      numexpanded=gal_binary_connected_components(workbin, &p->olabel,
                                                  workbin->ndim);

      /* Set all the input's blank pixels to blank in the labeled and
         binary arrays. */
      if( gal_blank_present(p->input, 1) )
        {
          b=workbin->array;
          i=p->input->array;
          of=(o=p->olabel->array)+p->olabel->size;
          do
            {
              if(isnan(*i++))
                {
                  *o=GAL_BLANK_INT32;
                  *b=GAL_BLANK_UINT8;
                }
              ++b;
            }
          while(++o<of);
        }

      /* Clean up. */
      gal_data_free(diffuseindexs);
    }

  /* Clean up and return */
  gal_data_free(p->expand_thresh);
  gal_data_free(p->exp_thresh_full);
  return numexpanded ? numexpanded : GAL_BLANK_SIZE_T;
}





/* The initial detection has been done, now we want to remove false
   detections. */
void
detection(struct noisechiselparams *p)
{
  char *msg;
  gal_data_t *workbin;
  struct timeval t0, t1;
  size_t num_true_initial, num_expanded;

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
      if( asprintf(&msg, "Pseudo-detection thresh (%.3f sigma) applied.",
                   p->dthresh)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* Find the real pseudo-detections. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  workbin=detection_pseudo_real(p);


  /* Only keep the initial detections that overlap with the real
     pseudo-detections. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  num_true_initial=detection_remove_false_initial(p, workbin);
  if(p->detectionname)
    {
      workbin->name="DETECTIONS-INIT-TRUE";
      gal_fits_img_write(workbin, p->detectionname, NULL,
                         PROGRAM_NAME);
      workbin->name=NULL;
    }
  if(!p->cp.quiet)
    {
      if( asprintf(&msg, "%zu false initial detections removed.",
                   p->numinitialdets - num_true_initial)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* If the user asked for dilation/expansion, then apply it and report the
     final number of detections. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  if(p->detgrowquant!=1.0f)
    {
      num_expanded=detection_quantile_expand(p, workbin);
      if(num_expanded!=GAL_BLANK_SIZE_T)
        num_true_initial=num_expanded;
    }


  /* Update the user on the progress, if necessary. */
  if(!p->cp.quiet && p->detgrowquant!=1.0f && num_expanded!=GAL_BLANK_SIZE_T )
    {
      /* If the user hasn't asked for a labeled image, then don't confuse
         them with the number of detections, just let them know that growth
         is complete. */
      if(p->label)
        {
          if( asprintf(&msg, "%zu detections after growth to %.3f "
                       "quantile.", num_true_initial, p->detgrowquant)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        {
          if( asprintf(&msg, "Growth to %.3f quantile complete.",
                       p->detgrowquant)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* When the final (grown or over-all object) detection's S/N is less than
     the pseudo-detection's S/N limit, the object is false. For a real
     detection, the actual object S/N should be higher than any of its
     pseudo-detection because it has a much larger area (and possibly more
     flux under it).  So when the final S/N is smaller than the minimum
     acceptable S/N threshold, we have a false pseudo-detection. */
  p->numdetections = ( p->cleangrowndet
                       ?  detection_final_remove_small_sn(p, workbin,
                                                          num_true_initial)
                       : num_true_initial );
  if(!p->cp.quiet)
    {
      /* If the user hasn't asked for a labeled image, then don't confuse
         them with the number of detections, just let them know that growth
         is complete. */
      if(p->label)
        {
          if( asprintf(&msg, "%zu final true detections.",
                       p->numdetections)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        gal_checkset_allocate_copy("Detection complete.", &msg);
      gal_timing_report(&t0, msg, 1);
      free(msg);
    }
  if(p->detectionname)
    {
      p->olabel->name="DETECTION-FINAL";
      gal_fits_img_write(p->olabel, p->detectionname, NULL,
                         PROGRAM_NAME);
      p->olabel->name=NULL;
    }


  /* p->binary was used to keep the initial pseudo-detection threshold. But
     we don't need it any more, so we'll just free it and put the 'workbin'
     array in its place. Note that 'workbin' has a map of all the detected
     objects, which is still necessary during NoiseChisel. */
  gal_data_free(p->binary);
  p->binary=workbin;


  /* The initial Sky and Sky STD values were only for detection. */
  gal_data_free(p->sky);
  gal_data_free(p->std);
  p->sky = p->std = NULL;


  /* If the user wanted to check the threshold and hasn't called
     'continueaftercheck', then stop NoiseChisel. */
  if(p->detectionname && !p->continueaftercheck)
    ui_abort_after_check(p, p->detectionname, NULL,
                         "showing all detection steps");
}
