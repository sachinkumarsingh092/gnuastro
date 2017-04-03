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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/statistics.h>
#include <gnuastro/interpolate.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "ui.h"
#include "threshold.h"










/**********************************************************************/
/***************        Apply a given threshold.      *****************/
/**********************************************************************/
void
threshold_apply(struct noisechiselparams *p, float *value1, float *value2,
                int type)
{
  size_t tid;
  void *tarray=NULL;
  gal_data_t *tile, *tblock=NULL;

  /* Clear the binary array (this is mainly because the input may contain
     blank values and we won't be doing the thresholding no those
     pixels. */
  memset(p->binary->array, 0, p->binary->size);

  /* Go over all the tiles. */
  for(tid=0; tid<p->cp.tl.tottiles; ++tid)
    {
      /* For easy reading. */
      tile=&p->cp.tl.tiles[tid];

      switch(type)
        {

        /* This is a quantile threshold. */
        case THRESHOLD_QUANTILES:
          /* Correct the tile's pointers to apply the threshold on the
             convolved image. */
          if(p->conv)
            {
              tarray=tile->array; tblock=tile->block;
              tile->array=gal_tile_block_relative_to_other(tile, p->conv);
              tile->block=p->conv;
            }

          /* Apply the threshold. */
          GAL_TILE_PARSE_OPERATE({
              *o = ( *i > value1[tid]
                     ? ( *i > value2[tid] ? THRESHOLD_NO_ERODE_VALUE : 1 )
                     : 0 );
            }, tile, p->binary, 1, 1);

          /* Revert the tile's pointers back to what they were. */
          if(p->conv) { tile->array=tarray; tile->block=tblock; }
          break;


        /* This is a Sky and Sky STD threshold. */
        case THRESHOLD_SKY_STD:

          /* The threshold is always low. So for the majority of non-NaN
             pixels in the image, the condition above will be true. If we
             come over a NaN pixel, then by definition of NaN, all
             conditionals will fail.

             If an image doesn't have any NaN pixels, only the pixels below
             the threshold have to be checked for a NaN which are by
             definition a very small fraction of the total pixels. And if
             there are NaN pixels in the image. */
          GAL_TILE_PARSE_OPERATE({
              *o = ( ( *i - value1[tid] > p->dthresh * value2[tid] )
                     ? 1 : *i==*i ? 0 : GAL_BLANK_UINT8 );
            }, tile, p->binary, 1, 1);
          break;


        default:
          error(EXIT_FAILURE, 0, "a bug! please contact us at %s so we can "
                "address the problem. For some reason a value of %d had "
                "been given to `type' in `threshold_apply'",
                PACKAGE_BUGREPORT, type);
        }
    }
}




















/**********************************************************************/
/***************            Write S/N values          *****************/
/**********************************************************************/
void
threshold_write_sn_table(struct noisechiselparams *p, gal_data_t *insn,
                         gal_data_t *inind, char *filename,
                         struct gal_linkedlist_stll *comments)
{
  gal_data_t *sn, *ind, *cols;

  /* Remove all blank elements. The index and sn values must have the same
     set of blank elements, but checking on the integer array is faster. */
  if( gal_blank_present(inind) )
    {
      ind=gal_data_copy(inind);
      sn=gal_data_copy(insn);
      gal_blank_remove(ind);
      gal_blank_remove(sn);
    }
  else
    {
      sn  = insn;
      ind = inind;
    }

  /* Set the columns. */
  cols       = ind;
  cols->next = sn;


  /* Prepare the comments. */
  gal_table_comments_add_intro(&comments, PROGRAM_STRING, &p->rawtime);


  /* write the table. */
  gal_table_write(cols, comments, p->cp.tableformat, filename, 1);


  /* Clean up (if necessary). */
  if(sn!=insn) gal_data_free(sn);
  if(ind==inind) ind->next=NULL; else gal_data_free(ind);
}



















/**********************************************************************/
/***************     Interpolation and smoothing      *****************/
/**********************************************************************/
/* Interpolate and smooth the values for each tile over the whole image. */
static void
threshold_interp_smooth(struct noisechiselparams *p, gal_data_t **first,
                        gal_data_t **second, char *filename)
{
  gal_data_t *tmp;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;

  /* A small sanity check. */
  if( (*first)->next )
    error(EXIT_FAILURE, 0, "the `first' argument to "
          "`threshold_interp_smooth' must not have any `next' pointer.");


  /* Do the interpolation of both arrays. */
  (*first)->next = *second;
  tmp=gal_interpolate_close_neighbors(*first, tl, cp->interpnumngb,
                                      cp->numthreads, cp->interponlyblank, 1);
  gal_data_free(*first);
  gal_data_free(*second);
  *first=tmp;
  *second=tmp->next;
  (*first)->next=(*second)->next=NULL;
  if(filename)
    {
      gal_tile_full_values_write(*first, tl, filename, PROGRAM_STRING);
      gal_tile_full_values_write(*second, tl, filename, PROGRAM_STRING);
    }


  /* Smooth the threshold if requested. */
  if(p->smoothwidth>1)
    {
      /* Smooth the first. */
      tmp=gal_tile_full_values_smooth(*first, tl, p->smoothwidth,
                                      p->cp.numthreads);
      gal_data_free(*first);
      *first=tmp;

      /* Smooth the second */
      tmp=gal_tile_full_values_smooth(*second, tl, p->smoothwidth,
                                      p->cp.numthreads);
      gal_data_free(*second);
      *second=tmp;

      /* Add them to the check image. */
      if(filename)
        {
          gal_tile_full_values_write(*first, tl, filename, PROGRAM_STRING);
          gal_tile_full_values_write(*second, tl, filename, PROGRAM_STRING);
        }
    }
}




















/****************************************************************
 ************           Quantile threshold           ************
 ****************************************************************/
struct qthreshparams
{
  gal_data_t        *erode_th;
  gal_data_t      *noerode_th;
  void                 *usage;
  struct noisechiselparams *p;
};





static void *
qthresh_on_tile(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct qthreshparams *qprm=(struct qthreshparams *)tprm->params;
  struct noisechiselparams *p=qprm->p;

  double *darr;
  void *tarray=NULL;
  int type=qprm->erode_th->type;
  gal_data_t *tile, *mode, *qvalue, *usage, *tblock=NULL;
  size_t i, tind, twidth=gal_type_sizeof(type), ndim=p->input->ndim;

  /* Put the temporary usage space for this thread into a data set for easy
     processing. */
  usage=gal_data_alloc(gal_data_ptr_increment(qprm->usage,
                                              tprm->id*p->maxtcontig, type),
                       type, ndim, p->maxtsize, NULL, 0, p->cp.minmapsize,
                       NULL, NULL, NULL);

  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    {
      /* Re-initialize the usage array's space (will be changed in
         `gal_data_copy_to_allocated' for each tile). */
      usage->ndim=ndim;
      usage->size=p->maxtcontig;
      memcpy(usage->dsize, p->maxtsize, ndim*sizeof *p->maxtsize);


      /* For easy reading. */
      tind = tprm->indexs[i];
      tile = &p->cp.tl.tiles[tind];


      /* If we have a convolved image, temporarily change the tile's
         pointers so we can do the work on the convolved image, then copy
         the desired contents into the already allocated `usage' array. */
      if(p->conv)
        {
          tarray=tile->array; tblock=tile->block;
          tile->array=gal_tile_block_relative_to_other(tile, p->conv);
          tile->block=p->conv;
        }
      gal_data_copy_to_allocated(tile, usage);
      if(p->conv) { tile->array=tarray; tile->block=tblock; }


      /* Find the mode on this dataset, note that we have set the `inplace'
         flag to `1'. So as a byproduct of finding the mode, `usage' is not
         going to have any blank elements and will be sorted (thus ready to
         be used by the quantile functions). */
      mode=gal_statistics_mode(usage, p->mirrordist, 1);


      /* Check the mode value. Note that if the mode is not accurate, then
         the contents of `darr' will be NaN and all conditions will
         fail. In such cases, the tile will be ignored. */
      darr=mode->array;
      if( fabs(darr[1]-0.5f) < p->modmedqdiff )
        {
          /* Get the erosion quantile for this tile and save it. Note that
             the type of `qvalue' is the same as the input dataset. */
          qvalue=gal_statistics_quantile(usage, p->qthresh, 1);
          memcpy(gal_data_ptr_increment(qprm->erode_th->array, tind, type),
                 qvalue->array, twidth);
          gal_data_free(qvalue);

          /* Do the same for the no-erode quantile. */
          qvalue=gal_statistics_quantile(usage, p->noerodequant, 1);
          memcpy(gal_data_ptr_increment(qprm->noerode_th->array, tind, type),
                 qvalue->array, twidth);
          gal_data_free(qvalue);
        }
      else
        {
          gal_blank_write(gal_data_ptr_increment(qprm->erode_th->array,
                                                 tind, type), type);
          gal_blank_write(gal_data_ptr_increment(qprm->noerode_th->array,
                                                 tind, type), type);
        }

      /* Clean up and fix the tile's pointers. */
      gal_data_free(mode);
    }

  /* Clean up and wait for the other threads to finish, then return. */
  usage->array=NULL;  /* Not allocated here. */
  gal_data_free(usage);
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





void
threshold_quantile_find_apply(struct noisechiselparams *p)
{
  char *msg;
  struct timeval t1;
  struct qthreshparams qprm;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;


  /* Get the starting time if necessary. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);


  /* Add image to check image if requested. If the user has asked for
     `oneelempertile', then the size of values is not going to be the same
     as the input, making it hard to inspect visually. So we'll only put
     the full input when `oneelempertile' isn't requested. */
  if(p->qthreshname && !tl->oneelempertile)
    gal_fits_img_write(p->conv ? p->conv : p->input, p->qthreshname, NULL,
                       PROGRAM_STRING);


  /* Allocate space for the quantile threshold values. */
  qprm.erode_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                               tl->numtiles, NULL, 0, cp->minmapsize,
                               "QTHRESH-ERODE", p->input->unit, NULL);
  qprm.noerode_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                                 tl->numtiles, NULL, 0, cp->minmapsize,
                                 "QTHRESH-NOERODE", p->input->unit, NULL);


  /* Allocate temporary space for processing in each tile. */
  qprm.usage=gal_data_malloc_array(p->input->type,
                                   cp->numthreads * p->maxtcontig);

  /* Find the threshold on each tile, then clean up the temporary space. */
  qprm.p=p;
  gal_threads_spin_off(qthresh_on_tile, &qprm, tl->tottiles, cp->numthreads);
  if(p->qthreshname)
    {
      gal_tile_full_values_write(qprm.erode_th, tl, p->qthreshname,
                                 PROGRAM_STRING);
      gal_tile_full_values_write(qprm.noerode_th, tl, p->qthreshname,
                                 PROGRAM_STRING);
    }
  free(qprm.usage);


  /* Interpolate and smooth the derived values. */
  threshold_interp_smooth(p, &qprm.erode_th, &qprm.noerode_th,
                          p->qthreshname);


  /* We now have a threshold for all tiles, apply it. */
  threshold_apply(p, qprm.erode_th->array, qprm.noerode_th->array,
                  THRESHOLD_QUANTILES);


  /* Write the binary image if check is requested. */
  if(p->qthreshname && !tl->oneelempertile)
    gal_fits_img_write(p->binary, p->qthreshname, NULL, PROGRAM_STRING);


  /* Clean up and report duration if necessary. */
  gal_data_free(qprm.erode_th);
  gal_data_free(qprm.noerode_th);
  if(!p->cp.quiet)
    {
      asprintf(&msg, "%.2f & %0.2f quantile thresholds applied.",
               p->qthresh, p->noerodequant);
      gal_timing_report(&t1, msg, 2);
    }


  /* If the user wanted to check the threshold and hasn't called
     `continueaftercheck', then stop NoiseChisel. */
  if(p->qthreshname && !p->continueaftercheck)
    ui_abort_after_check(p, p->qthreshname, "quantile threshold check");
}




















/****************************************************************
 ************       Mean and STD of undetected       ************
 ****************************************************************/
static void *
threshold_mean_std_undetected(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct noisechiselparams *p=(struct noisechiselparams *)tprm->params;

  double *darr, s, s2;
  int type=p->sky->type;
  size_t i, tind, numsky, dsize=2;
  gal_data_t *tile, *meanstd_d, *meanstd, *bintile;


  /* A dataset to keep the mean and STD in double type. */
  meanstd_d=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &dsize,
                           NULL, 0, -1, NULL, NULL, NULL);
  darr=meanstd_d->array;


  /* An empty dataset to replicate a tile on the binary array. */
  bintile=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &dsize,
                         NULL, 0, -1, NULL, NULL, NULL);
  free(bintile->array);
  free(bintile->dsize);
  bintile->block=p->binary;
  bintile->ndim=p->binary->ndim;


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    {
      /* Basic definitions */
      numsky=0;
      tind = tprm->indexs[i];
      tile = &p->cp.tl.tiles[tind];

      /* Correct the fake binary tile's properties to be the same as this
         one, then count the number of zero valued elements in it. */
      bintile->size=tile->size;
      bintile->dsize=tile->dsize;
      bintile->array=gal_tile_block_relative_to_other(tile, p->binary);
      GAL_TILE_PARSE_OPERATE({if(!*o) numsky++;}, tile, bintile, 1, 1);

      /* Only continue, if the fraction of Sky values are less than the
         requested fraction. */
      if( (float)(numsky)/(float)(tile->size) > p->minbfrac)
        {
          /* Calculate the mean and STD over this tile. */
          s=s2=0.0f;
          GAL_TILE_PARSE_OPERATE(
                                 {
                                   if(!*o)
                                     {
                                       s  += *i;
                                       s2 += *i * *i;
                                     }
                                 }, tile, bintile, 1, 1);
          darr[0]=s/numsky;
          darr[1]=sqrt( (s2-s*s/numsky)/numsky );

          /* Convert the mean and std into the same type as the sky and std
             arrays. */
          meanstd=gal_data_copy_to_new_type(meanstd_d, type);

          /* Copy the mean and STD to their respective places in the tile
             arrays. */
          memcpy(gal_data_ptr_increment(p->sky->array, tind, type),
                 meanstd->array, gal_type_sizeof(type));
          memcpy(gal_data_ptr_increment(p->std->array, tind, type),
                 gal_data_ptr_increment(meanstd->array, 1, type),
                 gal_type_sizeof(type));
        }
      else
        {
          gal_blank_write(gal_data_ptr_increment(p->sky->array, tind, type),
                          type);
          gal_blank_write(gal_data_ptr_increment(p->std->array, tind, type),
                          type);
        }
    }

  /* Clean up and wait for other threads to finish and abort. */
  bintile->array=NULL;
  bintile->dsize=NULL;
  gal_data_free(bintile);
  gal_data_free(meanstd_d);
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





void
threshold_sky_and_std(struct noisechiselparams *p)
{
  gal_data_t *tmp;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;


  /* When the check image has the same resolution as the input, write the
     binary array as a reference to help in the comparison. */
  if(p->detskyname && !tl->oneelempertile)
    gal_fits_img_write(p->binary, p->detskyname, NULL, PROGRAM_STRING);


  /* Allocate space for the mean and standard deviation. */
  p->sky=gal_data_alloc(NULL, p->input->type, p->input->ndim, tl->numtiles,
                        NULL, 0, cp->minmapsize, "SKY", p->input->unit, NULL);
  p->std=gal_data_alloc(NULL, p->input->type, p->input->ndim, tl->numtiles,
                        NULL, 0, cp->minmapsize, "STD", p->input->unit, NULL);


  /* Find the Sky and its STD on proper tiles. */
  gal_threads_spin_off(threshold_mean_std_undetected, p, tl->tottiles,
                       cp->numthreads);
  if(p->detskyname)
    {
      gal_tile_full_values_write(p->sky, tl, p->detskyname, PROGRAM_STRING);
      gal_tile_full_values_write(p->std, tl, p->detskyname, PROGRAM_STRING);
    }

  /* Get the basic information about the standard deviation
     distribution. */
  tmp=gal_statistics_median(p->std, 0);
  tmp=gal_data_copy_to_new_type(tmp, GAL_TYPE_FLOAT32);
  memcpy(&p->medstd, tmp->array, sizeof p->medstd);
  free(tmp);

  tmp=gal_statistics_minimum(p->std);
  tmp=gal_data_copy_to_new_type(tmp, GAL_TYPE_FLOAT32);
  memcpy(&p->minstd, tmp->array, sizeof p->minstd);
  free(tmp);

  tmp=gal_statistics_maximum(p->std);
  tmp=gal_data_copy_to_new_type(tmp, GAL_TYPE_FLOAT32);
  memcpy(&p->maxstd, tmp->array, sizeof p->maxstd);
  free(tmp);

  /* In case the image is in electrons or counts per second, the standard
     deviation of the noise will become smaller than unity, so we need to
     correct it in the S/N calculation. So, we'll calculate the correction
     factor here. */
  p->cpscorr = p->minstd>1 ? 1.0f : p->minstd;

  /* Interpolate and smooth the derived values. */
  threshold_interp_smooth(p, &p->sky, &p->std, p->detskyname);


  /* If a check was requested, abort NoiseChisel. */
  if(p->detskyname && !p->continueaftercheck)
    ui_abort_after_check(p, p->detskyname, "showing derivation of Sky value"
                         "and its standard deviation, or STD");
}
