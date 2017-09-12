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
struct threshold_apply_p
{
  float               *value1;
  float               *value2;
  int                    kind;
  struct noisechiselparams *p;
};




/* Apply the threshold on the tiles given to this thread. */
static void *
threshold_apply_on_thread(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct threshold_apply_p *taprm=(struct threshold_apply_p *)(tprm->params);
  struct noisechiselparams *p=taprm->p;

  size_t i, tid;
  void *tarray=NULL;
  gal_data_t *tile, *tblock=NULL;
  float *value1=taprm->value1, *value2=taprm->value2;

  /* Go over all the tiles assigned to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. */
      tid=tprm->indexs[i];
      tile=&p->cp.tl.tiles[tid];

      /* Based on the kind of threshold. */
      switch(taprm->kind)
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

          /* Apply the threshold: When the `>' comparison fails, it can be
             either because the pixel was actually smaller than the
             threshold, or that it was a NaN value. In the first case,
             return 0, in the second, return a blank value.

             We already know if a tile contains a blank value (which is a
             constant over the whole loop). So before checking if the value
             is blank, see if the tile actually has a blank value. This
             will help in efficiency, because the compiler can move this
             check out of the loop and only check for NaN values when we
             know the tile has blank pixels. */
          GAL_TILE_PO_OISET(float, uint8_t, tile, p->binary, 1, 0, {
              *o = ( *i > value1[tid]
                     ? ( *i > value2[tid] ? THRESHOLD_NO_ERODE_VALUE : 1 )
                     : ( (tile->flag & GAL_DATA_FLAG_HASBLANK) && !(*i==*i)
                         ? GAL_BLANK_UINT8 : 0 ) );
            });

          /* Revert the tile's pointers back to what they were. */
          if(p->conv) { tile->array=tarray; tile->block=tblock; }
          break;


        /* This is a Sky and Sky STD threshold. */
        case THRESHOLD_SKY_STD:

          /* See the explanation above the same step in the quantile
             threshold for an explanation. */
          GAL_TILE_PO_OISET(float, uint8_t, tile, p->binary, 1, 0, {
              *o = ( ( *i - value1[tid] > p->dthresh * value2[tid] )
                     ? 1
                     : ( (tile->flag & GAL_DATA_FLAG_HASBLANK) && !(*i==*i)
                         ? GAL_BLANK_UINT8 : 0 ) );
            });
          break;


        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we "
                "can address the problem. A value of %d had for "
                "`taprm->kind' is not valid", __func__, PACKAGE_BUGREPORT,
                taprm->kind);
        }
    }

  /* Wait until all the other threads finish. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* Apply a given threshold threshold on the tiles. */
void
threshold_apply(struct noisechiselparams *p, float *value1,
                float *value2, int kind)
{
  struct threshold_apply_p taprm={value1, value2, kind, p};
  gal_threads_spin_off(threshold_apply_on_thread, &taprm, p->cp.tl.tottiles,
                       p->cp.numthreads);
}



















/**********************************************************************/
/***************            Write S/N values          *****************/
/**********************************************************************/
void
threshold_write_sn_table(struct noisechiselparams *p, gal_data_t *insn,
                         gal_data_t *inind, char *filename,
                         gal_list_str_t *comments)
{
  gal_data_t *sn, *ind, *cols;

  /* Remove all blank elements. The index and sn values must have the same
     set of blank elements, but checking on the integer array is faster. */
  if( gal_blank_present(inind, 1) )
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
void
threshold_interp_smooth(struct noisechiselparams *p, gal_data_t **first,
                        gal_data_t **second, gal_data_t **third,
                        char *filename)
{
  gal_data_t *tmp;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;

  /* A small sanity check. */
  if( (*first)->next )
    error(EXIT_FAILURE, 0, "%s: `first' must not have any `next' pointer.",
          __func__);
  if( (*second)->next )
    error(EXIT_FAILURE, 0, "%s: `second' must not have any `next' pointer.",
          __func__);
  if( third && (*third)->next )
    error(EXIT_FAILURE, 0, "%s: `third' must not have any `next' pointer.",
          __func__);


  /* Do the interpolation of both arrays. */
  (*first)->next = *second;
  if(third) (*second)->next = *third;
  tmp=gal_interpolate_close_neighbors(*first, tl, cp->interpnumngb,
                                      cp->numthreads, cp->interponlyblank, 1);
  gal_data_free(*first);
  gal_data_free(*second);
  if(third) gal_data_free(*third);
  *first=tmp;
  *second=(*first)->next;
  if(third)
    {
      *third=(*second)->next;
      (*third)->next=NULL;
    }
  (*first)->next=(*second)->next=NULL;
  if(filename)
    {
      (*first)->name="THRESH1_INTERP";
      (*second)->name="THRESH2_INTERP";
      if(third) (*third)->name="THRESH3_INTERP";
      gal_tile_full_values_write(*first, tl, filename, NULL, PROGRAM_NAME);
      gal_tile_full_values_write(*second, tl, filename, NULL, PROGRAM_NAME);
      if(third)
        gal_tile_full_values_write(*third, tl, filename, NULL, PROGRAM_NAME);
      (*first)->name = (*second)->name = NULL;
      if(third) (*third)->name=NULL;
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

      /* Smooth the third */
      if(third)
        {
          tmp=gal_tile_full_values_smooth(*third, tl, p->smoothwidth,
                                          p->cp.numthreads);
          gal_data_free(*third);
          *third=tmp;
        }

      /* Add them to the check image. */
      if(filename)
        {
          (*first)->name="THRESH1_SMOOTH";
          (*second)->name="THRESH2_SMOOTH";
          if(third) (*third)->name="THRESH3_SMOOTH";
          gal_tile_full_values_write(*first, tl, filename, NULL,
                                     PROGRAM_NAME);
          gal_tile_full_values_write(*second, tl, filename, NULL,
                                     PROGRAM_NAME);
          if(third)
            gal_tile_full_values_write(*third, tl, filename, NULL,
                                       PROGRAM_NAME);
          (*first)->name = (*second)->name = NULL;
          if(third) (*third)->name=NULL;
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
  gal_data_t       *expand_th;
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
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
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

          /* Same for the no-erode quantile. */
          qvalue=gal_statistics_quantile(usage, p->noerodequant, 1);
          memcpy(gal_data_ptr_increment(qprm->noerode_th->array, tind, type),
                 qvalue->array, twidth);
          gal_data_free(qvalue);

          /* Same for the expansion quantile. */
          qvalue=gal_statistics_quantile(usage, p->detgrowquant, 1);
          memcpy(gal_data_ptr_increment(qprm->expand_th->array, tind, type),
                 qvalue->array, twidth);
          gal_data_free(qvalue);
        }
      else
        {
          gal_blank_write(gal_data_ptr_increment(qprm->erode_th->array,
                                                 tind, type), type);
          gal_blank_write(gal_data_ptr_increment(qprm->noerode_th->array,
                                                 tind, type), type);
          gal_blank_write(gal_data_ptr_increment(qprm->expand_th->array,
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
                       PROGRAM_NAME);


  /* Allocate space for the quantile threshold values. */
  qprm.erode_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                               tl->numtiles, NULL, 0, cp->minmapsize,
                               NULL, p->input->unit, NULL);
  qprm.noerode_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                                 tl->numtiles, NULL, 0, cp->minmapsize,
                                 NULL, p->input->unit, NULL);
  if(p->detgrowquant!=1.0f)
    qprm.expand_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                                  tl->numtiles, NULL, 0, cp->minmapsize,
                                  NULL, p->input->unit, NULL);
  else
    qprm.expand_th=NULL;


  /* Allocate temporary space for processing in each tile. */
  qprm.usage=gal_data_malloc_array(p->input->type,
                                   cp->numthreads * p->maxtcontig,
                                   __func__, "qprm.usage");


  /* Find the threshold on each tile, free the temporary processing space
     and set the blank flag on both. Since they have the same blank
     elements, it is only necessary to check one (with the `updateflag'
     value set to 1), then update the next. */
  qprm.p=p;
  gal_threads_spin_off(qthresh_on_tile, &qprm, tl->tottiles, cp->numthreads);
  free(qprm.usage);
  if( gal_blank_present(qprm.erode_th, 1) )
    {
      qprm.noerode_th->flag |= GAL_DATA_FLAG_HASBLANK;
      if(qprm.expand_th) qprm.expand_th->flag  |= GAL_DATA_FLAG_HASBLANK;
    }
  qprm.noerode_th->flag |= GAL_DATA_FLAG_BLANK_CH;
  qprm.expand_th->flag  |= GAL_DATA_FLAG_BLANK_CH;
  if(p->qthreshname)
    {
      qprm.erode_th->name="QTHRESH_ERODE";
      qprm.noerode_th->name="QTHRESH_NOERODE";
      gal_tile_full_values_write(qprm.erode_th, tl, p->qthreshname, NULL,
                                 PROGRAM_NAME);
      gal_tile_full_values_write(qprm.noerode_th, tl, p->qthreshname, NULL,
                                 PROGRAM_NAME);
      qprm.erode_th->name=qprm.noerode_th->name=NULL;

      if(qprm.expand_th)
        {
          qprm.expand_th->name="QTHRESH_EXPAND";
          gal_tile_full_values_write(qprm.expand_th, tl, p->qthreshname,
                                     NULL, PROGRAM_NAME);
          qprm.expand_th->name=NULL;
        }
    }


  /* Interpolate and smooth the derived values. */
  threshold_interp_smooth(p, &qprm.erode_th, &qprm.noerode_th,
                          &qprm.expand_th, p->qthreshname);


  /* We now have a threshold for all tiles, apply it. */
  threshold_apply(p, qprm.erode_th->array, qprm.noerode_th->array,
                  THRESHOLD_QUANTILES);


  /* Write the binary image if check is requested. */
  if(p->qthreshname && !tl->oneelempertile)
    gal_fits_img_write(p->binary, p->qthreshname, NULL, PROGRAM_NAME);


  /* Set the expansion quantile if necessary. */
  p->expand_thresh = qprm.expand_th ? qprm.expand_th : NULL;


  /* Clean up and report duration if necessary. */
  gal_data_free(qprm.erode_th);
  gal_data_free(qprm.noerode_th);
  if(!p->cp.quiet)
    {
      asprintf(&msg, "%.2f & %0.2f quantile thresholds applied.",
               p->qthresh, p->noerodequant);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* If the user wanted to check the threshold and hasn't called
     `continueaftercheck', then stop NoiseChisel. */
  if(p->qthreshname && !p->continueaftercheck)
    ui_abort_after_check(p, p->qthreshname, NULL, "quantile threshold check");
}
