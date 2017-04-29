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

#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/statistics.h>

#include "main.h"

#include "ui.h"
#include "threshold.h"








/****************************************************************
 ************            Estimate the Sky            ************
 ****************************************************************/
static void *
sky_mean_std_undetected(void *in_prm)
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
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
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
      if( (float)(numsky)/(float)(tile->size) > p->minskyfrac)
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

          /* Clean up. */
          gal_data_free(meanstd);
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
sky_and_std(struct noisechiselparams *p, char *checkname)
{
  gal_data_t *tmp;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;


  /* When the check image has the same resolution as the input, write the
     binary array as a reference to help in the comparison. */
  if(checkname && !tl->oneelempertile)
    gal_fits_img_write(p->binary, checkname, NULL, PROGRAM_STRING);


  /* Allocate space for the mean and standard deviation. */
  p->sky=gal_data_alloc(NULL, p->input->type, p->input->ndim, tl->numtiles,
                        NULL, 0, cp->minmapsize, "SKY", p->input->unit, NULL);
  p->std=gal_data_alloc(NULL, p->input->type, p->input->ndim, tl->numtiles,
                        NULL, 0, cp->minmapsize, "STD", p->input->unit, NULL);


  /* Find the Sky and its STD on proper tiles. */
  gal_threads_spin_off(sky_mean_std_undetected, p, tl->tottiles,
                       cp->numthreads);
  if(checkname)
    {
      gal_tile_full_values_write(p->sky, tl, checkname, NULL, PROGRAM_STRING);
      gal_tile_full_values_write(p->std, tl, checkname, NULL, PROGRAM_STRING);
    }

  /* Get the basic information about the standard deviation
     distribution. */
  tmp=gal_statistics_median(p->std, 0);
  tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
  memcpy(&p->medstd, tmp->array, sizeof p->medstd);
  gal_data_free(tmp);

  tmp=gal_statistics_minimum(p->std);
  tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
  memcpy(&p->minstd, tmp->array, sizeof p->minstd);
  gal_data_free(tmp);

  tmp=gal_statistics_maximum(p->std);
  tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
  memcpy(&p->maxstd, tmp->array, sizeof p->maxstd);
  gal_data_free(tmp);

  /* In case the image is in electrons or counts per second, the standard
     deviation of the noise will become smaller than unity, so we need to
     correct it in the S/N calculation. So, we'll calculate the correction
     factor here. */
  p->cpscorr = p->minstd>1 ? 1.0f : p->minstd;

  /* Interpolate and smooth the derived values. */
  threshold_interp_smooth(p, &p->sky, &p->std, checkname);


  /* If a check was requested, abort NoiseChisel. */
  if(checkname && !p->continueaftercheck)
    ui_abort_after_check(p, checkname, NULL, "showing derivation of Sky "
                         "value and its standard deviation, or STD");
}




















/****************************************************************
 ************            Subtract the Sky            ************
 ****************************************************************/
void
sky_subtract(struct noisechiselparams *p)
{
  size_t tid;
  void *tarray=NULL;
  float *sky=p->sky->array;
  gal_data_t *tile, *tblock=NULL;

  /* A small sanity check. */
  if(p->sky->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: only `float32' type is acceptable "
          "for sky values. but `p->sky' has type `%s'", __func__,
          gal_type_to_string(p->sky->type, 1));

  /* Go over all the tiles. */
  for(tid=0; tid<p->cp.tl.tottiles; ++tid)
    {
      /* For easy reading. */
      tile=&p->cp.tl.tiles[tid];

      /* First subtract the Sky value from the input image. */
      GAL_TILE_PARSE_OPERATE({*i-=sky[tid];}, tile, NULL, 0, 0);

      /* Change to the convolved image. */
      tarray=tile->array;
      tblock=tile->block;
      tile->array=gal_tile_block_relative_to_other(tile, p->conv);
      tile->block=p->conv;

      /* The threshold is always low. So for the majority of non-NaN
         pixels in the image, the condition above will be true. If we
         come over a NaN pixel, then by definition of NaN, all
         conditionals will fail.

         If an image doesn't have any NaN pixels, only the pixels below
         the threshold have to be checked for a NaN which are by
         definition a very small fraction of the total pixels. And if
         there are NaN pixels in the image. */
      GAL_TILE_PARSE_OPERATE({*i-=sky[tid];}, tile, NULL, 0, 0);

      /* Revert back to the original block. */
      tile->array=tarray;
      tile->block=tblock;
    }
}
