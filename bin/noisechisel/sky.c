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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/pointer.h>
#include <gnuastro/threads.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/tile-internal.h>

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

  int setblank, type=GAL_TYPE_FLOAT32;
  size_t i, tind, numsky, bdsize=2, ndim=p->sky->ndim;
  size_t refarea, twidth=gal_type_sizeof(GAL_TYPE_FLOAT32);
  gal_data_t *tile, *fusage, *busage, *bintile, *sigmaclip;


  /* Put the temporary usage space for this thread into a data set for easy
     processing. */
  fusage=gal_data_alloc(NULL, type, ndim, p->maxtsize, NULL, 0,
                        p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);
  busage=gal_data_alloc(NULL, GAL_TYPE_UINT8, ndim, p->maxtsize, NULL, 0,
                        p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);


  /* An empty dataset to replicate a tile on the binary array. */
  bintile=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &bdsize,
                         NULL, 0, -1, 1, NULL, NULL, NULL);
  bintile->ndim=ndim;
  free(bintile->array);
  free(bintile->dsize);
  bintile->block=p->binary;


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* Basic definitions */
      numsky=0;
      tind = tprm->indexs[i];
      tile = &p->cp.tl.tiles[tind];
      refarea = p->skyfracnoblank ? 0 : tile->size;

      /* Correct the fake binary tile's properties to be the same as this
         one, then count the number of zero valued elements in it. Note
         that the 'CHECK_BLANK' flag of 'GAL_TILE_PARSE_OPERATE' is set to
         1. So blank values in the input array are not counted. */
      bintile->size=tile->size;
      bintile->dsize=tile->dsize;
      bintile->array=gal_tile_block_relative_to_other(tile, p->binary);
      GAL_TILE_PARSE_OPERATE(tile, bintile, 1, 1, {
          if(p->skyfracnoblank) ++refarea;
          if(!*o)               ++numsky;
        });

      /* Only continue, if the fraction of Sky values is less than the
         requested fraction. */
      setblank=0;
      if( (float)(numsky)/(float)(refarea) > p->minskyfrac)
        {
          /* Re-initialize the usage array's size information (will be
             corrected to this tile's size by
             'gal_data_copy_to_allocated'). */
          busage->ndim = fusage->ndim = ndim;
          busage->size = fusage->size = p->maxtcontig;
          gal_data_copy_to_allocated(tile,    fusage);
          gal_data_copy_to_allocated(bintile, busage);


          /* Set all the non-zero pixels in 'busage' to NaN in 'fusage'. */
          busage->flag = fusage->flag = 0;
          gal_blank_flag_apply(fusage, busage);


          /* Do the sigma-clipping. */
          sigmaclip=gal_statistics_sigma_clip(fusage, p->sigmaclip[0],
                                              p->sigmaclip[1], 1, 1);

          /* When there are zero-valued pixels on the edges of the dataset
             (that have not been set to NaN/blank), given special
             conditions, the whole zero-valued region can get a binary
             value of 1 and so the Sky and its standard deviation can
             become zero. So, we need ignore such tiles. */
          if( ((float *)(sigmaclip->array))[3]==0.0 )
            setblank=1;
          else
            {
              /* Copy the sigma-clipped mean and STD to their respective
                 places in the tile arrays. But first, make sure
                 'sigmaclip' has the same type as the sky and std
                 arrays. */
              sigmaclip=gal_data_copy_to_new_type_free(sigmaclip, type);
              memcpy(gal_pointer_increment(p->sky->array, tind, type),
                     gal_pointer_increment(sigmaclip->array, 2, type),
                     twidth);
              memcpy(gal_pointer_increment(p->std->array, tind, type),
                     gal_pointer_increment(sigmaclip->array, 3, type),
                     twidth);
            }

          /* Clean up. */
          gal_data_free(sigmaclip);
        }
      else
        setblank=1;

      /* If the tile is marked as being blank, write blank values into
         it. */
      if(setblank==1)
        {
          gal_blank_write(gal_pointer_increment(p->sky->array, tind, type),
                          type);
          gal_blank_write(gal_pointer_increment(p->std->array, tind, type),
                          type);
        }
    }

  /* Clean up and wait for other threads to finish and abort. */
  bintile->array=NULL;
  bintile->dsize=NULL;
  gal_data_free(fusage);
  gal_data_free(busage);
  gal_data_free(bintile);
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
    {
      p->binary->name="DETECTED";
      gal_fits_img_write(p->binary, checkname, NULL, PROGRAM_NAME);
      p->binary->name=NULL;
    }


  /* Allocate space for the mean and standard deviation. */
  p->sky=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, p->input->ndim, tl->numtiles,
                        NULL, 0, cp->minmapsize, p->cp.quietmmap, NULL,
                        p->input->unit, NULL);
  p->std=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, p->input->ndim, tl->numtiles,
                        NULL, 0, cp->minmapsize, p->cp.quietmmap, NULL,
                        p->input->unit, NULL);


  /* Find the Sky and its STD on proper tiles. */
  gal_threads_spin_off(sky_mean_std_undetected, p, tl->tottiles,
                       cp->numthreads);
  if(checkname)
    {
      p->sky->name="SKY";
      p->std->name="STD";
      gal_tile_full_values_write(p->sky, tl, !p->ignoreblankintiles,
                                 checkname, NULL, PROGRAM_NAME);
      gal_tile_full_values_write(p->std, tl, !p->ignoreblankintiles,
                                 checkname, NULL, PROGRAM_NAME);
      p->sky->name=p->std->name=NULL;
    }


  /* Set the blank-checked bit of the arrays to zero so we are sure to
     check for blanks. */
  p->sky->flag &= ~GAL_DATA_FLAG_BLANK_CH;
  p->std->flag &= ~GAL_DATA_FLAG_BLANK_CH;


  /* Basic Sky standard deviation distribution measurements. */
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
  threshold_interp_smooth(p, &p->sky, &p->std, NULL, checkname);


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
  gal_data_t *tile;
  float *sky=p->sky->array;

  /* A small sanity check. */
  if(p->sky->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: only 'float32' type is acceptable "
          "for sky values. but 'p->sky' has type '%s'", __func__,
          gal_type_name(p->sky->type, 1));

  /* Go over all the tiles. */
  for(tid=0; tid<p->cp.tl.tottiles; ++tid)
    {
      /* For easy reading. */
      tile=&p->cp.tl.tiles[tid];

      /* Subtract the Sky value from the input image. */
      GAL_TILE_PARSE_OPERATE(tile, NULL, 0, 0, {*i-=sky[tid];});
    }
}
