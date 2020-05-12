/*********************************************************************
Common tile operations used by some Gnuastro programs, but too specific
to be in the general library.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2019-2020, Free Software Foundation, Inc.

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

#include <gnuastro/tile.h>
#include <gnuastro/pointer.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/tile-internal.h>


/* The main working function for 'threshold_no_outlier'. The main
   purpose/problem is this: when we have channels, the qthresh values for
   each channel should be treated independently. */
static void
tileinternal_no_outlier_work(gal_data_t *first, gal_data_t *second,
                             gal_data_t *third, size_t channelid,
                             size_t tottilesinch, double *outliersclip,
                             float outliersigma)
{
  gal_data_t *outlier, *nbs;
  size_t i, osize=first->size;
  size_t start=tottilesinch*channelid;
  float *oa1=NULL, *oa2=NULL, *oa3=NULL;
  float o, *arr1=NULL, *arr2=NULL, *arr3=NULL;

  /* A small sanity check. */
  if(first->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: datatype has to be float32", __func__);

  /* Correct the arrays (if necessary). IMPORTANT: The datasets are
     multi-dimensional. However, when estimating the quantile, their
     dimensionality doesn't matter (only the 'size' element is checked by
     'gal_statistics_quantile', not 'ndim' or `dsize'). So we just need to
     correct 'size' if channels are to be considered. */
  if(start || tottilesinch!=first->size)
    {
      /* Keep the original values for re-setting later. */
      oa1=first->array;
      oa2=second->array;
      if(third) oa3=third->array;

      /* Increment the array pointers. */
      first->array=gal_pointer_increment(first->array, start, first->type);
      second->array=gal_pointer_increment(second->array, start,
                                           second->type);
      if(third)
        third->array=gal_pointer_increment(third->array, start, third->type);

      /* Correct their sizes. */
      first->size=tottilesinch;
      second->size=tottilesinch;
      if(third) third->size=tottilesinch;
    }

  /* Find the quantile and remove all tiles that are more than it in the
     first array. */
  arr1=first->array;
  nbs=gal_statistics_no_blank_sorted(first, 0);
  outlier=gal_statistics_outlier_positive(nbs, nbs->size/2, outliersigma,
                                          outliersclip[0], outliersclip[1],
                                          0, 1);
  gal_data_free(nbs);
  if(outlier)
    {
      o = *((float *)(outlier->array));
      for(i=0;i<first->size;++i)
        /* Just note that we have blank (NaN) values, so to avoid doing a
           NaN check with 'isnan', we will check if the value is below the
           quantile, if it succeeds (isn't NaN and is below the quantile),
           then we'll put it's actual value, otherwise, a NaN. */
        arr1[i] = arr1[i]<=o ? arr1[i] : NAN;
      gal_data_free(outlier);
    }

  /* Second quantile threshold. We are finding the outliers independently
     on each dataset to later remove any tile that is blank in atleast one
     of them. */
  arr2=second->array;
  nbs=gal_statistics_no_blank_sorted(second, 0);
  outlier=gal_statistics_outlier_positive(nbs, nbs->size, outliersigma,
                                          outliersclip[0], outliersclip[1],
                                          0, 1);
  gal_data_free(nbs);
  if(outlier)
    {
      o = *((float *)(outlier->array));
      for(i=0;i<second->size;++i)
        arr2[i] = arr2[i]<=o ? arr2[i] : NAN;
      gal_data_free(outlier);
    }

  /* The third (if it exists). */
  if(third)
    {
      arr3=third->array;
      nbs=gal_statistics_no_blank_sorted(third, 0);
      outlier=gal_statistics_outlier_positive(nbs, nbs->size/2,
                                              outliersigma,
                                              outliersclip[0],
                                              outliersclip[1], 0, 1);
      gal_data_free(nbs);
      if(outlier)
        {
          o = *((float *)(outlier->array));
          for(i=0;i<third->size;++i)
            arr3[i] = arr3[i]<=o ? arr3[i] : NAN;
          gal_data_free(outlier);
        }
    }

  /* Make sure all three have the same NaN pixels. */
  for(i=0;i<first->size;++i)
    if( isnan(arr1[i]) || isnan(arr2[i]) || (third && isnan(arr3[i])) )
      {
        arr1[i] = arr2[i] = NAN;
        if(third) arr3[i] = NAN;
      }

  /* Correct the values, if they were changed. */
  if(start || tottilesinch!=osize)
    {
      first->array=oa1;
      second->array=oa2;
      first->size = second->size = osize;
      if(third) { third->array=oa3; third->size=osize; }
    }
}





/* Clean higher valued quantile thresholds: useful when the diffuse (almost
   flat) structures are much larger than the tile size. */
void
gal_tileinternal_no_outlier(gal_data_t *first, gal_data_t *second,
                            gal_data_t *third,
                            struct gal_tile_two_layer_params *tl,
                            double *outliersclip, float outliersigma,
                            char *filename)
{
  size_t i;

  /* A small sanity check: */
  if(first->size!=tl->tottiles)
    error(EXIT_FAILURE, 0, "%s: 'first->size' and 'tl->tottiles' must have "
          "the same value, but they don't: %zu, %zu", __func__, first->size,
          tl->tottiles);

  /* Do the work. */
  for(i=0;i<tl->totchannels;++i)
    tileinternal_no_outlier_work(first, second, third, i, tl->tottilesinch,
                                 outliersclip, outliersigma);

  /* If the user wants to see the steps. */
  if(filename)
    {
      first->name="VALUE1_NO_OUTLIER";
      second->name="VALUE2_NO_OUTLIER";
      gal_tile_full_values_write(first, tl, 1, filename, NULL, NULL);
      gal_tile_full_values_write(second, tl, 1, filename, NULL, NULL);
      first->name=second->name=NULL;
      if(third)
        {
          third->name="VALUE3_NO_OUTLIER";
          gal_tile_full_values_write(third, tl, 1, filename, NULL, NULL);
          third->name=NULL;
        }
    }
}
