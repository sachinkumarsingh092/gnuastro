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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/tile.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/statistics.h>
#include <gnuastro/interpolate.h>
#include <gnuastro/permutation.h>

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
  size_t i, osize=first->size;
  size_t start=tottilesinch*channelid;
  float *oa1=NULL, *oa2=NULL, *oa3=NULL;
  gal_data_t *nbs, *outlier_p, *outlier_n;
  float o_p, o_n, *arr1=NULL, *arr2=NULL, *arr3=NULL;

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
        third->array=gal_pointer_increment(third->array, start,
                                           third->type);

      /* Correct their sizes. */
      first->size=tottilesinch;
      second->size=tottilesinch;
      if(third) third->size=tottilesinch;
    }

  /* Find the quantile and remove all tiles that are more than it in the
     first array. */
  arr1=first->array;
  nbs=gal_statistics_no_blank_sorted(first, 0);
  outlier_p=gal_statistics_outlier_bydistance(1, nbs, nbs->size/2,
                                              outliersigma, outliersclip[0],
                                              outliersclip[1], 1, 1);
  outlier_n=gal_statistics_outlier_bydistance(0, nbs, nbs->size/2,
                                              outliersigma, outliersclip[0],
                                              outliersclip[1], 1, 1);
  /* For a check.
  {
    float *med;
    gal_data_t *median=gal_statistics_median(nbs, 1);
    float *out_n=outlier_n->array, *out_p=outlier_p->array;
    med=median->array;
    printf("vals: %f (out_n), %f (med), %f (out_p)\n",
           out_n[0], med[0], out_p[0]);
    exit(0);
  }
  */

  /* Clean up the temporary 'nbs' array. */
  gal_data_free(nbs);

  /* If outliers exist, then implement them. */
  if(outlier_p)
    {
      o_p = *((float *)(outlier_p->array));
      gal_data_free(outlier_p);
      if(outlier_n)
        {
          /* For easy reading, put the negative outlier into 'o_n'. */
          o_n = *((float *)(outlier_n->array));
          gal_data_free(outlier_n);

          /* See description below, it just includes a negative outlier. */
          for(i=0;i<first->size;++i)
            arr1[i] = arr1[i]<o_p ? (arr1[i]>o_n ? arr1[i] : NAN) : NAN;
        }
      else
          /* Just note that we have blank (NaN) values, so to avoid doing a
             NaN check with 'isnan', we will check if the value is below
             the quantile, if it succeeds (isn't NaN and is below the
             quantile), then we'll put it's actual value, otherwise, a
             NaN. */
        for(i=0;i<first->size;++i)
          arr1[i] = arr1[i]<o_p ? arr1[i] : NAN;
    }
  else
    if(outlier_n)
      {
        o_n = *((float *)(outlier_n->array));
        for(i=0;i<first->size;++i)
          arr1[i] = arr1[i]>o_n ? arr1[i] : NAN;
        gal_data_free(outlier_n);
      }

  /* Second quantile threshold. We are finding the outliers independently
     on each dataset to later remove any tile that is blank in atleast one
     of them. */
  arr2=second->array;
  nbs=gal_statistics_no_blank_sorted(second, 0);
  outlier_p=gal_statistics_outlier_bydistance(1, nbs, nbs->size,
                                              outliersigma, outliersclip[0],
                                              outliersclip[1], 1, 1);
  outlier_n=gal_statistics_outlier_bydistance(0, nbs, nbs->size,
                                              outliersigma, outliersclip[0],
                                              outliersclip[1], 1, 1);
  gal_data_free(nbs);
  if(outlier_p)
    {
      o_p = *((float *)(outlier_p->array));
      gal_data_free(outlier_p);
      if(outlier_n)
        {
          o_n = *((float *)(outlier_n->array));
          for(i=0;i<first->size;++i)
            arr2[i] = arr2[i]<o_p ? (arr2[i]>o_n ? arr2[i] : NAN) : NAN;
          gal_data_free(outlier_n);
        }
      else
        for(i=0;i<first->size;++i)
          arr2[i] = arr2[i]<o_p ? arr2[i] : NAN;
    }
  else
    if(outlier_n)
      {
        o_n = *((float *)(outlier_n->array));
        for(i=0;i<first->size;++i)
          arr2[i] = arr2[i]>o_n ? arr2[i] : NAN;
        gal_data_free(outlier_n);
      }

  /* The third (if it exists). */
  if(third)
    {
      arr3=third->array;
      nbs=gal_statistics_no_blank_sorted(third, 0);
      outlier_p=gal_statistics_outlier_bydistance(1, nbs, nbs->size/2,
                                                  outliersigma,
                                                  outliersclip[0],
                                                  outliersclip[1], 1, 1);
      outlier_n=gal_statistics_outlier_bydistance(0, nbs, nbs->size/2,
                                                  outliersigma,
                                                  outliersclip[0],
                                                  outliersclip[1], 1, 1);
      gal_data_free(nbs);
      if(outlier_p)
        {
          o_p = *((float *)(outlier_p->array));
          gal_data_free(outlier_p);
          if(outlier_n)
            {
              o_n = *((float *)(outlier_n->array));
              for(i=0;i<first->size;++i)
                arr3[i] = arr3[i]<o_p ? (arr3[i]>o_n ? arr3[i] : NAN) : NAN;
              gal_data_free(outlier_n);
            }
          else
            for(i=0;i<first->size;++i)
              arr3[i] = arr3[i]<o_p ? arr3[i] : NAN;
        }
      else
        if(outlier_n)
          {
            o_n = *((float *)(outlier_n->array));
            for(i=0;i<first->size;++i)
              arr3[i] = arr3[i]>o_n ? arr3[i] : NAN;
            gal_data_free(outlier_n);
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



















/*************************************************************/
/************        Local outlier removal        ************/
/*************************************************************/
#define TILEINTERNAL_OUTLIER_FLAGS_NO           0
#define TILEINTERNAL_OUTLIER_FLAGS_NGB_CHECKED  0x1
#define TILEINTERNAL_OUTLIER_FLAGS_BLANK        0x2

struct tileinternal_outlier_local
{
  gal_data_t                    *input;
  gal_data_t                  *measure;
  gal_data_t                   *blanks;
  size_t                  numneighbors;
  uint8_t                *thread_flags;
  gal_list_void_t            *ngb_vals;
  float (*metric)(size_t *, size_t *, size_t );

  struct gal_tile_two_layer_params *tl;
};





/* Run the interpolation on many threads. */
static void *
gal_tileinternal_no_outlier_local_on_thread(void *in_prm)
{
  /* Low-level variables that others depend on. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct tileinternal_outlier_local *prm=
    (struct tileinternal_outlier_local *)(tprm->params);

  /* Higher-level variables. */
  struct gal_tile_two_layer_params *tl=prm->tl;
  int correct_index=(tl && tl->totchannels>1 && !tl->workoverch);
  gal_data_t *input=prm->input;

  /* Rest of variables. */
  void *nv;
  uint8_t *b, *bf, *bb;
  gal_list_void_t *tvll;
  size_t ngb_counter, pind;
  gal_list_dosizet_t *lQ, *sQ;
  gal_data_t *tin, *tnear, *nearest=NULL;
  float dist, pdist, *tnarr, *marr=prm->measure->array;
  size_t i, index, fullind, chstart=0, ndim=input->ndim;
  size_t size = (correct_index ? tl->tottilesinch : input->size);
  size_t *dsize = (correct_index ? tl->numtilesinch : input->dsize);
  size_t *icoord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                      "icoord");
  size_t *ncoord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                      "ncoord");
  uint8_t *flag, *fullflag=&prm->thread_flags[tprm->id*input->size];

  /* Based on the above. */
  size_t *dinc=gal_dimension_increment(ndim, dsize);

  /* Initialize the flags array. We need two flags during this processing:
     1) to see if there are blanks. 2) to see if a neighbor has been
     checked. These are both binary (0 or 1). So to avoid wasting space, we
     will use bits to store them. We start with only setting the blank flag
     once for the whole thread. Then for each interpolated pixel, we reset
     the neighbor-check flag. */
  flag=fullflag;
  bb=prm->blanks->array;
  bf=(b=fullflag)+input->size;
  do *b = *bb++ ? TILEINTERNAL_OUTLIER_FLAGS_BLANK : 0; while(++b<bf);


  /* Put the allocated space to keep the neighbor values into a structure
     for easy processing. */
  tin=input;
  for(tvll=prm->ngb_vals; tvll!=NULL; tvll=tvll->next)
    {
      nv=gal_pointer_increment(tvll->v, tprm->id*prm->numneighbors,
                               input->type);
      gal_list_data_add_alloc(&nearest, nv, tin->type, 1,
                              &prm->numneighbors, NULL, 0, -1, 1,
                              NULL, NULL, NULL);
      tin=tin->next;
    }
  gal_list_data_reverse(&nearest);


  /* Go over all the points given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. */
      fullind=tprm->indexs[i];


      /* If we are on a blank element, then ignore this pixel. */
      if( (fullflag[fullind] & TILEINTERNAL_OUTLIER_FLAGS_BLANK) )
        { marr[fullind]=NAN; continue; }


      /* Correct the index (if necessary). When the values come from a
         tiled dataset, the caller might want to interpolate the values of
         each channel separately (not mix values from different
         channels). In such a case, the tiles of each channel (and their
         values in 'input' are contiguous. So we need to correct
         'tprm->indexs[i]' (which is the index over the whole tessellation,
         including all channels). */
      if(correct_index)
        {
          /* Index of this tile in its channel. */
          index = fullind % tl->tottilesinch;

          /* Index of the first tile in this channel. */
          chstart = (fullind / tl->tottilesinch) * tl->tottilesinch;

          /* Set the channel's starting pointer for the flags. */
          flag = gal_pointer_increment(fullflag, chstart, GAL_TYPE_UINT8);
        }
      else
        {
          chstart=0;
          index=fullind;
        }


      /* Reset all checked bits in the flags array to 0. */
      ngb_counter=0;
      bf=(b=flag)+size;
      do *b &= ~(TILEINTERNAL_OUTLIER_FLAGS_NGB_CHECKED); while(++b<bf);


      /* Get the coordinates of this pixel (to be interpolated). */
      gal_dimension_index_to_coord(index, ndim, dsize, icoord);


      /* Start parsing the neighbors. We will use a two-way ordered linked
         list structure. To start from the nearest and go out to the
         farthest. */
      lQ=sQ=NULL;
      gal_list_dosizet_add(&lQ, &sQ, index, 0.0f);
      while(sQ)
        {
          /* Pop-out (p) an index from the queue: */
          pind=gal_list_dosizet_pop_smallest(&lQ, &sQ, &pdist);

          /* If this isn't a blank value then add its values to the list of
             neighbor values. Note that we didn't check whether the values
             were blank or not when adding this pixel to the queue. */
          if( !(flag[pind] & TILEINTERNAL_OUTLIER_FLAGS_BLANK) )
            {
              tin=input;
              for(tnear=nearest; tnear!=NULL; tnear=tnear->next)
                {
                  memcpy(gal_pointer_increment(tnear->array, ngb_counter,
                                               tin->type),
                         gal_pointer_increment(tin->array, chstart+pind,
                                               tin->type),
                         gal_type_sizeof(tin->type));
                  tin=tin->next;
                }

              /* If we have filled all the elements clean up the linked
                 list and break out. */
              if(++ngb_counter>=prm->numneighbors)
                {
                  if(lQ) gal_list_dosizet_free(lQ);
                  break;
                }
            }

          /* Go over all the neighbors of this popped pixel and add them to
             the list of neighbors to be checked. */
          GAL_DIMENSION_NEIGHBOR_OP(pind, ndim, dsize, 1, dinc,
           {
             /* Only look at neighbors that have not been checked. VERY
                IMPORTANT: we must not check for blank values here,
                otherwise we won't be able to parse over extended blank
                regions. */
             if( !(flag[nind] & TILEINTERNAL_OUTLIER_FLAGS_NGB_CHECKED) )
               {
                 /* Get the coordinates of this neighbor. */
                 gal_dimension_index_to_coord(nind, ndim, dsize, ncoord);

                 /* Distance of this neighbor to the one to be filled. */
                 dist=prm->metric(icoord, ncoord, ndim);

                 /* Add this neighbor to the list. */
                 gal_list_dosizet_add(&lQ, &sQ, nind, dist);

                 /* Flag this neighbor as checked. */
                 flag[nind] |= TILEINTERNAL_OUTLIER_FLAGS_NGB_CHECKED;
               }
           } );

          /* If there are no more meshes to add to the queue, then this
             shows, there were not enough points for
             interpolation. Normally, this loop should only be exited
             through the 'currentnum>=numnearest' check above. */
          if(sQ==NULL)
            error(EXIT_FAILURE, 0, "%s: only %zu neighbors found while "
                  "you had asked to use %zu neighbors for close neighbor "
                  "interpolation", __func__, ngb_counter,
                  prm->numneighbors);
        }

      /* Calculate the desired statistic, and write it in the output. */
      for(tnear=nearest; tnear!=NULL; tnear=tnear->next)
        {
          /* First, reset the sorting flags (which remain from the last
             time). */
          tnear->flag &= ~(GAL_DATA_FLAG_SORT_CH | GAL_DATA_FLAG_BLANK_CH);

          /* For a check on the values.
          { size_t i; float *f=tnear->array;
            for(i=0;i<tnear->size;++i) printf("%f\n", f[i]); } */

          /* Sort the elements, then find the difference between the
             maximium and the value that is just after the minimum. We are
             doing this because the scatter in the minimum can be large. */
          tnarr=tnear->array;
          gal_statistics_sort_increasing(tnear);
          marr[fullind] = tnarr[tnear->size-1]-tnarr[1];
        }
    }


  /* Clean up. */
  for(tnear=nearest; tnear!=NULL; tnear=tnear->next) tnear->array=NULL;
  gal_list_data_free(nearest);
  free(icoord);
  free(ncoord);
  free(dinc);


  /* Wait for all the other threads to finish and return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





void
gal_tileinternal_no_outlier_local(gal_data_t *input, gal_data_t *second,
                                  gal_data_t *third,
                                  struct gal_tile_two_layer_params *tl,
                                  uint8_t metric, size_t numneighbors,
                                  size_t numthreads, double *outliersclip,
                                  double outliersigma, char *filename)
{
  gal_data_t *othresh;
  float *base, *f, *ff, thresh;
  struct tileinternal_outlier_local prm;
  size_t owindow, ngbvnum=numthreads*numneighbors;
  int permute=(tl && tl->totchannels>1 && tl->workoverch);


  /* Sanity checks. */
  if(numneighbors<=3)
    error(EXIT_FAILURE, 0, "interpnumngb has to be larger than 3, but "
          "is currently %zu", numneighbors);
  if(input->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
          "the problem. The input to this function (not NoiseChisel) "
          "should be in 32-bit floating point, but it is %s", __func__,
          PACKAGE_BUGREPORT, gal_type_name(input->type, 1));
  if(second && second->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
          "the problem. The 'second' argument to this function (not "
          "NoiseChisel) should be in 32-bit floating point, but it is "
          "%s", __func__, PACKAGE_BUGREPORT, gal_type_name(input->type, 1));
  if(third && third->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
          "the problem. The 'third' argument to this function (not "
          "NoiseChisel) should be in 32-bit floating point, but it is "
          "%s", __func__, PACKAGE_BUGREPORT, gal_type_name(input->type, 1));
  if(second && gal_dimension_is_different(input, second) )
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
          "the problem. The 'second' argument to this function (not "
          "NoiseChisel) doesn't have the same size as the input",
          __func__, PACKAGE_BUGREPORT);
  if(third && gal_dimension_is_different(input, third) )
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
          "the problem. The 'third' argument to this function (not "
          "NoiseChisel) doesn't have the same size as the input",
          __func__, PACKAGE_BUGREPORT);


  /* Initialize the constant parameters. */
  prm.tl           = tl;
  prm.ngb_vals     = NULL;
  prm.input        = input;
  prm.numneighbors = numneighbors;


  /* Set the distance metric. */
  switch(metric)
    {
    case GAL_INTERPOLATE_NEIGHBORS_METRIC_RADIAL:
      prm.metric=gal_dimension_dist_radial;
      break;
    case GAL_INTERPOLATE_NEIGHBORS_METRIC_MANHATTAN:
      prm.metric=gal_dimension_dist_manhattan;
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: %d is not a valid metric identifier",
            __func__, metric);
    }


  /* Flag the blank values. */
  prm.blanks=gal_blank_flag(input);


  /* If the input is from a tile structure and the user has asked to ignore
     channels, then re-order the values. */
  if(permute)
    {
      /* Prepare the permutation (if necessary/not already defined). */
      gal_tile_full_permutation(tl);

      /* Re-order values to ignore channels (if necessary). */
      gal_permutation_apply(input, tl->permutation);
      gal_permutation_apply(prm.blanks, tl->permutation);
    }


  /* Necessary allocations and basic checks. if we are given a list of
     datasets, make the necessary allocations. The reason we are doing this
     after a check of 'aslinkedlist' is that the 'input' might have a
     'next' element, but the caller might not have called
     'aslinkedlist'. */
  prm.measure=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, input->ndim,
                             input->dsize, input->wcs, 0, input->minmapsize,
                             input->quietmmap, NULL, input->unit, NULL);
  gal_list_void_add(&prm.ngb_vals,
                    gal_pointer_allocate(input->type, ngbvnum, 0,
                                         __func__, "prm.ngb_vals"));


  /* Allocate space for all the flag values of all the threads here (memory
     in each thread is limited) and this is cleaner. */
  prm.thread_flags=gal_pointer_allocate(GAL_TYPE_UINT8,
                                        numthreads*input->size, 0,
                                        __func__, "prm.thread_flags");


  /* Spin off the threads. */
  gal_threads_spin_off(gal_tileinternal_no_outlier_local_on_thread,
                       &prm, input->size, numthreads, input->minmapsize,
                       input->quietmmap);


  /* Find the outliers in the distribution, we will start from the first
     third of the cases to find the first outlier. Note that this should
     not be done in-place because we need the 'measure' arrray
     afterwards. */
  owindow=(prm.measure->size - gal_blank_number(prm.measure, 1))/3;
  othresh=gal_statistics_outlier_bydistance(1, prm.measure, owindow,
                                            outliersigma, outliersclip[0],
                                            outliersclip[1], 0, 1);


  /* If an outlier threshold was actually found, then mask all the tiles
     larger than that value. */
  if(othresh)
    {
      base=prm.measure->array;
      ff=(f=input->array)+input->size;
      thresh=((float *)(othresh->array))[0];
      do { *f = isnan(*f) ? *f : (*base>thresh ? NAN : *f); ++base; }
      while(++f<ff);
    }


  /* For a check.
  printf("measure-threshold: %f\n", thresh);
  if(permute)
    gal_permutation_apply_inverse(prm.measure, tl->permutation);
  gal_tile_full_values_write(prm.measure, tl, 1, "measure.fits",
                             NULL, NULL);
  */


  /* If the values were permuted for the interpolation, then re-order the
     values back to their original location (so they correspond to their
     tile indexs. */
  if(permute)
    gal_permutation_apply_inverse(input, tl->permutation);


  /* If the 'other' arrays are given, set all the blank elements here to
     blank there too. */
  if(second)
    {
      base=input->array; ff=(f=second->array)+second->size;
      do { *f = isnan(*base++) ? NAN : *f ;} while(++f<ff);
    }
  if(third)
    {
      base=input->array; ff=(f=third->array)+third->size;
      do { *f = isnan(*base++) ? NAN : *f ;} while(++f<ff);
    }


  /* Write the check images if necessary. */
  if(filename)
    {
      input->name="VALUE1_NO_OUTLIER";
      gal_tile_full_values_write(input, tl, 1, filename, NULL, NULL);
      input->name=NULL;
      if(second)
        {
          second->name="VALUE2_NO_OUTLIER";
          gal_tile_full_values_write(second, tl, 1, filename, NULL, NULL);
          second->name=NULL;
        }
      if(third)
        {
          third->name="VALUE3_NO_OUTLIER";
          gal_tile_full_values_write(third, tl, 1, filename, NULL, NULL);
          third->name=NULL;
        }
    }


  /* Clean up and return. */
  gal_data_free(othresh);
  free(prm.thread_flags);
  gal_data_free(prm.blanks);
  gal_data_free(prm.measure);
  gal_list_void_free(prm.ngb_vals, 1);
}
