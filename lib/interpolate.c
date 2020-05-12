/*********************************************************************
Interpolate - Fill blank values in a dataset
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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

#include <gnuastro/list.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/pointer.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>
#include <gnuastro/interpolate.h>
#include <gnuastro/permutation.h>

#include <gnuastro-internal/checkset.h>




/*********************************************************************/
/********************      Nearest neighbor       ********************/
/***************         (Dimension agnostic)         ****************/
/*********************************************************************/
/* These are bit-flags, so we're using hexadecimal notation. */
#define INTERPOLATE_FLAGS_NO       0
#define INTERPOLATE_FLAGS_CHECKED  0x1
#define INTERPOLATE_FLAGS_BLANK    0x2





/* Parameters for interpolation on threads. */
struct interpolate_params
{
  gal_data_t                    *input;
  size_t                           num;
  gal_data_t                      *out;
  gal_data_t                   *blanks;
  size_t                  numneighbors;
  uint8_t                *thread_flags;
  int                        onlyblank;
  gal_list_void_t            *ngb_vals;
  float (*metric)(size_t *, size_t *, size_t );

  struct gal_tile_two_layer_params *tl;
};





/* Run the interpolation on many threads. */
static void *
interpolate_close_neighbors_on_thread(void *in_prm)
{
  /* Variables that others depend on. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct interpolate_params *prm=(struct interpolate_params *)(tprm->params);
  struct gal_tile_two_layer_params *tl=prm->tl;
  int correct_index=(tl && tl->totchannels>1 && !tl->workoverch);
  gal_data_t *input=prm->input;

  /* Rest of variables. */
  void *nv;
  float dist, pdist;
  uint8_t *b, *bf, *bb;
  gal_list_void_t *tvll;
  gal_list_dosizet_t *lQ, *sQ;
  size_t ngb_counter, pind, *dinc;
  size_t i, index, fullind, chstart=0, ndim=input->ndim;
  gal_data_t *median, *tin, *tout, *tnear, *nearest=NULL;
  size_t size = (correct_index ? tl->tottilesinch : input->size);
  size_t *dsize = (correct_index ? tl->numtilesinch : input->dsize);
  size_t *icoord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                      "icoord");
  size_t *ncoord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                      "ncoord");
  uint8_t *fullflag=&prm->thread_flags[tprm->id*input->size], *flag=fullflag;


  /* Initializations. */
  bb=prm->blanks->array;
  bf=(b=fullflag)+input->size;
  dinc=gal_dimension_increment(ndim, dsize);
  do *b = *bb++ ? INTERPOLATE_FLAGS_BLANK : 0; while(++b<bf);


  /* Put the allocated space to keep the neighbor values into a structure
     for easy processing. */
  tin=input;
  for(tvll=prm->ngb_vals; tvll!=NULL; tvll=tvll->next)
    {
      nv=gal_pointer_increment(tvll->v, tprm->id*prm->numneighbors,
                               input->type);
      gal_list_data_add_alloc(&nearest, nv, tin->type, 1, &prm->numneighbors,
                              NULL, 0, -1, 1, NULL, NULL, NULL);
      tin=tin->next;
    }
  gal_list_data_reverse(&nearest);


  /* Go over all the points given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. */
      fullind=tprm->indexs[i];


      /* If the caller only wanted to interpolate over blank values and
         this value is not blank (we know from the flags), then just set
         the output value at this element to the input value and go to the
         next element. */
      if(prm->onlyblank && !(fullflag[fullind] & INTERPOLATE_FLAGS_BLANK) )
        {
          tin=input;
          for(tout=prm->out; tout!=NULL; tout=tout->next)
            {
              memcpy(gal_pointer_increment(tout->array, fullind, tin->type),
                     gal_pointer_increment(tin->array,  fullind, tin->type),
                     gal_type_sizeof(tin->type));
              tin=tin->next;
            }
          continue;
        }


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
      do *b &= ~(INTERPOLATE_FLAGS_CHECKED); while(++b<bf);


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
          if( !(flag[pind] & INTERPOLATE_FLAGS_BLANK) )
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
             if( !(flag[nind] & INTERPOLATE_FLAGS_CHECKED) )
               {
                 /* Get the coordinates of this neighbor. */
                 gal_dimension_index_to_coord(nind, ndim, dsize, ncoord);

                 /* Distance of this neighbor to the one to be filled. */
                 dist=prm->metric(icoord, ncoord, ndim);

                 /* Add this neighbor to the list. */
                 gal_list_dosizet_add(&lQ, &sQ, nind, dist);

                 /* Flag this neighbor as checked. */
                 flag[nind] |= INTERPOLATE_FLAGS_CHECKED;
               }
           } );

          /* If there are no more meshes to add to the queue, then this
             shows, there were not enough points for
             interpolation. Normally, this loop should only be exited
             through the 'currentnum>=numnearest' check above. */
          if(sQ==NULL)
            error(EXIT_FAILURE, 0, "%s: only %zu neighbors found while "
                  "you had asked to use %zu neighbors for close neighbor "
                  "interpolation", __func__, ngb_counter, prm->numneighbors);
        }

      /* Calculate the median of the values and write it in the output. */
      tout=prm->out;
      for(tnear=nearest; tnear!=NULL; tnear=tnear->next)
        {
          /* Find the median and copy it, but first, reset the flags (which
             remain from the last time). */
          tnear->flag &= ~(GAL_DATA_FLAG_SORT_CH | GAL_DATA_FLAG_BLANK_CH);
          median=gal_statistics_median(tnear, 1);
          memcpy(gal_pointer_increment(tout->array, fullind, tout->type),
                 median->array, gal_type_sizeof(tout->type));

          /* Clean up and go to next array. */
          gal_data_free(median);
          tout=tout->next;
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





/* When no interpolation is needed, then we can just copy the input into
   the output. */
static gal_data_t *
interpolate_copy_input(gal_data_t *input, int aslinkedlist)
{
  gal_data_t *tin, *tout;

  /* Make a copy of the first input. */
  tout=gal_data_copy(input);
  tout->next=NULL;

  /* If we have a linked list, copy each element. */
  if(aslinkedlist)
    {
      /* Note that we have already copied the first input. */
      for(tin=input->next; tin!=NULL; tin=tin->next)
        {
          /* Copy this dataset (will also copy flags). */
          tout->next=gal_data_copy(tin);
          tout=tout->next;
        }

      /* Output is the reverse of the input, so reverse it. */
      gal_list_data_reverse(&tout);
    }

  /* Return the copied list. */
  return tout;
}





/* Interpolate blank values in an array. If the 'tl!=NULL', then it is
   assumed that the tile values correspond to given tessellation. Such that
   'input[i]' corresponds to 'tiles[i]' in the tessellation. */
gal_data_t *
gal_interpolate_close_neighbors(gal_data_t *input,
                                struct gal_tile_two_layer_params *tl,
                                uint8_t metric, size_t numneighbors,
                                size_t numthreads, int onlyblank,
                                int aslinkedlist)
{
  gal_data_t *tin, *tout;
  struct interpolate_params prm;
  size_t ngbvnum=numthreads*numneighbors;
  int permute=(tl && tl->totchannels>1 && tl->workoverch);


  /* If there are no blank values in the array, AND we should only fill
     blank values, then simply copy the input and abort. */
  if( (input->flag | GAL_DATA_FLAG_BLANK_CH)     /* Zero bit is meaningful.*/
      && !(input->flag | GAL_DATA_FLAG_HASBLANK) /* There are no blanks.   */
      && onlyblank )                             /* Only interpolate blank.*/
    return interpolate_copy_input(input, aslinkedlist);


  /* Initialize the constant parameters. */
  prm.tl           = tl;
  prm.ngb_vals     = NULL;
  prm.input        = input;
  prm.onlyblank    = onlyblank;
  prm.numneighbors = numneighbors;
  prm.num          = aslinkedlist ? gal_list_data_number(input) : 1;


  /* Set the metric. */
  switch(metric)
    {
    case GAL_INTERPOLATE_CLOSE_METRIC_RADIAL:
      prm.metric=gal_dimension_dist_radial;
      break;
    case GAL_INTERPOLATE_CLOSE_METRIC_MANHATTAN:
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

      /* If this is a linked list, then permute remaining nodes. */
      if(aslinkedlist)
        for(tin=input->next; tin!=NULL; tin=tin->next)
          gal_permutation_apply(tin, tl->permutation);
    }


  /* Allocate space for the (first) output. */
  prm.out=gal_data_alloc(NULL, input->type, input->ndim, input->dsize,
                         input->wcs, 0, input->minmapsize,
                         input->quietmmap, NULL, input->unit, NULL);
  gal_list_void_add(&prm.ngb_vals,
                    gal_pointer_allocate(input->type, ngbvnum, 0, __func__,
                                         "prm.ngb_vals"));


  /* If we are given a list of datasets, make the necessary
     allocations. The reason we are doing this after a check of
     'aslinkedlist' is that the 'input' might have a 'next' element, but
     the caller might not have called 'aslinkedlist'. */
  prm.out->next=NULL;
  if(aslinkedlist)
    for(tin=input->next; tin!=NULL; tin=tin->next)
      {
        /* A small sanity check. */
        if( gal_dimension_is_different(input, tin) )
          error(EXIT_FAILURE, 0, "%s: all datasets in the list must have "
                "the same dimension and size", __func__);

        /* Allocate the output array for this node. */
        gal_list_data_add_alloc(&prm.out, NULL, tin->type, tin->ndim,
                                tin->dsize, tin->wcs, 0, tin->minmapsize,
                                tin->quietmmap, NULL, tin->unit, NULL);

        /* Allocate the space for the neighbor values of this input. */
        gal_list_void_add(&prm.ngb_vals,
                          gal_pointer_allocate(tin->type, ngbvnum, 0,
                                               __func__, "prm.ngb_vals"));
      }
  gal_list_data_reverse(&prm.out);
  gal_list_void_reverse(&prm.ngb_vals);


  /* Allocate space for all the flag values of all the threads here (memory
     in each thread is limited) and this is cleaner. */
  prm.thread_flags=gal_pointer_allocate(GAL_TYPE_UINT8,
                                        numthreads*input->size, 0, __func__,
                                        "prm.thread_flags");


  /* Spin off the threads. */
  gal_threads_spin_off(interpolate_close_neighbors_on_thread, &prm,
                       input->size, numthreads);


  /* If the values were permuted for the interpolation, then re-order the
     values back to their original location (so they correspond to their
     tile indexs. */
  if(permute)
    {
      gal_permutation_apply_inverse(input, tl->permutation);
      for(tout=prm.out; tout!=NULL; tout=tout->next)
        gal_permutation_apply_inverse(tout, tl->permutation);
    }


  /* The interpolated array doesn't have blank values. So set the blank
     flag to 0 and set the use-zero to 1. */
  for(tout=prm.out; tout!=NULL; tout=tout->next)
    {
      tout->flag |= GAL_DATA_FLAG_BLANK_CH;
      tout->flag &= ~GAL_DATA_FLAG_HASBLANK;
    }


  /* Clean up and return. */
  free(prm.thread_flags);
  gal_data_free(prm.blanks);
  gal_list_void_free(prm.ngb_vals, 1);
  return prm.out;
}




















/*********************************************************************/
/********************          1D on grid         ********************/
/*********************************************************************/
gsl_spline *
gal_interpolate_1d_make_gsl_spline(gal_data_t *X, gal_data_t *Y, int type_1d)
{
  size_t i, c;
  double *x, *y;
  gal_data_t *Xd, *Yd;
  gsl_spline *spline=NULL;
  const gsl_interp_type *itype=NULL;
  int Yhasblank=gal_blank_present(Y, 0);

  /* A small sanity check. */
  if(Y->ndim!=1)
    error(EXIT_FAILURE, 0, "%s: input dataset is not 1D (it is %zuD)",
          __func__, Y->ndim);
  if(X)
    {
      if( gal_dimension_is_different(X, Y) )
        error(EXIT_FAILURE, 0, "%s: when two inputs are given, they must "
              "have the same dimensions. X has %zu elements, while Y has "
              "%zu", __func__, X->size, Y->size);
      if(gal_blank_present(X, 0))
        error(EXIT_FAILURE, 0, "%s: the X dataset has blank elements",
              __func__);
    }

  /* Set the interpolation type. */
  switch(type_1d)
    {
    case GAL_INTERPOLATE_1D_LINEAR:
      itype=gsl_interp_linear;           break;
    case GAL_INTERPOLATE_1D_POLYNOMIAL:
      itype=gsl_interp_polynomial;       break;
    case GAL_INTERPOLATE_1D_CSPLINE:
      itype=gsl_interp_cspline;          break;
    case GAL_INTERPOLATE_1D_CSPLINE_PERIODIC:
      itype=gsl_interp_cspline_periodic; break;
    case GAL_INTERPOLATE_1D_AKIMA:
      itype=gsl_interp_akima;            break;
    case GAL_INTERPOLATE_1D_AKIMA_PERIODIC:
      itype=gsl_interp_akima_periodic;   break;
    case GAL_INTERPOLATE_1D_STEFFEN:
#if HAVE_DECL_GSL_INTERP_STEFFEN
      itype=gsl_interp_steffen;          break;
#else
      error(EXIT_FAILURE, 0, "%s: Steffen interpolation isn't available "
            "in the system's GNU Scientific Library (GSL). Please install "
            "a more recent GSL (version >= 2.0, released in October 2015) "
            "and rebuild Gnuastro", __func__);
#endif
    default:
      error(EXIT_FAILURE, 0, "%s: code %d not recognizable for the GSL "
            "interpolation type", __func__, type_1d);
    }

  /* Initializations. Note that if Y doesn't have any blank elements and is
     already in 'double' type, then we don't need to make a copy. */
  Yd = ( (Yhasblank || Y->type!=GAL_TYPE_FLOAT64)
         ? gal_data_copy_to_new_type(Y, GAL_TYPE_FLOAT64)
         : Y );
  Xd = ( X
         /* Has to be 'Yhasblank', we KNOW X doesn't have blank values. */
         ? ( (Yhasblank || X->type!=GAL_TYPE_FLOAT64)
             ? gal_data_copy_to_new_type(X, GAL_TYPE_FLOAT64)
             : X )
         : gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, Y->dsize, NULL,
                          0, -1, 1, NULL, NULL, NULL) );

  /* Fill in the X axis values while also removing NaN/blank elements. */
  c=0;
  x=Xd->array;
  y=Yd->array;
  for(i=0;i<Yd->size;++i)
    if( !isnan(y[i]) )
      {
        y[ c   ] = y[i];
	x[ c++ ] = X ? x[i] : i;
      }

  /* Make sure we have enough valid points for interpolation. */
  if( c>=gsl_interp_type_min_size(itype) )
    {
      spline=gsl_spline_alloc(itype, c);
      gsl_spline_init(spline, x, y, c);
    }
  else
    spline=NULL;

  /* Clean up and return. */
  if(Xd!=X) gal_data_free(Xd);
  if(Yd!=Y) gal_data_free(Yd);
  return spline;
}





/* Return 0 if all blanks were filled. */
static int
interpolate_1d_blank_write(gal_data_t *in, gsl_spline *spline,
			   gsl_interp_accel *acc)
{
  double tmp;
  int hasblank=0;
  uint8_t  *su8 =in->array, *u8 =in->array, *u8f =u8 +in->size;
  int8_t   *si8 =in->array, *i8 =in->array, *i8f =i8 +in->size;
  uint16_t *su16=in->array, *u16=in->array, *u16f=u16+in->size;
  int16_t  *si16=in->array, *i16=in->array, *i16f=i16+in->size;
  uint32_t *su32=in->array, *u32=in->array, *u32f=u32+in->size;
  int32_t  *si32=in->array, *i32=in->array, *i32f=i32+in->size;
  uint64_t *su64=in->array, *u64=in->array, *u64f=u64+in->size;
  int64_t  *si64=in->array, *i64=in->array, *i64f=i64+in->size;
  float    *sf32=in->array, *f32=in->array, *f32f=f32+in->size;
  double   *sf64=in->array, *f64=in->array, *f64f=f64+in->size;

  switch(in->type)
    {
    case GAL_TYPE_UINT8:
      do
        if(*u8==GAL_BLANK_UINT8)
          {
            /* If the evaluation is good, this function will return 0. */
            if( gsl_spline_eval_e(spline, u8-su8, acc, &tmp)==0 )
              *u8=tmp;
            else hasblank=1;
          }
      while(++u8<u8f);
      break;
    case GAL_TYPE_INT8:
      do
        if(*i8==GAL_BLANK_INT8)
          {
            if( gsl_spline_eval_e(spline, i8-si8, acc, &tmp)==0 )
              *u16=tmp;
            else hasblank=1;
          }
      while(++i8<i8f);
      break;
    case GAL_TYPE_UINT16:
      do
        if(*u16==GAL_BLANK_UINT16)
          {
            if( gsl_spline_eval_e(spline, u16-su16, acc, &tmp)==0 )
              *u16=tmp;
            else hasblank=1;
          }
      while(++u16<u16f);
      break;
    case GAL_TYPE_INT16:
      do
        if(*i16==GAL_BLANK_INT16)
          {
            if( gsl_spline_eval_e(spline, i16-si16, acc, &tmp)==0 )
              *i16=tmp;
            else hasblank=1;
          }
      while(++i16<i16f);
      break;
    case GAL_TYPE_UINT32:
      do
        if(*u32==GAL_BLANK_UINT32)
          {
            if( gsl_spline_eval_e(spline, u32-su32, acc, &tmp)==0 )
              *u32=tmp;
            else hasblank=1;
          }
      while(++u32<u32f);
      break;
    case GAL_TYPE_INT32:
      do
        if(*i32==GAL_BLANK_INT32)
          {
            if( gsl_spline_eval_e(spline, i32-si32, acc, &tmp)==0 )
              *i32=tmp;
            else hasblank=1;
          }
      while(++i32<i32f);
      break;
    case GAL_TYPE_UINT64:
      do
        if(*u64==GAL_BLANK_UINT64)
          {
            if( gsl_spline_eval_e(spline, u64-su64, acc, &tmp)==0 )
              *u64=tmp;
            else hasblank=1;
          }
      while(++u64<u64f);
      break;
    case GAL_TYPE_INT64:
      do
        if(*i64==GAL_BLANK_INT64)
          {
            if( gsl_spline_eval_e(spline, i64-si64, acc, &tmp)==0 )
              *i64=tmp;
            else hasblank=1;
          }
      while(++i64<i64f);
      break;
    case GAL_TYPE_FLOAT32:
      do
        if(isnan(*f32))
          {
            if( gsl_spline_eval_e(spline, f32-sf32, acc, &tmp)==0 )
              *f32=tmp;
            else hasblank=1;
          }
      while(++f32<f32f);
      break;
    case GAL_TYPE_FLOAT64:
      do
        if(isnan(*f64))
          {
            if( gsl_spline_eval_e(spline, f64-sf64, acc, f64) )
              hasblank=1;
          }
      while(++f64<f64f);
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: code %d is not a recognized data type",
	    __func__, in->type);
    }
  return hasblank;
}





void
gal_interpolate_1d_blank(gal_data_t *in, int type_1d)
{
  int hasblank;
  gsl_spline *spline;
  gsl_interp_accel *acc;

  /* If there are no blank elements, just return. */
  if(!gal_blank_present(in, 1)) return;

  /* Initialize the necessary structures. */
  spline=gal_interpolate_1d_make_gsl_spline(NULL, in, type_1d);

  /* If any interpolation structure was actually made. */
  if(spline)
    {
      /* Write the values in the blank elements. */
      acc=gsl_interp_accel_alloc();
      hasblank=interpolate_1d_blank_write(in, spline, acc);

      /* For a check.
      {
        size_t i;
        double *d;
        gal_data_t *check=gal_data_copy_to_new_type(in, GAL_TYPE_FLOAT64);
        d=check->array;
        for(i=0;i<check->size;++i)
          printf("%-10zu%f\n", i, d[i]);
        gal_data_free(check);
      }
      */

      /* Set the blank flags, note that 'GAL_DATA_FLAG_BLANK_CH' is already set
         by the top call to 'gal_blank_present'. */
      if(hasblank)
        in->flag |=  GAL_DATA_FLAG_HASBLANK;
      else
        in->flag &= ~GAL_DATA_FLAG_HASBLANK;

      /* Clean up. */
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
    }
}
