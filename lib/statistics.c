/*********************************************************************
Statistical functions.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/qsort.h>
#include <gnuastro/pointer.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/checkset.h>










/****************************************************************
 ********               Simple statistics                 *******
 ****************************************************************/
/* Return the number of non-blank elements in an array as a single element,
   'size_t' type data structure. */
gal_data_t *
gal_statistics_number(gal_data_t *input)
{
  size_t counter=0, dsize=1;
  gal_data_t *out=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &dsize,
                                 NULL, 1, -1, 1, NULL, NULL, NULL);

  /* If there is no blank values in the input, then the total number is
     just the size. */
  if(gal_blank_present(input, 0)) /* '{}' necessary for 'else'. */
    { GAL_TILE_PARSE_OPERATE(input, NULL, 0, 1, {++counter;}); }
  else
    counter = input->size;

  /* Write the value into memory. */
  *((size_t *)(out->array)) = counter;
  return out;
}





/* Return the minimum (non-blank) value of a dataset in the same type as
   the dataset. */
gal_data_t *
gal_statistics_minimum(gal_data_t *input)
{
  size_t dsize=1, n=0;
  gal_data_t *out=gal_data_alloc(NULL, gal_tile_block(input)->type, 1,
                                 &dsize, NULL, 1, -1, 1, NULL, NULL, NULL);

  /* See if the input actually has any elements. */
  if(input->size)
    {
      /* Initialize the output with the maximum possible value. */
      gal_type_max(out->type, out->array);

      /* Parse the full input. */
      GAL_TILE_PARSE_OPERATE( input, out, 0, 1,
                              {*o = *i < *o ? *i : *o; ++n;} );
    }

  /* If there were no usable elements, set the output to blank, then
     return. */
  if(n==0) gal_blank_write(out->array, out->type);
  return out;
}





/* Return the maximum (non-blank) value of a dataset in the same type as
   the dataset. */
gal_data_t *
gal_statistics_maximum(gal_data_t *input)
{
  size_t dsize=1, n=0;
  gal_data_t *out=gal_data_alloc(NULL, gal_tile_block(input)->type, 1,
                                 &dsize, NULL, 1, -1, 1, NULL, NULL, NULL);
  /* See if the input actually has any elements. */
  if(input->size)
    {
      /* Initialize the output with the minimum possible value. */
      gal_type_min(out->type, out->array);

      /* Parse the full input. */
      GAL_TILE_PARSE_OPERATE(input, out, 0, 1,
                             {*o = *i > *o ? *i : *o; ++n;});
    }

  /* If there were no usable elements, set the output to blank, then
     return. */
  if(n==0) gal_blank_write(out->array, out->type);
  return out;
}





/* Return the sum of the input dataset as a single element dataset of type
   float64. */
gal_data_t *
gal_statistics_sum(gal_data_t *input)
{
  size_t dsize=1, n=0;
  gal_data_t *out=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, 1, NULL, NULL, NULL);

  /* See if the input actually has any elements. */
  if(input->size)
    /* Parse the dataset. Note that in 'gal_data_alloc' we set the 'clear'
       flag to 1, so it will be 0.0f. */
    GAL_TILE_PARSE_OPERATE(input, out, 0, 1, {++n; *o += *i;});

  /* If there were no usable elements, set the output to blank, then
     return. */
  if(n==0) gal_blank_write(out->array, out->type);
  return out;
}





/* Return the mean of the input dataset as a float64 type single-element
   dataset. */
gal_data_t *
gal_statistics_mean(gal_data_t *input)
{
  size_t dsize=1, n=0;
  gal_data_t *out=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, 1, NULL, NULL, NULL);

  /* See if the input actually has any elements. */
  if(input->size)
    /* Parse the dataset. Note that in 'gal_data_alloc' we set the 'clear'
       flag to 1, so it will be 0.0f. */
    GAL_TILE_PARSE_OPERATE(input, out, 0, 1, {++n; *o += *i;});

  /* Above, we calculated the sum and number, so if there were any elements
     in the dataset ('n!=0'), divide the sum by the number, otherwise, put
     a blank value in the output. */
  if(n) *((double *)(out->array)) /= n;
  else gal_blank_write(out->array, out->type);
  return out;
}





/* Return the standard deviation of the input dataset as a single element
   dataset of type float64. */
gal_data_t *
gal_statistics_std(gal_data_t *input)
{
  size_t dsize=1, n=0;
  double s=0.0f, s2=0.0f, *o;
  gal_data_t *out=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, 1, NULL, NULL, NULL);

  /* See if the input actually has any elements. */
  o=out->array;
  switch(input->size)
    {
    /* No inputs. */
    case 0: o[0]=GAL_BLANK_FLOAT64; break;

    /* When we only have a single element, theoretically the standard
       deviation should be 0. But due to floating-point errors, it will
       probably not be. So we'll manually set it to zero. */
    case 1: o[0]=0; break;

    /* More than one element. */
    default:
      GAL_TILE_PARSE_OPERATE(input, out, 0, 1, {++n; s += *i; s2 += *i * *i;});
      o[0]=sqrt( (s2-s*s/n)/n );
      break;
    }

  /* Return the output dataset. */
  return out;
}





/* Return the mean and standard deviation of a dataset in one run in type
   float64. The output is a two element data structure, with the first
   value being the mean and the second value the standard deviation. */
gal_data_t *
gal_statistics_mean_std(gal_data_t *input)
{
  size_t dsize=2, n=0;
  double s=0.0f, s2=0.0f, *o;
  gal_data_t *out=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, 1, NULL, NULL, NULL);

  /* See if the input actually has any elements. */
  o=out->array;
  switch(input->size)
    {
    /* No inputs. */
    case 0: o[0]=o[1]=GAL_BLANK_FLOAT64; break;

    /* When we only have a single element, theoretically the standard
       deviation should be 0. But due to floating-point errors, it will
       probably not be. So we'll manually set it to zero. */
    case 1:
      GAL_TILE_PARSE_OPERATE(input, out, 0, 1, {s+=*i;});
      o[0]=s; o[1]=0;
      break;

    /* More than one element. */
    default:
      GAL_TILE_PARSE_OPERATE(input, out, 0, 1, {++n; s += *i; s2 += *i * *i;});
      o[0]=s/n;
      o[1]=sqrt( (s2-s*s/n)/n );
      break;
    }

  /* Return the output dataset. */
  return out;
}





/* The input is a sorted array with no blank values, we want the median
   value to be put inside the already allocated space which is pointed to
   by 'median'. It is in the same type as the input. */
#define MED_IN_SORTED(IT) {                                             \
    IT *a=sorted->array;                                                \
    *(IT *)median = n%2 ? a[n/2]  : (a[n/2]+a[n/2-1])/2;                \
  }
static void
statistics_median_in_sorted_no_blank(gal_data_t *sorted, void *median)
{
  size_t n=sorted->size;

  /* Do the processing if there are actually any elements. */
  if(sorted->size)
    switch(sorted->type)
      {
      case GAL_TYPE_UINT8:     MED_IN_SORTED( uint8_t  );    break;
      case GAL_TYPE_INT8:      MED_IN_SORTED( int8_t   );    break;
      case GAL_TYPE_UINT16:    MED_IN_SORTED( uint16_t );    break;
      case GAL_TYPE_INT16:     MED_IN_SORTED( int16_t  );    break;
      case GAL_TYPE_UINT32:    MED_IN_SORTED( uint32_t );    break;
      case GAL_TYPE_INT32:     MED_IN_SORTED( int32_t  );    break;
      case GAL_TYPE_UINT64:    MED_IN_SORTED( uint64_t );    break;
      case GAL_TYPE_INT64:     MED_IN_SORTED( int64_t  );    break;
      case GAL_TYPE_FLOAT32:   MED_IN_SORTED( float    );    break;
      case GAL_TYPE_FLOAT64:   MED_IN_SORTED( double   );    break;
      default:
        error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
              __func__, sorted->type);
      }
  else
    gal_blank_write(median, sorted->type);
}





/* Return the median value of the dataset in the same type as the input as
   a one element dataset. If the 'inplace' flag is set, the input data
   structure will be modified: it will have no blank values and will be
   sorted (increasing). */
gal_data_t *
gal_statistics_median(gal_data_t *input, int inplace)
{
  size_t dsize=1;
  gal_data_t *nbs=gal_statistics_no_blank_sorted(input, inplace);;
  gal_data_t *out=gal_data_alloc(NULL, nbs->type, 1, &dsize, NULL, 1, -1,
                                 1, NULL, NULL, NULL);

  /* Write the median. */
  if(nbs->size)
    statistics_median_in_sorted_no_blank(nbs, out->array);
  else
    gal_blank_write(out->array, out->type);

  /* Clean up (if necessary), then return the output */
  if(nbs!=input) gal_data_free(nbs);
  return out;
}




/* For a given size, return the index (starting from zero) that is at the
   given quantile.  */
size_t
gal_statistics_quantile_index(size_t size, double quantile)
{
  double floatindex;

  /* Some sanity checks. */
  if(size==0)
    {
      error(0, 0, "%s: 'size' is 0. The quantile is not defined for "
              "a zero-sized array\n", __func__);
      return GAL_BLANK_SIZE_T;
    }
  if(quantile<0.0f || quantile>1.0f)
    error(EXIT_FAILURE, 0, "%s: the input quantile should be between 0.0 "
          "and 1.0 (inclusive). You have asked for %g", __func__, quantile);

  /* Find the index of the quantile. */
  floatindex=(double)(size-1)*quantile;

  /*
  printf("quantile: %f, size: %zu, findex: %f\n", quantile, size, floatindex);
  */
  /* Note that in the conversion from float to size_t, the floor
     integer value of the float will be used. */
  if( floatindex - (int)floatindex > 0.5 )
    return floatindex+1;
  else
    return floatindex;
}





/* Return a single element dataset of the same type as input keeping the
   value that has the given quantile. */
gal_data_t *
gal_statistics_quantile(gal_data_t *input, double quantile, int inplace)
{
  void *blank;
  int increasing;
  size_t dsize=1, index;
  gal_data_t *nbs=gal_statistics_no_blank_sorted(input, inplace);
  gal_data_t *out=gal_data_alloc(NULL, nbs->type, 1, &dsize,
                                 NULL, 1, -1, 1, NULL, NULL, NULL);

  /* Only continue processing if there are non-blank elements. */
  if(nbs->size)
    {
      /* Set the increasing value. */
      increasing = nbs->flag & GAL_DATA_FLAG_SORTED_I;

      /* Find the index of the quantile, note that if it sorted in
         decreasing order, then we'll need to get the index of the inverse
         quantile. */
      index=gal_statistics_quantile_index(nbs->size,
                                          ( increasing
                                            ? quantile
                                            : (1.0f - quantile) ) );

      /* Write the value at this index into the output. */
      if(index==GAL_BLANK_SIZE_T)
        {
          blank=gal_pointer_allocate(nbs->type, 1, 0, __func__, "blank");
          memcpy(out->array, blank, gal_type_sizeof(nbs->type));
          free(blank);
        }
      else
        memcpy(out->array,
               gal_pointer_increment(nbs->array, index, nbs->type),
               gal_type_sizeof(nbs->type));
    }
  else
    gal_blank_write(out->array, out->type);

  /* Clean up and return. */
  if(nbs!=input) gal_data_free(nbs);
  return out;
}





/* Return the index of the (first) point in the sorted dataset that has the
   closest value to 'value' (which has to be the same type as the 'input'
   dataset). */
#define STATS_QFUNC_IND(IT) {                                           \
    IT *r, *a=nbs->array, *af=a+nbs->size, v=*((IT *)(value->array));   \
                                                                        \
    /* For a reference. Since we are comparing with the previous */     \
    /* element, we need to start with the second element.*/             \
    r=a++;                                                              \
                                                                        \
    /* Increasing array: */                                             \
    if( *a < *(a+1) )                                                   \
      {                                                                 \
        if( v>=*r )                                                     \
          {                                                             \
            do if(*a>v) { if( v - *(a-1) < *a - v ) --a; break; }       \
            while(++a<af);                                              \
            parsed=1;                                                   \
          }                                                             \
      }                                                                 \
                                                                        \
    /* Decreasing array. */                                             \
    else                                                                \
      {                                                                 \
        if(v<=*r)                                                       \
          {                                                             \
            do if(*a<v) { if( *(a-1) - v < v - *a ) --a; break; }       \
            while(++a<af);                                              \
            parsed=1;                                                   \
          }                                                             \
      }                                                                 \
                                                                        \
    /* Set the difference if the value is actually in the range. */     \
    if(parsed && a<af) index = a-r;                                     \
  }
size_t
gal_statistics_quantile_function_index(gal_data_t *input, gal_data_t *value,
                                       int inplace)
{
  int parsed=0;
  size_t index=GAL_BLANK_SIZE_T;
  gal_data_t *nbs=gal_statistics_no_blank_sorted(input, inplace);

  /* A sanity check. */
  if(nbs->type!=value->type)
    error(EXIT_FAILURE, 0, "%s: the types of the input dataset and requested "
          "value have to be the same", __func__);

  /* Only continue processing if we have non-blank elements. */
  if(nbs->size)
    /* Find the result: */
    switch(nbs->type)
      {
      case GAL_TYPE_UINT8:     STATS_QFUNC_IND( uint8_t  );     break;
      case GAL_TYPE_INT8:      STATS_QFUNC_IND( int8_t   );     break;
      case GAL_TYPE_UINT16:    STATS_QFUNC_IND( uint16_t );     break;
      case GAL_TYPE_INT16:     STATS_QFUNC_IND( int16_t  );     break;
      case GAL_TYPE_UINT32:    STATS_QFUNC_IND( uint32_t );     break;
      case GAL_TYPE_INT32:     STATS_QFUNC_IND( int32_t  );     break;
      case GAL_TYPE_UINT64:    STATS_QFUNC_IND( uint64_t );     break;
      case GAL_TYPE_INT64:     STATS_QFUNC_IND( int64_t  );     break;
      case GAL_TYPE_FLOAT32:   STATS_QFUNC_IND( float    );     break;
      case GAL_TYPE_FLOAT64:   STATS_QFUNC_IND( double   );     break;
      default:
        error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
              __func__, nbs->type);
      }
  else
    {
      error(0, 0, "%s: no non-blank elements. The quantile function is not "
            "defined for a zero-sized array\n", __func__);
      index=GAL_BLANK_SIZE_T;
    }

  /* Clean up and return. */
  if(nbs!=input) gal_data_free(nbs);
  return index;
}





/* Return the quantile function of the given value as float64. */
#define STATS_QFUNC(IT) {                                               \
    IT *a=nbs->array, v=*((IT *)(value->array));                        \
                                                                        \
    /* Increasing array: */                                             \
    if( *a < *(a+1) )                                                   \
      d[0] = v<*a ? -INFINITY : INFINITY;                               \
                                                                        \
    /* Decreasing array. */                                             \
    else                                                                \
      d[0] = v>*a ? INFINITY : -INFINITY;                               \
  }
gal_data_t *
gal_statistics_quantile_function(gal_data_t *input, gal_data_t *value,
                                 int inplace)
{
  double *d;
  size_t dsize=1;
  gal_data_t *nbs=gal_statistics_no_blank_sorted(input, inplace);
  gal_data_t *out=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, 1, NULL, NULL, NULL);
  size_t ind=gal_statistics_quantile_function_index(input, value, inplace);

  /* Only continue processing if there are non-blank values. */
  if(nbs->size)
    {
      /* Note that counting of the index starts from 0, so for the quantile
         we should divided by (size - 1). */
      d=out->array;
      if(ind==GAL_BLANK_SIZE_T)
        {
          /* See if the value is larger or smaller than the input's minimum
             or maximum. */
          switch(nbs->type)
            {
            case GAL_TYPE_UINT8:     STATS_QFUNC( uint8_t  );     break;
            case GAL_TYPE_INT8:      STATS_QFUNC( int8_t   );     break;
            case GAL_TYPE_UINT16:    STATS_QFUNC( uint16_t );     break;
            case GAL_TYPE_INT16:     STATS_QFUNC( int16_t  );     break;
            case GAL_TYPE_UINT32:    STATS_QFUNC( uint32_t );     break;
            case GAL_TYPE_INT32:     STATS_QFUNC( int32_t  );     break;
            case GAL_TYPE_UINT64:    STATS_QFUNC( uint64_t );     break;
            case GAL_TYPE_INT64:     STATS_QFUNC( int64_t  );     break;
            case GAL_TYPE_FLOAT32:   STATS_QFUNC( float    );     break;
            case GAL_TYPE_FLOAT64:   STATS_QFUNC( double   );     break;
            default:
              error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
                    __func__, nbs->type);
            }
        }
      else
        d[0] = (double)ind / ((double)(nbs->size - 1));
    }
  else
    gal_blank_write(out->array, out->type);

  /* Clean up and return. */
  if(nbs!=input) gal_data_free(nbs);
  return out;
}





/* Pull out unique elements */
#define UNIQUE_BYTYPE(TYPE) {                                           \
    size_t i, j;                                                        \
    TYPE *a=out->array, b;                                              \
                                                                        \
    /* Write the blank value for this type into 'b'. */                 \
    gal_blank_write(&b, out->type);                                     \
                                                                        \
    /* Go over the elements, and set the duplicates to blank. */        \
    /* Note that for integers and floats, the behavior of blank/NaN */  \
    /* differs: for floats (NaN), we can identify a blank using the  */ \
    /* fact that by definition, NaN!=NaN. */                            \
    if(b==b)                                                            \
      for(i=0;i<out->size;++i)                                          \
        { if(a[i]!=b)    for(j=i+1;j<out->size;++j) if(a[i]==a[j]) a[j]=b;} \
    else                                                                \
      for(i=0;i<out->size;++i)                                          \
        { if(a[i]==a[i]) for(j=i+1;j<out->size;++j) if(a[i]==a[j]) a[j]=b;} \
  }

gal_data_t *
gal_statistics_unique(gal_data_t *input, int inplace)
{
  gal_data_t *out = inplace ? input : gal_data_copy(input);

  /* Since we are replacing the repeated elements with blank, re-set the
     blank flags. */
  out->flag &= ~GAL_DATA_FLAG_BLANK_CH; /* Set bit to 0. */
  out->flag &= ~GAL_DATA_FLAG_HASBLANK; /* Set bit to 0. */

  /* Set all non-unique elements to blank. */
  switch(out->type)
    {
    case GAL_TYPE_UINT8:   UNIQUE_BYTYPE( uint8_t  ); break;
    case GAL_TYPE_INT8:    UNIQUE_BYTYPE( int8_t   ); break;
    case GAL_TYPE_UINT16:  UNIQUE_BYTYPE( uint16_t ); break;
    case GAL_TYPE_INT16:   UNIQUE_BYTYPE( int16_t  ); break;
    case GAL_TYPE_UINT32:  UNIQUE_BYTYPE( uint32_t ); break;
    case GAL_TYPE_INT32:   UNIQUE_BYTYPE( int32_t  ); break;
    case GAL_TYPE_UINT64:  UNIQUE_BYTYPE( uint64_t ); break;
    case GAL_TYPE_INT64:   UNIQUE_BYTYPE( int64_t  ); break;
    case GAL_TYPE_FLOAT32: UNIQUE_BYTYPE( float    ); break;
    case GAL_TYPE_FLOAT64: UNIQUE_BYTYPE( double   ); break;
    default:
      error(EXIT_FAILURE, 0, "the 'unique' operator doesn't support type "
            "code '%u'", out->type);
    }

  /* Remove all blank elements (note that 'gal_blank_remove' also corrects
     the size of the dataset and sets it to 1D). */
  gal_blank_remove_realloc(out);
  return out;
}




















/*********************************************************************/
/*****************              Mode           ***********************/
/*********************************************************************/
/* Main structure to keep mode parameters. */
struct statistics_mode_params
{
  gal_data_t   *data;   /* Sorted input dataset with no blank values. */
  size_t        lowi;   /* Lower quantile of interval.                */
  size_t        midi;   /* Index of the mid-interval point.           */
  size_t        midd;   /* Maximum CDF distance at the middle point.  */
  size_t       highi;   /* Higher quantile of interval.               */
  float    tolerance;   /* Tolerance level to terminate search.       */
  size_t    numcheck;   /* Number of pixels after mode to check.      */
  size_t    interval;   /* Interval to check pixels.                  */
  float   mirrordist;   /* Distance after mirror to check ( x STD).   */
};





/* Macros for the mode finding algorithm. */
#define MODE_MIN_Q        0.01f  /* Mode search lower interval quantile.  */
#define MODE_MAX_Q        0.55f  /* Mode search higher interval quantile. */
#define MODE_GOOD_LQ      0.02f  /* Least acceptable mode quantile.       */
#define MODE_SYM_LOW_Q    0.01f  /* Lower quantile to get symmetricity.   */
#define MODE_GOLDEN_RATIO 1.618034f /* Golden ratio: (1+sqrt(5))/2.       */
#define MODE_TWO_TAKE_GR  0.38197f  /* 2 - Golden ratio.                  */
#define MODE_MIRROR_ABOVE (size_t)(-1) /* Mirror is above the result.     */




/*
  Given a mirror point ('m'), return the maximum distance between the
  mirror distribution and the original distribution.

  The basic idea behind finding the mode is comparing the mirrored CDF
  (where the mirror is a test for the mode) with the original CDF for a
  given point. The job of this function is to return the maximum distance,
  given a mirror point. It takes the index of the mirror that is to be
  checked, it then finds the maximum difference between the mirrored CDF
  about the given point and the input CDF.

  'zf' keeps the value at the mirror (zero) point.  'i' is used to count
  the pixels after the mirror in the mirror distribution. So 'm+i' is the
  index of the mirrored distribution and mf=zf+(zf-a[m-i])=2*zf-a[m-i] is
  the mirrored flux at this point. Having found 'mf', we find the 'j' such
  that a[m+j] has the nearest flux to 'mf'.

  The desired difference between the input CDF and the mirrored one
  for each 'i' is then simply: 'j-i'.

  Once 'i' is incremented, 'mf' will increase, so to find the new 'j' we
  don't need to begin looking from 'j=0'. Remember that the array is
  sorted, so the desired 'j' is definitely larger than the previous
  'j'. So, if we keep the previous 'j' in 'prevj' then, all we have to do
  is to start incrementing 'j' from 'prevj'. This will really help in
  speeding up the job :-D. Only for the first element, 'prevj=0'. */
#define MIRR_MAX_DIFF(IT) {                                             \
    IT *a=p->data->array, zf=a[m], mf=2*zf-a[m-i];                      \
                                                                        \
    /* When a[m+j]>mf, we have reached the last pixel to check. Now, */ \
    /* we just have to see which one of a[m+j-1] or a[m+j] is closer */ \
    /* to 'mf'. We then change 'j' accordingly and break out of the  */ \
    /* 'j' loop. */                                                     \
    for(j=prevj;j<size-m;++j)                                           \
      if(a[m+j]>mf)                                                     \
        {                                                               \
          if( a[m+j]-mf < mf-a[m+j-1] )                                 \
            break;                                                      \
          else                                                          \
            {                                                           \
              j--;                                                      \
              break;                                                    \
            }                                                           \
        }                                                               \
  }

static size_t
mode_mirror_max_index_diff(struct statistics_mode_params *p, size_t m)
{
  /* The variables:
   i:        Index on mirror distribution.
   j:        Index on input distribution.
   prevj:    Index of previously checked point in the actual array.
   mf:       (in macro) Value that is approximately equal in both
             distributions.                                          */
  size_t i, j, absdiff, prevj=0, size=p->data->size;
  size_t  maxdiff=0, errordiff=p->mirrordist*sqrt(m);

  /*
  printf("###############\n###############\n");
  printf("### Mirror pixel: %zu (mirrordist: %f, sqrt(m): %f)\n", m,
         p->mirrordist, sqrt(m));
  printf("###############\n###############\n");
  */

  /* Go over the mirrored points. */
  for(i=1; i<p->numcheck && i<=m && m+i<size ;i+=p->interval)
    {
      /* Find 'j': the index of the closest point in the original
         distribution that has a value similar to the mirror
         distribution. */
      switch(p->data->type)
        {
        case GAL_TYPE_UINT8:     MIRR_MAX_DIFF( uint8_t  );   break;
        case GAL_TYPE_INT8:      MIRR_MAX_DIFF( int8_t   );   break;
        case GAL_TYPE_UINT16:    MIRR_MAX_DIFF( uint16_t );   break;
        case GAL_TYPE_INT16:     MIRR_MAX_DIFF( int16_t  );   break;
        case GAL_TYPE_UINT32:    MIRR_MAX_DIFF( uint32_t );   break;
        case GAL_TYPE_INT32:     MIRR_MAX_DIFF( int32_t  );   break;
        case GAL_TYPE_UINT64:    MIRR_MAX_DIFF( uint64_t );   break;
        case GAL_TYPE_INT64:     MIRR_MAX_DIFF( int64_t  );   break;
        case GAL_TYPE_FLOAT32:   MIRR_MAX_DIFF( float    );   break;
        case GAL_TYPE_FLOAT64:   MIRR_MAX_DIFF( double   );   break;
        default:
          error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
                __func__, p->data->type);
        }

      /*
      printf("i:%-5zu j:%-5zu diff:%-5d maxdiff: %zu\n",
             i, j, (int)j-(int)i, maxdiff);
      */

      /* The index of the actual CDF corresponding the the mirrored flux
         has been found. We want the mirrored distribution to be within the
         actual distribution, not beyond it, so the only acceptable results
         are when i<j. But we also have noise, so we can't simply use that
         as the criterion, small 'j's with 'i>j' are acceptable. So, only
         when 'i>j+errordiff' the result is not acceptable! */
      if(i>j+errordiff)
        {
          maxdiff = MODE_MIRROR_ABOVE;
          break;
        }
      absdiff  = i>j ? i-j : j-i;
      if(absdiff>maxdiff) maxdiff=absdiff;

      prevj=j;
    }

  /* Return the maximum difference  */
  return maxdiff;
}





/* Find the mode through the Golden-section search. It is assumed that
   'mode_mirror_max_index_diff' has one minimum (within the statistical
   errors) in the function. To find that minimum, the golden section search
   algorithm is going to used. Read the Wikipedia article for a very nice
   introduction.

   In summary we will constantly be finding middle points in the given
   interval and thus decreasing the interval until a certain tolerance is
   reached.

   If the input interval is on points 'a' and 'b', then the middle point
   (lets call it 'c', where c>a and c<b) to test should be positioned such
   that (b-c)/(c-a)=MODE_GOLDEN_RATIO. Once we open up this relation, we
   can find c using:

    c = ( b + MODE_GOLDEN_RATIO * a ) / ( 1 + MODE_GOLDEN_RATIO )

   We need a fourth point to be placed between. With this configuration,
   the probing point is located at: */
static size_t
mode_golden_section(struct statistics_mode_params *p)
{
  size_t di, dd;

  /* Find the probing point in the larger interval. */
  if(p->highi-p->midi > p->midi-p->lowi)
    di = p->midi + MODE_TWO_TAKE_GR * (float)(p->highi-p->midi);
  else
    di = p->midi - MODE_TWO_TAKE_GR * (float)(p->midi-p->lowi);

  /* Since these are all indexs (and positive) we don't need an absolute
     value, highi is also always larger than lowi! In some cases, the first
     (standard) condition might be satisfied, while highi-lowi<=2. In such
     cases, also jump out! */
  if( (p->highi - p->lowi) < p->tolerance*(p->midi+di)
      || (p->highi - p->lowi) <= 3)
    return (p->highi+p->lowi)/2;

  /* Find the maximum difference for this mirror point. */
  dd = mode_mirror_max_index_diff(p, di);

  /*------------------------------------------------------------------
  {
  static int counter=1;
  char outname[500], command[1000];
  char histsname[500], cfpsname[500];
  sprintf(outname, "%dcmp.pdf", counter);
  sprintf(cfpsname, "%dcfps.txt", counter);
  sprintf(histsname, "%dhists.txt", counter);
  gal_mode_make_mirror_plots(p->sorted, p->size, di, histsname, cfpsname);
  sprintf(command, "./plot.py %s %s %s", histsname, cfpsname, outname);
  system(command);
  }
  -------------------------------------------------------------------*/

  /*
  printf("lowi:%-5zu\tmidi:%-5zu(midd: %d)\thighi:%-5zu ----> "
         "dq: %-5zu di: %d\n",
         p->lowi, p->midi, (int)p->midd, p->highi,
         di, (int)dd);
  */

  /* +++++++++++++ Start of addition to the golden section search.

     The mirrored distribution's cumulative frequency plot has be lower
     than the actual's cfp. If it isn't, 'di' will be MODE_MIRROR_ABOVE. In
     this case, the normal golden section minimization is not going to give
     us what we want. So we have this modification. In such cases, we want
     the search to go to the lower interval. */
  if(dd==MODE_MIRROR_ABOVE)
    {
      if( p->midi < di )
        {
          p->highi=di;
          return mode_golden_section(p);
        }
      else
        {
          p->highi=p->midi;
          p->midi=di;
          p->midd=dd;
          return mode_golden_section(p);
        }
    }
  /* End of addition to the golden section search. +++++++++++++*/

  /* This is the standard golden section search: */
  if(dd<p->midd)
    {
      if(p->highi-p->midi > p->midi-p->lowi)
        {
          p->lowi  = p->midi;
          p->midi  = di;
          p->midd  = dd;
          return mode_golden_section(p);
        }
      else
        {
          p->highi = p->midi;
          p->midi  = di;
          p->midd  = dd;
          return mode_golden_section(p);
        }
    }
  else
    {
      if(p->highi-p->midi > p->midi-p->lowi)
        {
          p->highi = di;
          return mode_golden_section(p);
        }
      else
        {
          p->lowi  = di;
          return mode_golden_section(p);
        }
    }
}





/* Once the mode is found, we need to do a quality control. This quality
   control is the measure of its symmetricity. Let's assume the mode index
   is at 'm', since an index is just a count, from the Poisson
   distribution, the error in 'm' is sqrt(m).

   Now, let's take 'b' to be the first point that the difference between
   the cumulative distribution of the mirror and actual data deviate more
   than sqrt(m). For a scale parameter, lets assume that the index of 5% of
   'm' is 'a'. We could have taken the distribution minimum, but the
   scatter in the minimum can be too high!

   Now, the "symmetricity" of the mode can be defined as: (b-m)/(m-a). For
   a completly symmetric mode, this should be 1. Note that the search for
   'b' only goes to the 95% of the distribution.  */
#define MODE_SYM(IT) {                                                  \
    IT *a=p->data->array, af=0, bf=0, mf=0, fi;                         \
                                                                        \
    /* Set the values at the mirror and at 'a' (see above). */          \
    mf=a[m];                                                            \
    af=a[ gal_statistics_quantile_index(2*m+1, MODE_SYM_LOW_Q) ];       \
    if(mf<=af) return 0;                                                \
                                                                        \
    /* This loop is very similar to that of */                          \
    /* 'mode_mirror_max_index_diff'. It will find the index where the */\
    /* difference between the two cumulative frequency plots exceeds */ \
    /* that of the error in the mirror index.*/                         \
    for(i=1; i<topi-m ;i+=1)                                            \
      {                                                                 \
        fi=2*mf-a[m-i];                                                 \
                                                                        \
        for(j=prevj;j<size-m;++j)                                       \
          if(a[m+j]>fi)                                                 \
            {                                                           \
              if( a[m+j]-fi < fi-a[m+j-1] )                             \
                break;                                                  \
              else                                                      \
                {                                                       \
                  j--;                                                  \
                  break;                                                \
                }                                                       \
            }                                                           \
                                                                        \
        if(i>j+errdiff || j>i+errdiff)                                  \
          {                                                             \
            bi=m+i;                                                     \
            break;                                                      \
          }                                                             \
        prevj=j;                                                        \
      }                                                                 \
                                                                        \
    /* bi==0 shows that no point with a larger difference could be */   \
    /* found. So bi should be set to the end of the search region. */   \
    if(bi==0) bi=topi;                                                  \
                                                                        \
    bf = *(IT *)b_val = a[bi];                                          \
    /*printf("%zu: %f,%f,%f\n", m, (double)af, (double)mf, (double)bf);*/ \
                                                                        \
    /* For a bad result, return 0 (which will not output any mode). */  \
    return bf==af ? 0 : (double)(bf-mf)/(double)(mf-af);                \
  }
static double
mode_symmetricity(struct statistics_mode_params *p, size_t m, void *b_val)
{
  size_t i, j, bi=0, topi, errdiff, prevj=0, size=p->data->size;

  /* Set the basic constants. */
  topi = 2*m>size-1 ? size-1 : 2*m;
  errdiff = p->mirrordist * sqrt(m);

  /* Do the process. */
  switch(p->data->type)
    {
    case GAL_TYPE_UINT8:      MODE_SYM( uint8_t  );    break;
    case GAL_TYPE_INT8:       MODE_SYM( int8_t   );    break;
    case GAL_TYPE_UINT16:     MODE_SYM( uint16_t );    break;
    case GAL_TYPE_INT16:      MODE_SYM( int16_t  );    break;
    case GAL_TYPE_UINT32:     MODE_SYM( uint32_t );    break;
    case GAL_TYPE_INT32:      MODE_SYM( int32_t  );    break;
    case GAL_TYPE_UINT64:     MODE_SYM( uint64_t );    break;
    case GAL_TYPE_INT64:      MODE_SYM( int64_t  );    break;
    case GAL_TYPE_FLOAT32:    MODE_SYM( float    );    break;
    case GAL_TYPE_FLOAT64:    MODE_SYM( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, p->data->type);
    }

  /* Control shouldn't reach here! */
  error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s so we can "
        "address the problem. Control must not have reached the end of this "
        "function", __func__, PACKAGE_BUGREPORT);
  return NAN;
}





/* Return the mode and related parameters in a float64 'gal_data_t' with
   the following elements in its array, the array:

      array[0]: mode
      array[1]: mode quantile.
      array[2]: symmetricity.
      array[3]: value at the end of symmetricity.

  The inputs are:

    - 'input' is the input dataset, it doesn't have to be sorted and can
      have blank values.

    - 'mirrordist' is the maximum distance after the mirror point to check
      as a multiple of sigma.

    - 'inplace' is either 0 or 1. If it is 1 and the input array has blank
      values and is not sorted, then the removal of blank values and
      sorting will occur in-place (input will be modified): all blank
      elements in the input array will be removed and it will be sorted. */
gal_data_t *
gal_statistics_mode(gal_data_t *input, float mirrordist, int inplace)
{
  double *oa;
  size_t modeindex;
  size_t dsize=4, mdsize=1;
  struct statistics_mode_params p;
  int type=gal_tile_block(input)->type;
  gal_data_t *tmptype=gal_data_alloc(NULL, type, 1, &mdsize, NULL, 1, -1, 1,
                                     NULL, NULL, NULL);
  gal_data_t *b_val=gal_data_alloc(NULL, type, 1, &mdsize, NULL, 1, -1, 1,
                                   NULL, NULL, NULL);
  gal_data_t *out=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, 1, NULL, NULL, NULL);


  /* A small sanity check. */
  if(mirrordist<=0)
    error(EXIT_FAILURE, 0, "%s: %f not acceptable as a value to "
          "'mirrordist'. Only positive values can be given to it",
          __func__, mirrordist);


  /* Make sure the input doesn't have blank values and is sorted.  */
  p.data=gal_statistics_no_blank_sorted(input, inplace);


  /* It can happen that the whole array is blank. In such cases,
     'p.data->size==0', so set all output elements to NaN and return. */
  oa=out->array;
  if(p.data->size==0) { oa[0]=oa[1]=oa[2]=oa[3]=NAN; return out; }


  /* Basic constants. */
  p.tolerance    = 0.01;
  p.mirrordist   = mirrordist;
  p.numcheck     = p.data->size/2;


  /* Fill in the interval: Checking every single element is over-kill, so
     if the dataset is large enough, we'll set an interval to only check
     elements at an interval (so only 1000 elements are checked). */
  p.interval = p.numcheck>1000 ? p.numcheck/1000 : 1;


  /* Set the lower and higher acceptable indexes for the mode based on
     quantiles. */
  p.lowi  = gal_statistics_quantile_index(p.data->size, MODE_MIN_Q);
  p.highi = gal_statistics_quantile_index(p.data->size, MODE_MAX_Q);


  /* Having set the low and higher interval values, we will set the first
     middle point and also the maximum distance on that point. This is
     necessary to start the iteration. */
  p.midi = ( ( (float)p.highi + MODE_GOLDEN_RATIO * (float)p.lowi )
             / ( 1 + MODE_GOLDEN_RATIO ) );
  p.midd = mode_mirror_max_index_diff(&p, p.midi);


  /* Do the golden-section search iteration, read the mode value from the
     input array and save it in the 'tmptype' data structure that has the
     same type as the input. */
  modeindex = mode_golden_section(&p);
  memcpy( tmptype->array,
          gal_pointer_increment(p.data->array, modeindex, p.data->type),
          gal_type_sizeof(p.data->type) );


  /* Convert the mode (which is in the same type as the input at this
     stage) to float64. */
  tmptype=gal_data_copy_to_new_type_free(tmptype, GAL_TYPE_FLOAT64);


  /* Put the first three values into the output structure. */
  oa[0] = *((double *)(tmptype->array));
  oa[1] = ((double)modeindex) / ((double)(p.data->size-1));
  oa[2] = mode_symmetricity(&p, modeindex, b_val->array);


  /* If the symmetricity is good, put it in the output, otherwise set all
     output values to NaN. */
  if(oa[2]>GAL_STATISTICS_MODE_GOOD_SYM)
    {
      b_val=gal_data_copy_to_new_type_free(b_val, GAL_TYPE_FLOAT64);
      oa[3] = *((double *)(b_val->array));
    }
  else oa[0]=oa[1]=oa[2]=oa[3]=NAN;


  /* For a check:
  printf("mode: %g\nquantile: %g\nsymmetricity: %g\nsym value: %g\n",
         oa[0], oa[1], oa[2], oa[3]);
  */

  /* Clean up (if necessary), then return the output */
  if(p.data!=input) gal_data_free(p.data);
  gal_data_free(tmptype);
  gal_data_free(b_val);
  return out;
}





/* Make the mirror array. */
#define STATS_MKMIRROR(IT) {                                            \
    IT *a=noblank_sorted->array, *m=mirror->array;                      \
    IT zf=a[index];                                                     \
    *mirror_val=zf;                                                     \
    for(i=0;i<=index;++i) m[i]       = a[i];                            \
    for(i=1;i<=index;++i) m[index+i] = 2 * zf - m[index - i];           \
  }
static gal_data_t *
statistics_make_mirror(gal_data_t *noblank_sorted, size_t index,
                       double *mirror_val)
{
  size_t i, dsize = 2*index+1;
  gal_data_t *mirror=gal_data_alloc(NULL, noblank_sorted->type, 1, &dsize,
                                    NULL, 1, -1, 1, NULL, NULL, NULL);

  /* Make sure the index is less than or equal to the number of
     elements. */
  if( index >= noblank_sorted->size )
    error(EXIT_FAILURE, 0, "%s: the index value must be less than or equal "
          "to the number of elements in the input, but it isn't: index: "
          "%zu, size of input: %zu", __func__, index, noblank_sorted->size);

  /* Fill in the mirror array. */
  switch(noblank_sorted->type)
    {
    case GAL_TYPE_UINT8:     STATS_MKMIRROR( uint8_t  );     break;
    case GAL_TYPE_INT8:      STATS_MKMIRROR( int8_t   );     break;
    case GAL_TYPE_UINT16:    STATS_MKMIRROR( uint16_t );     break;
    case GAL_TYPE_INT16:     STATS_MKMIRROR( int16_t  );     break;
    case GAL_TYPE_UINT32:    STATS_MKMIRROR( uint32_t );     break;
    case GAL_TYPE_INT32:     STATS_MKMIRROR( int32_t  );     break;
    case GAL_TYPE_UINT64:    STATS_MKMIRROR( uint64_t );     break;
    case GAL_TYPE_INT64:     STATS_MKMIRROR( int64_t  );     break;
    case GAL_TYPE_FLOAT32:   STATS_MKMIRROR( float    );     break;
    case GAL_TYPE_FLOAT64:   STATS_MKMIRROR( double   );     break;
    }

  /* Return the mirrored distribution. */
  return mirror;
}





/* Make a mirrored histogram and cumulative frequency plot with the mirror
   distribution of the input with a value at 'value'.

   The output is a linked list of data structures: the first is the bins
   with one bin at the mirror point, the second is the histogram with a
   maximum of one and the third is the cumulative frequency plot. */
gal_data_t *
gal_statistics_mode_mirror_plots(gal_data_t *input, gal_data_t *value,
                                 size_t numbins, int inplace,
                                 double *mirror_val)
{
  gal_data_t *mirror, *bins, *hist, *cfp;
  gal_data_t *nbs=gal_statistics_no_blank_sorted(input, inplace);
  size_t ind=gal_statistics_quantile_function_index(nbs, value, inplace);

  /* Only continue if we actually have non-blank elements. */
  if(nbs->size==0) return NULL;

  /* If the given mirror was outside the range of the input, then index
     will be 0 (below the range) or -1 (above the range), in that case, we
     should return NULL. */
  if(ind==-1 || ind==0)
    return NULL;


  /* Make the mirror array. */
  mirror=statistics_make_mirror(nbs, ind, mirror_val);


  /* Set the bins for histogram and cdf. */
  bins=gal_statistics_regular_bins(mirror, NULL, numbins, *mirror_val);


  /* Make the histogram: set it's maximum value to 1 for a nice comparison
     with the CDF. */
  hist=gal_statistics_histogram(mirror, bins, 0, 1);


  /* Make the cumulative frequency plot. */
  cfp=gal_statistics_cfp(mirror, bins, 1);


  /* Set the pointers to make a table and return. */
  bins->next=hist;
  hist->next=cfp;
  return bins;
}



















/****************************************************************
 ********                      Sort                       *******
 ****************************************************************/
/* Check if the given dataset is sorted. */
enum is_sorted_return
{
  STATISTICS_IS_SORTED_NOT,                 /* ==0: by C standard. */
  STATISTICS_IS_SORTED_INCREASING,
  STATISTICS_IS_SORTED_DECREASING,
};

#define IS_SORTED(IT) {                                                 \
  IT *aa=input->array, *a=input->array, *af=a+input->size-1;            \
  if(a[1]>=a[0]) do if( *(a+1) < *a ) break; while(++a<af);             \
  else           do if( *(a+1) > *a ) break; while(++a<af);             \
  out=( a==af                   /* It reached the end of the array. */  \
          ? ( aa[1]>=aa[0]                                              \
                ? STATISTICS_IS_SORTED_INCREASING                       \
                : STATISTICS_IS_SORTED_DECREASING )                     \
          : STATISTICS_IS_SORTED_NOT );                                 \
  }

int
gal_statistics_is_sorted(gal_data_t *input, int updateflags)
{
  int out=GAL_BLANK_INT16; /* On some systems, int may be 16-bits wide. */

  /* If the flags are already set, don't bother going over the dataset. */
  if( input->flag & GAL_DATA_FLAG_SORT_CH )
    return ( input->flag & GAL_DATA_FLAG_SORTED_I
             ? STATISTICS_IS_SORTED_INCREASING
             : ( input->flag & GAL_DATA_FLAG_SORTED_D
                 ? STATISTICS_IS_SORTED_DECREASING
                 : STATISTICS_IS_SORTED_NOT ) );

  /* Parse the array (if necessary). */
  switch(input->size)
    {
    case 0:
      error(EXIT_FAILURE, 0, "%s: input dataset has 0 elements", __func__);

    /* A one-element dataset can be considered, sorted, so we'll say its
       increasing. */
    case 1:
      out=STATISTICS_IS_SORTED_INCREASING;
      break;

    /* Do the check when there is more than one element. */
    default:
      switch(input->type)
        {
        case GAL_TYPE_UINT8:     IS_SORTED( uint8_t  );    break;
        case GAL_TYPE_INT8:      IS_SORTED( int8_t   );    break;
        case GAL_TYPE_UINT16:    IS_SORTED( uint16_t );    break;
        case GAL_TYPE_INT16:     IS_SORTED( int16_t  );    break;
        case GAL_TYPE_UINT32:    IS_SORTED( uint32_t );    break;
        case GAL_TYPE_INT32:     IS_SORTED( int32_t  );    break;
        case GAL_TYPE_UINT64:    IS_SORTED( uint64_t );    break;
        case GAL_TYPE_INT64:     IS_SORTED( int64_t  );    break;
        case GAL_TYPE_FLOAT32:   IS_SORTED( float    );    break;
        case GAL_TYPE_FLOAT64:   IS_SORTED( double   );    break;
        default:
          error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
                __func__, input->type);
        }
    }

  /* Update the flags, if required. */
  if(updateflags)
    {
      input->flag |= GAL_DATA_FLAG_SORT_CH;
      switch(out)
        {
        case STATISTICS_IS_SORTED_NOT:
          input->flag &= ~GAL_DATA_FLAG_SORTED_I;
          input->flag &= ~GAL_DATA_FLAG_SORTED_D;
          break;

        case STATISTICS_IS_SORTED_INCREASING:
          input->flag |=  GAL_DATA_FLAG_SORTED_I;
          input->flag &= ~GAL_DATA_FLAG_SORTED_D;
          break;

        case STATISTICS_IS_SORTED_DECREASING:
          input->flag &= ~GAL_DATA_FLAG_SORTED_I;
          input->flag |=  GAL_DATA_FLAG_SORTED_D;
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
                "the problem. The value %d is not recognized for 'out'",
                __func__, PACKAGE_BUGREPORT, out);
        }
    }
  return out;
}





/* This function is ignorant to blank values, if you want to make sure
   there is no blank values, you can call 'gal_blank_remove' first. */
#define STATISTICS_SORT(QSORT_F) {                                      \
    qsort(input->array, input->size, gal_type_sizeof(input->type), QSORT_F); \
  }
void
gal_statistics_sort_increasing(gal_data_t *input)
{
  /* Do the sorting. */
  if(input->size)
    switch(input->type)
      {
      case GAL_TYPE_UINT8:
        STATISTICS_SORT(gal_qsort_uint8_i);    break;
      case GAL_TYPE_INT8:
        STATISTICS_SORT(gal_qsort_int8_i);     break;
      case GAL_TYPE_UINT16:
        STATISTICS_SORT(gal_qsort_uint16_i);   break;
      case GAL_TYPE_INT16:
        STATISTICS_SORT(gal_qsort_int16_i);    break;
      case GAL_TYPE_UINT32:
        STATISTICS_SORT(gal_qsort_uint32_i);   break;
      case GAL_TYPE_INT32:
        STATISTICS_SORT(gal_qsort_int32_i);    break;
      case GAL_TYPE_UINT64:
        STATISTICS_SORT(gal_qsort_uint64_i);   break;
      case GAL_TYPE_INT64:
        STATISTICS_SORT(gal_qsort_int64_i);    break;
      case GAL_TYPE_FLOAT32:
        STATISTICS_SORT(gal_qsort_float32_i);  break;
      case GAL_TYPE_FLOAT64:
        STATISTICS_SORT(gal_qsort_float64_i);  break;
      default:
        error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
              __func__, input->type);
      }

  /* Set the flags. */
  input->flag |=  GAL_DATA_FLAG_SORT_CH;
  input->flag |=  GAL_DATA_FLAG_SORTED_I;
  input->flag &= ~GAL_DATA_FLAG_SORTED_D;
}





/* See explanations above 'gal_statistics_sort_increasing'. */
void
gal_statistics_sort_decreasing(gal_data_t *input)
{
  /* Do the sorting. */
  if(input->size)
    switch(input->type)
      {
      case GAL_TYPE_UINT8:
        STATISTICS_SORT(gal_qsort_uint8_d);    break;
      case GAL_TYPE_INT8:
        STATISTICS_SORT(gal_qsort_int8_d);     break;
      case GAL_TYPE_UINT16:
        STATISTICS_SORT(gal_qsort_uint16_d);   break;
      case GAL_TYPE_INT16:
        STATISTICS_SORT(gal_qsort_int16_d);    break;
      case GAL_TYPE_UINT32:
        STATISTICS_SORT(gal_qsort_uint32_d);   break;
      case GAL_TYPE_INT32:
        STATISTICS_SORT(gal_qsort_int32_d);    break;
      case GAL_TYPE_UINT64:
        STATISTICS_SORT(gal_qsort_uint64_d);   break;
      case GAL_TYPE_INT64:
        STATISTICS_SORT(gal_qsort_int64_d);    break;
      case GAL_TYPE_FLOAT32:
        STATISTICS_SORT(gal_qsort_float32_d);  break;
      case GAL_TYPE_FLOAT64:
        STATISTICS_SORT(gal_qsort_float64_d);  break;
      default:
        error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
              __func__, input->type);
      }

  /* Set the flags. */
  input->flag |=  GAL_DATA_FLAG_SORT_CH;
  input->flag |=  GAL_DATA_FLAG_SORTED_D;
  input->flag &= ~GAL_DATA_FLAG_SORTED_I;
}





/* Return a dataset that doesn't have blank values and is sorted. If the
   'inplace' value is set to 1, then the input array will be modified,
   otherwise, a new array will be allocated with the desired properties. So
   if it is already sorted and has blank values, the 'inplace' variable is
   irrelevant.

   This function can also work on tiles, in that case, 'inplace' is
   useless, because a tile doesn't own its dataset and the dataset is not
   contiguous. */
gal_data_t *
gal_statistics_no_blank_sorted(gal_data_t *input, int inplace)
{
  gal_data_t *contig, *noblank, *sorted;

  /* We need to account for the case that there are no elements in the
     input. */
  if(input->size)
    {
      /* If this is a tile, then first we have to copy it into a contiguous
         piece of memory. After this step, we will only be dealing with
         'contig' (for a contiguous patch of memory). */
      if(input->block)
        {
          /* Copy the input into a contiguous patch of memory. */
          contig=gal_data_copy(input);

          /* When the data was a tile, we have already copied the array
             into a separate allocated space. So to avoid any further
             copying, we will just set the 'inplace' variable to 1. */
          inplace=1;
        }
      else contig=input;


      /* Make sure there are no blanks in the array that will be
         used. After this step, we won't be dealing with 'input' any more,
         but with 'noblank'. */
      if( gal_blank_present(contig, 1) )
        {
          /* See if we should allocate a new dataset to remove blanks or if
             we can use the actual contiguous patch of memory. */
          noblank = inplace ? contig : gal_data_copy(contig);
          gal_blank_remove(noblank);
        }
      else noblank=contig;

      /* Make sure the array is sorted. After this step, we won't be
         dealing with 'noblank' any more but with 'sorted'. */
      if(noblank->size)
        {
          if( gal_statistics_is_sorted(noblank, 1) )
            sorted = inplace ? noblank : gal_data_copy(noblank);
          else
            {
              if(inplace) sorted=noblank;
              else
                {
                  if(noblank!=input)   /* no-blank is already allocated. */
                    sorted=noblank;
                  else
                    sorted=gal_data_copy(noblank);
                }
              gal_statistics_sort_increasing(sorted);
            }
        }
      else
        sorted=noblank;
    }

  /* Input's size was zero. Note that we cannot simply copy the zero-sized
     input dataset, we'll have to allocate it here. */
  else
    sorted = ( inplace
               ? input
               : gal_data_alloc(NULL, input->type, 0, NULL, input->wcs, 0,
                                input->minmapsize, input->quietmmap,
                                NULL, NULL, NULL) );

  /* Set the blank and sorted flags if the dataset has zero-elements. Even
     if having blank values or being sorted is not defined on a
     zero-element dataset, it is up to different functions to choose what
     they will do with a zero-element dataset. The flags have to be set
     after this function any way. */
  if(sorted->size==0)
    {
      sorted->flag |= GAL_DATA_FLAG_SORT_CH;
      sorted->flag |= GAL_DATA_FLAG_BLANK_CH;
      sorted->flag |= GAL_DATA_FLAG_SORTED_I;
      sorted->flag &= ~GAL_DATA_FLAG_HASBLANK;
      sorted->flag &= ~GAL_DATA_FLAG_SORTED_D;
    }

  /* Return final array. */
  return sorted;
}




















/****************************************************************
 ********     Histogram and Cumulative Frequency Plot     *******
 ****************************************************************/
/* Generate an array of regularly spaced elements.

   Input arguments:

     * The 'input' set you want to apply the bins to. This is only
       necessary if the range argument is not complete, see below. If
       'range' has all the necessary information, you can pass a NULL
       pointer for 'input'.

     * The 'inrange' data structure keeps the desired range along each
       dimension of the input data structure, it has to be in float32
       type. Note that if

         - If you want the full range of the dataset (in any dimensions,
           then just set 'range' to NULL and the range will be specified
           from the minimum and maximum value of the dataset.

         - If there is one element for each dimension in range, then it is
           viewed as a quantile (Q), and the range will be: 'Q to 1-Q'.

         - If there are two elements for each dimension in range, then they
           are assumed to be your desired minimum and maximum values. When
           either of the two are NaN, the minimum and maximum will be
           calculated for it.

     * The number of bins: must be larger than 0.

     * 'onebinstart' A desired value for onebinstart. Note that with this
        option, the bins won't start and end exactly on the given range
        values, it will be slightly shifted to accommodate this
        request.

  The output is a 1D array (column) of type double, it has to be double to
  account for small differences on the bin edges.
*/
gal_data_t *
gal_statistics_regular_bins(gal_data_t *input, gal_data_t *inrange,
                            size_t numbins, double onebinstart)
{
  size_t i;
  gal_data_t *bins, *tmp, *range;
  double *b, *ra, min=NAN, max=NAN, hbw, diff, binwidth;


  /* Some sanity checks. */
  if(numbins==0)
    error(EXIT_FAILURE, 0, "%s: 'numbins' cannot be given a value of 0",
          __func__);
  if(input->size==0)
    error(EXIT_FAILURE, 0, "%s: input's size is 0", __func__);


  /* Set the minimum and maximum values. */
  if(inrange && inrange->size)
    {
      /* Make sure we are dealing with a double type range. */
      if(inrange->type==GAL_TYPE_FLOAT64)
        range=inrange;
      else
        range=gal_data_copy_to_new_type(inrange, GAL_TYPE_FLOAT64);

      /* Set the minimum and maximum of the bins. */
      ra=range->array;
      if( (range->size)%2 )
        error(EXIT_FAILURE, 0, "%s: quantile ranges are not implemented yet",
              __func__);
      else
        {
          /* If the minimum isn't set (is blank), find it. */
          if( isnan(ra[0]) )
            {
              tmp=gal_data_copy_to_new_type_free(
                            gal_statistics_minimum(input), GAL_TYPE_FLOAT64);
              min=*((double *)(tmp->array));
              gal_data_free(tmp);
            }
          else min=ra[0];

          /* For the maximum, when it isn't set, we'll add a very small
             value, so all points are included. */
          if( isnan(ra[1]) )
            {
              tmp=gal_data_copy_to_new_type_free(gal_statistics_maximum(input),
                                                 GAL_TYPE_FLOAT64);
              max=*((double *)(tmp->array));

              /* Clean up. */
              gal_data_free(tmp);
            }
          else max=ra[1];
        }

      /* Clean up: if 'range' was allocated. */
      if(range!=inrange) gal_data_free(range);
    }
  /* No range was given, find the minimum and maximum. */
  else
    {
      tmp=gal_data_copy_to_new_type_free(gal_statistics_minimum(input),
                                         GAL_TYPE_FLOAT64);
      min=*((double *)(tmp->array));
      gal_data_free(tmp);
      tmp=gal_data_copy_to_new_type_free(gal_statistics_maximum(input),
                                         GAL_TYPE_FLOAT64);
      max=*((double *)(tmp->array));

      /* Clean up. */
      gal_data_free(tmp);
    }


  /* Allocate the space for the bins. */
  bins=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numbins, NULL,
                      0, input->minmapsize, input->quietmmap, "bin_center",
                      input->unit, "Center value of each bin.");


  /* Set central bin values. */
  b=bins->array;
  hbw = ( binwidth=(max-min)/numbins )/2;
  for(i=0;i<numbins;++i) b[i] = min + i*binwidth + hbw;


  /* Go over all the bins and stop when the sign of the two sides
     of one bin are different. */
  if( !isnan(onebinstart) )
    {
      for(i=0;i<numbins-1;++i)
        if( (b[i]-hbw) < onebinstart && (b[i+1]-hbw) > onebinstart) break;
      if( i != numbins-1 )
        {
          diff = onebinstart - (b[i]-hbw);
          for(i=0;i<numbins;++i)
            b[i] += diff;
        }
    }

  /* For a check:
  printf("min: %g\n", min);
  printf("max: %g\n", max);
  printf("onebinstart: %.10f\n", onebinstart);
  printf("binwidth: %g\n", binwidth);
  for(i=0;i<numbins;++i)
    printf("%zu: %.4g\t(%g, %g)\n", i, b[i], b[i]-hbw, b[i]+hbw);
  */

  /* Set the status of the bins to regular and return. */
  bins->status=GAL_STATISTICS_BINS_REGULAR;
  return bins;
}





/* Make a histogram of all the elements in the given dataset with bin
   values that are defined in the 'inbins' structure (see
   'gal_statistics_regular_bins'). 'inbins' is not mandatory, if you pass a
   NULL pointer, the bins structure will be built within this function
   based on the 'numbins' input. As a result, when you have already defined
   the bins, 'numbins' is not used. */

#define HISTOGRAM_TYPESET(IT) {                                         \
    IT *a=input->array, *af=a+input->size;                              \
    do                                                                  \
      if(*a>=min && *a<=max)                                            \
        {                                                               \
          h_i=(*a-min)/binwidth;                                        \
          /* When '*a' is the largest element (within floating point */ \
          /* errors), 'h_i' can be one element larger than the       */ \
          /* number of bins. But since its in the dataset, we need   */ \
          /* to count it. So we'll put it in the last bin.           */ \
          ++h[ h_i - (h_i==hist->size ? 1 : 0) ];                       \
        }                                                               \
    while(++a<af);                                                      \
  }

gal_data_t *
gal_statistics_histogram(gal_data_t *input, gal_data_t *bins, int normalize,
                         int maxone)
{
  float *f, *ff;
  size_t *h, h_i;
  gal_data_t *hist;
  double *d, min, max, ref=NAN, binwidth;


  /* Check if the bins are regular or not. For irregular bins, we can
     either use the old implementation, or GSL's histogram
     functionality. */
  if(bins==NULL)
    error(EXIT_FAILURE, 0, "%s: 'bins' is NULL", __func__);
  if(bins->status!=GAL_STATISTICS_BINS_REGULAR)
    error(EXIT_FAILURE, 0, "%s: the input bins are not regular. Currently "
          "it is only implemented for regular bins", __func__);
  if(input->size==0)
    error(EXIT_FAILURE, 0, "%s: input's size is 0", __func__);


  /* Check if normalize and 'maxone' are not called together. */
  if(normalize && maxone)
    error(EXIT_FAILURE, 0, "%s: only one of 'normalize' and 'maxone' may "
          "be given", __func__);


  /* Allocate the histogram (note that we are clearning it so all values
     are zero. */
  hist=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, bins->ndim, bins->dsize,
                      NULL, 1, input->minmapsize, input->quietmmap,
                      "hist_number", "counts",
                      "Number of data points within each bin.");


  /* Set the minimum and maximum range of the histogram from the bins. */
  d=bins->array;
  binwidth=d[1]-d[0];
  min = d[ 0      ] - binwidth/2;
  max = d[ bins->size-1 ] + binwidth/2;


  /* Go through all the elements and find out which bin they belong to. */
  h=hist->array;
  switch(input->type)
    {
    case GAL_TYPE_UINT8:     HISTOGRAM_TYPESET(uint8_t);     break;
    case GAL_TYPE_INT8:      HISTOGRAM_TYPESET(int8_t);      break;
    case GAL_TYPE_UINT16:    HISTOGRAM_TYPESET(uint16_t);    break;
    case GAL_TYPE_INT16:     HISTOGRAM_TYPESET(int16_t);     break;
    case GAL_TYPE_UINT32:    HISTOGRAM_TYPESET(uint32_t);    break;
    case GAL_TYPE_INT32:     HISTOGRAM_TYPESET(int32_t);     break;
    case GAL_TYPE_UINT64:    HISTOGRAM_TYPESET(uint64_t);    break;
    case GAL_TYPE_INT64:     HISTOGRAM_TYPESET(int64_t);     break;
    case GAL_TYPE_FLOAT32:   HISTOGRAM_TYPESET(float);       break;
    case GAL_TYPE_FLOAT64:   HISTOGRAM_TYPESET(double);      break;
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, input->type);
    }


  /* For a check:
  {
    size_t i, *hh=hist->array;
    for(i=0;i<hist->size;++i) printf("%-10.4f%zu\n", f[i], hh[i]);
  }
  */


  /* Find the reference to correct the histogram if necessary. */
  if(normalize)
    {
      /* Set the reference. */
      ref=0.0f;
      hist=gal_data_copy_to_new_type_free(hist, GAL_TYPE_FLOAT32);
      ff=(f=hist->array)+hist->size; do ref += *f++;   while(f<ff);

      /* Correct the name, units and comments. */
      free(hist->name); free(hist->unit); free(hist->comment);
      gal_checkset_allocate_copy("hist_normalized", &hist->name);
      gal_checkset_allocate_copy("frac", &hist->unit);
      gal_checkset_allocate_copy("Normalized histogram value for this bin.",
                                 &hist->comment);
    }
  if(maxone)
    {
      /* Calculate the reference. */
      ref=-FLT_MAX;
      hist=gal_data_copy_to_new_type_free(hist, GAL_TYPE_FLOAT32);
      ff=(f=hist->array)+hist->size;
      do ref = *f>ref ? *f : ref; while(++f<ff);

      /* Correct the name, units and comments. */
      free(hist->name); free(hist->unit); free(hist->comment);
      gal_checkset_allocate_copy("hist_maxone", &hist->name);
      gal_checkset_allocate_copy("frac", &hist->unit);
      gal_checkset_allocate_copy("Fractional histogram value for this bin "
                                 "when maximum bin value is 1.0.",
                                 &hist->comment);
    }


  /* Correct the histogram if necessary. */
  if( !isnan(ref) )
    { ff=(f=hist->array)+hist->size; do *f++ /= ref;   while(f<ff); }


  /* Return the histogram. */
  return hist;
}





/* Make a cumulative frequency plot (CFP) of all the elements in the given
   dataset with bin values that are defined in the 'bins' structure (see
   'gal_statistics_regular_bins').

   The CFP is built from the histogram: in each bin, the value is the sum
   of all previous bins in the histogram. Thus, if you have already
   calculated the histogram before calling this function, you can pass it
   onto this function as the data structure in 'bins->next'. If
   'bin->next!=NULL', then it is assumed to be the histogram. If it is
   NULL, then the histogram will be calculated internally and freed after
   the job is finished.

   When a histogram is given and it is normalized, the CFP will also be
   normalized (even if the normalized flag is not set here): note that a
   normalized CFP's maximum value is 1. */
gal_data_t *
gal_statistics_cfp(gal_data_t *input, gal_data_t *bins, int normalize)
{
  double sum;
  float *f, *ff, *hf;
  gal_data_t *hist, *cfp;
  size_t *s, *sf, *hs, sums;


  /* Check if the bins are regular or not. For irregular bins, we can
     either use the old implementation, or GSL's histogram
     functionality. */
  if(bins->status!=GAL_STATISTICS_BINS_REGULAR)
    error(EXIT_FAILURE, 0, "%s: the input bins are not regular. Currently "
          "it is only implemented for regular bins", __func__);
  if(input->size==0)
    error(EXIT_FAILURE, 0, "%s: input's size is 0", __func__);


  /* Prepare the histogram. */
  hist = ( bins->next
           ? bins->next
           : gal_statistics_histogram(input, bins, 0, 0) );


  /* If the histogram has float32 type it was given by the user and is
     either normalized or its maximum was set to 1. We can only use it if
     it was normalized. If it isn't normalized, then we must ignore it and
     build the histogram here.*/
  if(hist->type==GAL_TYPE_FLOAT32)
    {
      sum=0.0f;
      ff=(f=hist->array)+hist->size; do sum += *f++;   while(f<ff);
      if(sum!=1.0f)
        hist=gal_statistics_histogram(input, bins, 0, 0);
    }


  /* Allocate the cumulative frequency plot's necessary space. */
  cfp=gal_data_alloc( NULL, hist->type, bins->ndim, bins->dsize,
                      NULL, 1, input->minmapsize, input->quietmmap,
                      ( hist->type==GAL_TYPE_FLOAT32
                        ? "cfp_normalized" : "cfp_number" ),
                      ( hist->type==GAL_TYPE_FLOAT32
                        ? "frac" : "count" ),
                      ( hist->type==GAL_TYPE_FLOAT32
                        ? "Fraction of data elements from the start to this "
                        "bin (inclusive)."
                        : "Number of data elements from the start to this "
                        "bin (inclusive).") );


  /* Fill in the cumulative frequency plot. */
  switch(hist->type)
    {
    case GAL_TYPE_SIZE_T:
      sums=0; hs=hist->array; sf=(s=cfp->array)+cfp->size;
      do sums = (*s += *hs++ + sums); while(++s<sf);
      break;

    case GAL_TYPE_FLOAT32:
      sum=0.0f; hf=hist->array; ff=(f=cfp->array)+cfp->size;
      do sum = (*f += *hf++ + sum);  while(++f<ff);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, cfp->type);
    }


  /* Normalize the CFP if the user asked for it and it wasn't normalized
     until now. */
  if(normalize && cfp->type==GAL_TYPE_SIZE_T)
    {
      /* Find the sum, then divide the plot by it. Note that the sum must
         come from the histogram, not the CFP!*/
      sums=0;
      cfp=gal_data_copy_to_new_type_free(cfp, GAL_TYPE_FLOAT32);
      sf=(s=hist->array)+hist->size; do sums += *s++;   while(s<sf);
      ff=(f=cfp->array)+cfp->size;   do *f++ /= sums;   while(f<ff);

      /* Correct the name, units and comments. */
      free(cfp->name); free(cfp->unit); free(cfp->comment);
      gal_checkset_allocate_copy("cfp_normalized", &cfp->name);
      gal_checkset_allocate_copy("frac", &cfp->unit);
      gal_checkset_allocate_copy("Fraction of data elements from the start "
                                 "to this bin (inclusive).", &cfp->comment);
    }

  /* If the histogram was allocated here, free it. */
  if(hist!=bins->next) gal_data_free(hist);
  return cfp;
}




















/****************************************************************
 *****************         Outliers          ********************
 ****************************************************************/
/* Sigma-cilp a given distribution:

   Inputs:

     - 'multip': multiple of the standard deviation,

     - 'param' must be positive and determines the type of clipping:

         - param<1.0: interpretted as a tolerance level to stop clipping.

         - param>=1.0 and an integer: a specific number of times to do the
           clippping.

   Output elements (type FLOAT32):

     - 0: Number of points used.
     - 1: Median.
     - 2: Mean.
     - 3: Standard deviation.

  The way this function works is very simple: first it will sort the input
  (if it isn't sorted). Afterwards, it will recursively change the starting
  point of the array and its size, calcluating the basic statistics in each
  round to define the new starting point and size.
*/
#define SIGCLIP(IT) {                                                   \
    IT *a  = nbs->array, *af = a  + nbs->size;                          \
    IT *bf = nbs->array, *b  = bf + nbs->size - 1;                      \
                                                                        \
    /* Remove all out-of-range elements from the start of the array. */ \
    if( nbs->flag & GAL_DATA_FLAG_SORTED_I )                            \
      do if( *a > (*med - (multip * *std)) )                            \
           { start=a; break; }                                          \
      while(++a<af);                                                    \
    else                                                                \
      do if( *a < (*med + (multip * *std)) )                            \
           { start=a; break; }                                          \
      while(++a<af);                                                    \
                                                                        \
    /* Remove all out-of-range elements from the end of the array. */   \
    if( nbs->flag & GAL_DATA_FLAG_SORTED_I )                            \
      do if( *b < (*med + (multip * *std)) )                            \
           { size=b-a+1; break; }                                       \
      while(--b>=bf);                                                   \
    else                                                                \
      do if( *b > (*med - (multip * *std)) )                            \
           { size=b-a+1; break; }                                       \
      while(--b>=bf);                                                   \
  }

gal_data_t *
gal_statistics_sigma_clip(gal_data_t *input, float multip, float param,
                          int inplace, int quiet)
{
  float *oa;
  void *start, *nbs_array;
  double *med, *mean, *std;
  uint8_t type=gal_tile_block(input)->type;
  uint8_t bytolerance = param>=1.0f ? 0 : 1;
  double oldmed=NAN, oldmean=NAN, oldstd=NAN;
  size_t num=0, one=1, four=4, size, oldsize;
  gal_data_t *fcopy, *median_i, *median_d, *out, *meanstd;
  gal_data_t *nbs=gal_statistics_no_blank_sorted(input, inplace);
  size_t maxnum = param>=1.0f ? param : GAL_STATISTICS_SIG_CLIP_MAX_CONVERGE;


  /* Some sanity checks. */
  if( multip<=0 )
    error(EXIT_FAILURE, 0, "%s: 'multip', must be greater than zero. The "
          "given value was %g", __func__, multip);
  if( param<=0 )
    error(EXIT_FAILURE, 0, "%s: 'param', must be greater than zero. The "
          "given value was %g", __func__, param);
  if( param >= 1.0f && ceil(param) != param )
    error(EXIT_FAILURE, 0, "%s: when 'param' is larger than 1.0, it is "
          "interpretted as an absolute number of clips. So it must be an "
          "integer. However, your given value %g", __func__, param);
  if( (nbs->flag & GAL_DATA_FLAG_SORT_CH)==0 )
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
          "problem. 'nbs->flag', doesn't have the 'GAL_DATA_FLAG_SORT_CH' "
          "bit activated", __func__, PACKAGE_BUGREPORT);
  if( (nbs->flag & GAL_DATA_FLAG_SORTED_I)==0
      && (nbs->flag & GAL_DATA_FLAG_SORTED_D)==0 )
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
          "problem. 'nbs' isn't sorted", __func__, PACKAGE_BUGREPORT);


  /* Allocate the necessary spaces. */
  out=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &four, NULL, 0,
                     input->minmapsize, input->quietmmap, NULL, NULL, NULL);
  median_i=gal_data_alloc(NULL, type, 1, &one, NULL, 0, input->minmapsize,
                          input->quietmmap, NULL, NULL, NULL);


  /* Only continue processing if we have non-blank elements. */
  oa=out->array;
  nbs_array=nbs->array;
  switch(nbs->size)
    {
    /* There was nothing in the input! */
    case 0:
      if(!quiet)
        printf("NO SIGMA-CLIPPING: all input elements are blank or input's "
               "size is zero.\n");
      oa[0] = oa[1] = oa[2] = oa[3] = NAN;
      break;

    /* Only one element, convert it to floating point and put it as the
       mean and median (the standard deviation will be zero by
       definition). */
    case 1:
      /* Write the values. */
      fcopy=gal_data_copy_to_new_type(nbs, GAL_TYPE_FLOAT32);
      oa[0] = 1;
      oa[1] = *((float *)(fcopy->array));
      oa[2] = *((float *)(fcopy->array));
      oa[3] = 0;
      gal_data_free(fcopy);

      /* Print the comments if requested. */
      if(!quiet)
        {
          printf("%-8s %-10s %-15s %-15s %-15s\n",
                 "round", "number", "median", "mean", "STD");
          printf("%-8d %-10.0f %-15g %-15g %-15g\n",
                 0, oa[0], oa[1], oa[2], oa[3]);
        }
      break;

    /* More than one element. */
    default:
      /* Print the comments if requested. */
      if(!quiet)
        printf("%-8s %-10s %-15s %-15s %-15s\n",
               "round", "number", "median", "mean", "STD");

      /* Do the clipping, but first initialize the values that will be
         changed during the clipping: the start of the array and the
         array's size. */
      size=nbs->size;
      start=nbs->array;
      while(num<maxnum && size)
        {
          /* Find the average and Standard deviation, note that both
             'start' and 'size' will be different in the next round. */
          nbs->array = start;
          nbs->size = oldsize = size;

          /* For a detailed check, just correct the type).
          if(!quiet)
            {
              size_t iii;
              printf("nbs->size: %zu\n", nbs->size);
              for(iii=0;iii<nbs->size;++iii)
                printf("%f\n", ((float *)(nbs->array))[iii]);
            }
          */

          /* Find the mean, median and standard deviation. */
          meanstd=gal_statistics_mean_std(nbs);
          statistics_median_in_sorted_no_blank(nbs, median_i->array);
          median_d=gal_data_copy_to_new_type(median_i, GAL_TYPE_FLOAT64);

          /* Put them in usable (with a type) pointers. */
          mean = meanstd->array;
          med  = median_d->array;
          std  = &((double *)(meanstd->array))[1];

          /* If the user wanted to view the steps, show it to them. */
          if(!quiet)
            printf("%-8zu %-10zu %-15g %-15g %-15g\n",
                   num+1, size, *med, *mean, *std);

          /* If we are to work by tolerance, then check if we should jump
             out of the loop. Normally, 'oldstd' should be larger than std,
             because the possible outliers have been removed. If it is not,
             it means that we have clipped too much and must stop anyway,
             so we don't need an absolute value on the difference!

             Note that when all the elements are identical after the clip,
             'std' will be zero. In this case we shouldn't calculate the
             tolerance (because it will be infinity and thus lager than the
             requested tolerance level value).*/
          if( bytolerance && num>0 )
            if( *std==0 || ((oldstd - *std) / *std) < param )
              {
                if(*std==0) {oldmed=*med; oldstd=*std; oldmean=*mean;}
                gal_data_free(meanstd);   gal_data_free(median_d);
                break;
              }

          /* Clip all the elements outside of the desired range: since the
             array is sorted, this means to just change the starting
             pointer and size of the array. */
          switch(type)
            {
            case GAL_TYPE_UINT8:     SIGCLIP( uint8_t  );   break;
            case GAL_TYPE_INT8:      SIGCLIP( int8_t   );   break;
            case GAL_TYPE_UINT16:    SIGCLIP( uint16_t );   break;
            case GAL_TYPE_INT16:     SIGCLIP( int16_t  );   break;
            case GAL_TYPE_UINT32:    SIGCLIP( uint32_t );   break;
            case GAL_TYPE_INT32:     SIGCLIP( int32_t  );   break;
            case GAL_TYPE_UINT64:    SIGCLIP( uint64_t );   break;
            case GAL_TYPE_INT64:     SIGCLIP( int64_t  );   break;
            case GAL_TYPE_FLOAT32:   SIGCLIP( float    );   break;
            case GAL_TYPE_FLOAT64:   SIGCLIP( double   );   break;
            default:
              error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
                    __func__, type);
            }

          /* Set the values from this round in the old elements, so the
             next round can compare with, and return then if necessary. */
          oldmed =  *med;
          oldstd  = *std;
          oldmean = *mean;
          ++num;

          /* Clean up: */
          gal_data_free(meanstd);
          gal_data_free(median_d);
        }

      /* If we were in tolerance mode and 'num' and 'maxnum' are equal (the
         loop didn't stop by tolerance), so the outputs should be NaN. */
      out->status=num;
      if( size==0 || (bytolerance && num==maxnum) )
        oa[0] = oa[1] = oa[2] = oa[3] = NAN;
      else
        {
          oa[0] = size;
          oa[1] = oldmed;
          oa[2] = oldmean;
          oa[3] = oldstd;
        }
    }

  /* Clean up and return. */
  nbs->array=nbs_array;
  gal_data_free(median_i);
  if(nbs!=input) gal_data_free(nbs);
  return out;
}





/* Find the first outlier in a distribution. */
#define OUTLIER_BYTYPE(IT) {                                            \
    IT *arr=nbs->array;                                                 \
    for(i=window_size;i<nbs->size;++i)                                  \
      {                                                                 \
        /* Fill in the distance array. */                               \
        for(j=0; j<wtakeone; ++j)                                       \
          darr[j] = arr[i-window_size+j+1] - arr[i-window_size+j];      \
                                                                        \
        /* Get the sigma-clipped information. */                        \
        sclip=gal_statistics_sigma_clip(dist, sigclip_multip,           \
                                        sigclip_param, 0, 1);           \
        sarr=sclip->array;                                              \
                                                                        \
        /* For a check */                                               \
        if(quiet==0)                                                    \
          printf("%f [%zu]: %f (%f, %f) %f\n", (float)(arr[i]), i,      \
                 (float)(arr[i]-arr[i-1]), sarr[1], sarr[3],            \
                 (((double)(arr[i]-arr[i-1])) - sarr[1])/sarr[3]);      \
                                                                        \
        /* Terminate the loop if the dist. is larger than requested. */ \
        /* This shows we have reached the first outlier's position. */  \
        if( (((double)(arr[i]-arr[i-1])) - sarr[1]) > sigma*sarr[3] )   \
          {                                                             \
            /* Allocate the output dataset. */                          \
            out=gal_data_alloc(NULL, input->type, 1, &one, NULL, 0, -1, \
                               1, NULL, NULL, NULL);                    \
                                                                        \
            /* Write the outlier, clean up and break. */                \
            *(IT *)(out->array)=arr[i-1];                               \
            gal_data_free(sclip);                                       \
            break;                                                      \
          }                                                             \
                                                                        \
        /* Clean up (if we get here). */                                \
        gal_data_free(sclip);                                           \
      }                                                                 \
  }
gal_data_t *
gal_statistics_outlier_positive(gal_data_t *input, size_t window_size,
                                float sigma, float sigclip_multip,
                                float sigclip_param, int inplace, int quiet)
{
  float *sarr;
  double *darr;
  size_t i, j, one=1, wtakeone;
  gal_data_t *dist, *sclip, *nbs, *out=NULL;

  /* Remove all blanks and sort the dataset. */
  nbs=gal_statistics_no_blank_sorted(input, inplace);

  /* If all elements are blank, simply return the default (NULL) output. */
  if(nbs->size==0) return out;

  /* Only continue if the window size is more than 2 elements (out
     "outlier" is hard to define on smaller datasets). */
  if(window_size>2)
    {
      /* For a check.
      if(nbs->type==GAL_TYPE_FLOAT32)
        {
          float *n=nbs->array;
          for(i=0;i<nbs->size;++i)
            printf("%f\n", n[i]);
          exit(0);
        }
      */


      /* Allocate space to keep the distances. */
      wtakeone=window_size-1;
      dist=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &wtakeone, NULL,
                          0, -1, 1, NULL, NULL, NULL);
      darr=dist->array;

      /* Find the outlier based on the type of the input dataset. */
      switch(input->type)
        {
        case GAL_TYPE_UINT8:     OUTLIER_BYTYPE( uint8_t  );   break;
        case GAL_TYPE_INT8:      OUTLIER_BYTYPE( int8_t   );   break;
        case GAL_TYPE_UINT16:    OUTLIER_BYTYPE( uint16_t );   break;
        case GAL_TYPE_INT16:     OUTLIER_BYTYPE( int16_t  );   break;
        case GAL_TYPE_UINT32:    OUTLIER_BYTYPE( uint32_t );   break;
        case GAL_TYPE_INT32:     OUTLIER_BYTYPE( int32_t  );   break;
        case GAL_TYPE_UINT64:    OUTLIER_BYTYPE( uint64_t );   break;
        case GAL_TYPE_INT64:     OUTLIER_BYTYPE( int64_t  );   break;
        case GAL_TYPE_FLOAT32:   OUTLIER_BYTYPE( float    );   break;
        case GAL_TYPE_FLOAT64:   OUTLIER_BYTYPE( double   );   break;
        default:
          error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
                __func__, input->type);
        }

      /* Clean up. */
      gal_data_free(dist);
    }

  /* Clean up and return. */
  if(nbs!=input) gal_data_free(nbs);
  return out;
}






/* Find the outliers using the average distance of the neighboring
   points. */
#define OUTLIER_FLAT_CFP_BYTYPE(IT) {                                   \
    IT diff, *pr=prev->array;                                           \
    IT *a=nbs->array, *p=a+d, *pp=a+nbs->size-d;                        \
                                                                        \
    do                                                                  \
      {                                                                 \
        diff=*(p+d)-*(p-d);                                             \
        if(p-a-d<numprev)                                               \
          {                                                             \
            pr[p-a-d]=diff;                                             \
            if(!quiet) printf("%-6zu%-15g%-15g\n", p-a, (float)(*p),    \
                              (float)diff);                             \
          }                                                             \
        else                                                            \
          {                                                             \
            /* Sigma-clipped median and std for a check. */             \
            prev->flag=0;                                               \
            prev->size=prev->dsize[0]=numprev;                          \
            sclip=gal_statistics_sigma_clip(prev, sigclip_multip,       \
                                            sigclip_param, 1, 1);       \
                                                                        \
            sarr=sclip->array;                                          \
            check = (diff - sarr[1]) / sarr[3];                         \
                                                                        \
            /* If requested, print the values. */                       \
            if(!quiet) printf("%-6zu%-15g%-15g%-15g (%g,%g)\n", p-a,    \
                              (float)(*p), (float)diff, check, sarr[1], \
                              sarr[3]);                                 \
                                                                        \
            /* When values are equal, std will be roughly zero */       \
            if(sarr[3]>1e-6 && check>thresh)                            \
              {                                                         \
                if(flatind==GAL_BLANK_SIZE_T)                           \
                  {                                                     \
                    ++counter;                                          \
                    flatind=p-a;                                        \
                  }                                                     \
                else                                                    \
                  {                                                     \
                    if(flatind==p-a-counter)                            \
                      { /* First element above thresh is 0, so for */   \
                        /* counting, when counting the number of */     \
                        /* contiguous elements, we have to add 1. */    \
                        if(counter+1==numcontig)                        \
                          {gal_data_free(sclip); break;}                \
                        else ++counter;                                 \
                      }                                                 \
                    else { flatind=GAL_BLANK_SIZE_T; counter=0; }       \
                  }                                                     \
              }                                                         \
            else { flatind=GAL_BLANK_SIZE_T; counter=0; }               \
            pr[(p-a-d)%numprev]=diff;                                   \
            gal_data_free(sclip);                                       \
          }                                                             \
      }                                                                 \
    while(++p<pp);                                                      \
    if(counter+1!=numcontig) flatind=GAL_BLANK_SIZE_T;                    \
  }

gal_data_t *
gal_statistics_outlier_flat_cfp(gal_data_t *input, size_t numprev,
                                float sigclip_multip, float sigclip_param,
                                float thresh, size_t numcontig, int inplace,
                                int quiet, size_t *index)
{
  float *sarr;
  double check;
  gal_data_t  *nbs, *prev, *out=NULL, *sclip;
  size_t d=2, counter=0, one=1, flatind=GAL_BLANK_SIZE_T;

  /* Sanity checks. */
  if(thresh<=0)
    error(EXIT_FAILURE, 0, "%s: the value of 'thresh' (%g) must be "
          "positive", __func__, thresh);
  if(numprev==0)
    error(EXIT_FAILURE, 0, "%s: 'numprev' (%zu) cannot be zero", __func__,
          numprev);

  /* Remove all blanks and sort the dataset. */
  nbs=gal_statistics_no_blank_sorted(input, inplace);

  /* Keep previous slopes. */
  prev=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numprev, NULL, 0, -1,
                      1, NULL, NULL, NULL);

  /* Find the index where the distribution becomes sufficiently flat. */
  switch(nbs->type)
    {
    case GAL_TYPE_UINT8:   OUTLIER_FLAT_CFP_BYTYPE( uint8_t  ); break;
    case GAL_TYPE_INT8:    OUTLIER_FLAT_CFP_BYTYPE( int8_t   ); break;
    case GAL_TYPE_UINT16:  OUTLIER_FLAT_CFP_BYTYPE( uint16_t ); break;
    case GAL_TYPE_INT16:   OUTLIER_FLAT_CFP_BYTYPE( int16_t  ); break;
    case GAL_TYPE_UINT32:  OUTLIER_FLAT_CFP_BYTYPE( uint32_t ); break;
    case GAL_TYPE_INT32:   OUTLIER_FLAT_CFP_BYTYPE( int32_t  ); break;
    case GAL_TYPE_UINT64:  OUTLIER_FLAT_CFP_BYTYPE( uint64_t ); break;
    case GAL_TYPE_INT64:   OUTLIER_FLAT_CFP_BYTYPE( int64_t  ); break;
    case GAL_TYPE_FLOAT32: OUTLIER_FLAT_CFP_BYTYPE( float    ); break;
    case GAL_TYPE_FLOAT64: OUTLIER_FLAT_CFP_BYTYPE( double   ); break;
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, nbs->type);
    }

  /* Write the output dataset: if no point flat part was found, return
     NULL. */
  if(flatind!=GAL_BLANK_SIZE_T)
    {
      out=gal_data_alloc(NULL, input->type, 1, &one, NULL, 0, -1, 1,
                         NULL, NULL, NULL);
      memcpy(out->array,
             gal_pointer_increment(nbs->array, flatind, nbs->type),
             gal_type_sizeof(nbs->type));
    }

  /* Clean up and return. */
  if(nbs!=input) gal_data_free(nbs);
  if(index) *index=flatind;
  gal_data_free(prev);
  return out;
}
