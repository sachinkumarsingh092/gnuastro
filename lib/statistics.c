/*********************************************************************
Statistical functions.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gnuastro/blank.h>
#include <gnuastro/qsort.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>

#include <checkset.h>

#include "mode.h"










/****************************************************************
 ********               Simple statistics                 *******
 ****************************************************************/
/* Return the number of non-blank elements in an array as a single element,
   uint64 type data structure. */
#define STATS_NUM(IT) {                                                 \
    IT b, *a=input->array, *af=a+input->size;                           \
    gal_blank_write(&b, input->type);                                   \
    if( gal_blank_present(input) )                                      \
      {                                                                 \
        if(b==b)   /* Blank can be found with the equal sign. */        \
          do if(*a!=b)       ++(*num);            while(++a<af);        \
        else       /* Blank is NaN (fails on equal).  */                \
          do if(*a==*a)     ++(*num);            while(++a<af);         \
      }                                                                 \
    else *num=input->size;                                              \
  }
gal_data_t *
gal_statistics_number(gal_data_t *input)
{
  uint64_t *num;
  size_t dsize=1;
  gal_data_t *out=gal_data_alloc(NULL, GAL_DATA_TYPE_UINT64, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);
  num=out->array;
  switch(input->type)
    {
    case GAL_DATA_TYPE_UINT8:     STATS_NUM( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      STATS_NUM( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    STATS_NUM( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     STATS_NUM( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    STATS_NUM( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     STATS_NUM( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    STATS_NUM( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     STATS_NUM( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   STATS_NUM( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   STATS_NUM( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_number'", input->type);
    }
  return out;
}





/* Return the minimum (non-blank) value of a dataset in the same type as
   the dataset. Note that a NaN (blank in floating point) will fail on any
   comparison. So when finding the minimum or maximum, when the blank value
   is NaN, we can safely assume there is no blank value at all. */
#define STATS_MIN(IT) {                                                 \
    IT b, *o=out->array, *a=input->array, *af=a+input->size;            \
    gal_blank_write(&b, input->type);                                   \
    if( b==b && gal_blank_present(input) )                              \
      do if(*a!=b) { *o = *a < *o ? *a : *o; ++n; } while(++a<af);      \
    else                                                                \
      do           { *o = *a < *o ? *a : *o; ++n; } while(++a<af);      \
                                                                        \
    /* If there were no usable elements, set the output to blank. */    \
    if(n==0) *o=b;                                                      \
  }
gal_data_t *
gal_statistics_minimum(gal_data_t *input)
{
  size_t dsize=1, n=0;
  gal_data_t *out=gal_data_alloc(NULL, input->type, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);
  gal_data_type_max(out->type, out->array);
  switch(input->type)
    {
    case GAL_DATA_TYPE_UINT8:     STATS_MIN( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      STATS_MIN( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    STATS_MIN( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     STATS_MIN( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    STATS_MIN( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     STATS_MIN( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    STATS_MIN( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     STATS_MIN( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   STATS_MIN( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   STATS_MIN( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_minimum'", input->type);
    }
  return out;
}





/* Return the maximum (non-blank) value of a dataset in the same type as
   the dataset. See explanations of `gal_statistics_minimum'. */
#define STATS_MAX(IT) {                                                 \
    IT b, *o=out->array, *a=input->array, *af=a+input->size;            \
    gal_blank_write(&b, input->type);                                   \
    if( b==b && gal_blank_present(input) )                              \
      do if(*a!=b) { *o = *a > *o ? *a : *o; ++n; } while(++a<af);      \
    else                                                                \
      do           { *o = *a > *o ? *a : *o; ++n; } while(++a<af);      \
                                                                        \
    /* If there were no usable elements, set the output to blank. */    \
    if(n==0) *o=b;                                                      \
  }
gal_data_t *
gal_statistics_maximum(gal_data_t *input)
{
  size_t dsize=1, n=0;
  gal_data_t *out=gal_data_alloc(NULL, input->type, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);
  gal_data_type_min(out->type, out->array);
  switch(input->type)
    {
    case GAL_DATA_TYPE_UINT8:     STATS_MAX( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      STATS_MAX( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    STATS_MAX( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     STATS_MAX( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    STATS_MAX( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     STATS_MAX( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    STATS_MAX( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     STATS_MAX( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   STATS_MAX( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   STATS_MAX( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_maximum'", input->type);
    }
  return out;
}





/* Return the sum of the input dataset as a single element dataset of type
   float64. */
#define STATS_SUM_NUM(IT) {                                             \
    IT b, *a=input->array, *af=a+input->size;                           \
    gal_blank_write(&b, input->type);                                   \
    if( gal_blank_present(input) )                                      \
      {                                                                 \
        if(b==b)       /* Blank value can be checked with an equals. */ \
          do if(*a!=b)       { ++n; sum += *a; }  while(++a<af);        \
        else           /* Blank value will fail an equal comparison. */ \
          do if(*a==*a)      { ++n; sum += *a; }  while(++a<af);        \
      }                                                                 \
    else                                                                \
      {                                                                 \
        do                          sum += *a;    while(++a<af);        \
        n = input->size;                                                \
      }                                                                 \
  }
gal_data_t *
gal_statistics_sum(gal_data_t *input)
{
  double sum=0.0f;
  size_t dsize=1, n=0;
  gal_data_t *out=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);
  switch(input->type)
    {
    case GAL_DATA_TYPE_UINT8:     STATS_SUM_NUM( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      STATS_SUM_NUM( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    STATS_SUM_NUM( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     STATS_SUM_NUM( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    STATS_SUM_NUM( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     STATS_SUM_NUM( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    STATS_SUM_NUM( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     STATS_SUM_NUM( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   STATS_SUM_NUM( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   STATS_SUM_NUM( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_maximum'", input->type);
    }
  *((double *)(out->array)) = n ? sum : GAL_BLANK_FLOAT64;
  return out;
}





/* Return the mean of the input dataset as a float64 type single-element
   dataset. */
gal_data_t *
gal_statistics_mean(gal_data_t *input)
{
  double sum=0.0f;
  size_t dsize=1, n=0;
  gal_data_t *out=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);
  switch(input->type)
    {
    case GAL_DATA_TYPE_UINT8:     STATS_SUM_NUM( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      STATS_SUM_NUM( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    STATS_SUM_NUM( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     STATS_SUM_NUM( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    STATS_SUM_NUM( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     STATS_SUM_NUM( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    STATS_SUM_NUM( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     STATS_SUM_NUM( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   STATS_SUM_NUM( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   STATS_SUM_NUM( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_maximum'", input->type);
    }
  *((double *)(out->array)) = n ? sum/n : GAL_BLANK_FLOAT64;
  return out;
}





/* Return the standard deviation of the input dataset as a single element
   dataset of type float64. */
#define STATS_NSS2(IT) {                                                \
    IT b, *a=input->array, *af=a+input->size;                           \
    gal_blank_write(&b, input->type);                                   \
    if( gal_blank_present(input) )                                      \
      {                                                                 \
        if(b==b)       /* Blank value can be checked with an equals. */ \
          do if(*a!=b) { ++n; s += *a; s2 += *a * *a; }  while(++a<af); \
        else           /* Blank value will fail an equal comparison. */ \
          do if(*a==*a){ ++n; s += *a; s2 += *a * *a; }  while(++a<af); \
      }                                                                 \
    else                                                                \
      {                                                                 \
        do             {      s += *a; s2 += *a * *a; }  while(++a<af); \
        n = input->size;                                                \
      }                                                                 \
  }
gal_data_t *
gal_statistics_std(gal_data_t *input)
{
  size_t dsize=1, n=0;
  double s=0.0f, s2=0.0f;
  gal_data_t *out=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);
  switch(input->type)
    {
    case GAL_DATA_TYPE_UINT8:     STATS_NSS2( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      STATS_NSS2( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    STATS_NSS2( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     STATS_NSS2( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    STATS_NSS2( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     STATS_NSS2( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    STATS_NSS2( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     STATS_NSS2( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   STATS_NSS2( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   STATS_NSS2( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_maximum'", input->type);
    }
  *((double *)(out->array)) = n ? sqrt( (s2-s*s/n)/n ) : GAL_BLANK_FLOAT64;
  return out;
}





/* Return the mean and standard deviation of a dataset in one run in type
   float64. The output is a two element data structure, with the first
   value being the mean and the second value the standard deviation. */
gal_data_t *
gal_statistics_mean_std(gal_data_t *input)
{
  size_t dsize=2, n=0;
  double s=0.0f, s2=0.0f;
  gal_data_t *out=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT64, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);
  switch(input->type)
    {
    case GAL_DATA_TYPE_UINT8:     STATS_NSS2( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      STATS_NSS2( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    STATS_NSS2( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     STATS_NSS2( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    STATS_NSS2( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     STATS_NSS2( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    STATS_NSS2( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     STATS_NSS2( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   STATS_NSS2( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   STATS_NSS2( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_maximum'", input->type);
    }
  ((double *)(out->array))[0] = n ? s/n                  : GAL_BLANK_FLOAT64;
  ((double *)(out->array))[1] = n ? sqrt( (s2-s*s/n)/n ) : GAL_BLANK_FLOAT64;
  return out;
}





/* The input is a sorted array with no blank values, we want the median
   value to be put inside the already allocated space which is pointed to
   by `median'. It is in the same type as the input. */
#define MED_IN_SORTED(IT) {                                             \
    IT *a=sorted->array;                                                \
    *(IT *)median = n%2 ? a[n/2]  : (a[n/2]+a[n/2-1])/2;                \
  }
static void
statistics_median_in_sorted_no_blank(gal_data_t *sorted, void *median)
{
  size_t n=sorted->size;

  switch(sorted->type)
    {
    case GAL_DATA_TYPE_UINT8:     MED_IN_SORTED( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      MED_IN_SORTED( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    MED_IN_SORTED( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     MED_IN_SORTED( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    MED_IN_SORTED( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     MED_IN_SORTED( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    MED_IN_SORTED( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     MED_IN_SORTED( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   MED_IN_SORTED( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   MED_IN_SORTED( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_median_in_sorted_no_blank'", sorted->type);
    }
}





/* Return the median value of the dataset in the same type as the input as
   a one element dataset. If the `inplace' flag is set, the input data
   structure will be modified: it will have no blank values and will be
   sorted (increasing). */
gal_data_t *
gal_statistics_median(gal_data_t *input, int inplace)
{
  size_t dsize=1;
  gal_data_t *nbs=gal_statistics_no_blank_sorted(input, inplace);
  gal_data_t *out=gal_data_alloc(NULL, input->type, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);

  /* Write the median. */
  statistics_median_in_sorted_no_blank(nbs, out->array);

  /* Clean up (if necessary), then return the output */
  if(nbs!=input) gal_data_free(nbs);
  return out;
}





/* Return a single element dataset of the same type as input keeping the
   value that has the given quantile. */
#define STATS_QUANT(IT) { IT *o=out->array, *a=nbs->array; *o = a[ind]; }
gal_data_t *
gal_statistsics_quantile(gal_data_t *input, float quantile, int inplace)
{
  size_t dsize=1, ind;
  gal_data_t *nbs=gal_statistics_no_blank_sorted(input, inplace);
  gal_data_t *out=gal_data_alloc(NULL, input->type, 1, &dsize,
                                 NULL, 1, -1, NULL, NULL, NULL);

  /* A small sanity check. */
  if(quantile<0 || quantile>1)
    error(EXIT_FAILURE, 0, "the `quantile' input to "
          "`gal_statistics_quantile' must be between 0 and 1 (inclusive)");

  /* Find the index of the quantile. */
  ind=gal_statistics_quantile_index(nbs->size, quantile);

  /* Set the value. */
  switch(input->type)
    {
    case GAL_DATA_TYPE_UINT8:     STATS_QUANT( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      STATS_QUANT( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    STATS_QUANT( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     STATS_QUANT( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    STATS_QUANT( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     STATS_QUANT( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    STATS_QUANT( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     STATS_QUANT( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   STATS_QUANT( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   STATS_QUANT( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_maximum'", input->type);
    }

  /* Clean up and return. */
  if(nbs!=input) gal_data_free(nbs);
  return out;
}





size_t
gal_statistics_quantile_index(size_t size, float quant)
{
  float floatindex;

  if(quant>1.0f)
    error(EXIT_FAILURE, 0, "the quantile in "
          "gal_statistics_index_from_quantile should be smaller than 1.0");

  /* Find the index of the quantile. */
  floatindex=(float)size*quant;

  /*
  printf("quant: %f, size: %zu, findex: %f\n", quant, size, floatindex);
  */
  /* Note that in the conversion from float to size_t, the floor
     integer value of the float will be used. */
  if( floatindex - (int)floatindex > 0.5 )
    return floatindex+1;
  else
    return floatindex;
}






















/****************************************************************
 ********                      Sort                       *******
 ****************************************************************/
/* Check if the given dataset is sorted. Output values are:

     - 0: Dataset is not sorted.
     - 1: Dataset is sorted and increasing or equal.
     - 2: dataset is sorted and decreasing.                  */

#define IS_SORTED(IT) {                                                 \
    IT *aa=data->array, *a=data->array, *af=a+data->size-1;             \
    if(a[1]>=a[0]) do if( *(a+1) < *a ) break; while(++a<af);           \
    else           do if( *(a+1) > *a ) break; while(++a<af);           \
    return ( a==af                                                      \
             ? ( aa[1]>=aa[0]                                           \
                 ? GAL_STATISTICS_SORTED_INCREASING                     \
                 : GAL_STATISTICS_SORTED_DECREASING )                   \
             : GAL_STATISTICS_SORTED_NOT );                             \
  }

int
gal_statistics_is_sorted(gal_data_t *data)
{
  /* A one-element dataset can be considered, sorted, so we'll just return
     1 (for sorted and increasing). */
  if(data->size==1) return GAL_STATISTICS_SORTED_INCREASING;

  /* Do the check. */
  switch(data->type)
    {
    case GAL_DATA_TYPE_UINT8:     IS_SORTED( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      IS_SORTED( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    IS_SORTED( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     IS_SORTED( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    IS_SORTED( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     IS_SORTED( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    IS_SORTED( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     IS_SORTED( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   IS_SORTED( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   IS_SORTED( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_is_sorted'", data->type);
    }

  /* Control shouldn't reach this point. */
  error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can fix the "
        "problem. For some reason, control has reached the end of "
        "`gal_statistics_is_sorted'", PACKAGE_BUGREPORT);
  return -1;
}





/* This function is ignorant to blank values, if you want to make sure
   there is no blank values, you can call `gal_blank_remove' first. */
#define STATISTICS_SORT(QSORT_F) {                                        \
    qsort(data->array, data->size, gal_data_sizeof(data->type), QSORT_F); \
  }
void
gal_statistics_sort_increasing(gal_data_t *data)
{
  switch(data->type)
    {
    case GAL_DATA_TYPE_UINT8:
      STATISTICS_SORT(gal_qsort_uint8_increasing);    break;
    case GAL_DATA_TYPE_INT8:
      STATISTICS_SORT(gal_qsort_int8_increasing);     break;
    case GAL_DATA_TYPE_UINT16:
      STATISTICS_SORT(gal_qsort_uint16_increasing);   break;
    case GAL_DATA_TYPE_INT16:
      STATISTICS_SORT(gal_qsort_int16_increasing);    break;
    case GAL_DATA_TYPE_UINT32:
      STATISTICS_SORT(gal_qsort_uint32_increasing);   break;
    case GAL_DATA_TYPE_INT32:
      STATISTICS_SORT(gal_qsort_int32_increasing);    break;
    case GAL_DATA_TYPE_UINT64:
      STATISTICS_SORT(gal_qsort_uint64_increasing);   break;
    case GAL_DATA_TYPE_INT64:
      STATISTICS_SORT(gal_qsort_int64_increasing);    break;
    case GAL_DATA_TYPE_FLOAT32:
      STATISTICS_SORT(gal_qsort_float32_increasing);  break;
    case GAL_DATA_TYPE_FLOAT64:
      STATISTICS_SORT(gal_qsort_float64_increasing);  break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_sort_increasing'", data->type);
    }
}





/* See explanations above `gal_statistics_sort_increasing'. */
void
gal_statistics_sort_decreasing(gal_data_t *data)
{
  switch(data->type)
    {
    case GAL_DATA_TYPE_UINT8:
      STATISTICS_SORT(gal_qsort_uint8_decreasing);    break;
    case GAL_DATA_TYPE_INT8:
      STATISTICS_SORT(gal_qsort_int8_decreasing);     break;
    case GAL_DATA_TYPE_UINT16:
      STATISTICS_SORT(gal_qsort_uint16_decreasing);   break;
    case GAL_DATA_TYPE_INT16:
      STATISTICS_SORT(gal_qsort_int16_decreasing);    break;
    case GAL_DATA_TYPE_UINT32:
      STATISTICS_SORT(gal_qsort_uint32_decreasing);   break;
    case GAL_DATA_TYPE_INT32:
      STATISTICS_SORT(gal_qsort_int32_decreasing);    break;
    case GAL_DATA_TYPE_UINT64:
      STATISTICS_SORT(gal_qsort_uint64_decreasing);   break;
    case GAL_DATA_TYPE_INT64:
      STATISTICS_SORT(gal_qsort_int64_decreasing);    break;
    case GAL_DATA_TYPE_FLOAT32:
      STATISTICS_SORT(gal_qsort_float32_decreasing);  break;
    case GAL_DATA_TYPE_FLOAT64:
      STATISTICS_SORT(gal_qsort_float64_decreasing);  break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_sort_decreasing'", data->type);
    }
}





/* Return a dataset that has doesn't have blank values and is sorted. If
   the `inplace' value is set to 1, then the input array will be modified,
   otherwise, a new array will be allocated with the desired properties. So
   if it is already sorted and has blank values, the `inplace' variable is
   irrelevant.*/
gal_data_t *
gal_statistics_no_blank_sorted(gal_data_t *input, int inplace)
{
  gal_data_t *noblank, *sorted;

  /* Make sure there is no blanks in the array that will be used. After
     this step, we won't be dealing with `input' any more, but with
     `noblank'.*/
  if( gal_blank_present(input) )
    {
      if(inplace)                     /* If we can free the input, then */
        {                             /* just remove the blank elements */
          gal_blank_remove(input);    /* in place. */
          noblank=input;
        }
      else
        {
          noblank=gal_data_copy(input);   /* We aren't allowed to touch */
          gal_blank_remove(input);        /* the input, so make a copy. */
        }
    }
  else noblank=input;

  /* Make sure the array is sorted. After this step, we won't be dealing
     with `noblank' any more but with `sorted'. */
  if( gal_statistics_is_sorted(noblank) ) sorted=noblank;
  else
    {
      if(inplace) sorted=noblank;
      else
        {
          if(noblank!=input)    /* no-blank has already been allocated. */
            sorted=noblank;
          else
            sorted=gal_data_copy(noblank);
        }
      gal_statistics_sort_increasing(sorted);
    }

  return sorted;
}




















/****************************************************************
 ********     Histogram and Cumulative Frequency Plot     *******
 ****************************************************************/
/* Generate an array of regularly spaced elements. For a 1D dataset, the
   output will be 1D, for 2D, it will be 2D.

   Input arguments:

     * The `data' set you want to apply the bins to. This is only necessary
       if the range argument is not complete, see below. If `range' has all
       the necessary information, you can pass a NULL pointer for `data'.

     * The `range' data structure keeps the desired range along each
       dimension of the input data structure, it has to be in float32
       type. Note that if

         - If you want the full range of the dataset (in any dimensions,
           then just set `range' to NULL and the range will be specified
           from the minimum and maximum value of the dataset.

         - If there is one element for each dimension in range, then it is
           viewed as a quantile (Q), and the range will be: `Q to 1-Q'.

         - If there are two elements for each dimension in range, then they
           are assumed to be your desired minimum and maximum values. When
           either of the two are NaN, the minimum and maximum will be
           calculated for it.

     * The number of bins: must be larger than 0.

     * `onebinstart' A desired value for onebinstart. Note that with this
        option, the bins won't start and end exactly on the given range
        values, it will be slightly shifted to accommodate this
        request.

  The output is a 1D array (column) of type double, it has to be double to
  account for small differences on the bin edges.
*/
gal_data_t *
gal_statistics_regular_bins(gal_data_t *data, gal_data_t *inrange,
                            size_t numbins, float onebinstart)
{
  size_t i;
  gal_data_t *bins, *tmp, *range;
  double *b, *ra, min, max, hbw, diff, binwidth;


  /* Some sanity checks. */
  if(numbins==0)
    error(EXIT_FAILURE, 0, "`numbins' in `gal_statistics_regular_bins' "
          "cannot be given a value of 0");


  /* Set the minimum and maximum values. */
  if(inrange && inrange->size)
    {
      /* Make sure we are dealing with a double type range. */
      if(inrange->type==GAL_DATA_TYPE_FLOAT64)
        range=inrange;
      else
        range=gal_data_copy_to_new_type(inrange, GAL_DATA_TYPE_FLOAT64);

      /* Set the minimum and maximum of the bins. */
      ra=range->array;
      if( (range->size)%2 )
        error(EXIT_FAILURE, 0, "Quantile ranges are not implemented in "
              "`gal_statistics_regular_bins' yet.");
      else
        {
          /* If the minimum isn't set (is blank), find it. */
          if( isnan(ra[0]) )
            {
              tmp=gal_data_copy_to_new_type_free(gal_statistics_minimum(data),
                                                 GAL_DATA_TYPE_FLOAT64);
              min=*((double *)(tmp->array));
              gal_data_free(tmp);
            }
          else min=ra[0];

          /* For the maximum, when it isn't set, we'll add a very small
             value, so all points are included. */
          if( isnan(ra[1]) )
            {
              tmp=gal_data_copy_to_new_type_free(gal_statistics_maximum(data),
                                                 GAL_DATA_TYPE_FLOAT64);
              max=*((double *)(tmp->array))+1e-6;
              gal_data_free(tmp);
            }
          else max=ra[1];
        }

      /* Clean up: if `range' was allocated. */
      if(range!=inrange) gal_data_free(range);
    }
  /* No range was given, find the minimum and maximum. */
  else
    {
      tmp=gal_data_copy_to_new_type_free(gal_statistics_minimum(data),
                                         GAL_DATA_TYPE_FLOAT64);
      min=*((double *)(tmp->array));
      gal_data_free(tmp);
      tmp=gal_data_copy_to_new_type_free(gal_statistics_maximum(data),
                                         GAL_DATA_TYPE_FLOAT64);
      max=*((double *)(tmp->array)) + 1e-6;
      gal_data_free(tmp);
    }


  /* Allocate the space for the bins. */
  bins=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT64, 1, &numbins, NULL,
                      0, data->minmapsize, "bin_center", data->unit,
                      "Center value of each bin.");


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
          diff=onebinstart-b[i];
          for(i=0;i<numbins;++i)
            b[i]+=diff;
        }
    }

  /* For a check
  printf("min: %g\n", min);
  printf("max: %g\n", max);
  printf("binwidth: %g\n", binwidth);
  for(i=0;i<numbins;++i)
    printf("%zu: %.4f\n", i, b[i]);
  */

  /* Set the status of the bins to regular and return. */
  bins->status=GAL_STATISTICS_BINS_REGULAR;
  return bins;
}





/* Make a histogram of all the elements in the given dataset with bin
   values that are defined in the `inbins' structure (see
   `gal_statistics_regular_bins'). `inbins' is not mandatory, if you pass a
   NULL pointer, the bins structure will be built within this function
   based on the `numbins' input. As a result, when you have already defined
   the bins, `numbins' is not used. */

#define HISTOGRAM_TYPESET(IT) {                                         \
    IT *a=data->array, *af=a+data->size;                                \
    do if( *a>=min && *a<max) ++h[ (size_t)( (*a-min)/binwidth ) ];     \
    while(++a<af);                                                      \
  }

gal_data_t *
gal_statistics_histogram(gal_data_t *data, gal_data_t *bins, int normalize,
                         int maxone)
{
  size_t *h;
  float *f, *ff;
  gal_data_t *hist;
  double *d, min, max, ref=NAN, binwidth;


  /* Check if the bins are regular or not. For irregular bins, we can
     either use the old implementation, or GSL's histogram
     functionality. */
  if(bins==NULL)
    error(EXIT_FAILURE, 0, "no `bins' in `gal_statistics_histogram");
  if(bins->status!=GAL_STATISTICS_BINS_REGULAR)
    error(EXIT_FAILURE, 0, "the input bins to `gal_statistics_histogram' "
          "are not regular. Currently it is only implemented for regular "
          "bins");


  /* Check if normalize and `maxone' are not called together. */
  if(normalize && maxone)
    error(EXIT_FAILURE, 0, "only one of `normalize' and `maxone' may "
          "be given to `gal_statistics_histogram'");


  /* Allocate the histogram (note that we are clearning it so all values
     are zero. */
  hist=gal_data_alloc(NULL, GAL_DATA_TYPE_SIZE_T, bins->ndim, bins->dsize,
                      NULL, 1, data->minmapsize, "hist_number", "counts",
                      "Number of data points within each bin.");


  /* Set the minimum and maximum range of the histogram from the bins. */
  d=bins->array;
  binwidth=d[1]-d[0];
  max = d[bins->size - 1] + binwidth/2;
  min = d[0]              - binwidth/2;


  /* Go through all the elements and find out which bin they belong to. */
  h=hist->array;
  switch(data->type)
    {
    case GAL_DATA_TYPE_UINT8:     HISTOGRAM_TYPESET(uint8_t);     break;
    case GAL_DATA_TYPE_INT8:      HISTOGRAM_TYPESET(int8_t);      break;
    case GAL_DATA_TYPE_UINT16:    HISTOGRAM_TYPESET(uint16_t);    break;
    case GAL_DATA_TYPE_INT16:     HISTOGRAM_TYPESET(int16_t);     break;
    case GAL_DATA_TYPE_UINT32:    HISTOGRAM_TYPESET(uint32_t);    break;
    case GAL_DATA_TYPE_INT32:     HISTOGRAM_TYPESET(int32_t);     break;
    case GAL_DATA_TYPE_UINT64:    HISTOGRAM_TYPESET(uint64_t);    break;
    case GAL_DATA_TYPE_INT64:     HISTOGRAM_TYPESET(int64_t);     break;
    case GAL_DATA_TYPE_FLOAT32:   HISTOGRAM_TYPESET(float);       break;
    case GAL_DATA_TYPE_FLOAT64:   HISTOGRAM_TYPESET(double);      break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_histogram'", data->type);
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
      hist=gal_data_copy_to_new_type_free(hist, GAL_DATA_TYPE_FLOAT32);
      ff=(f=hist->array)+hist->size; do ref += *f++;   while(f<ff);

      /* Correct the name, units and comments. */
      free(hist->name); free(hist->unit); free(hist->comment);
      gal_checkset_allocate_copy("hist_normalized", &hist->name);
      gal_checkset_allocate_copy("frac", &hist->unit);
      gal_checkset_allocate_copy("Normalized histogram value for this bin",
                                 &hist->comment);
    }
  if(maxone)
    {
      /* Calculate the reference. */
      ref=-FLT_MAX;
      hist=gal_data_copy_to_new_type_free(hist, GAL_DATA_TYPE_FLOAT32);
      ff=(f=hist->array)+hist->size;
      do ref = *f>ref ? *f : ref; while(++f<ff);

      /* Correct the name, units and comments. */
      free(hist->name); free(hist->unit); free(hist->comment);
      gal_checkset_allocate_copy("hist_maxone", &hist->name);
      gal_checkset_allocate_copy("frac", &hist->unit);
      gal_checkset_allocate_copy("Fractional histogram value for this bin "
                                 "when maximum bin value is 1.0",
                                 &hist->comment);
    }


  /* Correct the histogram if necessary. */
  if( !isnan(ref) )
    { ff=(f=hist->array)+hist->size; do *f++ /= ref;   while(f<ff); }


  /* Return the histogram. */
  return hist;
}





/* Make a cumulative frequency plot (CFP) of all the elements in the given
   dataset with bin values that are defined in the `bins' structure (see
   `gal_statistics_regular_bins').

   The CFP is built from the histogram: in each bin, the value is the sum
   of all previous bins in the histogram. Thus, if you have already
   calculated the histogram before calling this function, you can pass it
   onto this function as the data structure in `bins->next'. If
   `bin->next!=NULL', then it is assumed to be the histogram. If it is
   NULL, then the histogram will be calculated internally and freed after
   the job is finished.

   When a histogram is given and it is normalized, the CFP will also be
   normalized (even if the normalized flag is not set here): note that a
   normalized CFP's maximum value is 1.
*/
gal_data_t *
gal_statistics_cfp(gal_data_t *data, gal_data_t *bins, int normalize)
{
  double sum;
  float *f, *ff, *hf;
  gal_data_t *hist, *cfp;
  size_t *s, *sf, *hs, sums;


  /* Check if the bins are regular or not. For irregular bins, we can
     either use the old implementation, or GSL's histogram
     functionality. */
  if(bins->status!=GAL_STATISTICS_BINS_REGULAR)
    error(EXIT_FAILURE, 0, "the input bins to `gal_statistics_cfp' "
          "are not regular. Currently it is only implemented for regular "
          "bins");


  /* Prepare the histogram. */
  hist = ( bins->next
           ? bins->next
           : gal_statistics_histogram(data, bins, 0, 0) );


  /* If the histogram has float32 type it was given by the user and is
     either normalized or its maximum was set to 1. We can only use it if
     it was normalized. If it isn't normalized, then we must ignore it and
     build the histogram here.*/
  if(hist->type==GAL_DATA_TYPE_FLOAT32)
    {
      sum=0.0f;
      ff=(f=hist->array)+hist->size; do sum += *f++;   while(f<ff);
      if(sum!=1.0f)
        hist=gal_statistics_histogram(data, bins, 0, 0);
    }


  /* Allocate the cumulative frequency plot's necessary space. */
  cfp=gal_data_alloc( NULL, hist->type, bins->ndim, bins->dsize,
                      NULL, 1, data->minmapsize,
                      ( hist->type==GAL_DATA_TYPE_FLOAT32
                        ? "cfp_normalized" : "cfp_number" ),
                      ( hist->type==GAL_DATA_TYPE_FLOAT32
                        ? "frac" : "count" ),
                      ( hist->type==GAL_DATA_TYPE_FLOAT32
                        ? "Fraction of data elements from the start to this "
                        "bin (inclusive)."
                        : "Number of data elements from the start to this "
                        "bin (inclusive).") );


  /* Fill in the cumulative frequency plot. */
  switch(hist->type)
    {
    case GAL_DATA_TYPE_SIZE_T:
      sums=0; hs=hist->array; sf=(s=cfp->array)+cfp->size;
      do sums = (*s += *hs++ + sums); while(++s<sf);
      break;

    case GAL_DATA_TYPE_FLOAT32:
      sum=0.0f; hf=hist->array; ff=(f=cfp->array)+cfp->size;
      do sum = (*f += *hf++ + sum);  while(++f<ff);
      break;

    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_cfp'", cfp->type);
    }


  /* Normalize the CFP if the user asked for it and it wasn't normalized
     until now. */
  if(normalize && cfp->type==GAL_DATA_TYPE_SIZE_T)
    {
      /* Find the sum, then divide the plot by it. Note that the sum must
         come from the histogram, not the CFP!*/
      sums=0;
      cfp=gal_data_copy_to_new_type_free(cfp, GAL_DATA_TYPE_FLOAT32);
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
 *****************       Sigma clip          ********************
 ****************************************************************/
/* Return a data structure with an array of four values: the final number
   of points used, the median, average and standard deviation. The number
   of clips is put into the `status' element of the data structure.

   Inputs:

     - `multip': multiple of the standard deviation,

     - `param' must be positive and determines the type of clipping:

         - param<1.0: interpretted as a tolerance level to stop clipping.

         - param>=1.0 and an integer: a specific number of times to do the
           clippping.

  The way this function works is very simple: first it will sort the input
  (if it isn't sorted). Afterwards, it will recursively change the starting
  point of the array and its size, calcluating the basic statistics in each
  round to define the new starting point and size.
*/

#define SIGCLIP(IT) {                                                   \
    IT *a  = sorted->array, *af = a  + sorted->size;                    \
    IT *bf = sorted->array, *b  = bf + sorted->size - 1;                \
                                                                        \
    /* Remove all out-of-range elements from the start of the array. */ \
    if(sortstatus==GAL_STATISTICS_SORTED_INCREASING)                    \
      do if( *a > (*med - (multip * *std)) )                            \
           { start=a; break; }                                          \
      while(++a<af);                                                    \
    else                                                                \
      do if( *a < (*med + (multip * *std)) )                            \
           { start=a; break; }                                          \
      while(++a<af);                                                    \
                                                                        \
    /* Remove all out-of-range elements from the end of the array. */   \
    if(sortstatus==GAL_STATISTICS_SORTED_INCREASING)                    \
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
                          int quiet)
{
  int sortstatus;
  double *med, *mean, *std;
  void *start, *sorted_array;
  double oldmed, oldmean, oldstd;
  size_t num=0, dsize=4, size, oldsize;
  size_t maxnum = param>=1.0f ? param : 50;    /* for failing to converge */
  uint8_t bytolerance = param>=1.0f ? 0 : 1;
  gal_data_t *sorted, *median_i, *median_d, *out, *meanstd, *noblank;

  /* Some sanity checks. */
  if( multip<=0 )
    error(EXIT_FAILURE, 0, "`multip', must be greater than zero in "
          "`gal_statistics_sigma_clip'. The given value was %g", multip);
  if( param<=0 )
    error(EXIT_FAILURE, 0, "`param', must be greater than zero in "
          "`gal_statistics_sigma_clip'. The given value was %g", param);
  if( param >= 1.0f && ceil(param) != param )
    error(EXIT_FAILURE, 0, "when `param' is larger than 1.0, it is "
          "interpretted as an absolute number of clips in "
          "`gal_statistics_sigma_clip'. So it must be an integer. However, "
          "your given value %g", param);


  /* If there are blank elements, remove them (from a copied array). NOTE:
     we don't want to change the input. */
  if( gal_blank_present(input) )
    {
      noblank = gal_data_copy(input);
      gal_blank_remove(noblank);
    }
  else noblank=input;


  /* We want to have sorted (increasing) array. */
  sortstatus=gal_statistics_is_sorted(noblank);
  switch(sortstatus)
    {
    case GAL_STATISTICS_SORTED_NOT:
      /* Only copy (or allocate) space if it wasn't already allocated. */
      sorted = (noblank==input) ? gal_data_copy(input) : noblank;
      gal_statistics_sort_increasing(sorted);
      sortstatus=GAL_STATISTICS_SORTED_INCREASING;
      break;

    case GAL_STATISTICS_SORTED_DECREASING:
    case GAL_STATISTICS_SORTED_INCREASING:
      sorted=noblank;
      break;

    default:
      error(EXIT_FAILURE, 0, "sorted code %d not recognized in "
            "`gal_statistics_sigma_clip'", gal_statistics_is_sorted(input));
    }


  /* Allocate the necessary spaces. */
  out=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT32, 1, &dsize, NULL, 0,
                     input->minmapsize, NULL, NULL, NULL);
  median_i=gal_data_alloc(NULL, sorted->type, 1, &dsize, NULL, 0,
                          input->minmapsize, NULL, NULL, NULL);


  /* Print the comments. */
  if(!quiet)
    printf("%-8s %-10s %-15s %-15s %-15s\n",
           "round", "number", "median", "mean", "STD");


  /* Do the clipping, but first initialize the values that will be changed
     during the clipping: the start of the array and the array's size. */
  size=sorted->size;
  sorted_array=start=sorted->array;
  while(num<maxnum)
    {
      /* Find the median. */
      statistics_median_in_sorted_no_blank(sorted, median_i->array);
      median_d=gal_data_copy_to_new_type(median_i, GAL_DATA_TYPE_FLOAT64);

      /* Find the average and Standard deviation, note that both `start'
         and `size' will be different in the next round. */
      sorted->array = start;
      sorted->size = oldsize = size;
      meanstd=gal_statistics_mean_std(sorted);

      /* Put the three final values in usable (with a type) pointers. */
      med = median_d->array;
      mean = meanstd->array;
      std  = &((double *)(meanstd->array))[1];

      /* If the user wanted to view the steps, show it to them. */
      if(!quiet)
        printf("%-8zu %-10zu %-15g %-15g %-15g\n",
               num+1, size, *med, *mean, *std);

      /* If we are to work by tolerance, then check if we should jump out
         of the loop. Normally, `oldstd' should be larger than std, because
         the possible outliers have been removed. If it is not, it means
         that we have clipped too much and must stop anyway, so we don't
         need an absolute value on the difference! */
      if( bytolerance && num>0 && ((oldstd - *std) / *std) < param )
        break;

      /* Clip all the elements outside of the desired range: since the
         array is sorted, this means to just change the starting pointer
         and size of the array. */
      switch(sorted->type)
        {
        case GAL_DATA_TYPE_UINT8:     SIGCLIP( uint8_t  );   break;
        case GAL_DATA_TYPE_INT8:      SIGCLIP( int8_t   );   break;
        case GAL_DATA_TYPE_UINT16:    SIGCLIP( uint16_t );   break;
        case GAL_DATA_TYPE_INT16:     SIGCLIP( int16_t  );   break;
        case GAL_DATA_TYPE_UINT32:    SIGCLIP( uint32_t );   break;
        case GAL_DATA_TYPE_INT32:     SIGCLIP( int32_t  );   break;
        case GAL_DATA_TYPE_UINT64:    SIGCLIP( uint64_t );   break;
        case GAL_DATA_TYPE_INT64:     SIGCLIP( int64_t  );   break;
        case GAL_DATA_TYPE_FLOAT32:   SIGCLIP( float    );   break;
        case GAL_DATA_TYPE_FLOAT64:   SIGCLIP( double   );   break;
        default:
          error(EXIT_FAILURE, 0, "type code %d not recognized in "
                "`gal_statistics_sigma_clip'", sorted->type);
        }

      /* Set the values from this round in the old elements, so the next
         round can compare with, and return then if necessary. */
      oldmed =  *med;
      oldstd  = *std;
      oldmean = *mean;
      ++num;

      /* Clean up: */
      gal_data_free(meanstd);
      gal_data_free(median_d);
    }

  /* If we were in tolerance mode and `num' and `maxnum' are equal (the
     loop didn't stop by tolerance), so the outputs should be NaN. */
  out->status=num;
  if( bytolerance && num==maxnum )
    {
      ((float *)(out->array))[0] = NAN;
      ((float *)(out->array))[1] = NAN;
      ((float *)(out->array))[2] = NAN;
      ((float *)(out->array))[3] = NAN;
    }
  else
    {
      ((float *)(out->array))[0] = size;
      ((float *)(out->array))[1] = oldmed;
      ((float *)(out->array))[2] = oldmean;
      ((float *)(out->array))[3] = oldstd;
    }

  /* Clean up and return. */
  sorted->array=sorted_array;
  if(sorted!=input) gal_data_free(sorted);
  return out;
}




















/****************************************************************/
/*************         Identify outliers         ****************/
/****************************************************************/
/* Using the cumulative distribution function this funciton will
   remove outliers from a dataset. */
void
gal_statistics_remove_outliers_flat_cdf(float *sorted, size_t *outsize)
{
  printf("\n ... in gal_statistics_remove_outliers_flat_cdf ... \n");
  exit(1);
#if 0
  int firstfound=0;
  size_t size=*outsize, i, maxind;
  float *slopes, minslope, maxslope;

  /* Find a slopes array, think of the cumulative frequency plot when
     you want to think about slopes. */
  errno=0; slopes=malloc(size*sizeof *slopes);
  if(slopes==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for slopes in "
          "gal_statistics_remove_outliers_flat_cdf (statistics.c)",
          size*sizeof *slopes);

  /* Calcuate the slope of the CDF and put it in the slopes array. */
  for(i=1;i<size-1;++i)
    slopes[i]=2/(sorted[i+1]-sorted[i-1]);

  /* Find the position of the maximum slope, note that around the
     distribution mode, the difference between the values varies less,
     so two neighbouring elements have the closest values, hence the
     largest slope (when their difference is in the denominator). */
  gal_statistics_f_max_with_index(slopes+1, size-2, &maxslope, &maxind);

  /* Find the minimum slope from the second element (for the first the
     slope is not defined. NOTE; maxind is one smaller than it should
     be because the input array to find it began from the second
     element. */
  gal_statistics_float_second_min(slopes+1, maxind+1, &minslope);

  /* Find the second place where the slope falls below `minslope`
     after the maximum position. When found, add it with one to
     account for error. Note that incase there are no outliers by this
     definition, then the final size will be equal to size! */
  for(i=maxind+1;i<size-1;++i)
    if(slopes[i]<minslope)
      {
        if(firstfound)
          break;
        else
          firstfound=1;
      }
  *outsize=i+1;

  /*
  for(i=0;i<size;++i)
    printf("%zu\t%.3f\t%.3f\n", i, arr[i], slopes[i]);
  printf("\n\nPlace to cut off for outliers is: %zu\n\n", *outsize);
  */

  free(slopes);
#endif
}




















/****************************************************************/
/*************          Mode calculation         ****************/
/****************************************************************/
void
gal_statistics_mode_mirror_plots(float *sorted, size_t size, size_t mirrorindex,
                                 float min, float max, size_t numbins,
                                 char *histsname, char *cfpsname,
                                 float mirrorplotdist)
{
  printf("\n... in gal_statistics_mode_mirror_plots ...\n");
  exit(1);
#if 0
  FILE *fp;
  size_t i, msize;
  float *out, maxhist=-FLT_MAX, maxcfp, d;
  int normhist=0, maxhistone=0, normcfp=0;
  float *bins, *mirror, *actual, mf, onebinvalue=0.0f;


  /* Find the index of the mirror and make the mirror array: */
  mf=sorted[mirrorindex];
  gal_array_float_copy(sorted, size, &actual);
  makemirrored(sorted, mirrorindex, &mirror, &msize);


  /* Set the best range if asked for, such that the mirror is on the
     1/3 of the image scale. */
  if(mirrorplotdist!=0.0f)
    {
      min=actual[gal_statistics_index_from_quantile(size, 0.001f)];
      max=mf+mirrorplotdist*(mf-min);
    }


  /* set the mirror pixel to have the value of zero.*/
  min-=mf;
  max-=mf;
  gal_array_fsum_const(actual, size, -1.0f*mf);
  gal_array_fsum_const(mirror, msize, -1.0f*mf);


  /* Allocate space for the 3 column array keeping the histograms:*/
  errno=0;
  out=malloc(numbins*3*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for out in "
          "gal_mode_make_mirror_plots (mode.c)", numbins*3*sizeof *out);


  /* Define the bin sides: */
  gal_statistics_set_bins(actual, size, numbins, min,
                          max, onebinvalue, 0, &bins);


  /* Find the histogram of the actual data and put it in out. Note
     that maxhistone=0, because here we want to use one value for both
     histograms so they are comparable. */
  gal_statistics_histogram(actual, size, bins, numbins,
                           normhist, maxhistone);
  for(i=0;i<numbins;++i)
    if(bins[i*2+1]>maxhist) maxhist=bins[i*2+1];
  for(i=0;i<numbins;++i)
    { out[i*3]=bins[i*2]; out[i*3+1]=bins[i*2+1]/maxhist; bins[i*2+1]=0.0f;}
  bins[i*2+1]=0.0f; /* bins[] actually has numbins+1 elements. */
  d=(out[3]-out[0])/2;


  /* Find the histogram of the mirrired distribution and put it in
     out: */
  gal_statistics_histogram(mirror, msize, bins, numbins, normhist,
                           maxhistone);
  for(i=0;i<numbins;++i)
    { out[i*3+2]=bins[i*2+1]/maxhist; bins[i*2+1]=0.0f;}
  bins[i*2+1]=0.0f; /* bins[] actually has numbins+1 elements. */


  /* Print out the histogram: */
  errno=0;
  fp=fopen(histsname, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "could not open file %s", histsname);
  fprintf(fp, "# Histogram of actual and mirrored distributions.\n");
  fprintf(fp, "# Column 0: Value in the middle of this bin.\n");
  fprintf(fp, "# Column 1: Input data.\n");
  fprintf(fp, "# Column 2: Mirror distribution.\n");
  for(i=0;i<numbins;++i)
    fprintf(fp, "%-25.6f%-25.6f%-25.6f\n", out[i*3]+d,
            out[i*3+1], out[i*3+2]);
  fclose(fp);




  /* Find the cumulative frequency plots of the two distributions: */
  gal_statistics_cumulative_fp(actual, size, bins, numbins, normcfp);
  for(i=0;i<numbins;++i)
    { out[i*3+1]=bins[i*2+1]; bins[i*2+1]=0.0f; }
  bins[i*2+1]=0.0f; /* bins[] actually has numbins+1 elements. */
  gal_statistics_cumulative_fp(mirror, msize, bins, numbins, normcfp);
  for(i=0;i<numbins;++i)
    { out[i*3+2]=bins[i*2+1]; bins[i*2+1]=0.0f;}
  bins[i*2+1]=0.0f; /* bins[] actually has numbins+1 elements. */


  /* Since the two cumultiave frequency plots have to be on scale, see
     which one is larger and divided both CFPs by the size of the
     larger one. Then print the CFPs. */
  if(size>msize) maxcfp=size;
  else maxcfp=msize;
  errno=0;
  fp=fopen(cfpsname, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "could not open file %s", cfpsname);
  fprintf(fp, "# Cumulative frequency plot (average index in bin) of\n"
          "# Actual and mirrored distributions.\n");
  fprintf(fp, "# Column 0: Value in the middle of this bin.\n");
  fprintf(fp, "# Column 1: Actual data.\n");
  fprintf(fp, "# Column 2: Mirror distribution.\n");
  for(i=0;i<numbins;++i)
    fprintf(fp, "%-25.6f%-25.6f%-25.6f\n", out[i*3],
            out[i*3+1]/maxcfp, out[i*3+2]/maxcfp);
  fclose(fp);


  /* Clean up. */
  free(out);
  free(bins);
  free(mirror);
  free(actual);
#endif
}





/* It happens that you have the symmetricity and you want the flux
   value at that point, this function will do that job. In practice,
   it just finds bf from the equation to calculate symmetricity in
   modesymmetricity. */
float
gal_statistics_mode_value_from_sym(float *sorted, size_t size,
                                   size_t modeindex, float sym)
{
  float mf=sorted[modeindex];
  float af=
    sorted[gal_statistics_quantile_index(2*modeindex+1,
                           GAL_STATISTICS_MODE_SYMMETRICITY_LOW_QUANT)];
  return sym*(mf-af)+mf;
}





/* Find the quantile of the mode of a sorted distribution. The return
   value is either 0 (not accurate) or 1 (accurate). Accuracy is
   defined based on the difference between the maximum and minimum
   maxdiffs that were found during the golden section search.

   A good mode will have:

   modequant=(float)(modeindex)/(float)size;
   modesym>GAL_MODE_SYM_GOOD && modequant>GAL_MODE_LOW_QUANT_GOOD
*/
void
gal_statistics_mode_index_in_sorted(float *sorted, size_t size, float errorstdm,
                                    size_t *modeindex, float *modesym)
{
  struct gal_statistics_mode_params mp;

  /* Initialize the gal_mode_params structure: */
  mp.size=size;
  mp.sorted=sorted;
  mp.tolerance=0.01;
  mp.numcheck=size/2;
  mp.errorstdm=errorstdm;
  if(mp.numcheck>1000)
    mp.interval=mp.numcheck/1000;
  else mp.interval=1;
  mp.lowi  = gal_statistics_quantile_index(size,
                                 GAL_STATISTICS_MODE_LOW_QUANTILE);
  mp.highi = gal_statistics_quantile_index(size,
                                 GAL_STATISTICS_MODE_HIGH_QUANTILE);
  mp.midi  = (((float)mp.highi
               + MODE_GOLDEN_RATIO*(float)mp.lowi)
              /(1+MODE_GOLDEN_RATIO));
  mp.midd  = mirrormaxdiff(mp.sorted, mp.size, mp.midi, mp.numcheck,
                           mp.interval, mp.errorstdm);

  /* Do the golden section search and find the resulting
     symmetricity. */
  *modeindex=modegoldenselection(&mp);
  modesymmetricity(sorted, size, *modeindex, errorstdm, modesym);
}
