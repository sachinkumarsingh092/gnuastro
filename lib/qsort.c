/*********************************************************************
qsort -- Functions used by qsort to sort an array.
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
#include <stdlib.h>
#include <stdint.h>

#include <fitsio.h>

#include <gnuastro/qsort.h>


/*****************************************************************/
/**********                  Macros               ****************/
/*****************************************************************/
/* When one or both elements are NaN, the simple comparison, like '(tb >
   ta) - (tb < ta)', will give 0 (as if the elements are equal). However,
   some preference has to be given to the NaN element in a comparison,
   otherwise the output is not going to be reasonable. We also don't want
   to check NaNs on every comparison (it will slow down the processing).

   So we'll exploit the fact that when there comparison result doesn't
   equal zero, we don't have any NaNs and this 'COMPARE_FLOAT_POSTPROCESS'
   macro is called only when the comparison gives zero. Being larger or
   smaller isn't defined for NaNs, so we'll just put them in the end of the
   sorted list whether it is sorted by decreasing or increasing mode.*/
#define COMPARE_FLOAT_POSTPROCESS (                                     \
   isnan(ta) && isnan(tb)                                               \
   ? 0                                /* Both NaN, define as equal. */  \
   /* One is NaN, one isn't. */                                         \
   : ( isnan(ta)                                                        \
       ? 1                         /* First is NaN, set as smaller. */  \
       : ( isnan(tb)                                                    \
           ? -1                    /* Second is NaN, set as larger. */  \
           : 0 )                      /* None are NaN, set as equal.*/  \
       )                                                                \
)



/*****************************************************************/
/**********        Sorting of actual array        ****************/
/*****************************************************************/
int
gal_qsort_uint8_d(const void *a, const void *b)
{
  return ( *(uint8_t *)b - *(uint8_t *)a );
}

int
gal_qsort_uint8_i(const void *a, const void *b)
{
  return ( *(uint8_t *)a - *(uint8_t *)b );
}

int
gal_qsort_int8_d(const void *a, const void *b)
{
  return ( *(int8_t *)b - *(int8_t *)a );
}

int
gal_qsort_int8_i(const void *a, const void *b)
{
  return ( *(int8_t *)a - *(int8_t *)b );
}

int
gal_qsort_uint16_d(const void *a, const void *b)
{
  return ( *(uint16_t *)b - *(uint16_t *)a );
}

int
gal_qsort_uint16_i(const void *a, const void *b)
{
  return ( *(uint16_t *)a - *(uint16_t *)b );
}

int
gal_qsort_int16_d(const void *a, const void *b)
{
  return ( *(int16_t *)b - *(int16_t *)a );
}

int
gal_qsort_int16_i(const void *a, const void *b)
{
  return ( *(int16_t *)a - *(int16_t *)b );
}

int
gal_qsort_uint32_d(const void *a, const void *b)
{
  return ( *(uint32_t *)b - *(uint32_t *)a );
}

int
gal_qsort_uint32_i(const void *a, const void *b)
{
  return ( *(uint32_t *)a - *(uint32_t *)b );
}

int
gal_qsort_int32_d(const void *a, const void *b)
{
  return ( *(int32_t *)b - *(int32_t *)a );
}

int
gal_qsort_int32_i(const void *a, const void *b)
{
  return ( *(int32_t *)a - *(int32_t *)b );
}

int
gal_qsort_uint64_d(const void *a, const void *b)
{
  return ( *(uint64_t *)b - *(uint64_t *)a );
}

int
gal_qsort_uint64_i(const void *a, const void *b)
{
  return ( *(uint64_t *)a - *(uint64_t *)b );
}

int
gal_qsort_int64_d(const void *a, const void *b)
{
  return ( *(int64_t *)b - *(int64_t *)a );
}

int
gal_qsort_int64_i(const void *a, const void *b)
{
  return ( *(int64_t *)a - *(int64_t *)b );
}

int
gal_qsort_float32_d(const void *a, const void *b)
{
  float ta=*(float*)a;
  float tb=*(float*)b;
  int out=(tb > ta) - (tb < ta);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}

int
gal_qsort_float32_i(const void *a, const void *b)
{
  float ta=*(float*)a;
  float tb=*(float*)b;
  int out=(ta > tb) - (ta < tb);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}

int
gal_qsort_float64_d(const void *a, const void *b)
{
  double ta=*(double*)a;
  double tb=*(double*)b;
  int out=(tb > ta) - (tb < ta);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}

int
gal_qsort_float64_i(const void *a, const void *b)
{
  double ta=*(double*)a;
  double tb=*(double*)b;
  int out=(ta > tb) - (ta < tb);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}




















/*****************************************************************/
/***************          Sorting indexs        ******************/
/*****************************************************************/
/* Initialize the array for sorting indexs to NULL. */
void *gal_qsort_index_single=NULL;

int
gal_qsort_index_single_uint8_d(const void *a, const void *b)
{
  uint8_t ta=((uint8_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  uint8_t tb=((uint8_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_single_uint8_i(const void *a, const void *b)
{
  uint8_t ta=((uint8_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  uint8_t tb=((uint8_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_index_single_int8_d(const void *a, const void *b)
{
  int8_t ta=((int8_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  int8_t tb=((int8_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_single_int8_i(const void *a, const void *b)
{
  int8_t ta=((int8_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  int8_t tb=((int8_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_index_single_uint16_d(const void *a, const void *b)
{
  uint16_t ta=((uint16_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  uint16_t tb=((uint16_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_single_uint16_i(const void *a, const void *b)
{
  uint16_t ta=((uint16_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  uint16_t tb=((uint16_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_index_single_int16_d(const void *a, const void *b)
{
  int16_t ta=((int16_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  int16_t tb=((int16_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_single_int16_i(const void *a, const void *b)
{
  int16_t ta=((int16_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  int16_t tb=((int16_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_index_single_uint32_d(const void *a, const void *b)
{
  uint32_t ta=((uint32_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  uint32_t tb=((uint32_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_single_uint32_i(const void *a, const void *b)
{
  uint32_t ta=((uint32_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  uint32_t tb=((uint32_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_index_single_int32_d(const void *a, const void *b)
{
  int32_t ta=((int32_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  int32_t tb=((int32_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_single_int32_i(const void *a, const void *b)
{
  int32_t ta=((int32_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  int32_t tb=((int32_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_index_single_uint64_d(const void *a, const void *b)
{
  uint64_t ta=((uint64_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  uint64_t tb=((uint64_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_single_uint64_i(const void *a, const void *b)
{
  uint64_t ta=((uint64_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  uint64_t tb=((uint64_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_index_single_int64_d(const void *a, const void *b)
{
  int64_t ta=((int64_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  int64_t tb=((int64_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_single_int64_i(const void *a, const void *b)
{
  int64_t ta=((int64_t *)(gal_qsort_index_single))[ *(size_t *)a ];
  int64_t tb=((int64_t *)(gal_qsort_index_single))[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_index_single_float32_d(const void *a, const void *b)
{
  float ta=((float *)(gal_qsort_index_single))[ *(size_t *)a ];
  float tb=((float *)(gal_qsort_index_single))[ *(size_t *)b ];
  int out=(tb > ta) - (tb < ta);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}

int
gal_qsort_index_single_float32_i(const void *a, const void *b)
{
  float ta=((float *)(gal_qsort_index_single))[ *(size_t *)a ];
  float tb=((float *)(gal_qsort_index_single))[ *(size_t *)b ];
  int out=(ta > tb) - (ta < tb);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}

int
gal_qsort_index_single_float64_d(const void *a, const void *b)
{
  double ta=((double *)(gal_qsort_index_single))[ *(size_t *)a ];
  double tb=((double *)(gal_qsort_index_single))[ *(size_t *)b ];
  int out=(tb > ta) - (tb < ta);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}

int
gal_qsort_index_single_float64_i(const void *a, const void *b)
{
  double ta=((double *)(gal_qsort_index_single))[ *(size_t *)a ];
  double tb=((double *)(gal_qsort_index_single))[ *(size_t *)b ];
  int out=(ta > tb) - (ta < tb);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}

int
gal_qsort_index_multi_d(const void *a, const void *b)
{
  /* Define the structures. */
  struct gal_qsort_index_multi *A = (struct gal_qsort_index_multi *)a;
  struct gal_qsort_index_multi *B = (struct gal_qsort_index_multi *)b;

  /* For easy reading. */
  float ta=A->values[ A->index ];
  float tb=B->values[ B->index ];

  /* Return the result. */
  int out=(tb > ta) - (tb < ta);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}

int
gal_qsort_index_multi_i(const void *a, const void *b)
{
  /* Define the structures. */
  struct gal_qsort_index_multi *A = (struct gal_qsort_index_multi *)a;
  struct gal_qsort_index_multi *B = (struct gal_qsort_index_multi *)b;

  /* For easy reading. */
  float ta=A->values[ A->index ];
  float tb=B->values[ B->index ];

  /* Return the result. */
  int out=(ta > tb) - (ta < tb);
  return out ? out : COMPARE_FLOAT_POSTPROCESS;
}
