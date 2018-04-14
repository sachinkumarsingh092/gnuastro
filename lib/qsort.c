/*********************************************************************
forqsort -- Functions used by qsort to sort an array.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2018, Free Software Foundation, Inc.

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

#include <stdlib.h>
#include <stdint.h>

#include <fitsio.h>

#include <gnuastro/qsort.h>

/* Initialize the array for sorting indexs to NULL. */
float *gal_qsort_index_arr;

int
gal_qsort_index_float_decreasing(const void * a, const void * b)
{
  float ta=gal_qsort_index_arr[ *(size_t *)a ];
  float tb=gal_qsort_index_arr[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_index_float_increasing(const void * a, const void * b)
{
  float ta=gal_qsort_index_arr[ *(size_t *)a ];
  float tb=gal_qsort_index_arr[ *(size_t *)b ];
  return (ta > tb) - (ta < tb);
}










int
gal_qsort_uint8_decreasing(const void *a, const void *b)
{
  return ( *(uint8_t *)b - *(uint8_t *)a );
}

int
gal_qsort_uint8_increasing(const void *a, const void *b)
{
  return ( *(uint8_t *)a - *(uint8_t *)b );
}

int
gal_qsort_int8_decreasing(const void *a, const void *b)
{
  return ( *(int8_t *)b - *(int8_t *)a );
}

int
gal_qsort_int8_increasing(const void *a, const void *b)
{
  return ( *(int8_t *)a - *(int8_t *)b );
}

int
gal_qsort_uint16_decreasing(const void *a, const void *b)
{
  return ( *(uint16_t *)b - *(uint16_t *)a );
}

int
gal_qsort_uint16_increasing(const void *a, const void *b)
{
  return ( *(uint16_t *)a - *(uint16_t *)b );
}

int
gal_qsort_int16_decreasing(const void *a, const void *b)
{
  return ( *(int16_t *)b - *(int16_t *)a );
}

int
gal_qsort_int16_increasing(const void *a, const void *b)
{
  return ( *(int16_t *)a - *(int16_t *)b );
}

int
gal_qsort_uint32_decreasing(const void *a, const void *b)
{
  return ( *(uint32_t *)b - *(uint32_t *)a );
}

int
gal_qsort_uint32_increasing(const void *a, const void *b)
{
  return ( *(uint32_t *)a - *(uint32_t *)b );
}

int
gal_qsort_int32_decreasing(const void *a, const void *b)
{
  return ( *(int32_t *)b - *(int32_t *)a );
}

int
gal_qsort_int32_increasing(const void *a, const void *b)
{
  return ( *(int32_t *)a - *(int32_t *)b );
}

int
gal_qsort_uint64_decreasing(const void *a, const void *b)
{
  return ( *(uint64_t *)b - *(uint64_t *)a );
}

int
gal_qsort_uint64_increasing(const void *a, const void *b)
{
  return ( *(uint64_t *)a - *(uint64_t *)b );
}


int
gal_qsort_int64_decreasing(const void *a, const void *b)
{
  return ( *(int64_t *)b - *(int64_t *)a );
}

int
gal_qsort_int64_increasing(const void *a, const void *b)
{
  return ( *(int64_t *)a - *(int64_t *)b );
}

int
gal_qsort_float32_decreasing(const void *a, const void *b)
{
  float ta=*(float*)a;
  float tb=*(float*)b;
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_float32_increasing(const void *a, const void *b)
{
  float ta=*(float*)a;
  float tb=*(float*)b;
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_float64_decreasing(const void *a, const void *b)
{
  double ta=*(double*)a;
  double tb=*(double*)b;
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_float64_increasing(const void *a, const void *b)
{
  double ta=*(double*)a;
  double tb=*(double*)b;
  return (ta > tb) - (ta < tb);
}
