/*********************************************************************
forqsort -- Functions used by qsort to sort an array.
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

#include <stdlib.h>

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
gal_qsort_uchar_decreasing(const void * a, const void * b)
{
  return ( *(unsigned char*)b - *(unsigned char*)a );
}

int
gal_qsort_uchar_increasing(const void * a, const void * b)
{
  return ( *(unsigned char*)a - *(unsigned char*)b );
}

int
gal_qsort_char_decreasing(const void * a, const void * b)
{
  return ( *(char*)b - *(char*)a );
}

int
gal_qsort_char_increasing(const void * a, const void * b)
{
  return ( *(char*)a - *(char*)b );
}

int
gal_qsort_ushort_decreasing(const void * a, const void * b)
{
  return ( *(unsigned short*)b - *(unsigned short*)a );
}

int
gal_qsort_ushort_increasing(const void * a, const void * b)
{
  return ( *(unsigned short*)a - *(unsigned short*)b );
}

int
gal_qsort_short_decreasing(const void * a, const void * b)
{
  return ( *(short*)b - *(short*)a );
}

int
gal_qsort_short_increasing(const void * a, const void * b)
{
  return ( *(short*)a - *(short*)b );
}

int
gal_qsort_uint_decreasing(const void * a, const void * b)
{
  return ( *(unsigned int*)b - *(unsigned int*)a );
}

int
gal_qsort_uint_increasing(const void * a, const void * b)
{
  return ( *(unsigned int*)a - *(unsigned int*)b );
}

int
gal_qsort_int_decreasing(const void * a, const void * b)
{
  return ( *(int*)b - *(int*)a );
}

int
gal_qsort_int_increasing(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int
gal_qsort_ulong_decreasing(const void * a, const void * b)
{
  return ( *(unsigned long*)b - *(unsigned long*)a );
}

int
gal_qsort_ulong_increasing(const void * a, const void * b)
{
  return ( *(unsigned long*)a - *(unsigned long*)b );
}


int
gal_qsort_long_decreasing(const void * a, const void * b)
{
  return ( *(long*)b - *(long*)a );
}

int
gal_qsort_long_increasing(const void * a, const void * b)
{
  return ( *(long*)a - *(long*)b );
}

int
gal_qsort_longlong_decreasing(const void * a, const void * b)
{
  return ( *(LONGLONG*)b - *(LONGLONG*)a );
}

int
gal_qsort_longlong_increasing(const void * a, const void * b)
{
  return ( *(LONGLONG*)a - *(LONGLONG*)b );
}

int
gal_qsort_float_decreasing(const void * a, const void * b)
{
  float ta=*(float*)a;
  float tb=*(float*)b;
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_float_increasing(const void * a, const void * b)
{
  float ta=*(float*)a;
  float tb=*(float*)b;
  return (ta > tb) - (ta < tb);
}

int
gal_qsort_double_decreasing(const void * a, const void * b)
{
  double ta=*(double*)a;
  double tb=*(double*)b;
  return (tb > ta) - (tb < ta);
}

int
gal_qsort_double_increasing(const void * a, const void * b)
{
  double ta=*(double*)a;
  double tb=*(double*)b;
  return (ta > tb) - (ta < tb);
}
