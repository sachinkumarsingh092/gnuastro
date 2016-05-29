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
