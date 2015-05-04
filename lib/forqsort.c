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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdlib.h>

#include "forqsort.h"

/* Initialize the array for sorting indexs to NULL. */
float *forqsortindexarr;

int
indexfloatdecreasing(const void * a, const void * b)
{
  float ta=forqsortindexarr[ *(size_t *)a ];
  float tb=forqsortindexarr[ *(size_t *)b ];
  return (tb > ta) - (tb < ta);
}




int
intdecreasing(const void * a, const void * b)
{
  return ( *(int*)b - *(int*)a );
}

int
intincreasing(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int
floatdecreasing(const void * a, const void * b)
{
  float ta=*(float*)a;
  float tb=*(float*)b;
  return (tb > ta) - (tb < ta);
}

int
floatincreasing(const void * a, const void * b)
{
  float ta=*(float*)a;
  float tb=*(float*)b;
  return (ta > tb) - (ta < tb);
}

int
doubledecreasing(const void * a, const void * b)
{
  double ta=*(double*)a;
  double tb=*(double*)b;
  return (tb > ta) - (tb < ta);
}

int
doubleincreasing(const void * a, const void * b)
{
  double ta=*(double*)a;
  double tb=*(double*)b;
  return (ta > tb) - (ta < tb);
}
