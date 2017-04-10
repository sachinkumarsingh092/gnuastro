/*********************************************************************
dimension -- Functions for multi-dimensional operations.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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
#include <stdlib.h>

#include <gnuastro/dimension.h>





/************************************************************************/
/********************             Info             **********************/
/************************************************************************/
size_t
gal_dimension_total_size(size_t ndim, size_t *dsize)
{
  size_t i, num=1;
  for(i=0;i<ndim;++i) num *= dsize[i];
  return num;
}





/* Calculate the values necessary to increment/decrement along each
   dimension of a dataset with size `dsize'. */
size_t *
gal_dimension_increment(size_t ndim, size_t *dsize)
{
  int i;
  size_t *out=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);

  /* Along the fastest dimension, it is 1. */
  out[ndim-1]=1;

  /* For the rest of the dimensions, it is the multiple of the faster
     dimension's length and the value for the previous dimension. */
  if(ndim>1)
    for(i=ndim-2;i>=0;--i)
      out[i]=dsize[i+1]*out[i+1];

  /* Return the allocated array. */
  return out;
}





size_t
gal_dimension_num_neighbors(size_t ndim)
{
  if(ndim)
    return pow(3, ndim)-1;
  else
    error(EXIT_FAILURE, 0, "ndim cannot be zero in "
          "`gal_dimension_num_neighbors'");
  return 0;
}


















/************************************************************************/
/********************          Coordinates         **********************/
/************************************************************************/
void
gal_dimension_add_coords(size_t *c1, size_t *c2, size_t *out, size_t ndim)
{
  size_t *end=c1+ndim;
  do *out++ = *c1++ + *c2++; while(c1<end);
}





/* Return the index of an element from its coordinates. The index is the
   position in the contiguous array (assuming it is a 1D arrray). */
size_t
gal_dimension_coord_to_index(size_t ndim, size_t *dsize, size_t *coord)
{
  size_t i, d, ind=0, in_all_faster_dim;

  switch(ndim)
    {
    case 0:
      error(EXIT_FAILURE, 0, "`gal_dimension_coord_to_index' doesn't accept "
            "0 dimensional arrays");

    case 1:
      ind=coord[0];
      break;

    case 2:
      ind=coord[0]*dsize[1]+coord[1];
      break;

    default:
      for(d=0;d<ndim;++d)
        {
          /* First, find the number of elements in all dimensions faster
             than this one. */
          in_all_faster_dim=1;
          for(i=d+1;i<ndim;++i)
            in_all_faster_dim *= dsize[i];

          /* Multiply it by the coordinate value of this dimension and add
             to the index. */
          ind += coord[d] * in_all_faster_dim;
        }
    }

  /* Return the derived index. */
  return ind;
}





/* You know the index (`ind') of a point/tile in an n-dimensional (`ndim')
   array which has `dsize[i]' elements along dimension `i'. You want to
   know the coordinates of that point along each dimension. The output is
   not actually returned, it must be allocated (`ndim' elements) before
   calling this function. This function will just fill it. The reason for
   this is that this function will often be called with a loop and a single
   allocated space would be enough for the whole loop. */
void
gal_dimension_index_to_coord(size_t ind, size_t ndim, size_t *dsize,
                             size_t *coord)
{
  size_t d, *dinc;

  switch(ndim)
    {
    case 0:
      error(EXIT_FAILURE, 0, "a 0-dimensional dataset is not defined in "
            "`gal_dimension_ind_to_coord'");

    /* One dimensional dataset. */
    case 1:
      coord[0] = ind;
      break;

    /* 2D dataset. */
    case 2:
      coord[0] = ind / dsize[1];
      coord[1] = ind % dsize[1];
      break;

    /* Higher dimensional datasets. */
    default:
      /* Set the incrementation values for each dimension. */
      dinc=gal_dimension_increment(ndim, dsize);

      /* We start with the slowest dimension (first in the C standard) and
         continue until (but not including) the fastest dimension. This is
         because except for the fastest (coniguous) dimension, the other
         coordinates can be found by division. */
      for(d=0;d<ndim;++d)
        {
          /* Set the coordinate value for this dimension. */
          coord[d] = ind / dinc[d];

          /* Replace the index with its remainder with the number of
             elements in all faster dimensions. */
          ind  %= dinc[d];
        }

      /* Clean up. */
      free(dinc);
    }
}




















/************************************************************************/
/********************           Distances          **********************/
/************************************************************************/
size_t
gal_dimension_dist_manhattan(size_t *a, size_t *b, size_t ndim)
{
  size_t i, out=0;
  for(i=0;i<ndim;++i) out += labs(a[i]-b[i]);
  return out;
}
