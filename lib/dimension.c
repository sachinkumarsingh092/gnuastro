/*********************************************************************
dimension -- Functions for multi-dimensional operations.
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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/wcs.h>
#include <gnuastro/pointer.h>
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





int
gal_dimension_is_different(gal_data_t *first, gal_data_t *second)
{
  size_t i;

  /* First make sure that the dimensionality is the same. */
  if(first->ndim!=second->ndim)
    return 1;

  /* If the sizes are not zero, check if each dimension also has the same
     length. */
  if(first->size==0 && first->size==second->size)
    return 0;
  else
    for(i=0;i<first->ndim;++i)
      if( first->dsize[i] != second->dsize[i] )
        return 1;

  /* If it got to here, we know the dimensions have the same length. */
  return 0;
}





/* Calculate the values necessary to increment/decrement along each
   dimension of a dataset with size 'dsize'. */
size_t *
gal_dimension_increment(size_t ndim, size_t *dsize)
{
  int i;
  size_t *out=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__, "out");

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
    error(EXIT_FAILURE, 0, "%s: ndim cannot be zero", __func__);
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
      error(EXIT_FAILURE, 0, "%s: doesn't accept 0 dimensional arrays",
            __func__);

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





/* You know the index ('ind') of a point/tile in an n-dimensional ('ndim')
   array which has 'dsize[i]' elements along dimension 'i'. You want to
   know the coordinates of that point along each dimension. The output is
   not actually returned, it must be allocated ('ndim' elements) before
   calling this function. This function will just fill it. The reason for
   this is that this function will often be called with a loop and a single
   allocated space would be enough for the whole loop. */
void
gal_dimension_index_to_coord(size_t index, size_t ndim, size_t *dsize,
                             size_t *coord)
{
  size_t d, *dinc;

  switch(ndim)
    {
    case 0:
      error(EXIT_FAILURE, 0, "%s: a 0-dimensional dataset is not defined",
            __func__);

    /* One dimensional dataset. */
    case 1:
      coord[0] = index;
      break;

    /* 2D dataset. */
    case 2:
      coord[0] = index / dsize[1];
      coord[1] = index % dsize[1];
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
          coord[d] = index / dinc[d];

          /* Replace the index with its remainder with the number of
             elements in all faster dimensions. */
          index  %= dinc[d];
        }

      /* Clean up. */
      free(dinc);
    }
}




















/************************************************************************/
/********************           Distances          **********************/
/************************************************************************/
float
gal_dimension_dist_manhattan(size_t *a, size_t *b, size_t ndim)
{
  size_t i, out=0;
  for(i=0;i<ndim;++i) out += (a[i] > b[i]) ? (a[i]-b[i]) : (b[i]-a[i]);
  return out;
}





float
gal_dimension_dist_radial(size_t *a, size_t *b, size_t ndim)
{
  size_t i, out=0;
  for(i=0;i<ndim;++i) out += (a[i]-b[i])*(a[i]-b[i]);
  return sqrt(out);
}




















/************************************************************************/
/********************    Collapsing a dimension    **********************/
/************************************************************************/
enum dimension_collapse_operation
{
 DIMENSION_COLLAPSE_INVALID,    /* ==0 by C standard. */

 DIMENSION_COLLAPSE_SUM,
 DIMENSION_COLLAPSE_MAX,
 DIMENSION_COLLAPSE_MIN,
 DIMENSION_COLLAPSE_MEAN,
 DIMENSION_COLLAPSE_NUMBER,
};





static gal_data_t *
dimension_collapse_sanity_check(gal_data_t *in, gal_data_t *weight,
                                size_t c_dim, int hasblank, size_t *cnum,
                                double **warr)
{
  gal_data_t *wht=NULL;

  /* The requested dimension to collapse cannot be larger than the input's
     number of dimensions. */
  if( c_dim > (in->ndim-1) )
    error(EXIT_FAILURE, 0, "%s: the input has %zu dimension(s), but you have "
          "asked to collapse dimension %zu", __func__, in->ndim, c_dim);

  /* If there is no blank value, there is no point in calculating the
     number of points in each collapsed dataset (when necessary). In that
     case, 'cnum!=0'. */
  if(hasblank==0)
    *cnum=in->dsize[c_dim];

  /* Weight sanity checks */
  if(weight)
    {
      if( weight->ndim!=1 )
        error(EXIT_FAILURE, 0, "%s: the weight dataset has %zu dimensions, "
              "it must be one-dimensional", __func__, weight->ndim);
      if( in->dsize[c_dim]!=weight->size )
        error(EXIT_FAILURE, 0, "%s: the weight dataset has %zu elements, "
              "but the input dataset has %zu elements in dimension %zu",
              __func__, weight->size, in->dsize[c_dim], c_dim);
      wht = ( weight->type == GAL_TYPE_FLOAT64
              ? weight
              : gal_data_copy_to_new_type(weight, GAL_TYPE_FLOAT64) );
      *warr = wht->array;
    }

  /* Return the weight data structure. */
  return wht;
}




/* Set the collapsed output sizes. */
static void
dimension_collapse_sizes(gal_data_t *in, size_t c_dim, size_t *outndim,
                         size_t *outdsize)
{
  size_t i, a=0;

  if(in->ndim==1)
    *outndim=outdsize[0]=1;
  else
    {
      *outndim=in->ndim-1;
      for(i=0;i<in->ndim;++i)
        if(i!=c_dim) outdsize[a++]=in->dsize[i];
    }
}





/* Depending on the operator, write the result into the output. */
#define COLLAPSE_WRITE(OIND,IIND) {                                     \
    /* Sum */                                                           \
    if(farr)                                                            \
      farr[ OIND ] += (warr ? warr[w] : 1) * inarr[ IIND ];             \
                                                                        \
    /* Number */                                                        \
    if(iarr)                                                            \
      {                                                                 \
        if(num->type==GAL_TYPE_UINT8) iarr[ OIND ] = 1;                 \
        else                        ++iarr[ OIND ];                     \
      }                                                                 \
                                                                        \
    /* Sum of weights. */                                               \
    if(wsumarr) wsumarr[ OIND ] += warr[w];                             \
                                                                        \
    /* Minimum or maximum. */                                           \
    if(mmarr)                                                           \
      mmarr[OIND] = ( max1_min0                                         \
                      ? (mmarr[OIND]>inarr[IIND]?mmarr[OIND]:inarr[IIND]) \
                      : (mmarr[OIND]<inarr[IIND]?mmarr[OIND]:inarr[IIND]) ); \
  }


/* Deal properly with blanks. */
#define COLLAPSE_CHECKBLANK(OIND,IIND) {                                \
    if(hasblank)                                                        \
      {                                                                 \
        if(B==B) /* An integer type: blank can be checked with '=='. */ \
          {                                                             \
            if( inarr[IIND] != B )           COLLAPSE_WRITE(OIND,IIND); \
          }                                                             \
        else     /* A floating point type where NAN != NAN. */          \
          {                                                             \
            if( inarr[IIND] == inarr[IIND] ) COLLAPSE_WRITE(OIND,IIND); \
          }                                                             \
      }                                                                 \
    else                                     COLLAPSE_WRITE(OIND,IIND); \
  }



#define COLLAPSE_DIM(IT) {                                              \
    IT m, B, *inarr=in->array, *mmarr=minmax?minmax->array:NULL;        \
    if(hasblank) gal_blank_write(&B, in->type);                         \
                                                                        \
    /* Initialize the array for minimum or maximum. */                  \
    if(mmarr)                                                           \
      {                                                                 \
        if(max1_min0) gal_type_min(in->type, &m);                       \
        else          gal_type_max(in->type, &m);                       \
        for(i=0;i<minmax->size;++i) mmarr[i]=m;                         \
      }                                                                 \
                                                                        \
    /* Collapse the dataset. */                                         \
    switch(in->ndim)                                                    \
      {                                                                 \
      /* 1D input dataset. */                                           \
      case 1:                                                           \
        for(i=0;i<in->dsize[0];++i)                                     \
          {                                                             \
            if(weight) w=i;                                             \
            COLLAPSE_CHECKBLANK(0,i);                                   \
          }                                                             \
        break;                                                          \
                                                                        \
      /* 2D input dataset. */                                           \
      case 2:                                                           \
        for(i=0;i<in->dsize[0];++i)                                     \
          for(j=0;j<in->dsize[1];++j)                                   \
            {                                                           \
              /* In a more easy to understand format:                   \
                 dim==0 --> a=j;                                        \
                 dim==1 --> a=i; */                                     \
              a = c_dim==0 ? j : i;                                     \
              if(weight) w = c_dim == 0 ? i : j;                        \
              COLLAPSE_CHECKBLANK(a, i*in->dsize[1] + j);               \
            }                                                           \
        break;                                                          \
                                                                        \
      /* 3D input dataset. */                                           \
      case 3:                                                           \
        slice=in->dsize[1]*in->dsize[2];                                \
        for(i=0;i<in->dsize[0];++i)                                     \
          for(j=0;j<in->dsize[1];++j)                                   \
            for(k=0;k<in->dsize[2];++k)                                 \
              {                                                         \
                /* In a more easy to understand format:                 \
                   dim==0 --> a=j; b=k;                                 \
                   dim==1 --> a=i; b=k;                                 \
                   dim==2 --> a=i; b=j;   */                            \
                a = c_dim==0 ? j : i;                                   \
                b = c_dim==2 ? j : k;                                   \
                if(weight) w = c_dim==0 ? i : (c_dim==1 ? j : k);       \
                COLLAPSE_CHECKBLANK(a*outdsize[1]+b,                    \
                                    i*slice + j*in->dsize[2] + k);      \
              }                                                         \
        break;                                                          \
                                                                        \
        /* Input dataset's dimensionality not yet supported. */         \
      default:                                                          \
        error(EXIT_FAILURE, 0, "%s: %zu-dimensional datasets not yet "  \
              "supported, please contact us at %s to add this feature", \
              __func__, in->ndim, PACKAGE_BUGREPORT);                   \
      }                                                                 \
                                                                        \
    /* For minimum or maximum, elements with no input must be blank. */ \
    if(mmarr && iarr)                                                   \
      for(i=0;i<minmax->size;++i) if(iarr[i]==0) mmarr[i]=B;            \
  }





gal_data_t *
gal_dimension_collapse_sum(gal_data_t *in, size_t c_dim, gal_data_t *weight)
{
  int max1_min0=0;
  double *wsumarr=NULL;
  uint8_t *ii, *iarr=NULL;
  size_t a, b, i, j, k, w=-1, cnum=0;
  size_t outdsize[10], slice, outndim;
  int hasblank=gal_blank_present(in, 0);
  double *dd, *df, *warr=NULL, *farr=NULL;
  gal_data_t *sum=NULL, *wht=NULL, *num=NULL, *minmax=NULL;

  /* Basic sanity checks. */
  wht=dimension_collapse_sanity_check(in, weight, c_dim, hasblank,
                                      &cnum, &warr);

  /* Set the size of the collapsed output. */
  dimension_collapse_sizes(in, c_dim, &outndim, outdsize);

  /* Allocate the sum (output) dataset. */
  sum=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, outndim, outdsize, in->wcs,
                     1, in->minmapsize, in->quietmmap, NULL, NULL, NULL);

  /* The number dataset (when there are blank values).*/
  if(hasblank)
    num=gal_data_alloc(NULL, GAL_TYPE_INT8, outndim, outdsize, NULL,
                       1, in->minmapsize, in->quietmmap, NULL, NULL, NULL);

  /* Set the array pointers. */
  if(sum) farr=sum->array;
  if(num) iarr=num->array;

  /* Parse the dataset. */
  switch(in->type)
    {
    case GAL_TYPE_UINT8:     COLLAPSE_DIM( uint8_t  );   break;
    case GAL_TYPE_INT8:      COLLAPSE_DIM( int8_t   );   break;
    case GAL_TYPE_UINT16:    COLLAPSE_DIM( uint16_t );   break;
    case GAL_TYPE_INT16:     COLLAPSE_DIM( int16_t  );   break;
    case GAL_TYPE_UINT32:    COLLAPSE_DIM( uint32_t );   break;
    case GAL_TYPE_INT32:     COLLAPSE_DIM( int32_t  );   break;
    case GAL_TYPE_UINT64:    COLLAPSE_DIM( uint64_t );   break;
    case GAL_TYPE_INT64:     COLLAPSE_DIM( int64_t  );   break;
    case GAL_TYPE_FLOAT32:   COLLAPSE_DIM( float    );   break;
    case GAL_TYPE_FLOAT64:   COLLAPSE_DIM( double   );   break;
    default:
      error(EXIT_FAILURE, 0, "%s: type value (%d) not recognized",
            __func__, in->type);
    }

  /* If 'num' is zero on any element, set its sum to NaN. */
  if(num)
    {
      ii = num->array;
      df = (dd=sum->array) + sum->size;
      do if(*ii++==0) *dd=NAN; while(++dd<df);
    }

  /* Remove the respective dimension in the WCS structure also (if any
     exists). Note that 'sum->ndim' has already been changed. So we'll use
     'in->wcs'. */
  gal_wcs_remove_dimension(sum->wcs, in->ndim-c_dim);

  /* Clean up and return. */
  if(wht!=weight) gal_data_free(wht);
  if(num) gal_data_free(num);
  return sum;
}





gal_data_t *
gal_dimension_collapse_mean(gal_data_t *in, size_t c_dim,
                            gal_data_t *weight)
{
  int max1_min0=0;
  double wsum=NAN;
  double *wsumarr=NULL;
  int32_t *ii, *iarr=NULL;
  size_t a, b, i, j, k, w=-1, cnum=0;
  size_t outdsize[10], slice, outndim;
  int hasblank=gal_blank_present(in, 0);
  double *dd, *dw, *df, *warr=NULL, *farr=NULL;
  gal_data_t *sum=NULL, *wht=NULL, *num=NULL, *minmax=NULL;


  /* Basic sanity checks. */
  wht=dimension_collapse_sanity_check(in, weight, c_dim, hasblank,
                                      &cnum, &warr);

  /* Set the size of the collapsed output. */
  dimension_collapse_sizes(in, c_dim, &outndim, outdsize);

  /* The sum array. */
  sum=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, outndim, outdsize, in->wcs,
                     1, in->minmapsize, in->quietmmap, NULL, NULL, NULL);

  /* If a weighted mean is requested. */
  if( weight )
    {
      /* There are blank values, so we'll need to keep the sums of the
         weights for each collapsed dimension */
      if( hasblank )
        wsumarr=gal_pointer_allocate(GAL_TYPE_FLOAT64, sum->size, 1,
                                     __func__, "wsumarr");

      /* There aren't any blank values, so one summation over the
         weights is enough to calculate the weighted mean. */
      else
        {
          wsum=0.0f;
          df=(dd=weight->array)+weight->size;
          do wsum += *dd++; while(dd<df);
        }
    }
  /* No weight is given, so we'll need the number of elements. */
  else if( hasblank )
    num=gal_data_alloc(NULL, GAL_TYPE_INT32, outndim, outdsize, NULL,
                       1, in->minmapsize, in->quietmmap, NULL, NULL, NULL);

  /* Set the array pointers. */
  if(sum) farr=sum->array;
  if(num) iarr=num->array;

  /* Parse the dataset. */
  switch(in->type)
    {
    case GAL_TYPE_UINT8:     COLLAPSE_DIM( uint8_t  );   break;
    case GAL_TYPE_INT8:      COLLAPSE_DIM( int8_t   );   break;
    case GAL_TYPE_UINT16:    COLLAPSE_DIM( uint16_t );   break;
    case GAL_TYPE_INT16:     COLLAPSE_DIM( int16_t  );   break;
    case GAL_TYPE_UINT32:    COLLAPSE_DIM( uint32_t );   break;
    case GAL_TYPE_INT32:     COLLAPSE_DIM( int32_t  );   break;
    case GAL_TYPE_UINT64:    COLLAPSE_DIM( uint64_t );   break;
    case GAL_TYPE_INT64:     COLLAPSE_DIM( int64_t  );   break;
    case GAL_TYPE_FLOAT32:   COLLAPSE_DIM( float    );   break;
    case GAL_TYPE_FLOAT64:   COLLAPSE_DIM( double   );   break;
    default:
      error(EXIT_FAILURE, 0, "%s: type value (%d) not recognized",
            __func__, in->type);
    }

  /* If 'num' is zero on any element, set its sum to NaN. */
  if(num)
    {
      ii = num->array;
      df = (dd=sum->array) + sum->size;
      do if(*ii++==0) *dd=NAN; while(++dd<df);
    }

  /* Divide the sum by the number. */
  df = (dd=sum->array) + sum->size;
  if(weight)
    {
      if(hasblank) { dw=wsumarr;  do *dd /= *dw++; while(++dd<df); }
      else                        do *dd /= wsum;  while(++dd<df);
    }
  else
    if(num) { ii = num->array;    do *dd /= *ii++; while(++dd<df); }
    else                          do *dd /= cnum;  while(++dd<df);

  /* Correct the WCS, clean up and return. */
  gal_wcs_remove_dimension(sum->wcs, in->ndim-c_dim);
  if(wht!=weight) gal_data_free(wht);
  if(wsumarr) free(wsumarr);
  gal_data_free(num);
  return sum;
}





gal_data_t *
gal_dimension_collapse_number(gal_data_t *in, size_t c_dim)
{
  int max1_min0=0;
  double *wsumarr=NULL;
  double *warr=NULL, *farr=NULL;
  int32_t *ii, *iif, *iarr=NULL;
  size_t a, b, i, j, k, w, cnum=0;
  size_t outdsize[10], slice, outndim;
  int hasblank=gal_blank_present(in, 0);
  gal_data_t *weight=NULL, *wht=NULL, *num=NULL, *minmax=NULL;

  /* Basic sanity checks. */
  wht=dimension_collapse_sanity_check(in, weight, c_dim, hasblank,
                                      &cnum, &warr);

  /* Set the size of the collapsed output. */
  dimension_collapse_sizes(in, c_dim, &outndim, outdsize);

  /* The number dataset (when there are blank values).*/
  num=gal_data_alloc(NULL, GAL_TYPE_INT32, outndim, outdsize, in->wcs,
                     1, in->minmapsize, in->quietmmap, NULL, NULL, NULL);

  /* Set the array pointers. */
  iarr=num->array;

  /* Parse the input dataset (if necessary). */
  if(hasblank)
    switch(in->type)
      {
      case GAL_TYPE_UINT8:     COLLAPSE_DIM( uint8_t  );   break;
      case GAL_TYPE_INT8:      COLLAPSE_DIM( int8_t   );   break;
      case GAL_TYPE_UINT16:    COLLAPSE_DIM( uint16_t );   break;
      case GAL_TYPE_INT16:     COLLAPSE_DIM( int16_t  );   break;
      case GAL_TYPE_UINT32:    COLLAPSE_DIM( uint32_t );   break;
      case GAL_TYPE_INT32:     COLLAPSE_DIM( int32_t  );   break;
      case GAL_TYPE_UINT64:    COLLAPSE_DIM( uint64_t );   break;
      case GAL_TYPE_INT64:     COLLAPSE_DIM( int64_t  );   break;
      case GAL_TYPE_FLOAT32:   COLLAPSE_DIM( float    );   break;
      case GAL_TYPE_FLOAT64:   COLLAPSE_DIM( double   );   break;
      default:
        error(EXIT_FAILURE, 0, "%s: type value (%d) not recognized",
              __func__, in->type);
      }
  else
    {
      iif=(ii=num->array)+num->size;
      do *ii++ = cnum; while(ii<iif);
    }

  /* Remove the respective dimension in the WCS structure also (if any
     exists). Note that 'sum->ndim' has already been changed. So we'll use
     'in->wcs'. */
  gal_wcs_remove_dimension(num->wcs, in->ndim-c_dim);

  /* Return. */
  if(wht!=weight) gal_data_free(wht);
  return num;
}





gal_data_t *
gal_dimension_collapse_minmax(gal_data_t *in, size_t c_dim, int max1_min0)
{
  uint8_t *iarr=NULL;
  double *wsumarr=NULL;
  double *warr=NULL, *farr=NULL;
  size_t a, b, i, j, k, w, cnum=0;
  size_t outdsize[10], slice, outndim;
  int hasblank=gal_blank_present(in, 0);
  gal_data_t *weight=NULL, *wht=NULL, *num=NULL, *minmax=NULL;

  /* Basic sanity checks. */
  wht=dimension_collapse_sanity_check(in, weight, c_dim, hasblank,
                                      &cnum, &warr);

  /* Set the size of the collapsed output. */
  dimension_collapse_sizes(in, c_dim, &outndim, outdsize);

  /* Allocate the necessary datasets. If there are blank pixels, we'll need
     to count how many elements whent into the calculation so we can set
     them to blank. */
  minmax=gal_data_alloc(NULL, in->type, outndim, outdsize, in->wcs,
                        0, in->minmapsize, in->quietmmap, NULL, NULL, NULL);
  if(hasblank)
    {
      num=gal_data_alloc(NULL, GAL_TYPE_UINT8, outndim, outdsize, in->wcs,
                         1, in->minmapsize, in->quietmmap, NULL, NULL, NULL);
      iarr=num->array;
    }

  /* Parse the input dataset (if necessary). */
  switch(in->type)
    {
    case GAL_TYPE_UINT8:     COLLAPSE_DIM( uint8_t  );   break;
    case GAL_TYPE_INT8:      COLLAPSE_DIM( int8_t   );   break;
    case GAL_TYPE_UINT16:    COLLAPSE_DIM( uint16_t );   break;
    case GAL_TYPE_INT16:     COLLAPSE_DIM( int16_t  );   break;
    case GAL_TYPE_UINT32:    COLLAPSE_DIM( uint32_t );   break;
    case GAL_TYPE_INT32:     COLLAPSE_DIM( int32_t  );   break;
    case GAL_TYPE_UINT64:    COLLAPSE_DIM( uint64_t );   break;
    case GAL_TYPE_INT64:     COLLAPSE_DIM( int64_t  );   break;
    case GAL_TYPE_FLOAT32:   COLLAPSE_DIM( float    );   break;
    case GAL_TYPE_FLOAT64:   COLLAPSE_DIM( double   );   break;
    default:
      error(EXIT_FAILURE, 0, "%s: type value (%d) not recognized",
            __func__, in->type);
    }

  /* Remove the respective dimension in the WCS structure also (if any
     exists). Note that 'sum->ndim' has already been changed. So we'll use
     'in->wcs'. */
  gal_wcs_remove_dimension(minmax->wcs, in->ndim-c_dim);

  /* Clean up and return. */
  if(wht!=weight) gal_data_free(wht);
  if(num) gal_data_free(num);
  return minmax;
}




















/************************************************************************/
/********************             Other            **********************/
/************************************************************************/
size_t
gal_dimension_remove_extra(size_t ndim, size_t *dsize, struct wcsprm *wcs)
{
  size_t i, j;

  for(i=0;i<ndim;++i)
    if(dsize[i]==1)
      {
        /* Correct the WCS. */
        if(wcs) gal_wcs_remove_dimension(wcs, ndim-i);

        /* Shift all subsequent dimensions to replace this one. */
        for(j=i;j<ndim-1;++j) dsize[j]=dsize[j+1];

        /* Decrement the 'i' and the total number of dimension. */
        --i;
        --ndim;
      }

  /* Return the number of dimensions. */
  return ndim;
}
