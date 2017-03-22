/*********************************************************************
multidim -- Functions for multi-dimensional operations.
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
#ifndef __GAL_MULTIDIM_H__
#define __GAL_MULTIDIM_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>

/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */

/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */





/************************************************************************/
/********************             Info             **********************/
/************************************************************************/
size_t
gal_dimension_total_size(size_t ndim, size_t *dsize);

size_t *
gal_dimension_increment(size_t ndim, size_t *dsize);

size_t
gal_dimension_num_neighbors(size_t ndim);





/************************************************************************/
/********************          Coordinates         **********************/
/************************************************************************/
void
gal_dimension_add_coords(size_t *c1, size_t *c2, size_t *out, size_t ndim);

size_t
gal_dimension_coord_to_index(size_t ndim, size_t *dsize, size_t *coord);

void
gal_dimension_index_to_coord(size_t ind, size_t ndim, size_t *dsize,
                             size_t *coord);





/************************************************************************/
/********************           Distances          **********************/
/************************************************************************/
size_t
gal_dimension_dist_manhattan(size_t *a, size_t *b, size_t ndim);





/************************************************************************/
/********************          Neighbors           **********************/
/************************************************************************/
/* Purpose
   -------

   This macro will allow you to do a fixed operation on the neighbors of an
   eleemnt. The identifier for the element is its dimension-agnostic index
   (distance from start of array). It is defined as a macro (and not a
   function) it is often necessary to loop over a very large number of
   pixels/indexs and the number of neighbors differs (in different
   dimensions and on the edges of the image).

   Usage
   -----

   The inputs are:

       `index': The index of the desired point.

       `ndim':  The number of dimensions in the dataset.

       `dsize': The size of the dataset along each dimension (C order).

       `connectivity': The connectivity of the neighbors. This is a single
                integer with a value between `1' and `ndim' (it is
                irrelevant for a 1D dataset). Below, the case for 3D data
                is shown. Each `dn' corresponds to `dinc[n]'. See `dinc'
                below for the `dinc' values.

                   `1': At most ONE addition/subtraction. In 1D: only the
                        first line should be used. In 2D, this is 4
                        connectivity. In 3D all the elements that share an
                        2D-face with the cube at index `i'.

                                       i - d0     i + d0        (1D, 2D & 3D)
                                       i - d1     i + d1        (2D & 3D)
                                       i - d2     i + d2        (3D)

                   `2': At most TWO additions/subtractions. In 2D: 8
                        connectivity. In 3D: all elements that share a 1D
                        edge with the cube at index `i'.

                                  i - d0 - d1     i - d0 + d1   (2D & 3D)
                                  i + d0 - d1     i + d0 + d1   (2D & 3D)
                                  i - d0 - d2     i - d0 + d2   (3D)
                                  i + d0 - d2     i + d0 + d2   (3D)
                                  i - d1 - d2     i - d1 + d2   (3D)
                                  i + d1 - d2     i + d1 + d2   (3D)

                   `3': At most THREE additions/subtractions (only for 3D
                        and higher dimensions). All cubes that share a 0-D
                        vertex with the cube at index `i'.

                             i - d0 - d1 - d2     i - d0 - d1 + d2
                             i - d0 + d1 - d2     i - d0 + d1 + d2
                             i + d0 - d1 - d2     i + d0 - d1 + d2
                             i + d0 + d1 - d2     i + d0 + d1 + d2

        `dinc': An array keeping the length necessary to increment along
                each dimension. You can make this array with the following
                function:

                  size_t *dinc=gal_dimension_increment(ndim, dsize);

                Don't forget to free it afterwards.

        `operation': Any C operation. `nind' is a `size_t' type variable
                that is defined by this macro and will have the index of
                each neighbor. You can use this `nind' for any processing
                that you like on the neighbor. Note that `op' will be
                repeated the number of times there is a neighbor.


   Implementation
   --------------

   To be most efficient (avoid as many `if's as possible), we will start
   parsing the neighbors from the fastest dimension. When-ever the element
   is on the edge of the dataset in any dimension, we will store it in a
   bit-wise array (one bit for each dimension, to mark if it is on the edge
   (either side 1, or in the middle 0). Far many more elements will be in
   the middle, so there is no problem to check the edge. The good thing
   with a bit-array is that it can take one register on the CPU and stay
   there until the end. But if we want to have an array of values, multiple
   registers will be necessary.

   The bit information is in two two-byte spaces, so in theory, this works
   for 16 dimensions.
*/
#define GAL_DIMENSION_NEIGHBOR_OP(index, ndim, dsize, connectivity,     \
                                  dinc, operation) {                    \
    uint32_t bitstr=0;                                                  \
    size_t nind, ind=index;                                             \
    uint8_t D, *is_start, *is_end, *is_edge, one=1;                     \
                                                                        \
    /* Initialize the start/end. */                                     \
    is_start=(uint8_t *)(&bitstr);                                      \
    is_end=(uint8_t *)(&bitstr)+1;                                      \
                                                                        \
    /* Start with the slowest dimension and see if it is on the edge */ \
    /* or not, similar to `gal_dimension_index_to_coord'. In the */     \
    /* process, also fill the `connectivity==1' neighbors. */           \
    for(D=0;D<ndim;++D)                                                 \
      {                                                                 \
        /* If this dimension is only one element wide, no neighbors. */ \
        if( (dsize)[D] == 1 )                                           \
          {                                                             \
            *is_start |= 1<<D;                                          \
            *is_end   |= 1<<D;                                          \
          }                                                             \
        else                                                            \
          {                                                             \
            if( ind / (dinc)[D] )                                       \
              {                                                         \
                /* We are at the end of this dimension. */              \
                if( ind / (dinc)[D] == (dsize)[D]-1 )                   \
                  {                                                     \
                    *is_end |= one<<D;                                  \
                    nind = (index) - (dinc)[D]; {operation;};           \
                  }                                                     \
                                                                        \
                /* Middle of the dimension: both +1 and -1 possible. */ \
                else                                                    \
                  {                                                     \
                    nind = (index) - (dinc)[D]; {operation;};           \
                    nind = (index) + (dinc)[D]; {operation;};           \
                  }                                                     \
              }                                                         \
            else                                                        \
              {                                                         \
                *is_start |= one<<D;                                    \
                nind = (index) + (dinc)[D]; {operation;};               \
              }                                                         \
          }                                                             \
                                                                        \
        /* Change `ind' to the remainder of previous dimensions. */     \
        ind %= dinc[D];                                                 \
      }                                                                 \
                                                                        \
    /* We now know if the index is on the edge or not. During the */    \
    /* process above, we also finished the `connectivity==1' case. */   \
    /* So we'll just have to find the rest of the terms. */             \
    if(connectivity>1 && ndim>1)                                        \
      {                                                                 \
        /* Finalize `is_edge' (bit value 1 for respective dim.). */     \
        is_edge=(uint8_t *)(&bitstr)+2;                                 \
        *is_edge = *is_start | *is_end;                                 \
                                                                        \
        /* Shared between 2D and 3D datasets. */                        \
        if(*is_edge)                                                    \
          { /* NOTE: these are bitwise operators, not conditionals. */  \
            if( !( *is_start & ( one | one<<1 ) ) )                     \
              { nind=(index) - (dinc)[0] - (dinc)[1]; {operation;}; }   \
            if( !( *is_start & one ) && !( *is_end   & one<<1 ) )       \
              { nind=(index) - (dinc)[0] + (dinc)[1]; {operation;}; }   \
            if( !( *is_end   & one ) && !( *is_start & one<<1 ) )       \
              { nind=(index) + (dinc)[0] - (dinc)[1]; {operation;}; }   \
            if( !( *is_end   & ( one | one<<1 ) ) )                     \
              { nind=(index) + (dinc)[0] + (dinc)[1]; {operation;}; }   \
          }                                                             \
        else                                                            \
          {                                                             \
            nind=(index) - (dinc)[0] - (dinc)[1]; {operation;};         \
            nind=(index) - (dinc)[0] + (dinc)[1]; {operation;};         \
            nind=(index) + (dinc)[0] - (dinc)[1]; {operation;};         \
            nind=(index) + (dinc)[0] + (dinc)[1]; {operation;};         \
          }                                                             \
                                                                        \
        /* Only for 3D datasets. */                                     \
        if(ndim>2)                                                      \
          {                                                             \
            /* Connectivity == 2. */                                    \
            if(*is_edge)                                                \
              for(D=0;D<2;++D)                                          \
                {                                                       \
                  if( !( *is_start & ( one<<D | one<<2 ) ) )            \
                    { nind=(index) - (dinc)[D] - (dinc)[2]; {operation;}; } \
                  if( !( *is_start & one<<D ) && !( *is_end   & one<<2 ) ) \
                    { nind=(index) - (dinc)[D] + (dinc)[2]; {operation;}; } \
                  if( !( *is_end   & one<<D ) && !( *is_start & one<<2 ) ) \
                    { nind=(index) + (dinc)[D] - (dinc)[2]; {operation;}; } \
                  if( !( *is_end   & ( one<<D | one<<2 ) ) )            \
                    { nind=(index) + (dinc)[D] + (dinc)[2]; {operation;}; } \
                }                                                       \
            else                                                        \
              for(D=0;D<2;++D)                                          \
                {                                                       \
                  nind=(index) - (dinc)[D] - (dinc)[2]; {operation;};   \
                  nind=(index) - (dinc)[D] + (dinc)[2]; {operation;};   \
                  nind=(index) + (dinc)[D] - (dinc)[2]; {operation;};   \
                  nind=(index) + (dinc)[D] + (dinc)[2]; {operation;};   \
                }                                                       \
                                                                        \
            /* Connectivity == 3. */                                    \
            if(connectivity>2)                                          \
              {                                                         \
                if(*is_edge)                                            \
                  {                                                     \
                    if( !*is_start )                                    \
                      { nind=(index) - (dinc)[0] - (dinc)[1] - (dinc)[2]; \
                        {operation;}; }                                 \
                                                                        \
                    if( !(*is_start & one) && !(*is_start & one<<1)     \
                        && !(*is_end & one<<2))                         \
                      { nind=(index) - (dinc)[0] - (dinc)[1] + (dinc)[2]; \
                        {operation;}; }                                 \
                                                                        \
                    if( !(*is_start & one) && !(*is_end & one<<1)       \
                        && !(*is_start & one<<2))                       \
                      { nind=(index) - (dinc)[0] + (dinc)[1] - (dinc)[2]; \
                        {operation;}; }                                 \
                                                                        \
                    if( !(*is_start & one) && !(*is_end & one<<1)       \
                        && !(*is_end & one<<2))                         \
                      { nind=(index) - (dinc)[0] + (dinc)[1] + (dinc)[2]; \
                        {operation;}; }                                 \
                                                                        \
                    if( !(*is_end & one) && !(*is_start & one<<1)       \
                        && !(*is_start & one<<2))                       \
                      { nind=(index) + (dinc)[0] - (dinc)[1] - (dinc)[2]; \
                        {operation;}; }                                 \
                                                                        \
                    if( !(*is_end & one) && !(*is_start & one<<1)       \
                        && !(*is_end & one<<2))                         \
                      { nind=(index) + (dinc)[0] - (dinc)[1] + (dinc)[2]; \
                        {operation;}; }                                 \
                                                                        \
                    if( !(*is_end & one) && !(*is_end & one<<1)         \
                        && !(*is_start & one<<2))                       \
                      { nind=(index) + (dinc)[0] + (dinc)[1] - (dinc)[2]; \
                        {operation;}; }                                 \
                                                                        \
                    if( !*is_end )                                      \
                      { nind=(index) + (dinc)[0] + (dinc)[1] + (dinc)[2]; \
                        {operation;}; }                                 \
                  }                                                     \
                else                                                    \
                  {                                                     \
                    nind=(index) - (dinc)[0] - (dinc)[1] - (dinc)[2];   \
                    {operation;};                                       \
                    nind=(index) - (dinc)[0] - (dinc)[1] + (dinc)[2];   \
                    {operation;};                                       \
                    nind=(index) - (dinc)[0] + (dinc)[1] - (dinc)[2];   \
                    {operation;};                                       \
                    nind=(index) - (dinc)[0] + (dinc)[1] + (dinc)[2];   \
                    {operation;};                                       \
                    nind=(index) + (dinc)[0] - (dinc)[1] - (dinc)[2];   \
                    {operation;};                                       \
                    nind=(index) + (dinc)[0] - (dinc)[1] + (dinc)[2];   \
                    {operation;};                                       \
                    nind=(index) + (dinc)[0] + (dinc)[1] - (dinc)[2];   \
                    {operation;};                                       \
                    nind=(index) + (dinc)[0] + (dinc)[1] + (dinc)[2];   \
                    {operation;};                                       \
                  }                                                     \
              }                                                         \
          }                                                             \
      }                                                                 \
                                                                        \
    /* For a check. */                                                  \
    /* printf("\nEdge bit flags: "); */                                 \
    /* gal_data_bit_print_stream(&bitstr, 3); printf("\n"); */          \
  }


__END_C_DECLS    /* From C++ preparations */

#endif
