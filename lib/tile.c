/*********************************************************************
tile -- work with tesselations over a host dataset.
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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/fits.h>
#include <gnuastro/tile.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/convolve.h>
#include <gnuastro/dimension.h>
#include <gnuastro/interpolate.h>
#include <gnuastro/permutation.h>

#include <gnuastro-internal/checkset.h>









/***********************************************************************/
/**************              Single tile              ******************/
/***********************************************************************/
/* Calculate the starting coordinates of a tile in the allocated block of
   memory. */
void
gal_tile_start_coord(gal_data_t *tile, size_t *start_coord)
{
  size_t ind, ndim=tile->ndim;
  gal_data_t *block=gal_tile_block(tile);

  /* If the input tile is actually the same as the block, then the start is
     at 0 (in all dimensions). */
  if(block==tile)
    memset(start_coord, 0, ndim*gal_type_sizeof(GAL_TYPE_SIZE_T));
  else
    {
      /* Calculate the coordinates of the first pixel of the tile. */
      ind = gal_pointer_num_between(block->array, tile->array, block->type);
      gal_dimension_index_to_coord(ind, ndim, block->dsize, start_coord);
    }
}





/* Put the starting and ending (end point is not inclusive) coordinates of
   a tile into the 'start_end' array. It is assumed that a space of
   '2*tile->ndim' has been already allocated (static or dynamic) before
   this function is called.

   'rel_block' (or relative-to-block) is only relevant when the tile has an
   intermediate tile between it and the allocated space (like a channel,
   see 'gal_tile_full_two_layers'). If it doesn't ('tile->block' points the
   allocated dataset), then the value to 'rel_block' is irrelevant.

   However, when 'tile->block' is its self a larger block and 'rel_block'
   is set to 0, then the starting and ending positions will be based on the
   position within 'tile->block', not the allocated space. */
void
gal_tile_start_end_coord(gal_data_t *tile, size_t *start_end, int rel_block)
{
  size_t *s, *sf, *h;
  gal_data_t *block=gal_tile_block(tile);
  gal_data_t *host=rel_block ? block : tile->block;
  size_t *hcoord, start_ind, ndim=tile->ndim, *end=start_end+ndim;

  /* Get the starting index. Note that for the type we need the allocated
     block dataset and can't rely on the tiles. */
  start_ind=gal_pointer_num_between(block->array, tile->array, block->type);

  /* Get the coordinates of the starting point relative to the allocated
     block. */
  gal_dimension_index_to_coord(start_ind, ndim, block->dsize, start_end);

  /* When the host is different from the block, the tile's starting
     position needs to be corrected. */
  if(host!=block)
    {
      /* Get the host's starting coordinates. */
      start_ind=gal_pointer_num_between(block->array, host->array,
                                        block->type);

      /* Temporarily put the host's coordinates in the place held for the
         ending coordinates. */
      hcoord=end;
      gal_dimension_index_to_coord(start_ind, ndim, block->dsize, hcoord);
      sf=(s=start_end)+ndim; h=hcoord; do *s++ -= *h++; while(s<sf);
    }

  /* Add the dimensions of the tile to the starting coordinate. Note that
     the ending coordinates are stored immediately after the start.*/
  gal_dimension_add_coords(start_end, tile->dsize, end, ndim);
}





/* Put the indexs of the first/start and last/end pixels (inclusive) in a
   tile into the 'start_end' array (that has two elements). It will then
   return the pointer to the start of the tile in the 'work' data
   structure. */
void *
gal_tile_start_end_ind_inclusive(gal_data_t *tile, gal_data_t *work,
                                 size_t *start_end_inc)
{
  gal_data_t *block=gal_tile_block(tile);
  size_t ndim=tile->ndim, *s, *e, *l, *sf;
  size_t *start_coord = gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0,
                                             __func__, "start_coord");
  size_t *end_coord   = gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0,
                                             __func__, "end_coord");


  /* The starting index can be found from the distance of the 'tile->array'
     pointer and 'block->array' pointer. IMPORTANT: with the type of the
     block array.  */
  start_end_inc[0]=gal_pointer_num_between(block->array, tile->array,
                                           block->type);


  /* To find the end index, we need to know the coordinates of the starting
     point in the allocated block.  */
  gal_dimension_index_to_coord(start_end_inc[0], ndim, block->dsize,
                              start_coord);


  /* 'end_coord' is one unit ahead of the last element in the tile in every
     dimension. To have less potential for bugs, we will remove that extra
     value, so we get the coordinates of the last pixel in the tile
     (inclusive). We will finally, increment that value by one to get to
     the pixel immediately outside of the tile.*/
  e=end_coord;
  l=tile->dsize;
  sf=(s=start_coord)+ndim; do *e++ = *s + *l++ - 1; while(++s<sf);


  /* Convert the (inclusive) ending point's coordinates into an index. */
  start_end_inc[1]=gal_dimension_coord_to_index(ndim, block->dsize,
                                                end_coord);


  /* For a check:
  printf("\ntile_dsize: %zu, %zu, %zu\n", tile->dsize[0], tile->dsize[1],
         tile->dsize[2]);
  printf("start_coord: %zu, %zu, %zu\n", start_coord[0], start_coord[1],
         start_coord[2]);
  printf("end_coord: %zu, %zu, %zu\n", end_coord[0], end_coord[1],
         end_coord[2]);
  printf("start_index: %zu\n", start_end_inc[0]);
  printf("end_index: %zu\n", start_end_inc[1]);
  exit(1);
  */


  /* Clean up and return the pointer in the work array that the tile starts
     from. */
  free(end_coord);
  free(start_coord);
  return gal_pointer_increment(work->array, start_end_inc[0], work->type);
}




















/***********************************************************************/
/**************           Series of tiles             ******************/
/***********************************************************************/
/* Construct a list of tile(s) given positional minimum(s) and maximum(s).
   The output is an allocated an allocated array that can later be freed
   with 'gal_data_array_free'. The minimum and maximums are assumed to be
   inclusive.

   The array keeping the minmium and maximum coordinates for each label
   will have the following format:

       | min0_d0 | min0_d1 | max0_d0 | max0_d1 | ...

                       ... | minN_d0 | minN_d1 | maxN_d0 | maxN_d1 |   */
gal_data_t *
gal_tile_series_from_minmax(gal_data_t *block, size_t *minmax, size_t number)
{
  size_t ndim=block->ndim;

  size_t *min, *max;
  size_t i, d, ind, size, width=2*ndim;
  gal_data_t *tiles=gal_data_array_calloc(number);

  /* Fill the tile information.  */
  for(i=0;i<number;++i)
    {
      /* To make things more readable. */
      min = &minmax[ i * width        ];
      max = &minmax[ i * width + ndim ];

      /* Tile types should be invalid (we shouldn't use tiles directly),
         also se the other simple values. */
      tiles[i].flag  = 0;
      tiles[i].block = block;
      tiles[i].type  = GAL_TYPE_INVALID;
      tiles[i].next  = i==number-1 ? NULL : &tiles[i+1];

      /* Set the size related constants. */
      size = 1;
      tiles[i].ndim  = ndim;
      tiles[i].dsize = gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0,
                                            __func__, "tiles[i].dsize");
      for(d=0;d<ndim;++d) size *= tiles[i].dsize[d] = max[d] - min[d] + 1;
      tiles[i].size  = size;

      /* Tile's array pointer. */
      ind=gal_dimension_coord_to_index(ndim, block->dsize, min);
      tiles[i].array = gal_pointer_increment(block->array, ind, block->type);
    }

  /* For a check (put all the objects in an extension of a test file).
  {
    gal_data_t *copy;
    for(i=0;i<number;++i)
      {
        copy=gal_data_copy(&tiles[i]);
        gal_fits_img_write(copy, "tiles.fits", NULL, NULL);
      }
  }
  */

  /* Return the final pointer. */
  return tiles;
}




















/***********************************************************************/
/**************        Allocated block of memory      ******************/
/***********************************************************************/
/* When you are working on an array, it important to know the size of the
   allocated space in each dimension. This simple function will just follow
   the block pointer and return the 'dsize' element of lowest-level
   structure. */
gal_data_t *
gal_tile_block(gal_data_t *tile)
{
  while(tile->block!=NULL) tile=tile->block;
  return tile;
}





/* Return the increment necessary to start at the next series of contiguous
   memory (fastest dimension) associated with a tile.

   1D and 2D cases are simple and need no extra explanation, but the case
   for higher dimensions can be alittle more complicated, So we will go
   over some examples. The notations below are:

       'n'     number of dimensions (same in tile and block).
       't[]'   size of the tile in each dimension.
       'b[]'   size of the allocated block in each dimension.

   It is just important to see the output of this function as an increment
   from the the last patch of contiguous memory associated with the
   tile. So when the increment number is 't[n-1]' (the first 2D slice of
   the tile has been parsed), simply incrementing by 'b[n-2] * b[n-1]' will
   take us to the last row of

  num_increment      coord         increment
  -------------      -----         ---------
         1          (...0,0,0)     b[n-1]: fastest dimension of the block.
         2          (...0,1,0)     Similar to previous
         .              .               .
         .              .               .
       t[n-2]       (...1,0,0)     (b[n-2] * b[n-1]) - ( (t[n-2]-1) * b[n-1] )
      t[n-2] + 1    (...1,1,0)      b[n-1]
         .              .               .
         .              .               .
      2 * t[n-2]    (...2,0,0)     b[n-2] * b[n-1]
      t[n-2]+1      (...2,1,0)     b[n-1]
         .              .                .
         .              .                .
   t[n-3] * t[n-2]  (..1,0,0,0)    b[n-3] * b[n-2] * b[n-1]

 */
size_t
gal_tile_block_increment(gal_data_t *block, size_t *tsize,
                         size_t num_increment, size_t *coord)
{
  size_t n=block->ndim;
  size_t *b=block->dsize, *t=tsize;
  size_t increment=GAL_BLANK_SIZE_T;

  if(n>3)
    error(EXIT_FAILURE, 0, "%s: currently only implemented for at most 3 "
          "dimensions", __func__);

  switch(n)
    {
    /* A zero-dimensional dataset is not defined. */
    case 0:
      error(EXIT_FAILURE, 0, "%s: zero dimensional input is not acceptable",
            __func__);

    /* 1D: the increment is just the tile size. */
    case 1:
      increment=t[0];
      if(coord) coord[0]+=increment;
      break;

    /* 2D: the increment is the block's number of fastest axis pixels. */
    case 2:
      increment=b[1];
      if(coord) ++coord[0];
      break;

    /* 3D: The increment depends on which dimension we are reaching. */
    case 3:
      if(num_increment % t[1])
        {
          increment = b[2];
          if(coord) ++coord[1];
        }
      else
        {
          increment=(b[1] * b[2]) - ( (t[1]-1) * b[2] );
          if(coord) { ++coord[0]; coord[1] -= t[1]-1; coord[2]=0; }
        }
      break;
    }

  /* Return the final increment value. */
  return increment;
}





/* Write a constant value for each tile into each pixel covered by the
   input tiles in an array the size of the block and return it.

   Arguments
   ---------

     'tilevalues': This must be an array that has the same number of
        elements as that in 'tilesll' and in the same order that 'tilesll'
        elements are parsed (from first to last). As a result the
        dimensionality of this array is irrelevant. Note that unlike
        'tiles', 'tilevalues' must be an array.

     'tilesll': This will be parsed as a linked list (using the 'next'
        element). Internally, it might be stored as an array, but this
        function doesn't care! The position of the tile over its block will
        be determined according to the 'block' element and the pointer of
        its 'array' as fully described in 'gnuastro/data.h'. This function
        will not pop/free the list, it will only parse it from start to
        end.

     'initialize': Initialize the allocated space with blank values before
        writing in the constant values. This can be useful when the tiles
        don't cover the full allocated block. */
gal_data_t *
gal_tile_block_write_const_value(gal_data_t *tilevalues, gal_data_t *tilesll,
                                 int withblank, int initialize)
{
  void *in;
  int type=tilevalues->type;
  size_t tile_ind, nt=0, nv=tilevalues->size;
  gal_data_t *tofill, *tile, *block=gal_tile_block(tilesll);

  /* A small sanity check. */
  for(tile=tilesll; tile!=NULL; tile=tile->next) ++nt;
  if(nt!=nv)
    error(EXIT_FAILURE, 0, "%s: the number of elements in 'tilevalues' (%zu) "
          "and 'tilesll' (%zu) must be the same", __func__, nv, nt);

  /* Allocate the output array. */
  tofill=gal_data_alloc(NULL, type, block->ndim, block->dsize, block->wcs,
                        0, block->minmapsize, block->quietmmap,
                        tilevalues->name, tilevalues->unit,
                        tilevalues->comment);

  /* If requested, initialize 'tofill', otherwise it is assumed that the
     full area of the output is covered by the tiles. */
  if(withblank || initialize) gal_blank_initialize(tofill);
  else
    {
      /* Copy the flags. */
      tofill->flag=tilevalues->flag;

      /* If we have more than one dimension, then remove the possibly
         sorted flags. */
      if(block->ndim>1)
        {
          tofill->flag &= ~GAL_DATA_FLAG_SORTED_I;
          tofill->flag &= ~GAL_DATA_FLAG_SORTED_D;
        }
    }

  /* Go over the tiles and write the values in. Recall that 'tofill' has
     the same type as 'tilevalues'. So we are using memcopy. */
  tile_ind=0;
  for(tile=tilesll; tile!=NULL; tile=tile->next)
    {
      /* Set the pointer to use as input. The 'if(o)' statement is set
         because GCC 7.1.1 complained about the possiblity of the first
         argument of 'memcpy' being NULL. Recall that 'o' is a pointer. */
      in=gal_pointer_increment(tilevalues->array, tile_ind++, type);
      GAL_TILE_PARSE_OPERATE( tile, tofill, 1, withblank, {
          if(o) memcpy(o, in, gal_type_sizeof(type));
        } );
    }

  return tofill;
}





/* Make a copy of the memory block and fill it with the index of each tile
   in 'tilesll' (counting from 0). The non-filled areas will have blank
   values. The output dataset will have a type of 'GAL_TYPE_INT32'. */
gal_data_t *
gal_tile_block_check_tiles(gal_data_t *tilesll)
{
  int32_t *arr;
  size_t i, dsize=gal_list_data_number(tilesll);
  gal_data_t *ids, *out, *block=gal_tile_block(tilesll);

  /* Allocate the array to keep the IDs of each tile. */
  ids=gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &dsize,
                     NULL, 0, block->minmapsize, block->quietmmap,
                     NULL, NULL, NULL);

  /* Put the IDs into the array. */
  arr=ids->array; for(i=0;i<dsize;++i) arr[i]=i;

  /* Make the output. */
  out=gal_tile_block_write_const_value(ids, tilesll, 0, 1);

  /* Clean up and return. */
  gal_data_free(ids);
  return out;
}





/* Return the pointer corresponding to the tile in another data
   structure (can have another type). */
void *
gal_tile_block_relative_to_other(gal_data_t *tile, gal_data_t *other)
{
  gal_data_t *block=gal_tile_block(tile);
  return gal_pointer_increment(other->array,
                               gal_pointer_num_between(block->array,
                                                       tile->array,
                                                       block->type),
                               other->type);
}





/* To use within 'gal_tile_full_blank_flag'. */
static void *
tile_block_blank_flag(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  gal_data_t *tile_ll=(gal_data_t *)(tprm->params);

  size_t i;
  gal_data_t *tile;

  /* Check all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      tile=&tile_ll[ tprm->indexs[i] ];
      gal_blank_present(tile, 1);
    }

  /* Wait for all the other threads to finish. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* Update the blank flag on the tiles within the list of input tiles. */
void
gal_tile_block_blank_flag(gal_data_t *tile_ll, size_t numthreads)
{
  /* Go over all the tiles and update their blank flag. */
  gal_threads_spin_off(tile_block_blank_flag, tile_ll,
                       gal_list_data_number(tile_ll), numthreads);
}




















/***********************************************************************/
/**************           Tile full dataset         ********************/
/***********************************************************************/
/* The user's specified tile size might not be an exact multiple of the
   parent's size. This function is useful in such cases. It will give the
   starting tile's size along each dimension.

   The most simplistic way to manage the tiles is to put the regular tiles
   at the start. The line below can be the length along any dimension, and
   the tile size along that dimension.

        | tile size | tile size | tile size | tile size | remainder
        |           |           |           |           |       |
        ---------------------------------------------------------

   The remainder of the scenario above will always be smaller than 'tile
   size' (can be even 1-pixel wide). So, we will merge the first tile size
   with the remainder.  In this way, the size of the first tile will always
   be between between one and two times the size of the regular tile:

        | first tile        | tile size | tile size | tile size |
        |                   |           |           |           |
        ---------------------------------------------------------

   When there is only a small remainder (for example one or two pixels),
   then this layout is fine. But when the remainder is significant compared
   to the regular tile size (like the example above), then it will make
   more sense to cut the first tile into two halfs ('f-half' and 'l-half')
   and put them at the start and end of the full length:


        | f-half  | tile size | tile size | tile size | l-half  |
        |         |           |           |           |         |
        ---------------------------------------------------------

   So in any case, knowing the size of the first tile, will allow us to
   parse all the tiles. We just have to make sure we don't go over the full
   input's length. */
static void
gal_tile_full_regular_first(gal_data_t *parent, size_t *regular,
                            float remainderfrac, size_t *first, size_t *last,
                            size_t *tsize)
{
  size_t i, remainder, *dsize=parent->dsize;;

  /* For each dimension, set the size of the first tile. */
  for(i=0;i<parent->ndim;++i)
    {
      /* It might happen that the tile size is bigger than the parent size
         in a dimension, in that case the analysis in the comments above
         are useless and only one tile should cover this dimension with the
         size of the parent. */
      if( regular[i] >= dsize[i] )
        {
          tsize[i]=1;
          first[i]=last[i]=dsize[i];
        }
      else
        {
          /* Calculate the remainder in this dimension. */
          remainder=dsize[i] % regular[i];

          /* Depending on the remainder, set the first tile size and
             number. */
          if(remainder)
            {
              if( remainder > remainderfrac * regular[i] )
                {
                  first[i]  = ( remainder + regular[i] )/2;
                  tsize[i]  = dsize[i]/regular[i] + 1 ;

                  /* If we only have one tile along the dimension, then
                     'first[i]==dsize[i]'. In this case, the first and last
                     tiles are the same and must have the same size. */
                  last[i]   = ( first[i]==dsize[i]
                                ? first[i]
                                : ( dsize[i]
                                    - ( first[i] + regular[i]*(tsize[i]-2) ) ) );
                }
              else
                {
                  first[i]  = remainder + regular[i];
                  tsize[i]  = dsize[i]/regular[i];
                  last[i]   = first[i]==dsize[i] ? first[i] : regular[i];
                }
            }
          else
            {
              first[i]  = last[i] = regular[i];
              tsize[i] = dsize[i]/regular[i];
            }
        }
    }

  /* For a check:
  printf("%s: first: %zu, %zu\n", __func__, first[0], first[1]);
  printf("%s: last: %zu, %zu\n",  __func__, last[0],  last[1]);
  */
}





/* Cover the full dataset with (mostly) identical tiles. The regular tile
   size is determined from the 'size' array. If the input data's size is
   not an exact multiple of 'size' for each dimension, then the tiles
   touching the edges in that dimension will have a different size to fully
   cover every element of the input. For a full description of tiling in
   'gal_data_t', please see 'data.h'.

   Inputs
   ------

     'input' is the gal_data_t which you want to tile (only used for its
        sizes).

     'regular' is the size of the regular tiles along each of the input's
        dimensions. So it must have the same number of elements as the
        dimensions of 'input'.

     'remainderfrac' is the significant fraction of the remainder space if
        the width of the input isn't an exact multiple of the tile size
        along a dimension, see 'gal_tile_full_regular_first'.

     'out' is the pointer to the array of data structures that is to keep
        the tile parameters. If '*out==NULL', then the necessary space will
        be allocated. If it is not NULL, then all the tile information will
        be filled from the given element, see 'multiple' for more.

     'multiple': When the '*out' array is to be allocated, allocate
        'multiple' times the necessary space. This can be very useful when
        you have several more identically sized 'inputs', and you want all
        their tiles to be allocated (and thus indexed) together, even
        though they have different 'block' datasets (that then link to one
        allocated space).  See the 'gal_tile_full_two_layers' below.

     'firsttsize': The size of the first tile along every dimension. This
        is only different from the regular tile size when 'regular' is not
        an exact multiple of 'input''s length along every dimension. This
        array is allocated internally by this function.

   Output
   ------

     The returned output is an array of numbers (the same size as the input
        data structure's dimensions) keeping the number of tiles along each
        dimension.


   Implementation
   --------------

     In the most general case, to set the starting pointers for each tile
     we need the following sizes. If the input array has no parent/block,
     then both these sizes are equal to it's own size:

        1. block-size (or 'bsize'), which is the size of the allocated
           block in each dimension.

        2. parent-size (or 'psize') which is the size of the parent in each
           dimension (we don't want to go out of the paren't range). */
size_t *
gal_tile_full(gal_data_t *input, size_t *regular,
              float remainderfrac, gal_data_t **out, size_t multiple,
              size_t **firsttsize)
{
  size_t i, d, tind, numtiles, *start=NULL;
  gal_data_t *tiles, *block=gal_tile_block(input);
  size_t *last   = gal_pointer_allocate(GAL_TYPE_SIZE_T, input->ndim, 0,
                                      __func__, "last");
  size_t *first  = gal_pointer_allocate(GAL_TYPE_SIZE_T, input->ndim, 0,
                                      __func__, "first");
  size_t *coord  = gal_pointer_allocate(GAL_TYPE_SIZE_T, input->ndim, 0,
                                      __func__, "coord");
  size_t *tcoord = gal_pointer_allocate(GAL_TYPE_SIZE_T, input->ndim, 0,
                                      __func__, "tcoord");
  size_t *tsize  = gal_pointer_allocate(GAL_TYPE_SIZE_T, input->ndim+1, 0,
                                      __func__, "tsize");


  /* Set the first tile size and total number of tiles along each
     dimension, then allocate the array of tiles. */
  gal_tile_full_regular_first(input, regular, remainderfrac,
                              first, last, tsize);
  numtiles=gal_dimension_total_size(input->ndim, tsize);


  /* Allocate the necessary space for all the tiles (if necessary). */
  if(*out)        tiles = *out;
  else     *out = tiles = gal_data_array_calloc(numtiles*multiple);


  /* It is possible that the 'input' dataset is its-self a larger tile over
     a region of the allocated block. In that case, we need to account for
     the block's dimensions when calculating the position of this block. */
  if(input->block)
    {
      start=gal_pointer_allocate(GAL_TYPE_SIZE_T, input->ndim, 0, __func__,
                                 "start");
      gal_tile_start_coord(input, start);
    }


  /* Initialize each tile. */
  for(i=0;i<numtiles;++i)
    {
      /* Specify the coordinates of the tile between the other tiles. Note
         that we are dealing with tiles here, not pixels. */
      gal_dimension_index_to_coord(i, input->ndim, tsize, tcoord);

      /* The coordinates are currently in units of tiles, not
         pixels. Convert them to the coordinates of the first pixel in each
         tile. */
      for(d=0;d<input->ndim;++d)
        {
          /* Convert the tile coordinates to pixel coordinates within
             'input'. See the comments above 'gal_tile_full_regular_first':
             The first tile in every dimension can be different from the
             regular tile size. */
          coord[d] = tcoord[d] ? first[d] + (tcoord[d]-1)*regular[d] : 0;

          /* When the 'input' data structure (that is to be tiled here) was
             itself a tile over a larger allocated array, a 'start' array
             has been allocated to correct the coordinates so they refer to
             a physical position on the allocated block of memory. */
          if(start)
            coord[d] += start[d];
        }

      /* Convert the coordinates (that are now in element/pixel units on
         the allocated block of memory) into an index. */
      tind=gal_dimension_coord_to_index(block->ndim, block->dsize, coord);

      /* Now that we have the index of this tile's starting point compared
         to the allocated block, put it in to the tile's 'array'
         pointer. */
      tiles[i].array=gal_pointer_increment(block->array, tind, block->type);

      /* Set the sizes of the tile. */
      tiles[i].size=1; /* Just an initializer, will be changed. */
      tiles[i].ndim=input->ndim;
      tiles[i].minmapsize=input->minmapsize;
      tiles[i].dsize=gal_pointer_allocate(GAL_TYPE_SIZE_T,input->ndim, 0,
                                          __func__, "tiles[i].dsize");
      for(d=0;d<input->ndim;++d)
        {
          /* The size of the first and last tiles can be different from the
             majority of the 'regular' tiles that have the same size. When
             a tile is on the edge in one of the dimensions, then its
             'tcoord[d]' will be either 0 or the last. */
          if( first[d] != regular[d]
              && ( tcoord[d]==0 || tcoord[d]==tsize[d]-1 ) )
            {
              if( tcoord[d] == 0          ) tiles[i].dsize[d] = first[d];
              if( tcoord[d] == tsize[d]-1 ) tiles[i].dsize[d] = last[d];
            }
          else
            tiles[i].dsize[d]=regular[d];

          /* Set the size value. */
          tiles[i].size *= tiles[i].dsize[d];
        }

      /* Set the block structure for this tile to the 'input', and set the
         next pointer as the next tile. Note that only when we are dealing
         with the last tile should the 'next' pointer be set to NULL.*/
      tiles[i].flag  = 0;
      tiles[i].block = input;
      tiles[i].next  = i==numtiles-1 ? NULL : &tiles[i+1];

      /* For a check:
      printf("%zu:\n\tStart index: %zu\n\tsize: %zu x %zu\n", i, tind,
             tiles[i].dsize[1], tiles[i].dsize[0]);
      exit(0);
      */
    }


  /* Clean up and return. */
  free(last);
  free(coord);
  free(tcoord);
  *firsttsize=first;
  if(start) free(start);
  tsize[input->ndim]=-1; /* 'tsize' had ndim+1 values, we will mark the  */
  return tsize;          /* extra space with the largest possible value: */
}                        /* -1, see 'gal_tile_full_sanity_check'.        */





/* Make sure that the input parameters (in 'tl', short for two-layer) fit
   with the input dataset. The filename and HDU are only required for error
   messages. Also, allocate and fill the 'channelsize' array. */
void
gal_tile_full_sanity_check(char *filename, char *hdu, gal_data_t *input,
                           struct gal_tile_two_layer_params *tl)
{
  double d;
  size_t i, ndim=input->ndim;

  /* Check the tile's dimensions. */
  for(i=0;tl->tilesize[i]!=-1;++i)
    {
      /* Not equal to zero. */
      if(tl->tilesize[i]==0)
        error(EXIT_FAILURE, 0, "'--tilesize' must be larger than zero, "
              "the given value for dimension %zu was zero", ndim-i);

      /* If the tile size is larger than the dataset size in this
         dimension, then quietly change the tile size to the dataset size
         along that dimension. */
      if( tl->tilesize[i] > input->dsize[i] )
        tl->tilesize[i] = input->dsize[i];
    }


  /* Make sure the number of tile sizes (tile dimensions) are the same as
     the dataset's dimensions). */
  if(i!=ndim)
    error(EXIT_FAILURE, 0, "%s (hdu: %s): has %zu dimensions, but only %zu "
          "value(s) given for the tile size ('--tilesize' option).",
          filename, hdu, ndim, i);


  /* Check the channel's dimensions. */
  for(i=0; tl->numchannels[i]!=-1; ++i)
    if(tl->numchannels[i]==0)
      error(EXIT_FAILURE, 0, "the number of channels in all dimensions must "
            "be larger than zero. The number for dimension %zu was zero",
            i+1);
  if(i!=ndim)
    error(EXIT_FAILURE, 0, "%s (hdu: %s): has %zu dimensions, but only %zu "
          "value(s) given for the number of channels", filename, hdu, ndim,
          i);


  /* Allocate space for the channel sizes. */
  tl->channelsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                       "tl->channelsize");


  /* Check if the channels are exactly divisible by the input's size along
     each dimension and set the correct size. */
  for(i=0;i<ndim;++i)
    {
      /* Check if the number of channels is not more than the size of the
         image. Note that the reported dimension must be in FITS format.*/
      if( input->dsize[i] < tl->numchannels[i] )
        error(EXIT_FAILURE, 0, "the number of channels in dimension %zu "
              "(%zu) is more than the size of the '%s' (hdu: %s) in that "
              "dimension", ndim-i, tl->numchannels[i], filename, hdu);

      /* Also check the tile size. */
      if( input->dsize[i] < tl->tilesize[i] )
        error(EXIT_FAILURE, 0, "the tile size in dimension %zu (%zu) is "
              "more than the size of the '%s' (hdu: %su) in that dimension",
              ndim-i, tl->tilesize[i], filename, hdu);

      /* First check. */
      d=(double)input->dsize[i]/(double)(tl->numchannels[i]);
      if(ceil(d)!=d)
        error(EXIT_FAILURE, 0, "%zu (number of channels along dimension "
              "%zu) is not exactly divisible by %zu (the length of '%s' "
              "(hdu: %s) that dimension). The channels cover the input "
              "dataset, hence, they must be identical", tl->numchannels[i],
              ndim-i, input->dsize[i], filename, hdu);

      /* Put the channel size into the output. */
      tl->channelsize[i]=d;
    }
}





/* A dataset can be tiled with two layers that are related:

      Channels: A tesselation of larger tile sizes that all have the same
           size ('channel_size' must be an exact multiple of 'input's size
           along every dimension. In astronomy images, this can be seen as
           CCD amplifiers, that cover large parts of the image. If
           '*channels!=NULL' then it is assumed to be already present and
           will not be allocated.

      Tiles: A combined tesselation of each channel with smaller
           tiles. These tiles can be used to calculate things like
           gradients over each channel and thus over the whole image. */
void
gal_tile_full_two_layers(gal_data_t *input,
                         struct gal_tile_two_layer_params *tl)
{
  gal_data_t *t;
  size_t i, *junk, *junk2, ndim=tl->ndim=input->ndim;

  /* Initialize.  */
  tl->channels=tl->tiles=NULL;

  /* Initialize necessary values and do the channels tessellation. */
  junk = gal_tile_full(input, tl->channelsize, tl->remainderfrac,
                       &tl->channels, 1, &junk2);
  tl->totchannels = gal_dimension_total_size(ndim, tl->numchannels);
  for(i=0;i<ndim;++i)
    if(junk[i]!=tl->numchannels[i])
      error(EXIT_FAILURE, 0, "%s: the input and output number of channels "
            "don't match in dimension %zu: %zu and %zu respectively.",
            __func__, ndim-i, tl->numchannels[i], junk[i]);
  free(junk);
  free(junk2);

  /* Tile each channel. While tiling the first channel, we are also going
     to allocate the space for the other channels. Then pass those pointers
     when we want to fill in each tile of the other channels. */
  tl->numtilesinch = gal_tile_full(tl->channels, tl->tilesize,
                                   tl->remainderfrac, &tl->tiles,
                                   tl->totchannels, &tl->firsttsize);
  tl->tottilesinch = gal_dimension_total_size(ndim, tl->numtilesinch);
  for(i=1; i<tl->totchannels; ++i)
    {
      /* Set the first tile in this channel. Then use it it fill the 'next'
         pointer of the previous channel's tiles. Note that 'gal_tile_full'
         set this 'next' element to NULL. */
      t = tl->tiles + i * tl->tottilesinch;
      tl->tiles[ i * tl->tottilesinch - 1 ].next = t;

      /* Fill in the information for all the tiles in this channel. Note
         that we already have the returned value, so it isn't important.*/
      junk=gal_tile_full(&tl->channels[i], tl->tilesize, tl->remainderfrac,
                         &t, 1, &junk2);
      free(junk);
      free(junk2);
    }

  /* Multiply the number of tiles along each dimension OF ONE CHANNEL by
     the number of channels in each dimension to get the dimensionality of
     the full tile structure. */
  tl->numtiles = gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                      "tl->numtiles");
  for(i=0;i<ndim;++i)
    tl->numtiles[i] = tl->numtilesinch[i] * tl->numchannels[i];
  tl->tottiles = gal_dimension_total_size(ndim, tl->numtiles);
}





/* Usage
   -----

   Make a permutation to allow the conversion of tile location in memory to
   its location in the full input dataset and put it in the input's
   'permutation' element. If a permutation has already been defined for the
   tessellation, this function will not do anythin. If permutation won't be
   necessary, then this function will just return (the permutation must
   have been initialized to NULL). */
void
gal_tile_full_permutation(struct gal_tile_two_layer_params *tl)
{
  size_t *ch_coord, *tinch_coord;
  size_t i, p=0, t, ch, ind_in_all, ndim=tl->ndim;

  /* If the permutation has already been defined for this tessellation,
     then there is no need to do it again here. */
  if(tl->permutation) return;

  /* If there is only one channel or one dimension, return NULL. The
     permutation functions know that the input and output indexs are the
     same when the permutation is NULL. */
  if( ndim==1 || tl->totchannels==1) return;

  /* Allocate the space for the permutation and coordinates. */
  ch_coord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                "ch_coord");
  tinch_coord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                 "tinch_coord");
  tl->permutation=gal_pointer_allocate(GAL_TYPE_SIZE_T, tl->tottiles, 0,
                                       __func__, "tl->permutation");

  /* Fill in the permutation, we use the fact that the tiles are filled
     from the first channel to the last. */
  for(ch=0;ch<tl->totchannels;++ch)
    {
      /* Get the coordinates of this channel's first tile. */
      gal_dimension_index_to_coord(ch, ndim, tl->numchannels, ch_coord);
      for(i=0;i<ndim;++i) ch_coord[i] *= tl->numtilesinch[i];

      /* Go over all the tiles in this channel. */
      for(t=0;t<tl->tottilesinch;++t)
        {
          /* Convert its index to coordinates and add them to the channel's
             starting coordinates. */
          gal_dimension_index_to_coord(t, ndim, tl->numtilesinch,
                                       tinch_coord);
          for(i=0;i<ndim;++i) tinch_coord[i] += ch_coord[i];

          /* Convert the coordinates into an index. */
          ind_in_all = gal_dimension_coord_to_index(ndim, tl->numtiles,
                                                    tinch_coord);
          tl->permutation[ind_in_all] = p++;
        }
    }

  /* Clean up and return. */
  free(tinch_coord);
  free(ch_coord);
}





/* Write one value for each tile into a file.

   IMPORTANT: it is assumed that the values are in the same order as the
   tiles.

                      tile[i]  -->   tilevalues[i]                       */
void
gal_tile_full_values_write(gal_data_t *tilevalues,
                           struct gal_tile_two_layer_params *tl,
                           int withblank, char *filename,
                           gal_fits_list_key_t *keys, char *program_string)
{
  gal_data_t *disp;

  /* Make the dataset to be displayed. */
  if(tl->oneelempertile)
    {
      if(tl->ndim>1 && tl->totchannels>1)
        {
          /* A small sanity check. */
          if(tl->permutation==NULL)
            error(EXIT_FAILURE, 0, "%s: no permutation defined for the input "
                  "tessellation", __func__);

          /* Writing tile values to disk is not done for checking, not for
             efficiency. So to be safe (allow the caller to work on
             multiple threads), we will copy the tile values, then permute
             those. */
          disp = gal_data_copy(tilevalues);
          gal_permutation_apply(disp, tl->permutation);
        }
      else disp = tilevalues;
    }
  else
    disp=gal_tile_block_write_const_value(tilevalues, tl->tiles,
                                          withblank, 0);

  /* Write the array as a file and then clean up (if necessary). */
  gal_fits_img_write(disp, filename, keys, program_string);
  if(disp!=tilevalues) gal_data_free(disp);
}





/* Smooth the given values with a flat kernel of the given width. */
gal_data_t *
gal_tile_full_values_smooth(gal_data_t *tilevalues,
                            struct gal_tile_two_layer_params *tl,
                            size_t width, size_t numthreads)
{
  size_t *kdsize, knum, i;
  gal_data_t *kernel, *smoothed;
  struct gal_tile_two_layer_params ttl={0};
  int permute=tl->ndim>1 && tl->totchannels>1;


  /* Check if the width is odd. */
  if(width%2==0)
    error(EXIT_FAILURE, 0, "%s: %zu not acceptable as width. It has to be "
          "an odd number", __func__, width);


  /* Prepare the kernel size along every dimension. */
  kdsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, tl->ndim, 0, __func__,
                              "kdsize");
  for(i=0;i<tl->ndim;++i) kdsize[i]=width;


  /* Make the kernel. */
  kernel=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, tilevalues->ndim,
                        kdsize, NULL, 0, -1, 1, NULL, NULL, NULL);
  knum=gal_dimension_total_size(tl->ndim, kernel->dsize);
  for(i=0;i<knum;++i) ((float *)(kernel->array))[i]=1/((double)knum);

  /* Permute (if necessary). */
  if(permute)
    {
      gal_tile_full_permutation(tl);
      gal_permutation_apply(tilevalues, tl->permutation);
    }

  /* Do the smoothing. */
  if(tl->workoverch)
    smoothed=gal_convolve_spatial(tilevalues, kernel, numthreads, 1, 1);
  else
    {
      /* Create the tile structure. */
      ttl.tilesize=tl->numtilesinch;
      ttl.numchannels=tl->numchannels;
      gal_tile_full_sanity_check("IMPOSSIBLE", "IMP_HDU", tilevalues, &ttl);
      gal_tile_full_two_layers(tilevalues, &ttl);

      /* Do the convolution separately on each channel. */
      smoothed=gal_convolve_spatial(ttl.tiles, kernel, numthreads, 1, 0);

      /* Clean up. */
      ttl.tilesize=ttl.numchannels=NULL;
      gal_tile_full_free_contents(&ttl);
    }

  /* Reverse the permutation. */
  if(permute) gal_permutation_apply_inverse(smoothed, tl->permutation);

  /* Clean up and return; */
  free(kdsize);
  gal_data_free(kernel);
  return smoothed;
}





size_t
gal_tile_full_id_from_coord(struct gal_tile_two_layer_params *tl,
                            size_t *coord)
{
  /* This function only works for 10 dimensions. */
  size_t i, tr, chid, tile[10];


  /* Host channel's ID. */
  for(i=0;i<tl->ndim;++i)
    tile[i] = tl->totchannels == 1 ? 0 : coord[i] / tl->channelsize[i];
  chid=gal_dimension_coord_to_index(tl->ndim, tl->numchannels, tile);


  /* Find the tile within the channel. */
  for(i=0;i<tl->ndim;++i)
    {
      tr=coord[i] % tl->channelsize[i];
      if( tl->firsttsize[i] != tl->tilesize[i] )
        tile[i] = ( tr <= tl->firsttsize[i]
                    ? 0
                    : 1 + (tr - tl->firsttsize[i]) / tl->tilesize[i] );
      else
        tile[i] = tr / tl->tilesize[i];
    }


  /* Return the tile ID. */
  return ( chid * tl->tottilesinch
           + gal_dimension_coord_to_index(tl->ndim, tl->numtilesinch, tile) );
}






/* Clean up the allocated spaces in the parameters. */
void
gal_tile_full_free_contents(struct gal_tile_two_layer_params *tl)
{
  /* Free the simply allocated spaces. */
  if(tl->tilesize)      free(tl->tilesize);
  if(tl->numchannels)   free(tl->numchannels);
  if(tl->channelsize)   free(tl->channelsize);
  if(tl->numtiles)      free(tl->numtiles);
  if(tl->numtilesinch)  free(tl->numtilesinch);
  if(tl->tilecheckname) free(tl->tilecheckname);
  if(tl->permutation)   free(tl->permutation);
  if(tl->firsttsize)    free(tl->firsttsize);

  /* Free the arrays of 'gal_data_t' for each tile and channel. */
  if(tl->tiles)    gal_data_array_free(tl->tiles,    tl->tottiles,    0);
  if(tl->channels) gal_data_array_free(tl->channels, tl->totchannels, 0);
}
