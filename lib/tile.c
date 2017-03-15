/*********************************************************************
tile -- work with tesselations over a host dataset.
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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/fits.h>
#include <gnuastro/tile.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/convolve.h>
#include <gnuastro/multidim.h>

#include "checkset.h"









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
    memset(start_coord, 0, ndim*gal_data_sizeof(GAL_DATA_TYPE_SIZE_T));
  else
    {
      /* Calculate the coordinates of the first pixel of the tile. */
      ind = gal_data_ptr_dist(block->array, tile->array, block->type);
      gal_multidim_index_to_coord(ind, ndim, block->dsize, start_coord);
    }
}





/* Put the starting and ending (end point is not inclusive) coordinates of
   a tile into the `start_end' array. It is assumed that a space of
   `2*tile->ndim' has been already allocated (static or dynamic) before
   this function is called.

   `rel_block' (or relative-to-block) is only relevant when the tile has an
   intermediate tile between it and the allocated space (like a channel,
   see `gal_tile_full_two_layers'). If it doesn't (`tile->block' points the
   allocated dataset), then the value to `rel_block' is irrelevant.

   However, when `tile->block' is its self a larger block and `rel_block'
   is set to 0, then the starting and ending positions will be based on the
   position within `tile->block', not the allocated space. */
void
gal_tile_start_end_coord(gal_data_t *tile, size_t *start_end, int rel_block)
{
  size_t *s, *sf, *h;
  gal_data_t *block=gal_tile_block(tile);
  gal_data_t *host=rel_block ? block : tile->block;
  size_t *hcoord, start_ind, ndim=tile->ndim, *end=start_end+ndim;

  /* Get the starting index. Note that for the type we need the allocated
     block dataset and can't rely on the tiles. */
  start_ind=gal_data_ptr_dist(block->array, tile->array, block->type);

  /* Get the coordinates of the starting point relative to the allocated
     block. */
  gal_multidim_index_to_coord(start_ind, ndim, block->dsize, start_end);

  /* When the host is different from the block, the tile's starting
     position needs to be corrected. */
  if(host!=block)
    {
      /* Get the host's starting coordinates. */
      start_ind=gal_data_ptr_dist(block->array, host->array, block->type);

      /* Temporarily put the host's coordinates in the place held for the
         ending coordinates. */
      hcoord=end;
      gal_multidim_index_to_coord(start_ind, ndim, block->dsize, hcoord);
      sf=(s=start_end)+ndim; h=hcoord; do *s++ -= *h++; while(s<sf);
    }

  /* Add the dimensions of the tile to the starting coordinate. Note that
     the ending coordinates are stored immediately after the start.*/
  gal_multidim_add_coords(start_end, tile->dsize, end, ndim);
}





/* Put the indexs of the first/start and last/end pixels (inclusive) in a
   tile into the `start_end' array (that has two elements). It will then
   return the pointer to the start of the tile in the `work' data
   structure. */
void *
gal_tile_start_end_ind_inclusive(gal_data_t *tile, gal_data_t *work,
                                 size_t *start_end_inc)
{
  gal_data_t *block=gal_tile_block(tile);
  size_t ndim=tile->ndim, *s, *e, *l, *sf;
  size_t *start_coord = gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, ndim);
  size_t *end_coord   = gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, ndim);


  /* The starting index can be found from the distance of the `tile->array'
     pointer and `block->array' pointer. IMPORTANT: with the type of the
     block array.  */
  start_end_inc[0]=gal_data_ptr_dist(block->array, tile->array, block->type);


  /* To find the end index, we need to know the coordinates of the starting
     point in the allocated block.  */
  gal_multidim_index_to_coord(start_end_inc[0], ndim, block->dsize,
                              start_coord);


  /* `end_coord' is one unit ahead of the last element in the tile in every
     dimension. To have less potential for bugs, we will remove that extra
     value, so we get the coordinates of the last pixel in the tile
     (inclusive). We will finally, increment that value by one to get to
     the pixel immediately outside of the tile.*/
  e=end_coord;
  l=tile->dsize;
  sf=(s=start_coord)+ndim; do *e++ = *s + *l++ - 1; while(++s<sf);


  /* Convert the (inclusive) ending point's coordinates into an index. */
  start_end_inc[1]=gal_multidim_coord_to_index(ndim, block->dsize, end_coord);


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
  return gal_data_ptr_increment(work->array, start_end_inc[0], work->type);
}




















/***********************************************************************/
/**************        Allocated block of memory      ******************/
/***********************************************************************/
/* When you are working on an array, it important to know the size of the
   allocated space in each dimension. This simple function will just follow
   the block pointer and return the `dsize' element of lowest-level
   structure. */
gal_data_t *
gal_tile_block(gal_data_t *input)
{
  while(input->block!=NULL) input=input->block;
  return input;
}





/* Return the increment necessary to start at the next series of contiguous
   memory (fastest dimension) associated with a tile. See
   `gal_tile_block_check_tiles' as one example application of this
   function.

   1D and 2D cases are simple and need no extra explanation, but the case
   for higher dimensions can be alittle more complicated, So we will go
   over some examples. The notations below are:

       `n'     number of dimensions (same in tile and block).
       `t[]'   size of the tile in each dimension.
       `b[]'   size of the allocated block in each dimension.

   It is just important to see the output of this function as an increment
   from the the last patch of contiguous memory associated with the
   tile. So when the increment number is `t[n-1]' (the first 2D slice of
   the tile has been parsed), simply incrementing by `b[n-2] * b[n-1]' will
   take us to the last row of

  num_increment      coord         increment
  -------------      -----         ---------
         0          (...0,0,0)     b[n-1]: fastest dimension of the block.
         1          (...0,1,0)     Similar to previous
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
  size_t increment;
  size_t n=block->ndim;
  size_t *b=block->dsize, *t=tsize;

  if(n>3)
    error(EXIT_FAILURE, 0, "`gal_tile_block-increment' is currently only "
          "implemented for at most 3 dimensions");

  switch(n)
    {
    case 0:
      error(EXIT_FAILURE, 0, "zero dimensional input is not acceptable in "
            "`gal_tile_block_parse'");

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

    /* Higher dimensions. */
    default:
      if(num_increment % t[n-2])
        {
          increment=b[n-1];
          if(coord) ++coord[n-2];
        }
      else
        {
          increment=(b[n-2] * b[n-1]) - ( (t[n-2]-1) * b[n-1] );
          ++coord[n-3];
          if(coord) coord[n-2]=coord[n-1]=0;
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

     `tilevalues': This must be an array that has the same number of
        elements as that in `tilesll' and in the same order that `tilesll'
        elements are parsed (from first to last). As a result the
        dimensionality of this array is irrelevant. Note that unlike
        `tiles', `tilevalues' must be an array.

     `tiles': This will be parsed as a linked list (using the `next'
        element). Internally, it might be stored as an array, but this
        function doesn't care! The position of the tile over its block will
        be determined according to the `block' element and the pointer of
        its `array' as fully described in `gnuastro/data.h'. This function
        will not pop/free the list, it will only parse it from start to
        end.

     `initialize': Initialize the allocated space with blank values before
        writing in the constant values. This can be useful when the tiles
        don't cover the full allocated block. */
gal_data_t *
gal_tile_block_write_const_value(gal_data_t *tilevalues, gal_data_t *tilesll,
                                 int initialize)
{
  void *in;
  int type=tilevalues->type;
  size_t tile_ind, nt=0, nv=tilevalues->size;
  gal_data_t *tofill, *tile, *block=gal_tile_block(tilesll);

  /* A small sanity check. */
  for(tile=tilesll; tile!=NULL; tile=tile->next) ++nt;
  if(nt!=nv)
    error(EXIT_FAILURE, 0, "the number of elements in `tilevalues' (%zu) "
          "and `tilesll' (%zu) must be the same in "
          "`gal_tile_block_write_const_value'", nv, nt);

  /* Allocate the output array. */
  tofill=gal_data_alloc(NULL, type, block->ndim, block->dsize, block->wcs,
                        0, block->minmapsize, tilevalues->name,
                        tilevalues->unit, tilevalues->comment);

  /* If requested, initialize `tofill'. */
  if(initialize) gal_blank_initialize(tofill);

  /* Go over the tiles and write the values in. Note Recall that `tofill'
     has the same type as `tilevalues'. So we are using memcopy. */
  tile_ind=0;
  for(tile=tilesll; tile!=NULL; tile=tile->next)
    {
      /* Set the pointer to use as input. */
      in=gal_data_ptr_increment(tilevalues->array, tile_ind++, type);;
      GAL_TILE_PARSE_OPERATE({memcpy(o, in, gal_data_sizeof(type));},
                             tile, tofill, 1, 0);
    }

  return tofill;
}





/* Make a copy of the memory block in integer type and fill it with the ID
   of each tile, the non-filled areas have blank values. Finally, save the
   final array into a FITS file, specified with `filename'. This is done
   mainly for inspecting the positioning of tiles. We are using a signed
   32-bit type because this is the standard FITS standard type for
   integers. */
gal_data_t *
gal_tile_block_check_tiles(gal_data_t *tilesll)
{
  int32_t *arr;
  size_t i, dsize=gal_data_num_in_ll(tilesll);
  gal_data_t *ids, *out, *block=gal_tile_block(tilesll);

  /* Allocate the array to keep the IDs of each tile. */
  ids=gal_data_alloc(NULL, GAL_DATA_TYPE_INT32, 1, &dsize,
                     NULL, 0, block->minmapsize, NULL, NULL, NULL);

  /* Put the IDs into the array. */
  arr=ids->array; for(i=0;i<dsize;++i) arr[i]=i;

  /* Make the output. */
  out=gal_tile_block_write_const_value(ids, tilesll, 1);

  /* Clean up and return. */
  gal_data_free(ids);
  return out;
}



















/***********************************************************************/
/**************           Tile full dataset         ********************/
/***********************************************************************/
/* The user's specified tile size might not be an exact multiple of the
   parent's size. This function is useful in such cases. It will give the
   starting tile's size along each dimension. The line below can be the
   length along any dimension and the tile size along that dimension. You
   see that when we start with the tile size, we will end up with a last
   tile that contains the remainder elements.

        | tile size | tile size | tile size | tile size | remainder
        |           |           |           |           |       |
        ---------------------------------------------------------

   The remainder will always be smaller than `tile size'. So, we will merge
   the last tile size with the remainder and move that tile to the start.
   In this way, the size of the first tile will always be between between
   one and two times the size of the regular tile:

        | first tile        | tile size | tile size | tile size |
        |                   |           |           |           |
        ---------------------------------------------------------

   When there is only a small remainder (for example one or two pixels),
   then this layout is fine. But when the remainder is significant compared
   to the regular tile size (like the example above), then it will make
   more sense to cut the first tile into two halfs (`f-half' and `l-half')
   and put them at the start and end of the full length:


        | f-half  | tile size | tile size | tile size | l-half  |
        |         |           |           |           |         |
        ---------------------------------------------------------

   So in any case, knowing the size of the first tile, will allow us to
   parse all the tiles. We just have to make sure we don't go over the full
   input's length. */
static void
gal_tile_full_regular_first(gal_data_t *parent, size_t *regular,
                            float significance, size_t *first, size_t *last,
                            size_t *tsize)
{
  size_t i, remainder, *dsize=parent->dsize;;

  /* For each dimension, set the size of the first tile. */
  for(i=0;i<parent->ndim;++i)
    {
      /* Calculate the remainder in this dimension. */
      remainder=dsize[i] % regular[i];

      /* Depending on the remainder, set the first tile size and number. */
      if(remainder)
        {
          if( remainder > significance * regular[i] )
            {
              first[i]  = ( remainder + regular[i] )/2;
              tsize[i] = dsize[i]/regular[i] + 1 ;
              last[i]   = dsize[i] - ( first[i] + regular[i]*(tsize[i]-2) );
            }
          else
            {
              first[i]  = remainder + regular[i];
              tsize[i]  = dsize[i]/regular[i];
              last[i]   = regular[i];
            }
        }
      else
        {
          first[i]  = last[i] = regular[i];
          tsize[i] = dsize[i]/regular[i];
        }
    }
}





/* Cover the full dataset with (mostly) identical tiles. The regular tile
   size is determined from the `size' array. If the input data's size is
   not an exact multiple of `size' for each dimension, then the tiles
   touching the edges in that dimension will have a different size to fully
   cover every element of the input. For a full description of tiling in
   `gal_data_t', please see `data.h'.

   Inputs
   ------

     `input' is the gal_data_t which you want to tile (only used for its
        sizes).

     `regular' is the size of the regular tiles along each of the input's
        dimensions. So it must have the same number of elements as the
        dimensions of `input'.

     `remainderfrac' is the significant fraction of the remainder space if
        the width of the input isn't an exact multiple of the tile size
        along a dimension, see `gal_tile_full_regular_first'.

     `out' is the pointer to the array of data structures that is to keep
        the tile parameters. If `*out==NULL', then the necessary space will
        be allocated. If it is not NULL, then all the tile information will
        be filled from the given element, see `multiple' for more.

     `multiple': When the `*out' array is to be allocated, allocate
        `multiple' times the necessary space. This can be very useful when
        you have several more identically sized `inputs', and you want all
        their tiles to be allocated (and thus indexed) together, even
        though they have different `block' datasets (that then link to one
        allocated space).  See the `gal_tile_full_two_layers' below.

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

        1. block-size (or `bsize'), which is the size of the allocated
           block in each dimension.

        2. parent-size (or `psize') which is the size of the parent in each
           dimension (we don't want to go out of the paren't range).
*/
size_t *
gal_tile_full(gal_data_t *input, size_t *regular,
              float remainderfrac, gal_data_t **out, size_t multiple)
{
  size_t i, d, tind, numtiles, *start=NULL;
  gal_data_t *tiles, *block=gal_tile_block(input);
  size_t *last   = gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, input->ndim);
  size_t *tsize  = gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, input->ndim);
  size_t *first  = gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, input->ndim);
  size_t *coord  = gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, input->ndim);
  size_t *tcoord = gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, input->ndim);


  /* Set the first tile size and total number of tiles along each
     dimension, then allocate the array of tiles. */
  gal_tile_full_regular_first(input, regular, remainderfrac,
                             first, last, tsize);
  numtiles=gal_multidim_total_size(input->ndim, tsize);


  /* Allocate the necessary space for all the tiles (if necessary). */
  if(*out)        tiles = *out;
  else     *out = tiles = gal_data_array_calloc(numtiles*multiple);


  /* It is possible that the `input' dataset is its-self a larger tile over
     a region of the allocated block. In that case, we need to account for
     the block's dimensions when calculating the position of this block. */
  if(input->block)
    {
      start=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, input->ndim);
      gal_tile_start_coord(input, start);
    }


  /* Initialize each tile. */
  for(i=0;i<numtiles;++i)
    {
      /* Specify the coordinates of the tile between the other tiles. Note
         that we are dealing with tiles here, not pixels. */
      gal_multidim_index_to_coord(i, input->ndim, tsize, tcoord);

      /* The coordinates are currently in units of tiles, not
         pixels. Convert them to the coordinates of the first pixel in each
         tile. */
      for(d=0;d<input->ndim;++d)
        {
          /* Convert the tile coordinates to pixel coordinates within
             `input'. See the comments above `gal_tile_full_regular_first':
             The first tile in every dimension can be different from the
             regular tile size. */
          coord[d] = tcoord[d] ? first[d] + (tcoord[d]-1)*regular[d] : 0;

          /* When the `input' data structure (that is to be tiled here) was
             itself a tile over a larger allocated array, a `start' array
             has been allocated to correct the coordinates so they refer to
             a physical position on the allocated block of memory. */
          if(start)
            coord[d] += start[d];
        }

      /* Convert the coordinates (that are now in element/pixel units on
         the allocated block of memory) into an index. */
      tind=gal_multidim_coord_to_index(block->ndim, block->dsize, coord);

      /* Now that we have the index of this tile's starting point compared
         to the allocated block, put it in to the tile's `array'
         pointer. */
      tiles[i].array=gal_data_ptr_increment(block->array, tind, block->type);

      /* Set the sizes of the tile. */
      tiles[i].size=1;
      tiles[i].ndim=input->ndim;
      tiles[i].minmapsize=input->minmapsize;
      tiles[i].dsize=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T,input->ndim);
      for(d=0;d<input->ndim;++d)
        {
          /* The size of the first and last tiles can be different from the
             majority of the `regular' tiles that have the same size. When
             a tile is on the edge in one of the dimensions, then its
             `coord[d]' will be either 0 or the last. */
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

      /* Set the block structure for this tile to the `input', and set the
         next pointer as the next tile. Note that only when we are dealing
         with the last tile should the `next' pointer be set to NULL.*/
      tiles[i].block = input;
      tiles[i].next = i==numtiles-1 ? NULL : &tiles[i+1];

      /* For a check:
      printf("%zu:\n\tStart index: %zu\n\tsize: %zu x %zu\n", i, tind,
             tiles[i].dsize[1], tiles[i].dsize[0]);
      exit(0);
      */
    }


  /* Clean up and return. */
  free(last);
  free(first);
  free(coord);
  free(tcoord);
  if(start) free(start);
  return tsize;
}





/* Make sure that the input parameters (in `tl', short for two-layer) fit
   with the input dataset. The filename and HDU are only required for error
   messages. Also, allocate and fill the `channelsize' array. */
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
        error(EXIT_FAILURE, 0, "`--tilesize' must be larger than zero, "
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
          "value(s) given for the tile size (`--tilesize' option).",
          filename, hdu, ndim, i);


  /* Check the channel's dimensions. */
  for(i=0; tl->numchannels[i]!=-1; ++i)
    if(tl->numchannels[i]==0)
      error(EXIT_FAILURE, 0, "the number of channels in all dimensions must "
            "be larger than zero. The number for dimension %zu was zero",
            ndim);
  if(i!=ndim)
    error(EXIT_FAILURE, 0, "%s (hdu: %s): has %zu dimensions, but only %zu "
          "value(s) given for the number of channels", filename, hdu, ndim,
          i);


  /* Allocate space for the channel sizes. */
  tl->channelsize=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, ndim);


  /* Check if the channels are exactly divisible by the input's size along
     each dimension and set the correct size. */
  for(i=0;i<ndim;++i)
    {
      /* Check if the number of channels is not more than the size of the
         image. Note that the reported dimension must be in FITS format.*/
      if( input->dsize[i] < tl->numchannels[i] )
        error(EXIT_FAILURE, 0, "the number of channels in dimension %zu "
              "(%zu) is more than the size of the `%s' (hdu: %s) in that "
              "dimension", ndim-i, tl->numchannels[i], filename, hdu);

      /* Also check the tile size. */
      if( input->dsize[i] < tl->tilesize[i] )
        error(EXIT_FAILURE, 0, "the tile size in dimension %zu (%zu) is "
              "more than the size of the `%s' (hdu: %su) in that dimension",
              ndim-i, tl->tilesize[i], filename, hdu);

      /* First check. */
      d=(double)input->dsize[i]/(double)(tl->numchannels[i]);
      if(ceil(d)!=d)
        error(EXIT_FAILURE, 0, "%zu (number of channels along dimension "
              "%zu) is not exactly divisible by %zu (the length of `%s' "
              "(hdu: %s) that dimension). The channels cover the input "
              "dataset, hence, they must be identical", tl->numchannels[i],
              ndim-i, input->dsize[i], filename, hdu);

      /* Put the channel size into the output. */
      tl->channelsize[i]=d;
    }
}





/* A dataset can be tiled with two layers that are related:

      Channels: A tesselation of larger tile sizes that all have the same
           size (`channel_size' must be an exact multiple of `input's size
           along every dimension. In astronomy images, this can be seen as
           CCD amplifiers, that cover large parts of the image. If
           `*channels!=NULL' then it is assumed to be already present and
           will not be allocated.

      Tiles: A combined tesselation of each channel with smaller
           tiles. These tiles can be used to calculate things like
           gradients over each channel and thus over the whole image. */
void
gal_tile_full_two_layers(gal_data_t *input,
                         struct gal_tile_two_layer_params *tl)
{
  gal_data_t *t;
  size_t i, *junk, ndim=tl->ndim=input->ndim;

  /* Initialize. Note that `numchannels might have already been
     allocated. */
  tl->channels=tl->tiles=NULL;
  if(tl->numchannels) free(tl->numchannels);

  /* Initialize necessary values and do the channels tessellation. */
  tl->numchannels = gal_tile_full(input, tl->channelsize, tl->remainderfrac,
                                &tl->channels, 1);
  tl->totchannels = gal_multidim_total_size(ndim, tl->numchannels);


  /* Tile each channel. While tiling the first channel, we are also going
     to allocate the space for the other channels. Then pass those pointers
     when we want to fill in each tile of the other channels. */
  tl->numtilesinch = gal_tile_full(tl->channels, tl->tilesize,
                                   tl->remainderfrac, &tl->tiles,
                                   tl->totchannels);
  tl->tottilesinch = gal_multidim_total_size(ndim, tl->numtilesinch);
  for(i=1; i<tl->totchannels; ++i)
    {
      /* Set the first tile in this channel. Then use it it fill the `next'
         pointer of the previous channel's tiles. Note that `gal_tile_full'
         set this `next' element to NULL. */
      t = tl->tiles + i * tl->tottilesinch;
      tl->tiles[ i * tl->tottilesinch - 1 ].next = t;

      /* Fill in the information for all the tiles in this channel. Note
         that we already have the returned value, so it isn't important.*/
      junk=gal_tile_full(&tl->channels[i], tl->tilesize, tl->remainderfrac,
                         &t, 1);
      free(junk);
    }

  /* Multiply the number of tiles along each dimension OF ONE CHANNEL by
     the number of channels in each dimension to get the dimensionality of
     the full tile structure. */
  tl->numtiles = gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, ndim);
  for(i=0;i<ndim;++i)
    tl->numtiles[i] = tl->numtilesinch[i] * tl->numchannels[i];
  tl->tottiles = gal_multidim_total_size(ndim, tl->numtiles);
}





/* If necessary (see below), make a new Re-ordered array that has one
   pixel/element for each tile, such that it is ready for displaying. When
   it isn't necessary, this function will just return the input data
   structure. So you should free the output only if it is different from
   the input.

   When there is only one channel OR one dimension, you can just save the
   array and everything will be fine. However, when there is more than one
   channel AND more than one dimension this won't work and the values will
   not be in the places corresponding to the original dataset. This happens
   because the tiles (and thus the value you assign to each tile) in each
   channel are taken to be a contiguous patch of memory. But when you look
   at the whole tessellation, you want the tiles to be ordered based on
   their position in the final dataset/image.

   The example below may help clarify: assume you have a 6x6 tessellation
   with two channels in the horizontal and one in the vertical. On the left
   you can see how the tiles (and their values) are allocated in memory. On
   the right is how they will be displayed if you just save it as a FITS
   file with no modification (this will not correspond to your input
   dataset).

              Allocation map                  When displayed
              --------------                  --------------
              15 16 17 33 34 35               30 31 32 33 34 35
              12 13 14 30 31 32               24 25 26 27 28 29
              09 10 11 27 28 29               18 19 20 21 22 23
              06 07 08 24 25 26               12 13 14 15 16 17
              03 04 05 21 22 23               06 07 08 09 10 11
              00 01 02 18 19 20               00 01 02 03 04 05       */
gal_data_t *
gal_tile_full_values_for_disp(gal_data_t *tilevalues,
                              struct gal_tile_two_layer_params *tl)
{
  void *start;
  gal_data_t *out, *fake;
  size_t o_inc, num_increment, start_ind;
  size_t i, ch_ind, contig, i_inc=0, ndim=tl->ndim;
  size_t *chstart=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, ndim);
  size_t *s_e_inc=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, ndim);

  /* Make sure that the input and tessellation correspond to each
     other. This is important, because the user might mistakenly give the
     input dataset as input, not the dataset that has one value per tile.*/
  if(tilevalues->ndim!=tl->ndim)
    error(EXIT_FAILURE, 0, "the input (`tilevalues') and tessellation "
          "(`tl') in `gal_tile_full_values_for_disp' have %zu and %zu "
          "dimensions respectively! They must have the same number of "
          "dimensions", tilevalues->ndim, tl->ndim);
  for(i=0; i<ndim; ++i)
    if(tilevalues->dsize[i]!=tl->numtiles[i])
      error(EXIT_FAILURE, 0, "the total number of tiles (in all channels) "
            "of dimension %zu does not correspond between the input "
            "(`tilevalues': %zu) and tessellation (`tl': %zu)", ndim-i,
            tilevalues->dsize[i], tl->numtiles[i]);

  /* If there is only one dimension or one channel, then just return the
     input `tilevalues'. */
  if(tilevalues->ndim==1 || tl->totchannels==1)
    return tilevalues;

  /* Allocate the output. */
  out=gal_data_alloc(NULL, tilevalues->type, tilevalues->ndim,
                     tilevalues->dsize, tilevalues->wcs, 0,
                     tilevalues->minmapsize, tilevalues->name,
                     tilevalues->unit, tilevalues->comment);

  /* Allocate a fake tile (with a size equal to the number of tiles in one
     channel) for easy parsing with standard techniques. */
  fake=gal_data_alloc(NULL, GAL_DATA_TYPE_UINT8, tilevalues->ndim,
                      tl->numtilesinch, NULL, 0, tilevalues->minmapsize,
                      NULL, NULL, NULL);

  /* Initialize the fake array: we will set it's `block' to the output.*/
  fake->block=out;
  free(fake->array);
  fake->type=GAL_DATA_TYPE_INVALID;

  /* Put the values into the proper place. */
  contig=tl->numtilesinch[ndim-1];
  for(ch_ind=0; ch_ind<tl->totchannels; ++ch_ind)
    {
      /* Set the starting pointer of this channel for the `fake' tile that
         represents it. */
      gal_multidim_index_to_coord(ch_ind, ndim, tl->numchannels, chstart);
      for(i=0;i<ndim;++i) chstart[i]*=tl->numtilesinch[i];
      start_ind=gal_multidim_coord_to_index(ndim, tl->numtiles, chstart);
      fake->array=gal_data_ptr_increment(out->array, start_ind, out->type);

      /* Find the starting and ending tile indexs (inclusive in the output)
         of this channel. */
      start=gal_tile_start_end_ind_inclusive(fake, out, s_e_inc);

      /* Start parsing the channel and putting it in the output. */
      o_inc=0;
      num_increment=1;
      while(s_e_inc[0] + o_inc <= s_e_inc[1])
        {
          /* Copy the contiguous region from the `tilevalues' array to the
             output array. Note that since they have the same type, we can
             just use `memcpy' and don't actually have to parse it.*/
          memcpy(gal_data_ptr_increment(start,             o_inc, out->type),
                 gal_data_ptr_increment(tilevalues->array, i_inc, out->type),
                 gal_data_sizeof(out->type)*contig);

          /* Set the incrementation for the next contiguous patch. */
          o_inc += gal_tile_block_increment(out, fake->dsize,
                                            num_increment++, NULL);
          i_inc += contig;
        }
    }

  /* Clean up and return. */
  fake->array=NULL;
  gal_data_free(fake);
  free(chstart);
  free(s_e_inc);
  return out;
}





void
gal_tile_full_values_write(gal_data_t *tilevalues,
                           struct gal_tile_two_layer_params *tl,
                           char *filename, char *program_string)
{
  gal_data_t *disp;

  /* Make the dataset to be displayed. */
  disp = ( tl->oneelempertile
           ? gal_tile_full_values_for_disp(tilevalues, tl)
           : gal_tile_block_write_const_value(tilevalues, tl->tiles, 0) );

  /* Write the array as a file and then clean up. */
  gal_fits_img_write(disp, filename, NULL, program_string);
  gal_data_free(disp);
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

  /* Free the arrays of `gal_data_c' for each tile and channel. */
  if(tl->tiles)    gal_data_array_free(tl->tiles,    tl->tottiles,    0);
  if(tl->channels) gal_data_array_free(tl->channels, tl->totchannels, 0);
}



















/*********************************************************************/
/********************         On threads          ********************/
/*********************************************************************/
/* Run a given function on the given tiles. */
void
gal_tile_function_on_threads(gal_data_t *tiles, void *(*function)(void *),
                             size_t numthreads, void *caller_params)
{
  int err;
  pthread_t t;          /* All thread ids saved in this, not used. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct gal_tile_thread_param *prm;
  size_t i, *indexs, thrdcols, numbarriers;
  size_t numtiles=gal_data_num_in_ll(tiles);

  /* Allocate the array of parameters structure structures. */
  prm=malloc(numthreads*sizeof *prm);
  if(prm==NULL)
    {
      fprintf(stderr, "%zu bytes could not be allocated for prm.",
              numthreads*sizeof *prm);
      exit(EXIT_FAILURE);
    }

  /* Distribute the actions into the threads: */
  gal_threads_dist_in_threads(numtiles, numthreads, &indexs, &thrdcols);

  /* Do the job: when only one thread is necessary, there is no need to
     spin off one thread, just call the function directly (spinning off
     threads is expensive). This is for the generic thread spinner
     function, not this simple function where `numthreads' is a
     constant. */
  if(numthreads==1)
    {
      prm[0].id=0;
      prm[0].b=NULL;
      prm[0].indexs=indexs;
      prm[0].params=caller_params;
      function(&prm[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
         (that spinns off the nt threads) is also a thread, so the
         number the barriers should be one more than the number of
         threads spinned off. */
      numbarriers = (numtiles<numthreads ? numtiles : numthreads) + 1;
      gal_threads_attr_barrier_init(&attr, &b, numbarriers);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
        if(indexs[i*thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            prm[i].id=i;
            prm[i].b=&b;
            prm[i].params=caller_params;
            prm[i].indexs=&indexs[i*thrdcols];
            err=pthread_create(&t, &attr, function, &prm[i]);
            if(err)
              {
                fprintf(stderr, "can't create thread %zu", i);
                exit(EXIT_FAILURE);
              }
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }

  /* Clean up. */
  free(indexs);
}
