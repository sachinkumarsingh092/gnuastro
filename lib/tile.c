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

#include <gnuastro/tile.h>
#include <gnuastro/multidim.h>






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





/* Calculate the starting coordinates of a tile in the allocated block of
   memory. */
void
gal_tile_block_tile_start_coord(gal_data_t *tile, size_t *start_coord)
{
  size_t *s, *sf, ind;
  gal_data_t *block=gal_tile_block(tile);

  /* If the input tile is actually the same as the block, then the
     reference is all zeros. */
  if(block==tile)
    {
      sf = (s=start_coord) + tile->ndim;
      do *s++=0; while(s<sf);
      return;
    }

  /* Calculate the coordinates of the first pixel of the tile. */
  ind = tile->array - block->array;
  gal_multidim_index_to_coord(ind, tile->ndim, block->dsize, start_coord);
}



















/***********************************************************************/
/**************           Tile full dataset         ********************/
/***********************************************************************/
/* This option is mainly used by Gnuastro's programs to see if the user's
   given values for the tile size (`tile') and number of channels
   (`channels') corresponds to the input dataset. */
size_t *
gal_tile_all_sanity_check(char *filename, char *hdu, gal_data_t *input,
                          size_t *tile, size_t *numchannels)
{
  double d;
  size_t i, *out;

  /* Check the tile's dimensions. */
  for(i=0;tile[i]!=-1;++i)
    {
      /* Not equal to zero. */
      if(tile[i]==0)
        error(EXIT_FAILURE, 0, "the tile size must be larger than zero");

      /* If the tile size is larger than the dataset size in this
         dimension, then change the tile size to the dataset size. */
      if(tile[i]>input->dsize[i]) tile[i]=input->dsize[i];
    }


  /* Make sure the number of tile sizes (tile dimensions) are the same as
     the dataset's dimensions). */
  if(i!=input->ndim)
    error(EXIT_FAILURE, 0, "%s (hdu: %s): has %zu dimensions, but only %zu "
          "value(s) given for the tile size", filename, hdu, input->ndim, i);


  /* Check the channel's dimensions. */
  for(i=0;numchannels[i]!=-1;++i)
    if(tile[i]==0)
      error(EXIT_FAILURE, 0, "the number of channels must be larger than "
            "zero");
  if(i!=input->ndim)
    error(EXIT_FAILURE, 0, "%s (hdu: %s): has %zu dimensions, but only %zu "
          "value(s) given for the number of channels", filename, hdu,
          input->ndim, i);


  /* Allocate space for the channel sizes. */
  out=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, input->ndim);


  /* Check if the channels are exactly divisible by the input's size along
     each dimension and set the correct size. */
  for(i=0;i<input->ndim;++i)
    {
      /* Check if the number of channels is not more than the size of the
         image. Note that the reported dimension must be in FITS format.*/
      if(input->dsize[i]<numchannels[i])
        error(EXIT_FAILURE, 0, "the number of channels in dimension %zu "
              "(%zu) is more than the size of the `%s' (hdu: %s) in that "
              "dimension", input->ndim-i, numchannels[i], filename, hdu);

      /* Also check the tile size. */
      if(input->dsize[i]<tile[i])
        error(EXIT_FAILURE, 0, "the tile size in dimension %zu (%zu) is "
              "more than the size of the `%s' (hdu: %su) in that dimension",
              input->ndim-i, tile[i], filename, hdu);

      /* First check. */
      d=(double)input->dsize[i]/(double)numchannels[i];
      if(ceil(d)!=d)
        error(EXIT_FAILURE, 0, "%zu (number of channels along dimension "
              "%zu) is not exactly divisible by %zu (the length of `%s' "
              "(hdu: %s) that dimension). The channels cover the input "
              "dataset, hence, they must be identical", numchannels[i],
              input->ndim-i, input->dsize[i], filename, hdu);

      /* Put the channel size into the output. */
      out[i]=d;
    }

  /* Return the output array (size of channel along each dimension). */
  return out;
}





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
gal_tile_all_regular_first(gal_data_t *parent, size_t *regular,
                           float significance, size_t *first, size_t *last,
                           size_t *number)
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
              number[i] = dsize[i]/regular[i] + 1 ;
              last[i]   = dsize[i] - ( first[i] + regular[i]*(number[i]-2) );
            }
          else
            {
              first[i]  = remainder + regular[i];
              number[i] = dsize[i]/regular[i];
              last[i]   = regular[i];
            }
        }
      else
        {
          first[i]  = last[i] = regular[i];
          number[i] = dsize[i]/regular[i];
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

     `out' is the pointer to the array of data structures that is to keep
        the tile parameters. If `*out==NULL', then the necessary space will
        be allocated. If it is not NULL, then all the tile information will
        be filled from the given element, see `multiple' for more.

     `multiple': When the `*out' array is to be allocated, allocate
        `multiple' times the necessary space. This can be very useful when
        you have several more identically sized `inputs', and you want all
        their tiles to be allocated (and thus indexed) together, even
        though they have different `block' datasets (that then link to one
        allocated space).  See the `gal_tile_all_position_channel_tile'
        below.

   Output
   ------

     The returned output is an array of `gal_data_t' with `numtiles'
     elements, note that you have to pass the pointer to `numtiles' as one
     of the arguments.


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
size_t
gal_tile_all_position(gal_data_t *input, size_t *regular, gal_data_t **out,
                      size_t multiple)
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
  gal_tile_all_regular_first(input, regular, 0.3, first, last, tsize);
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
      gal_tile_block_tile_start_coord(input, start);
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
             `input'. See the comments above `gal_tile_all_regular_first':
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
      tiles[i].array=block->array+tind;

      /* Set the sizes of the tile. */
      tiles[i].ndim=input->ndim;
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
        }

      /* Set the block structure for this tile to the `input'. */
      tiles[i].block=input;

      /* For a check:
      printf("%zu:\n\tStart index: %zu\n\tsize: %zu x %zu\n", i, tind,
             tiles[i].dsize[1], tiles[i].dsize[0]);
      exit(0);
      */
    }


  /* Clean up and return. */
  free(last);
  free(tsize);
  free(first);
  free(coord);
  free(tcoord);
  if(start) free(start);
  return numtiles;
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
           gradients over each channel and thus over the whole image.  */
void
gal_tile_all_position_two_layers(gal_data_t *input, size_t *channel_size,
                                 size_t *tile_size, gal_data_t **channels,
                                 gal_data_t **tiles, size_t *numchannels,
                                 size_t *numtiles)
{
  gal_data_t *ch, *t;
  size_t i, nch, ntiles_in_ch;

  /* First allocate the channels. Note that the channels tesselation might
     have already been set. */
  if(*channels)
    {
      nch=1;
      for(i=0;i<input->ndim;++i)
        nch *= input->dsize[i]/channel_size[i];
    }
  else
    /* Note that the actual allocated input array will be the direct
       `block' of each channel. */
    nch=gal_tile_all_position(input, channel_size, channels, 1);


  /* Now, tile each channel. While tiling the first channel, we are also
     going to allocate the space for the other channels. Then pass those
     pointers. */
  *tiles=NULL;
  ch=*channels;
  ntiles_in_ch = gal_tile_all_position(ch, tile_size, tiles, nch);
  for(i=1;i<nch;++i)
    {
      t = *tiles + i*ntiles_in_ch;
      gal_tile_all_position(&ch[i], tile_size, &t, 1);
    }

  /* Return the total number of channels and tiles */
  *numchannels = nch;
  *numtiles    = nch * ntiles_in_ch;
}
