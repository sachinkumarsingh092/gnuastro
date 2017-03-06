/*********************************************************************
convolve -- Convolve a dataset with a given kernel.
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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/tile.h>
#include <gnuastro/threads.h>
#include <gnuastro/multidim.h>
#include <gnuastro/convolve.h>








/*********************************************************************/
/********************          Utilities          ********************/
/*********************************************************************/
/* See if the tile is on the edge of the hosted region or not. It doesn't
   matter if the host is the allocated block of memory or a region in it (a
   channel). */
static int
convolve_tile_is_on_edge(size_t *h, size_t *start_end_coord, size_t *k,
                         size_t ndim)
{
  size_t *s=start_end_coord, *e=start_end_coord+ndim, *kf=k+ndim;

  /* If the starting point of the tile is smaller than the half-kernel
     length along that dimension, then the tile is on the edge.

     If the end of the tile in this dimension added with the half-kernel
     length along that dimension is larger than the host's size, then the
     tile is on the edge. */
  do if( (*s++ < *k/2) || (*e++ + *k/2 > *h++) ) return 1; while(++k<kf);

  /* If control reaches here, then this is a central tile. */
  return 0;
}




















/*********************************************************************/
/********************     Spatial convolution     ********************/
/*********************************************************************/
struct spatial_params
{
  gal_data_t    *out;
  gal_data_t  *tiles;
  gal_data_t *kernel;
  int     convoverch;
};





/* Convolve over one tile that is not touching the edge. */
static void
convolve_spatial_tile(struct spatial_params *cprm, size_t tile_ind,
                      size_t *parse_coords, gal_data_t *overlap)
{
  gal_data_t *tile=&cprm->tiles[tile_ind];

  int on_edge;
  double sum, ksum;
  size_t i, ndim=tile->ndim, csize=tile->dsize[ndim-1];
  gal_data_t *block=gal_tile_block(tile), *kernel=cprm->kernel;

  size_t i_inc, i_ninc, i_st_en[2], o_inc, o_ninc, o_st_en[2];
  size_t *p, *pf, *k, *o, *e, o_st_inc, k_start, start_fastdim;

  /* These variables depend on the type of the input. */
  float *kv, *iv, *ivf, *i_start, *o_start;
  float *in_v, *in=block->array, *out=cprm->out->array;

  /* The `parse_coords' array was allocated once before this function for
     all the tiles that are given to a thread. It has the space for all the
     necessary coordinates. At first, `pix' will be the starting element of
     the tile, but then it will be incremented as we carpet the tile. So we
     aren't calling it `start'. */
  size_t *pix           = parse_coords;
  size_t *end           = parse_coords + ndim;
  size_t *overlap_start = parse_coords + ndim * 2;
  size_t *o_c           = parse_coords + ndim * 3;


  /* Set the starting and ending coordinates of this tile (recall that the
     start and end are the first two allocated spaces in
     parse_coords). When `convoverch' is set, we want to convolve over the
     whole allocated block, not just one channel. So in effect, it is the
     same as `rel_block' in `gal_tile_start_end_coord'. */
  gal_tile_start_end_coord(tile, parse_coords, cprm->convoverch);
  start_fastdim = pix[ndim-1];


  /* See if this tile is on the edge or not. */
  on_edge=convolve_tile_is_on_edge( ( cprm->convoverch
                                      ? block->dsize
                                      : tile->block->dsize ),
                                    parse_coords, kernel->dsize, ndim);


  /* Go over the tile, first find its limits, then loop over the fastest
     dimension. */
  i_inc=0; i_ninc=1;
  i_start=gal_tile_start_end_ind_inclusive(tile, block, i_st_en);
  while( i_st_en[0] + i_inc <= i_st_en[1] )
    {
      /* Initialize the value along the fastest dimension (it is not
         incremented during `gal_tile_block_increment'). */
      pix[ndim-1]=start_fastdim;

      /* Go over each pixel to convolve. */
      for(i=0;i<csize;++i)
        {
          /* Pointer to the pixel under consideration. */
          in_v = i_start + i_inc + i;

          /* If the input on this pixel is a NaN, then just set the output
             to NaN too and go onto the next pixel. `in_v' is the pointer
             on this pixel. */
          if( isnan(*in_v) )
            out[ in_v - in ]=NAN;
          else
            {
              /* Coordinate to start convolution for this pixel. */
              pf=(p=pix)+ndim; e=end;
              o=overlap_start; k=kernel->dsize;
              do
                {
                  if(on_edge)
                    {
                      return;
                    }
                  else { *o++ = *p - *k++/2; k_start=0; }
                }
              while(++p<pf);
              o_st_inc=gal_multidim_coord_to_index(ndim, block->dsize,
                                                   overlap_start);

              /* Set the overlap parameters, then set the region. */
              overlap->dsize=kernel->dsize;
              overlap->array=gal_data_ptr_increment(block->array, o_st_inc,
                                                    block->type);
              o_start=gal_tile_start_end_ind_inclusive(overlap, block,
                                                       o_st_en);

              /* Go over the kernel-overlap region. */
              kv = ( k_start
                     ? gal_data_ptr_increment(kernel->array, k_start,
                                              kernel->type)
                     : kernel->array );
              ksum=sum=0.0f; o_inc=0; o_ninc=1;
              while( o_st_en[0] + o_inc <= o_st_en[1] )
                {
                  /* Go over the contiguous region. */
                  ivf = ( iv = o_start + o_inc ) + overlap->dsize[ndim-1];
                  do
                    {
                      if( !isnan(*iv) )
                        { ksum += *kv; sum += *iv * *kv; }
                      ++kv;
                    }
                  while(++iv<ivf);

                  /* Update the incrementation to the next contiguous
                     region of memory over this tile. */
                  o_inc += gal_tile_block_increment(block, overlap->dsize,
                                                    o_ninc++, o_c);
                }

              /* Set the output value. */
              out[ in_v - in ] = ksum==0.0f ? NAN : sum/ksum;
            }

          /* Increment the last coordinate. */
          pix[ndim-1]++;
        }

      /* Increase the increment from the start of the tile for the next
         contiguous patch. */
      i_inc += gal_tile_block_increment(block, tile->dsize, i_ninc++, pix);
    }
}





/* Do spatial convolution on each mesh. */
static void *
convolve_spatial_on_thread(void *inparam)
{
  struct gal_tile_thread_param *tprm=(struct gal_tile_thread_param *)inparam;
  struct spatial_params *cprm=(struct spatial_params *)(tprm->params);

  size_t i;
  gal_data_t overlap, *block=gal_tile_block(cprm->tiles);
  size_t *parse_coords=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T,
                                             4*block->ndim);

  /* Initialize the necessary overlap parameters. */
  overlap.block=block;
  overlap.ndim=block->ndim;
  overlap.dsize=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, block->ndim);


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    convolve_spatial_tile(cprm, tprm->indexs[i], parse_coords, &overlap);


  /* Clean up, wait until all other threads finish, then return. In a
     single thread situation, `tprm->b==NULL'. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* Convolve a dataset with a given kernel in the spatial domain. Since
   spatial convolution can be very series of tiles arranged as an array. */
gal_data_t *
gal_convolve_spatial(gal_data_t *tiles, gal_data_t *kernel,
                     size_t numthreads, int convoverch)
{
  struct spatial_params params;
  gal_data_t *out, *block=gal_tile_block(tiles);

  /* Small sanity checks. */
  if(tiles->ndim!=kernel->ndim)
    error(EXIT_FAILURE, 0, "The number of dimensions between the kernel and "
          "input should be the same in `gal_convolve_spatial'");
  if( block->type!=GAL_DATA_TYPE_FLOAT32
      || kernel->type!=GAL_DATA_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "`gal_convolve_spatial' currently only works on "
          "`float32' type input and kernel");

  /* Allocate the output convolved dataset. */
  out=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT32, block->ndim, block->dsize,
                     block->wcs, 0, block->minmapsize, "CONVOLVED",
                     block->unit, NULL);
  { float *f=out->array, *ff=f+out->size; do *f=NAN; while(++f<ff); }

  /* Set the pointers in the parameters structure. */
  params.out=out;
  params.tiles=tiles;
  params.kernel=kernel;
  params.convoverch=convoverch;

  /* Do the spatial convolution on threads. */
  gal_tile_function_on_threads(tiles, convolve_spatial_on_thread,
                               numthreads, &params);

  /* Clean up and return the output array. */
  return out;
}
