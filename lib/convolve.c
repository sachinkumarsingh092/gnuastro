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
#include <gnuastro/convolve.h>
#include <gnuastro/dimension.h>

#include <gnuastro-internal/checkset.h>






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
     length along that dimension is equal to or larger than the host's
     size, then the tile is on the edge. */
  do if( (*s++ < *k/2) || (*e++ + *k/2 >= *h++) ) return 1; while(++k<kf);

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
  int edgecorrection;
};




/* Define the overlap of the kernel and image over this part of the image,
   the necessary input image parameters are stored in `overlap' (its
   `array' and `dsize' elements). The pointer to the kernel array that the
   overlap starts will be returned. */
static int
convolve_spatial_overlap(int on_edge, size_t *parse_coords, size_t *hsize,
                         gal_data_t *block, gal_data_t *kernel,
                         gal_data_t *overlap, size_t *k_start_inc,
                         size_t tile_ind)
{
  size_t overlap_inc, ndim=kernel->ndim;
  int full_overlap=1, dim_full_overlap, is_start, is_end;

  size_t *pix           = parse_coords; /* needs 2*ndim, also for end. */
  size_t *overlap_start = parse_coords + ndim * 3;
  size_t *kernel_start  = parse_coords + ndim * 4;
  size_t *host_start    = parse_coords + ndim * 5;

  size_t *p=pix, *pf=pix+ndim, *os=overlap_start, *h=hsize;
  size_t *k=kernel->dsize, *ks=kernel_start, *od=overlap->dsize;
  /*
  if(tile_ind==2053)
    printf("pix: %zu, %zu\n", pix[0], pix[1]);
  */
  /* Coordinate to start convolution for this pixel. */
  do
    {
      /* Initialize the overlap for this dimension (we'll assume it
         overlaps because this is the most common case usually). */
      dim_full_overlap=1;


      /* When the tile is on the edge, still some pixels in it can have
         full overlap. So using the `dim_full_overlap', we will do the same
         thing we do for the tiles that don't overlap for them. */
      if(on_edge)
        {
          /* See if this pixel is on the start and/or end of the dimension
             relative to the kernel. */
          is_start = *p < *k/2;
          is_end   = *p + *k/2 >= *h;
          if( is_start || is_end )
            {
              /* Overlapping with the outside.

                 With the start: assume that in this dimension, the pixel
                 is at position 2, while the kernel is 11 pixels wide (or 5
                 pixels in half-width). As seen below, the kernel should
                 start from pixel `5-2=3' in this dimension and the overlap
                 size should decrease by the same amount.

                    image:            0 1 2 3 4 5 6 7 8 9 ...
                    pixel:                p
                    kernel:     0 1 2 3 4 5 6 7 8 9 10

                 With the end: Similar to above, but assume the pixel is
                 two pixels away from the edge of a 100-pixel image. We are
                 no longer worried about the overlap or kernel starting
                 point, it is the width that we need to decrease it by:

                   97 + 5 - 100 + 1 : The `1' is because we want the pixel
                                      immediately after the end.

                    image:        ... 92 93 94 95 96 97 98 99 | 100 101 102
                    pixel:                            p
                    kernel:                 0 1 2 3 4 5 6 7 8 | 9   10    */
              *ks++ = is_start ? *k/2 - *p : 0;
              *os++ = is_start ? 0         : *p - *k/2;

              /* We will start with the full kernel width, then decrease it
                 if the pixel is too close to the start or end along this
                 dimension. Note that the host array/image might actually
                 be smaller than kernel, so both cases might occur. */
              *od = *k;
              if(is_start) *od -= *k/2 - *p;
              if(is_end)   *od -= *p + *k/2 - *h + 1;

              /* Increment and finalize. */
              ++h;
              ++k;
              ++od;
              full_overlap=0;
              dim_full_overlap=0;
            }
        }

      /* There is full overlap for this pixel or tile over this
         dimension.  */
      if(dim_full_overlap)
        {
          /* Set the values. */
          *ks++ = 0;
          *od++ = *k;
          *os++ = *p - *k/2;

          /* Increment. */
          ++h;
          ++k;
        }
    }
  while(++p<pf);
  /*
  if(tile_ind==2053)
    {
      printf("\tk/2: %zu, %zu\n", kernel->dsize[0]/2, kernel->dsize[1]/2);
      printf("\toverlap_start: %zu, %zu\n", overlap_start[0],
             overlap_start[1]);
      printf("\toverlap->dsize: %zu, %zu\n", overlap->dsize[0],
             overlap->dsize[1]);
      exit(0);
    }
  */
  /* Set the increment to start working on the kernel. */
  *k_start_inc = ( full_overlap
                   ? 0
                   : gal_dimension_coord_to_index(ndim, kernel->dsize,
                                                  kernel_start) );

  /* Add the host's starting location (necessary when convolution over the
     host/channel is treated independently). Until now we worked as if the
     the host/channel is the full image so the edges don't get mixed. But
     from now on we will be working over the allocated block to look at
     pixel values, so we need to convert the location to the proper place
     within the allocated array. */
  gal_dimension_add_coords(overlap_start, host_start, overlap_start, ndim);

  /* Set the increment to start working on the overlap region and use that
     to set the starting pointer of the overlap region. */
  overlap_inc=gal_dimension_coord_to_index(ndim, block->dsize, overlap_start);
  overlap->array=gal_data_ptr_increment(block->array, overlap_inc,
                                        block->type);
  return full_overlap;
}





/* Convolve over one tile that is not touching the edge. */
static void
convolve_spatial_tile(struct spatial_params *cprm, size_t tile_ind,
                      size_t *parse_coords, gal_data_t *overlap)
{
  gal_data_t *tile=&cprm->tiles[tile_ind];
  gal_data_t *block=gal_tile_block(tile), *kernel=cprm->kernel;

  double sum, ksum;
  int on_edge, full_overlap;
  size_t i, ndim=tile->ndim, csize=tile->dsize[ndim-1];
  gal_data_t *host=cprm->convoverch ? block : tile->block;

  /* Variables for scanning a tile (`i_*') and the region around every
     pixel of a tile (`o_*'). */
  size_t start_fastdim, k_start_inc;
  size_t k_inc, i_inc, i_ninc, i_st_en[2], o_inc, o_ninc, o_st_en[2];

  /* These variables depend on the type of the input. */
  float *kv, *iv, *ivf, *i_start, *o_start, *k_start;
  float *in_v, *in=block->array, *out=cprm->out->array;

  /* The `parse_coords' array was allocated once before this function for
     all the tiles that are given to a thread. It has the space for all the
     necessary coordinates. At first, `pix' will be the starting element of
     the tile, but then it will be incremented as we carpet the tile. So we
     aren't calling it `start'. */
  size_t *pix           = parse_coords; /* needs 2*ndim (also for end). */
  size_t *o_c           = parse_coords + ndim * 2;
  size_t *host_start    = parse_coords + ndim * 5;


  /* Set the starting and ending coordinates of this tile (recall that the
     start and end are the first two allocated spaces in
     parse_coords). When `convoverch' is set, we want to convolve over the
     whole allocated block, not just one channel. So in effect, it is the
     same as `rel_block' in `gal_tile_start_end_coord'. */
  gal_tile_start_end_coord(tile, parse_coords, cprm->convoverch);
  start_fastdim = pix[ndim-1];

  /* Set the starting coordinate of the host for this tile. */
  gal_tile_start_coord(host, host_start);

  /* See if this tile is on the edge or not. */
  on_edge=convolve_tile_is_on_edge(host->dsize, parse_coords, kernel->dsize,
                                   ndim);
  /*
  if(tile_ind==2053)
    {
      printf("\ntile %zu...\n", tile_ind);
      printf("\tpix: %zu, %zu\n", pix[0], pix[1]);
      exit(0);
    }
  */
  /* Go over the tile. */
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
              /* Define the overlap region. */
              full_overlap=convolve_spatial_overlap(on_edge, parse_coords,
                                                    host->dsize, block,
                                                    kernel, overlap,
                                                    &k_start_inc, tile_ind);

              /* Set the starting pixel over the image (`o_start'). */
              o_start=gal_tile_start_end_ind_inclusive(overlap, block,
                                                       o_st_en);

              /* Set the starting kernel pixel. Note that `kernel_array' is
                 `void *' (pointer arithmetic is not defined on it). So we
                 will first put it in `k_start, and then increment that. */
              k_start=kernel->array; k_start+=k_start_inc;

              /* Go over the kernel-overlap region. */
              ksum = cprm->edgecorrection ? 0.0f : 1.0f;
              sum=0.0f; k_inc=0; o_inc=0; o_ninc=1; kv=k_start;
              while( o_st_en[0] + o_inc <= o_st_en[1] )
                {
                  /* Go over the contiguous region. When there is full
                     overlap, we don't need to calculate incrementation
                     over the kernel, it is always a single
                     incrementation. But when we have partial overlap,
                     we'll need to calculate a different incrementation. */
                  ivf = ( iv = o_start + o_inc ) + overlap->dsize[ndim-1];
                  if(full_overlap==0) kv = k_start + k_inc;
                  do
                    {
                      if( !isnan(*iv) )
                        {
                          sum += *iv * *kv;
                          if(cprm->edgecorrection) ksum += *kv;
                        }
                      ++kv;
                    }
                  while(++iv<ivf);

                  /* Update the incrementation to the next contiguous
                     region of memory over this tile. Note that the
                     contents of `o_c' are irrelevant here.*/
                  o_inc += gal_tile_block_increment(block, overlap->dsize,
                                                    o_ninc++, o_c);
                  if(full_overlap==0)
                    k_inc += gal_tile_block_increment(kernel, overlap->dsize,
                                                      o_ninc-1, NULL);
                }

              /* Set the output value. */
              out[ in_v - in ] = ( ksum==0.0f
                                   ? NAN
                                   : sum/ksum );
            }

          /* Increment the last coordinate. */
          pix[ndim-1]++;
        }

      /* Increase the increment from the start of the tile for the next
         contiguous patch. */
      i_inc += gal_tile_block_increment(block, tile->dsize, i_ninc++, pix);
    }
  /*
  if(tile_ind==2053)
    printf("... done.\n");
  */
}





/* Do spatial convolution on each mesh. */
static void *
convolve_spatial_on_thread(void *inparam)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)inparam;
  struct spatial_params *cprm=(struct spatial_params *)(tprm->params);

  size_t i;
  gal_data_t overlap, *block=gal_tile_block(cprm->tiles);
  size_t *parse_coords=gal_data_malloc_array(GAL_TYPE_SIZE_T, 6*block->ndim);

  /* Initialize the necessary overlap parameters. */
  overlap.block=block;
  overlap.ndim=block->ndim;
  overlap.dsize=gal_data_malloc_array(GAL_TYPE_SIZE_T, block->ndim);


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    convolve_spatial_tile(cprm, tprm->indexs[i], parse_coords, &overlap);


  /* Clean up, wait until all other threads finish, then return. In a
     single thread situation, `tprm->b==NULL'. */
  free(parse_coords);
  free(overlap.dsize);
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* Convolve a dataset with a given kernel in the spatial domain. Spatial
   convolution can be greatly sped up if it is done on separate tiles over
   the image (on multiple threads). So as input, you can either give tile
   values or one full array. Just note that if you give a single array as
   input, the `next' element has to be `NULL'.*/
gal_data_t *
gal_convolve_spatial(gal_data_t *tiles, gal_data_t *kernel,
                     size_t numthreads, int edgecorrection, int convoverch)
{
  char *name;
  struct spatial_params params;
  gal_data_t *out, *block=gal_tile_block(tiles);

  /* Small sanity checks. */
  if(tiles->ndim!=kernel->ndim)
    error(EXIT_FAILURE, 0, "The number of dimensions between the kernel and "
          "input should be the same in `gal_convolve_spatial'");
  if( block->type!=GAL_TYPE_FLOAT32
      || kernel->type!=GAL_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "`gal_convolve_spatial' currently only works on "
          "`float32' type input and kernel");

  /* Allocate the output convolved dataset. */
  name = ( block->name
           ? gal_checkset_malloc_cat("CONVL_", block->name) : NULL );
  out=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, block->ndim, block->dsize,
                     block->wcs, 0, block->minmapsize, name,
                     block->unit, NULL);

  /* Set the pointers in the parameters structure. */
  params.out=out;
  params.tiles=tiles;
  params.kernel=kernel;
  params.convoverch=convoverch;
  params.edgecorrection=edgecorrection;

  /* Do the spatial convolution on threads. */
  gal_threads_spin_off(convolve_spatial_on_thread, &params,
                       gal_data_num_in_ll(tiles), numthreads);

  /* Clean up and return the output array. */
  if(name) free(name);
  return out;
}
