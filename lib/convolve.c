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
/* Each thread needs one of these structures. */
struct per_thread_spatial_prm
{
  /* Internally stored/used values. */
  gal_data_t      *tile;     /* Tile this thread is working on.          */
  gal_data_t   *overlap;     /* Overlap pointer and starting point.      */
  size_t *overlap_start;     /* Starting coordinate of kernel overlap.   */
  size_t  *kernel_start;     /* Kernel starting point.                   */
  size_t    *host_start;     /* Starting coordinate of host.             */
  size_t           *pix;     /* 2*ndim: starting and ending of tile,
                                Later, just the pixel being convolved.   */
  int           on_edge;     /* If the tile is on the edge or not.       */
  gal_data_t      *host;     /* Size of host (channel or block).         */
  size_t    k_start_inc;     /* Increment for kernel.                    */
  struct spatial_params *cprm; /* Link to main structure for all threads.*/
};




/* This is for the full list. */
struct spatial_params
{
  /* Main input/output parameters. */
  gal_data_t       *out;     /* Output data structure.                   */
  gal_data_t     *tiles;     /* Tiles over the input image.              */
  gal_data_t     *block;     /* Pointer to block for this tile.          */
  gal_data_t    *kernel;     /* Kernel to convolve with input.           */
  gal_data_t *tocorrect;     /* (possible) convolved image to correct.   */
  int        convoverch;     /* Ignore channel edges in convolution.     */
  int    edgecorrection;     /* Correct convolution's edge effects.      */
  struct per_thread_spatial_prm *pprm; /* Array of per-thread parameters.*/
};




/* Define the overlap of the kernel and image over this part of the image,
   the necessary input image parameters are stored in `overlap' (its
   `array' and `dsize' elements).  */
static int
convolve_spatial_overlap(struct per_thread_spatial_prm *pprm,
                         int tocorrect)
{
  struct spatial_params *cprm=pprm->cprm;
  gal_data_t *block=cprm->block, *kernel=cprm->kernel;
  size_t *dsize = tocorrect ? block->dsize : pprm->host->dsize;

  size_t *pp, *ppf, *hs;
  size_t overlap_inc, ndim=block->ndim;
  size_t *h=dsize, *os=pprm->overlap_start;
  size_t *k=kernel->dsize, *ks=pprm->kernel_start;
  int full_overlap=1, dim_full_overlap, is_start, is_end;
  size_t *p=pprm->pix, *pf=pprm->pix+ndim, *od=pprm->overlap->dsize;


  /* In to-correct mode, the pix position needs to be relative to the
     block. */
  if(tocorrect)
    {
      hs=pprm->host_start;
      ppf=(pp=pprm->pix)+ndim; do *pp += *hs++; while(++pp<ppf);
    }


  /* Coordinate to start convolution for this pixel. */
  do
    {
      /* Initialize the overlap for this dimension (we'll assume it
         overlaps because this is the most common case usually). */
      dim_full_overlap=1;

      /* When the tile is on the edge, still some pixels in it can have
         full overlap. So using the `dim_full_overlap', we will do the same
         thing we do for the tiles that don't overlap for them. When
         `tocorrect!=0', then only pixels that are on the edge of the tile
         will get to this point, so it must always be checked. */
      if( tocorrect ? 1 : pprm->on_edge )
        {
          /* See if this pixel is on the start and/or end of the dimension
             relative to the kernel. */
          is_start = *p          <  *k/2;
          is_end   = (*p + *k/2) >= *h;
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


  /* To check, add an `int check' argument to the function.
  if(check)
    {
      printf("pix (within %s): %zu, %zu\n", tocorrect ? "full" : "channel",
             pprm->pix[0], pprm->pix[1]);
      printf("\tk/2: %zu, %zu\n", kernel->dsize[0]/2, kernel->dsize[1]/2);
      printf("\th: %zu, %zu\n", dsize[0], dsize[1]);
      printf("\toverlap_start: %zu, %zu\n", pprm->overlap_start[0],
             pprm->overlap_start[1]);
      printf("\toverlap->dsize: %zu, %zu\n", pprm->overlap->dsize[0],
             pprm->overlap->dsize[1]);
      printf("\tfulloverlap: %d\n", full_overlap);
    }
  */

  /* Set the increment to start working on the kernel. */
  pprm->k_start_inc = ( full_overlap
                        ? 0
                        : gal_dimension_coord_to_index(ndim,
                                                       kernel->dsize,
                                                       pprm->kernel_start) );

  /* Make correction.

      Normal mode (when `tocorrect==0'): add the host's starting location
         (necessary when convolution over the host/channel is treated
         independently). In this mode, until now we were working as if the
         the host/channel is the full image so the edges don't get
         mixed. But from now on we will be working over the allocated block
         to look at pixel values, so we need to convert the location to the
         proper place within the allocated array.

      To-correct mode: The boundaries were calculated with respect to the
         block, so we don't need to correct `overlap_start'. But we need to
         correct the pixel position back to its original state (relative to
         the channel). */
  hs=pprm->host_start;
  if(tocorrect)
    { ppf=(pp=pprm->pix)+ndim; do *pp -= *hs++; while(++pp<ppf); }
  else
    { ppf=(pp=pprm->overlap_start)+ndim; do *pp += *hs++; while(++pp<ppf); }


  /* Set the increment to start working on the overlap region and use that
     to set the starting pointer of the overlap region. */
  overlap_inc=gal_dimension_coord_to_index(ndim, block->dsize,
                                           pprm->overlap_start);
  pprm->overlap->array=gal_data_ptr_increment(block->array, overlap_inc,
                                              block->type);
  return full_overlap;
}





/* Convolve over one tile that is not touching the edge. */
static void
convolve_spatial_tile(struct per_thread_spatial_prm *pprm)
{

  double sum, ksum;
  int full_overlap;
  struct spatial_params *cprm=pprm->cprm;
  gal_data_t *tile=pprm->tile, *overlap=pprm->overlap;
  gal_data_t *block=cprm->block, *kernel=cprm->kernel;
  size_t i, ndim=block->ndim, csize=tile->dsize[ndim-1];

  /* Variables for scanning a tile (`i_*') and the region around every
     pixel of a tile (`o_*'). */
  size_t start_fastdim;
  size_t k_inc, i_inc, i_ninc, i_st_en[2], o_inc, o_ninc, o_st_en[2];

  /* These variables depend on the type of the input. */
  float *kv, *iv, *ivf, *i_start, *o_start, *k_start;
  float *in_v, *in=block->array, *out=cprm->out->array;


  /* Starting pixel for this tile. Note that when we are in `tocorrect'
     mode, this position has to be  */
  pprm->host=cprm->convoverch ? block : tile->block;
  gal_tile_start_coord(pprm->host, pprm->host_start);


  /* Set the starting and ending coordinates of this tile (recall that the
     start and end are the first two allocated spaces in
     parse_coords). When `convoverch' is set, we want to convolve over the
     whole allocated block, not just one channel. So in effect, it is the
     same as `rel_block' in `gal_tile_start_end_coord'. */
  gal_tile_start_end_coord(tile, pprm->pix, cprm->convoverch);
  start_fastdim = pprm->pix[ndim-1];


  /* See if this tile is on the edge or not. */
  pprm->on_edge=convolve_tile_is_on_edge(pprm->host->dsize, pprm->pix,
                                         kernel->dsize, ndim);


  /* If it isn't on the edge and we are correcting an already convolved
     image (`tocorrect!=NULL'), then this tile can be ignored. */
  if(cprm->tocorrect && pprm->on_edge==0) return;

  /*
  if(tile_ind==2053)
    {
      printf("\ntile %zu...\n", tile_ind);
      printf("\tpix: %zu, %zu\n", pprm->pix[0], pprm->pix[1]);
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
      pprm->pix[ndim-1]=start_fastdim;

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
              full_overlap=convolve_spatial_overlap(pprm, 0);

              /* If tocorrect has been given and we have full overlap, then
                 just ignore this pixel. */
              if( !(cprm->tocorrect && full_overlap) )
                {

                  /* If we are in correct mode, then re-calculate the
                     full-overlap and all the other necessary paramters as
                     if the channels didn't exist. */
                  if(cprm->tocorrect)
                    full_overlap=convolve_spatial_overlap(pprm, 1);

                  /* Set the starting pixel over the image (`o_start'). */
                  o_start=gal_tile_start_end_ind_inclusive(overlap, block,
                                                           o_st_en);

                  /* Set the starting kernel pixel. Note that
                     `kernel_array' is `void *' (pointer arithmetic is not
                     defined on it). So we will first put it in `k_start,
                     and then increment that. */
                  k_start=kernel->array; k_start += pprm->k_start_inc;

                  /* Go over the kernel-overlap region. */
                  ksum = cprm->edgecorrection ? 0.0f : 1.0f;
                  sum=0.0f; k_inc=0; o_inc=0; o_ninc=1; kv=k_start;
                  while( o_st_en[0] + o_inc <= o_st_en[1] )
                    {
                      /* Go over the contiguous region. When there is full
                         overlap, we don't need to calculate incrementation
                         over the kernel, it is always a single
                         incrementation. But when we have partial overlap,
                         we'll need to calculate a different
                         incrementation. */
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
                         region of memory over this tile. */
                      o_inc += gal_tile_block_increment(block, overlap->dsize,
                                                        o_ninc++, NULL);
                      if(full_overlap==0)
                        k_inc += gal_tile_block_increment(kernel,
                                                          overlap->dsize,
                                                          o_ninc-1, NULL);
                    }

                  /* Set the output value. */
                  out[ in_v - in ] = ( ksum==0.0f
                                       ? NAN
                                       : sum/ksum );
                }
            }

          /* Increment the last coordinate. */
          pprm->pix[ndim-1]++;
        }

      /* Increase the increment from the start of the tile for the next
         contiguous patch. */
      i_inc += gal_tile_block_increment(block, tile->dsize, i_ninc++,
                                        pprm->pix);
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
  gal_data_t *block=cprm->block;

  size_t i;
  struct per_thread_spatial_prm *pprm=&cprm->pprm[tprm->id];
  size_t ndim=block->ndim, *dsize=gal_data_malloc_array(GAL_TYPE_SIZE_T,ndim);


  /* Set all dsize values to 1 (the values within `overlap->dsize' will be
     changed during convolution). */
  for(i=0;i<ndim;++i) dsize[i]=1;


  /* Initialize/Allocate necessary items for this thread. */
  pprm->cprm           = cprm;
  pprm->pix            = gal_data_malloc_array(GAL_TYPE_SIZE_T, 2*ndim);
  pprm->host_start     = gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);
  pprm->kernel_start   = gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);
  pprm->overlap_start  = gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);
  pprm->overlap        = gal_data_alloc(NULL, block->type, ndim, dsize,
                                        NULL, 0, -1, NULL, NULL, NULL);
  free(dsize);
  free(pprm->overlap->array);
  pprm->overlap->block = cprm->block;


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    {
      /* Set this tile's pointer into this thread's parameters. */
      pprm->tile = &cprm->tiles[ tprm->indexs[i] ];

      /* Do the convolution on this tile. */
      convolve_spatial_tile(pprm);
    }


  /* Set the overlap dataset's array to NULL, it was used to point to
     different parts of the image during convolution. */
  pprm->overlap->array=NULL;


  /* Clean up, wait until all other threads finish, then return. In a
     single thread situation, `tprm->b==NULL'. */
  free(pprm->pix);
  free(pprm->host_start);
  free(pprm->kernel_start);
  free(pprm->overlap_start);
  gal_data_free(pprm->overlap);
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* General spatial convolve function. This function is called by both
   `gal_convolve_spatial' and */
static gal_data_t *
gal_convolve_spatial_general(gal_data_t *tiles, gal_data_t *kernel,
                             size_t numthreads, int edgecorrection,
                             int convoverch, gal_data_t *tocorrect)
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


  /* Set the output datastructure.  */
  if(tocorrect) out=tocorrect;
  else
    {
      name = ( block->name
               ? gal_checkset_malloc_cat("CONVL_", block->name) : NULL );
      out=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, block->ndim, block->dsize,
                         block->wcs, 0, block->minmapsize, name,
                         block->unit, NULL);
      if(name) free(name);
    }


  /* Set the pointers in the parameters structure. */
  params.out=out;
  params.tiles=tiles;
  params.block=block;
  params.kernel=kernel;
  params.tocorrect=tocorrect;
  params.convoverch=convoverch;
  params.edgecorrection=edgecorrection;


  /* Allocate the per-thread parameters. */
  errno=0;
  params.pprm=malloc(numthreads * sizeof *params.pprm);
  if(params.pprm==NULL)
    error(EXIT_FAILURE, 0, "%zu bytes for `params.p' in "
          "`gal_convolve_spatial_general'", numthreads * sizeof *params.pprm);


  /* Do the spatial convolution on threads. */
  gal_threads_spin_off(convolve_spatial_on_thread, &params,
                       gal_data_num_in_ll(tiles), numthreads);


  /* Clean up and return the output array. */
  free(params.pprm);
  return out;
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
  /* Call the general function. */
  return gal_convolve_spatial_general(tiles, kernel, numthreads,
                                      edgecorrection, convoverch, NULL);
}





/* Correct the edges of channels in an already convolved image when it was
   initially convolved with `gal_convolve_spatial' with `convoverch==0'. In
   that case, strong boundaries exist on the tile edges. So if you later
   need to remove those boundaries, you can call this function, it will
   only do convolution on the tiles that are near the edge, not the full
   area, so it is much faster. */
void
gal_convolve_spatial_correct_ch_edge(gal_data_t *tiles, gal_data_t *kernel,
                                     size_t numthreads, int edgecorrection,
                                     gal_data_t *tocorrect)
{
  gal_data_t *block=gal_tile_block(tiles);

  /* Some small sanity checks. */
  if( gal_data_dsize_is_different(block, tocorrect) )
    error(EXIT_FAILURE, 0, "the `tocorrect' dataset has to have the same "
          "dimensions/size as the block of the `tiles' input in "
          "`gal_convolve_spatial_correct_ch_edge'");
  if( block->type != tocorrect->type )
    error(EXIT_FAILURE, 0, "the `tocorrect' dataset has to have the same "
          "type as the block of the `tiles' input in "
          "`gal_convolve_spatial_correct_ch_edge'. The given types are `%s' "
          "and `%s' respectively", gal_type_to_string(tocorrect->type, 1),
          gal_type_to_string(block->type, 1));

  /* Call the general function, which will do the correction. */
  gal_convolve_spatial_general(tiles, kernel, numthreads,
                               edgecorrection, 0, tocorrect);
}
