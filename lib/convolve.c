/*********************************************************************
Convolve -- Convolve a dataset with a given kernel.
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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/list.h>
#include <gnuastro/tile.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
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
  size_t             id;     /* ID of tile being operatred on.           */
  gal_data_t      *tile;     /* Tile this thread is working on.          */
  gal_data_t *i_overlap;     /* Overlap tile over input dataset.         */
  gal_data_t *k_overlap;     /* Overlap tile over kernel dataset.        */
  size_t *overlap_start;     /* Starting coordinate of kernel overlap.   */
  size_t  *kernel_start;     /* Kernel starting point.                   */
  size_t    *host_start;     /* Starting coordinate of host.             */
  size_t           *pix;     /* 2*ndim: starting and ending of tile,
                                Later, just the pixel being convolved.   */
  int           on_edge;     /* If the tile is on the edge or not.       */
  gal_data_t      *host;     /* Size of host (channel or block).         */
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
   the necessary input image parameters are stored in 'overlap' (its
   'array' and 'dsize' elements).  */
static int
convolve_spatial_overlap(struct per_thread_spatial_prm *pprm, int tocorrect)
{
  struct spatial_params *cprm=pprm->cprm;
  gal_data_t *block=cprm->block, *kernel=cprm->kernel;
  size_t *dsize = tocorrect ? block->dsize : pprm->host->dsize;
  size_t ndim=block->ndim;

  size_t *kd=pprm->k_overlap->dsize;
  size_t *pp, *ppf, *hs, increment, size=1;
  size_t *h=dsize, *os=pprm->overlap_start;
  size_t *p, *pf, *od=pprm->i_overlap->dsize;
  size_t *k=kernel->dsize, *ks=pprm->kernel_start;
  int full_overlap=1, dim_full_overlap, is_start, is_end;


  /* In to-correct mode, the pix position needs to be relative to the
     block. */
  if(tocorrect)
    {
      hs=pprm->host_start;
      ppf=(pp=pprm->pix)+ndim; do *pp += *hs++; while(++pp<ppf);
    }


  /* Coordinate to start convolution for this pixel. */
  pf = (p=pprm->pix) + ndim;
  do
    {
      /* Initialize the overlap for this dimension (we'll assume it
         overlaps because this is the most common case usually). */
      dim_full_overlap=1;

      /* When the tile is on the edge, some pixels in it can have full
         overlap. So using the 'dim_full_overlap', we will do the same
         thing we do for the tiles that don't overlap for them. When
         'tocorrect!=0', then only pixels that are on the edge of the tile
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
                 start from pixel '5-2=3' in this dimension and the overlap
                 size should decrease by the same amount.

                    image:            0 1 2 3 4 5 6 7 8 9 ...
                    pixel:                p
                    kernel:     0 1 2 3 4 5 6 7 8 9 10

                 With the end: Similar to above, but assume the pixel is
                 two pixels away from the edge of a 100-pixel image. We are
                 no longer worried about the overlap or kernel starting
                 point, it is the width that we need to decrease it by:

                   97 + 5 - 100 + 1 : The '1' is because we want the pixel
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

              /* Put the overlap size into the kernel's overlap 'dsize'
                 also and then use it to update the total size of the
                 overlap. */
              *kd++ = *od;
              size *= *od;

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
          size *= *k;
          *kd++ = *od = *k;
          *os++ = *p - *k/2;

          /* Increment. */
          ++h;
          ++k;
          ++od;
        }
    }
  while(++p<pf);


  /* Update the 'size' element of both overlap datasets. */
  pprm->i_overlap->size = pprm->k_overlap->size = size;


  /* Make correction.

      Normal mode (when 'tocorrect==0'): add the host's starting location
         (necessary when convolution over the host/channel is treated
         independently). In this mode, until now we were working as if the
         the host/channel is the full image so the edges don't get
         mixed. But from now on we will be working over the allocated block
         to look at pixel values, so we need to convert the location to the
         proper place within the allocated array.

      To-correct mode: The boundaries were calculated with respect to the
         block, so we don't need to correct 'overlap_start'. But we need to
         correct the pixel position back to its original state (relative to
         the channel). */
  hs=pprm->host_start;
  if(tocorrect)
    { ppf=(pp=pprm->pix)           + ndim; do *pp -= *hs++; while(++pp<ppf); }
  else
    { ppf=(pp=pprm->overlap_start) + ndim; do *pp += *hs++; while(++pp<ppf); }


  /* Set the starting point of the dataset overlap tile. */
  increment=gal_dimension_coord_to_index(ndim, block->dsize,
                                         pprm->overlap_start);
  pprm->i_overlap->array=gal_pointer_increment(block->array, increment,
                                               block->type);


  /* Set the starting point of the kernel overlap tile. */
  increment = ( full_overlap
                ? 0
                : gal_dimension_coord_to_index(ndim,
                                               kernel->dsize,
                                               pprm->kernel_start) );
  pprm->k_overlap->array=gal_pointer_increment(kernel->array, increment,
                                               kernel->type);
  return full_overlap;
}





/* Convolve over one tile that is not touching the edge. */
static void
convolve_spatial_tile(struct per_thread_spatial_prm *pprm)
{
  gal_data_t *tile=pprm->tile;

  int full_overlap;
  double sum, ksum;
  struct spatial_params *cprm=pprm->cprm;
  gal_data_t *block=cprm->block, *kernel=cprm->kernel;
  size_t j, ndim=block->ndim, csize=tile->dsize[ndim-1];
  gal_data_t *i_overlap=pprm->i_overlap, *k_overlap=pprm->k_overlap;

  /* Variables for scanning a tile ('i_*') and the region around every
     pixel of a tile ('o_*'). */
  size_t start_fastdim;
  size_t i_inc, i_ninc, i_st_en[2];

  /* These variables depend on the type of the input. */
  float *i_start;
  float *in_v, *in=block->array, *out=cprm->out->array;


  /* Starting pixel for the host of this tile. Note that when we are in
     'convoverch' mode, 'host' refers to the fully allocated block of
     memory. */
  pprm->host=cprm->convoverch ? block : tile->block;
  gal_tile_start_coord(pprm->host, pprm->host_start);


  /* Set the starting and ending coordinates of this tile (recall that the
     space for the start and end coordinates is stored in 'p->pix'). When
     'convoverch' is set, we want to convolve over the whole allocated
     block, not just one channel. So in effect, it is the same as
     'rel_block' in 'gal_tile_start_end_coord'. */
  gal_tile_start_end_coord(tile, pprm->pix, cprm->convoverch);
  start_fastdim = pprm->pix[ndim-1];


  /* See if this tile is on the edge or not. */
  pprm->on_edge=convolve_tile_is_on_edge(pprm->host->dsize, pprm->pix,
                                         kernel->dsize, ndim);


  /* If it isn't on the edge and we are correcting an already convolved
     image ('tocorrect!=NULL'), then this tile can be ignored. */
  if(cprm->tocorrect && pprm->on_edge==0) return;


  /* Parse over all the tile elements. */
  i_inc=0; i_ninc=1;
  i_start=gal_tile_start_end_ind_inclusive(tile, block, i_st_en);
  while( i_st_en[0] + i_inc <= i_st_en[1] )
    {
      /* Initialize the value along the fastest dimension (it is not
         incremented during 'gal_tile_block_increment'). */
      pprm->pix[ndim-1]=start_fastdim;

      /* Go over each pixel to convolve. */
      for(j=0;j<csize;++j)
        {
          /* Pointer to the pixel under consideration. */
          in_v = i_start + i_inc + j;

          /* If the input on this pixel is a NaN, then just set the output
             to NaN too and go onto the next pixel. 'in_v' is the pointer
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

                  /* Initialize the necessary values. */
                  sum  = 0.0L;
                  ksum = cprm->edgecorrection ? 0.0L : 1.0L;

                  /* Parse over both the overlap tiles. */
                  GAL_TILE_PO_OISET(float, float, i_overlap, k_overlap, 1, 0, {
                      if( !isnan(*i) )
                        {
                          sum += *i * *o;
                          if(cprm->edgecorrection) ksum += *o;
                        }
                    });

                  /* Set the output value. */
                  out[ in_v - in ] = ksum==0.0L ? NAN : sum/ksum;
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
  if(pprm->id==2053)
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
  size_t ndim=block->ndim;
  struct per_thread_spatial_prm *pprm=&cprm->pprm[tprm->id];
  size_t *dsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                     "dsize");


  /* Set all dsize values to 1 (the values within 'overlap->dsize' will be
     changed during convolution). */
  for(i=0;i<ndim;++i) dsize[i]=1;


  /* Initialize/Allocate necessary items for this thread. */
  pprm->cprm          = cprm;
  pprm->pix           = gal_pointer_allocate(GAL_TYPE_SIZE_T, 2*ndim, 0,
                                             __func__, "pprm->pix");
  pprm->host_start    = gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0,
                                             __func__, "pprm->host_start");
  pprm->kernel_start  = gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0,
                                             __func__, "pprm->kernel_start");
  pprm->overlap_start = gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0,
                                             __func__, "pprm->overlap_start");
  pprm->i_overlap     = gal_data_alloc(NULL, block->type, ndim, dsize,
                                       NULL, 0, -1, 1, NULL, NULL, NULL);
  pprm->k_overlap     = gal_data_alloc(NULL, cprm->kernel->type, ndim, dsize,
                                       NULL, 0, -1, 1, NULL, NULL, NULL);
  free(dsize);
  free(pprm->i_overlap->array);
  free(pprm->k_overlap->array);
  pprm->i_overlap->block = cprm->block;
  pprm->k_overlap->block = cprm->kernel;


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* Set this tile's pointer into this thread's parameters. */
      pprm->id   = tprm->indexs[i];
      pprm->tile = &cprm->tiles[ pprm->id ];

      /* Do the convolution on this tile. */
      convolve_spatial_tile(pprm);
    }


  /* Clean up, wait until all other threads finish, then return. In a
     single thread situation, 'tprm->b==NULL'. */
  free(pprm->pix);
  free(pprm->host_start);
  free(pprm->kernel_start);
  free(pprm->overlap_start);
  gal_data_free(pprm->i_overlap);
  gal_data_free(pprm->k_overlap);
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* General spatial convolve function. This function is called by both
   'gal_convolve_spatial' and */
static gal_data_t *
gal_convolve_spatial_general(gal_data_t *tiles, gal_data_t *kernel,
                             size_t numthreads, int edgecorrection,
                             int convoverch, gal_data_t *tocorrect)
{
  struct spatial_params params;
  gal_data_t *out, *block=gal_tile_block(tiles);


  /* Small sanity checks. */
  if(tiles->ndim!=kernel->ndim)
    error(EXIT_FAILURE, 0, "%s: The number of dimensions between the kernel "
          "and input should be the same", __func__);
  if( block->type!=GAL_TYPE_FLOAT32 || kernel->type!=GAL_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "%s: only accepts 'float32' type input and "
          "kernel currently", __func__);

  /* It may happen that an input dataset is part of a linked list, but it
     is not actually a tile structure (the user wants to convolve the whole
     dataset without using tiles)! In that case, this function should break
     beacuse a linked list is interpretted as a tile structure here.*/
  if( tiles->block==NULL && tiles->next && tiles->next->block==NULL )
    error(EXIT_FAILURE, 0, "%s: the input is a linked list but not a "
          "tessellation (a list of tiles). This function is optimized to "
          "work on a list of tiles. Please (temporarily) set the 'next' "
          "element of the input to 'NULL' and call this function again",
          __func__);


  /* Set the output datastructure.  */
  if(tocorrect) out=tocorrect;
  else
    {
      /* Allocate the space for the convolved image. */
      out=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, block->ndim, block->dsize,
                         block->wcs, 0, block->minmapsize, block->quietmmap,
                         NULL, block->unit, NULL);

      /* Spatial convolution won't change the blank bit-flag, so use the
         block structure's blank bit flag. */
      out->flag = ( block->flag
                    | ( GAL_DATA_FLAG_BLANK_CH | GAL_DATA_FLAG_HASBLANK ) );
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
    error(EXIT_FAILURE, 0, "%s: %zu bytes for 'params.pprm'",
          __func__, numthreads * sizeof *params.pprm);


  /* Do the spatial convolution on threads. */
  gal_threads_spin_off(convolve_spatial_on_thread, &params,
                       gal_list_data_number(tiles), numthreads);


  /* Clean up and return the output array. */
  free(params.pprm);
  return out;
}





/* Convolve a dataset with a given kernel in the spatial domain. Spatial
   convolution can be greatly sped up if it is done on separate tiles over
   the image (on multiple threads). So as input, you can either give tile
   values or one full array. Just note that if you give a single array as
   input, the 'next' element has to be 'NULL'.*/
gal_data_t *
gal_convolve_spatial(gal_data_t *tiles, gal_data_t *kernel,
                     size_t numthreads, int edgecorrection, int convoverch)
{
  /* When there isn't any tile structure, 'convoverch' must be set to
     one. Recall that the input can be a single full dataset also. */
  if(tiles->block==NULL) convoverch=1;

  /* Call the general function. */
  return gal_convolve_spatial_general(tiles, kernel, numthreads,
                                      edgecorrection, convoverch, NULL);
}





/* Correct the edges of channels in an already convolved image when it was
   initially convolved with 'gal_convolve_spatial' with 'convoverch==0'. In
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
  if( gal_dimension_is_different(block, tocorrect) )
    error(EXIT_FAILURE, 0, "%s: the 'tocorrect' dataset has to have the "
          "same dimensions/size as the block of the 'tiles' input", __func__);
  if( block->type != tocorrect->type )
    error(EXIT_FAILURE, 0, "%s: the 'tocorrect' dataset has to have the same "
          "type as the block of the 'tiles' input. The given types are '%s' "
          "and '%s' respectively", __func__,
          gal_type_name(tocorrect->type, 1), gal_type_name(block->type, 1));

  /* Call the general function, which will do the correction. */
  gal_convolve_spatial_general(tiles, kernel, numthreads,
                               edgecorrection, 0, tocorrect);
}
