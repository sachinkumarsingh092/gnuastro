/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

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
#include <float.h>
#include <stdlib.h>

#include <gnuastro/tile.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include "main.h"
#include "mkcatalog.h"



/*********************************************************************/
/*******************       Tiles for clumps       ********************/
/*********************************************************************/
static gal_data_t *
upperlimit_make_clump_tiles(struct mkcatalog_passparams *pp)
{
  gal_data_t *input=pp->p->input;
  size_t ndim=input->ndim, *dsize=input->dsize;

  int32_t *O, *C;
  gal_data_t *tiles=NULL;
  float *I, *II, *start=input->array;
  size_t increment=0, num_increment=1;
  size_t i, d, *min, *max, width=2*ndim;
  size_t *coord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                      "coord");
  size_t *minmax=gal_data_malloc_array(GAL_TYPE_SIZE_T,
                                       width*pp->clumpsinobj, __func__,
                                       "minmax");

  /* Initialize the minimum and maximum position for each tile/clump. So,
     we'll initialize the minimum coordinates to the maximum possible
     `size_t' value (in `GAL_BLANK_SIZE_T') and the maximums to zero. */
  for(i=0;i<pp->clumpsinobj;++i)
    for(d=0;d<ndim;++d)
      {
        minmax[ i * width +        d ] = GAL_BLANK_SIZE_T; /* Minimum. */
        minmax[ i * width + ndim + d ] = 0;                /* Maximum. */
      }

  /* Parse over the object and get the clump's minimum and maximum. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the pointers for this tile. */
      I = pp->st_i + increment;
      O = pp->st_o + increment;
      C = pp->st_c + increment;

      /* Go over the contiguous region. */
      II = I + dsize[ndim-1];
      do
        {
          /* Only consider clumps. */
          if( *O==pp->object && *C>0 )
            {
              /* Get the coordinates of this pixel. */
              gal_dimension_index_to_coord(I-start, ndim, dsize, coord);

              /* Check to see if this coordinate is the smallest/largest
                 found so far for this label. Note that labels start from
                 1, while indexs here start from zero. */
              min = &minmax[ (*C-1) * width        ];
              max = &minmax[ (*C-1) * width + ndim ];
              for(d=0;d<ndim;++d)
                {
                  if( coord[d] < min[d] ) min[d] = coord[d];
                  if( coord[d] > max[d] ) max[d] = coord[d];
                }
            }

          /* Increment the other pointers. */
          ++O; ++C;
        }
      while(++I<II);

      /* Increment to the next contiguous region. */
      increment += ( gal_tile_block_increment(input, dsize, num_increment++,
                                              NULL) );
    }

  /* For a check.
  for(i=0;i<pp->clumpsinobj;++i)
    printf("%zu: (%zu, %zu) --> (%zu, %zu)\n", i+1, minmax[i*width],
           minmax[i*width+1], minmax[i*width+2], minmax[i*width+3]);
  */

  /* Make the tiles. */
  tiles=gal_tile_series_from_minmax(input, minmax, pp->clumpsinobj);

  /* Cleanup and return. */
  free(coord);
  free(minmax);
  return tiles;
}




















/*********************************************************************/
/*******************         For one tile         ********************/
/*********************************************************************/
static double
upperlimit_one_tile(struct mkcatalog_passparams *pp, gal_data_t *tile,
                    unsigned long seed, int32_t clumplab)
{
  struct mkcatalogparams *p=pp->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  void *tarray;
  double sum, out;
  int continueparse;
  gal_data_t *sigclip;
  uint8_t *M=NULL, *st_m=NULL;
  float *uparr=pp->up_vals->array;
  size_t increment, num_increment;
  float *I, *II, *SK, *st_i, *st_sky;
  size_t d, tcounter=0, counter=0, se_inc[2];
  int32_t *O, *oO, *st_o, *st_oo, *st_oc, *oC=NULL;
  size_t maxcount = p->upnum * MKCATALOG_UPPERLIMIT_STOP_MULTIP;
  size_t *rcoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                       "rcoord");

  /* Initializations. */
  tarray=tile->array;
  gsl_rng_set(pp->rng, seed);


  /* `se_inc' is just used temporarily, the important thing here is
     `st_oo'. */
  st_oo = ( clumplab
            ? gal_tile_start_end_ind_inclusive(tile, p->objects, se_inc)
            : pp->st_o );
  st_oc = clumplab ? (int32_t *)(p->clumps->array) + se_inc[0] : NULL;


  /* Continue measuring magnitudes randomly until we get the desired
     number. */
  while(tcounter<maxcount && counter<p->upnum)
    {
      /* Get the random coordinates, note that `gsl_rng_uniform_int'
         returns an inclusive value. */
      for(d=0;d<ndim;++d)
        rcoord[d] = gsl_rng_uniform_int(pp->rng, dsize[d]-tile->dsize[d]-1);

      /* Set the tile's new starting pointer. */
      tile->array = gal_data_ptr_increment(p->input->array,
                          gal_dimension_coord_to_index(ndim, dsize, rcoord),
                                           p->input->type);

      /* Starting and ending coordinates for this random position, note
         that in `pp' we have the starting and ending coordinates of the
         actual tile. */
      increment     = 0;
      num_increment = 1;
      continueparse = 1;
      sum           = 0.0f;

      /* Starting pointers for the random tile. */
      st_i   = gal_tile_start_end_ind_inclusive(tile, p->input, se_inc);
      st_o               = (int32_t *)(p->objects->array) + se_inc[0];
      st_sky             = (float   *)(p->sky->array)     + se_inc[0];
      if(p->upmask) st_m = (uint8_t *)(p->upmask->array)  + se_inc[0];


      /* Starting pointers for the original tile.*/

      /* Parse over this object/clump. */
      while( se_inc[0] + increment <= se_inc[1] )
        {
          /* Set the pointers. */
          I                = st_i   + increment;    /* Random tile.   */
          SK               = st_sky + increment;    /* Random tile.   */
          O                = st_o   + increment;    /* Random tile.   */
          if(st_m) M       = st_m   + increment;    /* Random tile.   */
          oO               = st_oo  + increment;    /* Original tile. */
          if(clumplab) oC  = st_oc  + increment;    /* Original tile. */


          /* Parse over this contiguous region, similar to the first and
             second pass functions. */
          II = I + tile->dsize[ndim-1];
          do
            {
              /* Only use pixels over this object/clump. */
              if( *oO==pp->object
                  && ( oC==NULL || clumplab==0 || *oC==clumplab ) )
                {
                  if( *O || (M && *M) || ( p->hasblank && isnan(*I) ) )
                    continueparse=0;
                  else
                    sum += *I-*SK;
                }

              /* Increment the other pointers. */
              ++SK; ++O; ++oO; if(oC) ++oC;
            }
          while(continueparse && ++I<II);


          /* Increment to the next contiguous region of this tile. */
          if(continueparse)
            increment += ( gal_tile_block_increment(p->input, dsize,
                                                    num_increment++, NULL) );
          else break;
        }

      /* Further processing is only necessary if this random tile
         actually covered the sky region. */
      if(continueparse) uparr[ counter++ ] = sum;

      /* Increment the total-counter. */
      ++tcounter;
    }

  /* Calculate the standard deviation of this distribution. */
  if(counter==p->upnum)
    {
      sigclip=gal_statistics_sigma_clip(pp->up_vals, p->upsigmaclip[0],
                                        p->upsigmaclip[1], 1, 1);
      out = ((float *)(sigclip->array))[3] * p->upnsigma;
    }
  else out=NAN;

  /* Reset the tile's array pointer, clean up and return. */
  tile->array=tarray;
  free(rcoord);
  return out;
}




















/*********************************************************************/
/*******************     High level funciton      ********************/
/*********************************************************************/
void
upperlimit_calculate(struct mkcatalog_passparams *pp)
{
  size_t i;
  double *ci;
  unsigned long seed;
  gal_data_t *clumptiles;
  struct mkcatalogparams *p=pp->p;

  /* First find the upper limit magnitude for this object. */
  pp->oi[OCOL_UPPERLIMIT_B] = upperlimit_one_tile(pp, pp->tile,
                                                  p->seed+pp->object, 0);

  /* If a clumps image is present (a clump catalog is requested( and this
     object has clumps, then find the upper limit magnitude for the clumps
     within this object. */
  if(p->clumps && pp->clumpsinobj)
    {
      /* Make tiles covering the clumps. */
      clumptiles=upperlimit_make_clump_tiles(pp);

      /* Go over all the clumps. The random number generator seed for each
         clump/object has to be unique, but also reproducible (given the
         intial seed and identical inputs). So we have defined it based on
         the total number of objects and clumps and this object and clump's
         IDs. */
      for(i=0;i<pp->clumpsinobj;++i)
        {
          ci=&pp->ci[ i * CCOL_NUMCOLS ];
          seed = p->seed + p->numobjects + p->numclumps * pp->object + i;
          ci[CCOL_UPPERLIMIT_B] = upperlimit_one_tile(pp, &clumptiles[i],
                                                          seed, i+1);
        }

      /* Clean up the clump tiles. */
      gal_data_array_free(clumptiles, pp->clumpsinobj, 0);
    }
}
