/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2018, Free Software Foundation, Inc.

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

#include "ui.h"
#include "mkcatalog.h"



/*********************************************************************/
/*******************       Tiles for clumps       ********************/
/*********************************************************************/
static gal_data_t *
upperlimit_make_clump_tiles(struct mkcatalog_passparams *pp)
{
  gal_data_t *input=pp->p->input;
  size_t ndim=input->ndim, *tsize=pp->tile->dsize;

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

  /* Parse over the object and get the clump's minimum and maximum
     positions.*/
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the pointers for this tile. */
      I = pp->st_i + increment;
      O = pp->st_o + increment;
      C = pp->st_c + increment;

      /* Go over the contiguous region. */
      II = I + tsize[ndim-1];
      do
        {
          /* Only consider clumps. */
          if( *O==pp->object && *C>0 )
            {
              /* Get the coordinates of this pixel. */
              gal_dimension_index_to_coord(I-start, ndim, input->dsize,
                                           coord);

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
      increment += ( gal_tile_block_increment(input, tsize, num_increment++,
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
static void
upperlimit_random_range(struct mkcatalog_passparams *pp, gal_data_t *tile,
                        size_t *min, size_t *max, int32_t clumplab)
{
  struct mkcatalogparams *p=pp->p;
  size_t d, tstart, minext, maxext, coord[]={0,0};
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  /* Set the minimum and maximum acceptable value for the range.  */
  if(p->uprange)
    {
      tstart=gal_data_ptr_dist(tile->block->array, tile->array,
                               p->input->type);
      gal_dimension_index_to_coord(tstart, ndim, dsize, coord);
    }

  /* Go over the dimensions and set the range along each dimension. */
  for(d=0;d<ndim;++d)
    {
      /* If uprange is given and it is not zero, then use it, otherwise,
         just use the full possible range. */
      if( p->uprange && p->uprange[d] )
        {
          /* Set the minimum of the random range. Since `size_t' is always
             positive, to make sure the difference isn't negative, we need
             to convert them to integer first. */
          if( (int)coord[d] - ((int)p->uprange[d])/2 > 0 )
            {
              min[d] = coord[d]-p->uprange[d]/2;
              maxext = 0;
            }
          else
            {
              min[d] = 0;
              maxext = -1 * ((int)coord[d] - ((int)p->uprange[d])/2);
            }

          /* Set the maximum of the random range. */
          if( coord[d] + p->uprange[d]/2 < dsize[d] - tile->dsize[d] )
            {
              max[d] = coord[d] + p->uprange[d]/2;
              minext = 0;
            }
          else
            {
              max[d] = dsize[d] - tile->dsize[d] - 1;
              minext = ( (coord[d] + p->uprange[d]/2)
                         - (dsize[d] - tile->dsize[d]) );
            }

          /* `minadd' and `maxadd' were defined to account for the removed
             smaller range when an object is on the edge. Their role is to
             add to the other side of the range as much as possible when
             one side is decreased on an edge. */
          if(minext)
            min[d] = ((int)(min[d]) - (int)minext >= 0) ? (min[d]-minext) : 0;
          if(maxext)
            max[d] = ( (max[d] + maxext < dsize[d] - tile->dsize[d])
                       ? (max[d] + maxext) : (dsize[d]-tile->dsize[d]-1) );
        }
      else
        {
          min[d]=0;
          max[d]=dsize[d]-tile->dsize[d]-1;
        }

      /* A small sanity check. */
      if( max[d]-min[d] < 2*tile->dsize[d] )
        {
          if(clumplab)
            fprintf(stderr, "WARNING: object %d clump %d: range of random "
                    "positions (%zu) along dimension %zu for upper-limit "
                    "calculations is smaller than double of its size (%zu) "
                    "in this dimension.\n\n", pp->object, clumplab,
                    max[d]-min[d], ndim-d, 2*tile->dsize[d]);
          else
            fprintf(stderr, "WARNING: object %d: range of random "
                    "positions (%zu) along dimension %zu for upper-limit "
                    "calculations is smaller than double of its size (%zu) "
                    "in this dimension.\n\n", pp->object, max[d]-min[d],
                    ndim-d, 2*tile->dsize[d]);
        }
    }
}





/* Return a random position in the requested dimension. */
static size_t
upperlimit_random_position(struct mkcatalog_passparams *pp, gal_data_t *tile,
                           size_t dim, size_t *min, size_t *max)
{
  size_t r;
  struct mkcatalogparams *p=pp->p;

  /* `gsl_rng_get' returns an inclusive value between the minimum and
     maximum of the particular generator. It may happen that the labeled
     region extends the full range of a dimension. In that case, the only
     possible starting point would be 0. */
  if( (int)(p->input->dsize[dim]) - (int)(tile->dsize[dim]) > 0 )
    {
      r=gsl_rng_get(pp->rng); /* For easy reading. */
      return lrint( (float)(min[dim])
                    + ( (float)(r-p->rngmin)/(float)(p->rngdiff)
                        * (float)(max[dim] - min[dim]) ) );
    }
  else
    return 0;
}





/* Given the distribution of values, do the upper-limit calculations. */
static void
upperlimit_measure(struct mkcatalog_passparams *pp, int32_t clumplab,
                   int do_measurement)
{
  float *scarr;
  gal_data_t *column;
  size_t init_size, col, one=1;
  struct mkcatalogparams *p=pp->p;
  gal_data_t *sigclip=NULL, *sum, *qfunc;
  double *o = ( clumplab
                ? &pp->ci[ (clumplab-1) * CCOL_NUMCOLS ]
                : pp->oi );

  /* If the random distribution exsits, then fill it in. */
  if(do_measurement)
    {
      /* These columns are for both objects and clumps, so if they are
         requested in objects, they will also be written for clumps here
         (the order is irrelevant here). */
      for(column=p->objectcols; column!=NULL; column=column->next)
        {
          switch(column->status)
            {
            /* Columns that depend on the sigma of the distribution. */
            case UI_KEY_UPPERLIMIT:
            case UI_KEY_UPPERLIMITMAG:
            case UI_KEY_UPPERLIMITONESIGMA:

              /* We only need to do this once. */
              if(sigclip==NULL)
                {
                  /* Calculate the sigma-clipped standard deviation. Since
                     it is done in place, the size will change, so we'll
                     keep the size here and put it back after we are
                     done. */
                  init_size=pp->up_vals->size;
                  sigclip=gal_statistics_sigma_clip(pp->up_vals,
                                                    p->upsigmaclip[0],
                                                    p->upsigmaclip[1], 1, 1);
                  pp->up_vals->size=pp->up_vals->dsize[0]=init_size;
                  scarr=sigclip->array;

                  /* Write the raw sigma. */
                  col = clumplab ? CCOL_UPPERLIMIT_S : OCOL_UPPERLIMIT_S;
                  o[col] = scarr[3];

                  /* Write the multiple of `upnsigma'. */
                  col = clumplab ? CCOL_UPPERLIMIT_B : OCOL_UPPERLIMIT_B;
                  o[col] = scarr[3] * p->upnsigma;

                  /* Clean up. */
                  gal_data_free(sigclip);
                }
              break;

            /* Quantile column. */
            case UI_KEY_UPPERLIMITQUANTILE:

              /* Similar to the case for sigma-clipping, we'll need to keep
                 the size here also. */
              init_size=pp->up_vals->size;
              sum=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL, 0,
                                 -1, NULL, NULL, NULL);
              ((float *)(sum->array))[0]=o[clumplab ? CCOL_SUM : OCOL_SUM];
              qfunc=gal_statistics_quantile_function(pp->up_vals, sum, 1);

              /* Fill in the column. */
              col = clumplab ? CCOL_UPPERLIMIT_Q : OCOL_UPPERLIMIT_Q;
              pp->up_vals->size=pp->up_vals->dsize[0]=init_size;
              o[col] = ((double *)(qfunc->array))[0];
              break;
            }
        }
    }
  else
    {
      o[ clumplab ? CCOL_UPPERLIMIT_S : OCOL_UPPERLIMIT_S ] = NAN;
      o[ clumplab ? CCOL_UPPERLIMIT_B : OCOL_UPPERLIMIT_B ] = NAN;
      o[ clumplab ? CCOL_UPPERLIMIT_Q : OCOL_UPPERLIMIT_Q ] = NAN;
    }
}





static void
upperlimit_one_tile(struct mkcatalog_passparams *pp, gal_data_t *tile,
                    unsigned long seed, int32_t clumplab)
{
  struct mkcatalogparams *p=pp->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  double sum;
  void *tarray;
  int continueparse;
  uint8_t *M=NULL, *st_m=NULL;
  float *uparr=pp->up_vals->array;
  float *I, *II, *SK, *st_i, *st_sky;
  size_t d, tcounter=0, counter=0, se_inc[2];
  size_t min[2], max[2], increment, num_increment;
  int32_t *O, *oO, *st_o, *st_oo, *st_oc, *oC=NULL;
  size_t maxcount = p->upnum * MKCATALOG_UPPERLIMIT_STOP_MULTIP;
  size_t *rcoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                       "rcoord");

  /* Initializations. */
  tarray=tile->array;
  gsl_rng_set(pp->rng, seed);


  /* Set the range of random values for this tile. */
  upperlimit_random_range(pp, tile, min, max, clumplab);


  /* `se_inc' is just used temporarily, the important thing here is
     `st_oo'. */
  st_oo = ( clumplab
            ? gal_tile_start_end_ind_inclusive(tile, p->objects, se_inc)
            : pp->st_o );
  st_oc = clumplab ? (int32_t *)(p->clumps->array) + se_inc[0] : NULL;


  /* Continue measuring randomly until we get the desired total number. */
  while(tcounter<maxcount && counter<p->upnum)
    {
      /* Get the random coordinates. */
      for(d=0;d<ndim;++d)
        rcoord[d] = upperlimit_random_position(pp, tile, d, min, max);

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

  /* Do the measurement on the random distribution. */
  upperlimit_measure(pp, clumplab, counter==p->upnum);

  /* Reset the tile's array pointer, clean up and return. */
  tile->array=tarray;
  free(rcoord);
}




















/*********************************************************************/
/*******************     High level funciton      ********************/
/*********************************************************************/
void
upperlimit_calculate(struct mkcatalog_passparams *pp)
{
  size_t i;
  unsigned long seed;
  gal_data_t *clumptiles;
  struct mkcatalogparams *p=pp->p;

  /* First find the upper limit magnitude for this object. */
  upperlimit_one_tile(pp, pp->tile, p->seed+pp->object, 0);

  /* If a clumps image is present (a clump catalog is requested) and this
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
          seed = p->seed + p->numobjects + p->numclumps * pp->object + i;
          upperlimit_one_tile(pp, &clumptiles[i], seed, i+1);
        }

      /* Clean up the clump tiles. */
      gal_data_array_free(clumptiles, pp->clumpsinobj, 0);
    }
}
