/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2019, Free Software Foundation, Inc.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include "main.h"
#include "mkcatalog.h"

#include "parse.h"





/* Both passes are going to need their starting pointers set, so we'll do
   that here. */
void
parse_initialize(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;

  size_t i, ndim=p->objects->ndim;
  size_t *start_end=pp->start_end_inc;

  /* Initialize the number of clumps in this object. */
  pp->clumpsinobj=0;


  /* Initialize the intermediate values to zero. */
  memset(pp->oi, 0, OCOL_NUMCOLS * sizeof *pp->oi);


  /* Set the shifts in every dimension to avoid round-off errors in large
     numbers for the non-linear calculations. We are using the first pixel
     of each object's tile as the shift parameter to keep the mean
     (average) reasonably near to the standard deviation. Otherwise, when
     the object is far out in the image (large x and y positions), then
     roundoff errors are going to decrease the accuracy of the second order
     calculations. */
  if(pp->shift)
    {
      /* Get the coordinates of the tile's starting point. */
      gal_dimension_index_to_coord( ( (float *)(pp->tile->array)
                                      - (float *)(pp->tile->block->array) ),
                                    ndim, p->objects->dsize, pp->shift);

      /* Change their counting to start from 1, not zero, since we will be
         using them as FITS coordinates. */
      for(i=0;i<ndim;++i) ++pp->shift[i];
    }


  /* Set the starting and ending indexs of this tile/object on all (the
     possible) input arrays. */
  pp->st_o   = gal_tile_start_end_ind_inclusive(pp->tile, p->objects,
                                                start_end);
  pp->st_c   = (p->clumps
                ? (int32_t *)(p->clumps->array) + start_end[0] : NULL);
  pp->st_v   = (p->values
                ? (float *)(p->values->array)   + start_end[0] : NULL);
  pp->st_sky = ( p->sky
                 ? ( p->sky->size==p->objects->size
                     ? (float *)(p->sky->array) + start_end[0]
                     : NULL )
                 : NULL);
  pp->st_std = ( p->std
                 ? ( p->std->size==p->objects->size
                     ? (float *)(p->std->array) + start_end[0]
                     : NULL )
                 : NULL );
}





void
parse_objects(struct mkcatalog_passparams *pp)
{
  uint8_t *oif=pp->p->oiflag;
  struct mkcatalogparams *p=pp->p;
  size_t ndim=p->objects->ndim, *dsize=p->objects->dsize;

  double *oi=pp->oi;
  size_t *tsize=pp->tile->dsize;
  size_t d, increment=0, num_increment=1;
  float var, sval, *V=NULL, *SK=NULL, *ST=NULL;
  int32_t *O, *OO, *C=NULL, *objarr=p->objects->array;
  float *std=p->std?p->std->array:NULL, *sky=p->sky?p->sky->array:NULL;

  /* If tile processing isn't necessary, set `tid' to a blank value. */
  size_t tid = ( ( (p->sky     && p->sky->size>1 && pp->st_sky == NULL )
                   || ( p->std && p->std->size>1 && pp->st_std == NULL ) )
                 ? 0 : GAL_BLANK_SIZE_T );

  /* Coordinate shift. */
  size_t *sc = ( pp->shift
                 ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                        "sc")
                 : NULL );

  /* If any coordinate columns are requested. */
  size_t *c = (
               /* Coordinate-related columns. */
               ( oif[    OCOL_GX   ]
                 || oif[ OCOL_GY   ]
                 || oif[ OCOL_VX   ]
                 || oif[ OCOL_VY   ]
                 || oif[ OCOL_C_GX ]
                 || oif[ OCOL_C_GY ]
                 || sc
                 /* When the sky and its STD are tiles, we'll also need
                    the coordinate to find which tile a pixel belongs
                    to. */
                 || tid==GAL_BLANK_SIZE_T )
               ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__, "c")
               : NULL );


  /* Parse each contiguous patch of memory covered by this object. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the contiguous range to parse. The pixel-to-pixel counting
         along the fastest dimension will be done over the `O' pointer. */
      if( p->clumps            ) C  = pp->st_c   + increment;
      if( p->values            ) V  = pp->st_v   + increment;
      if( p->sky && pp->st_sky ) SK = pp->st_sky + increment;
      if( p->std && pp->st_std ) ST = pp->st_std + increment;
      OO = ( O = pp->st_o + increment ) + tsize[ndim-1];

      /* Parse the tile. */
      do
        {
          /* If this pixel belongs to the requested object then do the
             processing.  */
          if( *O==pp->object )
            {
              /* INTERNAL: Get the number of clumps in this object: it is
                 the largest clump ID over each object. */
              if( p->clumps && *C>0 )
                pp->clumpsinobj = *C > pp->clumpsinobj ? *C : pp->clumpsinobj;


              /* Add to the area of this object. */
              if(oif[ OCOL_NUMALL ]) oi[ OCOL_NUMALL ]++;


              /* Geometric coordinate measurements. */
              if(c)
                {
                  /* Convert the index to coordinate. */
                  gal_dimension_index_to_coord(O-objarr, ndim, dsize, c);

                  /* If we need tile-ID, get the tile ID now. */
                  if(tid!=GAL_BLANK_SIZE_T)
                    tid=gal_tile_full_id_from_coord(&p->cp.tl, c);

                  /* Do the general geometric (independent of pixel value)
                     calculations. */
                  if(oif[ OCOL_GX ]) oi[ OCOL_GX ] += c[1]+1;
                  if(oif[ OCOL_GY ]) oi[ OCOL_GY ] += c[0]+1;
                  if(pp->shift)
                    {
                      /* Calculate the shifted coordinates for second order
                         calculations. The coordinate is incremented because
                         from now on, the positions are in the FITS standard
                         (starting from one).  */
                      for(d=0;d<ndim;++d) sc[d] = c[d] + 1 - pp->shift[d];

                      /* Include the shifted values, note that the second
                         order moments are never needed independently, they
                         are used together to find the ellipticity
                         parameters. */
                      oi[ OCOL_GXX ] += sc[1] * sc[1];
                      oi[ OCOL_GYY ] += sc[0] * sc[0];
                      oi[ OCOL_GXY ] += sc[1] * sc[0];
                    }
                  if(p->clumps && *C>0)
                    {
                      if(oif[ OCOL_C_NUMALL ]) oi[ OCOL_C_NUMALL ]++;
                      if(oif[ OCOL_C_GX     ]) oi[ OCOL_C_GX     ] += c[1]+1;
                      if(oif[ OCOL_C_GY     ]) oi[ OCOL_C_GY     ] += c[0]+1;
                    }
                }


              /* Value related measurements. */
              if( p->values && !( p->hasblank && isnan(*V) ) )
                {
                  /* General flux summations. */
                  if(oif[ OCOL_NUM    ]) oi[ OCOL_NUM     ]++;
                  if(oif[ OCOL_SUM    ]) oi[ OCOL_SUM     ] += *V;

                  /* Get the necessary clump information. */
                  if(p->clumps && *C>0)
                    {
                      if(oif[ OCOL_C_NUM ]) oi[ OCOL_C_NUM ]++;
                      if(oif[ OCOL_C_SUM ]) oi[ OCOL_C_SUM ] += *V;
                    }

                  /* For flux weighted centers, we can only use positive
                     values, so do those measurements here. */
                  if( *V > 0.0f )
                    {
                      if(oif[ OCOL_NUMWHT ]) oi[ OCOL_NUMWHT ]++;
                      if(oif[ OCOL_SUMWHT ]) oi[ OCOL_SUMWHT ] += *V;
                      if(oif[ OCOL_VX     ]) oi[ OCOL_VX     ] += *V*(c[1]+1);
                      if(oif[ OCOL_VY     ]) oi[ OCOL_VY     ] += *V*(c[0]+1);
                      if(pp->shift)
                        {
                          oi[ OCOL_VXX    ] += *V * sc[1] * sc[1];
                          oi[ OCOL_VYY    ] += *V * sc[0] * sc[0];
                          oi[ OCOL_VXY    ] += *V * sc[1] * sc[0];
                        }
                      if(p->clumps && *C>0)
                        {
                          if(oif[ OCOL_C_NUMWHT ]) oi[ OCOL_C_NUMWHT ]++;
                          if(oif[ OCOL_C_SUMWHT ]) oi[ OCOL_C_SUMWHT ] += *V;
                          if(oif[ OCOL_C_VX ])
                            oi[   OCOL_C_VX ] += *V * (c[1]+1);
                          if(oif[ OCOL_C_VY ])
                            oi[   OCOL_C_VY ] += *V * (c[0]+1);
                        }
                    }
                }


              /* Sky value based measurements. */
              if(p->sky)
                if(oif[ OCOL_SUMSKY ])
                  oi[ OCOL_SUMSKY  ] += ( pp->st_sky
                                          ? *SK             /* Full array   */
                                          : ( p->sky->size>1
                                              ? sky[tid]    /* Tile         */
                                              : sky[0] ) ); /* Single value */


              /* Sky standard deviation based measurements.*/
              if(p->std)
                {
                  sval = pp->st_std ? *ST : (p->std->size>1?std[tid]:std[0]);
                  var = p->variance ? sval : sval*sval;
                  if(oif[ OCOL_SUMVAR ]) oi[ OCOL_SUMVAR  ] += var;
                  /* For each pixel, we have a sky contribution to the
                     counts and the signal's contribution. The standard
                     deviation in the sky is simply `sval', but the
                     standard deviation of the signal (independent of the
                     sky) is `sqrt(*V)'. Therefore the total variance of
                     this pixel is the variance of the sky added with the
                     absolute value of its sky-subtracted flux. We use the
                     absolute value, because especially as the signal gets
                     noisy there will be negative values, and we don't want
                     them to decrease the variance. */
                  if(oif[ OCOL_SUM_VAR ])
                    oi[ OCOL_SUM_VAR ] += ( (p->variance ? var : sval)
                                            + fabs(*V) );
                }
            }

          /* Increment the other pointers. */
          if( p->values            ) ++V;
          if( p->clumps            ) ++C;
          if( p->sky && pp->st_sky ) ++SK;
          if( p->std && pp->st_std ) ++ST;
        }
      while(++O<OO);

      /* Increment to the next contiguous region of this tile. */
      increment += ( gal_tile_block_increment(p->objects, tsize,
                                              num_increment++, NULL) );
    }

  /* Clean up. */
  if(c)  free(c);
  if(sc) free(sc);
}





/* Macro to help in finding the minimum and maximum coordinates. */
#define CMIN(COL, DIM) ( ci[ CCOL_NUMALL ]==1.0f                        \
                         ? (c[ DIM ]+1)                                 \
                         : ( (c[ DIM ]+1) < ci[ COL ]                   \
                             ? (c[ DIM ]+1) : ci[ COL ] ) )
#define CMAX(COL, DIM) ( ci[ CCOL_NUMALL ]==1.0f                        \
                         ? (c[ DIM ]+1)                                 \
                         : ( (c[ DIM ]+1) > ci[ COL ]                   \
                             ? (c[ DIM ]+1) : ci[ COL ] ) )




/* Parse over the clumps within an object.  */
void
parse_clumps(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;
  size_t ndim=p->objects->ndim, *dsize=p->objects->dsize;

  double *ci, *cir;
  uint8_t *cif=p->ciflag;
  size_t *tsize=pp->tile->dsize;
  int32_t *O, *OO, *C=NULL, nlab;
  float var, sval, *V=NULL, *SK=NULL, *ST=NULL;
  size_t i, ii, d, increment=0, num_increment=1;
  size_t nngb=gal_dimension_num_neighbors(ndim);
  int32_t *objects=p->objects->array, *clumps=p->clumps->array;
  float *std=p->std?p->std->array:NULL, *sky=p->sky?p->sky->array:NULL;

  /* If tile processing isn't necessary, set `tid' to a blank value. */
  size_t tid = ( ( (p->sky     && p->sky->size>1 && pp->st_sky == NULL )
                   || ( p->std && p->std->size>1 && pp->st_std == NULL ) )
                 ? 0 : GAL_BLANK_SIZE_T );

  /* Coordinate shift. */
  size_t *sc = ( pp->shift
                 ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                        "sc")
                 : NULL );

  /* If any coordinate columns are requested. */
  size_t *c = ( ( cif[    CCOL_GX   ]
                  || cif[ CCOL_GY   ]
                  || cif[ CCOL_VX   ]
                  || cif[ CCOL_VY   ]
                  || cif[ CCOL_MINX ]
                  || cif[ CCOL_MAXX ]
                  || cif[ CCOL_MINY ]
                  || cif[ CCOL_MAXY ]
                  || sc
                  || tid==GAL_BLANK_SIZE_T )
                ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                       "c")
                : NULL );

  /* Preparations for neighbor parsing. */
  int32_t *ngblabs=( ( cif[    CCOL_RIV_NUM     ]
                       || cif[ CCOL_RIV_SUM     ]
                       || cif[ CCOL_RIV_SUM_VAR ] )
                     ? gal_pointer_allocate(GAL_TYPE_INT32, nngb, 0,
                                             __func__, "ngblabs")
                     : NULL );
  size_t *dinc = ngblabs ? gal_dimension_increment(ndim, dsize) : NULL;


  /* Parse each contiguous patch of memory covered by this object. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the contiguous range to parse. The pixel-to-pixel counting
         along the fastest dimension will be done over the `O' pointer. */
      C = pp->st_c + increment;
      if( p->values            ) V  = pp->st_v   + increment;
      if( p->sky && pp->st_sky ) SK = pp->st_sky + increment;
      if( p->std && pp->st_std ) ST = pp->st_std + increment;
      OO = ( O = pp->st_o + increment ) + tsize[ndim-1];

      /* Parse the tile */
      do
        {
          /* If this pixel belongs to the requested object then do the
             processing. */
          if( *O==pp->object )
            {
              /* We are on a clump. */
              if(p->clumps && *C>0)
                {
                  /* Pointer to make things easier. Note that the clump
                     labels start from 1, but the array indexs from 0.*/
                  ci=&pp->ci[ (*C-1) * CCOL_NUMCOLS ];

                  /* Add to the area of this object. */
                  if( cif[ CCOL_NUMALL ]
                      || cif[ CCOL_MINX ] || cif[ CCOL_MAXX ]
                      || cif[ CCOL_MINY ] || cif[ CCOL_MAXY ] )
                    ci[ CCOL_NUMALL ]++;

                  /* Raw-position related measurements. */
                  if(c)
                    {
                      /* Get "C" the coordinates of this point. */
                      gal_dimension_index_to_coord(O-objects, ndim, dsize, c);

                      /* Position extrema measurements. */
                      if(cif[ CCOL_MINX ]) ci[CCOL_MINX]=CMIN(CCOL_MINX, 1);
                      if(cif[ CCOL_MAXX ]) ci[CCOL_MAXX]=CMAX(CCOL_MAXX, 1);
                      if(cif[ CCOL_MINY ]) ci[CCOL_MINY]=CMIN(CCOL_MINY, 0);
                      if(cif[ CCOL_MAXY ]) ci[CCOL_MAXY]=CMAX(CCOL_MAXY, 0);

                      /* If we need tile-ID, get the tile ID now. */
                      if(tid!=GAL_BLANK_SIZE_T)
                        tid=gal_tile_full_id_from_coord(&p->cp.tl, c);

                      /* General geometric (independent of pixel value)
                         calculations. */
                      if(cif[ CCOL_GX ]) ci[ CCOL_GX ] += c[1]+1;
                      if(cif[ CCOL_GY ]) ci[ CCOL_GY ] += c[0]+1;
                      if(pp->shift)
                        {
                          /* Shifted coordinates for second order moments,
                             see explanations in the first pass.*/
                          for(d=0;d<ndim;++d) sc[d] = c[d] + 1 - pp->shift[d];

                          /* Raw second-order measurements. */
                          ci[ CCOL_GXX ] += sc[1] * sc[1];
                          ci[ CCOL_GYY ] += sc[0] * sc[0];
                          ci[ CCOL_GXY ] += sc[1] * sc[0];
                        }
                    }

                  /* Value related measurements, see `parse_objects' for
                     comments. */
                  if( p->values && !( p->hasblank && isnan(*V) ) )
                    {
                      /* Fill in the necessary information. */
                      if(cif[ CCOL_NUM ]) ci[ CCOL_NUM ]++;
                      if(cif[ CCOL_SUM ]) ci[ CCOL_SUM ] += *V;
                      if( *V > 0.0f )
                        {
                          if(cif[ CCOL_NUMWHT ]) ci[ CCOL_NUMWHT ]++;
                          if(cif[ CCOL_SUMWHT ]) ci[ CCOL_SUMWHT ] += *V;
                          if(cif[ CCOL_VX     ]) ci[ CCOL_VX  ]+=*V*(c[1]+1);
                          if(cif[ CCOL_VY     ]) ci[ CCOL_VY  ]+=*V*(c[0]+1);
                          if(pp->shift)
                            {
                              ci[ CCOL_VXX ] += *V * sc[1] * sc[1];
                              ci[ CCOL_VYY ] += *V * sc[0] * sc[0];
                              ci[ CCOL_VXY ] += *V * sc[1] * sc[0];
                            }
                        }
                    }

                  /* Sky based measurements. */
                  if(p->sky)
                    if(cif[ CCOL_SUMSKY ])
                      ci[ CCOL_SUMSKY  ] += ( pp->st_sky
                                              ? *SK             /* Full */
                                              : ( p->sky->size>1
                                                  ? sky[tid]    /* Tile */
                                                  : sky[0] ) ); /* 1 value */

                  /* Sky Standard deviation based measurements, see
                     `parse_objects' for comments. */
                  if(p->std)
                    {
                      sval = ( pp->st_std
                               ? *ST
                               : (p->std->size>1 ? std[tid] : std[0]) );
                      var = p->variance ? sval : sval*sval;
                      if(cif[ CCOL_SUMVAR  ]) ci[ CCOL_SUMVAR ] += var;
                      if(cif[ CCOL_SUM_VAR ])
                        ci[ CCOL_SUM_VAR ] += ( (p->variance ? var : sval)
                                                + fabs(*V) );
                    }
                }

              /* This pixel is on the diffuse region (and the object
                 actually has clumps). If any river-based measurements are
                 necessary check to see if it is touching a clump or not,
                 but only if this object actually has any clumps. */
              else if(ngblabs && pp->clumpsinobj)
                {
                  /* We are on a diffuse (possibly a river) pixel. So the
                     value of this pixel has to be added to any of the
                     clumps in touches. But since it might touch a labeled
                     region more than once, we use `ngblabs' to keep track
                     of which label we have already added its value
                     to. `ii' is the number of different labels this river
                     pixel has already been considered for. `ngblabs' will
                     keep the list labels. */
                  ii=0;
                  memset(ngblabs, 0, nngb*sizeof *ngblabs);

                  /* Go over the neighbors and see if this pixel is
                     touching a clump or not. */
                  GAL_DIMENSION_NEIGHBOR_OP(O-objects, ndim, dsize, ndim,
                                            dinc,
                     {
                       /* Neighbor's label (mainly for easy reading). */
                       nlab=clumps[nind];

                       /* We only want neighbors that are a clump and part
                          of this object and part of the same object. */
                       if( nlab>0 && objects[nind]==pp->object)
                         {
                           /* Go over all already checked labels and make
                              sure this clump hasn't already been
                              considered. */
                           for(i=0;i<ii;++i) if(ngblabs[i]==nlab) break;

                           /* It hasn't been considered yet: */
                           if(i==ii)
                             {
                               /* Make sure it won't be considered any
                                  more. */
                               ngblabs[ii++] = nlab;

                               /* To help in reading. */
                               cir=&pp->ci[ (nlab-1) * CCOL_NUMCOLS ];

                               /* Write in the necessary values. */
                               if(cif[ CCOL_RIV_NUM  ])
                                 cir[ CCOL_RIV_NUM ]++;

                               if(cif[ CCOL_RIV_SUM  ])
                                 cir[ CCOL_RIV_SUM ] += *V;

                               if(cif[ CCOL_RIV_SUM_VAR  ])
                                 {
                                   sval = ( pp->st_std
                                            ? *ST
                                            : ( p->std->size>1
                                                ? std[tid]
                                                : std[0] )     );
                                   cir[ CCOL_RIV_SUM_VAR ] += fabs(*V)
                                     + (p->variance ? sval : sval*sval);
                                 }
                             }
                         }
                     });
                }
            }

          /* Increment the other pointers. */
          ++C;
          if( p->values            ) ++V;
          if( p->sky && pp->st_sky ) ++SK;
          if( p->std && pp->st_std ) ++ST;
        }
      while(++O<OO);

      /* Increment to the next contiguous region of this tile. */
      increment += ( gal_tile_block_increment(p->objects, tsize,
                                              num_increment++, NULL) );
    }

  /* Clean up. */
  if(c)       free(c);
  if(sc)      free(sc);
  if(dinc)    free(dinc);
  if(ngblabs) free(ngblabs);
}





void
parse_median(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;

  float *V;
  double *ci;
  gal_data_t *median;
  int32_t *O, *OO, *C=NULL;
  gal_data_t **clumpsmed=NULL;
  size_t i, increment=0, num_increment=1;
  size_t *tsize=pp->tile->dsize, ndim=p->objects->ndim;
  size_t counter=0, *ccounter=NULL, tmpsize=pp->oi[OCOL_NUM];
  gal_data_t *objmed=gal_data_alloc(NULL, p->values->type, 1, &tmpsize, NULL,
                                    0, p->cp.minmapsize, NULL, NULL, NULL);

  /* Allocate space for the clump medians. */
  if(p->clumps)
    {
      errno=0;
      clumpsmed=malloc(pp->clumpsinobj * sizeof *clumpsmed);
      if(clumpsmed==NULL)
        error(EXIT_FAILURE, errno, "%s: couldn't allocate `clumpsmed' for "
              "%zu clumps", __func__, pp->clumpsinobj);


      /* Allocate the array necessary to keep the values of each clump. */
      ccounter=gal_pointer_allocate(GAL_TYPE_SIZE_T, pp->clumpsinobj, 1,
                                    __func__, "ccounter");
      for(i=0;i<pp->clumpsinobj;++i)
        {
          tmpsize=pp->ci[ i * CCOL_NUMCOLS + CCOL_NUM ];
          clumpsmed[i]=gal_data_alloc(NULL, p->values->type, 1, &tmpsize,
                                      NULL, 0, p->cp.minmapsize, NULL, NULL,
                                      NULL);
        }
    }


  /* Parse each contiguous patch of memory covered by this object. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the contiguous range to parse. The pixel-to-pixel counting
         along the fastest dimension will be done over the `O' pointer. */
      V = pp->st_v + increment;
      if(p->clumps) C = pp->st_c + increment;
      OO = ( O = pp->st_o + increment ) + tsize[ndim-1];

      /* Parse the next contiguous region of this tile. */
      do
        {
          /* If this pixel belongs to the requested object, then do the
             processing. `hasblank' is constant, so when the values doesn't
             have any blank values, the `isnan' will never be checked. */
          if( *O==pp->object && !( p->hasblank && isnan(*V) ) )
            {
              /* Copy the value for the whole object. */
              memcpy( gal_pointer_increment(objmed->array, counter++,
                                             p->values->type), V,
                      gal_type_sizeof(p->values->type) );

              /* We are also on a clump. */
              if(p->clumps && *C>0)
                memcpy( gal_pointer_increment(clumpsmed[*C-1]->array,
                                               ccounter[*C-1]++,
                                               p->values->type), V,
                        gal_type_sizeof(p->values->type) );
            }

          /* Increment the other pointers. */
          ++V;
          if(p->clumps) ++C;
        }
      while(++O<OO);

      /* Increment to the next contiguous region of this tile. */
      increment += ( gal_tile_block_increment(p->objects, tsize,
                                              num_increment++, NULL) );
    }


  /* Calculate the final medians for objects. */
  median=gal_data_copy_to_new_type_free(gal_statistics_median(objmed, 1),
                                        GAL_TYPE_FLOAT64);
  pp->oi[OCOL_MEDIAN]=*((double *)(median->array));
  gal_data_free(objmed);
  gal_data_free(median);


  /* Calculate the median for clumps. */
  if(p->clumps)
    {
      for(i=0;i<pp->clumpsinobj;++i)
        {
          ci=&pp->ci[ i * CCOL_NUMCOLS ];
          median=gal_statistics_median(clumpsmed[i], 1);
          median=gal_data_copy_to_new_type_free(median, GAL_TYPE_FLOAT64);
          ci[ CCOL_MEDIAN ] = ( *((double *)(median->array))
                                - (ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]) );
          gal_data_free(clumpsmed[i]);
          gal_data_free(median);
        }
      free(clumpsmed);
      free(ccounter);
    }
}
