/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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





static size_t *
parse_spectrum_pepare(struct mkcatalog_passparams *pp, size_t *start_end_inc,
                      int32_t **st_o, float **st_v, float **st_std)
{
  size_t *tsize;
  gal_data_t *spectile;
  struct mkcatalogparams *p=pp->p;
  size_t coord[3], minmax[6], numslices=p->objects->dsize[0];
  gal_data_t *area, *sum, *esum, *proj, *eproj, *oarea, *osum, *oesum;

  /* Get the coordinates of the spectral tile's starting element, then make
     the tile. */
  gal_dimension_index_to_coord(gal_pointer_num_between(p->objects->array,
                                                       pp->tile->array,
                                                       p->objects->type),
                               p->objects->ndim, p->objects->dsize, coord);
  minmax[0]=0;                                   /* Changed to first slice.*/
  minmax[1]=coord[1];
  minmax[2]=coord[2];
  minmax[3]=p->objects->dsize[0]-1;              /* Changed to last slice. */
  minmax[4]=coord[1]+pp->tile->dsize[1]-1;
  minmax[5]=coord[2]+pp->tile->dsize[2]-1;
  spectile=gal_tile_series_from_minmax(p->objects, minmax, 1);

  /* Find the starting (and ending) pointers on each of the datasets. */
  *st_o   = gal_tile_start_end_ind_inclusive(spectile, p->objects,
                                             start_end_inc);
  *st_v   = (float *)(p->values->array) + start_end_inc[0];
  *st_std = ( p->std
                 ? ( p->std->size==p->objects->size
                     ? (float *)(p->std->array) + start_end_inc[0]
                     : NULL )
                 : NULL );

  /* Allocate the columns. */
  area  = gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &numslices, NULL, 1,
                         p->cp.minmapsize, p->cp.quietmmap, "AREA",
                         "counter", "Area of object in a slice");
  sum   = gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 1,
                         p->cp.minmapsize, p->cp.quietmmap, "SUM",
                         p->values->unit, "Sum of values with this label.");
  esum  = gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 1,
                         p->cp.minmapsize, p->cp.quietmmap, "SUM_ERR",
                         p->values->unit, "Error in SUM column.");
  proj  = gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 1,
                         p->cp.minmapsize, p->cp.quietmmap, "SUM_PROJECTED",
                         p->values->unit, "Sum of full projected 2D area on "
                         "a slice.");
  eproj = gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 1,
                         p->cp.minmapsize, p->cp.quietmmap, "SUM_PROJECTED_ERR",
                         p->values->unit, "Error in SUM_PROJECTED column.");
  oarea = gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &numslices, NULL, 1,
                         p->cp.minmapsize, p->cp.quietmmap, "AREA_OTHER",
                         "counter", "Area covered by other labels in a slice.");
  osum  = gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 1,
                         p->cp.minmapsize, p->cp.quietmmap, "SUM_OTHER",
                         p->values->unit, "Sum of values in other labels on "
                         "a slice.");
  oesum = gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 1,
                         p->cp.minmapsize, p->cp.quietmmap, "SUM_OTHER_ERR",
                         p->values->unit, "Error in SUM_OTHER column.");

  /* Fill up the contents of the first element (note that the first
     'gal_data_t' is actually in an array, so the skeleton is already
     allocated, we just have to allocate its contents. */
  gal_data_initialize(pp->spectrum, NULL, p->specsliceinfo->type, 1,
                      &numslices, NULL, 0, p->cp.minmapsize,
                      p->cp.quietmmap, NULL, NULL, NULL);
  gal_data_copy_to_allocated(p->specsliceinfo, pp->spectrum);
  pp->spectrum->next=gal_data_copy(p->specsliceinfo->next);


  /* Add all the other columns in the final spectrum table. */
  pp->spectrum->next->next                       = area;
  area->next                                     = sum;
  area->next->next                               = esum;
  area->next->next->next                         = proj;
  area->next->next->next->next                   = eproj;
  area->next->next->next->next->next             = oarea;
  area->next->next->next->next->next->next       = osum;
  area->next->next->next->next->next->next->next = oesum;

  /* Clean up and return. */
  tsize=spectile->dsize;
  spectile->dsize=NULL;
  gal_data_free(spectile);
  return tsize;
}





/* Since spectra will be a large table for many objects, it is very
   important to not consume too much space in for columns that don't need
   it. This function will check the integer columns and if they are smaller
   than the maximum values of smaller types, */
static void
parse_spectrum_uint32_to_best_type(gal_data_t **input)
{
  gal_data_t *tmp=gal_statistics_maximum(*input);

  /* If maximum is smaller than UINT8_MAX, convert it to uint8_t */
  if( *(uint32_t *)(tmp->array) < UINT8_MAX )
    *input = gal_data_copy_to_new_type_free(*input, GAL_TYPE_UINT8);

  /* otherwise, if it is smaller than UINT16_MAX, convert it to uint16_t */
  else if( *(uint32_t *)(tmp->array)      < UINT16_MAX )
    *input = gal_data_copy_to_new_type_free(*input, GAL_TYPE_UINT16);

  /* Clean up. */
  gal_data_free(tmp);
}





static void
parse_spectrum_end(struct mkcatalog_passparams *pp, gal_data_t *xybin)
{
  size_t i;
  double *searr, *pearr, *osearr;
  struct mkcatalogparams *p=pp->p;

  /* The datasets and their pointers. */
  gal_data_t *area  = pp->spectrum->next->next;
  gal_data_t *sum   = area->next;
  gal_data_t *esum  = area->next->next;
  gal_data_t *proj  = area->next->next->next;
  gal_data_t *eproj = area->next->next->next->next;
  gal_data_t *oarea = area->next->next->next->next->next;
  gal_data_t *osum  = area->next->next->next->next->next->next;
  gal_data_t *oesum = area->next->next->next->next->next->next->next;

  /* Apply corrections to the columns that need it. */
  searr  = esum->array;
  pearr  = eproj->array;
  osearr = oesum->array;
  for(i=0; i<p->objects->dsize[0]; ++i)
    {
      searr[i]  = sqrt( searr[i]  );
      pearr[i]  = sqrt( pearr[i]  );
      osearr[i] = sqrt( osearr[i] );
    }

  /* Convert the 'double' type columns to 'float'. The extra precision of
     'double' was necessary when we were summing values in each slice. But
     afterwards, it is not necessary at all (the measurement error is much
     larger than a double-precision floating point number (15
     decimals). But the extra space gained (double) is very useful in not
     wasting too much memory and hard-disk space or online transfer time.*/
  sum   = gal_data_copy_to_new_type_free(sum,   GAL_TYPE_FLOAT32);
  esum  = gal_data_copy_to_new_type_free(esum,  GAL_TYPE_FLOAT32);
  proj  = gal_data_copy_to_new_type_free(proj,  GAL_TYPE_FLOAT32);
  eproj = gal_data_copy_to_new_type_free(eproj, GAL_TYPE_FLOAT32);
  osum  = gal_data_copy_to_new_type_free(osum,  GAL_TYPE_FLOAT32);
  oesum = gal_data_copy_to_new_type_free(oesum, GAL_TYPE_FLOAT32);

  /* For the two area columns, find their maximum value and convert the
     dataset to the smallest type that can hold them. */
  parse_spectrum_uint32_to_best_type(&area);
  parse_spectrum_uint32_to_best_type(&oarea);

  /* List the datasets and write them into the pointer for this object
     (exact copy of the statement in 'parse_spectrum_pepare'). */
  pp->spectrum->next->next                       = area;
  area->next                                     = sum;
  area->next->next                               = esum;
  area->next->next->next                         = proj;
  area->next->next->next->next                   = eproj;
  area->next->next->next->next->next             = oarea;
  area->next->next->next->next->next->next       = osum;
  area->next->next->next->next->next->next->next = oesum;
}




/* Each spectrum is a multi-column table (note that the slice counter and
   wavelength are written in the end):

     Column 3:  Number of object pixels.
     Column 4:  Sum of object pixel values.
     Column 5:  Error in Column 2.
     Column 6:  Sum over all 2D projection over whole specturm.
     Column 7:  Error in Column 4.
     Column 8:  Area of other labels in this slice.
     Column 9:  Flux by other objects in projected area.
     Column 10: Error in Column 9.
 */
static void
parse_spectrum(struct mkcatalog_passparams *pp, gal_data_t *xybin)
{
  struct mkcatalogparams *p=pp->p;

  gal_data_t *area;
  float *st_v, *st_std;
  uint32_t *narr, *oarr;
  size_t nproj=0, *tsize, start_end_inc[2];
  uint8_t *xybinarr = xybin ? xybin->array : NULL;
  int32_t *O, *OO, *st_o, *objarr=p->objects->array;
  size_t tid, *dsize=p->objects->dsize, num_increment=1;
  double var, *sarr, *searr, *parr, *pearr, *osarr, *osearr;
  size_t increment=0, pind=0, sind=0, ndim=p->objects->ndim, c[3];
  float st, sval, *V=NULL, *ST=NULL, *std=p->std?p->std->array:NULL;

  /* Prepare the columns to write in. */
  tsize  = parse_spectrum_pepare(pp, start_end_inc, &st_o, &st_v, &st_std);
  area   = pp->spectrum->next->next;
  narr   = area->array;
  sarr   = area->next->array;
  searr  = area->next->next->array;
  parr   = area->next->next->next->array;
  pearr  = area->next->next->next->next->array;
  oarr   = area->next->next->next->next->next->array;
  osarr  = area->next->next->next->next->next->next->array;
  osearr = area->next->next->next->next->next->next->next->array;

  /* If tile-id isn't necessary, set 'tid' to a blank value. */
  tid = (p->std && p->std->size>1 && st_std == NULL) ? 0 : GAL_BLANK_SIZE_T;

  /* Parse each contiguous patch of memory covered by this object. */
  while( start_end_inc[0] + increment <= start_end_inc[1] )
    {
      /* Set the contiguous range to parse. The pixel-to-pixel counting
         along the fastest dimension will be done over the 'O' pointer. */
      if( p->values        ) V  = st_v   + increment;
      if( p->std && st_std ) ST = st_std + increment;
      OO = ( O = st_o + increment ) + pp->tile->dsize[ndim-1];

      /* Parse the tile. */
      do
        {
          /* Only continue if this voxel is useful: it isn't NaN, or its
             covered by the projected area or object's label. */
          if( !isnan(*V) && (xybin && xybinarr[pind]==2) )
            {
              /* Get the error in measuing this pixel's flux. */
              if(p->std)
                {
                  /* If the standard deviation is given on a tile
                     structure, estimate the tile ID. */
                  if(tid != GAL_BLANK_SIZE_T)
                    {
                      gal_dimension_index_to_coord(O-objarr, ndim, dsize, c);
                      tid=gal_tile_full_id_from_coord(&p->cp.tl, c);
                    }

                  /* Get the error associated with this voxel. Note that if
                     we are given a variance dataset already, there is no
                     need to use 'st*st', we can directly use 'sval'. */
                  sval = st_std ? *ST : (p->std->size>1?std[tid]:std[0]);
                  st = p->variance ? sqrt(sval) : sval;
                  var = (p->variance ? sval : st*st) + fabs(*V);
                }
              else var = NAN;


              /* Projected spectra: see if we have a value of '2' in the
                 'xybin' array (showing that there is atleast one non-blank
                 element there over the whole spectrum.  */
              ++nproj;
              parr [ sind ] += *V;
              pearr[ sind ] += var;

              /* Calculate the number of labeled/detected pixels that
                 don't belong to this object. */
              if(*O>0)
                {
                  if(*O==pp->object)
                    {
                      ++narr[ sind ];
                      sarr  [ sind ] += *V;
                      searr [ sind ] += var;
                    }
                  else
                    {
                      ++oarr [ sind ];
                      osarr  [ sind ] += *V;
                      osearr [ sind ] += var;
                    }
                }
            }

          /* Increment the pointers. */
          if( xybin            ) ++pind;
          if( p->values        ) ++V;
          if( p->std && st_std ) ++ST;
        }
      while(++O<OO);

      /* Increment to the next contiguous region of this tile. */
      increment += ( gal_tile_block_increment(p->objects, tsize,
                                              num_increment++, NULL) );

      /* Increment the slice number, 'sind', and reset the projection (2D)
         index 'pind' if we have just finished parsing a slice. */
      if( (num_increment-1)%pp->tile->dsize[1]==0 )
        {
          /* If there was no measurement, set NaN for the values and their
             errors (zero is meaningful). */
          if( nproj      ==0 ) parr[sind]  = pearr[sind]  = NAN;
          if( narr[sind] ==0 ) sarr[sind]  = searr[sind]  = NAN;
          if( oarr[sind] ==0 ) osarr[sind] = osearr[sind] = NAN;

          nproj=pind=0;
          ++sind;
        }
    }

  /* Finalize the spectrum generation and clean up. */
  parse_spectrum_end(pp, xybin);
  free(tsize);

  /* For a check.
  gal_table_write(pp->spectrum, NULL, GAL_TABLE_FORMAT_BFITS,
                  "spectrum.fits", "SPECTRUM", 0);
  */
}





void
parse_objects(struct mkcatalog_passparams *pp)
{
  uint8_t *oif=pp->p->oiflag;
  struct mkcatalogparams *p=pp->p;
  size_t ndim=p->objects->ndim, *dsize=p->objects->dsize;

  double *oi=pp->oi;
  gal_data_t *xybin=NULL;
  size_t *tsize=pp->tile->dsize;
  uint8_t *u, *uf, goodvalue, *xybinarr=NULL;
  size_t d, pind=0, increment=0, num_increment=1;
  int32_t *O, *OO, *C=NULL, *objarr=p->objects->array;
  float var, sval, varval, skyval, *V=NULL, *SK=NULL, *ST=NULL;
  float *std=p->std?p->std->array:NULL, *sky=p->sky?p->sky->array:NULL;

  /* If tile processing isn't necessary, set 'tid' to a blank value. */
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
                 || oif[ OCOL_GZ   ]
                 || oif[ OCOL_VX   ]
                 || oif[ OCOL_VY   ]
                 || oif[ OCOL_VZ   ]
                 || oif[ OCOL_C_GX ]
                 || oif[ OCOL_C_GY ]
                 || oif[ OCOL_C_GZ ]
                 || sc
                 /* When the sky and its STD are tiles, we'll also need
                    the coordinate to find which tile a pixel belongs
                    to. */
                 || tid==GAL_BLANK_SIZE_T )
               ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__, "c")
               : NULL );


  /* If an XY projection area is necessary, we'll need to allocate an array
     to keep the projected space. */
  if( p->spectrum
      || oif[ OCOL_NUMALLXY ]
      || oif[ OCOL_NUMXY    ] )
    {
      xybin=gal_data_alloc(NULL, GAL_TYPE_UINT8, 2, &tsize[1], NULL,
                           1, p->cp.minmapsize, p->cp.quietmmap,
                           NULL, NULL, NULL);
      xybinarr=xybin->array;
    }


  /* Parse each contiguous patch of memory covered by this object. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the contiguous range to parse. The pixel-to-pixel counting
         along the fastest dimension will be done over the 'O' pointer. */
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
              if(xybin) xybinarr[ pind ]=1;
              if(oif[ OCOL_NUMALL   ]) oi[ OCOL_NUMALL ]++;


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
                  if(oif[ OCOL_GX ]) oi[ OCOL_GX ] += c[ ndim-1 ]+1;
                  if(oif[ OCOL_GY ]) oi[ OCOL_GY ] += c[ ndim-2 ]+1;
                  if(oif[ OCOL_GZ ]) oi[ OCOL_GZ ] += c[ ndim-3 ]+1;
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
                      if(oif[ OCOL_C_GX ]) oi[ OCOL_C_GX ] += c[ ndim-1 ]+1;
                      if(oif[ OCOL_C_GY ]) oi[ OCOL_C_GY ] += c[ ndim-2 ]+1;
                      if(oif[ OCOL_C_GZ ]) oi[ OCOL_C_GZ ] += c[ ndim-3 ]+1;
                    }
                }


              /* Value related measurements. */
              goodvalue=0;
              if( p->values && !( p->hasblank && isnan(*V) ) )
                {
                  /* For the standard-deviation measurements later. */
                  goodvalue=1;

                  /* General flux summations. */
                  if(xybin) xybinarr[ pind ]=2;
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
                      if(oif[ OCOL_VX ]) oi[ OCOL_VX ] += *V*(c[ ndim-1 ]+1);
                      if(oif[ OCOL_VY ]) oi[ OCOL_VY ] += *V*(c[ ndim-2 ]+1);
                      if(oif[ OCOL_VZ ]) oi[ OCOL_VZ ] += *V*(c[ ndim-3 ]+1);
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
                            oi[   OCOL_C_VX ] += *V * (c[ ndim-1 ]+1);
                          if(oif[ OCOL_C_VY ])
                            oi[   OCOL_C_VY ] += *V * (c[ ndim-2 ]+1);
                          if(oif[ OCOL_C_VZ ])
                            oi[   OCOL_C_VZ ] += *V * (c[ ndim-3 ]+1);
                        }
                    }
                }


              /* Sky value based measurements. */
              if(p->sky && oif[ OCOL_SUMSKY ])
                {
                  skyval = ( pp->st_sky
                             ? (isnan(*SK)?0:*SK)               /* Full array  */
                             : ( p->sky->size>1
                                 ? (isnan(sky[tid])?0:sky[tid]) /* Tile        */
                                 : sky[0] ) );                  /* Single value*/
                  if(!isnan(skyval))
                    {
                      oi[ OCOL_NUMSKY  ]++;
                      oi[ OCOL_SUMSKY  ] += skyval;
                    }
                }


              /* Sky standard deviation based measurements.*/
              if(p->std)
                {
                  sval = pp->st_std ? *ST : (p->std->size>1?std[tid]:std[0]);
                  var = p->variance ? sval : sval*sval;
                  if(oif[ OCOL_SUMVAR ] && (!isnan(var)))
                    {
                      oi[ OCOL_NUMVAR  ]++;
                      oi[ OCOL_SUMVAR  ] += var;
                    }
                  /* For each pixel, we have a sky contribution to the
                     counts and the signal's contribution. The standard
                     deviation in the sky is simply 'sval', but the
                     standard deviation of the signal (independent of the
                     sky) is 'sqrt(*V)'. Therefore the total variance of
                     this pixel is the variance of the sky added with the
                     absolute value of its sky-subtracted flux. We use the
                     absolute value, because especially as the signal gets
                     noisy there will be negative values, and we don't want
                     them to decrease the variance. */
                  if(oif[ OCOL_SUM_VAR ] && goodvalue)
                    {
                      varval=p->variance ? var : sval;
                      if(!isnan(varval)) oi[ OCOL_SUM_VAR ] += varval + fabs(*V);
                    }
                }
            }

          /* Increment the other pointers. */
          if( xybin                ) ++pind;
          if( p->values            ) ++V;
          if( p->clumps            ) ++C;
          if( p->sky && pp->st_sky ) ++SK;
          if( p->std && pp->st_std ) ++ST;
        }
      while(++O<OO);

      /* Increment to the next contiguous region of this tile. */
      increment += ( gal_tile_block_increment(p->objects, tsize,
                                              num_increment++, NULL) );

      /* If a 2D projection is requested, see if we should initialize (set
         to zero) the projection-index ('pind') not. */
      if(xybin && (num_increment-1)%tsize[1]==0 )
        pind=0;
    }

  /* Write the projected area columns. */
  if(xybin)
    {
      /* Any non-zero pixel must be set for NUMALLXY. */
      uf=(u=xybin->array)+xybin->size;
      do
        if(*u)
          {
            if(oif[ OCOL_NUMALLXY ]          ) oi[ OCOL_NUMALLXY ]++;
            if(oif[ OCOL_NUMXY    ] && *u==2 ) oi[ OCOL_NUMXY    ]++;
          }
      while(++u<uf);

      /* For a check on the projected 2D areas.
      if(xybin && pp->object==2)
        {
          gal_fits_img_write(xybin, "xybin.fits", NULL, NULL);
          exit(0);
        }
      */
    }

  /* Generate the Spectrum. */
  if(p->spectrum)
    parse_spectrum(pp, xybin);

  /* Clean up. */
  if(c)     free(c);
  if(sc)    free(sc);
  if(xybin) gal_data_free(xybin);
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
  gal_data_t *xybin=NULL;
  size_t *tsize=pp->tile->dsize;
  int32_t *O, *OO, *C=NULL, nlab;
  uint8_t *u, *uf, goodvalue, *cif=p->ciflag;
  size_t nngb=gal_dimension_num_neighbors(ndim);
  size_t i, ii, d, pind=0, increment=0, num_increment=1;
  float var, sval, varval, skyval, *V=NULL, *SK=NULL, *ST=NULL;
  int32_t *objects=p->objects->array, *clumps=p->clumps->array;
  float *std=p->std?p->std->array:NULL, *sky=p->sky?p->sky->array:NULL;

  /* If tile processing isn't necessary, set 'tid' to a blank value. */
  size_t tid = ( ( (p->sky     && p->sky->size>1 && pp->st_sky == NULL )
                   || ( p->std && p->std->size>1 && pp->st_std == NULL ) )
                 ? 0 : GAL_BLANK_SIZE_T );

  /* Coordinate shift. */
  size_t *sc = ( pp->shift
                 ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                        "sc")
                 : NULL );

  /* If any coordinate columns are requested. */
  size_t *c = ( ( cif[    CCOL_GX ]
                  || cif[ CCOL_GY ]
                  || cif[ CCOL_GZ ]
                  || cif[ CCOL_VX ]
                  || cif[ CCOL_VY ]
                  || cif[ CCOL_VZ ]
                  || cif[ CCOL_MINX ]
                  || cif[ CCOL_MAXX ]
                  || cif[ CCOL_MINY ]
                  || cif[ CCOL_MAXY ]
                  || cif[ CCOL_MINZ ]
                  || cif[ CCOL_MAXZ ]
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

  /* If an XY projection area is requested, we'll need to allocate an array
     to keep the projected space.*/
  if( cif[    CCOL_NUMALLXY ]
      || cif[ CCOL_NUMXY    ] )
    {
      xybin=gal_data_array_calloc(pp->clumpsinobj);
      for(i=0;i<pp->clumpsinobj;++i)
        gal_data_initialize(&xybin[i], NULL, GAL_TYPE_UINT8, 2, &tsize[1],
                            NULL, 1, p->cp.minmapsize, p->cp.quietmmap,
                            NULL, NULL, NULL);
    }


  /* Parse each contiguous patch of memory covered by this object. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the contiguous range to parse. The pixel-to-pixel counting
         along the fastest dimension will be done over the 'O' pointer. */
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
                      || cif[ CCOL_MINY ] || cif[ CCOL_MAXY ]
                      || cif[ CCOL_MINZ ] || cif[ CCOL_MAXZ ] )
                    ci[ CCOL_NUMALL ]++;
                  if(cif[ CCOL_NUMALLXY ])
                    ((uint8_t *)(xybin[*C-1].array))[ pind ] = 1;

                  /* Raw-position related measurements. */
                  if(c)
                    {
                      /* Get "C" the coordinates of this point. */
                      gal_dimension_index_to_coord(O-objects, ndim, dsize, c);

                      /* Position extrema measurements. */
                      if(cif[ CCOL_MINX ])
                        ci[CCOL_MINX]=CMIN(CCOL_MINX, ndim-1);
                      if(cif[ CCOL_MAXX ])
                        ci[CCOL_MAXX]=CMAX(CCOL_MAXX, ndim-1);
                      if(cif[ CCOL_MINY ])
                        ci[CCOL_MINY]=CMIN(CCOL_MINY, ndim-2);
                      if(cif[ CCOL_MAXY ])
                        ci[CCOL_MAXY]=CMAX(CCOL_MAXY, ndim-2);
                      if(cif[ CCOL_MINZ ])
                        ci[CCOL_MINZ]=CMIN(CCOL_MINZ, ndim-3);
                      if(cif[ CCOL_MAXZ ])
                        ci[CCOL_MAXZ]=CMAX(CCOL_MAXZ, ndim-3);

                      /* If we need tile-ID, get the tile ID now. */
                      if(tid!=GAL_BLANK_SIZE_T)
                        tid=gal_tile_full_id_from_coord(&p->cp.tl, c);

                      /* General geometric (independent of pixel value)
                         calculations. */
                      if(cif[ CCOL_GX ]) ci[ CCOL_GX ] += c[ ndim-1 ]+1;
                      if(cif[ CCOL_GY ]) ci[ CCOL_GY ] += c[ ndim-2 ]+1;
                      if(cif[ CCOL_GZ ]) ci[ CCOL_GZ ] += c[ ndim-3 ]+1;
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

                  /* Value related measurements, see 'parse_objects' for
                     comments. */
                  goodvalue=0;
                  if( p->values && !( p->hasblank && isnan(*V) ) )
                    {
                      /* For the standard-deviation measurement. */
                      goodvalue=1;

                      /* Fill in the necessary information. */
                      if(cif[ CCOL_NUM   ]) ci[ CCOL_NUM ]++;
                      if(cif[ CCOL_SUM   ]) ci[ CCOL_SUM ] += *V;
                      if(cif[ CCOL_NUMXY ])
                        ((uint8_t *)(xybin[*C-1].array))[ pind ] = 2;
                      if( *V > 0.0f )
                        {
                          if(cif[ CCOL_NUMWHT ]) ci[ CCOL_NUMWHT ]++;
                          if(cif[ CCOL_SUMWHT ]) ci[ CCOL_SUMWHT ] += *V;
                          if(cif[ CCOL_VX ])
                            ci[   CCOL_VX ] += *V * (c[ ndim-1 ]+1);
                          if(cif[ CCOL_VY ])
                            ci[   CCOL_VY ] += *V * (c[ ndim-2 ]+1);
                          if(cif[ CCOL_VZ ])
                            ci[   CCOL_VZ ] += *V * (c[ ndim-3 ]+1);
                          if(pp->shift)
                            {
                              ci[ CCOL_VXX ] += *V * sc[1] * sc[1];
                              ci[ CCOL_VYY ] += *V * sc[0] * sc[0];
                              ci[ CCOL_VXY ] += *V * sc[1] * sc[0];
                            }
                        }
                    }

                  /* Sky based measurements. */
                  if(p->sky && cif[ CCOL_SUMSKY ])
                    {
                      skyval = ( pp->st_sky
                                 ? *SK             /* Full. */
                                 : ( p->sky->size>1
                                     ? sky[tid]    /* Tile. */
                                     : sky[0] ) ); /* 1 value. */
                      if(!isnan(skyval))
                        {
                          ci[ CCOL_NUMSKY  ]++;
                          ci[ CCOL_SUMSKY  ] += skyval;
                        }
                    }

                  /* Sky Standard deviation based measurements, see
                     'parse_objects' for comments. */
                  if(p->std)
                    {
                      sval = ( pp->st_std
                               ? *ST
                               : (p->std->size>1 ? std[tid] : std[0]) );
                      var = p->variance ? sval : sval*sval;
                      if(cif[ CCOL_SUMVAR  ] && (!isnan(var)))
                        {
                          ci[ CCOL_NUMVAR ]++;
                          ci[ CCOL_SUMVAR ] += var;
                        }
                      if(cif[ CCOL_SUM_VAR ] && goodvalue)
                        {
                          varval=p->variance ? var : sval;
                          if(!isnan(varval))
                            ci[ CCOL_SUM_VAR ] += varval + fabs(*V);
                        }
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
                     region more than once, we use 'ngblabs' to keep track
                     of which label we have already added its value
                     to. 'ii' is the number of different labels this river
                     pixel has already been considered for. 'ngblabs' will
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
          if( xybin                ) ++pind;
          if( p->values            ) ++V;
          if( p->sky && pp->st_sky ) ++SK;
          if( p->std && pp->st_std ) ++ST;
        }
      while(++O<OO);

      /* Increment to the next contiguous region of this tile. */
      increment += ( gal_tile_block_increment(p->objects, tsize,
                                              num_increment++, NULL) );

      /* If a 2D projection is requested, see if we should initialize (set
         to zero) the projection-index ('pind') not. */
      if(xybin && (num_increment-1) % tsize[1]==0 )
        pind=0;
    }


  /* Write the projected area columns. */
  if(xybin)
    for(i=0;i<pp->clumpsinobj;++i)
      {
        /* Pointer to make things easier. */
        ci=&pp->ci[ i * CCOL_NUMCOLS ];

        /* Any non-zero pixel must be set for NUMALLXY. */
        uf=(u=xybin[i].array)+xybin[i].size;
        do
          if(*u)
            {
              if(cif[ CCOL_NUMALLXY ]          ) ci[ CCOL_NUMALLXY ]++;
              if(cif[ CCOL_NUMXY    ] && *u==2 ) ci[ CCOL_NUMXY    ]++;
            }
        while(++u<uf);

        /* For a check on the projected 2D areas. */
        if(xybin && pp->object==2)
          gal_fits_img_write(&xybin[i], "xybin.fits", NULL, NULL);

      }


  /* Clean up. */
  if(c)       free(c);
  if(sc)      free(sc);
  if(dinc)    free(dinc);
  if(ngblabs) free(ngblabs);
  if(xybin)   gal_data_array_free(xybin, pp->clumpsinobj, 1);
}





void
parse_order_based(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;

  float *V;
  double *ci;
  float *sigcliparr;
  gal_data_t *result;
  int32_t *O, *OO, *C=NULL;
  gal_data_t **clumpsvals=NULL;
  size_t i, increment=0, num_increment=1;
  size_t *tsize=pp->tile->dsize, ndim=p->objects->ndim;
  size_t counter=0, *ccounter=NULL, tmpsize=pp->oi[OCOL_NUM];
  gal_data_t *objvals=gal_data_alloc(NULL, p->values->type, 1, &tmpsize, NULL,
                                     0, p->cp.minmapsize, p->cp.quietmmap,
                                     NULL, NULL, NULL);


  /* Allocate space for the clump values. */
  if(p->clumps)
    {
      errno=0;
      clumpsvals=malloc(pp->clumpsinobj * sizeof *clumpsvals);
      if(clumpsvals==NULL)
        error(EXIT_FAILURE, errno, "%s: couldn't allocate 'clumpsvals' for "
              "%zu clumps", __func__, pp->clumpsinobj);


      /* Allocate the array necessary to keep the values of each clump. */
      ccounter=gal_pointer_allocate(GAL_TYPE_SIZE_T, pp->clumpsinobj, 1,
                                    __func__, "ccounter");
      for(i=0;i<pp->clumpsinobj;++i)
        {
          tmpsize=pp->ci[ i * CCOL_NUMCOLS + CCOL_NUM ];
          clumpsvals[i]=gal_data_alloc(NULL, p->values->type, 1, &tmpsize,
                                       NULL, 0, p->cp.minmapsize,
                                       p->cp.quietmmap, NULL, NULL, NULL);
        }
    }


  /* Parse each contiguous patch of memory covered by this object. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the contiguous range to parse. The pixel-to-pixel counting
         along the fastest dimension will be done over the 'O' pointer. */
      V = pp->st_v + increment;
      if(p->clumps) C = pp->st_c + increment;
      OO = ( O = pp->st_o + increment ) + tsize[ndim-1];

      /* Parse the next contiguous region of this tile. */
      do
        {
          /* If this pixel belongs to the requested object, then do the
             processing. 'hasblank' is constant, so when the values doesn't
             have any blank values, the 'isnan' will never be checked. */
          if( *O==pp->object && !( p->hasblank && isnan(*V) ) )
            {
              /* Copy the value for the whole object. */
              memcpy( gal_pointer_increment(objvals->array, counter++,
                                             p->values->type), V,
                      gal_type_sizeof(p->values->type) );

              /* We are also on a clump. */
              if(p->clumps && *C>0)
                memcpy( gal_pointer_increment(clumpsvals[*C-1]->array,
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


  /* Calculate the necessary values for the objects. */
  if(p->oiflag[ OCOL_MEDIAN ])
    {
      result=gal_data_copy_to_new_type_free(gal_statistics_median(objvals, 1),
                                            GAL_TYPE_FLOAT64);
      pp->oi[OCOL_MEDIAN]=*((double *)(result->array));
      gal_data_free(result);
    }
  if(p->oiflag[ OCOL_SIGCLIPNUM ]
     || p->oiflag[ OCOL_SIGCLIPSTD ]
     || p->oiflag[ OCOL_SIGCLIPMEAN ]
     || p->oiflag[ OCOL_SIGCLIPMEDIAN ])
    {
      /* Calculate the sigma-clipped results and write them in any
         requested column. */
      result=gal_statistics_sigma_clip(objvals, p->sigmaclip[0],
                                       p->sigmaclip[1], 1, 1);
      sigcliparr=result->array;
      if(p->oiflag[ OCOL_SIGCLIPNUM ])
        pp->oi[OCOL_SIGCLIPNUM]=sigcliparr[0];
      if(p->oiflag[ OCOL_SIGCLIPSTD ])
        pp->oi[OCOL_SIGCLIPSTD]=sigcliparr[3];
      if(p->oiflag[ OCOL_SIGCLIPMEAN ])
        pp->oi[OCOL_SIGCLIPMEAN]=sigcliparr[2];
      if(p->oiflag[ OCOL_SIGCLIPMEDIAN ])
        pp->oi[OCOL_SIGCLIPMEDIAN]=sigcliparr[1];

      /* Clean up the sigma-clipped values. */
      gal_data_free(result);
    }


  /* Clean up the object values. */
  gal_data_free(objvals);


  /* Calculate the necessary value for clumps. */
  if(p->clumps)
    {
      for(i=0;i<pp->clumpsinobj;++i)
        {
          /* Set the main row to fill. */
          ci=&pp->ci[ i * CCOL_NUMCOLS ];

          /* Do the necessary calculation. */
          if(p->ciflag[ CCOL_MEDIAN ])
            {
              result=gal_statistics_median(clumpsvals[i], 1);
              result=gal_data_copy_to_new_type_free(result, GAL_TYPE_FLOAT64);
              ci[ CCOL_MEDIAN ] = ( *((double *)(result->array))
                                    - (ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]) );
              gal_data_free(result);
            }
          if(p->ciflag[ CCOL_SIGCLIPNUM ]
             || p->ciflag[ CCOL_SIGCLIPSTD ]
             || p->ciflag[ CCOL_SIGCLIPMEAN ]
             || p->ciflag[ CCOL_SIGCLIPMEDIAN ])
            {
              /* Calculate the sigma-clipped results and write them in any
                 requested column. */
              result=gal_statistics_sigma_clip(clumpsvals[i], p->sigmaclip[0],
                                               p->sigmaclip[1], 1, 1);
              sigcliparr=result->array;
              if(p->ciflag[ CCOL_SIGCLIPNUM ])
                pp->ci[CCOL_SIGCLIPNUM]=sigcliparr[0];
              if(p->ciflag[ CCOL_SIGCLIPSTD ])
                pp->ci[CCOL_SIGCLIPSTD]=( sigcliparr[3]
                                          - (ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]));
              if(p->ciflag[ CCOL_SIGCLIPMEAN ])
                pp->ci[CCOL_SIGCLIPMEAN]=( sigcliparr[2]
                                           - (ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]));
              if(p->ciflag[ CCOL_SIGCLIPMEDIAN ])
                pp->ci[CCOL_SIGCLIPMEDIAN]=( sigcliparr[1]
                                             - (ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]));

              /* Clean up the sigma-clipped values. */
              gal_data_free(result);
            }

          /* Clean up this clump's values. */
          gal_data_free(clumpsvals[i]);
        }
      free(clumpsvals);
      free(ccounter);
    }
}
