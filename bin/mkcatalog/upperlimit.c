/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <inttypes.h>

#include <gnuastro/tile.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
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
  gal_data_t *objects=pp->p->objects;
  size_t ndim=objects->ndim, *tsize=pp->tile->dsize;

  gal_data_t *tiles=NULL;
  size_t increment=0, num_increment=1;
  size_t i, d, *min, *max, width=2*ndim;
  int32_t *O, *OO, *C, *start=objects->array;
  size_t *coord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                     "coord");
  size_t *minmax=gal_pointer_allocate(GAL_TYPE_SIZE_T,
                                      width*pp->clumpsinobj, 0, __func__,
                                      "minmax");

  /* Initialize the minimum and maximum position for each tile/clump. So,
     we'll initialize the minimum coordinates to the maximum possible
     'size_t' value (in 'GAL_BLANK_SIZE_T') and the maximums to zero. */
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
      C  = pp->st_c + increment;
      OO = ( O = pp->st_o + increment ) + tsize[ndim-1];

      /* Go over the contiguous region. */
      do
        {
          /* Only consider clumps. */
          if( *O==pp->object && *C>0 )
            {
              /* Get the coordinates of this pixel. */
              gal_dimension_index_to_coord(O-start, ndim, objects->dsize,
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
          ++C;
        }
      while(++O<OO);

      /* Increment to the next contiguous region. */
      increment += ( gal_tile_block_increment(objects, tsize, num_increment++,
                                              NULL) );
    }

  /* For a check.
  for(i=0;i<pp->clumpsinobj;++i)
    printf("%zu: (%zu, %zu) --> (%zu, %zu)\n", i+1, minmax[i*width],
           minmax[i*width+1], minmax[i*width+2], minmax[i*width+3]);
  */

  /* Make the tiles. */
  tiles=gal_tile_series_from_minmax(objects, minmax, pp->clumpsinobj);

  /* Cleanup and return. */
  free(coord);
  free(minmax);
  return tiles;
}




















/*********************************************************************/
/*******************         For one tile         ********************/
/*********************************************************************/
/* Set the minimum and maximum possible range to place the FIRST pixel of
   the object/clump tile over the dataset. */
static void
upperlimit_random_range(struct mkcatalog_passparams *pp, gal_data_t *tile,
                        size_t *min, size_t *max, int32_t clumplab)
{
  struct mkcatalogparams *p=pp->p;
  size_t d, tstart, minext, maxext, coord[]={0,0};
  size_t ndim=p->objects->ndim, *dsize=p->objects->dsize;

  /* Set the minimum and maximum acceptable value for the range.  */
  if(p->uprange)
    {
      tstart=gal_pointer_num_between(tile->block->array, tile->array,
                                     p->objects->type);
      gal_dimension_index_to_coord(tstart, ndim, dsize, coord);
    }

  /* Go over the dimensions and set the range along each dimension. */
  for(d=0;d<ndim;++d)
    {
      /* If uprange is given and it is not zero, then use it, otherwise,
         just use the full possible range. */
      if( p->uprange && p->uprange[d] )
        {
          /* Set the minimum of the random range. Since 'size_t' is always
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

          /* 'minadd' and 'maxadd' were defined to account for the removed
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
          /* We are positioning the FIRST pixel of the tile, not the
             center. So, the minimum possible value is zero, and in order
             to not push out of the image, the maximum is the
             'tile->dsize[d]' away from the edge. */
          min[d]=0;
          max[d]=dsize[d]-tile->dsize[d]-1;
        }

      /* A small warning to the user if the range isn't large enough. */
      if( max[d]-min[d] < 2*tile->dsize[d] )
        {
          p->uprangewarning=1;
          if(clumplab)
            fprintf(stderr, "WARNING-UPPERLIMIT: object %d clump %d, "
                    "dimension %zu: range (%zu) < 2*size (%zu).\n",
                    pp->object, clumplab, ndim-d, max[d]-min[d],
                    2*tile->dsize[d]);
          else
            fprintf(stderr, "WARNING-UPPERLIMIT: object %d, dimension %zu: "
                    "range (%zu) < 2*size (%zu).\n", pp->object, ndim-d,
                    max[d]-min[d], 2*tile->dsize[d]);
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

  /* 'gsl_rng_get' returns an inclusive value between the minimum and
     maximum of the particular generator. It may happen that the labeled
     region extends the full range of a dimension. In that case, the only
     possible starting point would be 0. */
  if( (int)(p->objects->dsize[dim]) - (int)(tile->dsize[dim]) > 0 )
    {
      r=gsl_rng_get(pp->rng); /* For easy reading. */
      return lrint( (float)(min[dim])
                    + ( (float)(r-p->rngmin)/(float)(p->rngdiff)
                        * (float)(max[dim] - min[dim]) ) );
    }
  else
    return 0;
}





/* It is necessary to write the upperlimit parameters into the output
   tables. The same set of information will thus be necessary both in the
   upperlimit check table and also the final output. This function will do
   the job in both cases.

   Note that in the check output, the sigma-clipping information is not
   used/necessary, so to avoid confusion, we won't write it.
*/
void
upperlimit_write_comments(struct mkcatalogparams *p,
                          gal_list_str_t **comments, int withsigclip)
{
  char *str;

  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      if(asprintf(&str, "--------- Upper-limit measurement ---------")<0)
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  if( asprintf(&str, "Number of usable random samples: %zu", p->upnum)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_str_add(comments, str, 0);

  if(p->uprange)
    {
      switch(p->objects->ndim)
        {
        case 2:
          if( asprintf(&str, "Range of random samples about target: "
                       "%zu, %zu", p->uprange[1], p->uprange[0])<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          break;
        case 3:
          if( asprintf(&str, "Range of random samples about target: %zu, "
                       "%zu, %zu", p->uprange[2], p->uprange[1],
                       p->uprange[0])<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                "address the problem. The value %zu is not recognized for "
                "'p->input->ndim'", __func__, PACKAGE_BUGREPORT,
                p->objects->ndim);
        }
      gal_list_str_add(comments, str, 0);
    }

  if( asprintf(&str, "Random number generator name: %s", p->rng_name)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_str_add(comments, str, 0);

  if( asprintf(&str, "Random number generator seed: %lu", p->rng_seed)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_str_add(comments, str, 0);

  if(withsigclip)
    {
      if( asprintf(&str, "Multiple of STD used for sigma-clipping: %.3f",
                   p->upsigmaclip[0])<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);

      if(p->upsigmaclip[1]>=1.0f)
        {
          if( asprintf(&str, "Number of clips for sigma-clipping: %.0f",
                       p->upsigmaclip[1])<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        {
          if( asprintf(&str, "Tolerance level to sigma-clipping: %.3f",
                       p->upsigmaclip[1])<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      gal_list_str_add(comments, str, 0);

      if( p->oiflag[ OCOL_UPPERLIMIT_B ] )
        {
          if( asprintf(&str, "Multiple of sigma-clipped STD for upper-limit: "
                       "%.3f", p->upnsigma)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_list_str_add(comments, str, 0);
        }
    }
}





/* Write the values into a table for the user */
static void
upperlimit_write_check(struct mkcatalogparams *p, gal_list_sizet_t *check_x,
                       gal_list_sizet_t *check_y, gal_list_sizet_t *check_z,
                       gal_list_f32_t *check_s)
{
  float *sarr;
  char *tmp=NULL, *tmp2=NULL;
  gal_list_str_t *comments=NULL;
  size_t *xarr, *yarr, *zarr=NULL, tnum, ttnum, num;
  gal_data_t *x=NULL, *y=NULL, *z=NULL, *s=NULL; /* To avoid warnings. */


  /* Convert the lists to an array. */
  xarr=gal_list_sizet_to_array(check_x, 1, &num);
  yarr=gal_list_sizet_to_array(check_y, 1, &tnum);
  if(check_z) zarr=gal_list_sizet_to_array(check_z, 1, &ttnum);
  if(tnum!=num || (check_z && ttnum!=num) )
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
          "problem. For some reason the size of the input lists don't "
          "match (%zu, %zu)", __func__, PACKAGE_BUGREPORT, tnum, num);
  sarr=gal_list_f32_to_array(check_s, 1, &tnum);
  if(tnum!=num)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
          "problem. For some reason the size of the input lists don't "
          "match (%zu, %zu)", __func__, PACKAGE_BUGREPORT, tnum, num);


  /* Put the arrays into a data container. */
  x=gal_data_alloc(xarr, GAL_TYPE_SIZE_T, 1, &num, NULL, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "RANDOM_X", "pixel",
                   "X-axis position of random footprint's first pixel.");
  y=gal_data_alloc(yarr, GAL_TYPE_SIZE_T, 1, &num, NULL, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "RANDOM_Y", "pixel",
                   "Y-axis position of random footprint's first pixel.");
  if(check_z)
    z=gal_data_alloc(zarr, GAL_TYPE_SIZE_T, 1, &num, NULL, 0,
                     p->cp.minmapsize, p->cp.quietmmap, "RANDOM_Z", "pixel",
                     "Z-axis position of random footprint's first pixel.");
  s=gal_data_alloc(sarr, GAL_TYPE_FLOAT32, 1, &num, NULL, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "RANDOM_SUM",
                   p->values->unit ? p->values->unit : MKCATALOG_NO_UNIT,
                   "Sum of pixel values over random footprint.");


  /* If 'size_t' isn't 32-bit on this system, then convert the unsigned
     64-bit values to 32-bit because the FITS table format doesn't
     recognize 64-bit integers.*/
  if( GAL_TYPE_SIZE_T != GAL_TYPE_UINT32 )
    {
      x=gal_data_copy_to_new_type_free(   x, GAL_TYPE_UINT32);
      y=gal_data_copy_to_new_type_free(   y, GAL_TYPE_UINT32);
      if(check_z)
        z=gal_data_copy_to_new_type_free( z, GAL_TYPE_UINT32);
    }


  /* Write exactly what object/clump this table is for. */
  if( p->checkuplim[1]!=GAL_BLANK_INT32 )
    if( asprintf(&tmp2, ", Clump %d", p->checkuplim[1]) <0 )
      error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  if( asprintf(&tmp, "Upperlimit distribution for Object %d%s",
               p->checkuplim[0],
               ( p->checkuplim[1]==GAL_BLANK_INT32
                 ? "" : tmp2) ) <0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_str_add(&comments, tmp, 0);
  if(tmp2) {free(tmp2); tmp2=NULL;}


  /* Write the basic info, and conclude the comments. */
  mkcatalog_write_inputs_in_comments(p, &comments, 0, 0);
  upperlimit_write_comments(p, &comments, 0);
  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      if( asprintf(&tmp, "--------- Table columns ---------")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, tmp, 0);
    }


  /* Define a list from the containers and write them into a table. */
  x->next=y;
  if(check_z) { y->next=z; z->next=s; }
  else        { y->next=s;            }
  gal_list_str_reverse(&comments);
  gal_table_write(x, comments, p->cp.tableformat, p->upcheckout,
                  "UPPERLIMIT_CHECK", 0);

  /* Inform the user. */
  if(!p->cp.quiet)
    printf("  - Upperlimit check table: %s\n", p->upcheckout);

  /* Clean up. */
  gal_data_free(x);
  gal_data_free(y);
  gal_data_free(s);
  if(check_z) gal_data_free(z);
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
  gal_data_t *sum, *qfunc=NULL, *sigclip=NULL;
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
            case UI_KEY_UPPERLIMITSIGMA:
            case UI_KEY_UPPERLIMITONESIGMA:

              /* We only need to do this once, but the columns can be
                 requested in any order. */
              if(sigclip==NULL)
                {
                  /* Calculate the sigma-clipped standard deviation. Since
                     it is done in place, the size will change, so we'll
                     keep the size here and put it back after we are
                     done. */
                  init_size=pp->up_vals->size;
                  sigclip=gal_statistics_sigma_clip(pp->up_vals,
                                                    p->upsigmaclip[0],
                                                    p->upsigmaclip[1],1,1);
                  pp->up_vals->size=pp->up_vals->dsize[0]=init_size;
                  scarr=sigclip->array;

                  /* 1-sigma. */
                  col = clumplab ? CCOL_UPPERLIMIT_S : OCOL_UPPERLIMIT_S;
                  o[col] = scarr[3];

                  /* sigma multiplied by 'upnsigma'. */
                  col = clumplab ? CCOL_UPPERLIMIT_B : OCOL_UPPERLIMIT_B;
                  o[col] = scarr[3] * p->upnsigma;

                  /* Nonparametric skewness [ (Mean-Median)/STD ]. */
                  col = clumplab?CCOL_UPPERLIMIT_SKEW:OCOL_UPPERLIMIT_SKEW;
                  o[col] = ( scarr[2] - scarr[1] ) / scarr[3];

                  /* Clean up. */
                  gal_data_free(sigclip);
                }
              break;

            /* Quantile column. */
            case UI_KEY_UPPERLIMITQUANTILE:

              /* Also only necessary once (if requested multiple times). */
              if(qfunc==NULL)
                {
                  /* Similar to the case for sigma-clipping, we'll need to
                     keep the size here also. */
                  init_size=pp->up_vals->size;
                  sum=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL, 0,
                                     -1, 1, NULL, NULL, NULL);
                  ((float *)(sum->array))[0]=o[clumplab?CCOL_SUM:OCOL_SUM];
                  qfunc=gal_statistics_quantile_function(pp->up_vals, sum, 1);

                  /* Fill in the column. */
                  col = clumplab ? CCOL_UPPERLIMIT_Q : OCOL_UPPERLIMIT_Q;
                  pp->up_vals->size=pp->up_vals->dsize[0]=init_size;
                  o[col] = ((double *)(qfunc->array))[0];

                  /* Clean up. */
                  gal_data_free(sum);
                  gal_data_free(qfunc);
                }
              break;
            }
        }
    }
  else
    {
      o[ clumplab ? CCOL_UPPERLIMIT_B : OCOL_UPPERLIMIT_B ] = NAN;
      o[ clumplab ? CCOL_UPPERLIMIT_S : OCOL_UPPERLIMIT_S ] = NAN;
      o[ clumplab ? CCOL_UPPERLIMIT_Q : OCOL_UPPERLIMIT_Q ] = NAN;
    }
}





static void
upperlimit_one_tile(struct mkcatalog_passparams *pp, gal_data_t *tile,
                    unsigned long seed, int32_t clumplab)
{
  struct mkcatalogparams *p=pp->p;
  size_t ndim=p->objects->ndim, *dsize=p->objects->dsize;

  double sum;
  void *tarray;
  uint8_t *M=NULL, *st_m=NULL;
  int continueparse, writecheck=0;
  struct gal_list_f32_t *check_s=NULL;
  size_t d, counter=0, se_inc[2], nfailed=0;
  float *V, *st_v, *uparr=pp->up_vals->array;
  size_t min[3], max[3], increment, num_increment;
  int32_t *O, *OO, *oO, *st_o, *st_oo, *st_oc, *oC=NULL;
  size_t maxfails = p->upnum * MKCATALOG_UPPERLIMIT_MAXFAILS_MULTIP;
  struct gal_list_sizet_t *check_x=NULL, *check_y=NULL, *check_z=NULL;
  size_t *rcoord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                      "rcoord");

  /* See if a check table must be created for this distribution. */
  if( p->checkuplim[0]==pp->object )
    {
      /* We are on a clump */
      if( clumplab )
        {
          if( p->checkuplim[1]==clumplab )
            writecheck=1;
        }
      else
        if( p->checkuplim[1]==GAL_BLANK_INT32 )
          writecheck=1;
    }


  /* Initializations. */
  tarray=tile->array;
  gsl_rng_set(pp->rng, seed);
  pp->up_vals->flag &= ~GAL_DATA_FLAG_SORT_CH;


  /* Set the range of random values for this tile. */
  upperlimit_random_range(pp, tile, min, max, clumplab);


  /* 'se_inc' is just used temporarily, the important thing here is
     'st_oo'. */
  st_oo = ( clumplab
            ? gal_tile_start_end_ind_inclusive(tile, p->objects, se_inc)
            : pp->st_o );
  st_oc = clumplab ? (int32_t *)(p->clumps->array) + se_inc[0] : NULL;


  /* Continue measuring randomly until we get the desired total number. */
  while(nfailed<maxfails && counter<p->upnum)
    {
      /* Get the random coordinates. */
      for(d=0;d<ndim;++d)
        rcoord[d] = upperlimit_random_position(pp, tile, d, min, max);

      /* Set the tile's new starting pointer. */
      tile->array = gal_pointer_increment(p->objects->array,
                          gal_dimension_coord_to_index(ndim, dsize, rcoord),
                                           p->objects->type);

      /* Starting and ending coordinates for this random position, note
         that in 'pp' we have the starting and ending coordinates of the
         actual tile. */
      increment     = 0;
      num_increment = 1;
      continueparse = 1;
      sum           = 0.0f;

      /* Starting pointers for the random tile. */
      st_v   = gal_tile_start_end_ind_inclusive(tile, p->values, se_inc);
      st_o               = (int32_t *)(p->objects->array) + se_inc[0];
      if(p->upmask) st_m = (uint8_t *)(p->upmask->array)  + se_inc[0];

      /* Parse over this object/clump. */
      while( se_inc[0] + increment <= se_inc[1] )
        {
          /* Set the pointers. */
          V               = st_v  + increment;    /* Random tile.   */
          O               = st_o  + increment;    /* Random tile.   */
          if(st_m) M      = st_m  + increment;    /* Random tile.   */
          oO              = st_oo + increment;    /* Original tile. */
          if(clumplab) oC = st_oc + increment;    /* Original tile. */


          /* Parse over this contiguous region, similar to the first and
             second pass functions. */
          OO = O + tile->dsize[ndim-1];
          do
            {
              /* Only use pixels over this object/clump. */
              if( *oO==pp->object && ( oC==NULL || *oC==clumplab ) )
                {
                  /* If this pixel is a non-zero object code, or is masked,
                     or has a blank value, then stop parsing. */
                  if( *O || (M && *M) || ( p->hasblank && isnan(*V) ) )
                    continueparse=0;
                  else
                    sum += *V;
                }

              /* Increment the other pointers. */
              ++V;
              ++oO;
              if(M) ++M;
              if(oC) ++oC;
            }
          while(continueparse && ++O<OO);


          /* Increment to the next contiguous region of this tile. */
          if(continueparse)
            increment += ( gal_tile_block_increment(p->objects, dsize,
                                                    num_increment++, NULL) );
          else break;
        }


      /* Further processing is only necessary if this random tile was fully
         parsed. If it was, we must reset 'nfailed' to zero again. */
      if(continueparse)
        {
          nfailed=0;
          uparr[ counter++ ] = sum;
        }
      else ++nfailed;


      /* If a check is necessary, write in the values (in FITS
         coordinates). */
      if(writecheck)
        {
          switch(ndim)
            {
            case 2:
              gal_list_sizet_add(&check_x, rcoord[1]+1);
              gal_list_sizet_add(&check_y, rcoord[0]+1);
              break;

            case 3:
              gal_list_sizet_add(&check_x, rcoord[2]+1);
              gal_list_sizet_add(&check_y, rcoord[1]+1);
              gal_list_sizet_add(&check_z, rcoord[0]+1);
              break;

            default:
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s "
                    "to fix the problem. 'ndim' value of %zu is not "
                    "recognized", __func__, PACKAGE_BUGREPORT, ndim);
            }
          gal_list_f32_add(&check_s, continueparse ? sum : NAN);
        }
    }

  /* If a check is necessary, then write the values. */
  if(writecheck)
    upperlimit_write_check(p, check_x, check_y, check_z, check_s);

  /* Do the measurement on the random distribution. */
  upperlimit_measure(pp, clumplab, counter==p->upnum);

  /* Reset the tile's array pointer, clean up and return. */
  free(rcoord);
  tile->array=tarray;
  gal_list_f32_free(check_s);
  gal_list_sizet_free(check_x);
  gal_list_sizet_free(check_y);
}




















/*********************************************************************/
/*******************     High level function      ********************/
/*********************************************************************/
void
upperlimit_calculate(struct mkcatalog_passparams *pp)
{
  size_t i;
  unsigned long seed;
  gal_data_t *clumptiles;
  struct mkcatalogparams *p=pp->p;

  /* First find the upper limit magnitude for this object. */
  upperlimit_one_tile(pp, pp->tile, p->rng_seed+pp->object, 0);

  /* If a clumps image is present (a clump catalog is requested) and this
     object has clumps, then find the upper limit magnitude for the clumps
     within this object. */
  if(p->clumps && pp->clumpsinobj)
    {
      /* If an upper-limit check image is requested, then make sure that
         the clump label is not more than the number of clumps in this
         object. */
      if( p->checkuplim[0] == pp->object
          && p->checkuplim[1] != GAL_BLANK_INT32
          && p->checkuplim[1] > pp->clumpsinobj )
        error(EXIT_FAILURE, 0, "object %d has %zu clumps, but an upperlimit "
              "check table (using the '--checkuplim' option) has been "
              "requested for clump %d", pp->object, pp->clumpsinobj,
              p->checkuplim[1]);

      /* Make tiles covering the clumps. */
      clumptiles=upperlimit_make_clump_tiles(pp);

      /* Go over all the clumps. The random number generator seed for each
         clump/object has to be unique, but also reproducible (given the
         intial seed and identical inputs). So we have defined it based on
         the total number of objects and clumps and this object and clump's
         IDs. */
      for(i=0;i<pp->clumpsinobj;++i)
        {
          seed = p->rng_seed + p->numobjects + p->numclumps * pp->object + i;
          upperlimit_one_tile(pp, &clumptiles[i], seed, i+1);
        }

      /* Clean up the clump tiles. */
      gal_data_array_free(clumptiles, pp->clumpsinobj, 0);
    }
}
