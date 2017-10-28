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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include <gnuastro/git.h>
#include <gnuastro/wcs.h>
#include <gnuastro/data.h>
#include <gnuastro/fits.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>

#include <gnuastro-internal/timing.h>

#include "main.h"
#include "mkcatalog.h"

#include "ui.h"
#include "columns.h"
#include "upperlimit.h"





/*********************************************************************/
/*************      Definitions and initialization     ***************/
/*********************************************************************/
/* Both passes are going to need their starting pointers set, so we'll do
   that here. */
static void
mkcatalog_initialize_params(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;

  /* Initialize the number of clumps in this object. */
  pp->clumpsinobj=0;


  /* Initialize the intermediate values. */
  memset(pp->oi, 0, OCOL_NUMCOLS * sizeof *pp->oi);


  /* Set the shifts in every dimension to avoid round-off errors in large
     numbers for the non-linear calculations. We are using the first pixel
     of each object's tile as the shift parameter to keep the mean
     (average) reasonably near to the standard deviation. Otherwise, when
     the object is far out in the image (large x and y positions), then
     roundoff errors are going to decrease the accuracy of the second order
     calculations. */
  gal_dimension_index_to_coord( ( (float *)(pp->tile->array)
                                  - (float *)(pp->tile->block->array) ),
                                p->input->ndim, p->input->dsize, pp->shift);


  /* Set the starting and ending indexs of this tile/object.  */
  pp->st_i   = gal_tile_start_end_ind_inclusive(pp->tile, p->input,
                                                pp->start_end_inc);
  pp->st_sky = (float *)(p->sky->array)       + pp->start_end_inc[0];
  pp->st_std = (float *)(p->std->array)       + pp->start_end_inc[0];
  pp->st_o   = (int32_t *)(p->objects->array) + pp->start_end_inc[0];
  pp->st_c   = ( p->clumps
                 ? (int32_t *)(p->clumps->array)  + pp->start_end_inc[0]
                 : NULL );
}




















/*********************************************************************/
/*************         First and second passes         ***************/
/*********************************************************************/
static void
mkcatalog_first_pass(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  double *oi=pp->oi;
  int32_t *O, *C=NULL;
  size_t d, increment=0, num_increment=1;
  float ss, *I, *II, *SK, *ST, *input=p->input->array;
  size_t *c=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__, "c");
  size_t *sc=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__, "sc");


  /* Parse each contiguous patch of memory covered by this object. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the contiguous range to parse, we will check the count
         over the `I' pointer and just increment the rest. */
      O  = pp->st_o   + increment;
      SK = pp->st_sky + increment;
      ST = pp->st_std + increment;
      if(p->clumps) C = pp->st_c + increment;
      II = ( I = pp->st_i + increment ) + pp->tile->dsize[ndim-1];

      /* Parse the tile. */
      do
        {
          /* If this pixel belongs to the requested object then do the
             processing.  */
          if( *O==pp->object )
            {
              /* Get the number of clumps in this object: the largest clump
                 ID over each object. */
              if( p->clumps && *C>0 )
                pp->clumpsinobj = *C > pp->clumpsinobj ? *C : pp->clumpsinobj;


              /* Get the coordinates of this point. */
              gal_dimension_index_to_coord(I-input, ndim, dsize, c);


              /* Calculate the shifted coordinates for second order
                 calculations. The coordinate is incremented because from
                 now on, the positions are in the FITS standard (starting
                 from one).

                 IMPORTANT NOTE: this is a postfix increment, so after the
                 expression (difference) is evaluated, the coordinate is
                 going to change. This is necessary because `shift' is also
                 starting from zero.  */
              for(d=0;d<ndim;++d) sc[d] = c[d]++ - pp->shift[d];


              /* Do the general geometric (independent of pixel value)
                 calculations. */
              oi[ OCOL_NUMALL ]++;
              oi[ OCOL_GX     ] += c[1];
              oi[ OCOL_GY     ] += c[0];
              oi[ OCOL_GXX    ] += sc[1] * sc[1];
              oi[ OCOL_GYY    ] += sc[0] * sc[0];
              oi[ OCOL_GXY    ] += sc[1] * sc[0];
              if(p->clumps && *C>0)
                {
                  oi[ OCOL_C_GX    ] += c[1];
                  oi[ OCOL_C_GY    ] += c[0];
                }


              /* Start the pixel value related parameters.

                 ABOUT THE CHECK: The reason this condition is given
                 like this is that the `threshold' value is optional
                 and we don't want to do multiple checks.

                 The basic idea is this: when the user doesn't want any
                 thresholds applied, then `p->threshold==NAN' and any
                 conditional that involves a NaN will fail, so its logical
                 negation will be positive and the calculations below will
                 be done. However, if the user does specify a threhold and
                 the pixel is above the threshold, then (`ss < p->threshold
                 * *ST') will be false and its logical negation will be
                 positive, so the pixel will be included. */
              if( !( p->hasblank && isnan(*I) )
                  && !( (ss = *I - *SK) < p->threshold * *ST ) )
                {
                  /* General flux summations. */
                  oi[ OCOL_NUM    ]++;
                  oi[ OCOL_SUM    ] += ss;
                  oi[ OCOL_SUMSKY ] += *SK;
                  oi[ OCOL_SUMSTD ] += *ST;
                  if(p->clumps && *C>0)
                    {
                      oi[ OCOL_C_NUM ]++;
                      oi[ OCOL_C_SUM ] += ss;
                    }

                  /* For flux weighted centers, we can only use positive
                     values, so do those measurements here. */
                  if( ss > 0.0f )
                    {
                      oi[ OCOL_NUMWHT ]++;
                      oi[ OCOL_SUMWHT ] += ss;
                      oi[ OCOL_VX     ] += ss * c[1];
                      oi[ OCOL_VY     ] += ss * c[0];
                      oi[ OCOL_VXX    ] += ss * sc[1] * sc[1];
                      oi[ OCOL_VYY    ] += ss * sc[0] * sc[0];
                      oi[ OCOL_VXY    ] += ss * sc[1] * sc[0];
                      if(p->clumps && *C>0)
                        {
                          oi[ OCOL_C_NUMWHT ]++;
                          oi[ OCOL_C_SUMWHT ] += ss;
                          oi[ OCOL_C_VX     ] += ss * c[1];
                          oi[ OCOL_C_VY     ] += ss * c[0];
                        }
                    }
                }
            }

          /* Increment the other pointers. */
          ++O; ++SK; ++ST; if(p->clumps) ++C;
        }
      while(++I<II);

      /* Increment to the next contiguous region of this tile. */
      increment += ( gal_tile_block_increment(p->input, dsize,
                                              num_increment++, NULL) );
    }

  /* Clean up. */
  free(c);
  free(sc);
}





/* Do the second pass  */
static void
mkcatalog_second_pass(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  double *ci;
  int32_t *O, *C=NULL, nlab, *ngblabs;
  size_t i, ii, d, increment=0, num_increment=1;
  size_t nngb=gal_dimension_num_neighbors(ndim);
  size_t *dinc=gal_dimension_increment(ndim, dsize);
  float ss, *I, *II, *SK, *ST, *input=p->input->array;
  int32_t *objects=p->objects->array, *clumps=p->clumps->array;
  size_t *c=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__, "c");
  size_t *sc=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__, "sc");

  /* Allocate array to keep the neighbor labels. */
  ngblabs=gal_data_malloc_array(GAL_TYPE_INT32, nngb, __func__, "ngblabs");

  /* Parse each contiguous patch of memory covered by this object. */
  while( pp->start_end_inc[0] + increment <= pp->start_end_inc[1] )
    {
      /* Set the contiguous range to parse, we will check the count
         over the `I' pointer and just increment the rest. */
      O  = pp->st_o   + increment;
      SK = pp->st_sky + increment;
      ST = pp->st_std + increment;
      if(p->clumps) C = pp->st_c + increment;
      II = ( I = pp->st_i + increment ) + pp->tile->dsize[ndim-1];

      /* Parse the next contiguous region of this tile. */
      do
        {
          /* If this pixel belongs to the requested object, is a clumps and
             isn't NAN, then do the processing. `hasblank' is constant, so
             when the input doesn't have any blank values, the `isnan' will
             never be checked. */
          if( *O==pp->object )
            {
              /* We are on a clump. */
              if(p->clumps && *C>0)
                {
                  /* Pointer to make things easier. Note that the clump
                     labels start from 1, but the array indexs from 0.*/
                  ci=&pp->ci[ (*C-1) * CCOL_NUMCOLS ];

                  /* Get the coordinates of this point. */
                  gal_dimension_index_to_coord(I-input, ndim, dsize, c);

                  /* Shifted coordinates for second order moments, see
                     explanations in the first pass.*/
                  for(d=0;d<ndim;++d) sc[d] = c[d]++ - pp->shift[d];

                  /* Geometric measurements (independent of pixel value). */
                  ci[ CCOL_NUMALL ]++;
                  ci[ CCOL_GX     ] += c[1];
                  ci[ CCOL_GY     ] += c[0];
                  ci[ CCOL_GXX    ] += sc[1] * sc[1];
                  ci[ CCOL_GYY    ] += sc[0] * sc[0];
                  ci[ CCOL_GXY    ] += sc[1] * sc[0];

                  /* Only use pixels above the threshold, see explanations in
                     first pass for an explanation. */
                  if( !( p->hasblank && isnan(*I) )
                      && !( (ss = *I - *SK) < p->threshold * *ST ) )
                    {
                      /* Fill in the necessary information. */
                      ci[ CCOL_NUM    ]++;
                      ci[ CCOL_SUM    ] += ss;
                      ci[ CCOL_SUMSKY ] += *SK;
                      ci[ CCOL_SUMSTD ] += *ST;
                      if( ss > 0.0f )
                        {
                          ci[ CCOL_NUMWHT ]++;
                          ci[ CCOL_SUMWHT ] += ss;
                          ci[ CCOL_VX     ] += ss * c[1];
                          ci[ CCOL_VY     ] += ss * c[0];
                          ci[ CCOL_VXX    ] += ss * sc[1] * sc[1];
                          ci[ CCOL_VYY    ] += ss * sc[0] * sc[0];
                          ci[ CCOL_VXY    ] += ss * sc[1] * sc[0];
                        }
                    }
                }

              /* This pixel is on the diffuse region, check to see if it is
                 touching a clump or not, but only if this object actually
                 has any clumps. */
              else if(pp->clumpsinobj)
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
                  GAL_DIMENSION_NEIGHBOR_OP(I-input, ndim, dsize, ndim, dinc,
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
                               ngblabs[ii++] = nlab;

                               ++(pp->ci)[ (nlab-1) * CCOL_NUMCOLS
                                           + CCOL_RIV_NUM ];
                               pp->ci[ (nlab-1) * CCOL_NUMCOLS
                                       + CCOL_RIV_SUM ] += *I-*SK;
                             }
                         }
                     });
                }
            }

          /* Increment the other pointers. */
          ++O; ++SK; ++ST; if(p->clumps) ++C;
        }
      while(++I<II);

      /* Increment to the next contiguous region of this tile. */
      increment += ( gal_tile_block_increment(p->input, dsize,
                                              num_increment++, NULL) );
    }

  /* Clean up. */
  free(c);
  free(sc);
  free(dinc);
  free(ngblabs);
}




















/*********************************************************************/
/*****************       High-level funcitons      *******************/
/*********************************************************************/
static void
mkcatalog_clump_starting_index(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;

  /* Lock the mutex if we are working on more than one thread. NOTE: it is
     very important to keep the number of operations within the mutex to a
     minimum so other threads don't get delayed. */
  if(p->cp.numthreads>1)
    pthread_mutex_lock(&p->mutex);

  /* Put the current total number of rows filled into the output, then
     increment the total number by the number of clumps. */
  pp->clumpstartindex = p->clumprowsfilled;
  p->clumprowsfilled += pp->clumpsinobj;

  /* Unlock the mutex (if it was locked). */
  if(p->cp.numthreads>1)
    pthread_mutex_unlock(&p->mutex);
}





/* Each thread will call this function once. It will go over all the
   objects that are assigned to it. */
static void *
mkcatalog_single_object(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct mkcatalogparams *p=(struct mkcatalogparams *)(tprm->params);
  size_t ndim=p->input->ndim;

  size_t i;
  struct mkcatalog_passparams pp;

  /* Initialize the mkcatalog_passparams elements. */
  pp.p               = p;
  pp.clumpstartindex = 0;
  pp.rng             = p->rng ? gsl_rng_clone(p->rng) : NULL;
  pp.shift           = gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                             "pp.shift");
  pp.oi              = gal_data_malloc_array(GAL_TYPE_FLOAT64, OCOL_NUMCOLS,
                                             __func__, "pp.oi");

  /* If we have upper-limit mode, then allocate the container to keep the
     values to calculate the standard deviation. */
  pp.up_vals = p->upperlimit ? gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1,
                                              &p->upnum, NULL, 0,
                                              p->cp.minmapsize, NULL, NULL,
                                              NULL) : NULL;

  /* Fill the desired columns for all the objects given to this thread. */
  for(i=0; tprm->indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading, Note that the object IDs start from one while
         the array positions start from 0. */
      pp.ci=NULL;
      pp.object = tprm->indexs[i] + 1;
      pp.tile   = &p->tiles[ tprm->indexs[i] ];

      /* Initialize the parameters for this object/tile. */
      mkcatalog_initialize_params(&pp);

      /* Get the first pass information. */
      mkcatalog_first_pass(&pp);

      /* Currently the second pass is only necessary when there is a clumps
         image. */
      if(p->clumps)
        {
          /* Allocate space for the properties of each clump. */
          pp.ci = gal_data_calloc_array(GAL_TYPE_FLOAT64,
                                        pp.clumpsinobj * CCOL_NUMCOLS,
                                        __func__, "pp.ci");

          /* Get the starting row of this object's clumps in the final
             catalog. This index is also necessary for the unique random
             number generator seeds of each clump. */
          mkcatalog_clump_starting_index(&pp);

          /* Get the second pass information. */
          mkcatalog_second_pass(&pp);
        }

      /* Calculate the upper limit magnitude (if necessary). */
      if(p->upperlimit) upperlimit_calculate(&pp);

      /* Write the pass information into the columns. */
      columns_fill(&pp);

      /* Clean up for this object. */
      if(pp.ci) free(pp.ci);
    }

  /* Clean up. */
  free(pp.oi);
  free(pp.shift);
  gal_data_free(pp.up_vals);
  if(pp.rng) gsl_rng_free(pp.rng);

  /* Wait until all the threads finish and return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}




















/*********************************************************************/
/********         Processing after threads finish        *************/
/*********************************************************************/
/* Convert internal image coordinates to WCS for table.

   Note that from the beginning (during the passing steps), we saved FITS
   coordinates. Also note that we are doing the conversion in place. */
static void
mkcatalog_wcs_conversion(struct mkcatalogparams *p)
{
  gal_data_t *c;
  gal_data_t *column;

  /* Flux weighted center positions for clumps and objects. */
  if(p->wcs_vo)
    {
      gal_wcs_img_to_world(p->wcs_vo, p->input->wcs, 1);
      if(p->wcs_vc)
        gal_wcs_img_to_world(p->wcs_vc, p->input->wcs, 1);
    }


  /* Geometric center positions for clumps and objects. */
  if(p->wcs_go)
    {
      gal_wcs_img_to_world(p->wcs_go, p->input->wcs, 1);
      if(p->wcs_gc)
        gal_wcs_img_to_world(p->wcs_gc, p->input->wcs, 1);
    }


  /* All clumps flux weighted center. */
  if(p->wcs_vcc)
    gal_wcs_img_to_world(p->wcs_vcc, p->input->wcs, 1);


  /* All clumps geometric center. */
  if(p->wcs_gcc)
    gal_wcs_img_to_world(p->wcs_gcc, p->input->wcs, 1);


  /* Go over all the object columns and fill in the values. */
  for(column=p->objectcols; column!=NULL; column=column->next)
    {
      /* Definitions */
      c=NULL;

      /* Set `c' for the columns that must be corrected. Note that this
         `switch' statement doesn't need any `default', because there are
         probably columns that don't need any correction. */
      switch(column->status)
        {
        case UI_KEY_W1:           c=p->wcs_vo;                break;
        case UI_KEY_W2:           c=p->wcs_vo->next;          break;
        case UI_KEY_GEOW1:        c=p->wcs_go;                break;
        case UI_KEY_GEOW2:        c=p->wcs_go->next;          break;
        case UI_KEY_CLUMPSW1:     c=p->wcs_vcc;               break;
        case UI_KEY_CLUMPSW2:     c=p->wcs_vcc->next;         break;
        case UI_KEY_CLUMPSGEOW1:  c=p->wcs_gcc;               break;
        case UI_KEY_CLUMPSGEOW2:  c=p->wcs_gcc->next;         break;
        }

      /* Copy the elements into the output column. */
      if(c)
        memcpy(column->array, c->array,
               column->size*gal_type_sizeof(c->type));
    }


  /* Go over all the clump columns and fill in the values. */
  for(column=p->clumpcols; column!=NULL; column=column->next)
    {
      /* Definitions */
      c=NULL;

      /* Set `c' for the columns that must be corrected. Note that this
         `switch' statement doesn't need any `default', because there are
         probably columns that don't need any correction. */
      switch(column->status)
        {
        case UI_KEY_W1:           c=p->wcs_vc;                break;
        case UI_KEY_W2:           c=p->wcs_vc->next;          break;
        case UI_KEY_GEOW1:        c=p->wcs_gc;                break;
        case UI_KEY_GEOW2:        c=p->wcs_gc->next;          break;
        }

      /* Copy the elements into the output column. */
      if(c)
        memcpy(column->array, c->array,
               column->size*gal_type_sizeof(c->type));
    }
}





/* Write the similar information. */
static gal_list_str_t *
mkcatalog_outputs_same_start(struct mkcatalogparams *p, int o0c1,
                             char *ObjClump)
{
  float snlim;
  char *str, *tstr;
  double pixarea=NAN;
  gal_list_str_t *comments=NULL;
  char *skyfile=p->skyfile ? p->skyfile : p->inputname;
  char *stdfile=p->stdfile ? p->stdfile : p->inputname;
  char *clumpsfile=p->clumpsfile ? p->clumpsfile : p->inputname;
  char *objectsfile=p->objectsfile ? p->objectsfile : p->inputname;

  asprintf(&str, "%s catalog of %s", o0c1 ? "Object" : "Clump",
           PROGRAM_STRING);
  gal_list_str_add(&comments, str, 0);

  /* If in a Git controlled directory and output isn't a FITS file (in
     FITS, this will be automatically included). */
  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT && gal_git_describe())
    {
      asprintf(&str, "Working directory commit %s", gal_git_describe());
      gal_list_str_add(&comments, str, 0);
    }

  /* Write the date. However, `ctime' is going to put a new-line character
     in the end of its string, so we are going to remove it manually. */
  asprintf(&str, "%s started on %s", PROGRAM_NAME, ctime(&p->rawtime));
  str[strlen(str)-1]='\0';
  gal_list_str_add(&comments, str, 0);


  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      asprintf(&str, "--------- Input files ---------");
      gal_list_str_add(&comments, str, 0);
    }

  asprintf(&str, "Values:  %s (hdu: %s).", p->inputname, p->cp.hdu);
  gal_list_str_add(&comments, str, 0);

  asprintf(&str, "Objects: %s (hdu: %s).", objectsfile, p->objectshdu);
  gal_list_str_add(&comments, str, 0);

  if(p->clumps)
    {
      asprintf(&str, "Clumps:  %s (hdu: %s).", clumpsfile, p->clumpshdu);
      gal_list_str_add(&comments, str, 0);
    }

  asprintf(&str, "Sky:     %s (hdu: %s).", skyfile, p->skyhdu);
  gal_list_str_add(&comments, str, 0);

  asprintf(&str, "Sky STD: %s (hdu: %s).", stdfile, p->stdhdu);
  gal_list_str_add(&comments, str, 0);

  if(p->upmaskfile)
    {
      asprintf(&str, "Upperlimit mask: %s (hdu: %s).", p->upmaskfile,
               p->upmaskhdu);
      gal_list_str_add(&comments, str, 0);
    }

  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      asprintf(&str, "--------- Supplimentary information ---------");
      gal_list_str_add(&comments, str, 0);
    }

  if(p->input->wcs)
    {
      pixarea=gal_wcs_pixel_area_arcsec2(p->input->wcs);
      if( isnan(pixarea)==0 )
        {
          asprintf(&str, "Pixel area (arcsec^2): %g", pixarea);
          gal_list_str_add(&comments, str, 0);
        }
    }

  if(p->hasmag)
    {
      asprintf(&str, "Zeropoint magnitude: %.4f", p->zeropoint);
      gal_list_str_add(&comments, str, 0);
    }

  /* Print surface brightness limits. */
  if( !isnan(p->zeropoint) &&  !isnan(p->sfmagnsigma) )
    {
      /* Per pixel. */
      asprintf(&str, "%g sigma surface brightness (magnitude/pixel): %.3f",
               p->sfmagnsigma, ( -2.5f
                                 *log10( p->sfmagnsigma
                                         * p->medstd )
                                 + p->zeropoint ) );
      gal_list_str_add(&comments, str, 0);

      /* Per requested projected area: if a pixel area could be measured (a
         WCS was given), then also estimate the surface brightness over one
         arcsecond^2. From the pixel area, we know how many pixels are
         necessary to fill the requested projected area (in
         arcsecond^2). We also know that as the number of samples (pixels)
         increases (to N), the noise increases by sqrt(N), see the full
         discussion in the book. */
      if(!isnan(pixarea) && !isnan(p->sfmagarea))
        {
          /* Prepare the comment/information. */
          if(p->sfmagarea==1.0f) tstr=NULL;
          else                   asprintf(&tstr, "%g-", p->sfmagarea);
          asprintf(&str, "%g sigma surface brightness (magnitude/%sarcsec^2): "
                   "%.3f", p->sfmagnsigma, tstr ? tstr : "",
                   ( -2.5f * log10( p->sfmagnsigma
                                    * p->medstd
                                    * sqrt( p->sfmagarea / pixarea) )
                     + p->zeropoint ) );

          /* Add the final string/line to the catalog comments. */
          gal_list_str_add(&comments, str, 0);

          /* Clean up (if necessary). */
          if (tstr)
            {
              free(tstr);
              tstr=NULL;
            }
        }

      /* Notice: */
      asprintf(&str, "Pixel STD for surface brightness calculation%s: %f",
               (!isnan(pixarea) && !isnan(p->sfmagarea))?"s":"", p->medstd);
      gal_list_str_add(&comments, str, 0);
    }

  snlim = o0c1 ? p->clumpsn : p->detsn;
  if( !isnan(snlim) )
    {
      asprintf(&str, "%s limiting signal-to-noise ratio: %.3f", ObjClump,
               snlim);
      gal_list_str_add(&comments, str, 0);
    }

  if(o0c1==0)
    {
      asprintf(&str, "(NOTE: S/N limit above is for pseudo-detections, "
               "not objects.)");
      gal_list_str_add(&comments, str, 0);
    }

  if(p->cpscorr>1.0f)
    {
      asprintf(&str, "Counts-per-second correction: %.3f", p->cpscorr);
      gal_list_str_add(&comments, str, 0);
    }

  if( !isnan(p->threshold) )
    {
      asprintf(&str, "**IMPORTANT** Pixel threshold (multiple of local "
               "std): %.3f", p->threshold);
      gal_list_str_add(&comments, str, 0);
    }


  if(p->upperlimit)
    {
      if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
        {
          asprintf(&str, "--------- Upper-limit measurement ---------");
          gal_list_str_add(&comments, str, 0);
        }

      asprintf(&str, "Number of random samples: %zu", p->upnum);
      gal_list_str_add(&comments, str, 0);

      asprintf(&str, "Random number generator name: %s", p->rngname);
      gal_list_str_add(&comments, str, 0);

      asprintf(&str, "Random number generator seed: %"PRIu64, p->seed);
      gal_list_str_add(&comments, str, 0);

      asprintf(&str, "Multiple of STD used for sigma-clipping: %.3f",
               p->upsigmaclip[0]);
      gal_list_str_add(&comments, str, 0);

      if(p->upsigmaclip[1]>=1.0f)
        asprintf(&str, "Number of clips for sigma-clipping: %.0f",
                 p->upsigmaclip[1]);
      else
        asprintf(&str, "Tolerance level to sigma-clipping: %.3f",
                 p->upsigmaclip[1]);
      gal_list_str_add(&comments, str, 0);

      asprintf(&str, "Multiple of sigma-clipped STD for upper-limit: %.3f",
               p->upnsigma);
      gal_list_str_add(&comments, str, 0);
    }



  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      asprintf(&str, "--------- Table columns ---------");
      gal_list_str_add(&comments, str, 0);
    }

  /* Return the comments. */
  return comments;
}





/* Write the produced columns into the output */
static void
mkcatalog_write_outputs(struct mkcatalogparams *p)
{
  /*char *str;*/
  gal_list_str_t *comments;


  /* OBJECT CATALOG
     ============== */
  comments=mkcatalog_outputs_same_start(p, 0, "Detection");


  /* Write objects catalog
     ---------------------

     Reverse the comments list (so it is printed in the same order here),
     write the objects catalog and free the comments. */
  gal_list_str_reverse(&comments);
  gal_table_write(p->objectcols, comments, p->cp.tableformat, p->objectsout,
                  "OBJECTS");
  gal_list_str_free(comments, 1);



  /* CLUMPS CATALOG
     ============== */
  if(p->clumps)
    {
      comments=mkcatalog_outputs_same_start(p, 1, "Clumps");



      /* Write objects catalog
         ---------------------

         Reverse the comments list (so it is printed in the same order here),
         write the objects catalog and free the comments. */
      gal_list_str_reverse(&comments);
      gal_table_write(p->clumpcols, comments, p->cp.tableformat, p->clumpsout,
                      "CLUMPS");
      gal_list_str_free(comments, 1);
    }
}




















/*********************************************************************/
/*****************       Top-level function        *******************/
/*********************************************************************/
void
mkcatalog(struct mkcatalogparams *p)
{
  /* When more than one thread is to be used, initialize the mutex: we need
     it to assign a column to the clumps in the final catalog. */
  if( p->cp.numthreads > 1 ) pthread_mutex_init(&p->mutex, NULL);


  /* Do the processing on each thread. */
  gal_threads_spin_off(mkcatalog_single_object, p, p->numobjects,
                       p->cp.numthreads);


  /* Post-thread processing, for example to convert image coordinates to RA
     and Dec. */
  mkcatalog_wcs_conversion(p);


  /* Write the filled columns into the output. */
  mkcatalog_write_outputs(p);


  /* Destroy the mutex. */
  if( p->cp.numthreads>1 ) pthread_mutex_destroy(&p->mutex);
}
