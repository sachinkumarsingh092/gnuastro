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
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>

#include "main.h"
#include "mkcatalog.h"

#include "ui.h"
#include "parse.h"
#include "columns.h"
#include "upperlimit.h"










/*********************************************************************/
/**************       Manage a single object       *******************/
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
  size_t ndim=p->objects->ndim;

  size_t i;
  uint8_t *oif=p->oiflag;
  struct mkcatalog_passparams pp;


  /* Initialize the mkcatalog_passparams elements. */
  pp.p               = p;
  pp.clumpstartindex = 0;
  pp.rng             = p->rng ? gsl_rng_clone(p->rng) : NULL;
  pp.oi              = gal_data_malloc_array(GAL_TYPE_FLOAT64, OCOL_NUMCOLS,
                                             __func__, "pp.oi");

  /* If we have second order measurements, allocate the array keeping the
     temporary shift values for each object of this thread. Note that the
     clumps catalog (if requested), will have the same measurements, so its
     just enough to check the objects. */
  pp.shift = ( ( oif[    OCOL_GXX ]
                 || oif[ OCOL_GYY ]
                 || oif[ OCOL_GXY ]
                 || oif[ OCOL_VXX ]
                 || oif[ OCOL_VYY ]
                 || oif[ OCOL_VXY ] )
               ? gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                       "pp.shift")
               : NULL );

  /* If we have upper-limit mode, then allocate the container to keep the
     values to calculate the standard deviation. */
  pp.up_vals = ( p->upperlimit
                 ? gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &p->upnum,
                                  NULL, 0, p->cp.minmapsize, NULL, NULL,
                                  NULL)
                 : NULL );

  /* Fill the desired columns for all the objects given to this thread. */
  for(i=0; tprm->indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. Note that the object IDs start from one while
         the array positions start from 0. */
      pp.ci=NULL;
      pp.object = tprm->indexs[i] + 1;
      pp.tile   = &p->tiles[ tprm->indexs[i] ];

      /* Initialize the parameters for this object/tile. */
      parse_initialize(&pp);

      /* Get the first pass information. */
      parse_objects(&pp);

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
          parse_clumps(&pp);
        }

      /* If the median is requested, another pass is necessary. */
      if( p->oiflag[ OCOL_MEDIAN ] )
        parse_median(&pp);

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
      gal_wcs_img_to_world(p->wcs_vo, p->objects->wcs, 1);
      if(p->wcs_vc)
        gal_wcs_img_to_world(p->wcs_vc, p->objects->wcs, 1);
    }


  /* Geometric center positions for clumps and objects. */
  if(p->wcs_go)
    {
      gal_wcs_img_to_world(p->wcs_go, p->objects->wcs, 1);
      if(p->wcs_gc)
        gal_wcs_img_to_world(p->wcs_gc, p->objects->wcs, 1);
    }


  /* All clumps flux weighted center. */
  if(p->wcs_vcc)
    gal_wcs_img_to_world(p->wcs_vcc, p->objects->wcs, 1);


  /* All clumps geometric center. */
  if(p->wcs_gcc)
    gal_wcs_img_to_world(p->wcs_gcc, p->objects->wcs, 1);


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





void
mkcatalog_write_inputs_in_comments(struct mkcatalogparams *p,
                                   gal_list_str_t **comments, int withsky,
                                   int withstd)
{
  char *str;

  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      if( asprintf(&str, "--------- Input files ---------")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  if( asprintf(&str, "Objects: %s (hdu: %s).", p->objectsfile, p->cp.hdu)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_str_add(comments, str, 0);

  if(p->clumps)
    {
      if(asprintf(&str, "Clumps:  %s (hdu: %s).", p->usedclumpsfile,
                  p->clumpshdu)<0)
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  if(p->values)
    {
      if( asprintf(&str, "Values:  %s (hdu: %s).", p->usedvaluesfile,
                   p->valueshdu)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  if(withsky && p->sky)
    {
      if( asprintf(&str, "Sky:     %s (hdu: %s).", p->usedskyfile,
                   p->skyhdu)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  if(withstd && p->std)
    {
      if( asprintf(&str, "Sky STD: %s (hdu: %s).", p->usedstdfile,
                   p->stdhdu)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  if(p->upmaskfile)
    {
      if( asprintf(&str, "Upperlimit mask: %s (hdu: %s).", p->upmaskfile,
                   p->upmaskhdu)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }
}





/* Write the similar information. */
static gal_list_str_t *
mkcatalog_outputs_same_start(struct mkcatalogparams *p, int o0c1,
                             char *ObjClump)
{
  char *str, *tstr;
  double pixarea=NAN;
  gal_list_str_t *comments=NULL;

  if( asprintf(&str, "%s catalog of %s", o0c1 ? "Object" : "Clump",
               PROGRAM_STRING)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_str_add(&comments, str, 0);

  /* If in a Git controlled directory and output isn't a FITS file (in
     FITS, this will be automatically included). */
  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT && gal_git_describe())
    {
      if(asprintf(&str, "Working directory commit %s", gal_git_describe())<0)
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  /* Write the date. However, `ctime' is going to put a new-line character
     in the end of its string, so we are going to remove it manually. */
  if( asprintf(&str, "%s started on %s", PROGRAM_NAME, ctime(&p->rawtime))<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  str[strlen(str)-1]='\0';
  gal_list_str_add(&comments, str, 0);


  /* Write the basic information. */
  mkcatalog_write_inputs_in_comments(p, &comments, 1, 1);


  /* Write other supplimentary information. */
  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      if( asprintf(&str, "--------- Supplimentary information ---------")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  if(p->objects->wcs)
    {
      pixarea=gal_wcs_pixel_area_arcsec2(p->objects->wcs);
      if( isnan(pixarea)==0 )
        {
          if( asprintf(&str, "Pixel area (arcsec^2): %g", pixarea)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_list_str_add(&comments, str, 0);
        }
    }

  if(p->hasmag)
    {
      if( asprintf(&str, "Zeropoint magnitude: %.4f", p->zeropoint)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  /* Print surface brightness limits. */
  if( !isnan(p->medstd) && !isnan(p->zeropoint) &&  !isnan(p->sfmagnsigma) )
    {
      /* Per pixel. */
      if( asprintf(&str, "%g sigma surface brightness (magnitude/pixel): "
                   "%.3f", p->sfmagnsigma, ( -2.5f
                                             *log10( p->sfmagnsigma
                                                     * p->medstd )
                                             + p->zeropoint ) )<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);

      /* Requested projected area: if a pixel area could be measured (a WCS
         was given), then also estimate the surface brightness over one
         arcsecond^2. From the pixel area, we know how many pixels are
         necessary to fill the requested projected area (in
         arcsecond^2). We also know that as the number of samples (pixels)
         increases (to N), the noise increases by sqrt(N), see the full
         discussion in the book. */
      if(!isnan(pixarea) && !isnan(p->sfmagarea))
        {
          /* Prepare the comment/information. */
          if(p->sfmagarea==1.0f)
            tstr=NULL;
          else
            if( asprintf(&tstr, "%g-", p->sfmagarea)<0 )
              error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          if( asprintf(&str, "%g sigma surface brightness "
                       "(magnitude/%sarcsec^2): %.3f", p->sfmagnsigma,
                       tstr ? tstr : "",
                       ( -2.5f * log10( p->sfmagnsigma
                                        * p->medstd
                                        * sqrt( p->sfmagarea / pixarea) )
                         + p->zeropoint ) )<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);

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
      if( asprintf(&str, "Pixel STD for surface brightness calculation%s: %f",
                   (!isnan(pixarea) && !isnan(p->sfmagarea))?"s":"",
                   p->medstd)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  if(p->cpscorr>1.0f)
    {
      if( asprintf(&str, "Counts-per-second correction: %.3f", p->cpscorr)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  if(p->upperlimit)
    upperlimit_write_comments(p, &comments, 1);



  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      if( asprintf(&str, "--------- Table columns ---------")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
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

  /* Inform the user. */
  if(!p->cp.quiet)
    {
      if(p->clumpsout==p->objectsout)
        printf("  - Output catalog: %s\n", p->objectsout);
      else
        {
          printf("  - Output objects catalog: %s\n", p->objectsout);
          if(p->clumps)
            printf("  - Output clumps catalog: %s\n", p->clumpsout);
        }
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
