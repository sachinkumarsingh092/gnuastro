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
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>
#include <gnuastro/permutation.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

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
  pp.oi              = gal_pointer_allocate(GAL_TYPE_FLOAT64, OCOL_NUMCOLS,
                                            0, __func__, "pp.oi");

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
               ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                       "pp.shift")
               : NULL );

  /* If we have upper-limit mode, then allocate the container to keep the
     values to calculate the standard deviation. */
  if(p->upperlimit)
    {
      /* Allocate the space to keep the upper-limit values. */
      pp.up_vals = gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &p->upnum,
                                  NULL, 0, p->cp.minmapsize, p->cp.quietmmap,
                                  NULL, NULL, NULL);

      /* Set the blank checked flag to 1. By definition, this dataset won't
         have any blank values. Also 'flag' is initialized to '0'. So we
         just have to set the checked flag ('GAL_DATA_FLAG_BLANK_CH') to
         one to inform later steps that there are no blank values. */
      pp.up_vals->flag |= GAL_DATA_FLAG_BLANK_CH;
    }
  else
    pp.up_vals=NULL;


  /* Fill the desired columns for all the objects given to this thread. */
  for(i=0; tprm->indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. Note that the object IDs start from one while
         the array positions start from 0. */
      pp.ci       = NULL;
      pp.object   = ( p->outlabs
                      ? p->outlabs[tprm->indexs[i]]
                      : tprm->indexs[i] + 1 );
      pp.tile     = &p->tiles[   tprm->indexs[i] ];
      pp.spectrum = &p->spectra[ tprm->indexs[i] ];

      /* Initialize the parameters for this object/tile. */
      parse_initialize(&pp);

      /* Get the first pass information. */
      parse_objects(&pp);

      /* Currently the second pass is only necessary when there is a clumps
         image. */
      if(p->clumps)
        {
          /* Allocate space for the properties of each clump. */
          pp.ci = gal_pointer_allocate(GAL_TYPE_FLOAT64,
                                       pp.clumpsinobj * CCOL_NUMCOLS, 1,
                                       __func__, "pp.ci");

          /* Get the starting row of this object's clumps in the final
             catalog. This index is also necessary for the unique random
             number generator seeds of each clump. */
          mkcatalog_clump_starting_index(&pp);

          /* Get the second pass information. */
          parse_clumps(&pp);
        }

      /* If an order-based calculation is requested, another pass is
         necessary. */
      if( p->oiflag[ OCOL_MEDIAN ]
          || p->oiflag[ OCOL_SIGCLIPNUM ]
          || p->oiflag[ OCOL_SIGCLIPSTD ]
          || p->oiflag[ OCOL_SIGCLIPMEAN ]
          || p->oiflag[ OCOL_SIGCLIPMEDIAN ])
        parse_order_based(&pp);

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

      /* Set 'c' for the columns that must be corrected. Note that this
         'switch' statement doesn't need any 'default', because there are
         probably columns that don't need any correction. */
      switch(column->status)
        {
        case UI_KEY_W1:           c=p->wcs_vo;                break;
        case UI_KEY_W2:           c=p->wcs_vo->next;          break;
        case UI_KEY_W3:           c=p->wcs_vo->next->next;    break;
        case UI_KEY_GEOW1:        c=p->wcs_go;                break;
        case UI_KEY_GEOW2:        c=p->wcs_go->next;          break;
        case UI_KEY_GEOW3:        c=p->wcs_go->next->next;    break;
        case UI_KEY_CLUMPSW1:     c=p->wcs_vcc;               break;
        case UI_KEY_CLUMPSW2:     c=p->wcs_vcc->next;         break;
        case UI_KEY_CLUMPSW3:     c=p->wcs_vcc->next->next;   break;
        case UI_KEY_CLUMPSGEOW1:  c=p->wcs_gcc;               break;
        case UI_KEY_CLUMPSGEOW2:  c=p->wcs_gcc->next;         break;
        case UI_KEY_CLUMPSGEOW3:  c=p->wcs_gcc->next->next;   break;
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

      /* Set 'c' for the columns that must be corrected. Note that this
         'switch' statement doesn't need any 'default', because there are
         probably columns that don't need any correction. */
      switch(column->status)
        {
        case UI_KEY_W1:           c=p->wcs_vc;                break;
        case UI_KEY_W2:           c=p->wcs_vc->next;          break;
        case UI_KEY_W3:           c=p->wcs_vc->next->next;    break;
        case UI_KEY_GEOW1:        c=p->wcs_gc;                break;
        case UI_KEY_GEOW2:        c=p->wcs_gc->next;          break;
        case UI_KEY_GEOW3:        c=p->wcs_gc->next->next;    break;
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
  char *tmp, *str;

  /* Basic classifiers for plain text outputs. */
  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      if( asprintf(&str, "--------- Input files ---------")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  /* Object labels. */
  if( asprintf(&str, "Objects: %s (hdu: %s).", p->objectsfile, p->cp.hdu)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_str_add(comments, str, 0);

  /* Clump labels. */
  if(p->clumps)
    {
      if(asprintf(&str, "Clumps:  %s (hdu: %s).", p->usedclumpsfile,
                  p->clumpshdu)<0)
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  /* Values dataset. */
  if(p->values)
    {
      if( asprintf(&str, "Values:  %s (hdu: %s).", p->usedvaluesfile,
                   p->valueshdu)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(comments, str, 0);
    }

  /* Sky dataset. */
  if(withsky && p->sky)
    {
      if(p->sky->size==1)
        {
          if( asprintf(&str, "Sky:     %g.", *((float *)(p->sky->array)) )<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        {
          if( asprintf(&str, "Sky:     %s (hdu: %s).", p->usedskyfile,
                       p->skyhdu)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      gal_list_str_add(comments, str, 0);
    }

  /* Sky standard deviation dataset. */
  tmp = p->variance ? "VAR" : "STD";
  if(withstd && p->std)
    {
      if(p->std->size==1)
        {
          if( asprintf(&str, "Sky %s: %g.", tmp,
                       *((float *)(p->std->array)) )<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        {
          if( asprintf(&str, "Sky %s: %s (hdu: %s).", tmp, p->usedstdfile,
                       p->stdhdu)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      gal_list_str_add(comments, str, 0);
    }

  /* Upper limit mask. */
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

  /* Write the date. However, 'ctime' is going to put a new-line character
     in the end of its string, so we are going to remove it manually. */
  if( asprintf(&str, "%s started on %s", PROGRAM_NAME, ctime(&p->rawtime))<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  str[strlen(str)-1]='\0';
  gal_list_str_add(&comments, str, 0);


  /* Write the basic information. */
  mkcatalog_write_inputs_in_comments(p, &comments, 1, 1);


  /* Write other supplementary information. */
  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      if( asprintf(&str, "--------- Supplementary information ---------")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  /* Pixel area. */
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

  /* Zeropoint magnitude */
  if(p->hasmag)
    {
      if( asprintf(&str, "Zeropoint magnitude: %.4f", p->zeropoint)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  /* Print surface brightness limits. */
  if( !isnan(p->medstd) && !isnan(p->sfmagnsigma) )
    {
      /* Only print magnitudes if a zeropoint is given. */
      if( !isnan(p->zeropoint) )
        {
          /* Per pixel. */
          if( asprintf(&str, "%g sigma surface brightness (magnitude/pixel): "
                       "%.3f", p->sfmagnsigma, ( -2.5f
                                                 *log10( p->sfmagnsigma
                                                         * p->medstd )
                                                 + p->zeropoint ) )<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_list_str_add(&comments, str, 0);

          /* Requested projected area: if a pixel area could be measured (a
             WCS was given), then also estimate the surface brightness over
             one arcsecond^2. From the pixel area, we know how many pixels
             are necessary to fill the requested projected area (in
             arcsecond^2). We also know that as the number of samples
             (pixels) increases (to N), the noise increases by sqrt(N), see
             the full discussion in the book. */
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
        }

      /* Notice: */
      if( asprintf(&str, "Pixel STD for surface brightness calculation%s: %f",
                   (!isnan(pixarea) && !isnan(p->sfmagarea))?"s":"",
                   p->medstd)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }
  else
    {
      gal_checkset_allocate_copy("No surface brightness calcuations "
                                 "because no STD image used.", &str);
      gal_list_str_add(&comments, str, 0);
      gal_checkset_allocate_copy("Ask for column that uses the STD image, "
                                 "or '--forcereadstd'.", &str);
      gal_list_str_add(&comments, str, 0);
    }

  /* The count-per-second correction. */
  if(p->cpscorr>1.0f)
    {
      if( asprintf(&str, "Counts-per-second correction: %.3f", p->cpscorr)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  /* Print upper-limit parameters. */
  if(p->upperlimit)
    upperlimit_write_comments(p, &comments, 1);

  /* Start column metadata. */
  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      if( asprintf(&str, "--------- Table columns ---------")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, str, 0);
    }

  /* Return the comments. */
  return comments;
}






/* Since all the measurements were done in parallel (and we didn't know the
   number of clumps per object a-priori), the clumps informtion is just
   written in as they are measured. Here, we'll sort the clump columns by
   object ID. There is an option to disable this. */
static void
sort_clumps_by_objid(struct mkcatalogparams *p)
{
  gal_data_t *col;
  size_t o, i, j, *permute, *rowstart;

  /* Make sure everything is fine. */
  if(p->hostobjid_c==NULL || p->numclumps_c==NULL)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
          "problem. 'p->hostobjid_c' and 'p->numclumps_c' must not be "
          "NULL.", __func__, PACKAGE_BUGREPORT);


  /* Allocate the necessary arrays. */
  rowstart=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->numobjects, 0, __func__,
                                 "rowstart");
  permute=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->numclumps, 0, __func__,
                                "permute");


  /* The objects array is already sorted by object ID. So we should just
     add up the number of clumps to find the row where each object's clumps
     should start from in the final sorted clumps catalog. */
  rowstart[0] = 0;
  for(i=1;i<p->numobjects;++i)
    rowstart[i] = p->numclumps_c[i-1] + rowstart[i-1];

  /* Fill the permutation array. Note that WE KNOW that all the objects for
     one clump are after each other.*/
  i=0;
  while(i<p->numclumps)
    {
      o=p->hostobjid_c[i]-1;
      for(j=0; j<p->numclumps_c[o]; ++j)
        permute[i++] = rowstart[o] + j;
    }

  /* Permute all the clump columns. */
  for(col=p->clumpcols; col!=NULL; col=col->next)
    gal_permutation_apply_inverse(col, permute);

  /* Clean up */
  free(permute);
  free(rowstart);
}





/* Write the produced columns into the output */
static void
mkcatalog_write_outputs(struct mkcatalogparams *p)
{
  size_t i, scounter;
  char str[200], *fname;
  gal_list_str_t *comments;
  int outisfits=gal_fits_name_is_fits(p->objectsout);

  /* If a catalog is to be generated. */
  if(p->objectcols)
    {
      /* OBJECT catalog */
      comments=mkcatalog_outputs_same_start(p, 0, "Detection");

      /* Reverse the comments list (so it is printed in the same order
         here), write the objects catalog and free the comments. */
      gal_list_str_reverse(&comments);
      gal_table_write(p->objectcols, comments, p->cp.tableformat,
                      p->objectsout, "OBJECTS", 0);
      gal_list_str_free(comments, 1);


      /* CLUMPS catalog */
      if(p->clumps)
        {
          /* Make the comments. */
          comments=mkcatalog_outputs_same_start(p, 1, "Clumps");

          /* Write objects catalog
             ---------------------

             Reverse the comments list (so it is printed in the same order
             here), write the objects catalog and free the comments. */
          gal_list_str_reverse(&comments);
          gal_table_write(p->clumpcols, comments, p->cp.tableformat,
                          p->clumpsout, "CLUMPS", 0);
          gal_list_str_free(comments, 1);
        }
    }

  /* Spectra. */
  if(p->spectra)
    {
      /* Inform the user (Writing many FITS extensions can take
         long). */
      if(p->objectcols && outisfits)
        printf("  - Catalog(s) complete, writing spectra.\n");

      /* Start counting and writing the files. Note that due to some
         conditions (for example in debugging), a 'p->spectra[i]' may not
         actually contain any data. So we'll also count the number of
         spectra that are written. */
      scounter=0;
      for(i=0;i<p->numobjects;++i)
        if(p->spectra[i].ndim>0)
          {
            /* Increment the written spectra-counter. */
            ++scounter;

            /* Write the spectra based on the requested format. */
            if(outisfits)
              {
                /* Write the table. */
                sprintf(str, "SPECTRUM_%zu", i+1);
                gal_table_write(&p->spectra[i], NULL, GAL_TABLE_FORMAT_BFITS,
                                p->objectsout, str, 0);
              }
            else
              {
                sprintf(str, "-spec-%zu.txt", i+1);
                fname=gal_checkset_automatic_output(&p->cp, p->objectsout,
                                                    str);
                gal_table_write(&p->spectra[i], NULL, GAL_TABLE_FORMAT_TXT,
                                fname, NULL, 0);
                free(fname);
              }
          }
    }

  /* Configuration information. */
  if(outisfits)
    {
      gal_fits_key_write_filename("input", p->objectsfile, &p->cp.okeys, 1);
      gal_fits_key_write_config(&p->cp.okeys, "MakeCatalog configuration",
                                "MKCATALOG-CONFIG", p->objectsout, "0");
    }


  /* Inform the user */
  if(!p->cp.quiet)
    {
      if(p->objectcols)
        {
          if(p->clumpsout && strcmp(p->clumpsout,p->objectsout))
            {
              printf("  - Output objects catalog: %s\n", p->objectsout);
              if(p->clumps)
                printf("  - Output clumps catalog: %s\n", p->clumpsout);
            }
          else
            printf("  - Catalog written to %s\n", p->objectsout);

        }

      if(p->spectra)
        {
          if(outisfits)
            {
              if(p->objectcols)
                printf("  - Spectra in %zu extensions named 'SPECTRUM_NN'.\n",
                       p->numobjects);
              else
                printf("  - Output: %s (Spectra in %zu extensions named "
                       "'SPECTRUM_NN').\n)", p->objectsout, p->numobjects);
            }
          else
            printf("  - Spectra in %zu files with '-spec-NN.txt' suffix.\n",
                   p->numobjects);
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

  /* If the columns need to be sorted (by object ID), then some adjustments
     need to be made (possibly to both the objects and clumps catalogs). */
  if(p->hostobjid_c)
    sort_clumps_by_objid(p);

  /* Write the filled columns into the output. */
  mkcatalog_write_outputs(p);

  /* Destroy the mutex. */
  if( p->cp.numthreads>1 ) pthread_mutex_destroy(&p->mutex);
}
