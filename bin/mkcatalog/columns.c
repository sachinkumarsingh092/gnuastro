/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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
#include <stdlib.h>
#include <pthread.h>

#include "main.h"
#include "mkcatalog.h"

#include "ui.h"
#include "columns.h"




/******************************************************************/
/*******************     Intermediate arrays     ******************/
/******************************************************************/
/* Allocate RA-DEC internal arrays. These arrays are defined to keep all
   the positions in one place and do the RA-DEC conversion once in the
   end. They are all allocated together, but we don't know if RA is
   requested first or Dec or if they are requested multiple times. So
   before the allocation, we'll check the first one.

   The space that is allocated in `columns_define_alloc' is for the final
   values that are written in the output file. */
static void
columns_alloc_radec(struct mkcatalogparams *p)
{
  if(p->rd_vo==NULL)
    {
      /* Allocate the space for all dimensions. */
      errno=0;
      p->rd_vo = malloc(p->input->ndim * sizeof *p->rd_vo);
      if(p->rd_vo==NULL)
        error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for p->rd_vo",
              __func__, p->input->ndim * sizeof *p->rd_vo );

      /* Space for each dimension. */
      p->rd_vo[0] = gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numobjects,
                                          __func__, "p->rd_vo[0]");
      p->rd_vo[1] = gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numobjects,
                                          __func__, "p->rd_vo[1]");
      if(p->clumps)
        {
          /* Allocate the space for all dimensions. */
          errno=0;
          p->rd_vc = malloc(p->input->ndim * sizeof *p->rd_vc);
          if(p->rd_vc==NULL)
            error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for p->rd_vo",
                  __func__, p->input->ndim * sizeof *p->rd_vc );

          /* Space for each dimension. */
          p->rd_vc[0]=gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numclumps,
                                            __func__, "p->rd_vc[0]");
          p->rd_vc[1]=gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numclumps,
                                            __func__, "p->rd_vc[1]");
        }
    }
}





/* Similar to `columns_alloc_radec'. */
static void
columns_alloc_georadec(struct mkcatalogparams *p)
{
  if(p->rd_go==NULL)
    {
      /* Allocate the space for all dimensions. */
      errno=0;
      p->rd_go = malloc(p->input->ndim * sizeof *p->rd_go);
      if(p->rd_go==NULL)
        error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for `p->rd_go'",
              __func__, p->input->ndim * sizeof *p->rd_go );

      /* Space for each dimension. */
      p->rd_go[0] = gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numobjects,
                                          __func__, "p->rd_go[0]");
      p->rd_go[1] = gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numobjects,
                                          __func__, "p->rd_go[1]");
      if(p->clumps)
        {
          /* Allocate the space for all dimensions. */
          errno=0;
          p->rd_gc = malloc(p->input->ndim * sizeof *p->rd_gc);
          if(p->rd_gc==NULL)
            error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for `p->rd_go'",
                  __func__, p->input->ndim * sizeof *p->rd_gc );

          /* Space for each dimension. */
          p->rd_gc[0]=gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numclumps,
                                            __func__, "p->rd_gc[0]");
          p->rd_gc[1]=gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numclumps,
                                            __func__, "p->rd_gc[1]");
        }
    }
}





/* Similar to `columns_alloc_radec'. */
static void
columns_alloc_clumpsradec(struct mkcatalogparams *p)
{
  if(p->rd_vcc==NULL)
    {
      /* Allocate the space for all dimensions. */
      errno=0;
      p->rd_vcc = malloc(p->input->ndim * sizeof *p->rd_vcc);
      if(p->rd_vcc==NULL)
        error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for `p->rd_vcc'",
              __func__, p->input->ndim * sizeof *p->rd_vcc );

      /* Space for each dimension. */
      p->rd_vcc[0] = gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numobjects,
                                           __func__, "p->rd_vcc[0]");
      p->rd_vcc[1] = gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numobjects,
                                           __func__, "p->rd_vcc[1]");
    }
}





/* Similar to `columns_alloc_radec'. */
static void
columns_alloc_clumpsgeoradec(struct mkcatalogparams *p)
{
  if(p->rd_gcc==NULL)
    {
      /* Allocate the space for all dimensions. */
      errno=0;
      p->rd_gcc = malloc(p->input->ndim * sizeof *p->rd_gcc);
      if(p->rd_gcc==NULL)
        error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for `p->rd_gcc'",
              __func__, p->input->ndim * sizeof *p->rd_gcc );

      /* Space for each dimension. */
      p->rd_gcc[0] = gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numobjects,
                                           __func__, "p->rd_gcc[0]");
      p->rd_gcc[1] = gal_data_malloc_array(GAL_TYPE_FLOAT64, p->numobjects,
                                           __func__, "p->rd_gcc[1]");
    }
}




















/******************************************************************/
/**********       Column definition/allocation      ***************/
/******************************************************************/
/* Set the necessary parameters for each output column and allocate the
   space necessary to keep the values. */
void
columns_define_alloc(struct mkcatalogparams *p)
{
  gal_list_i32_t *colcode;
  gal_list_str_t *strtmp, *noclumpimg=NULL;
  int disp_fmt=0, disp_width=0, disp_precision=0;
  char *name=NULL, *unit=NULL, *ocomment=NULL, *ccomment=NULL;
  uint8_t otype=GAL_TYPE_INVALID, ctype=GAL_TYPE_INVALID, *oiflag, *ciflag;

  /* Allocate the array for which intermediate parameters are
     necessary. The basic issue is that higher-level calculations require a
     smaller domain of raw measurements. So to avoid having to calculate
     something multiple times, each parameter will flag the intermediate
     parameters it requires in these arrays. */
  oiflag = p->oiflag = gal_data_malloc_array(GAL_TYPE_UINT8, OCOL_NUMCOLS,
                                             __func__, "oiflag");
  ciflag = p->ciflag = gal_data_malloc_array(GAL_TYPE_UINT8, CCOL_NUMCOLS,
                                             __func__, "ciflag");

  /* Allocate the columns. */
  for(colcode=p->columnids; colcode!=NULL; colcode=colcode->next)
    {
      /* Set the column-specific parameters, please follow the same order
         as `args.h'. IMPORTANT: we want the names to be the same as the
         option names. Note that zero `disp_' variables will be
         automatically determined.*/
      switch(colcode->v)
        {
        case UI_KEY_OBJID:
          name           = "OBJ_ID";
          unit           = "counter";
          ocomment       = "Object identifier.";
          ccomment       = NULL;
          otype          = GAL_TYPE_INT32;  /* Same type as clumps image. */
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          /* Is an internal parameter. */
          break;

        case UI_KEY_HOSTOBJID:
          name           = "HOST_OBJ_ID";
          unit           = "counter";
          ocomment       = NULL;
          ccomment       = "Object identifier hosting this clump.";
          otype          = GAL_TYPE_INVALID;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          /* Is an internal parameter. */
          break;

        case UI_KEY_IDINHOSTOBJ:
          name           = "ID_IN_HOST_OBJ";
          unit           = "counter";
          ocomment       = NULL;
          ccomment       = "ID of clump in its host object.";
          otype          = GAL_TYPE_INVALID;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          /* Is an internal parameter. */
          break;

        case UI_KEY_NUMCLUMPS:
          name           = "NUM_CLUMPS";
          unit           = "counter";
          ocomment       = "Number of clumps in this object.";
          ccomment       = NULL;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = 0;
          disp_width     = 5;
          disp_precision = 0;
          /* Is an internal parameter. */
          break;

        case UI_KEY_AREA:
          name           = "AREA";
          unit           = "counter";
          ocomment       = "Number of pixels covered.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 5;
          disp_precision = 0;
          oiflag[ OCOL_NUM ] = 1;
          ciflag[ CCOL_NUM ] = 1;
          break;

        case UI_KEY_CLUMPSAREA:
          name           = "CLUMPS_AREA";
          unit           = "counter";
          ocomment       = "Total number of clump pixels in object.";
          ccomment       = NULL;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = 0;
          disp_width     = 5;
          disp_precision = 0;
          oiflag[ OCOL_C_NUM ] = 1;
          break;

        case UI_KEY_X:
          name           = "X";
          unit           = "position";
          ocomment       = "Flux weighted center (FITS axis 1).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_VX ] = 1;
          ciflag[ CCOL_VX ] = 1;
          break;

        case UI_KEY_Y:
          name           = "Y";
          unit           = "position";
          ocomment       = "Flux weighted center (FITS axis 2).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_VY ] = 1;
          ciflag[ CCOL_VY ] = 1;
          break;

        case UI_KEY_GEOX:
          name           = "GEO_X";
          unit           = "position";
          ocomment       = "Geometric center (FITS axis 1).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_GX ] = 1;
          ciflag[ CCOL_GX ] = 1;
          break;

        case UI_KEY_GEOY:
          name           = "GEO_Y";
          unit           = "position";
          ocomment       = "Geometric center (FITS axis 2).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_GY ] = 1;
          ciflag[ CCOL_GY ] = 1;
          break;

        case UI_KEY_CLUMPSX:
          name           = "CLUMPS_X";
          unit           = "position";
          ocomment       = "Flux weighted center of clumps (FITS axis 1).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_VX ] = 1;
          break;

        case UI_KEY_CLUMPSY:
          name           = "CLUMPS_Y";
          unit           = "position";
          ocomment       = "Flux weighted center of clumps (FITS axis 2).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_VY ] = 1;
          break;

        case UI_KEY_CLUMPSGEOX:
          name           = "CLUMPS_GEO_X";
          unit           = "position";
          ocomment       = "Geometric center of clumps (FITS axis 1).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_GX ] = 1;
          break;

        case UI_KEY_CLUMPSGEOY:
          name           = "CLUMPS_GEO_Y";
          unit           = "position";
          ocomment       = "Geometric center of clumps (FITS axis 2).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_GY ] = 1;
          break;

        case UI_KEY_RA:
          name           = "RA";
          unit           = "degrees";
          ocomment       = "Flux weighted center right ascension.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          oiflag[ OCOL_C_VY ] = 1;
          oiflag[ OCOL_C_VY ] = 1;
          columns_alloc_radec(p);
          break;

        case UI_KEY_DEC:
          name           = "DEC";
          unit           = "degrees";
          ocomment       = "Flux weighted center declination.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          oiflag[ OCOL_C_VY ] = 1;
          oiflag[ OCOL_C_VY ] = 1;
          columns_alloc_radec(p);
          break;

        case UI_KEY_GEORA:
          name           = "GEO_RA";
          unit           = "degrees";
          ocomment       = "Geometric center right ascension.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          oiflag[ OCOL_GX ] = 1;
          oiflag[ OCOL_GY ] = 1;
          ciflag[ CCOL_GX ] = 1;
          ciflag[ CCOL_GY ] = 1;
          columns_alloc_georadec(p);
          break;

        case UI_KEY_GEODEC:
          name           = "GEO_DEC";
          unit           = "degrees";
          ocomment       = "Geometric center declination.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          oiflag[ OCOL_GX ] = 1;
          oiflag[ OCOL_GY ] = 1;
          ciflag[ CCOL_GX ] = 1;
          ciflag[ CCOL_GY ] = 1;
          columns_alloc_georadec(p);
          break;

        case UI_KEY_CLUMPSRA:
          name           = "CLUMPS_RA";
          unit           = "degrees";
          ocomment       = "RA of all clumps flux weighted center.";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          oiflag[ OCOL_C_VX ] = 1;
          oiflag[ OCOL_C_VY ] = 1;
          columns_alloc_clumpsradec(p);
          break;

        case UI_KEY_CLUMPSDEC:
          name           = "CLUMPS_DEC";
          unit           = "degrees";
          ocomment       = "Declination of all clumps flux weighted center.";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 15;
          disp_precision = 7;
          oiflag[ OCOL_C_VX ] = 1;
          oiflag[ OCOL_C_VY ] = 1;
          columns_alloc_clumpsradec(p);
          break;

        case UI_KEY_CLUMPSGEORA:
          name           = "CLUMPS_RA";
          unit           = "degrees";
          ocomment       = "RA of all clumps geometric center.";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          oiflag[ OCOL_C_GX ] = 1;
          oiflag[ OCOL_C_GY ] = 1;
          columns_alloc_clumpsgeoradec(p);
          break;

        case UI_KEY_CLUMPSGEODEC:
          name           = "CLUMPS_DEC";
          unit           = "degrees";
          ocomment       = "Declination of all clumps geometric center.";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          oiflag[ OCOL_C_GX ] = 1;
          oiflag[ OCOL_C_GY ] = 1;
          columns_alloc_clumpsgeoradec(p);
          break;

        case UI_KEY_BRIGHTNESS:
          name           = "BRIGHTNESS";
          unit           = p->input->unit ? p->input->unit : "pixelunit";
          ocomment       = "Brightness (sum of sky subtracted values).";
          ccomment       = "Brightness (sum of pixels subtracted by rivers).";
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 4;
          oiflag[ OCOL_SUM ]     = 1;
          ciflag[ CCOL_SUM ]     = 1;
          ciflag[ CCOL_RIV_NUM ] = 1;
          ciflag[ CCOL_RIV_SUM ] = 1;
          break;

        case UI_KEY_CLUMPSBRIGHTNESS:
          name           = "CLUMPS_BRIGHTNESS";
          unit           = p->input->unit ? p->input->unit : "pixelunit";
          ocomment       = "Brightness (sum of pixel values) in clumps.";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 4;
          oiflag[ OCOL_C_SUM ] = 1;
          break;

        case UI_KEY_NORIVERBRIGHTNESS:
          name           = "NO_RIVER_BRIGHTNESS";
          unit           = p->input->unit ? p->input->unit : "pixelunit";
          ocomment       = NULL;
          ccomment       = "Brightness (sum of sky subtracted values).";
          otype          = GAL_TYPE_INVALID;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 4;
          ciflag[ CCOL_SUM ] = 1;
          break;

        case UI_KEY_MAGNITUDE:
          name           = "MAGNITUDE";
          unit           = "log";
          ocomment       = "Magnitude.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 8;
          disp_precision = 3;
          oiflag[ OCOL_SUM ] = 1;
          ciflag[ CCOL_SUM ] = 1;
          p->hasmag      = 1;
          break;

        case UI_KEY_MAGNITUDEERR:
          name           = "MAGNITUDE_ERR";
          unit           = "log";
          ocomment       = "Error in measuring magnitude.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 8;
          disp_precision = 3;
          oiflag[ OCOL_SUMSTD ] = 1;
          oiflag[ OCOL_NUM    ] = 1;
          oiflag[ OCOL_SUM    ] = 1;
          ciflag[ CCOL_SUMSTD ] = 1;
          ciflag[ CCOL_NUM    ] = 1;
          ciflag[ CCOL_SUM    ] = 1;
          break;

        case UI_KEY_CLUMPSMAGNITUDE:
          name           = "CLUMPS_MAGNITUDE";
          unit           = "log";
          ocomment       = "Magnitude in all clumps.";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 8;
          disp_precision = 3;
          oiflag[ OCOL_C_SUM ] = 1;
          p->hasmag      = 1;
          break;

        case UI_KEY_UPPERLIMIT:
          name           = "UPPERLIMIT";
          unit           = p->input->unit;
          ocomment       = "Upper limit value (random positionings).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 8;
          disp_precision = 3;
          p->upperlimit  = 1;
          /* Upper limit measurement doesn't need per-pixel calculations. */
          break;

        case UI_KEY_UPPERLIMITMAG:
          name           = "UPPERLIMIT_MAG";
          unit           = "log";
          ocomment       = "Upper limit magnitude (random positionings).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 8;
          disp_precision = 3;
          p->upperlimit  = 1;
          p->hasmag      = 1;
          /* Upper limit magnitude doesn't need per-pixel calculations. */
          break;

        case UI_KEY_RIVERAVE:
          name           = "RIVER_AVE";
          unit           = p->input->unit ? p->input->unit : "pixelunit";
          ocomment       = NULL;
          ccomment       = "Average river value surrounding this clump.";
          otype          = GAL_TYPE_INVALID;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 4;
          ciflag[ CCOL_RIV_NUM ] = 1;
          ciflag[ CCOL_RIV_SUM ] = 1;
          break;

        case UI_KEY_RIVERNUM:
          name           = "RIVER_NUM";
          unit           = "counter";
          ocomment       = NULL;
          ccomment       = "Number of river pixels around this clump.";
          otype          = GAL_TYPE_INVALID;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 5;
          disp_precision = 0;
          ciflag[ CCOL_RIV_NUM ] = 1;
          break;

        case UI_KEY_SN:
          name           = "SN";
          unit           = "ratio";
          ocomment       = "Signal to noise ratio.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_SUMSTD ] = 1;
          oiflag[ OCOL_NUM    ] = 1;
          oiflag[ OCOL_SUM    ] = 1;
          ciflag[ CCOL_SUMSTD ] = 1;
          ciflag[ CCOL_NUM    ] = 1;
          ciflag[ CCOL_SUM    ] = 1;
          break;

        case UI_KEY_SKY:
          name           = "SKY";
          unit           = p->input->unit ? p->input->unit : "pixelunit";
          ocomment       = "Average input sky value.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 4;
          oiflag[ OCOL_NUM    ] = 1;
          oiflag[ OCOL_SUMSKY ] = 1;
          ciflag[ CCOL_NUM    ] = 1;
          ciflag[ CCOL_SUMSKY ] = 1;
          break;

        case UI_KEY_STD:
          name           = "STD";
          unit           = p->input->unit ? p->input->unit : "pixelunit";
          ocomment       = "Average of input standard deviation.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 4;
          oiflag[ OCOL_NUM    ] = 1;
          oiflag[ OCOL_SUMSTD ] = 1;
          ciflag[ CCOL_NUM    ] = 1;
          ciflag[ CCOL_SUMSTD ] = 1;
          break;

        case UI_KEY_SEMIMAJOR:
          name           = "SEMI_MAJOR";
          unit           = "pixel";
          ocomment       = "Flux weighted semi-major axis.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_VXX ] = 1;
          oiflag[ OCOL_VYY ] = 1;
          oiflag[ OCOL_VXY ] = 1;
          ciflag[ CCOL_VXX ] = 1;
          ciflag[ CCOL_VYY ] = 1;
          ciflag[ CCOL_VXY ] = 1;
          break;

        case UI_KEY_SEMIMINOR:
          name           = "SEMI_MINOR";
          unit           = "pixel";
          ocomment       = "Flux weighted semi-minor axis.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_VXX ] = 1;
          oiflag[ OCOL_VYY ] = 1;
          oiflag[ OCOL_VXY ] = 1;
          ciflag[ CCOL_VXX ] = 1;
          ciflag[ CCOL_VYY ] = 1;
          ciflag[ CCOL_VXY ] = 1;
          break;

        case UI_KEY_AXISRATIO:
          name           = "AXIS_RATIO";
          unit           = "ratio";
          ocomment       = "Flux weighted axis ratio (minor/major).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 7;
          disp_precision = 3;
          oiflag[ OCOL_VXX ] = 1;
          oiflag[ OCOL_VYY ] = 1;
          oiflag[ OCOL_VXY ] = 1;
          ciflag[ CCOL_VXX ] = 1;
          ciflag[ CCOL_VYY ] = 1;
          ciflag[ CCOL_VXY ] = 1;
          break;

        case UI_KEY_POSITIONANGLE:
          name           = "POSITION_ANGLE";
          unit           = "degrees";
          ocomment       = "Position angle.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_VXX ] = 1;
          oiflag[ OCOL_VYY ] = 1;
          oiflag[ OCOL_VXY ] = 1;
          ciflag[ CCOL_VXX ] = 1;
          ciflag[ CCOL_VYY ] = 1;
          ciflag[ CCOL_VXY ] = 1;
          break;

        case UI_KEY_GEOSEMIMAJOR:
          name           = "GEO_SEMI_MAJOR";
          unit           = "pixel";
          ocomment       = "Geometric semi-major axis.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_GXX ] = 1;
          oiflag[ OCOL_GYY ] = 1;
          oiflag[ OCOL_GXY ] = 1;
          ciflag[ CCOL_GXX ] = 1;
          ciflag[ CCOL_GYY ] = 1;
          ciflag[ CCOL_GXY ] = 1;
          break;

        case UI_KEY_GEOSEMIMINOR:
          name           = "GEO_SEMI_MINOR";
          unit           = "pixel";
          ocomment       = "Geometric semi-minor axis.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_GXX ] = 1;
          oiflag[ OCOL_GYY ] = 1;
          oiflag[ OCOL_GXY ] = 1;
          ciflag[ CCOL_GXX ] = 1;
          ciflag[ CCOL_GYY ] = 1;
          ciflag[ CCOL_GXY ] = 1;
          break;

        case UI_KEY_GEOAXISRATIO:
          name           = "GEO_AXIS_RATIO";
          unit           = "ratio";
          ocomment       = "Geometric axis ratio (minor/major).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 7;
          disp_precision = 3;
          oiflag[ OCOL_VXX ] = 1;
          oiflag[ OCOL_VYY ] = 1;
          oiflag[ OCOL_VXY ] = 1;
          ciflag[ CCOL_VXX ] = 1;
          ciflag[ CCOL_VYY ] = 1;
          ciflag[ CCOL_VXY ] = 1;
          break;

        case UI_KEY_GEOPOSITIONANGLE:
          name           = "GEO_POSITION_ANGLE";
          unit           = "degrees";
          ocomment       = "Geometric Position angle.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_GXX ] = 1;
          oiflag[ OCOL_GYY ] = 1;
          oiflag[ OCOL_GXY ] = 1;
          ciflag[ CCOL_GXX ] = 1;
          ciflag[ CCOL_GYY ] = 1;
          ciflag[ CCOL_GXY ] = 1;
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to fix "
                "the problem. The code %d is not an internally recognized "
                "column code", __func__, PACKAGE_BUGREPORT, colcode->v);
        }

      /* If this is an objects column, add it to the list of columns. We
         will be using the `status' element to keep the MakeCatalog code
         for the columns. */
      if(otype!=GAL_TYPE_INVALID)
        {
          gal_list_data_add_alloc(&p->objectcols, NULL, otype, 1,
                                  &p->numobjects, NULL, 0, p->cp.minmapsize,
                                  name, unit, ocomment);
          p->objectcols->status         = colcode->v;
          p->objectcols->disp_fmt       = disp_fmt;
          p->objectcols->disp_width     = disp_width;
          p->objectcols->disp_precision = disp_precision;
        }


      /* Similar to the objects column above but for clumps, but since the
         clumps image is optional, we need a further check before actually
         allocating the column. */
      if(ctype!=GAL_TYPE_INVALID)
        {
          /* A clumps image has been given, so allocate space for this
             column. */
          if(p->clumps)
            {
              gal_list_data_add_alloc(&p->clumpcols, NULL, ctype, 1,
                                      &p->numclumps, NULL, 0,
                                      p->cp.minmapsize, name, unit, ccomment);
              p->clumpcols->status         = colcode->v;
              p->clumpcols->disp_fmt       = disp_fmt;
              p->clumpcols->disp_width     = disp_width;
              p->clumpcols->disp_precision = disp_precision;
            }


          /* If this is a clumps-only column and no clumps image was
             given. Add the column to the list of similar columns to inform
             the user. */
          else if(otype==GAL_TYPE_INVALID)
            gal_list_str_add(&noclumpimg, name, 1);
        }
    }


  /* If a warning for clumps columns and no clumps image is necessary make
     the warning. */
  if(noclumpimg)
    {
      gal_list_str_reverse(&noclumpimg);
      fprintf(stderr, "\n-------\n"
              "WARNING: the following column(s) are unique to "
              "clumps (not objects), but the objects image doesn't have "
              " `WCLUMPS' keyword. So these requested columns will be "
              "ignored.\n\n");
      for(strtmp=noclumpimg; strtmp!=NULL; strtmp=strtmp->next)
        fprintf(stderr, "\t%s\n", strtmp->v);
      gal_list_str_free(noclumpimg, 1);
      fprintf(stderr, "\n-------\n");
    }


  /* Free the general columns information because it is no longe needed,
     we'll set it back to NULL afterwards so it is not mistakenly used. */
  gal_list_i32_free(p->columnids);
  p->columnids=NULL;
}




















/******************************************************************/
/**********            Column calculation           ***************/
/******************************************************************/
#define MKC_RATIO(TOP,BOT) ( (BOT)!=0.0f ? (TOP)/(BOT) : NAN )
#define MKC_MAG(B)         ( ((B)>0) ? -2.5f * log10(B) + p->zeropoint : NAN )





/* Calculate the Signal to noise ratio for the object. */
static double
columns_sn(struct mkcatalogparams *p, double *row, int o0c1)
{
  double var, sn, std, Ni, I, O;

  /* Get all the values as averages (per pixel). */
  Ni  = row[ o0c1 ? CCOL_NUM : OCOL_NUM ];
  I   = MKC_RATIO( row[ o0c1 ? CCOL_SUM    : OCOL_SUM ],    Ni );
  std = MKC_RATIO( row[ o0c1 ? CCOL_SUMSTD : OCOL_SUMSTD ], Ni );
  var = (p->skysubtracted ? 2.0f : 1.0f) * std * std;

  /* Calculate the S/N. Note that when grown clumps are requested from
     NoiseChisel, some "clumps" will completely cover their objects and
     there will be no rivers. So if this is a clump, and the river area is
     0, we should treat the S/N as a an object. */
  if( o0c1 && row[ CCOL_RIV_NUM ] )
    {
      /* If the Sky is already subtracted, the varience should be counted
         two times. */
      O   = row[ CCOL_RIV_SUM ] / row[ CCOL_RIV_NUM ];  /* Outside.  */
      sn  = ( sqrt(Ni/p->cpscorr) * (I-O)
              / sqrt( (I>0?I:-1*I) + (O>0?O:-1*O) + var ) );
    }
  else
    sn  = sqrt(Ni/p->cpscorr) * I / sqrt( (I>0?I:-1*I) + var );

  /* Return the derived value. */
  return sn;
}





/* Do the second order calculations, see "Measuring elliptical parameters"
   section of the book/manual for a thorough explanation of the
   derivation. */
static double
columns_second_order(struct mkcatalog_passparams *pp, double *row,
                     int key, int o0c1)
{
  double x=NAN, y=NAN, xx=NAN, yy=NAN, xy=NAN;
  double denom, kx=pp->shift[1]+1, ky=pp->shift[0]+1;

  /* Preparations. */
  switch(key)
    {
    /* Brightness weighted. */
    case UI_KEY_SEMIMAJOR:
    case UI_KEY_SEMIMINOR:
    case UI_KEY_POSITIONANGLE:

      /* Denominator (to be divided). */
      denom = row[ o0c1 ? CCOL_SUMPOS : OCOL_SUMPOS ];

      /* First order. */
      x  = MKC_RATIO( row[ o0c1 ? CCOL_VX     : OCOL_VX     ], denom );
      y  = MKC_RATIO( row[ o0c1 ? CCOL_VY     : OCOL_VY     ], denom );

      /* Second order. */
      xx = ( MKC_RATIO( row[ o0c1 ? CCOL_VXX    : OCOL_VXX    ], denom )
             - (x-kx) * (x-kx) );
      yy = ( MKC_RATIO( row[ o0c1 ? CCOL_VYY    : OCOL_VYY    ], denom )
             - (y-ky) * (y-ky) );
      xy = ( MKC_RATIO( row[ o0c1 ? CCOL_VXY    : OCOL_VXY    ], denom )
             - (x-kx) * (y-ky) );
      break;

    /* Geometric. */
    case UI_KEY_GEOSEMIMAJOR:
    case UI_KEY_GEOSEMIMINOR:
    case UI_KEY_GEOPOSITIONANGLE:

      /* Denominator (to be divided). */
      denom = row[ o0c1 ? CCOL_NUM : OCOL_NUM ];

      /* First order. */
      x  = MKC_RATIO( row[ o0c1 ? CCOL_GX  : OCOL_GX  ], denom );
      y  = MKC_RATIO( row[ o0c1 ? CCOL_GY  : OCOL_GY  ], denom );

      /* Second order. */
      xx = ( MKC_RATIO( row[ o0c1 ? CCOL_GXX : OCOL_GXX ], denom )
             - (x-kx) * (x-kx) );
      yy = ( MKC_RATIO( row[ o0c1 ? CCOL_GYY : OCOL_GYY ], denom )
             - (y-ky) * (y-ky) );
      xy = ( MKC_RATIO( row[ o0c1 ? CCOL_GXY : OCOL_GXY ], denom )
             - (x-kx) * (y-ky) );
      break;

    /* Error. */
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Code %d not a recognized key",
            __func__, key);
    }

  /* Return the output. */
  switch(key)
    {
    /* Semi-major axis. */
    case UI_KEY_SEMIMAJOR:
    case UI_KEY_GEOSEMIMAJOR:
      return sqrt( ( xx + yy ) / 2
                   + sqrt( (xx - yy)/2 * (xx - yy)/2 + xy * xy ) );

    /* Semi-minor axis. */
    case UI_KEY_SEMIMINOR:
    case UI_KEY_GEOSEMIMINOR:
      /*printf("\nhere\n");*/
      return sqrt( ( xx + yy )/2
                   - sqrt( (xx - yy)/2 * (xx - yy)/2 + xy * xy ) );

    /* Position angle. */
    case UI_KEY_POSITIONANGLE:
    case UI_KEY_GEOPOSITIONANGLE:
      return 0.5f * atan2(2 * xy, xx - yy) * 180/M_PI;
    }


  /* Control should not reach here! If it does, its a bug, so abort and let
     the user know. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s, so we can "
        "address the problem. Control should not have reached the end of "
        "this function", __func__, PACKAGE_BUGREPORT);
  return NAN;
}





/* The magnitude error is directly derivable from the S/N:

   To derive the error in measuring the magnitude from the S/N, let's take
   `F' as the flux, `Z' is the zeropoint, `M' is the magnitude, `S' is the
   S/N, and `D' to stand for capital delta (or error in a value) then from

      `M = -2.5*log10(F) + Z'

   we get the following equation after calculating the derivative with
   respect to F.

      `dM/df = -2.5 * ( 1 / ( F * ln(10) ) )'

   From the Tailor series, `DM' can be written as:

      `DM = dM/dF * DF'

   So

      `DM = |-2.5/ln(10)| * DF/F'

   But `DF/F' is just the inverse of the Signal to noise ratio, or
  `1/S'. So

      `DM = 2.5 / ( S * ln(10) )'               */
#define MAG_ERROR(P,ROW,O0C1) (2.5f/(columns_sn((P),(ROW),(O0C1)) * log(10)))






/* All the raw first and second pass information has been collected, now
   write them into the output columns. The list of columns here is in the
   same order as `columns_alloc_set_out_cols', see there for the type of
   each column. */
void
columns_fill(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;

  int key;
  double tmp;
  void *colarr;
  gal_data_t *column;
  double *ci, *oi=pp->oi;
  size_t sr=pp->clumpstartindex, cind, coind;
  size_t oind=pp->object-1; /* IDs start from 1, indexs from 0. */

  /* Go over all the object columns and fill in the information. */
  for(column=p->objectcols; column!=NULL; column=column->next)
    {
      /* For easy reading. */
      key=column->status;
      colarr=column->array;

      /* Go over all the columns. */
      switch(key)
        {
        case UI_KEY_OBJID:
          ((int32_t *)colarr)[oind] = pp->object;
          break;

        case UI_KEY_NUMCLUMPS:
          ((int32_t *)colarr)[oind] = pp->clumpsinobj;
          break;

        case UI_KEY_AREA:
          ((int32_t *)colarr)[oind] = oi[OCOL_NUM];
          break;

        case UI_KEY_CLUMPSAREA:
          ((int32_t *)colarr)[oind] = oi[OCOL_C_NUM];
          break;

        case UI_KEY_X:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_VX], oi[OCOL_SUMPOS] );
          break;

        case UI_KEY_Y:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_VY], oi[OCOL_SUMPOS] );
          break;

        case UI_KEY_GEOX:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_GX], oi[OCOL_NUM] );
          break;

        case UI_KEY_GEOY:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_GY], oi[OCOL_NUM] );
          break;

        case UI_KEY_CLUMPSX:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_C_VX],
                                               oi[OCOL_C_SUMPOS] );
          break;

        case UI_KEY_CLUMPSY:
          ((float *)colarr)[oind] = MKC_RATIO(oi[OCOL_C_VY],
                                              oi[OCOL_C_SUMPOS] );
          break;

        case UI_KEY_CLUMPSGEOX:
          ((float *)colarr)[oind] = MKC_RATIO(oi[OCOL_C_GX], oi[OCOL_C_NUM]);
          break;

        case UI_KEY_CLUMPSGEOY:
          ((float *)colarr)[oind] = MKC_RATIO(oi[OCOL_C_GY], oi[OCOL_C_NUM]);
          break;

        case UI_KEY_RA:
        case UI_KEY_DEC:
          p->rd_vo[0][oind] = MKC_RATIO( oi[OCOL_VX], oi[OCOL_SUMPOS] );
          p->rd_vo[1][oind] = MKC_RATIO( oi[OCOL_VY], oi[OCOL_SUMPOS] );
          break;

        case UI_KEY_GEORA:
        case UI_KEY_GEODEC:
          p->rd_go[0][oind] = MKC_RATIO( oi[OCOL_GX], oi[OCOL_NUM] );
          p->rd_go[1][oind] = MKC_RATIO( oi[OCOL_GY], oi[OCOL_NUM] );
          break;

        case UI_KEY_CLUMPSRA:
        case UI_KEY_CLUMPSDEC:
          p->rd_vcc[0][oind] = MKC_RATIO( oi[OCOL_C_VX], oi[OCOL_C_SUMPOS] );
          p->rd_vcc[1][oind] = MKC_RATIO( oi[OCOL_C_VY], oi[OCOL_C_SUMPOS] );
          break;

        case UI_KEY_CLUMPSGEORA:
        case UI_KEY_CLUMPSGEODEC:
          p->rd_gcc[0][oind] = MKC_RATIO( oi[OCOL_C_GX], oi[OCOL_C_NUM] );
          p->rd_gcc[1][oind] = MKC_RATIO( oi[OCOL_C_GY], oi[OCOL_C_NUM] );
          break;

        case UI_KEY_BRIGHTNESS:
          ((float *)colarr)[oind] = oi[ OCOL_SUM ];
          break;

        case UI_KEY_CLUMPSBRIGHTNESS:
          ((float *)colarr)[oind] = oi[ OCOL_C_SUM ];
          break;

        case UI_KEY_MAGNITUDE:
          ((float *)colarr)[oind] = MKC_MAG(oi[ OCOL_SUM ]);
          break;

        case UI_KEY_MAGNITUDEERR:
          ((float *)colarr)[oind] = MAG_ERROR(p, oi, 0);
          break;

        case UI_KEY_CLUMPSMAGNITUDE:
          ((float *)colarr)[oind] = MKC_MAG(oi[ OCOL_C_SUM ]);
          break;

        case UI_KEY_UPPERLIMIT:
          ((float *)colarr)[oind] = oi[ OCOL_UPPERLIMIT_B ];
          break;

        case UI_KEY_UPPERLIMITMAG:
          ((float *)colarr)[oind] = MKC_MAG(oi[ OCOL_UPPERLIMIT_B ]);
          break;

        case UI_KEY_SN:
          ((float *)colarr)[oind] = columns_sn(p, oi, 0);
          break;

        case UI_KEY_SKY:
          ((float *)colarr)[oind] = MKC_RATIO(oi[OCOL_SUMSKY], oi[OCOL_NUM]);
          break;

        case UI_KEY_STD:
          ((float *)colarr)[oind] = MKC_RATIO(oi[OCOL_SUMSTD], oi[OCOL_NUM]);
          break;

        case UI_KEY_SEMIMAJOR:
          ((float *)colarr)[oind] = columns_second_order(pp, oi, key, 0);
          break;

        case UI_KEY_SEMIMINOR:
          ((float *)colarr)[oind] = columns_second_order(pp, oi, key, 0);
          break;

        case UI_KEY_AXISRATIO:
          ((float *)colarr)[oind]
            = ( columns_second_order(pp, oi, UI_KEY_SEMIMINOR, 0)
                / columns_second_order(pp, oi, UI_KEY_SEMIMAJOR, 0) );
          break;

        case UI_KEY_POSITIONANGLE:
          ((float *)colarr)[oind] = columns_second_order(pp, oi, key, 0);
          break;

        case UI_KEY_GEOSEMIMAJOR:
          ((float *)colarr)[oind] = columns_second_order(pp, oi, key, 0);
          break;

        case UI_KEY_GEOSEMIMINOR:
          ((float *)colarr)[oind] = columns_second_order(pp, oi, key, 0);
          break;

        case UI_KEY_GEOAXISRATIO:
          ((float *)colarr)[oind]
            = ( columns_second_order(pp, oi, UI_KEY_GEOSEMIMINOR, 0)
                / columns_second_order(pp, oi, UI_KEY_GEOSEMIMAJOR, 0) );
          break;

        case UI_KEY_GEOPOSITIONANGLE:
          ((float *)colarr)[oind] = columns_second_order(pp, oi, key, 0);
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                "solve the problem. the output column code %d not recognized "
                "(for objects). ", __func__, PACKAGE_BUGREPORT, key);
        }
    }

  /* Go over the clump columns and fill the information. */
  for(column=p->clumpcols; column!=NULL; column=column->next)
    for(coind=0;coind<pp->clumpsinobj;++coind)
      {
        /* `coind': clump-in-object-index.
           `cind': clump-index (over all the catalog). */
        cind   = sr + coind;
        colarr = column->array;
        key    = column->status;
        ci     = &pp->ci[ coind * CCOL_NUMCOLS ];

        /* Parse columns */
        switch(key)
          {
          case UI_KEY_HOSTOBJID:
            ((int32_t *)colarr)[cind]=pp->object;
            break;

          case UI_KEY_IDINHOSTOBJ:
            ((int32_t *)colarr)[cind]=coind+1;
            break;

          case UI_KEY_AREA:
            ((int32_t *)colarr)[cind]=ci[CCOL_NUM];
            break;

          case UI_KEY_X:
            ((float *)colarr)[cind] = MKC_RATIO( ci[CCOL_VX],
                                                 ci[CCOL_SUMPOS] );
            break;

          case UI_KEY_Y:
            ((float *)colarr)[cind] = MKC_RATIO( ci[CCOL_VY],
                                                 ci[CCOL_SUMPOS] );
            break;

          case UI_KEY_GEOX:
            ((float *)colarr)[cind] = MKC_RATIO( ci[CCOL_GX], ci[CCOL_NUM] );
            break;

          case UI_KEY_GEOY:
            ((float *)colarr)[cind] = MKC_RATIO( ci[CCOL_GY], ci[CCOL_NUM] );
            break;

          case UI_KEY_RA:
          case UI_KEY_DEC:
            p->rd_vc[0][cind] = MKC_RATIO( ci[CCOL_VX], ci[CCOL_SUMPOS] );
            p->rd_vc[1][cind] = MKC_RATIO( ci[CCOL_VY], ci[CCOL_SUMPOS] );
            break;

          case UI_KEY_GEORA:
          case UI_KEY_GEODEC:
            p->rd_gc[0][cind] = MKC_RATIO( ci[CCOL_GX], ci[CCOL_NUM] );
            p->rd_gc[1][cind] = MKC_RATIO( ci[CCOL_GY], ci[CCOL_NUM] );
            break;

          case UI_KEY_BRIGHTNESS:
            /* Calculate the river flux over the clump area. */
            tmp = ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]*ci[ CCOL_NUM ];

            /* Subtract it from the clump's brightness. */
            ((float *)colarr)[cind] = ci[ CCOL_SUM ] - tmp;
            break;

          case UI_KEY_NORIVERBRIGHTNESS:
            ((float *)colarr)[cind] = ci[ CCOL_SUM ];
            break;

          case UI_KEY_MAGNITUDE: /* Similar: brightness for clumps */
            tmp = ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]*ci[ CCOL_NUM ];
            ((float *)colarr)[cind] = MKC_MAG(ci[ CCOL_SUM ]-tmp);
            break;

          case UI_KEY_MAGNITUDEERR:
            ((float *)colarr)[cind] = MAG_ERROR(p, ci, 1);
            break;

          case UI_KEY_UPPERLIMIT:
            ((float *)colarr)[cind] = ci[ CCOL_UPPERLIMIT_B ];
            break;

          case UI_KEY_UPPERLIMITMAG:
            ((float *)colarr)[cind] = MKC_MAG(ci[ CCOL_UPPERLIMIT_B ]);
            break;

          case UI_KEY_RIVERAVE:
            ((float *)colarr)[cind] = ( ci[ CCOL_RIV_SUM]
                                        / ci[ CCOL_RIV_NUM] );
            break;

          case UI_KEY_RIVERNUM:
            ((int32_t *)colarr)[cind] = ci[ CCOL_RIV_NUM ];
            break;

          case UI_KEY_SN:
            ((float *)colarr)[cind] = columns_sn(p, ci, 1);
            break;

          case UI_KEY_SKY:
            ((float *)colarr)[cind] = MKC_RATIO( ci[ CCOL_SUMSKY],
                                                 ci[ CCOL_NUM] );
            break;

          case UI_KEY_STD:
            ((float *)colarr)[cind] = MKC_RATIO( ci[ CCOL_SUMSTD ],
                                                 ci[ CCOL_NUM ] );
            break;

          case UI_KEY_SEMIMAJOR:
            ((float *)colarr)[cind] = columns_second_order(pp, ci, key, 1);
            break;

          case UI_KEY_SEMIMINOR:
            ((float *)colarr)[cind] = columns_second_order(pp, ci, key, 1);
            break;

          case UI_KEY_AXISRATIO:
            ((float *)colarr)[cind]
              = ( columns_second_order(pp, ci, UI_KEY_SEMIMINOR, 1)
                  / columns_second_order(pp, ci, UI_KEY_SEMIMAJOR, 1) );
          break;

          case UI_KEY_POSITIONANGLE:
            ((float *)colarr)[cind] = columns_second_order(pp, ci, key, 1);
            break;

          case UI_KEY_GEOSEMIMAJOR:
            ((float *)colarr)[cind] = columns_second_order(pp, ci, key, 1);
            break;

          case UI_KEY_GEOSEMIMINOR:
            ((float *)colarr)[cind] = columns_second_order(pp, ci, key, 1);
            break;

          case UI_KEY_GEOAXISRATIO:
            ((float *)colarr)[cind]
              = ( columns_second_order(pp, ci, UI_KEY_GEOSEMIMINOR, 1)
                  / columns_second_order(pp, ci, UI_KEY_GEOSEMIMAJOR, 1) );
          break;

          case UI_KEY_GEOPOSITIONANGLE:
            ((float *)colarr)[cind] = columns_second_order(pp, ci, key, 1);
            break;

          default:
            error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                  "solve the problem. The output column code %d not "
                  "recognized (for clumps). ", __func__, PACKAGE_BUGREPORT,
                  key);
          }
      }
}
