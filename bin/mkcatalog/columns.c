/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <string.h>
#include <pthread.h>

#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>

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

   The space that is allocated in 'columns_define_alloc' is for the final
   values that are written in the output file. */
static void
columns_alloc_radec(struct mkcatalogparams *p)
{
  size_t i;

  /* For objects. */
  if(p->wcs_vo==NULL)
    for(i=0;i<p->objects->ndim;++i)
      gal_list_data_add_alloc(&p->wcs_vo, NULL, GAL_TYPE_FLOAT64, 1,
                              &p->numobjects, NULL, 0, p->cp.minmapsize,
                              p->cp.quietmmap, NULL, NULL, NULL);

  /* For clumps */
  if(p->clumps && p->wcs_vc==NULL)
    for(i=0;i<p->objects->ndim;++i)
      gal_list_data_add_alloc(&p->wcs_vc, NULL, GAL_TYPE_FLOAT64, 1,
                              &p->numclumps, NULL, 0, p->cp.minmapsize,
                              p->cp.quietmmap, NULL, NULL, NULL);
}





/* Similar to 'columns_alloc_radec'. */
static void
columns_alloc_georadec(struct mkcatalogparams *p)
{
  size_t i;

  /* For objects. */
  if(p->wcs_go==NULL)
    for(i=0;i<p->objects->ndim;++i)
      gal_list_data_add_alloc(&p->wcs_go, NULL, GAL_TYPE_FLOAT64, 1,
                              &p->numobjects, NULL, 0, p->cp.minmapsize,
                              p->cp.quietmmap, NULL, NULL, NULL);

  /* For clumps */
  if(p->clumps && p->wcs_gc==NULL)
    for(i=0;i<p->objects->ndim;++i)
      gal_list_data_add_alloc(&p->wcs_gc, NULL, GAL_TYPE_FLOAT64, 1,
                              &p->numclumps, NULL, 0, p->cp.minmapsize,
                              p->cp.quietmmap, NULL, NULL, NULL);
}





/* Similar to 'columns_alloc_radec'. */
static void
columns_alloc_clumpsradec(struct mkcatalogparams *p)
{
  size_t i;

  if(p->wcs_vcc==NULL)
    for(i=0;i<p->objects->ndim;++i)
      gal_list_data_add_alloc(&p->wcs_vcc, NULL, GAL_TYPE_FLOAT64, 1,
                              &p->numobjects, NULL, 0, p->cp.minmapsize,
                              p->cp.quietmmap, NULL, NULL, NULL);
}





/* Similar to 'columns_alloc_radec'. */
static void
columns_alloc_clumpsgeoradec(struct mkcatalogparams *p)
{
  size_t i;

  if(p->wcs_gcc==NULL)
    for(i=0;i<p->objects->ndim;++i)
      gal_list_data_add_alloc(&p->wcs_gcc, NULL, GAL_TYPE_FLOAT64, 1,
                              &p->numobjects, NULL, 0, p->cp.minmapsize,
                              p->cp.quietmmap, NULL, NULL, NULL);
}





/* Set pointers to fascilitate filling in the values. */
#define SET_WCS_PREPARE(ARR, LIST, ARRNAME) {                           \
    d=0;                                                                \
    errno=0;                                                            \
    (ARR)=malloc(p->objects->ndim * sizeof (ARR) );                     \
    if( (ARR)==NULL )                                                   \
      error(EXIT_FAILURE, 0, "%s: %zu bytes for %s", __func__,          \
            p->objects->ndim * sizeof (ARR), ARRNAME);                  \
    for(tmp=(LIST);tmp!=NULL;tmp=tmp->next) (ARR)[d++]=tmp->array;      \
  }

static void
columns_set_wcs_pointers(struct mkcatalogparams *p, double ***vo,
                         double ***vc, double ***go, double ***gc,
                         double ***vcc, double ***gcc)
{
  size_t d;
  gal_data_t *tmp;

  if(p->wcs_vo)  SET_WCS_PREPARE(*vo,  p->wcs_vo,  "vo" );
  if(p->wcs_vc)  SET_WCS_PREPARE(*vc,  p->wcs_vc,  "vc" );
  if(p->wcs_go)  SET_WCS_PREPARE(*go,  p->wcs_go,  "go" );
  if(p->wcs_gc)  SET_WCS_PREPARE(*gc,  p->wcs_gc,  "gc" );
  if(p->wcs_vcc) SET_WCS_PREPARE(*vcc, p->wcs_vcc, "vcc");
  if(p->wcs_gcc) SET_WCS_PREPARE(*gcc, p->wcs_gcc, "gcc");
}



















/******************************************************************/
/**********       Column definition/allocation      ***************/
/******************************************************************/
static void
columns_wcs_preparation(struct mkcatalogparams *p)
{
  size_t i;
  gal_list_i32_t *colcode;
  int continue_wcs_check=1;

  /* Make sure a WCS structure is present if we need it. */
  for(colcode=p->columnids; colcode!=NULL; colcode=colcode->next)
    {
      if(continue_wcs_check)
        {
          switch(colcode->v)
            {
            /* High-level. */
            case UI_KEY_RA:
            case UI_KEY_DEC:

            /* Low-level. */
            case UI_KEY_W1:
            case UI_KEY_W2:
            case UI_KEY_GEOW1:
            case UI_KEY_GEOW2:
            case UI_KEY_CLUMPSW1:
            case UI_KEY_CLUMPSW2:
            case UI_KEY_CLUMPSGEOW1:
            case UI_KEY_CLUMPSGEOW2:
              if(p->objects->wcs)
                continue_wcs_check=0;
              else
                error(EXIT_FAILURE, 0, "%s (hdu: %s): no WCS meta-data "
                      "found by WCSLIB. Atleast one of the requested columns "
                      "requires world coordinate system meta-data",
                      p->objectsfile, p->cp.hdu);
              break;
            }
        }
      else
        break;
    }

  /* Convert the high-level WCS columns to low-level ones. */
  for(colcode=p->columnids; colcode!=NULL; colcode=colcode->next)
    switch(colcode->v)
      {
      case UI_KEY_RA:
      case UI_KEY_DEC:
        /* Check all the CTYPES. */
        for(i=0;i<p->objects->ndim;++i)
          if( !strcmp(p->ctype[i], colcode->v==UI_KEY_RA ? "RA" : "DEC") )
            {
              colcode->v = i==0 ? UI_KEY_W1 : UI_KEY_W2;
              break;
            }

        /* Make sure it actually existed. */
        if(i==p->objects->ndim)
          error(EXIT_FAILURE, 0, "%s (hdu: %s): %s not present in any of "
                "the WCS axis types (CTYPE)", p->objectsfile, p->cp.hdu,
                colcode->v==UI_KEY_RA ? "RA" : "DEC");
        break;
      }
}





static void
columns_sanity_check(struct mkcatalogparams *p)
{
  gal_list_i32_t *colcode;

  /* If there is any columns that need WCS, the input image needs to have a
     WCS in its headers. So before anything, we need to check if a WCS is
     present or not. This can't be done after the initial setting of column
     properties because the WCS-related columns use information that is
     based on it (for units and names). */
  columns_wcs_preparation(p);

  /* If sigma-clipping measurements are requested, make sure the necessary
     parameters are provided. */
  for(colcode=p->columnids; colcode!=NULL; colcode=colcode->next)
    switch(colcode->v)
      {
      case UI_KEY_SIGCLIPSTD:
      case UI_KEY_SIGCLIPMEAN:
      case UI_KEY_SIGCLIPNUMBER:
      case UI_KEY_SIGCLIPMEDIAN:
        if(isnan(p->sigmaclip[0]) || isnan(p->sigmaclip[1]))
          error(EXIT_FAILURE, 0, "no sigma-clip defined! When any of the "
                "sigma-clipping columns are requested, it is necessary to "
                "specify the necessary sigma-clipping parameters with the "
                "'--sigmaclip' option");
        break;
      }

  /* Check for dimension-specific columns. */
  switch(p->objects->ndim)
    {
    case 2:
      for(colcode=p->columnids; colcode!=NULL; colcode=colcode->next)
        switch(colcode->v)
          {
          case UI_KEY_AREAXY:
          case UI_KEY_GEOAREAXY:
          case UI_KEY_Z:
          case UI_KEY_MINZ:
          case UI_KEY_MAXZ:
          case UI_KEY_GEOZ:
          case UI_KEY_CLUMPSZ:
          case UI_KEY_CLUMPSGEOZ:
          case UI_KEY_W3:
          case UI_KEY_GEOW3:
          case UI_KEY_CLUMPSW3:
          case UI_KEY_CLUMPSGEOW3:
            error(EXIT_FAILURE, 0, "%s (hdu %s) is a 2D dataset, so columns "
                  "relating to a third dimension cannot be requested",
                  p->objectsfile, p->cp.hdu);
          }
      break;

    case 3:
      for(colcode=p->columnids; colcode!=NULL; colcode=colcode->next)
        switch(colcode->v)
          {
          case UI_KEY_SEMIMAJOR:
          case UI_KEY_SEMIMINOR:
          case UI_KEY_AXISRATIO:
          case UI_KEY_POSITIONANGLE:
          case UI_KEY_GEOSEMIMAJOR:
          case UI_KEY_GEOSEMIMINOR:
          case UI_KEY_GEOAXISRATIO:
          case UI_KEY_GEOPOSITIONANGLE:
            error(EXIT_FAILURE, 0, "columns requiring second moment "
                  "calculations (like semi-major, semi-minor, axis ratio "
                  "and position angle) are not yet supported for 3D "
                  "inputs");
          }
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
            "fix the problem. MakeCatalog should not reach this point "
            "when the input has %zu dimensions", __func__,
            PACKAGE_BUGREPORT, p->objects->ndim);
    }
}





/* Set the necessary parameters for each output column and allocate the
   space necessary to keep the values. */
void
columns_define_alloc(struct mkcatalogparams *p)
{
  gal_list_i32_t *colcode;
  size_t ndim=p->objects->ndim;
  gal_list_str_t *strtmp, *noclumpimg=NULL;
  int disp_fmt=0, disp_width=0, disp_precision=0;
  char *name=NULL, *unit=NULL, *ocomment=NULL, *ccomment=NULL;
  uint8_t otype=GAL_TYPE_INVALID, ctype=GAL_TYPE_INVALID, *oiflag, *ciflag;

  /* Do a sanity check on the columns given the input dataset. */
  columns_sanity_check(p);

  /* Allocate the array for which intermediate parameters are
     necessary. The basic issue is that higher-level calculations require a
     smaller domain of raw measurements. So to avoid having to calculate
     something multiple times, each parameter will flag the intermediate
     parameters it requires in these arrays. */
  oiflag = p->oiflag = gal_pointer_allocate(GAL_TYPE_UINT8, OCOL_NUMCOLS, 1,
                                            __func__, "oiflag");
  ciflag = p->ciflag = gal_pointer_allocate(GAL_TYPE_UINT8, CCOL_NUMCOLS, 1,
                                            __func__, "ciflag");

  /* Allocate the columns. */
  for(colcode=p->columnids; colcode!=NULL; colcode=colcode->next)
    {
      /* Set the column-specific parameters, please follow the same order
         as 'args.h'. IMPORTANT: we want the names to be the same as the
         option names. Note that zero 'disp_' variables will be
         automatically determined.*/
      switch(colcode->v)
        {
        case UI_KEY_OBJID:
          name           = "OBJ_ID";
          unit           = "counter";
          ocomment       = "Object identifier.";
          ccomment       = NULL;
          otype          = GAL_TYPE_INT32;
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
          ocomment       = "Number of non-blank pixels.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          oiflag[ OCOL_NUM ] = ciflag[ CCOL_NUM ] = 1;
          break;

        case UI_KEY_AREAXY:
          name           = "AREAXY";
          unit           = "counter";
          ocomment       = "Projected valued pixels in first two dimensions.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          oiflag[ OCOL_NUMXY ] = ciflag[ CCOL_NUMXY ] = 1;
          break;

        case UI_KEY_CLUMPSAREA:
          name           = "AREA_CLUMPS";
          unit           = "counter";
          ocomment       = "Total number of clump pixels in object.";
          ccomment       = NULL;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          oiflag[ OCOL_C_NUM ] = 1;
          break;

        case UI_KEY_WEIGHTAREA:
          name           = "AREA_WEIGHT";
          unit           = "counter";
          ocomment       = "Area used for flux-weighted positions.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          oiflag[ OCOL_NUMWHT ] = ciflag[ CCOL_NUMWHT ] = 1;
          break;

        case UI_KEY_GEOAREA:
          name           = "AREA_FULL";
          unit           = "counter";
          ocomment       = "Full area of label (irrespective of values).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_GEOAREAXY:
          name           = "AREA_FULL_XY";
          unit           = "counter";
          ocomment       = "Project number in first two dimensions.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          oiflag[ OCOL_NUMALLXY ] = ciflag[ CCOL_NUMALLXY ] = 1;
          break;

        case UI_KEY_X:
          name           = "X";
          unit           = "pixel";
          ocomment       = "Flux weighted center (FITS axis 1).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_VX     ] = ciflag[ CCOL_VX     ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[ OCOL_NUMWHT ] = ciflag[ CCOL_NUMWHT ] = 1;
          break;

        case UI_KEY_Y:
          name           = "Y";
          unit           = "pixel";
          ocomment       = "Flux weighted center (FITS axis 2).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_VY     ] = ciflag[ CCOL_VY     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[ OCOL_NUMWHT ] = ciflag[ CCOL_NUMWHT ] = 1;
          break;

        case UI_KEY_Z:
          name           = "Z";
          unit           = "pixel";
          ocomment       = "Flux weighted center (FITS axis 3).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_VZ ] = 1;
          ciflag[ CCOL_VZ ] = 1;
          break;

        case UI_KEY_GEOX:
          name           = "GEO_X";
          unit           = "pixel";
          ocomment       = "Geometric center (FITS axis 1).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_GEOY:
          name           = "GEO_Y";
          unit           = "pixel";
          ocomment       = "Geometric center (FITS axis 2).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_GEOZ:
          name           = "GEO_Z";
          unit           = "pixel";
          ocomment       = "Geometric center (FITS axis 3).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_GZ ] = 1;
          ciflag[ CCOL_GZ ] = 1;
          break;

        case UI_KEY_CLUMPSX:
          name           = "CLUMPS_X";
          unit           = "pixel";
          ocomment       = "Flux weighted center of clumps (FITS axis 1).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_VX     ] = 1;
          oiflag[ OCOL_C_GX     ] = 1;
          oiflag[ OCOL_C_SUMWHT ] = 1;
          oiflag[ OCOL_C_NUMWHT ] = 1;
          break;

        case UI_KEY_CLUMPSY:
          name           = "CLUMPS_Y";
          unit           = "pixel";
          ocomment       = "Flux weighted center of clumps (FITS axis 2).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_VY     ] = 1;
          oiflag[ OCOL_C_GY     ] = 1;
          oiflag[ OCOL_C_SUMWHT ] = 1;
          oiflag[ OCOL_C_NUMWHT ] = 1;
          break;

        case UI_KEY_CLUMPSZ:
          name           = "CLUMPS_Z";
          unit           = "pixel";
          ocomment       = "Flux weighted center of clumps (FITS axis 3).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_VZ ] = 1;
          break;

        case UI_KEY_CLUMPSGEOX:
          name           = "CLUMPS_GEO_X";
          unit           = "pixel";
          ocomment       = "Geometric center of clumps (FITS axis 1).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_GX     ] = 1;
          oiflag[ OCOL_C_NUMALL ] = 1;
          break;

        case UI_KEY_CLUMPSGEOY:
          name           = "CLUMPS_GEO_Y";
          unit           = "pixel";
          ocomment       = "Geometric center of clumps (FITS axis 2).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_GY     ] = 1;
          oiflag[ OCOL_C_NUMALL ] = 1;
          break;

        case UI_KEY_CLUMPSGEOZ:
          name           = "CLUMPS_GEO_Z";
          unit           = "pixel";
          ocomment       = "Geometric center of clumps (FITS axis 3).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 10;
          disp_precision = 3;
          oiflag[ OCOL_C_GZ ] = 1;

        case UI_KEY_MINX:
          name           = "MIN_X";
          unit           = "pixel";
          ocomment       = "Minimum X axis pixel position.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_UINT32;
          ctype          = GAL_TYPE_UINT32;
          disp_fmt       = 0;
          disp_width     = 10;
          disp_precision = 0;
          ciflag[ CCOL_MINX ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_MAXX:
          name           = "MAX_X";
          unit           = "pixel";
          ocomment       = "Maximum X axis pixel position.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_UINT32;
          ctype          = GAL_TYPE_UINT32;
          disp_fmt       = 0;
          disp_width     = 10;
          disp_precision = 0;
          ciflag[ CCOL_MAXX ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_MINY:
          name           = "MIN_Y";
          unit           = "pixel";
          ocomment       = "Minimum Y axis pixel position.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_UINT32;
          ctype          = GAL_TYPE_UINT32;
          disp_fmt       = 0;
          disp_width     = 10;
          disp_precision = 0;
          ciflag[ CCOL_MINY ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_MAXY:
          name           = "MAX_Y";
          unit           = "pixel";
          ocomment       = "Maximum Y axis pixel position.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_UINT32;
          ctype          = GAL_TYPE_UINT32;
          disp_fmt       = 0;
          disp_width     = 10;
          disp_precision = 0;
          ciflag[ CCOL_MAXY ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_MINZ:
          name           = "MIN_Z";
          unit           = "pixel";
          ocomment       = "Minimum Z axis pixel position.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_UINT32;
          ctype          = GAL_TYPE_UINT32;
          disp_fmt       = 0;
          disp_width     = 10;
          disp_precision = 0;
          ciflag[ CCOL_MINZ ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_MAXZ:
          name           = "MAX_Z";
          unit           = "pixel";
          ocomment       = "Maximum Z axis pixel position.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_UINT32;
          ctype          = GAL_TYPE_UINT32;
          disp_fmt       = 0;
          disp_width     = 10;
          disp_precision = 0;
          ciflag[ CCOL_MAXZ ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_W1:
          name           = p->ctype[0];
          unit           = p->objects->wcs->cunit[0];
          ocomment       = "Flux weighted center (WCS axis 1).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_radec(p);
          oiflag[     OCOL_VX     ] = ciflag[ CCOL_VX     ] = 1;
          oiflag[     OCOL_VY     ] = ciflag[ CCOL_VY     ] = 1;
          oiflag[     OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[     OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[     OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[     OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          if(ndim==3)
            {
              oiflag[ OCOL_VZ ] = ciflag[ CCOL_VZ ] = 1;
              oiflag[ OCOL_GZ ] = ciflag[ CCOL_GZ ] = 1;
            }
          break;

        case UI_KEY_W2:
          name           = p->ctype[1];
          unit           = p->objects->wcs->cunit[1];
          ocomment       = "Flux weighted center (WCS axis 2).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_radec(p);
          oiflag[     OCOL_VX     ] = ciflag[ CCOL_VX     ] = 1;
          oiflag[     OCOL_VY     ] = ciflag[ CCOL_VY     ] = 1;
          oiflag[     OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[     OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[     OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[     OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          if(ndim==3)
            {
              oiflag[ OCOL_VZ     ] = ciflag[ CCOL_VZ     ] = 1;
              oiflag[ OCOL_GZ     ] = ciflag[ CCOL_GZ     ] = 1;
            }
          break;

        case UI_KEY_W3:
          name           = p->ctype[2];
          unit           = p->objects->wcs->cunit[2];
          ocomment       = "Flux weighted center (WCS axis 3).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_radec(p);
          oiflag[ OCOL_VX     ] = ciflag[ CCOL_VX     ] = 1;
          oiflag[ OCOL_VY     ] = ciflag[ CCOL_VY     ] = 1;
          oiflag[ OCOL_VZ     ] = ciflag[ CCOL_VZ     ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GZ     ] = ciflag[ CCOL_GZ     ] = 1;
          oiflag[ OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_GEOW1:
          name           = gal_checkset_malloc_cat("GEO_", p->ctype[0]);
          unit           = p->objects->wcs->cunit[0];
          ocomment       = "Geometric center (WCS axis 1).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_georadec(p);
          oiflag[   OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[   OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[   OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          if(ndim==3)
            oiflag[ OCOL_GZ     ] = ciflag[ CCOL_GZ     ] = 1;
          break;

        case UI_KEY_GEOW2:
          name           = gal_checkset_malloc_cat("GEO_", p->ctype[1]);
          unit           = p->objects->wcs->cunit[1];
          ocomment       = "Geometric center (WCS axis 2).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_georadec(p);
          oiflag[   OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[   OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[   OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          if(ndim==3)
            oiflag[ OCOL_GZ     ] = ciflag[ CCOL_GZ     ] = 1;
          break;

        case UI_KEY_GEOW3:
          name           = gal_checkset_malloc_cat("GEO_", p->ctype[2]);
          unit           = p->objects->wcs->cunit[2];
          ocomment       = "Geometric center (WCS axis 3).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_FLOAT64;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_georadec(p);
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GZ     ] = ciflag[ CCOL_GZ     ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          break;

        case UI_KEY_CLUMPSW1:
          name           = gal_checkset_malloc_cat("CLUMPS_", p->ctype[0]);
          unit           = p->objects->wcs->cunit[0];
          ocomment       = "Flux.wht center of all clumps (WCS axis 1).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_clumpsradec(p);
          oiflag[     OCOL_C_VX     ] = 1;
          oiflag[     OCOL_C_VY     ] = 1;
          oiflag[     OCOL_C_GX     ] = 1;
          oiflag[     OCOL_C_GY     ] = 1;
          oiflag[     OCOL_C_SUMWHT ] = 1;
          oiflag[     OCOL_C_NUMALL ] = 1;
          if(ndim==3)
            {
              oiflag[ OCOL_C_VZ     ] = 1;
              oiflag[ OCOL_C_GZ     ] = 1;
            }
          break;

        case UI_KEY_CLUMPSW2:
          name           = gal_checkset_malloc_cat("CLUMPS_", p->ctype[1]);
          unit           = p->objects->wcs->cunit[1];
          ocomment       = "Flux.wht center of all clumps (WCS axis 2).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 15;
          disp_precision = 7;
          columns_alloc_clumpsradec(p);
          oiflag[     OCOL_C_VX     ] = 1;
          oiflag[     OCOL_C_VY     ] = 1;
          oiflag[     OCOL_C_GX     ] = 1;
          oiflag[     OCOL_C_GY     ] = 1;
          oiflag[     OCOL_C_SUMWHT ] = 1;
          oiflag[     OCOL_C_NUMALL ] = 1;
          if(ndim==3)
            {
              oiflag[ OCOL_C_VZ     ] = 1;
              oiflag[ OCOL_C_GZ     ] = 1;
            }
          break;

        case UI_KEY_CLUMPSW3:
          name           = gal_checkset_malloc_cat("CLUMPS_", p->ctype[2]);
          unit           = p->objects->wcs->cunit[2];
          ocomment       = "Flux.wht center of all clumps (WCS axis 3).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 15;
          disp_precision = 7;
          columns_alloc_clumpsradec(p);
          oiflag[ OCOL_C_VX     ] = 1;
          oiflag[ OCOL_C_VY     ] = 1;
          oiflag[ OCOL_C_VZ     ] = 1;
          oiflag[ OCOL_C_GX     ] = 1;
          oiflag[ OCOL_C_GY     ] = 1;
          oiflag[ OCOL_C_GZ     ] = 1;
          oiflag[ OCOL_C_SUMWHT ] = 1;
          oiflag[ OCOL_C_NUMALL ] = 1;
          if(ndim==3)
            {
            }
          break;

        case UI_KEY_CLUMPSGEOW1:
          name           = gal_checkset_malloc_cat("CLUMPS_GEO", p->ctype[0]);
          unit           = p->objects->wcs->cunit[0];
          ocomment       = "Geometric center of all clumps (WCS axis 1).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_clumpsgeoradec(p);
          oiflag[   OCOL_C_GX     ] = 1;
          oiflag[   OCOL_C_GY     ] = 1;
          oiflag[   OCOL_C_NUMALL ] = 1;
          if(ndim==3)
            oiflag[ OCOL_C_GZ     ] = 1;
          break;

        case UI_KEY_CLUMPSGEOW2:
          name           = gal_checkset_malloc_cat("CLUMPS_GEO", p->ctype[1]);
          unit           = p->objects->wcs->cunit[1];
          ocomment       = "Geometric center of all clumps (WCS axis 2).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_clumpsgeoradec(p);
          oiflag[   OCOL_C_GX     ] = 1;
          oiflag[   OCOL_C_GY     ] = 1;
          oiflag[   OCOL_C_NUMALL ] = 1;
          if(ndim==3)
            oiflag[ OCOL_C_GZ     ] = 1;
          break;

        case UI_KEY_CLUMPSGEOW3:
          name           = gal_checkset_malloc_cat("CLUMPS_GEO", p->ctype[2]);
          unit           = p->objects->wcs->cunit[2];
          ocomment       = "Geometric center of all clumps (WCS axis 3).";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT64;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 13;
          disp_precision = 7;
          columns_alloc_clumpsgeoradec(p);
          oiflag[ OCOL_C_GX     ] = 1;
          oiflag[ OCOL_C_GY     ] = 1;
          oiflag[ OCOL_C_GZ     ] = 1;
          oiflag[ OCOL_C_NUMALL ] = 1;
          break;

        case UI_KEY_BRIGHTNESS:
          name           = "BRIGHTNESS";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Brightness (sum of sky subtracted values).";
          ccomment       = "Brightness (sum of pixels subtracted by rivers).";
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          oiflag[ OCOL_NUM     ] = ciflag[ CCOL_NUM     ] = 1;
          oiflag[ OCOL_SUM     ] = ciflag[ CCOL_SUM     ] = 1;
                                   ciflag[ CCOL_RIV_NUM ] = 1;
                                   ciflag[ CCOL_RIV_SUM ] = 1;
          break;

        case UI_KEY_BRIGHTNESSERR:
          name           = "BRIGHTNESS_ERROR";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Error (1-sigma) in measuring brightness.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          oiflag[ OCOL_NUM     ] = ciflag[ CCOL_NUM         ] = 1;
          oiflag[ OCOL_SUM_VAR ] = ciflag[ CCOL_SUM_VAR     ] = 1;
                                   ciflag[ CCOL_RIV_NUM     ] = 1;
                                   ciflag[ CCOL_RIV_SUM_VAR ] = 1;
          break;

        case UI_KEY_CLUMPSBRIGHTNESS:
          name           = "CLUMPS_BRIGHTNESS";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Brightness (sum of pixel values) in clumps.";
          ccomment       = NULL;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_INVALID;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          oiflag[ OCOL_C_NUM ] = 1;
          oiflag[ OCOL_C_SUM ] = 1;
          break;

        case UI_KEY_BRIGHTNESSNORIVER:
          name           = "NO_RIVER_BRIGHTNESS";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = NULL;
          ccomment       = "Brightness (sum of sky subtracted values).";
          otype          = GAL_TYPE_INVALID;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          ciflag[ CCOL_NUM ] = 1;
          ciflag[ CCOL_SUM ] = 1;
          break;

        case UI_KEY_MEAN:
          name           = "MEAN";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Mean of sky subtracted values.";
          ccomment       = "Mean of pixels subtracted by rivers.";
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          oiflag[ OCOL_NUM     ] = ciflag[ CCOL_NUM     ] = 1;
          oiflag[ OCOL_SUM     ] = ciflag[ CCOL_SUM     ] = 1;
                                   ciflag[ CCOL_RIV_NUM ] = 1;
                                   ciflag[ CCOL_RIV_SUM ] = 1;
          break;

        case UI_KEY_MEDIAN:
          name           = "MEDIAN";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Median of sky subtracted values.";
          ccomment       = "Median of pixels subtracted by rivers.";
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          oiflag[ OCOL_NUM     ] = ciflag[ CCOL_NUM     ] = 1;
          oiflag[ OCOL_MEDIAN  ] = ciflag[ CCOL_MEDIAN  ] = 1;
                                   ciflag[ CCOL_RIV_NUM ] = 1;
                                   ciflag[ CCOL_RIV_SUM ] = 1;
          break;

        case UI_KEY_SIGCLIPNUMBER:
          name           = "SIGCLIP-NUMBER";
          unit           = "counter";
          ocomment       = "Number of pixels in Sigma-clipped object";
          ccomment       = "Number of pixels in Sigma-clipped clump.";
          otype          = GAL_TYPE_INT32;
          ctype          = GAL_TYPE_INT32;
          disp_fmt       = 0;
          disp_width     = 6;
          disp_precision = 0;
          oiflag[ OCOL_NUM        ] = ciflag[ CCOL_NUM        ] = 1;
          oiflag[ OCOL_SIGCLIPNUM ] = ciflag[ CCOL_SIGCLIPNUM ] = 1;
                                      ciflag[ CCOL_RIV_NUM    ] = 1;
                                      ciflag[ CCOL_RIV_SUM    ] = 1;
          break;

        case UI_KEY_SIGCLIPMEDIAN:
          name           = "SIGCLIP-MEDIAN";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Sigma-clipped median of object pixels.";
          ccomment       = "Sigma-clipped median of clump pixels.";
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          oiflag[ OCOL_NUM           ] = ciflag[ CCOL_NUM           ] = 1;
          oiflag[ OCOL_SIGCLIPMEDIAN ] = ciflag[ CCOL_SIGCLIPMEDIAN ] = 1;
                                         ciflag[ CCOL_RIV_NUM       ] = 1;
                                         ciflag[ CCOL_RIV_SUM       ] = 1;
          break;

        case UI_KEY_SIGCLIPMEAN:
          name           = "SIGCLIP-MEAN";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Sigma-clipped mean of object pixels.";
          ccomment       = "Sigma-clipped mean of clump pixels.";
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          oiflag[ OCOL_NUM         ] = ciflag[ CCOL_NUM         ] = 1;
          oiflag[ OCOL_SIGCLIPMEAN ] = ciflag[ CCOL_SIGCLIPMEAN ] = 1;
                                       ciflag[ CCOL_RIV_NUM     ] = 1;
                                       ciflag[ CCOL_RIV_SUM     ] = 1;
          break;

        case UI_KEY_SIGCLIPSTD:
          name           = "SIGCLIP-STD";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Sigma-clipped standard deviation of object pixels.";
          ccomment       = "Sigma-clipped standard deviation of clump pixels.";
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          oiflag[ OCOL_NUM        ] = ciflag[ CCOL_NUM        ] = 1;
          oiflag[ OCOL_SIGCLIPSTD ] = ciflag[ CCOL_SIGCLIPSTD ] = 1;
                                      ciflag[ CCOL_RIV_NUM   ] = 1;
                                      ciflag[ CCOL_RIV_SUM   ] = 1;
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
          p->hasmag      = 1;
          oiflag[ OCOL_NUM     ] = ciflag[ CCOL_NUM     ] = 1;
          oiflag[ OCOL_SUM     ] = ciflag[ CCOL_SUM     ] = 1;
                                   ciflag[ CCOL_SUM     ] = 1;
                                   ciflag[ CCOL_RIV_NUM ] = 1;
                                   ciflag[ CCOL_RIV_SUM ] = 1;
                                   ciflag[ CCOL_RIV_NUM ] = 1;
          break;

        case UI_KEY_MAGNITUDEERR:
          name           = "MAGNITUDE_ERROR";
          unit           = "log";
          ocomment       = "Error in measuring magnitude.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 8;
          disp_precision = 3;
          oiflag[ OCOL_SUM         ] = ciflag[ CCOL_SUM         ] = 1;
          oiflag[ OCOL_SUM_VAR     ] = ciflag[ CCOL_SUM_VAR     ] = 1;
                                       ciflag[ CCOL_RIV_SUM     ] = 1;
                                       ciflag[ CCOL_RIV_SUM_VAR ] = 1;
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
          p->hasmag      = 1;
          oiflag[ OCOL_SUM         ] = ciflag[ CCOL_SUM         ] = 1;
          oiflag[ OCOL_SUM_VAR     ] = ciflag[ CCOL_SUM_VAR     ] = 1;
                                       ciflag[ CCOL_RIV_NUM     ] = 1;
                                       ciflag[ CCOL_RIV_SUM     ] = 1;
                                       ciflag[ CCOL_RIV_SUM_VAR ] = 1;
          break;

        case UI_KEY_UPPERLIMIT:
          name           = "UPPERLIMIT";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Upper limit value (random positionings).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          p->upperlimit  = 1;
          oiflag[ OCOL_UPPERLIMIT_B ] = ciflag[ CCOL_UPPERLIMIT_B ] = 1;

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
          p->hasmag      = 1;
          p->upperlimit  = 1;
          oiflag[ OCOL_UPPERLIMIT_B ] = ciflag[ CCOL_UPPERLIMIT_B ] = 1;
          break;

        case UI_KEY_UPPERLIMITONESIGMA:
          name           = "UPPERLIMIT_ONE_SIGMA";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "One sigma value of all random measurements.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          p->upperlimit  = 1;
          oiflag[ OCOL_UPPERLIMIT_S ] = ciflag[ CCOL_UPPERLIMIT_S ] = 1;
          break;

        case UI_KEY_UPPERLIMITSIGMA:
          name           = "UPPERLIMIT_SIGMA";
          unit           = "frac";
          ocomment       = "Place in upperlimit distribution (sigma "
                           "multiple).";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          p->upperlimit  = 1;
          oiflag[ OCOL_NUM          ] = ciflag[ CCOL_NUM          ] = 1;
          oiflag[ OCOL_SUM          ] = ciflag[ CCOL_SUM          ] = 1;
          oiflag[ OCOL_UPPERLIMIT_S ] = ciflag[ CCOL_UPPERLIMIT_S ] = 1;
                                        ciflag[ CCOL_RIV_NUM      ] = 1;
                                        ciflag[ CCOL_RIV_SUM      ] = 1;
          break;

        case UI_KEY_UPPERLIMITQUANTILE:
          name           = "UPPERLIMIT_QUANTILE";
          unit           = "quantile";
          ocomment       = "Quantile of brightness in random distribution.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          p->upperlimit  = 1;
          oiflag[ OCOL_UPPERLIMIT_Q ] = ciflag[ CCOL_UPPERLIMIT_Q ] = 1;
          break;

        case UI_KEY_UPPERLIMITSKEW:
          name           = "UPPERLIMIT_SKEW";
          unit           = "frac";
          ocomment       = "(Mean-Median)/STD of random distribution.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_FLOAT;
          disp_width     = 8;
          disp_precision = 3;
          p->upperlimit  = 1;
          oiflag[ OCOL_UPPERLIMIT_SKEW ] = oiflag[ CCOL_UPPERLIMIT_SKEW ] = 1;
          break;

        case UI_KEY_RIVERAVE:
          name           = "RIVER_AVE";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = NULL;
          ccomment       = "Average river value surrounding this clump.";
          otype          = GAL_TYPE_INVALID;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 5;
          ciflag[ CCOL_RIV_NUM ] = ciflag[ CCOL_RIV_SUM ] = 1;
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
          oiflag[ OCOL_SUM         ] = ciflag[ CCOL_SUM         ] = 1;
          oiflag[ OCOL_SUM_VAR     ] = ciflag[ CCOL_SUM_VAR     ] = 1;
                                       ciflag[ CCOL_NUM         ] = 1;
                                       ciflag[ CCOL_RIV_NUM     ] = 1;
                                       ciflag[ CCOL_RIV_SUM     ] = 1;
                                       ciflag[ CCOL_RIV_SUM_VAR ] = 1;
          break;

        case UI_KEY_SKY:
          name           = "SKY";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Average input sky value.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 4;
          oiflag[ OCOL_NUM    ] = ciflag[ CCOL_NUM    ] = 1;
          oiflag[ OCOL_SUMSKY ] = ciflag[ CCOL_SUMSKY ] = 1;
          break;

        case UI_KEY_STD:
          name           = "STD";
          unit           = MKCATALOG_NO_UNIT;
          ocomment       = "Sky standard deviation under object.";
          ccomment       = ocomment;
          otype          = GAL_TYPE_FLOAT32;
          ctype          = GAL_TYPE_FLOAT32;
          disp_fmt       = GAL_TABLE_DISPLAY_FMT_GENERAL;
          disp_width     = 10;
          disp_precision = 4;
          oiflag[ OCOL_NUM    ] = ciflag[ CCOL_NUM    ] = 1;
          oiflag[ OCOL_SUMVAR ] = ciflag[ CCOL_SUMVAR ] = 1;
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
          oiflag[ OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[ OCOL_VX     ] = ciflag[ CCOL_VX     ] = 1;
          oiflag[ OCOL_VY     ] = ciflag[ CCOL_VY     ] = 1;
          oiflag[ OCOL_VXX    ] = ciflag[ CCOL_VXX    ] = 1;
          oiflag[ OCOL_VYY    ] = ciflag[ CCOL_VYY    ] = 1;
          oiflag[ OCOL_VXY    ] = ciflag[ CCOL_VXY    ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GXX    ] = ciflag[ CCOL_GXX    ] = 1;
          oiflag[ OCOL_GYY    ] = ciflag[ CCOL_GYY    ] = 1;
          oiflag[ OCOL_GXY    ] = ciflag[ CCOL_GXY    ] = 1;
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
          oiflag[ OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[ OCOL_VX     ] = ciflag[ CCOL_VX     ] = 1;
          oiflag[ OCOL_VY     ] = ciflag[ CCOL_VY     ] = 1;
          oiflag[ OCOL_VXX    ] = ciflag[ CCOL_VXX    ] = 1;
          oiflag[ OCOL_VYY    ] = ciflag[ CCOL_VYY    ] = 1;
          oiflag[ OCOL_VXY    ] = ciflag[ CCOL_VXY    ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GXX    ] = ciflag[ CCOL_GXX    ] = 1;
          oiflag[ OCOL_GYY    ] = ciflag[ CCOL_GYY    ] = 1;
          oiflag[ OCOL_GXY    ] = ciflag[ CCOL_GXY    ] = 1;
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
          oiflag[ OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[ OCOL_VX     ] = ciflag[ CCOL_VX     ] = 1;
          oiflag[ OCOL_VY     ] = ciflag[ CCOL_VY     ] = 1;
          oiflag[ OCOL_VXX    ] = ciflag[ CCOL_VXX    ] = 1;
          oiflag[ OCOL_VYY    ] = ciflag[ CCOL_VYY    ] = 1;
          oiflag[ OCOL_VXY    ] = ciflag[ CCOL_VXY    ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GXX    ] = ciflag[ CCOL_GXX    ] = 1;
          oiflag[ OCOL_GYY    ] = ciflag[ CCOL_GYY    ] = 1;
          oiflag[ OCOL_GXY    ] = ciflag[ CCOL_GXY    ] = 1;
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
          oiflag[ OCOL_SUMWHT ] = ciflag[ CCOL_SUMWHT ] = 1;
          oiflag[ OCOL_VX     ] = ciflag[ CCOL_VX     ] = 1;
          oiflag[ OCOL_VY     ] = ciflag[ CCOL_VY     ] = 1;
          oiflag[ OCOL_VXX    ] = ciflag[ CCOL_VXX    ] = 1;
          oiflag[ OCOL_VYY    ] = ciflag[ CCOL_VYY    ] = 1;
          oiflag[ OCOL_VXY    ] = ciflag[ CCOL_VXY    ] = 1;
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GXX    ] = ciflag[ CCOL_GXX    ] = 1;
          oiflag[ OCOL_GYY    ] = ciflag[ CCOL_GYY    ] = 1;
          oiflag[ OCOL_GXY    ] = ciflag[ CCOL_GXY    ] = 1;
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
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GXX    ] = ciflag[ CCOL_GXX    ] = 1;
          oiflag[ OCOL_GYY    ] = ciflag[ CCOL_GYY    ] = 1;
          oiflag[ OCOL_GXY    ] = ciflag[ CCOL_GXY    ] = 1;
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
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GXX    ] = ciflag[ CCOL_GXX    ] = 1;
          oiflag[ OCOL_GYY    ] = ciflag[ CCOL_GYY    ] = 1;
          oiflag[ OCOL_GXY    ] = ciflag[ CCOL_GXY    ] = 1;
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
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GXX    ] = ciflag[ CCOL_GXX    ] = 1;
          oiflag[ OCOL_GYY    ] = ciflag[ CCOL_GYY    ] = 1;
          oiflag[ OCOL_GXY    ] = ciflag[ CCOL_GXY    ] = 1;
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
          oiflag[ OCOL_NUMALL ] = ciflag[ CCOL_NUMALL ] = 1;
          oiflag[ OCOL_GX     ] = ciflag[ CCOL_GX     ] = 1;
          oiflag[ OCOL_GY     ] = ciflag[ CCOL_GY     ] = 1;
          oiflag[ OCOL_GXX    ] = ciflag[ CCOL_GXX    ] = 1;
          oiflag[ OCOL_GYY    ] = ciflag[ CCOL_GYY    ] = 1;
          oiflag[ OCOL_GXY    ] = ciflag[ CCOL_GXY    ] = 1;
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to fix "
                "the problem. The code %d is not an internally recognized "
                "column code", __func__, PACKAGE_BUGREPORT, colcode->v);
        }


      /* If this is an objects column, add it to the list of columns. We
         will be using the 'status' element to keep the MakeCatalog code
         for the columns. */
      if(otype!=GAL_TYPE_INVALID)
        {
          gal_list_data_add_alloc(&p->objectcols, NULL, otype, 1,
                                  &p->numobjects, NULL, 0, p->cp.minmapsize,
                                  p->cp.quietmmap, name, unit, ocomment);
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
          /* If a clumps labeled image, add this column for the output. */
          if(p->clumps)
            {
              gal_list_data_add_alloc(&p->clumpcols, NULL, ctype, 1,
                                      &p->numclumps, NULL, 0,
                                      p->cp.minmapsize, p->cp.quietmmap,
                                      name, unit, ccomment);
              p->clumpcols->status         = colcode->v;
              p->clumpcols->disp_fmt       = disp_fmt;
              p->clumpcols->disp_width     = disp_width;
              p->clumpcols->disp_precision = disp_precision;
            }


          /* If this is a clumps-only column and no clumps image was
             given. Add the column to the list of similar columns to inform
             the user.

             We'll just ignore the clump-specific ID-related columns,
             because the '--ids' (generic for both objects and clumps) is a
             simple generic solution for identifiers. Also, ultimately,
             they aren't measurements. */
          else if( otype==GAL_TYPE_INVALID
                   && colcode->v!=UI_KEY_HOSTOBJID
                   && colcode->v!=UI_KEY_IDINHOSTOBJ )
            gal_list_str_add(&noclumpimg, name, 1);
        }
    }


  /* If the user has asked for clump-only columns, but no clumps catalog is
     to be created (the '--clumpscat' option was not given or there were no
     clumps in the specified image), then print an informative message that
     the columns in question will be ignored. */
  if(noclumpimg)
    {
      gal_list_str_reverse(&noclumpimg);
      fprintf(stderr, "WARNING: the following column(s) are unique to "
              "clumps (not objects), but the '--clumpscat' option has not "
              "been called, or there were no clumps in the clumps labeled "
              "image. Hence, these columns will be ignored in the "
              "output.\n\n");
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





/* Calculate the error in brightness. */
static double
columns_brightness_error(struct mkcatalogparams *p, double *row, int o0c1)
{
  double V = row[ o0c1 ? CCOL_SUM_VAR : OCOL_SUM_VAR ];
  double OV = (o0c1 && row[ CCOL_RIV_NUM ]) ? row[ CCOL_RIV_SUM_VAR ] : 0.0;
  return sqrt(V+OV);
}





/* Calculate the Signal to noise ratio for the object. */
static double
columns_sn(struct mkcatalogparams *p, double *row, int o0c1)
{
  double I = row[ o0c1 ? CCOL_SUM     : OCOL_SUM     ];

  /* When grown clumps are requested from NoiseChisel, some "clumps" will
     completely cover their objects and there will be no rivers. So if this
     is a clump, and the river area is 0, we should treat the S/N as a an
     object (and set the outer flux to 0.0). */
  double O = ( (o0c1 && row[ CCOL_RIV_NUM ])
               ? (row[ CCOL_NUM ]*row[ CCOL_RIV_SUM ]/row[ CCOL_RIV_NUM ])
               : 0.0 );

  /* Return the derived value. */
  return sqrt(1/p->cpscorr) * (I-O) / columns_brightness_error(p, row, o0c1);
}





/* Do the second order calculations, see "Measuring elliptical parameters"
   section of the book/manual for a thorough explanation of the
   derivation. */
static double
columns_second_order(struct mkcatalog_passparams *pp, double *row,
                     int key, int o0c1)
{
  double x=NAN, y=NAN, xx=NAN, yy=NAN, xy=NAN;
  double denom, kx=pp->shift[1], ky=pp->shift[0];

  /* Preparations. */
  switch(key)
    {
    /* Brightness weighted. */
    case UI_KEY_SEMIMAJOR:
    case UI_KEY_SEMIMINOR:
    case UI_KEY_POSITIONANGLE:

      /* Denominator (to be divided). */
      denom = row[ o0c1 ? CCOL_SUMWHT : OCOL_SUMWHT ];

      /* First order. */
      x  = MKC_RATIO( row[ o0c1 ? CCOL_VX : OCOL_VX ], denom );
      y  = MKC_RATIO( row[ o0c1 ? CCOL_VY : OCOL_VY ], denom );

      /* Second order. */
      xx = ( MKC_RATIO( row[ o0c1 ? CCOL_VXX : OCOL_VXX ], denom )
             - (x-kx) * (x-kx) );
      yy = ( MKC_RATIO( row[ o0c1 ? CCOL_VYY : OCOL_VYY ], denom )
             - (y-ky) * (y-ky) );
      xy = ( MKC_RATIO( row[ o0c1 ? CCOL_VXY : OCOL_VXY ], denom )
             - (x-kx) * (y-ky) );
      break;

    /* Geometric. */
    case UI_KEY_GEOSEMIMAJOR:
    case UI_KEY_GEOSEMIMINOR:
    case UI_KEY_GEOPOSITIONANGLE:

      /* Denominator (to be divided). */
      denom = row[ o0c1 ? CCOL_NUMALL : OCOL_NUMALL ];

      /* First order. */
      x  = MKC_RATIO( row[ o0c1 ? CCOL_GX : OCOL_GX  ], denom );
      y  = MKC_RATIO( row[ o0c1 ? CCOL_GY : OCOL_GY  ], denom );

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
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
            "problem. %d is not a recognized key", __func__, PACKAGE_BUGREPORT,
            key);
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





/* The clump brightness is needed in several places, so we've defined this
   function to have an easier code. */
static double
columns_clump_brightness(double *ci)
{
  double tmp;
  /* Calculate the river flux over the clump area. But only when rivers are
     present. When grown clumps are requested, the clumps can fully cover a
     detection (that has one or no clumps). */
  tmp = ( ci[ CCOL_RIV_NUM ]>0.0f
          ? ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]*ci[ CCOL_NUM ]
          : 0 );

  /* Subtract it from the clump's brightness. */
  return ci[ CCOL_NUM ]>0.0f ? (ci[ CCOL_SUM ] - tmp) : NAN;
}





/* Measure the minimum and maximum positions. */
static uint32_t
columns_xy_extrema(struct mkcatalog_passparams *pp, double *oi, size_t *coord, int key)
{
  size_t ndim=pp->tile->ndim;
  gal_data_t *tile=pp->tile, *block=tile->block;

  /* We only want to do the coordinate estimation once: in 'columns_fill',
     we initialized the coordinates with 'GAL_BLANK_SIZE_T'. When the
     coordinate has already been measured already, it won't have this value
     any more. */
  if(coord[0]==GAL_BLANK_SIZE_T)
    gal_dimension_index_to_coord(gal_pointer_num_between(block->array,
                                                         tile->array,
                                                         block->type),
                                 block->ndim, block->dsize, coord);

  /* Return the proper value: note that 'coord' is in C standard: starting
     from the slowest dimension and counting from zero. */
  if(oi[OCOL_NUMALL])
    switch(key)
      {
      case UI_KEY_MINX: return coord[ndim-1] + 1;                   break;
      case UI_KEY_MAXX: return coord[ndim-1] + tile->dsize[ndim-1]; break;
      case UI_KEY_MINY: return coord[ndim-2] + 1;                   break;
      case UI_KEY_MAXY: return coord[ndim-2] + tile->dsize[ndim-2]; break;
      case UI_KEY_MINZ: return coord[ndim-3] + 1;                   break;
      case UI_KEY_MAXZ: return coord[ndim-3] + tile->dsize[ndim-3]; break;
      default:
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
              "problem. The value %d is not a recognized value", __func__,
              PACKAGE_BUGREPORT, key);
      }
  else
    return 0;


  /* Control should not reach here. */
  error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return GAL_BLANK_UINT32;
}





/* The magnitude error is directly derivable from the S/N:

   To derive the error in measuring the magnitude from the S/N, let's take
   'F' as the flux, 'Z' is the zeropoint, 'M' is the magnitude, 'S' is the
   S/N, and 'D' to stand for capital delta (or error in a value) then from

      'M = -2.5*log10(F) + Z'

   we get the following equation after calculating the derivative with
   respect to F.

      'dM/df = -2.5 * ( 1 / ( F * ln(10) ) )'

   From the Tailor series, 'DM' can be written as:

      'DM = dM/dF * DF'

   So

      'DM = |-2.5/ln(10)| * DF/F'

   But 'DF/F' is just the inverse of the Signal to noise ratio, or
  '1/S'. So

      'DM = 2.5 / ( S * ln(10) )'               */
#define MAG_ERROR(P,ROW,O0C1) ( 2.5f                                    \
                                / ( ( columns_sn((P),(ROW),(O0C1)) > 0  \
                                      ? columns_sn((P),(ROW),(O0C1))    \
                                      : NAN )                           \
                                    * log(10) ) )






/* All the raw first and second pass information has been collected, now
   write them into the output columns. The list of columns here is in the
   same order as 'columns_alloc_set_out_cols', see there for the type of
   each column. */
#define POS_V_G(ARRAY, SUMWHT_COL, NUMALL_COL, V_COL, G_COL)            \
  ( (ARRAY)[ SUMWHT_COL ]>0                                             \
    ? MKC_RATIO( (ARRAY)[ V_COL ], (ARRAY)[ SUMWHT_COL ] )              \
    : MKC_RATIO( (ARRAY)[ G_COL ], (ARRAY)[ NUMALL_COL ] ) )
void
columns_fill(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;

  int key;
  double tmp;
  void *colarr;
  gal_data_t *column;
  double *ci, *oi=pp->oi;
  size_t coord[3]={GAL_BLANK_SIZE_T, GAL_BLANK_SIZE_T, GAL_BLANK_SIZE_T};

  size_t i, cind, coind, sr=pp->clumpstartindex, oind=GAL_BLANK_SIZE_T;
  double **vo=NULL, **vc=NULL, **go=NULL, **gc=NULL, **vcc=NULL, **gcc=NULL;

  /* Find the object's index in final catalog. */
  if(p->outlabs)
    {
      for(i=0;i<p->numobjects;++i)
        if(p->outlabs[i]==pp->object)
          { oind=i; break; }
    }
  else oind=pp->object-1;

  /* If a WCS column is requested (check will be done inside the function),
     then set the pointers. */
  columns_set_wcs_pointers(p, &vo, &vc, &go, &gc, &vcc, &gcc);

  /* Go over all the object columns and fill in the information. */
  for(column=p->objectcols; column!=NULL; column=column->next)
    {
      /* For easy reading. */
      key=column->status;
      colarr=column->array;

      /* Put the number of clumps in the internal array which we will need
         later to order the clump table by object ID. */
      if(p->numclumps_c)
        p->numclumps_c[oind]=pp->clumpsinobj;

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

        case UI_KEY_AREAXY:
          ((int32_t *)colarr)[oind] = oi[OCOL_NUMXY];
          break;

        case UI_KEY_CLUMPSAREA:
          ((int32_t *)colarr)[oind] = oi[OCOL_C_NUM];
          break;

        case UI_KEY_WEIGHTAREA:
          ((int32_t *)colarr)[oind] = oi[OCOL_NUMWHT];
          break;

        case UI_KEY_GEOAREA:
          ((int32_t *)colarr)[oind] = oi[OCOL_NUMALL];
          break;

        case UI_KEY_GEOAREAXY:
          ((int32_t *)colarr)[oind] = oi[OCOL_NUMALLXY];
          break;

        case UI_KEY_X:
          ((float *)colarr)[oind] = POS_V_G(oi, OCOL_SUMWHT, OCOL_NUMWHT,
                                            OCOL_VX, OCOL_GX);
          break;

        case UI_KEY_Y:
          ((float *)colarr)[oind] = POS_V_G(oi, OCOL_SUMWHT, OCOL_NUMWHT,
                                            OCOL_VY, OCOL_GY);
          break;

        case UI_KEY_Z:
          ((float *)colarr)[oind] = POS_V_G(oi, OCOL_SUMWHT, OCOL_NUMALL,
                                            OCOL_VZ, OCOL_GZ);
          break;

        case UI_KEY_GEOX:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_GX], oi[OCOL_NUMALL] );
          break;

        case UI_KEY_GEOY:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_GY], oi[OCOL_NUMALL] );
          break;

        case UI_KEY_GEOZ:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_GZ], oi[OCOL_NUMALL] );
          break;

        case UI_KEY_CLUMPSX:
          ((float *)colarr)[oind] = POS_V_G(oi, OCOL_C_SUMWHT, OCOL_C_NUMWHT,
                                            OCOL_C_VX, OCOL_C_GX);
          break;

        case UI_KEY_CLUMPSY:
          ((float *)colarr)[oind] = POS_V_G(oi, OCOL_C_SUMWHT, OCOL_C_NUMWHT,
                                            OCOL_C_VY, OCOL_C_GY);
          break;

        case UI_KEY_CLUMPSZ:
          ((float *)colarr)[oind] = POS_V_G(oi, OCOL_C_SUMWHT, OCOL_C_NUMALL,
                                            OCOL_C_VZ, OCOL_C_GZ);
          break;

        case UI_KEY_CLUMPSGEOX:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_C_GX],
                                               oi[OCOL_C_NUMALL] );
          break;

        case UI_KEY_CLUMPSGEOY:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_C_GY],
                                               oi[OCOL_C_NUMALL] );
          break;

        case UI_KEY_CLUMPSGEOZ:
          ((float *)colarr)[oind] = MKC_RATIO( oi[OCOL_C_GZ],
                                               oi[OCOL_C_NUMALL] );

        case UI_KEY_MINX:
        case UI_KEY_MAXX:
        case UI_KEY_MINY:
        case UI_KEY_MAXY:
        case UI_KEY_MINZ:
        case UI_KEY_MAXZ:
          ((uint32_t *)colarr)[oind]=columns_xy_extrema(pp, oi, coord, key);
          break;

        case UI_KEY_W1:
        case UI_KEY_W2:
        case UI_KEY_W3:
          vo[0][oind] = POS_V_G(oi, OCOL_SUMWHT, OCOL_NUMALL, OCOL_VX,
                                OCOL_GX);
          vo[1][oind] = POS_V_G(oi, OCOL_SUMWHT, OCOL_NUMALL, OCOL_VY,
                                OCOL_GY);
          if(p->objects->ndim==3)
            vo[2][oind] = POS_V_G(oi, OCOL_SUMWHT, OCOL_NUMALL, OCOL_VZ,
                                  OCOL_GZ);
          break;

        case UI_KEY_GEOW1:
        case UI_KEY_GEOW2:
        case UI_KEY_GEOW3:
          go[0][oind] = MKC_RATIO( oi[OCOL_GX], oi[OCOL_NUMALL] );
          go[1][oind] = MKC_RATIO( oi[OCOL_GY], oi[OCOL_NUMALL] );
          if(p->objects->ndim==3)
            go[2][oind] = MKC_RATIO( oi[OCOL_GZ], oi[OCOL_NUMALL] );
          break;

        case UI_KEY_CLUMPSW1:
        case UI_KEY_CLUMPSW2:
        case UI_KEY_CLUMPSW3:
          vcc[0][oind] = POS_V_G(oi, OCOL_C_SUMWHT, OCOL_C_NUMALL, OCOL_C_VX,
                                 OCOL_C_GX);
          vcc[1][oind] = POS_V_G(oi, OCOL_C_SUMWHT, OCOL_C_NUMALL, OCOL_C_VY,
                                 OCOL_C_GY);
          if(p->objects->ndim==3)
            vcc[2][oind] = POS_V_G(oi, OCOL_C_SUMWHT, OCOL_C_NUMALL,
                                   OCOL_C_VZ, OCOL_C_GZ);
          break;

        case UI_KEY_CLUMPSGEOW1:
        case UI_KEY_CLUMPSGEOW2:
        case UI_KEY_CLUMPSGEOW3:
          gcc[0][oind] = MKC_RATIO( oi[OCOL_C_GX], oi[OCOL_C_NUMALL] );
          gcc[1][oind] = MKC_RATIO( oi[OCOL_C_GY], oi[OCOL_C_NUMALL] );
          if(p->objects->ndim==3)
            gcc[2][oind] = MKC_RATIO( oi[OCOL_C_GZ], oi[OCOL_C_NUMALL] );
          break;

        case UI_KEY_BRIGHTNESS:
          ((float *)colarr)[oind] = ( oi[ OCOL_NUM ]>0.0f
                                      ? oi[ OCOL_SUM ]
                                      : NAN );
          break;

        case UI_KEY_BRIGHTNESSERR:
          ((float *)colarr)[oind] = ( oi[ OCOL_NUM ]>0.0f
                                      ? columns_brightness_error(p, oi, 0)
                                      : NAN );
          break;

        case UI_KEY_CLUMPSBRIGHTNESS:
          ((float *)colarr)[oind] = ( oi[ OCOL_C_NUM ]>0.0f
                                      ? oi[ OCOL_C_SUM ]
                                      : NAN );
          break;

        case UI_KEY_MEAN:
          ((float *)colarr)[oind] = ( oi[ OCOL_NUM ]>0.0f
                                      ? oi[ OCOL_SUM ] / oi[ OCOL_NUM ]
                                      : NAN );
          break;

        case UI_KEY_MEDIAN:
          ((float *)colarr)[oind] = ( oi[ OCOL_NUM ]>0.0f
                                      ? oi[ OCOL_MEDIAN ]
                                      : NAN );
          break;

        case UI_KEY_SIGCLIPNUMBER:
          ((int32_t *)colarr)[oind] = oi[ OCOL_SIGCLIPNUM ];
          break;

        case UI_KEY_SIGCLIPMEDIAN:
          ((float *)colarr)[oind] = oi[ OCOL_SIGCLIPMEDIAN ];
          break;

        case UI_KEY_SIGCLIPMEAN:
          ((float *)colarr)[oind] = oi[ OCOL_SIGCLIPMEAN ];
          break;

        case UI_KEY_SIGCLIPSTD:
          ((float *)colarr)[oind] = oi[ OCOL_SIGCLIPSTD ];
          break;

        case UI_KEY_MAGNITUDE:
          ((float *)colarr)[oind] = ( oi[ OCOL_NUM ]>0.0f
                                      ? MKC_MAG(oi[ OCOL_SUM ])
                                      : NAN );
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

        case UI_KEY_UPPERLIMITONESIGMA:
          ((float *)colarr)[oind] = oi[ OCOL_UPPERLIMIT_S ];
          break;

        case UI_KEY_UPPERLIMITSIGMA:
          ((float *)colarr)[oind] = ( ( oi[ OCOL_NUM ]>0.0f
                                        ? oi[ OCOL_SUM ] : NAN )
                                      / oi[ OCOL_UPPERLIMIT_S ] );
          break;

        case UI_KEY_UPPERLIMITQUANTILE:
          ((float *)colarr)[oind] = oi[ OCOL_UPPERLIMIT_Q ];
          break;

        case UI_KEY_UPPERLIMITSKEW:
          ((float *)colarr)[oind] = oi[ OCOL_UPPERLIMIT_SKEW ];
          break;

        case UI_KEY_SN:
          ((float *)colarr)[oind] = columns_sn(p, oi, 0);
          break;

        case UI_KEY_SKY:
          ((float *)colarr)[oind] = MKC_RATIO(oi[OCOL_SUMSKY], oi[OCOL_NUMSKY]);
          break;

        case UI_KEY_STD:
          ((float *)colarr)[oind] = sqrt( MKC_RATIO( oi[ OCOL_SUMVAR ],
                                                     oi[ OCOL_NUMVAR ]) );
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
        /* 'coind': clump-in-object-index.
           'cind': clump-index (over all the catalog). */
        cind   = sr + coind;
        colarr = column->array;
        key    = column->status;
        ci     = &pp->ci[ coind * CCOL_NUMCOLS ];

        /* Put the object ID of this clump into the temporary array that we
           will later need to sort the final clumps catalog. */
        if(p->hostobjid_c)
          p->hostobjid_c[cind]=pp->object;

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

          case UI_KEY_AREAXY:
            ((int32_t *)colarr)[cind]=ci[CCOL_NUMXY];
            break;

          case UI_KEY_WEIGHTAREA:
            ((int32_t *)colarr)[cind]=ci[CCOL_NUMWHT];
            break;

          case UI_KEY_GEOAREA:
            ((int32_t *)colarr)[cind]=ci[CCOL_NUMALL];
            break;

          case UI_KEY_GEOAREAXY:
            ((int32_t *)colarr)[cind]=ci[CCOL_NUMALLXY];
            break;

          case UI_KEY_X:
            ((float *)colarr)[cind] = POS_V_G(ci, CCOL_SUMWHT, CCOL_NUMALL,
                                              CCOL_VX, CCOL_GX);
            break;

          case UI_KEY_Y:
            ((float *)colarr)[cind] = POS_V_G(ci, CCOL_SUMWHT, CCOL_NUMALL,
                                              CCOL_VY, CCOL_GY);
            break;

          case UI_KEY_Z:
            ((float *)colarr)[cind] = POS_V_G(ci, CCOL_SUMWHT, CCOL_NUMALL,
                                              CCOL_VZ, CCOL_GZ);
            break;

          case UI_KEY_GEOX:
            ((float *)colarr)[cind] = MKC_RATIO( ci[CCOL_GX],
                                                 ci[CCOL_NUMALL] );
            break;

          case UI_KEY_GEOY:
            ((float *)colarr)[cind] = MKC_RATIO( ci[CCOL_GY],
                                                 ci[CCOL_NUMALL] );
            break;

          case UI_KEY_GEOZ:
            ((float *)colarr)[cind] = MKC_RATIO( ci[CCOL_GZ],
                                                 ci[CCOL_NUMALL] );

          case UI_KEY_MINX:
            ((uint32_t *)colarr)[cind] = ci[CCOL_MINX];
            break;

          case UI_KEY_MAXX:
            ((uint32_t *)colarr)[cind] = ci[CCOL_MAXX];
            break;

          case UI_KEY_MINY:
            ((uint32_t *)colarr)[cind] = ci[CCOL_MINY];
            break;

          case UI_KEY_MAXY:
            ((uint32_t *)colarr)[cind] = ci[CCOL_MAXY];
            break;

          case UI_KEY_MINZ:
            ((uint32_t *)colarr)[cind] = ci[CCOL_MINZ];
            break;

          case UI_KEY_MAXZ:
            ((uint32_t *)colarr)[cind] = ci[CCOL_MAXZ];
            break;

          case UI_KEY_W1:
          case UI_KEY_W2:
          case UI_KEY_W3:
            vc[0][cind] = POS_V_G(ci, CCOL_SUMWHT, CCOL_NUMALL, CCOL_VX,
                                  CCOL_GX);
            vc[1][cind] = POS_V_G(ci, CCOL_SUMWHT, CCOL_NUMALL, CCOL_VY,
                                  CCOL_GY);
            if(p->objects->ndim==3)
              vc[2][cind] = POS_V_G(ci, CCOL_SUMWHT, CCOL_NUMALL, CCOL_VZ,
                                    CCOL_GZ);
            break;

          case UI_KEY_GEOW1:
          case UI_KEY_GEOW2:
          case UI_KEY_GEOW3:
            gc[0][cind] = MKC_RATIO( ci[CCOL_GX], ci[CCOL_NUMALL] );
            gc[1][cind] = MKC_RATIO( ci[CCOL_GY], ci[CCOL_NUMALL] );
            if(p->objects->ndim==3)
              gc[2][cind] = MKC_RATIO( ci[CCOL_GZ], ci[CCOL_NUMALL] );
            break;

          case UI_KEY_BRIGHTNESS:
            ((float *)colarr)[cind] = columns_clump_brightness(ci);
            break;

          case UI_KEY_BRIGHTNESSERR:
            ((float *)colarr)[cind] = ( ci[ CCOL_NUM ]>0.0f
                                        ? columns_brightness_error(p, ci, 1)
                                        : NAN );
            break;

          case UI_KEY_BRIGHTNESSNORIVER:
            ((float *)colarr)[cind] = ( ci[ CCOL_NUM ]>0.0f
                                        ? ci[ CCOL_SUM ] : NAN );
            break;

          case UI_KEY_MEAN:
            ((float *)colarr)[cind] = ( columns_clump_brightness(ci)
                                        /ci[CCOL_NUM] );
            break;

          case UI_KEY_MEDIAN:
            ((float *)colarr)[cind] = ( ci[ CCOL_NUM ]>0.0f
                                        ? ci[ CCOL_MEDIAN ] : NAN );
            break;

          case UI_KEY_SIGCLIPNUMBER:
            ((int32_t *)colarr)[cind] = ci[ CCOL_SIGCLIPNUM ];
            break;

          case UI_KEY_SIGCLIPMEDIAN:
            ((float *)colarr)[cind] = ci[ CCOL_SIGCLIPMEDIAN ];
            break;

          case UI_KEY_SIGCLIPMEAN:
            ((float *)colarr)[cind] = ci[ CCOL_SIGCLIPMEAN ];
            break;

          case UI_KEY_SIGCLIPSTD:
            ((float *)colarr)[cind] = ci[ CCOL_SIGCLIPSTD ];
            break;

          case UI_KEY_MAGNITUDE: /* Similar: brightness for clumps */
            tmp = ( ci[ CCOL_RIV_NUM ]
                    ? ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM ]*ci[ CCOL_NUM ]
                    : 0 );
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

          case UI_KEY_UPPERLIMITONESIGMA:
            ((float *)colarr)[cind] = ci[ CCOL_UPPERLIMIT_S ];
            break;

          case UI_KEY_UPPERLIMITSIGMA:
            ((float *)colarr)[cind] = ( columns_clump_brightness(ci)
                                        / ci[ CCOL_UPPERLIMIT_S ] );
            break;

          case UI_KEY_UPPERLIMITQUANTILE:
            ((float *)colarr)[cind] = ci[ CCOL_UPPERLIMIT_Q ];
            break;

          case UI_KEY_UPPERLIMITSKEW:
            ((float *)colarr)[cind] = ci[ CCOL_UPPERLIMIT_SKEW ];
            break;

          case UI_KEY_RIVERAVE:
            ((float *)colarr)[cind] = ( ci[ CCOL_RIV_NUM]
                                        ? ci[ CCOL_RIV_SUM ]/ci[ CCOL_RIV_NUM]
                                        : NAN );
            break;

          case UI_KEY_RIVERNUM:
            ((int32_t *)colarr)[cind] = ci[ CCOL_RIV_NUM ];
            break;

          case UI_KEY_SN:
            ((float *)colarr)[cind] = columns_sn(p, ci, 1);
            break;

          case UI_KEY_SKY:
            ((float *)colarr)[cind] = MKC_RATIO( ci[ CCOL_SUMSKY],
                                                 ci[ CCOL_NUMSKY] );
            break;

          case UI_KEY_STD:
            ((float *)colarr)[cind] = sqrt( MKC_RATIO( ci[ CCOL_SUMVAR ],
                                                       ci[ CCOL_NUMVAR ] ));
            break;

          case UI_KEY_SEMIMAJOR:
            ((float *)colarr)[cind] = columns_second_order(pp, ci, key, 1);
            break;

          case UI_KEY_SEMIMINOR:
            ((float *)colarr)[cind] = columns_second_order(pp, ci, key, 1);
            break;

          case UI_KEY_AXISRATIO:
            ((float *)colarr)[cind]
              = ( columns_second_order(  pp, ci, UI_KEY_SEMIMINOR, 1)
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
              = ( columns_second_order(  pp, ci, UI_KEY_GEOSEMIMINOR, 1)
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

  /* Clean up. */
  if(vo)  free(vo);
  if(vc)  free(vc);
  if(go)  free(go);
  if(gc)  free(gc);
  if(vcc) free(vcc);
  if(gcc) free(gcc);
}
