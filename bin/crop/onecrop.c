/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/box.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/polygon.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "onecrop.h"
#include "wcsmode.h"























/*******************************************************************/
/************     Set/correct first and last pixel    **************/
/*******************************************************************/
/* Read the section string and set the starting and ending pixels
   based on that. */
void
onecrop_parse_section(struct cropparams *p, size_t *dsize,
                      long *fpixel, long *lpixel)
{
  int add;
  long read;
  char *tailptr;
  long naxes[MAXDIM];
  char forl='f', *pt=p->section;
  size_t i, ndim=p->imgs->ndim, dim=0;


  /* When the user asks for a section of the dataset, then the cropped
     region is not defined by its center. So it makes no sense to later
     check if the center is blank or not. Hence, we will over-write it with
     zero. */
  p->checkcenter=0;


  /* Initialize the fpixel and lpixel arrays (note that 'section' is only
     defined in image mode, so there will only be one element in 'imgs'. */
  for(i=0;i<ndim;++i)
    {
      fpixel[i] = 1;
      lpixel[i] = naxes[i] = p->imgs->dsize[ ndim - i - 1 ];
    }


  /* Parse the string: 'forl': "first-or-last". */
  while(*pt!='\0')
    {
      add=0;
      switch(*pt)
        {
        case ',':
          ++dim;
          if(dim>=ndim)
            error(EXIT_FAILURE, 0, "Extra ',' in '%s'", p->section);
          forl='f';
          ++pt;
          break;
        case ':':
          forl='l';
          ++pt;
          break;
        case '.':
          error(EXIT_FAILURE, 0, "the numbers in the argument to "
                "'--section' ('-s') have to be integers. You input "
                "includes a float number: %s", p->section);
          break;
        case ' ': case '\t':
          ++pt;
          break;
        case '*':
          add=1;                /* If it is an asterisk, then add the */
          ++pt;                 /* given value to the maximum size of */
          break;                /* the image. */

        /* Numerical characters signify the start of a number, so we don't
           need to increment the pointer and can just break out. */
        case '0': case '1': case '2': case '3': case '4': case '5':
        case '6': case '7': case '8': case '9': case '-':
          break;

        /* An un-recognized character should crash the program. */
        default:
          error(EXIT_FAILURE, 0, "value to '--section' must only contain "
                "integer numbers and these special characters between them: "
                "',', ':', '*' when necessary. But it is '%s' (the first "
                "non-acceptable character is '%c').\n\n"
                "Please run the command below to learn more about this "
                "option in Gnuastro's Crop program:\n\n"
                "    $ info gnuastro \"Crop section syntax\"\n", p->section,
                *pt);
          break;
        }

      /* Read the number: */
      read=strtol(pt, &tailptr, 0);

      /* Check if there actually was a number.
      printf("\n\n------\n%c: %ld (%s)\n", *pt, read, tailptr);
      */

      /* Make sure if a number was read at all? */
      if(tailptr==pt)           /* No number was read!                 */
        {
          if(add) read=0;       /* We have a * followed by ':' or ','. */
          else    continue;
        }

      /* Put it in the correct array. */
      if(forl=='f')
        fpixel[dim] = add ? naxes[dim]+read : read;
      else
        lpixel[dim] = add ? naxes[dim]+read : read;
      pt=tailptr;
    }


  /* Make sure the first pixel is located before/below the last pixel. */
  for(i=0;i<ndim;++i)
    if(fpixel[i]>lpixel[i])
      error(EXIT_FAILURE, 0, "the bottom left corner coordinates "
            "cannot be larger or equal to the top right's! Your section "
            "string (%s) has been read as: bottom left coordinate "
            "(%ld, %ld) to top right coordinate (%ld, %ld)",
            p->section, fpixel[0], fpixel[1], lpixel[0], lpixel[1]);

  /* For a check:
  printf("\n%s\n", p->section);
  printf("fpixel: ("); for(i=0;i<ndim;++i) printf("%ld, ", fpixel[i]);
  printf("\b\b)\n");
  printf("lpixel: ("); for(i=0;i<ndim;++i) printf("%ld, ", lpixel[i]);
  printf("\b\b)\n\n");
  exit(0);
  */
}





static void
onecrop_ipolygon_fl(double *ipolygon, size_t nvertices, long *fpixel,
                    long *lpixel)
{
  size_t i;
  double minx=FLT_MAX, miny=FLT_MAX;
  double maxx=-FLT_MAX, maxy=-FLT_MAX;

  /* Find their minimum and maximum values. */
  for(i=0;i<nvertices;++i)
    {
      if(ipolygon[i*2]>maxx) maxx=ipolygon[i*2];
      if(ipolygon[i*2]<minx) minx=ipolygon[i*2];
      if(ipolygon[i*2+1]>maxy) maxy=ipolygon[i*2+1];
      if(ipolygon[i*2+1]<miny) miny=ipolygon[i*2+1];
    }

  /* Set the first and last pixel. */
  fpixel[0] = minx - (int)minx >=0.5 ? (int)minx + 1 : (int)minx;
  fpixel[1] = miny - (int)miny >=0.5 ? (int)miny + 1 : (int)miny;
  lpixel[0] = maxx - (int)maxx >=0.5 ? (int)maxx + 1 : (int)maxx;
  lpixel[1] = maxy - (int)maxy >=0.5 ? (int)maxy + 1 : (int)maxy;
}





#define POLYGON_MASK(CTYPE) {                                           \
    CTYPE *ba=array, *bb=gal_blank_alloc_write(type);                   \
    for(i=0;i<size;++i)                                                 \
      {                                                                 \
        point[0]=i%s1+1; point[1]=i/s1+1;                               \
        if((*isinside)(ipolygon, point, nvertices)==polygonout)         \
          ba[i]=*bb;                                                    \
      }                                                                 \
    free(bb);                                                           \
  }


void
polygonmask(struct onecropparams *crp, void *array, long *fpixel_i,
            size_t s0, size_t s1)
{
  int type=crp->p->type;
  double *ipolygon, point[2];
  int polygonout=crp->p->polygonout;
  int (*isinside)(double *, double *, size_t);
  size_t i, *ordinds, size=s0*s1, nvertices=crp->p->nvertices;

  /* First of all, allocate enough space to put a copy of the input
     coordinates (we will be using that after sorting in an
     anti-clickwise manner.) */
  errno=0; ipolygon=malloc(2*nvertices*sizeof *ipolygon);
  if(ipolygon==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for ipolygon",
          __func__, 2*nvertices*sizeof *ipolygon);
  errno=0; ordinds=malloc(nvertices*sizeof *ordinds);
  if(ordinds==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for ordinds",
          __func__, nvertices*sizeof *ordinds);

  /* If the user wants to sort the edges, do it. If not, make sure its in
     counter-clockwise orientation. */
  if(crp->p->polygonsort)
    gal_polygon_vertices_sort(crp->ipolygon, nvertices, ordinds);
  else
    {
      /* Keep the original order of the vertices, just make it
         counter-clockwise. */
      for(i=0;i<nvertices;++i) ordinds[i]=i;
      gal_polygon_to_counterclockwise(crp->ipolygon, nvertices);
    }

  /* Fill the final polygon vertice positions within 'ipolygon' and also
     the fpixel_i coordinates from all the vertices to bring them into the
     crop image coordinates. */
  for(i=0;i<nvertices;++i)
    {
      ipolygon[i*2  ] = crp->ipolygon[ordinds[i]*2]   - fpixel_i[0];
      ipolygon[i*2+1] = crp->ipolygon[ordinds[i]*2+1] - fpixel_i[1];
    }

  /* Print a warning if we done the sorting, _and_ the sorted polygon is
     concave, _and_ the user hasn't activated the quiet mode. Note that we
     could'n do this immediately after calling 'gal_polygon_vertices_sort'
     because that function doesn't touch the actual vertice values, it only
     fills 'ordinds'. */
  if(crp->p->polygonsort
     && !crp->p->cp.quiet
     && !gal_polygon_is_convex(ipolygon, nvertices) )
    error(0, 0, "%s: warning: the given vertices belong to a concave "
          "polygon, but there is no unique solution to the sorting of "
          "concave polygons. Please check the cropped image, if it doesn't "
          "correspond to your desired polygon, sort the vertices yourself "
          "in counter-clockwise order _and_ don't use the '--polygonsort' "
          "option", __func__);

  /* Set the function for checking if a point is inside the polygon. For
     concave polygons the process is more complex and thus
     slower. Therefore when the polygon is convex, its better to use the
     simpler/faster function. */
  isinside = ( gal_polygon_is_convex(ipolygon, nvertices)
               ? gal_polygon_is_inside_convex
               : gal_polygon_is_inside );

  /* Go over all the pixels in the image and if they are within the
     polygon keep them if the user has asked for it.*/
  switch(type)
    {
    case GAL_TYPE_UINT8:    POLYGON_MASK(uint8_t);  break;
    case GAL_TYPE_INT8:     POLYGON_MASK(int8_t);   break;
    case GAL_TYPE_UINT16:   POLYGON_MASK(uint16_t); break;
    case GAL_TYPE_INT16:    POLYGON_MASK(int16_t);  break;
    case GAL_TYPE_UINT32:   POLYGON_MASK(uint32_t); break;
    case GAL_TYPE_INT32:    POLYGON_MASK(int32_t);  break;
    case GAL_TYPE_INT64:    POLYGON_MASK(int64_t);  break;
    case GAL_TYPE_FLOAT32:  POLYGON_MASK(float);    break;
    case GAL_TYPE_FLOAT64:  POLYGON_MASK(double);   break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s, so we "
            "can fix the problem. Type code %d is not recognized",
            __func__, PACKAGE_BUGREPORT, type);
    }

  /* Clean up: */
  free(ordinds);
  free(ipolygon);
}



















/*******************************************************************/
/******************          One crop.         *********************/
/*******************************************************************/
static void
onecrop_zero_to_nan(void *array, size_t size, int type)
{
  float *fp, *ffp;
  double *dp, *fdp;

  switch(type)
    {
    case GAL_TYPE_FLOAT32:
      ffp=(fp=array)+size;
      do
        if(*fp==0.0f) *fp=NAN;
      while(++fp<ffp);
      break;

    case GAL_TYPE_FLOAT64:
      fdp=(dp=array)+size;
      do
        if(*dp==0.0f) *dp=NAN;
      while(++dp<fdp);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: %d is not a recognized type",
            __func__, type);
    }
}






/* Set the output name and image output sizes. */
void
onecrop_name(struct onecropparams *crp)
{
  char **strarr;
  struct cropparams *p=crp->p;
  struct gal_options_common_params *cp=&p->cp;

  /* Set the output name and crop sides: */
  if(p->catname)
    {
      /* If a name column was set, use it, otherwise, use the ID of the
         profile. */
      if(p->name)
        {
          strarr=p->name;
          if( asprintf(&crp->name, "%s%s%s", cp->output, strarr[crp->out_ind],
                       p->suffix)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        if( asprintf(&crp->name, "%s%zu%s", cp->output, crp->out_ind+1,
                     p->suffix)<0 )
          error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);

      /* Make sure the file doesn't exist. */
      gal_checkset_writable_remove(crp->name, 0, cp->dontdelete);
    }
  else
    {
      /* Set the output name. */
      if(p->outnameisfile)            /* An output file was specified. */
        {
          crp->name=cp->output;
          gal_checkset_writable_remove(crp->name, 0, cp->dontdelete);
        }
      else          /* The output was a directory, use automatic output. */
        crp->name=gal_checkset_automatic_output(cp,
                                                p->imgs[crp->in_ind].name,
                                                p->suffix);
    }
}





/* Find the first and last pixel of a crop. */
static void
onecrop_flpixel(struct onecropparams *crp)
{
  struct cropparams *p=crp->p;
  size_t ndim=p->imgs->ndim;

  double center[MAXDIM];
  int ncoord=1, status=0;
  size_t i, *dsize=p->imgs[crp->in_ind].dsize;
  long *fpixel=crp->fpixel, *lpixel=crp->lpixel;
  double pixcrd[MAXDIM], imgcrd[MAXDIM], phi[1], theta[1];

  switch(p->mode)
    {

    case IMGCROP_MODE_IMG:
      if(p->section)            /* Defined by section.    */
        onecrop_parse_section(p, dsize, fpixel, lpixel);
      else if(p->polygon)       /* Defined by a polygon.  */
        {
          if(p->polygonout==0)
            onecrop_ipolygon_fl(p->ipolygon, p->nvertices, fpixel, lpixel);
        }
      else                      /* Defined by its center. */
        {
          for(i=0;i<ndim;++i) center[i] = p->centercoords[i][crp->out_ind];
          gal_box_border_from_center(center, ndim, p->iwidth, fpixel, lpixel);
        }
      break;


    case IMGCROP_MODE_WCS: /* In wcsmode, crp->world is already filled.   */
      if(p->polygon)       /* Note: p->iwidth was set based on p->wwidth. */
        {
          /* Fill crp->ipolygon in wcspolygonpixel, then set flpixel*/
          fillcrpipolygon(crp);
          if(p->polygonout==0)
            onecrop_ipolygon_fl(crp->ipolygon, p->nvertices, fpixel, lpixel);
        }
      else
        {
          /* Convert 'crp->world' (in WCS) into 'pixcrd' (image coord). */
          if(wcss2p(p->imgs[crp->in_ind].wcs, ncoord, ndim, crp->world,
                    phi, theta, imgcrd, pixcrd, &status) )
            if(status)
              error(EXIT_FAILURE, 0, "%s: wcss2p error %d: %s", __func__,
                    status, wcs_errmsg[status]);

          /* Find the first and last pixels of this crop. */
          gal_box_border_from_center(pixcrd, ndim, p->iwidth, fpixel, lpixel);
        }
      break;


    default:
      error(EXIT_FAILURE, 0, "%s: a bug! The domain (WCS or image) are not "
            "set. Please contact us at %s so we can see how it got to this "
            "impossible place", __func__, PACKAGE_BUGREPORT);
    }


  /* If the user only wants regions outside to the polygon, then set
     the fpixel and lpixel to cover the full input image. */
  if(p->polygon && p->polygonout)
    {
      crp->fpixel[0]=crp->fpixel[1]=1;
      crp->lpixel[0]=dsize[1];
      crp->lpixel[1]=dsize[0];
    }
}






/* Find the size of the final FITS image (irrespective of how many
   crops will be needed for it) and make the image to keep the
   data.

   NOTE: The fpixel and lpixel in crp keep the first and last pixel of
   the total image for this crop, irrespective of the final keeping
   blank areas or not. While the fpixel_i and lpixel_i arrays keep the
   first and last pixels after the blank pixels have been removed.
*/
static void
onecrop_make_array(struct onecropparams *crp, long *fpixel_i,
                   long *lpixel_i, long *fpixel_c, long *lpixel_c)
{
  double crpix;
  fitsfile *ofp;
  long naxes[MAXDIM];
  char *outname=crp->name;
  int status=0, type=crp->p->type;
  char **strarr, cpname[FLEN_KEYWORD];
  gal_data_t *rkey=gal_data_array_calloc(1);
  size_t i, ndim=crp->p->imgs->ndim, totsize;
  char *cp, *cpf, blankrec[80], titlerec[80];
  struct inputimgs *img=&crp->p->imgs[crp->in_ind];


  /* Set the last element of the blank array. */
  cpf=blankrec+79;
  *cpf='\0';
  titlerec[79]='\0';
  cp=blankrec; do *cp=' '; while(++cp<cpf);


  /* Set the size of the output, in WCS mode, noblank==0. */
  if(crp->p->noblank && crp->p->mode==IMGCROP_MODE_IMG)
    for(i=0;i<ndim;++i)
      {
        fpixel_c[i] = 1;
        lpixel_c[i] = naxes[i] = lpixel_i[i]-fpixel_i[i]+1;
      }
  else
    for(i=0;i<ndim;++i)
      naxes[i] = crp->lpixel[i]-crp->fpixel[i]+1;


  /* Create the FITS file with a blank first extension, then close it, so
     with the next 'fits_open_file', we build the image in the second
     extension. This way, atleast for Gnuastro's outputs, we can
     consistently use '-h1' (something like how you count columns, or
     generally everything from 1). */
  if(fits_create_file(&ofp, outname, &status))
    gal_fits_io_error(status, "creating file");
  fits_create_img(ofp, SHORT_IMG, 0, naxes, &status);
  fits_close_file(ofp, &status);


  /* Create the output crop image. */
  fits_open_file(&crp->outfits, outname, READWRITE, &status);
  fits_create_img(crp->outfits, gal_fits_type_to_bitpix(type),
                  ndim, naxes, &status);
  gal_fits_io_error(status, "creating image");
  ofp=crp->outfits;


  /* When CFITSIO creates a FITS extension it adds two comments linking to
     the FITS paper. Since we are mentioning the version of CFITSIO and
     only use its routines to read/write from/to FITS files, this is
     redundant. If 'status!=0', then 'gal_fits_io_error' will abort, but in
     case CFITSIO doesn't write the comments, status will become
     non-zero. So we are resetting it to zero after these (because not
     being able to delete them isn't an error). */
  fits_delete_key(ofp, "COMMENT", &status);
  fits_delete_key(ofp, "COMMENT", &status);
  status=0;


  /* Read the units of the input dataset and store them in the output. */
  rkey->next=NULL;
  rkey->name="BUNIT";
  rkey->type=GAL_TYPE_STRING;
  gal_fits_key_read_from_ptr(crp->infits, rkey, 1, 1);
  if(rkey->status==0)           /* The BUNIT keyword was read. */
    {
      strarr=rkey->array;
      fits_update_key(ofp, TSTRING, "BUNIT", strarr[0], "physical units",
                      &status);
      gal_fits_io_error(status, "writing BUNIT");
    }
  rkey->name=NULL;              /* 'name' wasn't allocated. */
  gal_data_free(rkey);


  /* Write the blank value as a FITS keyword if necessary. */
  if( type!=GAL_TYPE_FLOAT32 && type!=GAL_TYPE_FLOAT64 )
    if(fits_write_key(ofp, gal_fits_type_to_datatype(crp->p->type), "BLANK",
                      crp->p->blankptrwrite, "Pixels with no data.",
                      &status) )
      gal_fits_io_error(status, "adding Blank");
  totsize = naxes[0]*naxes[1] * (ndim==3?naxes[2]:1);
  if(fits_write_null_img(ofp, 1, totsize, &status))
    gal_fits_io_error(status, "writing null array");


  /* Write the WCS header keywords in the output FITS image, then
     update the header keywords. */
  if(img->wcs)
    {
      /* Write the WCS title and common WCS information. */
      if(fits_write_record(ofp, blankrec, &status))
        gal_fits_io_error(status, NULL);
      sprintf(titlerec, "%sWCS information", GAL_FITS_KEY_TITLE_START);
      for(i=strlen(titlerec);i<79;++i)
        titlerec[i]=' ';
      fits_write_record(ofp, titlerec, &status);
      for(i=0;i<img->nwcskeys-1;++i)
        fits_write_record(ofp, &img->wcstxt[i*80], &status);
      gal_fits_io_error(status, NULL);

      /* Correct the CRPIX keywords. */
      for(i=0;i<ndim;++i)
        {
          sprintf(cpname, "CRPIX%zu", i+1);
          crpix = img->wcs->crpix[i] - (fpixel_i[i]-1) + (fpixel_c[i]-1);
          fits_update_key(ofp, TDOUBLE, cpname, &crpix, NULL, &status);
          gal_fits_io_error(status, NULL);
        }
    }


  /* Add the Crop information. */
  sprintf(titlerec, "%sCrop information", GAL_FITS_KEY_TITLE_START);
  for(i=strlen(titlerec);i<79;++i)
    titlerec[i]=' ';
  if(fits_write_record(ofp, titlerec, &status))
    gal_fits_io_error(status, NULL);
}





/* The starting and ending points are set in the onecropparams structure
   for one crop from one image. Crop that region out of the input.

   On'basekeyname': To be safe, GCC 8.1 (and persumably later versions)
   assumes that we are writing the full statically allocated space into
   'regioinkey'! So it prints a warning that you may be writing outside the
   allocated space! With these variables, we are ultimately just writing
   the file counters, so we can never (with current techologies!!!) exceed
   'FLEN_KEYWORD' (which is 75 characters). To avoid compiler warnings, we
   are just removing a few characters ('FLEN_KEYWORD-5') to allow the
   suffix and remove the warnings. */
void
onecrop(struct onecropparams *crp)
{
  struct cropparams *p=crp->p;
  struct inputimgs *img=&p->imgs[crp->in_ind];

  void *array;
  int status=0, anynul=0;
  fitsfile *ifp=crp->infits, *ofp;
  char basekeyname[FLEN_KEYWORD-5];     /* '-5': avoid gcc 8.1+ warnings! */
  gal_fits_list_key_t *headers=NULL;    /* See above comment for more.    */
  size_t i, j, cropsize=1, ndim=img->ndim;
  char region[FLEN_VALUE], regionkey[FLEN_KEYWORD];
  long fpixel_o[MAXDIM], lpixel_o[MAXDIM], inc[MAXDIM];
  long naxes[MAXDIM], fpixel_i[MAXDIM], lpixel_i[MAXDIM];

  /* Fill the 'naxes' and 'inc' arrays. */
  for(i=0;i<ndim;++i)
    {
      inc[ i ]   = 1;
      naxes[ i ] = img->dsize[ ndim - i - 1 ];
    }


  /* Find the first and last pixel of this crop box from this input
     image. Then copy the first and last pixels into the '_i' arrays.*/
  onecrop_flpixel(crp);
  memcpy(fpixel_i, crp->fpixel, ndim*sizeof *fpixel_i);
  memcpy(lpixel_i, crp->lpixel, ndim*sizeof *lpixel_i);


  /* Find the overlap and apply it if there is any overlap. */
  if( gal_box_overlap(naxes, fpixel_i, lpixel_i, fpixel_o, lpixel_o, ndim) )
    {
      /* Make the output FITS image and initialize it with an array of
         NaN or BLANK values. */
      if(crp->outfits==NULL)
        onecrop_make_array(crp, fpixel_i, lpixel_i, fpixel_o, lpixel_o);
      ofp=crp->outfits;


      /* Allocate an array to keep the desired crop region, then read
         the desired pixels into it. */
      status=0;
      for(i=0;i<ndim;++i) cropsize *= ( lpixel_i[i] - fpixel_i[i] + 1 );
      array=gal_pointer_allocate(p->type, cropsize, 0, __func__, "array");
      if(fits_read_subset(ifp, gal_fits_type_to_datatype(p->type),
                          fpixel_i, lpixel_i, inc, p->blankptrread, array,
                          &anynul, &status))
        gal_fits_io_error(status, NULL);


      /* If we have a floating point or double image, pixels with zero
         value should actually be a NaN. Unless the user specificly
         asks for it, make the conversion.*/
      if(p->zeroisnotblank==0
         && (p->type==GAL_TYPE_FLOAT32
             || p->type==GAL_TYPE_FLOAT64) )
        onecrop_zero_to_nan(array, cropsize, p->type);


      /* If a polygon is given, remove all the pixels within or
         outside of it.*/
      if(p->polygon)
        {
          /* A small sanity check until this part supports 3D. */
          if(ndim!=2)
            error(EXIT_FAILURE, 0, "%s: polygons only implemented in 2D",
                  __func__);

          /* In WCS mode, crp->ipolygon was allocated and filled in
             wcspolygonflpixel (wcsmode.c). */
          if(p->mode==IMGCROP_MODE_IMG) crp->ipolygon=p->ipolygon;
          polygonmask(crp, array, fpixel_i, lpixel_i[1]-fpixel_i[1]+1,
                      lpixel_i[0]-fpixel_i[0]+1);
          if(p->mode==IMGCROP_MODE_WCS) free(crp->ipolygon);
        }


      /* Write the array into the image. */
      status=0;
      if( fits_write_subset(ofp, gal_fits_type_to_datatype(p->type),
                            fpixel_o, lpixel_o, array, &status) )
        gal_fits_io_error(status, NULL);


      /* Write the selected region of this image as a string to include as
         a FITS keyword. Then we want to delete the last coma ','.*/
      j=0;
      for(i=0;i<ndim;++i)
        j += sprintf(&region[j], "%ld:%ld,", fpixel_i[i], lpixel_i[i]);
      region[j-1]='\0';


      /* A section has been added to the cropped image from this input
         image, so save the information of this image. */
      sprintf(basekeyname, "ICF%zu", crp->numimg);
      gal_fits_key_write_filename(basekeyname, img->name, &headers, 0);
      sprintf(regionkey, "%sPIX", basekeyname);
      gal_fits_key_list_add_end(&headers, GAL_TYPE_STRING, regionkey,
                                0, region, 0, "Range of pixels used for "
                                "this output.", 0, NULL);
      gal_fits_key_write_in_ptr(&headers, ofp);


      /* Free the allocated array. */
      free(array);
    }
  else
    if(p->polygon && p->polygonout==0 && p->mode==IMGCROP_MODE_WCS)
      free(crp->ipolygon);

  /* The crop is complete. */
  return;
}




















/*******************************************************************/
/******************        Check center        *********************/
/*******************************************************************/
int
onecrop_center_filled(struct onecropparams *crp)
{
  struct cropparams *p=crp->p;

  void *array;
  size_t size, ndim, *dsize;
  fitsfile *ofp=crp->outfits;
  int status=0, anynul=0, type;
  long checkcenter=p->checkcenter;
  long naxes[3], fpixel[3], lpixel[3], inc[3]={1,1,1};

  /* If checkcenter is zero, then don't check. */
  if(checkcenter==0) return GAL_BLANK_UINT8;

  /* Get the final size of the output image. */
  gal_fits_img_info(ofp, &type, &ndim, &dsize, NULL, NULL);
  if(ndim==2)
    {
      naxes[0]=dsize[1];
      naxes[1]=dsize[0];
    }
  else
    {
      naxes[0]=dsize[2];
      naxes[1]=dsize[1];
      naxes[2]=dsize[0];
    }

  /* Get the size and range of the central region to check. The +1 is
     because in FITS, counting begins from 1, not zero. It might happen
     that the image is actually smaller than the width to check the center
     (for example 1 or 2 pixels wide). In that case, we'll just use the
     full image to check. */
  size = ( (naxes[0]>checkcenter ? checkcenter : naxes[0])
           * (naxes[1]>checkcenter ? checkcenter : naxes[1]) );
  fpixel[0] = naxes[0]>checkcenter ? ((naxes[0]/2+1)-checkcenter/2) : 1;
  fpixel[1] = naxes[1]>checkcenter ? ((naxes[1]/2+1)-checkcenter/2) : 1;
  lpixel[0] = ( naxes[0]>checkcenter
                ? ((naxes[0]/2+1)+checkcenter/2) : naxes[0] );
  lpixel[1] = ( naxes[1]>checkcenter
                ? ((naxes[1]/2+1)+checkcenter/2) : naxes[1] );


  /* For the third dimension. */
  if(ndim==3)
    {
      size *= (naxes[2]>checkcenter ? checkcenter : naxes[2]);
      fpixel[2] = naxes[2]>checkcenter ? ((naxes[2]/2+1)-checkcenter/2) : 1;
      lpixel[2] = ( naxes[2]>checkcenter
                    ? ((naxes[2]/2+1)+checkcenter/2) : naxes[2] );
    }


  /* For a check:
  printf("naxes: %ld, %ld\nfpixel: (%ld, %ld)\nlpixel: (%ld, %ld)\n"
         "size: %zu\n", naxes[0], naxes[1], fpixel[0], fpixel[1],
         lpixel[0], lpixel[1], size);
  */

  /* Allocate the array and read in the pixels. */
  array=gal_pointer_allocate(type, size, 0, __func__, "array");
  if( fits_read_subset(ofp, gal_fits_type_to_datatype(type), fpixel, lpixel,
                       inc, p->blankptrread, array, &anynul, &status) )
    gal_fits_io_error(status, NULL);
  free(array);

  /* CFITSIO already checks if there are any blank pixels. If there are,
     then 'anynul' will be 1, if there aren't it will be 0. So the output
     of this function is just the inverse of that number. */
  return !anynul;
}
