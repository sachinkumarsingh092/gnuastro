/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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
sectionparser(struct cropparams *p, size_t *dsize,
              long *fpixel, long *lpixel)
{
  int add;
  long read;
  size_t dim=0;
  char *tailptr;
  char forl='f', *pt=p->section;
  long naxes[2]={dsize[1], dsize[0]};

  /* If control reached here, then the cropped region is not defined by its
     center. So it makes no sense to check if the center is blank. */
  p->checkcenter=0;

  /* Initialize the fpixel and lpixel arrays: */
  lpixel[0]=naxes[0];
  lpixel[1]=naxes[1];
  fpixel[0]=fpixel[1]=1;

  /* Parse the string. */
  while(*pt!='\0')
    {
      add=0;
      switch(*pt)
        {
        case ',':
          ++dim;
          if(dim==2)
            error(EXIT_FAILURE, 0, "Extra `,` in `%s`", p->section);
          forl='f';
          ++pt;
          break;
        case ':':
          forl='l';
          ++pt;
          break;
        case '.':
          error(EXIT_FAILURE, 0, "the numbers in the argument to "
                "`--section` (`-s') have to be integers. You input "
                "includes a float number: %s", p->section);
          break;
        case ' ': case '\t':
          ++pt;
          break;
        case '*':
          add=1;                /* If it is an asterisk, then add the */
          ++pt;                 /* given value to the maximum size of */
          break;                /* the image. */
        default:
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
          if(add) read=0;       /* We have a * followed by `:' or `,'. */
          else    continue;
        }

      /* Put it in the correct array. Note that for the last point, we
         don't want to include the given pixel. Unlike CFITSIO, in Crop,
         the given intervals are not inclusive. fpixel and lpixel will be
         directly passed to CFITSIO. So we have to correct his here.*/
      if(forl=='f')
        fpixel[dim] = add ? naxes[dim]+read : read;
      else
        lpixel[dim] = add ? naxes[dim]+read : read;
      pt=tailptr;

      /* For a check:
      printf("\n\n[%ld, %ld]: fpixel=(%ld, %ld), lpixel=(%ld, %ld)\n\n",
             naxes[0], naxes[1],
             fpixel[0], fpixel[1], lpixel[0], lpixel[1]);
      */
    }

  /* Make sure the first pixel is located before/below the last pixel. */
  if(fpixel[0]>lpixel[0] || fpixel[1]>lpixel[1])
    error(EXIT_FAILURE, 0, "the bottom left corner coordinates "
          "cannot be larger or equal to the top right's! Your section "
          "string (%s) has been read as: bottom left coordinate "
          "(%ld, %ld) to top right coordinate (%ld, %ld)",
          p->section, fpixel[0], fpixel[1], lpixel[0], lpixel[1]);

  /*
  printf("\n%s\nfpixel=(%ld, %ld), lpixel=(%ld, %ld)\n\n", section,
         fpixel[0], fpixel[1], lpixel[0], lpixel[1]);
  exit(0);
  */
}






void
crop_polygonparser(struct cropparams *p)
{
  size_t dim=0;
  char *tailptr;
  char *pt=p->polygon;
  double read, *array;
  gal_list_f64_t *vertices=NULL;

  /* If control reached here, then the cropped region is not defined by its
     center. So it makes no sense to check if the center is blank. */
  p->checkcenter=0;

  /* Parse the string. */
  while(*pt!='\0')
    {
      switch(*pt)
        {
        case ',':
          ++dim;
          if(dim==2)
            error(EXIT_FAILURE, 0, "Extra `,` in `%s`", p->polygon);
          ++pt;
          break;
        case ':':
          if(dim==0)
            error(EXIT_FAILURE, 0, "not enough coordinates for at least "
                  "one polygon vertex (in %s)", p->polygon);
          dim=0;
          ++pt;
          break;
        default:
          break;
        }

      /* strtod will skip white spaces if they are before a number,
         but not when they are before a : or ,. So we need to remove
         all white spaces. White spaces are usually put beside each
         other, so if one is encountered, go along the string until
         the white space characters finish.  */
      if(isspace(*pt))
        ++pt;
      else
        {
          /* Read the number: */
          read=strtod(pt, &tailptr);

          /* Check if there actually was a number.
          printf("\n\n------\n%zu: %f (%s)\n", dim, read, tailptr);
          */

          /* Make sure if a number was read at all? */
          if(tailptr==pt) /* No number was read! */
            error(EXIT_FAILURE, 0, "%s could not be parsed as a floating "
                  "point number", tailptr);

          /* Check if there are no extra characters in the number, for
             example we don't have a case like `1.00132.17', or
             1.01i:2.0. Such errors are not uncommon when typing large
             numbers, and if ignored, they can lead to unpredictable
             results, so its best to abort and inform the user. */
          if( *tailptr!='\0'
              && !isspace(*tailptr)
              && strchr(":,", *tailptr)==NULL )
            error(EXIT_FAILURE, 0, "'%s' is an invalid floating point number "
                  "sequence in the value to the `--polygon' option, error "
                  "detected at '%s'", pt, tailptr);

          /* Add the read coordinate to the list of coordinates. */
          gal_list_f64_add(&vertices, read);

          /* The job here is done, start from tailptr */
          pt=tailptr;
        }
    }

  /* Put the coordinates into an array while reversing their order so they
     correspond to the user's order, then put it in the right place.*/
  array=gal_list_f64_to_array(vertices, 1, &p->nvertices);
  if(p->mode==IMGCROP_MODE_IMG) { p->ipolygon=array; p->wpolygon=NULL;  }
  else                          { p->ipolygon=NULL;  p->wpolygon=array; }

  /* The number of vertices is actually the number of nodes in the list
     divided by the dimension of the dataset (note that we were counting
     the dimension from 0. */
  p->nvertices/=(dim+1);

  /* For a check:
  {
    size_t i;
    double *polygon=p->mode==IMGCROP_MODE_IMG?p->ipolygon:p->wpolygon;
    for(i=0;i<p->nvertices;++i)
      printf("(%f, %f)\n", polygon[i*2], polygon[i*2+1]);
  }
  exit(0);
  */

  /* Clean up: */
  gal_list_f64_free(vertices);
}





void
imgpolygonflpixel(double *ipolygon, size_t nvertices, long *fpixel,
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
        if(gal_polygon_pin(ipolygon, point, nvertices)==outpolygon)     \
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
  int outpolygon=crp->p->outpolygon;
  size_t i, *ordinds, size=s0*s1, nvertices=crp->p->nvertices;


  /* First of all, allocate enough space to put a copy of the input
     coordinates (we will be using that after sorting in an
     anti-clickwise manner.) */
  errno=0; ipolygon=malloc(2*nvertices*sizeof *ipolygon);
  if(ipolygon==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for ipolygon in polygonmask "
          "(crop.c)", 2*nvertices*sizeof *ipolygon);
  errno=0; ordinds=malloc(nvertices*sizeof *ordinds);
  if(ordinds==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for ordinds in polygonmask "
          "(crop.c)", nvertices*sizeof *ordinds);


  /* Find the order of the polygons and put the elements in the proper
     order. Also subtract the fpixel_i coordinates from all the
     vertices to bring them into the crop image coordinates.*/
  gal_polygon_ordered_corners(crp->ipolygon, crp->p->nvertices, ordinds);
  for(i=0;i<crp->p->nvertices;++i)
    {
      ipolygon[i*2  ] = crp->ipolygon[ordinds[i]*2]   - fpixel_i[0];
      ipolygon[i*2+1] = crp->ipolygon[ordinds[i]*2+1] - fpixel_i[1];
    }


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
      error(EXIT_FAILURE, 0, "a bug! Please contact us at %s, so we "
            "can fix the problem. For some reason, an unrecognized "
            "type value (%d) has been seen in polygonmask (crop.c)",
            PACKAGE_BUGREPORT, type);
    }

  /* Clean up: */
  free(ordinds);
  free(ipolygon);
}



















/*******************************************************************/
/******************          One crop.         *********************/
/*******************************************************************/
void
changezerotonan(void *array, size_t size, int type)
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
      error(EXIT_FAILURE, 0, "%d is not a recognized type in "
            "`changezerotonan'", type);
    }
}






/* Set the output name and image output sizes. */
void
cropname(struct onecropparams *crp)
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
          asprintf(&crp->name, "%s%s%s", cp->output, strarr[crp->out_ind],
                   p->suffix);
        }
      else
        asprintf(&crp->name, "%s%zu%s", cp->output, crp->out_ind+1,
                 p->suffix);

      /* Make sure the file doesn't exist. */
      gal_checkset_check_remove_file(crp->name, 0, cp->dontdelete);
    }
  else
    {
      /* Set the output name. */
      if(p->outnameisfile)            /* An output file was specified. */
        {
          crp->name=cp->output;
          gal_checkset_check_remove_file(crp->name, 0, cp->dontdelete);
        }
      else          /* The output was a directory, use automatic output. */
        crp->name=gal_checkset_automatic_output(cp,
                                                p->imgs[crp->in_ind].name,
                                                p->suffix);
    }
}





/* Find the first and last pixel of a crop from its center point (in
   image mode or WCS mode). */
void
cropflpixel(struct onecropparams *crp)
{
  struct cropparams *p=crp->p;
  int ncoord=1, nelem=2, status[2]={0,0};
  size_t *dsize=p->imgs[crp->in_ind].dsize;
  double pixcrd[2], imgcrd[2], phi[1], theta[1];
  long *fpixel=crp->fpixel, *lpixel=crp->lpixel;

  switch(p->mode)
    {

    case IMGCROP_MODE_IMG:
      if(p->catname)
        gal_box_border_from_center(p->c1[crp->out_ind], p->c2[crp->out_ind],
                                   p->iwidth, fpixel, lpixel);
      else if(!isnan(p->xc))
        gal_box_border_from_center(p->xc, p->yc, p->iwidth, fpixel, lpixel);
      else if(p->section)
        sectionparser(p, dsize, fpixel, lpixel);
      else if(p->polygon)
        {
          if(p->outpolygon==0)
            imgpolygonflpixel(p->ipolygon, p->nvertices, fpixel, lpixel);
        }
      else
        error(EXIT_FAILURE, 0, "a bug! In image mode, neither of the "
              "following has been set: a catalog, a central pixel, "
              "a section or a polygon in the image. Please contact us "
              "to see how it got to this impossible place! You should "
              "have been warned of this condition long before Crop "
              "reaches this point");
      break;

    case IMGCROP_MODE_WCS: /* In wcsmode, crp->world is already filled.   */
      if(p->polygon)       /* Note: p->iwidth was set based on p->wwidth. */
        {
          /* Fill crp->ipolygon in wcspolygonpixel, then set flpixel*/
          fillcrpipolygon(crp);
          if(p->outpolygon==0)
            imgpolygonflpixel(crp->ipolygon, p->nvertices, fpixel, lpixel);
        }
      else
        {
          if(wcss2p(p->imgs[crp->in_ind].wcs, ncoord, nelem, crp->world,
                    phi, theta, imgcrd, pixcrd, status) )
            if(status[0] || status[1])
              error(EXIT_FAILURE, 0, "wcss2p error %d: %s",
                    status[0] ? status[0] : status[1],
                    wcs_errmsg[status[0] ? status[0] : status[1]]);
          gal_box_border_from_center(pixcrd[0], pixcrd[1], p->iwidth, fpixel,
                                     lpixel);
          /*
            printf("\n(%f, %f): (%ld, %ld) -- (%ld, %ld)\n\n", pixcrd[0],
                   pixcrd[1], fpixel[0], fpixel[1], lpixel[0], lpixel[1]);
          */
        }
      break;

    default:
      error(EXIT_FAILURE, 0, "a bug! in cropflpixel (crop.c), "
            "neither imgmode or wcsmode are set. Please contact us so "
            "we can see how it got to this impossible place");
    }

  /* If the user only wants regions outside to the polygon, then set
     the fpixel and lpixel to cover the full input image. */
  if(p->polygon && p->outpolygon)
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
void
firstcropmakearray(struct onecropparams *crp, long *fpixel_i,
                   long *lpixel_i, long *fpixel_c, long *lpixel_c)
{
  size_t i;
  fitsfile *ofp;
  long naxes[2];
  int type=crp->p->type;
  double crpix0, crpix1;
  int naxis=2, status=0;
  char *outname=crp->name;
  char *cp, *cpf, blankrec[80], titlerec[80];
  char startblank[]="                      / ";
  struct inputimgs *img=&crp->p->imgs[crp->in_ind];


  /* Set the last element of the blank array. */
  cpf=blankrec+79;
  *cpf='\0';
  titlerec[79]='\0';
  cp=blankrec; do *cp=' '; while(++cp<cpf);


  /* Set the size of the output, in WCS mode, noblank==0. */
  if(crp->p->noblank && crp->p->mode==IMGCROP_MODE_IMG)
    {
      fpixel_c[0]=fpixel_c[1]=1;
      lpixel_c[0]=naxes[0]=lpixel_i[0]-fpixel_i[0]+1;
      lpixel_c[1]=naxes[1]=lpixel_i[1]-fpixel_i[1]+1;
    }
  else
    {
      naxes[0]=crp->lpixel[0]-crp->fpixel[0]+1;
      naxes[1]=crp->lpixel[1]-crp->fpixel[1]+1;
    }


  /* Create the FITS file with a blank first extension, then close it, so
     with the next `fits_open_file', we build the image in the second
     extension. This way, atleast for Gnuastro's outputs, we can
     consistently use `-h1' (something like how you count columns, or
     generally everything from 1). */
  if(fits_create_file(&ofp, outname, &status))
    gal_fits_io_error(status, "creating file");
  fits_create_img(ofp, SHORT_IMG, 0, naxes, &status);
  fits_close_file(ofp, &status);

  /* Create the output crop image. */
  fits_open_file(&crp->outfits, outname, READWRITE, &status);
  fits_create_img(crp->outfits, gal_fits_type_to_bitpix(type),
                  naxis, naxes, &status);
  gal_fits_io_error(status, "creating image");
  ofp=crp->outfits;

  /* When CFITSIO creates a FITS extension it adds two comments linking to
     the FITS paper. Since we are mentioning the version of CFITSIO and
     only use its ruitines to read/write from/to FITS files, this is
     redundant. If `status!=0', then `gal_fits_io_error' will abort, but in
     case CFITSIO doesn't write the comments, status will become
     non-zero. So we are resetting it to zero after these (because not
     being able to delete them isn't an error). */
  fits_delete_key(ofp, "COMMENT", &status);
  fits_delete_key(ofp, "COMMENT", &status);
  status=0;

  /* Write the blank value if necessary. */
  if( type!=GAL_TYPE_FLOAT32 && type!=GAL_TYPE_FLOAT64 )
    if(fits_write_key(ofp, gal_fits_type_to_datatype(crp->p->type), "BLANK",
                      crp->p->bitnul, "pixels with no data", &status) )
      gal_fits_io_error(status, "adding Blank");
  if(fits_write_null_img(ofp, 1, naxes[0]*naxes[1], &status))
    gal_fits_io_error(status, "writing null array");


  /* Write the WCS header keywords in the output FITS image, then
     update the header keywords. */
  if(img->wcs)
    {
      crpix0 = img->wcs->crpix[0] - (fpixel_i[0]-1) + (fpixel_c[0]-1);
      crpix1 = img->wcs->crpix[1] - (fpixel_i[1]-1) + (fpixel_c[1]-1);
      if(fits_write_record(ofp, blankrec, &status))
        gal_fits_io_error(status, NULL);
      sprintf(titlerec, "%sWCS information", startblank);
      for(i=strlen(titlerec);i<79;++i)
        titlerec[i]=' ';
      fits_write_record(ofp, titlerec, &status);
      for(i=0;i<img->nwcskeys-1;++i)
        fits_write_record(ofp, &img->wcstxt[i*80], &status);
      fits_update_key(ofp, TDOUBLE, "CRPIX1", &crpix0, NULL, &status);
      fits_update_key(ofp, TDOUBLE, "CRPIX2", &crpix1, NULL, &status);
      gal_fits_io_error(status, NULL);
    }


  /* Add the Crop information. */
  if(fits_write_record(ofp, blankrec, &status))
    gal_fits_io_error(status, NULL);
  sprintf(titlerec, "%sCrop information", startblank);
  for(i=strlen(titlerec);i<79;++i)
    titlerec[i]=' ';
  if(fits_write_record(ofp, titlerec, &status))
    gal_fits_io_error(status, NULL);
}





/* The starting and ending points are set in the onecropparams structure
   for one crop from one image. Crop that region out.

   return values are:
   0: No crop was made (not in range of image).
   1: The input image covered at least part of the crop image.
 */
void
onecrop(struct onecropparams *crp)
{
  struct cropparams *p=crp->p;
  struct inputimgs *img=&p->imgs[crp->in_ind];

  void *array;
  size_t cropsize;
  int status=0, anynul=0;
  char basename[FLEN_KEYWORD];
  fitsfile *ifp=crp->infits, *ofp;
  struct gal_fits_key_ll *headers=NULL;
  long fpixel_o[2], lpixel_o[2], inc[2]={1,1};
  char region[FLEN_VALUE], regionkey[FLEN_KEYWORD];
  long naxes[]={img->dsize[1], img->dsize[0]}, fpixel_i[2] , lpixel_i[2];


  /* Find the first and last pixel of this crop box from this input
     image. If the outer polygon region is to be kept, then set the
     sides to the image sides.*/
  cropflpixel(crp);
  fpixel_i[0]=crp->fpixel[0];      fpixel_i[1]=crp->fpixel[1];
  lpixel_i[0]=crp->lpixel[0];      lpixel_i[1]=crp->lpixel[1];


  /* Find the overlap and apply it if there is any overlap. */
  if( gal_box_overlap(naxes, fpixel_i, lpixel_i, fpixel_o, lpixel_o) )
    {
      /* Make the output FITS image and initialize it with an array of
         NaN or BLANK values. Note that for FLOAT_IMG and DOUBLE_IMG,
         it will automatically fill them with the NaN value.*/
      if(crp->outfits==NULL)
        firstcropmakearray(crp, fpixel_i, lpixel_i, fpixel_o, lpixel_o);
      ofp=crp->outfits;


      /* Read the desired part of the image, then write it into this
         array. */
      status=0;
      cropsize=(lpixel_i[0]-fpixel_i[0]+1)*(lpixel_i[1]-fpixel_i[1]+1);
      array=gal_data_malloc_array(p->type, cropsize);
      if(fits_read_subset(ifp, gal_fits_type_to_datatype(p->type),
                          fpixel_i, lpixel_i, inc, p->bitnul, array,
                          &anynul, &status))
        gal_fits_io_error(status, NULL);

      /* If we have a floating point or double image, pixels with zero
         value should actually be a NaN. Unless the user specificly
         asks for it, make the conversion.*/
      if(p->zeroisnotblank==0
         && (p->type==GAL_TYPE_FLOAT32
             || p->type==GAL_TYPE_FLOAT64) )
        changezerotonan(array, cropsize, p->type);


      /* If a polygon is given, remove all the pixels within or
         outside of it.*/
      if(p->polygon)
        {
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


      /* A section has been added to the cropped image from this input
         image, so increment crp->imgcount and save the information of
         this image. */
      sprintf(basename, "ICF%zu", ++crp->numimg);
      gal_fits_key_write_filename(basename, img->name, &headers);
      sprintf(regionkey, "%sPIX", basename);
      sprintf(region, "%ld:%ld,%ld:%ld", fpixel_i[0], lpixel_i[0],
              fpixel_i[1], lpixel_i[1]);
      gal_fits_key_add_to_ll_end(&headers, GAL_TYPE_STRING, regionkey,
                                 0, region, 0,
                                 "Range of pixels used for this output.", 0,
                                 NULL);
      gal_fits_key_write(ofp, &headers);

      /* Free the allocated array. */
      free(array);
    }
  else
    if(p->polygon && p->outpolygon==0 && p->mode==IMGCROP_MODE_WCS)
      free(crp->ipolygon);


  return;
}




















/*******************************************************************/
/******************        Check center        *********************/
/*******************************************************************/
int
iscenterfilled(struct onecropparams *crp)
{
  struct cropparams *p=crp->p;

  void *array;
  size_t size, ndim, *dsize;
  fitsfile *ofp=crp->outfits;
  int status=0, anynul=0, type;
  long checkcenter=p->checkcenter;
  long naxes[2], fpixel[2], lpixel[2], inc[2]={1,1};

  /* If checkcenter is zero, then don't check. */
  if(checkcenter==0) return GAL_BLANK_UINT8;

  /* Get the final size of the output image. */
  gal_fits_img_info(ofp, &type, &ndim, &dsize);
  naxes[0]=dsize[1];
  naxes[1]=dsize[0];

  /* Get the size and range of the central region to check. The +1 is
     because in FITS, counting begins from 1, not zero. It might happen
     that the image is actually smaller than the width to check the center
     (for example 1 or 2 pixels wide). In that case, we'll just use the
     full image to check. */
  size = ( (naxes[0]>checkcenter ? checkcenter : naxes[0])
           * (naxes[1]>checkcenter ? checkcenter : naxes[1]) );
  fpixel[0] = naxes[0]>checkcenter ? (naxes[0]/2+1)-checkcenter/2 : 1;
  fpixel[1] = naxes[1]>checkcenter ? (naxes[1]/2+1)-checkcenter/2 : 1;
  lpixel[0] = naxes[0]>checkcenter ? (naxes[0]/2+1)+checkcenter/2 : naxes[0];
  lpixel[1] = naxes[1]>checkcenter ? (naxes[1]/2+1)+checkcenter/2 : naxes[1];

  /* For a check:
  printf("naxes: %ld, %ld\nfpixel: (%ld, %ld)\nlpixel: (%ld, %ld)\n"
         "size: %zu\n", naxes[0], naxes[1], fpixel[0], fpixel[1],
         lpixel[0], lpixel[1], size);
  */

  /* Allocate the array and read in the pixels. */
  array=gal_data_malloc_array(type, size);
  if( fits_read_subset(ofp, gal_fits_type_to_datatype(type), fpixel, lpixel,
                       inc, p->bitnul, array, &anynul, &status) )
    gal_fits_io_error(status, NULL);
  free(array);

  /* CFITSIO already checks if there are any blank pixels. If there are,
     then `anynul' will be 1, if there aren't it will be 0. So the output
     of this function is just the inverse of that number. */
  return !anynul;
}
