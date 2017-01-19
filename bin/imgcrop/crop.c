/*********************************************************************
ImageCrop - Crop a given size from one or multiple images.
ImageCrop is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <gnuastro/polygon.h>

#include <timing.h>
#include <checkset.h>

#include "main.h"

#include "crop.h"
#include "wcsmode.h"























/*******************************************************************/
/************     Set/correct first and last pixel    **************/
/*******************************************************************/
/* Read the section string and set the starting and ending pixels
   based on that. */
void
sectionparser(char *section, long *naxes, long *fpixel, long *lpixel)
{
  int add;
  long read;
  size_t dim=0;
  char *tailptr;
  char forl='f', *pt=section;

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
            error(EXIT_FAILURE, 0, "Extra `,` in `%s`", section);
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
                "includes a float number: %s",
                section);
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

      /* Put it in the correct array. Note that for the last
         point, we don't want to include the given pixel. Unlike
         CFITSIO, in ImageCrop, the given intervals are not
         inclusive. fpixel and lpixel will be directly passed to
         CFITSIO. So we have to correct his here.*/
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

  if(fpixel[0]>lpixel[0] || fpixel[1]>lpixel[1])
    error(EXIT_FAILURE, 0, "the bottom left corner coordinates "
          "cannot be larger or equal to the top right's! Your section "
          "string (%s) has been read as: bottom left coordinate "
          "(%ld, %ld) to top right coordinate (%ld, %ld)",
          section, fpixel[0], fpixel[1], lpixel[0], lpixel[1]);

  /*
  printf("\n%s\nfpixel=(%ld, %ld), lpixel=(%ld, %ld)\n\n", section,
         fpixel[0], fpixel[1], lpixel[0], lpixel[1]);
  exit(0);
  */
}






void
polygonparser(struct imgcropparams *p)
{
  size_t dim=0;
  char *tailptr;
  double read[2], *array;
  struct gal_linkedlist_tdll *gal_linkedlist_tdll=NULL;
  char *pt=p->up.polygon;

  /* Parse the string. */
  while(*pt!='\0')
    {
      switch(*pt)
        {
        case ',':
          ++dim;
          if(dim==2)
            error(EXIT_FAILURE, 0, "Extra `,` in `%s`", p->up.polygon);
          ++pt;
          break;
        case ':':
          if(dim==0)
            error(EXIT_FAILURE, 0, "not enough coordinates for at least "
                  "one polygon vertex (in %s)", p->up.polygon);
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
          read[dim]=strtod(pt, &tailptr);

          /* Check if there actually was a number.
          printf("\n\n------\n%zu: %f (%s)\n", dim, read[dim], tailptr);
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

          /* If this was the second dimension, then put the values
             into the linked list: */
          if(dim==1)
            gal_linkedlist_add_to_tdll(&gal_linkedlist_tdll, read[0],
                                       read[1]);

          /* The job here is done, start from tailptr */
          pt=tailptr;
        }
    }

  /* Convert the linked list to an array: */
  gal_linkedlist_tdll_to_array_inv(gal_linkedlist_tdll, &array,
                                   &p->nvertices);
  if(p->imgmode) { p->ipolygon=array; p->wpolygon=NULL;  }
  else           { p->ipolygon=NULL;  p->wpolygon=array; }

  /* Put them in the proper order in WCS mode: */

  /* For a check:
  {
    size_t i;
    double *polygon=p->imgmode?p->ipolygon:p->wpolygon;
    for(i=0;i<p->nvertices;++i)
      printf("(%f, %f)\n", polygon[i*2], polygon[i*2+1]);
  }
  */

  /* Clean up: */
  gal_linkedlist_free_tdll(gal_linkedlist_tdll);
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





void
polygonmask(struct cropparams *crp, void *array, long *fpixel_i,
            size_t s0, size_t s1)
{
  long *lb, *la=array;
  short *sb, *sa=array;
  float *fb, *fa=array;
  LONGLONG *Lb, *La=array;
  unsigned char *bb, *ba=array;
  double *db, *ipolygon, point[2], *da=array;
  size_t i, *ordinds, size=s0*s1, nvertices=crp->p->nvertices;
  int outpolygon=crp->p->outpolygon, datatype=crp->p->datatype;


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
  switch(datatype)
    {
    case TBYTE:
      bb=gal_fits_datatype_blank(datatype);
      for(i=0;i<size;++i)
        {
          point[0]=i%s1+1; point[1]=i/s1+1;
          if(gal_polygon_pin(ipolygon, point, nvertices)==outpolygon)
            ba[i]=*bb;
        }
      free(bb);
      break;
    case TSHORT:
      sb=gal_fits_datatype_blank(datatype);
      for(i=0;i<size;++i)
        {
          point[0]=i%s1+1; point[1]=i/s1+1;
          if(gal_polygon_pin(ipolygon, point, nvertices)==outpolygon)
            sa[i]=*sb;
        }
      free(sb);
      break;
    case TLONG:
      lb=gal_fits_datatype_blank(datatype);
      for(i=0;i<size;++i)
        {
          point[0]=i%s1+1; point[1]=i/s1+1;
          if(gal_polygon_pin(ipolygon, point, nvertices)==outpolygon)
            la[i]=*lb;
        }
      free(lb);
      break;
    case TLONGLONG:
      Lb=gal_fits_datatype_blank(datatype);
      for(i=0;i<size;++i)
        {
          point[0]=i%s1+1; point[1]=i/s1+1;
          if(gal_polygon_pin(ipolygon, point, nvertices)==outpolygon)
            La[i]=*Lb;
        }
      free(Lb);
      break;
    case TFLOAT:
      fb=gal_fits_datatype_blank(datatype);
      for(i=0;i<size;++i)
        {
          point[0]=i%s1+1; point[1]=i/s1+1;
          if(gal_polygon_pin(ipolygon, point, nvertices)==outpolygon)
            fa[i]=*fb;
        }
      free(fb);
      break;
    case TDOUBLE:
      db=gal_fits_datatype_blank(datatype);
      for(i=0;i<size;++i)
        {
          point[0]=i%s1+1; point[1]=i/s1+1;
          if(gal_polygon_pin(ipolygon, point, nvertices)==outpolygon)
            da[i]=*db;
        }
      free(db);
      break;
    default:
      error(EXIT_FAILURE, 0, "a bug! Please contact us at %s, so we "
            "can fix the problem. For some reason, an unrecognized "
            "datatype value (%d) has been seen in polygonmask (crop.c)",
            PACKAGE_BUGREPORT, datatype);
    }

  /* Clean up: */
  free(ordinds);
  free(ipolygon);
}



















/*******************************************************************/
/******************          One crop.         *********************/
/*******************************************************************/
void
changezerotonan(void *array, size_t size, int bitpix)
{
  float *fp, *ffp;
  double *dp, *fdp;

  switch(bitpix)
    {
    case FLOAT_IMG:
      ffp=(fp=array)+size;
      do
        if(*fp==0.0f) *fp=NAN;
      while(++fp<ffp);
      break;

    case DOUBLE_IMG:
      fdp=(dp=array)+size;
      do
        if(*dp==0.0f) *dp=NAN;
      while(++dp<fdp);
      break;

    default:
      error(EXIT_FAILURE, 0, "in changezerotonan, bitpix is not "
            "recognized! This is out of the users control and is a bug, "
            "please report it to us so we see how it was caused and fix "
            "it");
    }
}






/* Set the output name and image output sizes. */
void
cropname(struct cropparams *crp)
{
  struct imgcropparams *p=crp->p;
  struct gal_commonparams *cp=&p->cp;
  struct imgcroplog *log=&crp->p->log[crp->outindex];

  /* Set the output name and crop sides: */
  if(p->up.catset)
    {
      errno=0;
      log->name=malloc(crp->outlen);
      if(log->name==NULL)
        error(EXIT_FAILURE, errno, "imgmode.c, %zu bytes on "
              "imgcroponthreads", crp->outlen);
      sprintf(log->name, "%s%zu%s", cp->output, crp->outindex+1,
              p->suffix);
      gal_checkset_check_remove_file(log->name, cp->dontdelete);
    }
  else
    {
      /* Set the output name. */
      if(p->outnameisfile)            /* An output file was specified. */
        {
          log->name=cp->output;
          gal_checkset_check_remove_file(log->name, cp->dontdelete);
        }
      else          /* The output was a directory, use automatic output. */
        gal_checkset_automatic_output(p->imgs[crp->imgindex].name,
                                      p->suffix, cp->removedirinfo,
                                      cp->dontdelete, &log->name);
    }
}





/* Find the first and last pixel of a crop from its center point (in
   image mode or WCS mode). */
void
cropflpixel(struct cropparams *crp)
{
  struct imgcropparams *p=crp->p;
  int ncoord=1, nelem=2, status[2]={0,0};
  long *naxes=p->imgs[crp->imgindex].naxes;
  double pixcrd[2], imgcrd[2], phi[1], theta[1];
  long *fpixel=crp->fpixel, *lpixel=crp->lpixel;

  if(p->imgmode)
    {
      if(p->up.catset)
        gal_box_border_from_center(p->cat[crp->outindex*p->cs1+p->xcol],
                                   p->cat[crp->outindex*p->cs1+p->ycol],
                                   p->iwidth, fpixel, lpixel);
      else if(p->up.xcset)
        gal_box_border_from_center(p->xc, p->yc, p->iwidth, fpixel, lpixel);
      else if(p->up.sectionset)
        sectionparser(p->section, naxes, fpixel, lpixel);
      else if(p->up.polygonset)
        {
          if(p->outpolygon==0)
            imgpolygonflpixel(p->ipolygon, p->nvertices, fpixel, lpixel);
        }
      else
        error(EXIT_FAILURE, 0, "a bug! In image mode, neither of the "
              "following has been set: a catalog, a central pixel, "
              "a section or a polygon in the image. Please contact us "
              "to see how it got to this impossible place! You should "
              "have been warned of this condition long before ImageCrop "
              "reaches this point");
    }
  else if(p->wcsmode) /* In wcsmode, crp->world is already filled.       */
    {                 /* Note that p->iwidth was set based on p->wwidth. */
      if(p->up.polygonset)
        { /* Fill crp->ipolygon in wcspolygonpixel, then set flpixel*/
          fillcrpipolygon(crp);
          if(p->outpolygon==0)
            imgpolygonflpixel(crp->ipolygon, p->nvertices, fpixel, lpixel);
        }
      else
        {
          if(wcss2p(p->imgs[crp->imgindex].wcs, ncoord, nelem, crp->world,
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
    }
  else
    error(EXIT_FAILURE, 0, "a bug! in cropflpixel (crop.c), "
          "neither imgmode or wcsmode are set. Please contact us so "
          "we can see how it got to this impossible place");

  /* If the user only wants regions outside to the polygon, then set
     the fpixel and lpixel to cover the full input image. */
  if(p->up.polygonset && p->outpolygon)
    {
      crp->lpixel[0]=naxes[0];
      crp->lpixel[1]=naxes[1];
      crp->fpixel[0]=crp->fpixel[1]=1;
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
firstcropmakearray(struct cropparams *crp, long *fpixel_i,
                   long *lpixel_i, long *fpixel_c, long *lpixel_c)
{
  size_t i;
  fitsfile *ofp;
  long naxes[2];
  double crpix0, crpix1;
  int naxis=2, status=0;
  int bitpix=crp->p->bitpix;
  char *cp, *cpf, blankrec[80], titlerec[80];
  char startblank[]="                      / ";
  char *outname=crp->p->log[crp->outindex].name;
  struct inputimgs *img=&crp->p->imgs[crp->imgindex];


  /* Set the last element of the blank array. */
  cpf=blankrec+79;
  *cpf='\0';
  titlerec[79]='\0';
  cp=blankrec; do *cp=' '; while(++cp<cpf);


  /* Set the size of the output, in WCS mode, noblank==0. */
  if(crp->p->noblank && crp->p->wcsmode==0)
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


  /* Create the FITS image extension and array and fill it with null
     values. */
  if(fits_create_file(&crp->outfits, outname, &status))
    gal_fits_io_error(status, "creating file");
  ofp=crp->outfits;
  if(fits_create_img(ofp, bitpix, naxis, naxes, &status))
    gal_fits_io_error(status, "creating image");
  if(bitpix==BYTE_IMG || bitpix==SHORT_IMG
     || bitpix==LONG_IMG || bitpix==LONGLONG_IMG)
    if(fits_write_key(ofp, crp->p->datatype, "BLANK",
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





/* The starting and ending points are set in the cropparams structure
   for one crop from one image. Crop that region out.

   return values are:
   0: No crop was made (not in range of image).
   1: The input image covered at least part of the crop image.
 */
void
onecrop(struct cropparams *crp)
{
  struct imgcropparams *p=crp->p;
  struct inputimgs *img=&p->imgs[crp->imgindex];

  void *array;
  size_t cropsize;
  char basename[FLEN_KEYWORD];
  long fpixel_i[2] , lpixel_i[2];
  fitsfile *ifp=crp->infits, *ofp;
  struct gal_fits_key_ll *headers=NULL;
  int status=0, anynul=0, bitpix=p->bitpix;
  long fpixel_o[2], lpixel_o[2], inc[2]={1,1};
  char region[FLEN_VALUE], regionkey[FLEN_KEYWORD];


  /* Find the first and last pixel of this crop box from this input
     image. If the outer polygon region is to be kept, then set the
     sides to the image sides.*/
  cropflpixel(crp);
  fpixel_i[0]=crp->fpixel[0];      fpixel_i[1]=crp->fpixel[1];
  lpixel_i[0]=crp->lpixel[0];      lpixel_i[1]=crp->lpixel[1];


  /* Find the overlap and apply it if there is any overlap. */
  if( gal_box_overlap(img->naxes, fpixel_i, lpixel_i, fpixel_o, lpixel_o) )
    {
      /* Make the output FITS image and initialize it with an array of
         NaN or BLANK values. Note that for FLOAT_IMG and DOUBLE_IMG,
         it will automatically fill them with the NaN value.*/
      if(crp->outfits==NULL)
        firstcropmakearray(crp, fpixel_i, lpixel_i, fpixel_o, lpixel_o);
      ofp=crp->outfits;


      /* Read the desired part of the image, then write it into this
         array. */
      cropsize=(lpixel_i[0]-fpixel_i[0]+1)*(lpixel_i[1]-fpixel_i[1]+1);
      array=gal_fits_datatype_alloc(cropsize, p->datatype);
      status=0;
      if(fits_read_subset(ifp, p->datatype, fpixel_i, lpixel_i, inc,
                          p->bitnul, array, &anynul, &status))
        gal_fits_io_error(status, NULL);


      /* If we have a floating point or double image, pixels with zero
         value should actually be a NaN. Unless the user specificly
         asks for it, make the conversion.*/
      if(p->zeroisnotblank==0
         && (bitpix==FLOAT_IMG || bitpix==DOUBLE_IMG) )
        changezerotonan(array, cropsize, bitpix);


      /* If a polygon is given, remove all the pixels within or
         outside of it.*/
      if(p->up.polygonset)
        {
          /* In WCS mode, crp->ipolygon was allocated and filled in
             wcspolygonflpixel (wcsmode.c). */
          if(p->imgmode) crp->ipolygon=p->ipolygon;
          polygonmask(crp, array, fpixel_i, lpixel_i[1]-fpixel_i[1]+1,
                      lpixel_i[0]-fpixel_i[0]+1);
          if(p->wcsmode) free(crp->ipolygon);
        }


      /* Write the array into the image. */
      status=0;
      if( fits_write_subset(ofp, p->datatype, fpixel_o, lpixel_o,
                            array, &status) )
        gal_fits_io_error(status, NULL);


      /* A section has been added to the cropped image from this input
         image, so increment crp->imgcount and save the information of
         this image. */
      sprintf(basename, "ICF%zu", ++p->log[crp->outindex].numimg);
      gal_fits_file_name_in_keywords(basename, img->name, &headers);
      sprintf(regionkey, "%sPIX", basename);
      sprintf(region, "%ld:%ld,%ld:%ld", fpixel_i[0], lpixel_i[0],
              fpixel_i[1], lpixel_i[1]);
      gal_fits_add_to_key_ll_end(&headers, TSTRING, regionkey, 0, region, 0,
                                 "Range of pixels used for this output.", 0,
                                 NULL);
      gal_fits_update_keys(ofp, &headers);


      /* Free the allocated array. */
      free(array);
    }
  else
    if(p->up.polygonset && p->outpolygon==0 && p->wcsmode)
      free(crp->ipolygon);


  return;
}




















/*******************************************************************/
/******************        Check center        *********************/
/*******************************************************************/
int
iscenterfilled(struct cropparams *crp)
{
  struct imgcropparams *p=crp->p;

  void *array;
  int bitpix=p->bitpix;
  size_t size, nulcount;
  fitsfile *ofp=crp->outfits;
  int status=0, maxdim=10, anynul;
  long checkcenter=p->checkcenter;
  long naxes[2], fpixel[2], lpixel[2], inc[2]={1,1};

  uint8_t *b, *fb, *nb;
  int16_t *s, *fs, *ns;
  int32_t *l, *fl, *nl;
  int64_t *L, *fL, *nL;
  float   *f, *ff; /* isnan will check. */
  double  *d, *fd; /* isnan will check */

  /* Get the final size of the output image. */
  if( fits_get_img_size(ofp, maxdim, naxes, &status) )
    gal_fits_io_error(status, NULL);

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

  /* Allocate the array and read in the pixels. */
  array=gal_fits_datatype_alloc(size, gal_fits_bitpix_to_datatype(bitpix) );
  if( fits_read_subset(ofp, p->datatype, fpixel, lpixel, inc,
                       p->bitnul, array, &anynul, &status) )
    gal_fits_io_error(status, NULL);

  /* Depending on bitpix, check the central pixels of the image. */
  nulcount=0;
  switch(bitpix)
    {
    case BYTE_IMG:
      nb=p->bitnul;
      fb=(b=array)+size;
      do if(*b==*nb) ++nulcount; while(++b<fb);
      break;

    case SHORT_IMG:
      ns=p->bitnul;
      fs=(s=array)+size;
      do if(*s==*ns) ++nulcount; while(++s<fs);
      break;

    case LONG_IMG:
      nl=p->bitnul;
      fl=(l=array)+size;
      do if(*l==*nl) ++nulcount; while(++l<fl);
      break;

    case LONGLONG_IMG:
      nL=p->bitnul;
      fL=(L=array)+size;
      do if(*L==*nL) ++nulcount; while(++L<fL);
      break;

    case FLOAT_IMG:
      ff=(f=array)+size;
      do if(isnan(*f)) ++nulcount; while(++f<ff);
      break;

    case DOUBLE_IMG:
      fd=(d=array)+size;
      do if(isnan(*d)) ++nulcount; while(++d<fd);
      break;

    default:
      error(EXIT_FAILURE, 0, "in iscenterfilled, the bitbix is not "
            "recognized! This is not possible by the user, so it is a "
            "a bug. Please contact us so we can correct it");
    }
  free(array);

  if(nulcount==size)
    return 0;
  else
    return 1;
}



















/*******************************************************************/
/******************          Log file          *********************/
/*******************************************************************/
void
printlog(struct imgcropparams *p)
{
  size_t i;
  FILE *logfile;
  char msg[GAL_TIMING_VERB_MSG_LENGTH_V];
  struct imgcroplog *log=p->log;
  size_t numfiles=0, numcentfilled=0, numstitched=0;

  /* Only for a catalog are these statistics worth it! */
  if(p->up.catset && p->cp.verb)
    for(i=0;log[i].name;++i)
      if(log[i].numimg)
        {
          if(log[i].centerfilled || p->keepblankcenter)
            {
              ++numfiles;
              if(log[i].numimg>1)
                ++numstitched;
            }
          if(log[i].centerfilled)
            ++numcentfilled;
        }

  /* Check to see if the file exists and remove if if it is ok. */
  gal_checkset_check_remove_file(LOGFILENAME, p->cp.dontdelete);

  /* Make the file and print the top comments. If the file can't be
     opened for write mode, there is no problem, this is a log file,
     the user can set it on verbose mode and the same information will
     be printed. */
  errno=0;
  logfile=fopen(LOGFILENAME, "w");
  if(logfile)
    {
      /* First print the comments to the file. */
      fprintf(logfile,
              "# "SPACK_STRING" log file.\n"
              "# "SPACK_NAME" was run on %s#\n",
              ctime(&p->rawtime));
      if(p->keepblankcenter==0)
        fprintf(logfile, "# NOTE: by default images with a blank "
                "center are deleted.\n# To keep such images, run again "
                "with `--keepblankcenter`.\n#\n");
      fprintf(logfile,
              "# Column numbers below start from zero.\n"
              "# 0: Output file name.\n"
              "# 1: Number of images used in this cropped image.\n"
              "# 2: Are the central %zu pixels filled? (1: yes, 0: no)\n",
              p->checkcenter);

      /* Then print each output's information. */
      for(i=0;log[i].name;++i)
        fprintf(logfile, "%s     %-8zu%-2d\n", log[i].name,
                log[i].numimg, log[i].centerfilled);

      /* Report Summary: */
      if(p->cp.verb && p->up.catset)
        {
          sprintf(msg, "%zu images created.", numfiles);
          gal_timing_report(NULL, msg, 1);
          sprintf(msg, "%zu filled in the center.",
                  numcentfilled);
          gal_timing_report(NULL, msg, 1);
          if(numstitched)
            {
              sprintf(msg, "%zu used more than one input.",
                      numstitched);
              gal_timing_report(NULL, msg, 1);
            }
        }

      /* Close the file. */
      errno=0;
      if(fclose(logfile))
        error(EXIT_FAILURE, errno, LOGFILENAME" could not be closed");
    }
}
