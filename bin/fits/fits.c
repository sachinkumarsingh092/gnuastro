/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2021, Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>
#include <string.h>
#include <unistd.h>

#include <gnuastro/wcs.h>
#include <gnuastro/list.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/pointer.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "fits.h"
#include "keywords.h"



int
fits_has_error(struct fitsparams *p, int actioncode, char *string, int status)
{
  char *action=NULL;
  int r=EXIT_SUCCESS;

  switch(actioncode)
    {
    case FITS_ACTION_DELETE:        action="deleted";      break;
    case FITS_ACTION_RENAME:        action="renamed";      break;
    case FITS_ACTION_UPDATE:        action="updated";      break;
    case FITS_ACTION_WRITE:         action="written";      break;
    case FITS_ACTION_COPY:          action="copied";       break;
    case FITS_ACTION_REMOVE:        action="removed";      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' so we "
            "can fix this problem. The value of 'actioncode' must not be %d",
            __func__, PACKAGE_BUGREPORT, actioncode);
    }

  if(p->quitonerror)
    {
      fits_report_error(stderr, status);
      error(EXIT_FAILURE, 0, "%s: %s: not %s\n", __func__, string, action);
    }
  else
    {
      fprintf(stderr, "%s: Not %s.\n", string, action);
      r=EXIT_FAILURE;
    }
  return r;
}





/* Print all the extension informations. */
void
fits_print_extension_info(struct fitsparams *p)
{
  uint16_t *ui16;
  fitsfile *fptr;
  size_t i, numext, *dsize, ndim;
  int hascomments=0, hasblankname=0;
  char **tstra, **estra, **sstra, **cmstra;
  int j, nc, numhdu, hdutype, status=0, type;
  gal_data_t *tmp, *cols=NULL, *sizecol, *commentscol;
  char *msg, *tstr=NULL, sstr[1000], extname[FLEN_VALUE];


  /* Open the FITS file and read the first extension type, upon moving to
     the next extension, we will read its type, so for the first we will
     need to do it explicitly. */
  fptr=gal_fits_hdu_open(p->filename, "0", READONLY);
  if (fits_get_hdu_type(fptr, &hdutype, &status) )
    gal_fits_io_error(status, "reading first extension");


  /* Get the number of HDUs. */
  if( fits_get_num_hdus(fptr, &numhdu, &status) )
    gal_fits_io_error(status, "finding number of HDUs");
  numext=numhdu;


  /* Allocate all the columns (in reverse order, since this is a simple
     linked list). */
  gal_list_data_add_alloc(&cols, NULL, GAL_TYPE_STRING, 1, &numext, NULL, 1,
                          -1, 1, "HDU_COMMENT", "note", "Possible comment");
  gal_list_data_add_alloc(&cols, NULL, GAL_TYPE_STRING, 1, &numext, NULL, 1,
                          -1, 1, "HDU_SIZE", "name", "Size of image or table "
                          "number of rows and columns.");
  gal_list_data_add_alloc(&cols, NULL, GAL_TYPE_STRING, 1, &numext, NULL, 1,
                          -1, 1, "HDU_TYPE", "name", "Image data type or "
                          "'table' format (ASCII or binary).");
  gal_list_data_add_alloc(&cols, NULL, GAL_TYPE_STRING, 1, &numext, NULL, 1,
                          -1, 1, "EXTNAME", "name", "Extension name of this "
                          "HDU (EXTNAME in FITS).");
  gal_list_data_add_alloc(&cols, NULL, GAL_TYPE_UINT16, 1, &numext, NULL, 1,
                          -1, 1, "HDU_INDEX", "count", "Index (starting from "
                          "zero) of each HDU (extension).");


  /* Keep pointers to the array of each column for easy writing. */
  ui16  = cols->array;
  estra = cols->next->array;
  tstra = cols->next->next->array;

  sizecol = cols->next->next->next;
  sstra = sizecol->array;

  commentscol = cols->next->next->next->next;
  cmstra = commentscol->array;

  cols->next->disp_width=15;
  cols->next->next->disp_width=15;


  /* Fill in each column. */
  for(i=0;i<numext;++i)
    {
      /* Work based on the type of the extension. */
      switch(hdutype)
        {
        case IMAGE_HDU:
          gal_fits_img_info(fptr, &type, &ndim, &dsize, NULL, NULL);
          tstr = ndim==0 ? "no-data" : gal_type_name(type , 1);
          break;

        case ASCII_TBL:
        case BINARY_TBL:
          ndim=2;
          tstr = hdutype==ASCII_TBL ? "table_ascii" : "table_binary";
          dsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, 2, 0, __func__,
                                     "dsize");
          gal_fits_tab_size(fptr, dsize+1, dsize);
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! the 'hdutype' code %d not "
                "recognized", __func__, hdutype);
        }


      /* Read the extension name*/
      fits_read_keyword(fptr, "EXTNAME", extname, NULL, &status);
      switch(status)
        {
        case 0:
          gal_fits_key_clean_str_value(extname);
          break;

        case KEY_NO_EXIST:
          sprintf(extname, "%s", GAL_BLANK_STRING);
          hasblankname=1;
          status=0;
          break;

        default:
          gal_fits_io_error(status, "reading EXTNAME keyword");
        }
      status=0;


      /* Check if its a healpix grid. */
      if( gal_fits_hdu_is_healpix(fptr) )
        {
          hascomments=1;
          gal_checkset_allocate_copy("HEALpix", cmstra+i);
        }
      else
        gal_checkset_allocate_copy(GAL_BLANK_STRING, cmstra+i);


      /* Write the size into a string. 'sprintf' returns the number of
         written characters (excluding the '\0'). So for each dimension's
         size that is written, we add to 'nc' (the number of
         characters). Note that FITS allows blank extensions, in those
         cases, return "0". */
      if(ndim>0)
        {
          nc=0;
          for(j=ndim-1;j>=0;--j)
            nc += sprintf(sstr+nc, "%zux", dsize[j]);
          sstr[nc-1]='\0';
          free(dsize);
        }
      else
        {
          sstr[0]='0';
          sstr[1]='\0';
        }


      /* Write the strings into the columns. */
      j=0;
      for(tmp=cols; tmp!=NULL; tmp=tmp->next)
        {
          switch(j)
            {
            case 0: ui16[i]=i;                                    break;
            case 1: gal_checkset_allocate_copy(extname, estra+i); break;
            case 2: gal_checkset_allocate_copy(tstr, tstra+i);    break;
            case 3: gal_checkset_allocate_copy(sstr, sstra+i);    break;
            }
          ++j;
        }


      /* Move to the next extension if we aren't on the last extension. */
      if( i!=numext-1 && fits_movrel_hdu(fptr, 1, &hdutype, &status) )
        {
          if( asprintf(&msg, "moving to hdu %zu", i+1)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_fits_io_error(status, msg);
        }
    }

  /* Close the file. */
  fits_close_file(fptr, &status);

  /* If there weren't any comments, don't print the comment column. */
  if(hascomments==0)
    {
      sizecol->next=NULL;
      gal_data_free(commentscol);
    }

  /* Print the results. */
  if(!p->cp.quiet)
    {
      printf("%s\nRun on %s-----\n", PROGRAM_STRING, ctime(&p->rawtime));
      printf("HDU (extension) information: '%s'.\n", p->filename);
      printf(" Column 1: Index (counting from 0, usable with '--hdu').\n");
      printf(" Column 2: Name ('EXTNAME' in FITS standard, usable with "
             "'--hdu').\n");
      if(hasblankname)
        printf("           ('%s' means that no name is specified for this "
               "HDU)\n", GAL_BLANK_STRING);
      printf(" Column 3: Image data type or 'table' format (ASCII or "
             "binary).\n");
      printf(" Column 4: Size of data in HDU.\n");
      if(hascomments)
        printf(" Column 5: Comments about the HDU (e.g., if its HEALpix, or "
               "etc).\n");
      printf("-----\n");
    }
  gal_table_write(cols, NULL, NULL, GAL_TABLE_FORMAT_TXT, NULL, NULL, 0);
  gal_list_data_free(cols);
}





static void
fits_hdu_number(struct fitsparams *p)
{
  fitsfile *fptr;
  int numhdu, status=0;

  /* Read the first extension (necessary for reading the rest). */
  fptr=gal_fits_hdu_open(p->filename, "0", READONLY);

  /* Get the number of HDUs. */
  if( fits_get_num_hdus(fptr, &numhdu, &status) )
    gal_fits_io_error(status, "finding number of HDUs");

  /* Close the file. */
  fits_close_file(fptr, &status);

  /* Print the result. */
  printf("%d\n", numhdu);
}





static void
fits_datasum(struct fitsparams *p)
{
  printf("%ld\n", gal_fits_hdu_datasum(p->filename, p->cp.hdu));
}





static void
fits_pixelscale(struct fitsparams *p)
{
  int nwcs=0;
  size_t i, ndim=0;
  struct wcsprm *wcs;
  double multip, *pixelscale;

  /* Read the desired WCS. */
  wcs=gal_wcs_read(p->filename, p->cp.hdu, 0, 0, &nwcs);

  /* If a WCS doesn't exist, let the user know and return. */
  if(wcs)
    ndim=wcs->naxis;
  else
    error(EXIT_FAILURE, 0, "%s (hdu %s): no WCS could be read by WCSLIB, "
          "hence the pixel-scale cannot be determined", p->filename,
          p->cp.hdu);

  /* Calculate the pixel-scale in each dimension. */
  pixelscale=gal_wcs_pixel_scale(wcs);

  /* If not in quiet-mode, print some extra information. We don't want the
     last number to have a space after it, so we'll write the last one
     outside the loop.*/
  if(p->cp.quiet==0)
    {
      printf("Basic information for --pixelscale (remove extra info "
             "with '--quiet' or '-q')\n");
      printf("  Input: %s (hdu %s) has %zu dimensions.\n", p->filename,
             p->cp.hdu, ndim);
      printf("  Pixel scale in each FITS dimension:\n");
      for(i=0;i<ndim;++i)
        {
          if( !strcmp(wcs->cunit[i], "deg") )
            printf("    %zu: %g (%s/pixel) = %g (arcsec/pixel)\n", i+1,
                   pixelscale[i], wcs->cunit[i], pixelscale[i]*3600);
          else
            printf("    %zu: %g (%s/slice)\n", i+1,
                   pixelscale[i], wcs->cunit[i]);
        }

      /* Pixel area/volume. */
      if(ndim>=2)
        {
          /* Multiply the values in each dimension. */
          multip=pixelscale[0]*pixelscale[1];

          /* Pixel scale (applicable to 2 or 3 dimensions). */
          printf("  Pixel area%s:\n", ndim==2?"":" (on each 2D slice) ");
          if( strcmp(wcs->cunit[0], wcs->cunit[1]) )
            printf("    %g (%s*%s)\n", multip, wcs->cunit[0],
                   wcs->cunit[1]);
          else
            if( strcmp(wcs->cunit[0], "deg") )
              printf("    %g (%s^2)\n", multip, wcs->cunit[0]);
            else
              printf("    %g (deg^2) = %g (arcsec^2)\n",
                     multip, multip*3600*3600 );

          /* For a 3 dimensional dataset, we need need extra info. */
          if(ndim>=3)
            {
              multip*=pixelscale[2];
              printf("  Voxel volume:\n");
              if( strcmp(wcs->cunit[0], wcs->cunit[1]) )
                printf("    %g (%s*%s*%s)\n", multip, wcs->cunit[0],
                       wcs->cunit[1], wcs->cunit[2]);
              else
                if( strcmp(wcs->cunit[0], "deg") )
                  printf("    %g (%s^2*%s)\n", multip, wcs->cunit[0],
                         wcs->cunit[2]);
                else
                  {
                    if( strcmp(wcs->cunit[2], "m") )
                      printf("    %g (deg^2*%s) = %g (arcsec^2*%s)\n",
                             multip, wcs->cunit[2], multip*3600*3600,
                             wcs->cunit[2]);
                    else
                      printf("    %g (deg^2*m) = %g (arcsec^2*m) = "
                             "%g (arcsec^2*A)\n", multip, multip*3600*3600,
                             multip*3600*3600*1e10);
                  }
            }
        }
    }
  else
    {
      multip=1;
      for(i=0;i<ndim;++i)
        {
          multip*=pixelscale[i];
          printf("%g ", pixelscale[i]);
        }
      switch(ndim)
        {
        case 2: printf("%g\n", multip); break;
        case 3: printf("%g %g\n", pixelscale[0]*pixelscale[1], multip); break;
        }
    }

  /* Clean up. */
  wcsfree(wcs);
  free(pixelscale);
}





static void
fits_skycoverage(struct fitsparams *p)
{
  int nwcs=0;
  size_t i, ndim;
  struct wcsprm *wcs;
  double *center, *width, *min, *max;

  /* Find the coverage. */
  if( gal_wcs_coverage(p->filename, p->cp.hdu, &ndim,
                       &center, &width, &min, &max)==0 )
    error(EXIT_FAILURE, 0, "%s (hdu %s): is not usable for finding "
          "sky coverage (either doesn't have a WCS, or isn't an image "
          "or cube HDU with 2 or 3 dimensions", p->filename, p->cp.hdu);

  /* Inform the user. */
  if(p->cp.quiet)
    {
      /* First print the center and full-width. */
      for(i=0;i<ndim;++i) printf("%-15.10g ", center[i]);
      for(i=0;i<ndim;++i) printf("%-15.10g ", width[i]);
      printf("\n");

      /* Then print the range in coordinates. */
      for(i=0;i<ndim;++i) printf("%-15.10g %-15.10g ", min[i], max[i]);
      printf("\n");
    }
  else
    {
      printf("Input file: %s (hdu: %s)\n", p->filename, p->cp.hdu);
      printf("\nSky coverage by center and (full) width:\n");
      switch(ndim)
        {
        case 2:
          printf("  Center: %-15.10g%-15.10g\n", center[0], center[1]);
          printf("  Width:  %-15.10g%-15.10g\n", width[0],  width[1]);
          break;
        case 3:
          printf("  Center: %-15.10g%-15.10g%-15.10g\n", center[0],
                 center[1], center[2]);
          printf("  width:  %-15.10g%-15.10g%-15.10g\n", width[0],
                 width[1], width[2]);
          break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
                "the problem. 'ndim' value %zu is not recognized", __func__,
                PACKAGE_BUGREPORT, ndim);
        }

      /* For the range type of coverage. */
      wcs=gal_wcs_read(p->filename, p->cp.hdu, 0, 0, &nwcs);
      printf("\nSky coverage by range along dimensions:\n");
      for(i=0;i<ndim;++i)
        printf("  %-8s %-15.10g%-15.10g\n", gal_wcs_dimension_name(wcs, i),
               min[i], max[i]);
      wcsfree(wcs);
    }

  /* Clean up. */
  free(min);
  free(max);
  free(width);
  free(center);
}





static void
fits_hdu_remove(struct fitsparams *p, int *r)
{
  char *hdu;
  fitsfile *fptr;
  int status=0, hdutype;

  while(p->remove)
    {
      /* Pop-out the top element. */
      hdu=gal_list_str_pop(&p->remove);

      /* Open the FITS file at the specified HDU. */
      fptr=gal_fits_hdu_open(p->filename, hdu, READWRITE);

      /* Delete the extension. */
      if( fits_delete_hdu(fptr, &hdutype, &status) )
        *r=fits_has_error(p, FITS_ACTION_REMOVE, hdu, status);
      status=0;

      /* Close the file. */
      fits_close_file(fptr, &status);
    }
}





/* This is similar to the library's 'gal_fits_open_to_write', except that
   it won't create an empty first extension.*/
fitsfile *
fits_open_to_write_no_blank(char *filename)
{
  int status=0;
  fitsfile *fptr;

  /* When the file exists just open it. Otherwise, create the file. But we
     want to leave the first extension as a blank extension and put the
     image in the next extension to be consistent between tables and
     images. */
  if(access(filename,F_OK) == -1 )
    {
      /* Create the file. */
      if( fits_create_file(&fptr, filename, &status) )
        gal_fits_io_error(status, NULL);
    }

  /* Open the file, ready for later steps. */
  if( fits_open_file(&fptr, filename, READWRITE, &status) )
    gal_fits_io_error(status, NULL);

  /* Return the pointer. */
  return fptr;
}





static void
fits_hdu_copy(struct fitsparams *p, int cut1_copy0, int *r)
{
  char *hdu;
  int status=0, hdutype;
  fitsfile *in, *out=NULL;
  gal_list_str_t *list = cut1_copy0 ? p->cut : p->copy;

  /* Copy all the given extensions. */
  while(list)
    {
      /* Pop-out the top element. */
      hdu=gal_list_str_pop(&list);

      /* Open the FITS file at the specified HDU. */
      in=gal_fits_hdu_open(p->filename, hdu,
                           cut1_copy0 ? READWRITE : READONLY);

      /* If the output isn't opened yet, open it.  */
      if(out==NULL)
        out = ( ( gal_fits_hdu_format(p->filename, hdu)==IMAGE_HDU
                  && p->primaryimghdu )
                ? fits_open_to_write_no_blank(p->cp.output)
                : gal_fits_open_to_write(p->cp.output) );


      /* Copy to the extension. */
      if( fits_copy_hdu(in, out, 0, &status) )
        *r=fits_has_error(p, FITS_ACTION_COPY, hdu, status);
      status=0;

      /* If this is a 'cut' operation, then remove the extension. */
      if(cut1_copy0)
        {
          if( fits_delete_hdu(in, &hdutype, &status) )
            *r=fits_has_error(p, FITS_ACTION_REMOVE, hdu, status);
          status=0;
        }

      /* Close the input file. */
      fits_close_file(in, &status);
    }

  /* Close the output file. */
  if(out) fits_close_file(out, &status);
}





int
fits(struct fitsparams *p)
{
  int r=EXIT_SUCCESS, printhduinfo=1;

  switch(p->mode)
    {
    /* Keywords, we have a separate set of functions in 'keywords.c'. */
    case FITS_MODE_KEY:
      r=keywords(p);
      break;

    /* HDU, functions defined here. */
    case FITS_MODE_HDU:

      /* Options that must be called alone. */
      if(p->numhdus) fits_hdu_number(p);
      else if(p->datasum) fits_datasum(p);
      else if(p->pixelscale) fits_pixelscale(p);
      else if(p->skycoverage) fits_skycoverage(p);

      /* Options that can be called together. */
      else
        {
          if(p->copy)
            {
              fits_hdu_copy(p, 0, &r);
              printhduinfo=0;
            }
          if(p->cut)
            {
              fits_hdu_copy(p, 1, &r);
              printhduinfo=0;
            }
          if(p->remove)
            {
              fits_hdu_remove(p, &r);
              printhduinfo=0;
            }

          if(printhduinfo)
            fits_print_extension_info(p);
        }
      break;

    /* Not recognized. */
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to address "
            "the problem. The code %d is not recognized for p->mode",
            __func__, PACKAGE_BUGREPORT, p->mode);
    }

  return r;
}
