/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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

#include <gnuastro/fits.h>
#include <gnuastro/blank.h>

#include "timing.h"
#include "checkset.h"

#include "main.h"

#include "header.h"



/* Print all the extension informations. */
void
fits_print_extension_info(struct fitsparams *p)
{
  uint16_t *ui16;
  fitsfile *fptr;
  gal_data_t *cols=NULL, *tmp;
  char **tstra, **estra, **sstra;
  size_t i, numext, *dsize, ndim;
  int j, nc, numhdu, hdutype, status=0, type;
  char *msg, *tstr, *estr, sstr[1000], extname[FLEN_VALUE];


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
  gal_data_add_to_ll(&cols, NULL, GAL_DATA_TYPE_STRING, 1, &numext, NULL, 1,
                     p->cp.minmapsize, "HDU_SIZE", "name",
                     "Size of image or table number of rows and columns.");
  gal_data_add_to_ll(&cols, NULL, GAL_DATA_TYPE_STRING, 1, &numext, NULL, 1,
                     p->cp.minmapsize, "HDU_TYPE", "name",
                     "Image data type or `table' format (ASCII or binary).");
  gal_data_add_to_ll(&cols, NULL, GAL_DATA_TYPE_STRING, 1, &numext, NULL, 1,
                     p->cp.minmapsize, "EXTNAME", "name",
                     "Extension name of this HDU (EXTNAME in FITS).");
  gal_data_add_to_ll(&cols, NULL, GAL_DATA_TYPE_UINT16, 1, &numext, NULL, 1,
                     p->cp.minmapsize, "HDU_INDEX", "count",
                     "Index (starting from zero) of each HDU (extension).");


  /* Keep pointers to the array of each column for easy writing. */
  ui16  = cols->array;
  estra = cols->next->array;
  tstra = cols->next->next->array;
  sstra = cols->next->next->next->array;

  cols->next->disp_width=15;
  cols->next->next->disp_width=15;


  /* Fill in each column. */
  for(i=0;i<numext;++i)
    {
      /* Work based on the type of the extension. */
      switch(hdutype)
        {
        case IMAGE_HDU:
          gal_fits_img_info(fptr, &type, &ndim, &dsize);
          tstr=gal_data_type_as_string(type , 1);
          break;

        case ASCII_TBL:
        case BINARY_TBL:
          ndim=2;
          tstr = hdutype==ASCII_TBL ? "table_ascii" : "table_binary";
          dsize=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, 2);
          gal_fits_tab_size(fptr, dsize+1, dsize);
          break;

        default:
          error(EXIT_FAILURE, 0, "a bug! the `hdutype' code %d not "
                "recognized", hdutype);
        }


      /* Read the extension name*/
      fits_read_keyword(fptr, "EXTNAME", extname, NULL, &status);
      switch(status)
        {
        case 0:
          estr=gal_fits_key_clean_str_value(extname);
          break;

        case KEY_NO_EXIST:
          sprintf(extname, "%s", GAL_BLANK_STRING);
          estr=extname;
          status=0;
          break;

        default:
          gal_fits_io_error(status, "reading EXTNAME keyword");
        }
      status=0;


      /* Write the size into a string. `sprintf' returns the number of
         written characters (excluding the `\0'). So for each dimension's
         size that is written, we add to `nc' (the number of
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
            case 0: ui16[i]=i;                                   break;
            case 1: gal_checkset_allocate_copy(estr, estra+i);   break;
            case 2: gal_checkset_allocate_copy(tstr, tstra+i);   break;
            case 3: gal_checkset_allocate_copy(sstr, sstra+i);   break;
            }
          ++j;
        }


      /* Move to the next extension if we aren't on the last extension. */
      if( i!=numext-1 && fits_movrel_hdu(fptr, 1, &hdutype, &status) )
        {
          asprintf(&msg, "moving to hdu %zu", i+1);
          gal_fits_io_error(status, msg);
        }
    }


  /* Print the resutls. */
  if(!p->cp.quiet)
    {
      printf("%s\nRun on %s-----\n", PACKAGE_STRING, ctime(&p->rawtime));
      printf("HDU (extension) information: `%s'.\n", p->filename);
      printf(" Column 1: Index (counting from 0).\n");
      printf(" Column 2: Name (`EXTNAME' in FITS standard).\n");
      printf(" Column 3: Image data type or `table' format (ASCII or "
             "binary).\n");
      printf(" Column 4: Size of data in HDU.\n");
      printf("-----\n");
    }
  gal_table_write(cols, NULL, GAL_TABLE_FORMAT_TXT, NULL, 0);
  if(!p->cp.quiet) printf("-----\n");
  gal_data_free_ll(cols);
}





int
fits(struct fitsparams *p)
{
  int r=EXIT_SUCCESS;
  if(p->mode==FITS_MODE_KEY)
    r=header(p);
  else
    fits_print_extension_info(p);
  return r;
}
