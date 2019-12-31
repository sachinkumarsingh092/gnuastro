/*********************************************************************
array -- Functions for I/O on arrays (images or cubes)
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/txt.h>
#include <gnuastro/fits.h>
#include <gnuastro/jpeg.h>
#include <gnuastro/tiff.h>
#include <gnuastro/array.h>










/*********************************************************************/
/*****************        High-level functions        ****************/
/*********************************************************************/
int
gal_array_name_recognized(char *name)
{
  if( gal_array_name_recognized_multiext(name) ) return 1;
  else if ( gal_jpeg_name_is_jpeg(name)        ) return 1;
  else                                           return 0;

  /* Control should not get to here, but just to avoid compiler warnings,
     we'll return a NULL. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to solve the "
        "problem. Control must not reach the end of this function", __func__,
        PACKAGE_BUGREPORT);
  return 0;
}





int
gal_array_name_recognized_multiext(char *name)
{
  if(       gal_fits_name_is_fits(name) ) return 1;
  else if ( gal_tiff_name_is_tiff(name) ) return 1;
  else                                    return 0;

  /* Control should not get to here, but just to avoid compiler warnings,
     we'll return a NULL. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to solve the "
        "problem. Control must not reach the end of this function", __func__,
        PACKAGE_BUGREPORT);
  return 0;
}





/* Read (all the possibly existing) color channels within each
   extension/dir of the given file. */
gal_data_t *
gal_array_read(char *filename, char *extension, gal_list_str_t *lines,
               size_t minmapsize, int quietmmap)
{
  size_t ext;

  /* FITS  */
  if( gal_fits_name_is_fits(filename) )
    return gal_fits_img_read(filename, extension, minmapsize, quietmmap);

  /* TIFF */
  else if ( gal_tiff_name_is_tiff(filename) )
    {
      ext=gal_tiff_dir_string_read(extension);
      return gal_tiff_read(filename, ext, minmapsize, quietmmap);
    }

  /* JPEG */
  else if ( gal_jpeg_name_is_jpeg(filename) )
    return gal_jpeg_read(filename, minmapsize, quietmmap);

  /* Default: plain text. */
  else
    return gal_txt_image_read(filename, lines, minmapsize, quietmmap);

  /* Control should not get to here, but just to avoid compiler warnings,
     we'll return a NULL. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to solve the "
        "problem. Control must not reach the end of this function", __func__,
        PACKAGE_BUGREPORT);
  return NULL;
}





/* Read the contents of the given file/extension to a specific type. */
gal_data_t *
gal_array_read_to_type(char *filename, char *extension,
                       gal_list_str_t *lines, uint8_t type,
                       size_t minmapsize, int quietmmap)
{
  gal_data_t *out=NULL;
  gal_data_t *next, *in=gal_array_read(filename, extension, lines,
                                       minmapsize, quietmmap);

  /* Go over all the channels. */
  while(in)
    {
      next=in->next;
      in->next=NULL;
      gal_list_data_add(&out, gal_data_copy_to_new_type_free(in, type));
      in=next;
    }

  /* Invert the reverse list and return. */
  gal_list_data_reverse(&out);
  return out;
}





/* Read the input array and make sure it is only one channel. */
gal_data_t *
gal_array_read_one_ch(char *filename, char *extension, gal_list_str_t *lines,
                      size_t minmapsize, int quietmmap)
{
  char *fname;
  gal_data_t *out;
  out=gal_array_read(filename, extension, lines, minmapsize, quietmmap);

  if(out->next)
    {
      if(extension)
        {
          if( asprintf(&fname, "%s (hdu %s)", filename, extension)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation error", __func__);
        }
      else
        fname=filename;

      error(EXIT_FAILURE, 0, "%s: contains %zu channels (it isn't "
            "monochrome).\n\n"
            "You can use Gnuastro's ConvertType program to separate the "
            "(color) channels into separate extensions of a FITS file, with "
            "a command like this:\n\n"
            "    $ astconvertt %s -h%s --output=sep-ch.fits",
            fname, gal_list_data_number(out), filename, extension);
    }

  return out;
}





/* Read a single-channel dataset into a specific type. */
gal_data_t *
gal_array_read_one_ch_to_type(char *filename, char *extension,
                              gal_list_str_t *lines, uint8_t type,
                              size_t minmapsize, int quietmmap)
{
  gal_data_t *out=gal_array_read_one_ch(filename, extension, lines,
                                        minmapsize, quietmmap);

  return gal_data_copy_to_new_type_free(out, type);
}
