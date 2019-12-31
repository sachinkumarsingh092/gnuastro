/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef MAIN_H
#define MAIN_H

/* Include necessary headers */
#include <gnuastro/data.h>

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "ConvertType"   /* Program full name.       */
#define PROGRAM_EXEC   "astconvertt"   /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION


/* Constants/macros. */
#define BLANK_CHANNEL_NAME "blank"


enum output_formats
{
  OUT_FORMAT_INVALID,           /* ==0 by C standard */

  OUT_FORMAT_TXT,
  OUT_FORMAT_EPS,
  OUT_FORMAT_PDF,
  OUT_FORMAT_FITS,
  OUT_FORMAT_JPEG,
};

enum colorspace_names
{
  COLOR_INVALID,           /* ==0 by C standard */

  COLOR_RGB,
  COLOR_HSV,
  COLOR_SLS,
  COLOR_VIRIDIS,
  COLOR_GRAY,
};





/* This is used in converting certain values in the array. */
struct change
{
  gal_data_t    *from;
  gal_data_t      *to;
  struct change *next;
};





struct converttparams
{
  /* From command-line */
  struct gal_options_common_params cp;   /* Common parameters.         */
  gal_list_str_t          *inputnames;   /* The names of input files.  */
  gal_list_str_t                *hdus;   /* The names of input hdus.   */
  uint8_t           quality;  /* Quality of JPEG image.                */
  float           widthincm;  /* Width in centimeters.                 */
  uint32_t      borderwidth;  /* Width of border in PostScript points. */
  uint8_t               hex;  /* Use hexadecimal not ASCII85 encoding. */
  gal_data_t      *colormap;  /* Color space to use for single/mono.   */
  uint8_t          rgbtohsv;  /* Convert input RGB channels to HSV.    */
  char          *fluxlowstr;  /* Lower flux truncation value.          */
  char         *fluxhighstr;  /* Higher flux truncation value.         */
  uint8_t           maxbyte;  /* Maximum byte value.                   */
  uint8_t          forcemin;  /* fluxlow is minimum.                   */
  uint8_t          forcemax;  /* fluxhigh is maximum.                  */
  char           *changestr;  /* String of change values.              */
  uint8_t  changeaftertrunc;  /* First convert, then truncate.         */
  uint8_t            invert;  /* ==1: invert the output image.         */

  /* Internal */
  struct change     *change;  /* The value conversion string.          */
  gal_data_t       *fluxlow;  /* The lower flux truncation.            */
  gal_data_t      *fluxhigh;  /* The higher flux truncation.           */
  time_t            rawtime;  /* Starting time of the program.         */
  int             outformat;  /* The format of the output file.        */
  size_t              numch;  /* Current Channel.                      */
  gal_data_t          *chll;  /* Linked list of color channels.        */
};

#endif
