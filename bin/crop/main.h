/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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
#define PROGRAM_NAME   "Crop"     /* Program full name.       */
#define PROGRAM_EXEC   "astcrop"    /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION





/* Macros */
#define LOGFILENAME             PROGRAM_EXEC".log"
#define FILENAME_BUFFER_IN_VERB 30
#define MAXDIM                  3


/* Modes to interpret coordinates. */
enum crop_modes
{
  IMGCROP_MODE_INVALID,         /* For sanity checks.     */

  IMGCROP_MODE_IMG,             /* Use image coordinates. */
  IMGCROP_MODE_WCS,             /* Use WCS coordinates.   */
};




/* The sides of the image keep the celestial coordinates of the four
   sides of this image. With respect to the pixels they are. */
struct inputimgs
{
  char             *name;  /* File name of input image.                   */
  size_t            ndim;  /* Number of dimensions of this image.         */
  size_t          *dsize;  /* Size of the image.                          */
  int               nwcs;  /* Number of WCS in each input image.          */
  struct wcsprm     *wcs;  /* WCS structure of each input image.          */
  char           *wcstxt;  /* Text output of each WCS.                    */
  int           nwcskeys;  /* Number of keywords in the header WCS.       */
  double     corners[24];  /* WCS of corners (24: for 3D, 8: for 2D).     */
  double   sized[MAXDIM];  /* Width and height of image in degrees.       */
  double  equatorcorr[2];  /* If image crosses the equator, see wcsmode.c.*/
};





/* Main program parameters: */
struct cropparams
{
  /* Directly from command-line */
  struct gal_options_common_params cp;  /* Common parameters.             */
  gal_list_str_t       *inputs;  /* All input FITS files.                 */
  size_t             hstartwcs;  /* Header keyword No. to start read WCS. */
  size_t               hendwcs;  /* Header keyword No. to end read WCS.   */
  int                     mode;  /* Image or WCS mode.                    */
  uint8_t       zeroisnotblank;  /* ==1: In float or double, keep 0.0.    */
  uint8_t              noblank;  /* ==1: no blank (out of image) pixels.  */
  char                 *suffix;  /* Ending of output file name.           */
  gal_data_t    *incheckcenter;  /* Value given to '--checkcenter'.       */
  gal_data_t           *center;  /* Center position of crop.              */
  gal_data_t            *width;  /* Width of crop when defined by center. */
  char                *catname;  /* Name of input catalog.                */
  char                 *cathdu;  /* HDU of catalog if its a FITS file.    */
  char                *namecol;  /* Filename (without suffix) of crop col.*/
  gal_list_str_t     *coordcol;  /* Column in cat containing coordinates. */
  char                *section;  /* Section string.                       */
  gal_data_t          *polygon;  /* Input string of polygon vertices.     */
  uint8_t           polygonout;  /* ==1: Keep the inner polygon region.   */
  uint8_t          polygonsort;  /* Don't sort polygon vertices.          */

  /* Internal */
  size_t                 numin;  /* Number of input images.               */
  size_t                numout;  /* Number of output images.              */
  double        **centercoords;  /* A 1D array for the center position.   */
  size_t           checkcenter;  /* width of a box to check for zeros     */
  char                  **name;  /* filename of crop in row.              */
  double             *wpolygon;  /* Array of WCS polygon vertices.        */
  double             *ipolygon;  /* Array of image polygon vertices.      */
  size_t             nvertices;  /* Number of polygon vertices.           */
  long          iwidth[MAXDIM];  /* Image mode width (in pixels).         */
  double             *pixscale;  /* Raw resolution in each dimension.     */
  time_t               rawtime;  /* Starting time of the program.         */
  int            outnameisfile;  /* Output filename is a directory.       */
  int                     type;  /* Type of output(s).                    */
  void           *blankptrread;  /* Null value for reading of output type.*/
  void          *blankptrwrite;  /* Null value for writing of output type.*/
  struct inputimgs       *imgs;  /* WCS and size information for inputs.  */
  gal_data_t              *log;  /* Log file contents.                    */
};

#endif
