/*********************************************************************
ImageCrop - Crop a given size from one or multiple images.
ImageCrop is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef MAIN_H
#define MAIN_H

/* Include necessary headers */
#include <gnuastro/data.h>

#include <options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "ImageCrop"     /* Program full name.       */
#define PROGRAM_EXEC   "astimgcrop"    /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION





/* Macros */
#define LOGFILENAME             PROGRAM_EXEC".log"
#define FILENAME_BUFFER_IN_VERB 30



/* Modes of operation. */
enum imgcrop_modes
{
  IMGCROP_MODE_INVALID,         /* For sanity checks.     */

  IMGCROP_MODE_IMG,             /* Use image coordinates. */
  IMGCROP_MODE_WCS,             /* Use WCS coordinates.   */
};




/* The sides of the image keep the celestial coordinates of the four
   sides of this image. With respect to the pixels they are.
*/
struct inputimgs
{
  char             *name;  /* File name of input image.                   */
  size_t            ndim;  /* Number of dimensions of this image.         */
  size_t          *dsize;  /* Size of the image.                          */
  int               nwcs;  /* Number of WCS in each input image.          */
  struct wcsprm     *wcs;  /* WCS structure of each input image.          */
  char           *wcstxt;  /* Text output of each WCS.                    */
  int           nwcskeys;  /* Number of keywords in the header WCS.       */
  double      corners[8];  /* RA and Dec of this image corners (within).  */
  double        sized[2];  /* Width and height of image in degrees.       */
  double  equatorcorr[2];  /* If image crosses the equator, see wcsmode.c.*/
};





/* Main program parameters: */
struct imgcropparams
{
  /* Directly from command-line */
  struct gal_options_common_params cp;  /* Common parameters.             */
  struct gal_linkedlist_stll  *inputs;  /* All input FITS files.          */
  size_t             hstartwcs;  /* Header keyword No. to start read WCS. */
  size_t               hendwcs;  /* Header keyword No. to end read WCS.   */
  uint8_t       zeroisnotblank;  /* ==1: In float or double, keep 0.0.    */
  uint8_t              noblank;  /* ==1: no blank (out of image) pixels.  */
  char                 *suffix;  /* Ending of output file name.           */
  size_t           checkcenter;  /* width of a box to check for zeros     */
  size_t              iwidthin;  /* Image mode width (in pixels).         */
  double                wwidth;  /* WCS mode width (in arcseconds).       */
  double                    ra;  /* RA of one crop box center.            */
  double                   dec;  /* Dec of one crop box center.           */
  double                    xc;  /* Center point, one crop (FITS stnrd).  */
  double                    yc;  /* Center point, one crop (FITS stnrd).  */
  char                *catname;  /* Name of input catalog.                */
  char                 *cathdu;  /* HDU of catalog if its a FITS file.    */
  char                *namecol;  /* Filename (without suffix) of crop col.*/
  char                  *racol;  /* Catalog RA column                     */
  char                 *deccol;  /* Catalog Dec column                    */
  char                   *xcol;  /* Catalog X column                      */
  char                   *ycol;  /* Catalog Y column                      */
  char                *section;  /* Section string.                       */
  char                *polygon;  /* Input string of polygon vertices.     */
  uint8_t           outpolygon;  /* ==1: Keep the inner polygon region.   */
  char                *modestr;  /* ==1: will use X and Y coordiates.     */

  /* Internal */
  int                     mode;  /* Image or WCS mode.                    */
  size_t                 numin;  /* Number of input images.               */
  size_t                numout;  /* Number of output images.              */
  double                   *c1;  /* First coordinate from catalog.        */
  double                   *c2;  /* Second coordinate from catalog.       */
  char                  **name;  /* filename of crop in row.              */
  double             *wpolygon;  /* Array of WCS polygon vertices.        */
  double             *ipolygon;  /* Array of image polygon vertices.      */
  size_t             nvertices;  /* Number of polygon vertices.           */
  long               iwidth[2];  /* Image mode width (in pixels).         */
  double                   res;  /* Resolution in arcseconds              */
  time_t               rawtime;  /* Starting time of the program.         */
  int            outnameisfile;  /* Output filename is a directory.       */
  int                     type;  /* Type of output(s).                    */
  void                 *bitnul;  /* Null value for this data-type.        */
  struct inputimgs       *imgs;  /* WCS and size information for inputs.  */
  gal_data_t              *log;  /* Log file contents.                    */
};

#endif
