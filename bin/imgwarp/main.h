/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef MAIN_H
#define MAIN_H

#include <pthread.h>

#include <commonparams.h>





/* Progarm name macros: */
#define SPACK           "astimgwarp" /* Subpackage executable name. */
#define SPACK_NAME      "ImageWarp"  /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION
#define LOGFILENAME     SPACK".log"







struct uiparams
{
  char        *inputname;  /* Name of input file.                      */
  char       *matrixname;  /* Name of transform file.                  */
  char     *matrixstring;  /* String containing transform elements.    */

  int              align;  /* ==1: Align the image.                    */
  int           alignset;

  int    matrixstringset;
  int    maxblankfracset;
  int       hstartwcsset;
  int         hendwcsset;
};





struct imgwarpparams
{
  /* Other structures: */
  struct uiparams     up;  /* User interface parameters.                 */
  struct gal_commonparams cp; /* Common parameters.                      */

  /* Input: */
  double           *input;  /* Name of input FITS file.                  */
  double          *matrix;  /* Warp/Transformation matrix.               */
  size_t              is0;  /* Number of rows in input image.            */
  size_t              is1;  /* Number of columns in input image.         */
  size_t              ms0;  /* Matrix number of rows.                    */
  size_t              ms1;  /* Matrix number of columns.                 */
  int         inputbitpix;  /* The type of the input array.              */
  int                nwcs;  /* Number of WCS structures.                 */
  struct wcsprm      *wcs;  /* Pointer to WCS structures.                */
  size_t        hstartwcs;  /* Header keyword No. to start reading WCS.  */
  size_t          hendwcs;  /* Header keyword No. to end reading WCS.    */

  /* Output: */
  size_t           numnul;  /* Number of blank pixels in output.         */
  int          correctwcs;  /* Wrap the warped/transfomed pixels.        */
  int          doubletype;  /* Save output in double not in input type.  */
  int      zerofornoinput;  /* Set the pixels with no input to zero.     */
  float      maxblankfrac;  /* Maximum fraction of blank pixel in out.   */

  /* Internal parameters: */
  double          *output;  /* Warped image array.                       */
  size_t        onaxes[2];  /* Output image size                         */
  size_t        knaxes[2];  /* Output image size                         */
  double         *inverse;  /* Inverse of the input matrix.              */
  time_t          rawtime;  /* Starting time of the program.             */
  size_t       extinds[4];  /* Indexs of the minimum and maximum values. */
  size_t       ordinds[4];  /* Indexs of anticlockwise vertices.         */
  double    outfpixval[2];  /* Pixel value of first output pixel.        */
  double         opixarea;  /* Area of output pix in units of input pix. */
};

#endif
