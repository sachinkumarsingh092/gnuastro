/*********************************************************************
Warp - Warp images using projective mapping.
Warp is part of GNU Astronomy Utilities (Gnuastro) package.

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
#define PROGRAM_NAME   "Warp"        /* Program full name.       */
#define PROGRAM_EXEC   "astwarp"     /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION






struct optionwarpsll
{
  int                   type;
  double                  v1;
  double                  v2;
  struct optionwarpsll *next;
};





struct warpparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.             */
  struct optionwarpsll *owll;    /* List of 2D rotations.                */
  char         *inputname;  /* Name of input file.                       */
  size_t        hstartwcs;  /* Header keyword No. to start reading WCS.  */
  size_t          hendwcs;  /* Header keyword No. to end reading WCS.    */
  unsigned char keepinputwcs;  /* Wrap the warped/transfomed pixels.     */
  float      maxblankfrac;  /* Maximum fraction of blank pixel in out.   */
  char           *typestr;  /* Type of output image.                     */
  unsigned char     align;  /* Align the image with celestial coord.     */
  char         *rotatestr;  /* String given for rotation.                */
  char          *scalestr;  /* String given for scaling.                 */
  char           *flipstr;  /* String given for flipping.                */
  char          *shearstr;  /* String given for shearing.                */
  char      *translatestr;  /* String given for translation.             */
  char        *projectstr;  /* String given for projection.              */
  char         *matrixstr;  /* String containing transform elements.     */

  /* Internal parameters: */
  gal_data_t       *input;  /* Name of input FITS file.                  */
  int             outtype;  /* Output type.                              */
  double        matrix[9];  /* Warp/Transformation matrix.               */
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
