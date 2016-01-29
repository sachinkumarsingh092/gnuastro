/*********************************************************************
ImageArithmetic - Do arithmetic operations on images.
ImageArithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include "linkedlist.h"
#include "fitsarrayvv.h"
#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astarithmetic"   /* Subpackage executable name. */
#define SPACK_NAME      "Arithmetic"      /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION


/* Constants: */
#define NEGDASHREPLACE 11  /* A vertical tab (ASCII=11) for negative dash */
#define NOOPTARRAY     NULL
#define NOOPTNUMBER    NAN
#define MAXNUMIMAGES   10
#define NOOPTFILENAME  ""


/* The operands can be in three ways:

      1. A file name which is not read yet. When inactive, NOOPTFILENAME.
      2. A number which is not used yet. When inactive, NOOPTNUMBER.
      3. An array (operator output). When inactive, NOOPTARRAY.

   In every node of the operand linked list, only one of these should
   be active. The other two should be inactive. Otherwise it will be a
   bug and will cause problems. All the operands operate on this
   premise.
*/
struct operand
{
  char  *filename;
  char       *hdu;
  double   number;
  double   *array;
  struct operand *next;
};


struct uiparams
{
  char        *maskname;    /* Name of mask image.                     */
  char            *mhdu;    /* Mask image HDU.                         */
  int       masknameset;
  int masknameallocated;
  int           mhduset;
  char *hdus[MAXNUMIMAGES];  /* Array of pointers to HDU strings.      */
};


struct imgarithparams
{
  /* Other structures: */
  struct uiparams       up;  /* User interface parameters.              */
  struct commonparams   cp;  /* Common parameters.                      */

  /* Input: */
  struct stll      *tokens;  /* Tokens to do arithmetic.                */
  char          *firstname;  /* Name of first input for output name.    */
  float             *array;  /* Main array to keep results.             */
  float               *tmp;  /* Secondary array for temporary reading.  */
  size_t        addcounter;  /* The number of FITS images added.        */
  size_t        popcounter;  /* The number of FITS images popped.       */
  size_t                s0;  /* Length of image along first C axis.     */
  size_t                s1;  /* Length of image along second C axis.    */
  int              obitpix;  /* The type of the output image.           */
  int                 nwcs;  /* The number of WCS coordinates.          */
  struct wcsprm       *wcs;  /* The WCS structure.                      */
  int             anyblank;  /* If there are blank pixels in the image. */

  /* Output: */


  /* Operating mode: */

  /* Internal: */
  struct operand *operands;  /* The operands linked list.               */
  time_t           rawtime;  /* Starting time of the program.           */
};

#endif
