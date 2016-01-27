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
#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astimgarith"       /* Subpackage executable name. */
#define SPACK_NAME      "ImageArithmetic"  /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION


/* Constants: */
#define MAXNUMIMAGES       10


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
  float             *array;  /* Main array to keep results.             */
  float               *tmp;  /* Secondary array for temporary reading.  */
  char          *firstname;  /* First image name.                       */
  size_t                s0;  /* Length of image along first C axis.     */
  size_t                s1;  /* Length of image along second C axis.    */

  /* Output: */


  /* Operating mode: */

  /* Internal: */
  time_t           rawtime;  /* Starting time of the program.           */
};

#endif
