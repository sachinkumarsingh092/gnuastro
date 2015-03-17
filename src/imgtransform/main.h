/*********************************************************************
ImageTransform - Transform images (e.g., rotate, scale, sheer and ...)
ImageTransform is part of GNU Astronomy Utilities (gnuastro) package.

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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef MAIN_H
#define MAIN_H

#include <pthread.h>

#include "commonparams.h"


/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astimgtransform" /* Subpackage executable name. */
#define SPACK_NAME      "ImageTransform"  /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION
#define LOGFILENAME     SPACK".log"







struct uiparams
{
  char        *inputname;  /* Name of input file.                      */
  char    *transformname;  /* Name of transform file.                  */
  char  *transformstring;  /* String containing transform elements.    */

  int transformstringset;
};





struct imgtransformparams
{
  /* Other structures: */
  struct uiparams     up;  /* User interface parameters.               */
  struct commonparams cp;  /* Common parameters.                       */

  /* Input: */
  float            *input;  /* Name of input FITS file.                */
  double       *transform;  /* Name of input transformation.           */
  size_t              ts0;  /* Transform number of rows.               */
  size_t              ts1;  /* Transform number of columns.            */

  /* Internal parameters: */
  time_t         rawtime;  /* Starting time of the program.            */
};

#endif
