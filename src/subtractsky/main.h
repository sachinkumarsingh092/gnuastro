/*********************************************************************
SubtractSky - Find and subtract the sky value from an image.
SubtractSky is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include "fitsarrayvv.h"
#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astsubtractsky" /* Subpackage executable name. */
#define SPACK_NAME      "SubtractSky"    /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION






struct uiparams
{
  char      *inputname;  /* Name of input file.                 */
  int    numnearestset;
  int     mininterpset;
  int   kernelwidthset;
};





struct subtractskyparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.          */
  struct commonparams      cp;  /* Common parameters.                  */

  /* Input: */
  float                  *img;  /* The input image array.              */
  int                    nwcs;  /* Number of WCS structures.           */
  struct wcsprm          *wcs;  /* Pointer to WCS structures.          */
  size_t           numnearest;  /* Number of nearest for interpolation.*/
  size_t            mininterp;  /* Minimum number for interpolation.   */
  size_t          kernelwidth;  /* Width of smoothing kernel.          */

  /* output: */

  /* Operating mode: */

  /* Internal: */
  time_t              rawtime;  /* Starting time of the program.       */
};

#endif
