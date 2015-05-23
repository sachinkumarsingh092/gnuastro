/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include "mesh.h"
#include "fitsarrayvv.h"
#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astnoisechisel" /* Subpackage executable name. */
#define SPACK_NAME      "NoiseChisel"    /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION




/* Macros: */
#define MINACCEPTABLENEAREST 3



struct uiparams
{
  char         *inputname;  /* Name of input file.                 */
  char          *maskname;  /* Name of mask image file.            */
  char              *mhdu;  /* Name of mask image header name.     */
  int         masknameset;
  int   masknameallocated;
  int             mhduset;
  int       numnearestset;
  int      smoothwidthset;
  int       mirrordistset;
  int         minmodeqset;
  int    sigclipmultipset;
  int sigcliptoleranceset;
  int         meshsizeset;
  int             nch1set;
  int             nch2set;
  int     lastmeshfracset;
};





struct noisechiselparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.          */
  struct commonparams      cp;  /* Common parameters.                  */
  struct meshparams        mp;  /* Mesh grid of input image.           */

  /* Input: */
  int                    nwcs;  /* Number of WCS structures.           */
  struct wcsprm          *wcs;  /* Pointer to WCS structures.          */
  int                  bitpix;  /* Input image bitpix value.           */
  size_t             numblank;  /* Number of blank pixels in image.    */

  /* output: */
  int                checkstd;  /* ==1: include the sky STD in checks. */
  char              *meshname;  /* Name of --checkmesh output.         */

  /* Operating mode: */

  /* Internal: */
  time_t              rawtime;  /* Starting time of the program.       */
};

#endif
