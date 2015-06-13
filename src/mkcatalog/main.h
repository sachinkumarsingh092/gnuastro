/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

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
#define SPACK           "astmkcatalog" /* Subpackage executable name. */
#define SPACK_NAME      "MakeCatalog"  /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION






struct uiparams
{
  char             *inputname;  /* Name of input file.               */
  char              *maskname;  /* Name of masked file.              */
  char                  *mhdu;  /* HDU of mask image.                */
  char           *objlabsname;  /* Name of object labels file.       */
  char                *objhdu;  /* HDU of object labels image.       */
  char         *clumplabsname;  /* Name of clump labels file .       */
  char              *clumphdu;  /* HDU of clump labels image.        */
  char               *skyname;  /* Sky value image file name.        */
  char                *skyhdu;  /* Sky HDU name.                     */
  char               *stdname;  /* Sky STD value image file name.    */
  char                *stdhdu;  /* Sky STD HDU name.                 */

  int             masknameset;
  int                 mhduset;
  int          objlabsnameset;
  int               objhduset;
  int        clumplabsnameset;
  int             clumphduset;
  int              skynameset;
  int               skyhduset;
  int              stdnameset;
  int               stdhduset;
};





struct mkcatalogparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.         */
  struct commonparams      cp;  /* Common parameters.                 */

  /* Input: */
  float                  *img;  /* Input image.                       */
  long               *objects;  /* Object labels on each pixel.       */
  long                *clumps;  /* Clump labels on each pixel.        */
  float                  *sky;  /* Sky value on each pixel.           */
  float                  *std;  /* Sky STD value on each pixel.       */
  int                    nwcs;  /* Number of WCS structures.          */
  struct wcsprm          *wcs;  /* Pointer to WCS structures.         */
  size_t                   s0;  /* Size of input (first C axis).      */
  size_t                   s1;  /* Size of input (second C axis).     */

  /* Output: */

  /* Operating mode: */

  /* Internal: */
  time_t              rawtime;  /* Starting time of the program.      */
};

#endif
