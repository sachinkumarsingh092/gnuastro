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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef MAIN_H
#define MAIN_H

#include <gnuastro/mesh.h>
#include <gnuastro/fits.h>

#include <commonparams.h>

/* Progarm name macros: */
#define SPACK           "astsubtractsky" /* Subpackage executable name. */
#define SPACK_NAME      "SubtractSky"    /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION





struct uiparams
{
  char          *inputname;  /* Name of input file.                 */
  char           *maskname;  /* Name of mask image file.            */
  char               *mhdu;  /* Name of mask image header name.     */
  char         *kernelname;
  char               *khdu;
  int          masknameset;
  int              mhduset;
  int        kernelnameset;
  int              khduset;
  int        numnearestset;
  int       smoothwidthset;
  int        mirrordistset;
  int          minmodeqset;
  int   fullconvolutionset;
  int fullinterpolationset;
  int        fullsmoothset;
  int     sigclipmultipset;
  int  sigcliptoleranceset;
  int          meshsizeset;
  int              nch1set;
  int              nch2set;
  int      lastmeshfracset;
};





struct subtractskyparams
{
  /* Other structures: */
  struct uiparams         up;  /* User interface parameters.             */
  struct gal_commonparams cp; /* Common parameters.                      */
  struct gal_mesh_params  mp; /* Mesh grid of input image.               */

  /* Input: */
  int               nwcs;  /* Number of WCS structures.                  */
  struct wcsprm     *wcs;  /* Pointer to WCS structures.                 */
  int             bitpix;  /* Input image bitpix value.                  */
  int           anyblank;  /* ==1: thereare blank pixels in input image. */

  /* output: */
  int           checkstd;  /* ==1: include the sky STD in checks.        */
  char         *meshname;  /* Name of --checkmesh output.                */
  char         *convname;  /* Name of --checkconvolution output.         */
  char          *skyname;  /* Name of sky and its STD image.             */

  /* Statistics: */
  float    sigclipmultip; /* Multiple of standard deviation, sigma clip. */
  float sigcliptolerance; /* Tolerance in sigma clip.                    */

  /* Operating mode: */

  /* Internal: */
  float            *conv;  /* Convolved input image.                     */
  time_t         rawtime;  /* Starting time of the program.              */
};

#endif
