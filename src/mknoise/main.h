/*********************************************************************
MakeNoise - Add noise to a dataset.
MakeNoise is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <gsl/gsl_rng.h>

#include "fitsarrayvv.h"
#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astmknoise" /* Subpackage executable name. */
#define SPACK_NAME      "MakeNoise"     /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION






struct uiparams
{
  char       *inputname;   /* Name of input file.                      */

  int     backgroundset;
  int      zeropointset;
  int         stdaddset;
  int       randseedset;
};





struct mknoiseparams
{
  /* Other structures: */
  struct uiparams     up;  /* User interface parameters.               */
  struct gal_commonparams cp; /* Common parameters.                    */

  /* Input: */
  int            envseed;  /* ==1, generate a random seed.             */
  double          *input;  /* Input image data in double precision.    */
  int        inputbitpix;  /* Input BITPIX header keyword value.       */
  size_t             is0;  /* The number of rows in the input image.   */
  size_t             is1;  /* The number of columns in the input image.*/
  int           anyblank;  /* ==1: There are blank pixels in input.    */
  int               nwcs;  /* Number of WCS structures.                */
  struct wcsprm     *wcs;  /* Pointer to WCS structures.               */
  double      background;  /* Mean of noise probability distribution.  */
  double     mbackground;  /* Background in magnitudes.                */
  double          stdadd;  /* Standard deviation constants.            */
  double       zeropoint;  /* Zeropoint magnitude of image.            */

  /* Random number generator */
  gsl_rng           *rng;  /* Main instance of random number generator.*/

  /* Output: */
  int         doubletype;  /* Save the output in double type.          */

  /* Internal: */
  time_t          rawtime; /* Starting time of the program.            */
  char rng_type[FLEN_VALUE];   /* The type of the Random number gen.  */
  long           rng_seed; /* Seed of Random number generator.         */
};

#endif
