/*********************************************************************
MakeNoise - Add noise to a dataset.
MakeNoise is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <gsl/gsl_rng.h>

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "MakeNoise"     /* Program full name.       */
#define PROGRAM_EXEC   "astmknoise"    /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION





/* Main program parameters structure */
struct mknoiseparams
{
  /* From command-line */
  struct gal_options_common_params cp;   /* Common parameters.           */
  char        *inputname;    /* Input filename.                          */
  double           sigma;    /* Total noise sigma (ignoring others).     */
  double    instrumental;    /* Standard deviation constants.            */
  double       zeropoint;    /* Zeropoint magnitude of image.            */
  double  background_mag;    /* Background in magnitudes.                */
  uint8_t        envseed;    /* ==1, generate a random seed.             */

  /* Internal */
  gal_data_t      *input;    /* Input image data in double precision.    */
  double      background;    /* Background in units of brightness.       */
  gsl_rng           *rng;    /* Main instance of random number generator.*/
  const char   *rng_name;    /* The type/name of the Random number gen.  */
  unsigned long rng_seed;    /* Seed of Random number generator.         */
  time_t         rawtime;    /* Starting time of the program.            */
};

#endif
