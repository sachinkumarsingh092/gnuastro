/*********************************************************************
CosmicCalculator - Calculate cosmological parameters
CosmicCalculator is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "CosmicCalculator" /* Program full name.       */
#define PROGRAM_EXEC   "astcosmiccal"     /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION







/* Main program parameters structure */
struct cosmiccalparams
{
  /* Other structures: */
  struct gal_options_common_params cp;  /* Common parameters.           */

  /* Input: */
  double              redshift; /* Redshift of interest.                */
  double                    H0; /* Current expansion rate (km/sec/Mpc). */
  double               olambda; /* Current cosmological constant dens.  */
  double               omatter; /* Current matter density.              */
  double            oradiation; /* Current radiation density.           */
  double            solidangle; /* Solid angle for volume (in stradian).*/

  /* Output: */
  uint8_t           onlyvolume; /* Only print the volume in Mpc^3.      */
  uint8_t       onlyabsmagconv; /* Only print conversion to abs. mag.   */

  /* Internal: */
  time_t               rawtime; /* Starting time of the program.        */
};

#endif
