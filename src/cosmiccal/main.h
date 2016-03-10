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

#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astcosmiccal" /* Subpackage executable name. */
#define SPACK_NAME      "CosmicCalculator"  /* Subpackage full name.  */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION






struct uiparams
{
  int             redshiftset;
  int            curvatureset;
  int                   H0set;
  int              olambdaset;
  int              omatterset;
  int           oradiationset;
};





struct cosmiccalparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.         */
  struct commonparams      cp;  /* Common parameters.                 */

  /* Input: */
  double             redshift;  /* Redshift of interest.              */
  double            curvature;  /* Curvature of the universe.         */
  double                   H0;  /* Current expansion rate.            */
  double              olambda;  /* Current cosmological constant dens.*/
  double              omatter;  /* Current matter density.            */
  double           oradiation;  /* Current radiation density.         */

  /* Output: */

  /* Operating mode: */

  /* Internal: */
  double                ocurv;  /* Curvature density today.           */

  time_t              rawtime;  /* Starting time of the program.      */
};

#endif
