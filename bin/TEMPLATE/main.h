/*********************************************************************
TEMPLATE - A minimal set of files and functions to define a program.
TEMPLATE is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2019, Free Software Foundation, Inc.

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
#define PROGRAM_NAME   "TEMPLATE"    /* Program full name.       */
#define PROGRAM_EXEC   "astTEMPLATE" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION







/* Main program parameters structure */
struct TEMPLATEparams
{
  /* From command-line */
  struct gal_options_common_params     cp; /* Common parameters.           */
  char             *inputname;  /* Input filename.                         */
  gal_list_str_t  *multivalue;  /* List of values given to "multivalue"    */
  uint8_t               onoff;  /* ==1 if option is called, ==0 otherwise. */

  /* Output: */
  time_t              rawtime;  /* Starting time of the program.           */
};

#endif
