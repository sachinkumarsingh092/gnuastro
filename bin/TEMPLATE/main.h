/*********************************************************************
TEMPLATE - A minimal set of files and functions to define a program.
TEMPLATE is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Your Name <your@email>
Contributing author(s):
Copyright (C) YYYY, Free Software Foundation, Inc.

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
  struct gal_options_common_params     cp; /* Common parameters.        */
  struct gal_linkedlist_strll *multivalue; /* Pointer to multivalue.    */
  char             *inputname;  /* Input filename.                      */
  uint8_t              *onoff;  /* How to store on/off options.         */

  /* Output: */
  time_t              rawtime;  /* Starting time of the program.        */
};

#endif
