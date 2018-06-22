/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2018, Free Software Foundation, Inc.

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
#define PROGRAM_NAME   "Table"         /* Program full name.       */
#define PROGRAM_EXEC   "asttable"      /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION







/* Main program parameters structure */
struct tableparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.            */
  char              *filename;  /* Input filename.                      */
  gal_list_str_t     *columns;  /* List of given columns.               */
  uint8_t         information;  /* ==1: only print FITS information.    */
  uint8_t     colinfoinstdout;  /* ==1: print column metadata in CL.    */

  /* Output: */
  gal_data_t           *table;  /* Linked list of output table columns. */
  gal_data_t      *allcolinfo;  /* Information of all the columns.      */
  time_t              rawtime;  /* Starting time of the program.        */
};

#endif
