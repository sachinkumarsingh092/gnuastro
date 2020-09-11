/*********************************************************************
Query - Retreive data from a remote data server.
Query is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2020, Free Software Foundation, Inc.

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
#define PROGRAM_NAME   "query"    /* Program full name.       */
#define PROGRAM_EXEC   "astquery" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION







/* Main program parameters structure */
struct queryparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.           */
  int                 database;  /* ID of database to use.             */
  char             *datasetstr;  /* ID of dataset in database to use.  */
  gal_data_t           *center;  /* Center position of query.          */
  double                radius;  /* Radius around center.              */
  char                  *query;  /* Raw query string.                  */
  gal_list_str_t      *columns;  /* Columns to extract from database.  */

  /* Output: */
  time_t               rawtime;  /* Starting time of the program.      */
};

#endif
