/*********************************************************************
Match - A program to match catalogs and WCS warps
Match is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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
#define PROGRAM_NAME   "Match"    /* Program full name.       */
#define PROGRAM_EXEC   "astmatch" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION



enum match_modes
{
  MATCH_MODE_INVALID,           /* ==0 by default. */
  MATCH_MODE_WCS,
  MATCH_MODE_CATALOG,
};




/* Main program parameters structure */
struct matchparams
{
  /* From command-line */
  struct gal_options_common_params     cp; /* Common parameters.        */
  char            *input1name;  /* First input filename.                */
  char            *input2name;  /* Second input filename.               */
  char                  *hdu2;  /* Second input's HDU.                  */
  gal_data_t           *ccol1;  /* Array of first input column names.   */
  gal_data_t           *ccol2;  /* Array of second input column names.  */
  gal_data_t           *coord;  /* Array of manual coordinate values.   */
  gal_data_t         *outcols;  /* Array of second input column names.  */
  gal_data_t        *aperture;  /* Acceptable matching aperture.        */
  uint8_t         logasoutput;  /* Don't rearrange inputs, out is log.  */
  uint8_t          notmatched;  /* Output is rows that don't match.     */

  /* Internal */
  int                    mode;  /* Mode of operation: image or catalog. */
  gal_data_t           *cols1;  /* Column values of first input.        */
  gal_data_t           *cols2;  /* Column values of second input.       */
  gal_list_str_t       *acols;  /* Output columns from first input.     */
  gal_list_str_t       *bcols;  /* Output columns from second input.    */
  size_t                 anum;  /* Number of columns in first input.    */
  size_t                 bnum;  /* Number of columns in second input.   */
  char               *logname;  /* Name of log file.                    */
  char              *out1name;  /* Name of first matched output.        */
  char              *out2name;  /* Name of second matched output.       */
  gal_list_str_t  *stdinlines;  /* Lines given by Standard input.       */

  /* Output: */
  time_t              rawtime;  /* Starting time of the program.        */
};

#endif
