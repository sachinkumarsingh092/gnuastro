/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <gnuastro/fits.h>

#include <options.h>

/* Progarm name macros: */
#define PROGRAM_NAME  "Fits"        /* Program full name.       */
#define PROGRAM_EXEC  "astfits"     /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION





enum fits_mode
  {
    FITS_MODE_INVALID,          /* ==0, by C standard. */

    FITS_MODE_HDU,
    FITS_MODE_KEY,
  };





/* Main program's structure */
struct fitsparams
{
  /* From the environment. */
  struct gal_options_common_params cp;  /* Common parameters.              */
  char                      *filename;  /* Name of input file.             */
  uint8_t                    printall;  /* Print all the header keywords.  */
  uint8_t                        date;  /* Set DATE to current time.       */
  char                       *comment;  /* COMMENT value.                  */
  char                       *history;  /* HISTORY value.                  */
  struct gal_linkedlist_stll    *asis;  /* Strings to write asis.          */
  struct gal_linkedlist_stll  *delete;  /* Keywords to remove.             */
  struct gal_linkedlist_stll  *rename;  /* Rename a keyword.               */
  struct gal_linkedlist_stll  *update;  /* For keywords to update.         */
  struct gal_linkedlist_stll   *write;  /* Full arg. for keywords to add.  */
  uint8_t                 quitonerror;  /* Quit if an error occurs.        */

  /* Internal: */
  int                            mode;
  struct gal_fits_key_ll  *write_keys;  /* Keys to write in the header.    */
  struct gal_fits_key_ll *update_keys;  /* Keys to update in the header.   */
  time_t                      rawtime;  /* Starting time of the program.   */
};

#endif
