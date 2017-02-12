/*********************************************************************
Header - View and manipulate a data file header
Header is part of GNU Astronomy Utilities (Gnuastro) package.

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
#define PROGRAM_NAME  "Header"        /* Program full name.       */
#define PROGRAM_EXEC  "astheader"     /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION







struct headerparams
{
  /* Common parameters */
  struct gal_options_common_params  cp;  /* Common parameters.        */

  /* Input: */
  int                            nwcs;  /* Number of WCS structures.       */
  fitsfile                      *fptr;  /* FITS file pointer.              */
  struct wcsprm                  *wcs;  /* Pointer to WCS structures.      */

  /* Output: */
  uint8_t                        date;  /* Set DATE to current time.       */
  char                       *comment;  /* COMMENT value.                  */
  char                       *history;  /* HISTORY value.                  */
  struct gal_linkedlist_stll    *asis;  /* Strings to write asis.          */
  struct gal_linkedlist_stll  *delete;  /* Keywords to remove.             */
  struct gal_linkedlist_stll *renamefrom; /* Initial value of the keyword. */
  struct gal_linkedlist_stll *renameto; /* The final value of the keyword. */
  struct gal_fits_key_ll      *update;  /* Keywords to update. */
  struct gal_fits_key_ll       *write;  /* Keywords to add.                */

  /* Operating mode: */
  uint8_t                 quitonerror;  /* Quit if an error occurs.        */

  /* Internal: */
  int                        onlyview;
  time_t                      rawtime;  /* Starting time of the program.   */
  char                      *filename;  /* Name of input file.             */
  struct gal_linkedlist_stll  *rename;  /* Rename a keyword.               */
  struct gal_linkedlist_stll *updatestr;/* For keywords to update.         */
  struct gal_linkedlist_stll *writestr; /* Full arg. for keywords to add.  */
};

#endif
