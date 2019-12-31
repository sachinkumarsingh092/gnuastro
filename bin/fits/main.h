/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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

#include <gnuastro/list.h>
#include <gnuastro/fits.h>

#include <gnuastro-internal/options.h>

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
  struct gal_options_common_params cp;  /* Common parameters.           */
  int  hdu_in_commandline;     /* HDU wasn't given in config. file.     */
  char          *filename;     /* Name of input file.                   */
  char            *outhdu;     /* HDU of output (only when necessary).  */
  gal_list_str_t  *remove;     /* Remove extensions from a file.        */
  gal_list_str_t    *copy;     /* Copy extensions to output.            */
  gal_list_str_t     *cut;     /* Copy ext. to output and remove.       */
  uint8_t         numhdus;     /* Print number of HDUs in FITS file.    */
  uint8_t         datasum;     /* Calculate and print HDU's datasum.    */
  uint8_t   primaryimghdu;     /* Copy/cut HDU into primary HDU.        */
  uint8_t    printallkeys;     /* Print all the header keywords.        */
  uint8_t            date;     /* Set DATE to current time.             */
  gal_list_str_t    *asis;     /* Strings to write asis.                */
  gal_list_str_t  *delete;     /* Keywords to remove.                   */
  gal_list_str_t  *rename;     /* Rename a keyword.                     */
  gal_list_str_t  *update;     /* For keywords to update.               */
  gal_list_str_t   *write;     /* Full arg. for keywords to add.        */
  gal_list_str_t *history;     /* HISTORY value.                        */
  gal_list_str_t *comment;     /* COMMENT value.                        */
  uint8_t         *verify;     /* Verify the CHECKSUM and DATASUM keys. */
  char          *copykeys;     /* Range of keywords to copy in output.  */
  char         *datetosec;     /* Convert FITS date to seconds.         */
  uint8_t     quitonerror;     /* Quit if an error occurs.              */

  /* Internal: */
  int                         mode;  /* Operating on HDUs or keywords.  */
  long            copykeysrange[2];  /* Start and end of copy.          */
  gal_fits_list_key_t  *write_keys;  /* Keys to write in the header.    */
  gal_fits_list_key_t *update_keys;  /* Keys to update in the header.   */
  time_t                   rawtime;  /* Starting time of the program.   */
};

#endif
