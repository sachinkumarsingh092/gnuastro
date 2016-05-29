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

#include "fitsarrayvv.h"
#include "commonparams.h"

/* Progarm name macros: */
#define SPACK           "astheader" /* Subpackage executable name. */
#define SPACK_NAME      "Header"    /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION






struct uiparams
{
  char                     *inputname;  /* Name of input file.             */
  struct gal_linkedlist_stll  *rename;  /* Rename a keyword.               */
  struct gal_linkedlist_stll  *update;  /* For keywords to update.         */
  struct gal_linkedlist_stll   *write;  /* Full arg. for keywords to add.  */
};





struct headerparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.         */
  struct gal_commonparams  cp;  /* Common parameters.                 */

  /* Input: */
  int                    nwcs;  /* Number of WCS structures.          */
  fitsfile              *fptr;  /* FITS file pointer.                 */
  struct wcsprm          *wcs;  /* Pointer to WCS structures.         */

  /* Output: */
  int                    date;  /* Set DATE to current time.          */
  char               *comment;  /* COMMENT value.                     */
  char               *history;  /* HISTORY value.                     */
  struct gal_linkedlist_stll *asis; /* Strings to write as is into header.*/
  struct gal_linkedlist_stll *delete; /* Keywords to remove.          */
  struct gal_linkedlist_stll *renamefrom; /* The initial value of the keyword. */
  struct gal_linkedlist_stll *renameto; /* The final value of the keyword. */
  struct gal_fitsarray_header_ll *update; /* Linked list of keywords to update. */
  struct gal_fitsarray_header_ll *write; /* Linked list of keywords to add. */

  /* Operating mode: */
  int             quitonerror;  /* Quit if an error occurs.           */

  /* Internal: */
  int                onlyview;
  time_t              rawtime;  /* Starting time of the program.      */
};

#endif
