/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
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
#define PROGRAM_NAME   "Table"         /* Program full name.       */
#define PROGRAM_EXEC   "asttable"      /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION





struct list_range
{
  gal_data_t           *v;
  struct list_range *next;
};





/* Main program parameters structure */
struct tableparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.            */
  char              *filename;  /* Input filename.                      */
  char               *wcsfile;  /* File with WCS.                       */
  char                *wcshdu;  /* HDU in file with WCS.                */
  gal_list_str_t     *columns;  /* List of given columns.               */
  uint8_t         information;  /* ==1: only print FITS information.    */
  uint8_t     colinfoinstdout;  /* ==1: print column metadata in CL.    */
  gal_data_t           *range;  /* Range to limit output.               */
  char                  *sort;  /* Column name or number for sorting.   */
  uint8_t          descending;  /* Sort columns in descending order.    */
  size_t                 head;  /* Output only the no. of top rows.     */
  size_t                 tail;  /* Output only the no. of bottom rows.  */

  /* Internal. */
  gal_data_t           *table;  /* Linked list of output table columns. */
  struct wcsprm          *wcs;  /* WCS structure for conversion.        */
  int                    nwcs;  /* Number of WCS structures.            */
  gal_data_t      *allcolinfo;  /* Information of all the columns.      */
  gal_data_t         *sortcol;  /* Column to define a sorting.          */
  struct list_range *rangecol;  /* Column to define a range.            */
  uint8_t            freesort;  /* If the sort column should be freed.  */
  uint8_t          *freerange;  /* If the range column should be freed. */
  uint8_t              sortin;  /* If the sort column is in the output. */
  time_t              rawtime;  /* Starting time of the program.        */

  /* For arithmetic operators. */
  gal_list_str_t  *wcstoimg_p;  /* Pointer to the node.                 */
  gal_list_str_t  *imgtowcs_p;  /* Pointer to the node.                 */
  size_t             wcstoimg;  /* Output column no, for conversion.    */
  size_t             imgtowcs;  /* Output column no, for conversion.    */
};

#endif
