/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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

/* Row selection types. */
enum select_types
{
 /* Different types of row-selection */
 SELECT_TYPE_RANGE,             /* 0 by C standard */
 SELECT_TYPE_INPOLYGON,
 SELECT_TYPE_OUTPOLYGON,
 SELECT_TYPE_EQUAL,
 SELECT_TYPE_NOTEQUAL,

 /* This marks the total number of row-selection criteria. */
 SELECT_TYPE_NUMBER,
};



/* Basic structure. */
struct list_select
{
  gal_data_t          *col;
  int                 type;
  struct list_select *next;
};

struct arithmetic_token
{
  /* First layer elements. */
  int            operator;  /* OPERATOR: Code of operator.                */
  size_t     num_operands;  /* OPERATOR: Number of required operands.     */
  size_t            index;  /* OPERAND: Index in requested columns.       */
  gal_data_t    *constant;  /* OPERAND: a constant/single number.         */
  struct arithmetic_token *next;  /* Pointer to next token.               */
};

struct column_pack
{
  size_t                    start; /* Starting ind. in requested columns. */
  size_t                numsimple; /* Number of simple columns.           */
  struct arithmetic_token *tokens; /* Arithmetic tokens to use.           */
  struct column_pack        *next; /* Next output column.                 */
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
  gal_data_t       *inpolygon;  /* Columns to check if inside polygon.  */
  gal_data_t      *outpolygon;  /* Columns to check if outside polygon. */
  gal_data_t         *polygon;  /* Values of vertices of the polygon.   */
  gal_data_t           *equal;  /* Values to keep in output.            */
  gal_data_t        *notequal;  /* Values to not include in output.     */
  char                  *sort;  /* Column name or number for sorting.   */
  uint8_t          descending;  /* Sort columns in descending order.    */
  size_t                 head;  /* Output only the no. of top rows.     */
  size_t                 tail;  /* Output only the no. of bottom rows.  */
  gal_list_str_t   *catcolumn;  /* Filename to concat column wise.      */
  gal_list_str_t   *catcolhdu;  /* HDU/extension for the catcolumn.     */

  /* Internal. */
  struct column_pack *outcols;  /* Output column packages.              */
  gal_data_t           *table;  /* Linked list of output table columns. */
  struct wcsprm          *wcs;  /* WCS structure for conversion.        */
  int                    nwcs;  /* Number of WCS structures.            */
  gal_data_t      *allcolinfo;  /* Information of all the columns.      */
  gal_data_t         *sortcol;  /* Column to define a sorting.          */
  int               selection;  /* Any row-selection is requested.      */
  gal_data_t          *select;  /* Select rows for output.              */
  struct list_select *selectcol; /* Column to define a range.           */
  uint8_t            freesort;  /* If the sort column should be freed.  */
  uint8_t         *freeselect;  /* If selection columns should be freed.*/
  uint8_t              sortin;  /* If the sort column is in the output. */
  time_t              rawtime;  /* Starting time of the program.        */
  gal_data_t       **colarray;  /* Array of columns, with arithmetic.   */
  size_t          numcolarray;  /* Number of elements in 'colarray'.    */

  /* For arithmetic operators. */
  gal_list_str_t  *wcstoimg_p;  /* Pointer to the node.                 */
  gal_list_str_t  *imgtowcs_p;  /* Pointer to the node.                 */
  size_t             wcstoimg;  /* Output column no, for conversion.    */
  size_t             imgtowcs;  /* Output column no, for conversion.    */
};

#endif
