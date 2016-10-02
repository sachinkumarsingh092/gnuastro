/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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

#include <commonparams.h>

/* Progarm name macros: */
#define SPACK           "asttable" /* Subpackage executable name. */
#define SPACK_NAME      "Table"    /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION

#define MAX_COL_FORMAT_LENGTH 20



/* Structure to keep the information for each output column */
struct outcolumn
{
  size_t inindex;               /* Row index (from 0) in input array.   */
  int   datatype;               /* Type of data (from CFITSIO macros).  */
  int     anynul;               /* If there is any blank characters.    */
  void   *nulval;               /* The blank value for this column.     */
  void     *data;               /* Array keeping the column data.       */
  char fmt[MAX_COL_FORMAT_LENGTH];  /* format to use in printf.         */
};





/* User interface structure. */
struct uiparams
{
  int             information;  /* ==1, only print FITS information.    */
  char              *fitsname;  /* Name of input FITS file.             */
  char               *txtname;  /* Name of input text file.             */
  int              ignorecase;  /* Ignore case matching column names.   */

  /* Input table parameters. */
  size_t                ncols;  /* Number of columns in table.          */
  int               *datatype;  /* Type of data in column.              */
  char                **ttstr;  /* TFORM (another format for type).     */
  char                **tname;  /* Column name (one word).              */
  char                **tunit;  /* Unit of values in column.            */
  double            *txtarray;  /* Array keeping text file values.      */
  int         infitstabletype;  /* Input table is ASCII or binary.      */

  /* Print parameters: */
  int                     feg;  /* format of floating points.           */
  size_t            sintwidth;  /* Full width for short integers.       */
  size_t            lintwidth;  /* Full width for short integers.       */
  size_t           floatwidth;  /* Full width for all floats.           */
  size_t          doublewidth;  /* Full width for all doubles.          */
  size_t             strwidth;  /* Full width for all floats.           */
  size_t       floatprecision;  /* Number of decimals for floats.       */
  size_t      doubleprecision;  /* Number of decimals for doubles.      */

  /* If values are set: */
  int                inputset;
  int          informationset;
  int           ignorecaseset;
  int                  fegset;
  int            sintwidthset;
  int            lintwidthset;
  int           floatwidthset;
  int          doublewidthset;
  int             strwidthset;
  int       floatprecisionset;
  int      doubleprecisionset;
  int        fitstabletypeset;


  struct gal_linkedlist_stll *columns;
};





/* Main program parameters structure */
struct tableparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.           */
  struct gal_commonparams  cp;  /* Common parameters.                   */

  /* Input: */
  fitsfile           *fitsptr;  /* FITS pointer (input or output).      */

  /* Output: */
  size_t                nrows;  /* Number of rows in table.             */
  size_t               nocols;  /* Number of output columns.            */
  struct outcolumn     *ocols;  /* Array of output column informatio.   */
  int            outputtofits;  /* ==1: output is a FITS file.          */
  int             outputtotxt;  /* ==1: output is a text file.          */
  int          outputtostdout;  /* ==1: output is the standard output.  */
  int           fitstabletype;  /* ASCII, or binary table CFITSIO macro.*/

  /* Internal: */
  int                onlyview;
  time_t              rawtime;  /* Starting time of the program.        */
};

#endif
