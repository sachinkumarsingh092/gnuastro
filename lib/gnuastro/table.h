/*********************************************************************
table -- functions for table input and output.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef __GAL_TABLE_H__
#define __GAL_TABLE_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */

#include <gnuastro/fits.h> /* Includes `gnuastro/data.h' and `fitsio.h' */
#include <gnuastro/linkedlist.h>



/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */





/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */





/* `printf' default formattings. Note that the string type shouldn't have
   any precision and for the width,  */
#define GAL_TABLE_DEF_STR_WIDTH       6
#define GAL_TABLE_DEF_INT_WIDTH       6
#define GAL_TABLE_DEF_LINT_WIDTH      10
#define GAL_TABLE_DEF_FLT_WIDTH       12
#define GAL_TABLE_DEF_DBL_WIDTH       18

#define GAL_TABLE_DEF_INT_PRECISION   0
#define GAL_TABLE_DEF_FLT_PRECISION   6
#define GAL_TABLE_DEF_DBL_PRECISION   14



/* Formats of table storage for input or output, as strings and
   integers. */
#define GAL_TABLE_STRING_FORMAT_TXT   "txt"
#define GAL_TABLE_STRING_FORMAT_AFITS "fits-ascii"
#define GAL_TABLE_STRING_FORMAT_BFITS "fits-binary"
enum gal_table_types
{
  GAL_TABLE_FORMAT_INVALID,       /* Invalid (=0 by C standard).       */

  GAL_TABLE_FORMAT_TXT,           /* Plain text table.                 */
  GAL_TABLE_FORMAT_AFITS,         /* FITS ASCII table.                 */
  GAL_TABLE_FORMAT_BFITS,         /* FITS binary table.                */
};





/* When the desired column is not a number, should the string match or
   regular expression search be in the name, units or comments of the
   columns? */
#define GAL_TABLE_STRING_SEARCH_NAME    "name"
#define GAL_TABLE_STRING_SEARCH_UNIT    "unit"
#define GAL_TABLE_STRING_SEARCH_COMMENT "comment"
enum gal_table_where_to_search
{
  GAL_TABLE_SEARCH_INVALID,       /* Invalid (=0 by C standard).     */

  GAL_TABLE_SEARCH_NAME,          /* Match/search in names.          */
  GAL_TABLE_SEARCH_UNIT,          /* Match/search in units.          */
  GAL_TABLE_SEARCH_COMMENT,       /* Match/search in comments.       */
};





/* Particular display formats for different types: Integers or floating
   point numbers can be printed in various display formats. The values here
   are stored in `gal_data_t' when necessary to help in printing of the
   data.*/
enum gal_table_diplay_formats
{
  GAL_TABLE_DISPLAY_FMT_INVALID,        /* Invalid (=0 by C standard). */

  GAL_TABLE_DISPLAY_FMT_STRING,         /* String/Character.           */
  GAL_TABLE_DISPLAY_FMT_DECIMAL,        /* Integers: signed decimal.   */
  GAL_TABLE_DISPLAY_FMT_UDECIMAL,       /* Integers: unsigned decimal. */
  GAL_TABLE_DISPLAY_FMT_OCTAL,          /* Integers: octal.            */
  GAL_TABLE_DISPLAY_FMT_HEX,            /* Integers: hexadecimal.      */
  GAL_TABLE_DISPLAY_FMT_FLOAT,          /* Floats: with decimal point. */
  GAL_TABLE_DISPLAY_FMT_EXP,            /* Floats: as exponential.     */
  GAL_TABLE_DISPLAY_FMT_GENERAL,        /* Floats: general (%g in C).  */
};





/* Functions */
int
gal_table_string_to_format(char *string);

void
gal_table_check_fits_format(char *filename, int tableformat);

int
gal_table_string_to_searchin(char *string);

void
gal_table_too_many_columns(char *filename);

void
gal_table_print_info(gal_data_t *allcols, size_t numcols, size_t numrows);

void
gal_table_col_print_info(gal_data_t *col, int tableformat,
                         char *fmt, char *lng);

void
gal_table_read_blank(gal_data_t *col, char *blank);

gal_data_t *
gal_table_info(char *filename, char *hdu, size_t *numcols,
               size_t *numrows, int *tableformat);

gal_data_t *
gal_table_read(char *filename, char *hdu, struct gal_linkedlist_stll *cols,
               int searchin, int ignorecase, int minmapsize);

void
gal_table_write(gal_data_t *cols, char *comments, int tableformat,
                char *filename, int dontdelete);





__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_TABLE_H__ */
