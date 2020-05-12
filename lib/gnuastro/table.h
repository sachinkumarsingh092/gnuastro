/*********************************************************************
table -- functions for table input and output.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef __GAL_TABLE_H__
#define __GAL_TABLE_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */

#include <gnuastro/list.h>
#include <gnuastro/fits.h>



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





/* 'printf' default formattings. Note that the string type shouldn't have
   any precision and for the width,  */
#define GAL_TABLE_DEF_WIDTH_STR       6
#define GAL_TABLE_DEF_WIDTH_INT       6
#define GAL_TABLE_DEF_WIDTH_LINT      10
#define GAL_TABLE_DEF_WIDTH_FLT       13
#define GAL_TABLE_DEF_WIDTH_DBL       18

#define GAL_TABLE_DEF_PRECISION_INT   0
#define GAL_TABLE_DEF_PRECISION_FLT   6
#define GAL_TABLE_DEF_PRECISION_DBL   14





/* Particular display formats for different types: Integers or floating
   point numbers can be printed in various display formats. The values here
   are stored in 'gal_data_t' when necessary to help in printing of the
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





/* Formats of table storage for input or output, as strings and
   integers. */
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
enum gal_table_where_to_search
{
  GAL_TABLE_SEARCH_INVALID,       /* Invalid (=0 by C standard).     */

  GAL_TABLE_SEARCH_NAME,          /* Match/search in names.          */
  GAL_TABLE_SEARCH_UNIT,          /* Match/search in units.          */
  GAL_TABLE_SEARCH_COMMENT,       /* Match/search in comments.       */
};



/************************************************************************/
/***************         Information about a table        ***************/
/************************************************************************/
gal_data_t *
gal_table_info(char *filename, char *hdu, gal_list_str_t *lines,
               size_t *numcols, size_t *numrows, int *tableformat);

void
gal_table_print_info(gal_data_t *allcols, size_t numcols, size_t numrows);



/************************************************************************/
/***************               Read a table               ***************/
/************************************************************************/
gal_data_t *
gal_table_read(char *filename, char *hdu, gal_list_str_t *lines,
               gal_list_str_t *cols, int searchin, int ignorecase,
               size_t minmapsize, int quietmmap, size_t *colmatch);

gal_list_sizet_t *
gal_table_list_of_indexs(gal_list_str_t *cols, gal_data_t *allcols,
                         size_t numcols, int searchin, int ignorecase,
                         char *filename, char *hdu, size_t *colmatch);



/************************************************************************/
/***************              Write a table               ***************/
/************************************************************************/
void
gal_table_comments_add_intro(gal_list_str_t **comments,
                             char *program_string, time_t *rawtime);

void
gal_table_write(gal_data_t *cols, gal_list_str_t *comments,
                int tableformat, char *filename, char *extname,
                uint8_t colinfoinstdout);

void
gal_table_write_log(gal_data_t *logll, char *program_string,
                    time_t *rawtime, gal_list_str_t *comments,
                    char *filename, int quiet);





__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_TABLE_H__ */
