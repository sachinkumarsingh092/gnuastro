/*********************************************************************
tableintern -- Internalfunctions for table input and output.
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
#ifndef __GAL_TABLEINTERN_H__
#define __GAL_TABLEINTERN_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */

#include <gnuastro/fits.h> /* Includes 'gnuastro/data.h' and 'fitsio.h' */
#include <gnuastro/list.h>



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





/************************************************************************/
/***************              Error messages              ***************/
/************************************************************************/
void
gal_tableintern_error_col_selection(char *filename, char *hdu,
                                    char *errorstring);




/************************************************************************/
/***************                 Formats                  ***************/
/************************************************************************/
uint8_t
gal_tableintern_string_to_format(char *string);

char *
gal_tableintern_format_as_string(uint8_t tableformat);

uint8_t
gal_tableintern_string_to_searchin(char *string);

char *
gal_tableintern_searchin_as_string(uint8_t searchin);

void
gal_tableintern_check_fits_format(char *filename, int tableformat);



/************************************************************************/
/***************          Printing information            ***************/
/************************************************************************/
void
gal_tableintern_col_print_info(gal_data_t *col, int tableformat,
                               char *fmt, char *lng);

void
gal_tableintern_read_blank(gal_data_t *col, char *blank);




__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_TABLE_H__ */
