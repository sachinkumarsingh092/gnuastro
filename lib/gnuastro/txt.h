/*********************************************************************
txt -- functions to deal with plain text files.
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
#ifndef __GAL_TXT_H__
#define __GAL_TXT_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */

#include <gnuastro/data.h>
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




/* Macros.

 Note on `GAL_TXT_STRING_BLANK': Within Gnuastro, a blank string value is
 just a NULL pointer, but here, we are dealing with the outside world. In
 particular we want to write the column to the outside world. So this macro
 acts as a place holder for the blank string columns ONLY in ASCII
 tables.*/
#define GAL_TXT_DELIMITERS     " ,\t\f\v"
#define GAL_TXT_MAX_FMT_LENGTH 20





/* Status of a line: */
enum gal_txt_line_status_enums
{
  GAL_TXT_LINESTAT_INVALID,

  GAL_TXT_LINESTAT_BLANK,
  GAL_TXT_LINESTAT_COMMENT,
  GAL_TXT_LINESTAT_DATAROW,
};





/* Functions */
int
gal_txt_line_stat(char *line);

gal_data_t *
gal_txt_table_info(char *filename, size_t *numcols, size_t *numrows);

gal_data_t *
gal_txt_table_read(char *filename, size_t numrows, gal_data_t *colinfo,
                   struct gal_linkedlist_sll *indexll, size_t minmapsize);

gal_data_t *
gal_txt_image_read(char *filename, size_t minmapsize);

void
gal_txt_write(gal_data_t *cols, struct gal_linkedlist_stll *comment,
              char *filename, int dontdelete);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_TXT_H__ */
