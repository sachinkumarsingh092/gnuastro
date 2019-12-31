/*********************************************************************
array - Functions to read/write arrays from/to files.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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
#ifndef __GAL_ARRAY_H__
#define __GAL_ARRAY_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>


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





/* Functions */
int
gal_array_name_recognized(char *name);

int
gal_array_name_recognized_multiext(char *name);

gal_data_t *
gal_array_read(char *filename, char *extension, gal_list_str_t *lines,
               size_t minmapsize, int quietmmap);

gal_data_t *
gal_array_read_to_type(char *filename, char *extension,
                       gal_list_str_t *lines, uint8_t type,
                       size_t minmapsize, int quietmmap);

gal_data_t *
gal_array_read_one_ch(char *filename, char *extension, gal_list_str_t *lines,
                      size_t minmapsize, int quietmmap);

gal_data_t *
gal_array_read_one_ch_to_type(char *filename, char *extension,
                              gal_list_str_t *lines, uint8_t type,
                              size_t minmapsize, int quietmmap);


__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_TIFF_H__ */
