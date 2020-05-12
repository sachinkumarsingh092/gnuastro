/*********************************************************************
binary -- Work on binary (0 and 1 valued) datasets.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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
#ifndef __GAL_BINARY_H__
#define __GAL_BINARY_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>
#include <gnuastro/blank.h>

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



/* The binary functions will be working on a 'uint8_t' type dataset with
   values of 1 or 0 (no other pixel will be touched). However, in some
   cases, it is necessary to put temporary values in each element during
   the processing of the functions. So if your input datasets have values
   other than 0 and 1 that you don't want these functions to work on, be
   sure they are not equal to this value. It is chosen as the immediate
   value before the maximum value for this type (which is the blank value
   for this type), so blank values will also not be touched by this
   function. */
#define GAL_BINARY_TMP_VALUE GAL_BLANK_UINT8-1






/*********************************************************************/
/*****************      Erosion and dilation      ********************/
/*********************************************************************/
gal_data_t *
gal_binary_erode(gal_data_t *input, size_t num, int connectivity,
                 int inplace);

gal_data_t *
gal_binary_dilate(gal_data_t *input, size_t num, int connectivity,
                  int inplace);

gal_data_t *
gal_binary_open(gal_data_t *input, size_t num, int connectivity,
                int inplace);





/*********************************************************************/
/*****************      Connected components      ********************/
/*********************************************************************/
size_t
gal_binary_connected_components(gal_data_t *binary, gal_data_t **out,
                                int connectivity);

gal_data_t *
gal_binary_connected_indexs(gal_data_t *binary, int connectivity);

gal_data_t *
gal_binary_connected_adjacency_matrix(gal_data_t *adjacency,
                                      size_t *numnewlabs);


/*********************************************************************/
/*****************            Fill holes          ********************/
/*********************************************************************/
gal_data_t *
gal_binary_holes_label(gal_data_t *input, int connectivity,
                       size_t *numholes);

void
gal_binary_holes_fill(gal_data_t *input, int connectivity, size_t maxsize);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_BINARY_H__ */
