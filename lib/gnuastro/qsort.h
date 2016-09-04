/*********************************************************************
forqsort -- Functions used by qsort to sort an array.
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
#ifndef __GAL_QSORT_H__
#define __GAL_QSORT_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */



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




/* Pointer used to sort the indexs of an array based on their flux
   (value in this array). */
extern float *gal_qsort_index_arr;

int
gal_qsort_index_float_decreasing(const void * a, const void * b);

int
gal_qsort_int_decreasing(const void * a, const void * b);

int
gal_qsort_int_increasing(const void * a, const void * b);

int
gal_qsort_float_decreasing(const void * a, const void * b);

int
gal_qsort_float_increasing(const void * a, const void * b);

int
gal_qsort_double_decreasing(const void * a, const void * b);

int
gal_qsort_double_increasing(const void * a, const void * b);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_QSORT_H__ */
