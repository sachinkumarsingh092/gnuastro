/*********************************************************************
qsort -- Functions used by qsort to sort an array.
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






/*****************************************************************/
/**********        Sorting of actual array        ****************/
/*****************************************************************/
int
gal_qsort_uint8_d(const void *a, const void *b);

int
gal_qsort_uint8_i(const void *a, const void *b);

int
gal_qsort_int8_d(const void *a, const void *b);

int
gal_qsort_int8_i(const void *a, const void *b);

int
gal_qsort_uint16_d(const void *a, const void *b);

int
gal_qsort_uint16_i(const void *a, const void *b);

int
gal_qsort_int16_d(const void *a, const void *b);

int
gal_qsort_int16_i(const void *a, const void *b);

int
gal_qsort_uint32_d(const void *a, const void *b);

int
gal_qsort_uint32_i(const void *a, const void *b);

int
gal_qsort_int32_d(const void *a, const void *b);

int
gal_qsort_int32_i(const void *a, const void *b);

int
gal_qsort_uint64_d(const void *a, const void *b);

int
gal_qsort_uint64_i(const void *a, const void *b);

int
gal_qsort_int64_d(const void *a, const void *b);

int
gal_qsort_int64_i(const void *a, const void *b);

int
gal_qsort_float32_d(const void *a, const void *b);

int
gal_qsort_float32_i(const void *a, const void *b);

int
gal_qsort_float64_d(const void *a, const void *b);

int
gal_qsort_float64_i(const void *a, const void *b);





/*****************************************************************/
/***************          Sorting indexs        ******************/
/*****************************************************************/
/* Pointer used to sort the indexs of an array based on their flux (value
   in this array). Note: when EACH THREAD USES A DIFFERENT ARRAY, this is
   not thread-safe . */
extern void *gal_qsort_index_single;


/* When each thread is working on a different array, we'll need to keep the
   pointer to the array in question for every index. */
struct gal_qsort_index_multi
{
  float *values; /* Array of values (pointer, so original is not touched). */
                 /* This should be identical in all elements.              */
  size_t  index; /* Index of each element to be sorted.                    */
};

int
gal_qsort_index_single_uint8_d(const void *a, const void *b);

int
gal_qsort_index_single_uint8_i(const void *a, const void *b);

int
gal_qsort_index_single_int8_d(const void *a, const void *b);

int
gal_qsort_index_single_int8_i(const void *a, const void *b);

int
gal_qsort_index_single_uint16_d(const void *a, const void *b);

int
gal_qsort_index_single_uint16_i(const void *a, const void *b);

int
gal_qsort_index_single_int16_d(const void *a, const void *b);

int
gal_qsort_index_single_int16_i(const void *a, const void *b);

int
gal_qsort_index_single_uint32_d(const void *a, const void *b);

int
gal_qsort_index_single_uint32_i(const void *a, const void *b);

int
gal_qsort_index_single_int32_d(const void *a, const void *b);

int
gal_qsort_index_single_int32_i(const void *a, const void *b);

int
gal_qsort_index_single_uint64_d(const void *a, const void *b);

int
gal_qsort_index_single_uint64_i(const void *a, const void *b);

int
gal_qsort_index_single_int64_d(const void *a, const void *b);

int
gal_qsort_index_single_int64_i(const void *a, const void *b);

int
gal_qsort_index_single_float32_d(const void *a, const void *b);

int
gal_qsort_index_single_float32_i(const void *a, const void *b);

int
gal_qsort_index_single_float64_d(const void *a, const void *b);

int
gal_qsort_index_single_float64_i(const void *a, const void *b);

int
gal_qsort_index_multi_d(const void *a, const void *b);

int
gal_qsort_index_multi_i(const void *a, const void *b);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_QSORT_H__ */
