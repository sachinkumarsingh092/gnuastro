/*********************************************************************
Functions to manipulate arrays.
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
#ifndef __GAL_ARRAY_H__
#define __GAL_ARRAY_H__

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



void
gal_array_uchar_init_on_region(unsigned char *in, const unsigned char v,
                                    size_t start, size_t s0, size_t s1,
                                    size_t is1);

void
gal_array_long_init(long *in, size_t size, const long v);

void
gal_array_long_init_on_region(long *in, const long v, size_t start,
                              size_t s0, size_t s1, size_t is1);

void
gal_array_uchar_copy(unsigned char *in, size_t size, unsigned char **out);

void
gal_array_float_copy(float *in, size_t size, float **out);

void
gal_array_float_copy_values(float *in, size_t size, float **out);

void
gal_array_fset_const(float *in, size_t size, float a);

void
gal_array_freplace_value(float *in, size_t size, float from, float to);

void
gal_array_freplace_nonnans(float *in, size_t size, float to);

void
gal_array_no_nans(float *in, size_t *size);

void
gal_array_no_nans_double(double *in, size_t *size);

void
gal_array_uchar_replace(unsigned char *in, size_t size,
                        unsigned char from, unsigned char to);

void
gal_array_fmultip_const(float *in, size_t size, float a);

void
gal_array_fsum_const(float *in, size_t size, float a);

float *
gal_array_fsum_arrays(float *in1, float *in2, size_t size);


void
gal_array_dmultip_const(double *in, size_t size, double a);

void
gal_array_dmultip_arrays(double *in1, double *in2, size_t size);


void
gal_array_ddivide_const(double *in, size_t size, double a);

void
gal_array_dconst_divide(double *in, size_t size, double a);

void
gal_array_ddivide_arrays(double *in1, double *in2, size_t size);


void
gal_array_dsum_const(double *in, size_t size, double a);

void
gal_array_dsum_arrays(double *in1, double *in2, size_t size);


void
gal_array_dsubtract_const(double *in, size_t size, double a);

void
gal_array_dconst_subtract(double *in, size_t size, double a);

void
gal_array_dsubtract_arrays(double *in1, double *in2, size_t size);


void
gal_array_dpower_const(double *in, size_t size, double a);

void
gal_array_dconst_power(double *in, size_t size, double a);

void
gal_array_dpower_arrays(double *in1, double *in2, size_t size);


void
gal_array_dlog_array(double *in1, size_t size);

void
gal_array_dlog10_array(double *in1, size_t size);

void
gal_array_dabs_array(double *in1, size_t size);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_ARRAY_H__ */
