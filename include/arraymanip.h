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
#ifndef ARRAYMANIP_H
#define ARRAYMANIP_H

void
gal_arraymanip_uchar_init_on_region(unsigned char *in, const unsigned char v,
                                    size_t start, size_t s0, size_t s1,
                                    size_t is1);

void
gal_arraymanip_long_init(long *in, size_t size, const long v);

void
gal_arraymanip_long_init_on_region(long *in, const long v, size_t start,
                                   size_t s0, size_t s1, size_t is1);

void
gal_arraymanip_uchar_copy(unsigned char *in, size_t size, unsigned char **out);

void
gal_arraymanip_float_copy(float *in, size_t size, float **out);

void
gal_arraymanip_float_copy_values(float *in, size_t size, float **out);

void
gal_arraymanip_fset_const(float *in, size_t size, float a);

void
gal_arraymanip_freplace_value(float *in, size_t size, float from, float to);

void
gal_arraymanip_no_nans(float *in, size_t *size);

void
gal_arraymanip_fmultip_const(float *in, size_t size, float a);

void
gal_arraymanip_fsum_const(float *in, size_t size, float a);

float *
gal_arraymanip_fsum_arrays(float *in1, float *in2, size_t size);

#endif
