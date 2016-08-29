/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
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
#ifndef BINARY_H
#define BINARY_H





/* For the GAL_FITSARRAY_BYTE_BLANK value: */
#include <gnuastro/fits.h>





/* Special values:

     BINARYNOOP:   Value that no binary operation should be preformed on.
     BINARYTMP:    Temporary value to use within one function.

   Note that through the `setbytblank' function, practically we are using
   the value `GAL_FITS_BYTE_BLANK' (from `gnuastro/fits.h') for blank
   binary values. Recall that due to the nature of the CPU (which operates
   on 8-bits), in practice it is much more efficient to work on a byte (or
   8-bits) rather than each bit. So in practice we can use 256 values for
   meta-data analyais (like blank values, or temporary values at etc),
   eventhough the main values we are working with are 0 and 1.
*/
#define BINARYNOOP  2
#define BINARYTMP   3





void
setbytblank(float *img, unsigned char *byt, size_t size);

void
count_f_b_onregion(unsigned char *byt, size_t startind, size_t s0,
                   size_t s1, size_t is1, size_t *numf, size_t *numb,
                   int *anyblank);

void
index_f_b_onregion(unsigned char *byt, size_t startind, size_t s0,
                   size_t s1, size_t is1, size_t *inds,
                   unsigned char b0f1);

void
dilate0_erode1_4con(unsigned char *byt, size_t nr, size_t nc,
                    unsigned char b0_f1);

void
dilate0_erode1_8con(unsigned char *byt, size_t nr, size_t nc,
                    unsigned char b0_f1);

void
opening(unsigned char *byt, size_t s0, size_t s1,
        size_t depth, int con_type);

void
fillboundedholes(unsigned char *in, size_t s0, size_t s1, int anyblank);

void
maskbackorforeground(float *in, size_t size, unsigned char *byt,
                     unsigned char b0f1);

#endif
