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
#include "fitsarrayvv.h"





/* Special values: */
#define BINARYTMP 2





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
