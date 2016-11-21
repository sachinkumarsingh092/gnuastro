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
#ifndef LABEL_H
#define LABEL_H





/* For the GAL_FITS_LONG_BLANK and GAL_FITS_BYTE_BLANK macros: */
#include <gnuastro/fits.h>


/* In case the blank value for the long type is negative, the only
   necessary check to see if we are on a label or not is to see if its
   value is positive. But since we don't want to impose this condition on
   the fits.h macro (which is used by all Gnuastro programs) we allow for
   it to be positive through this C preprocessor check. The compiler has no
   idea of our convention to label things with positive indices, so we
   can't rely on the compiler to optimize this. */
#if GAL_FITS_LONG_BLANK<0
#define ISINDEXABLELABEL (*lab>0)
#else
#define ISINDEXABLELABEL (*lab && *lab!=GAL_DATA_BLANK_LONG)
#endif



size_t
BF_concmp(unsigned char *byt, long *lab, size_t s0, size_t s1,
          int anyblank, const size_t connectivity);

size_t
BF_concomp_AdjMatrix(int *adj, size_t numside, long **outnewlabs);

void
removesmallarea_relabel(long *in, unsigned char *byt, size_t size,
                        size_t *numlabs, size_t minarea);

void
labindexs(long *lab, size_t size, size_t numlabs, size_t **outareas,
          size_t ***outlabinds);

#endif
