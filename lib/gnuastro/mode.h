/*********************************************************************
mode -- Find the mode of a distribution.
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
#ifndef __MODE_H__
#define __MODE_H__

#include <gnuastro/statistics.h>

#define MODE_GOLDEN_RATIO           1.618034f
#define MODE_TWO_TAKE_GOLDEN_RATIO  0.38197f
#define MODE_MIRROR_IS_ABOVE_RESULT (size_t)(-1)

size_t
mirrormaxdiff(float *a, size_t size, size_t m,
              size_t numcheck, size_t interval, size_t stdm);

void
modesymmetricity(float *a, size_t size, size_t mi, float errorstdm,
                 float *sym);

size_t
modegoldenselection(struct gal_statistics_mode_params *mp);

void
makemirrored(float *in, size_t mi, float **outmirror, size_t *outsize);

#endif
