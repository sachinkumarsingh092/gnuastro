/*********************************************************************
SpatialConvolve - Convolve an image in the spatial domain.
This is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef SPATIALCONVOLVE_H
#define SPATIALCONVOLVE_H

int
spatialconvolve(float *input, size_t is0, size_t is1,
                float *kernel, size_t kso, size_t ks1,
                size_t numthreads, float **out);

#endif
