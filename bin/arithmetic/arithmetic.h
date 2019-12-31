/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef IMGARITH_H
#define IMGARITH_H


#include <gnuastro/arithmetic.h>


/* These are operator codes for functions that aren't in the arithmetic
   library. */
enum arithmetic_prog_operators
{
  ARITHMETIC_OP_FILTER_MEDIAN = GAL_ARITHMETIC_OP_LAST_CODE,
  ARITHMETIC_OP_FILTER_MEAN,
  ARITHMETIC_OP_FILTER_SIGCLIP_MEAN,
  ARITHMETIC_OP_FILTER_SIGCLIP_MEDIAN,
  ARITHMETIC_OP_ERODE,
  ARITHMETIC_OP_DILATE,
  ARITHMETIC_OP_CONNECTED_COMPONENTS,
  ARITHMETIC_OP_FILL_HOLES,
  ARITHMETIC_OP_INVERT,
  ARITHMETIC_OP_INTERPOLATE_MEDIANNGB,
  ARITHMETIC_OP_COLLAPSE_SUM,
  ARITHMETIC_OP_COLLAPSE_MIN,
  ARITHMETIC_OP_COLLAPSE_MAX,
  ARITHMETIC_OP_COLLAPSE_MEAN,
  ARITHMETIC_OP_COLLAPSE_NUMBER,
  ARITHMETIC_OP_UNIQUE,
  ARITHMETIC_OP_ADD_DIMENSION,
};




void
arithmetic(struct arithmeticparams *p);

#endif
