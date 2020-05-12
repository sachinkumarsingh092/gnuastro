/*********************************************************************
Arithmetic operations on data structures.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/blank.h>

#include <gnuastro-internal/arithmetic-binary.h>






/* 'BINARY_SET_LT_INT' is defined in 'arithmetic-binary.h'. As you see
   there, this is a deep macro (calls other macros) to deal with different
   types. This allows efficiency in processing (after compilation), but
   compilation will be very slow. Therefore, for each operator we have
   defined a separate '.c' file so they are built separately and when built
   in parallel can be much faster than having them all in a single file. */
void
arithmetic_bitor(gal_data_t *l, gal_data_t *r, gal_data_t *o)
{
  /* A small sanity check. */
  if( l->type==GAL_TYPE_FLOAT32 || l->type==GAL_TYPE_FLOAT64
      || r->type==GAL_TYPE_FLOAT32 || r->type==GAL_TYPE_FLOAT64 )
      error(EXIT_FAILURE, 0, "%s: the bitor operator can only work on "
            "integer type operands", __func__);

  /* Do the operation. */
  BINARY_SET_LT_INT( ( o->type==l->type
                       ? ARITHMETIC_BINARY_OUT_TYPE_LEFT
                       : ARITHMETIC_BINARY_OUT_TYPE_RIGHT ), | );
}
