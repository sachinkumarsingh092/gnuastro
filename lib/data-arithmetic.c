/*********************************************************************
Arithmetic operations on data structures.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>
#include <stdarg.h>
#include <stdlib.h>

#include <gnuastro/data.h>

#include <data-arithmetic-binary.h>





gal_data_t *
gal_data_arithmetic(int operator, unsigned char flags, ...)
{
  va_list va;
  int out_type;
  size_t out_size;
  gal_data_t *o=NULL;

  /* Prepare the variable arguments (starting after the flags argument). */
  va_start(va, flags);

  /* Depending on the operator do the job: */
  switch(operator)
    {
    case GAL_DATA_OPERATOR_PLUS:     BINARY_INTERNAL(+, 0); break;
    case GAL_DATA_OPERATOR_MINUS:    BINARY_INTERNAL(-,  0); break;
    case GAL_DATA_OPERATOR_MULTIPLY: BINARY_INTERNAL(*,  0); break;
    case GAL_DATA_OPERATOR_DIVIDE:   BINARY_INTERNAL(/,  0); break;

    case GAL_DATA_OPERATOR_LT:  BINARY_INTERNAL(<,  GAL_DATA_TYPE_UCHAR); break;
    case GAL_DATA_OPERATOR_LE:  BINARY_INTERNAL(<=, GAL_DATA_TYPE_UCHAR); break;
    case GAL_DATA_OPERATOR_GT:  BINARY_INTERNAL(>,  GAL_DATA_TYPE_UCHAR); break;
    case GAL_DATA_OPERATOR_GE:  BINARY_INTERNAL(>=, GAL_DATA_TYPE_UCHAR); break;
    case GAL_DATA_OPERATOR_EQ:  BINARY_INTERNAL(==, GAL_DATA_TYPE_UCHAR); break;
    case GAL_DATA_OPERATOR_NE:  BINARY_INTERNAL(!=, GAL_DATA_TYPE_UCHAR); break;
    case GAL_DATA_OPERATOR_AND: BINARY_INTERNAL(&&, GAL_DATA_TYPE_UCHAR); break;
    case GAL_DATA_OPERATOR_OR:  BINARY_INTERNAL(||, GAL_DATA_TYPE_UCHAR); break;

#if 0
  else if(!strcmp(operator, "abs"))       takeabs(p);
  else if(!strcmp(operator, "pow"))       topower(p, NULL);
  else if(!strcmp(operator, "sqrt"))      takesqrt(p);
  else if(!strcmp(operator, "log"))       takelog(p);
  else if(!strcmp(operator, "log10"))     takelog10(p);
  else if(!strcmp(operator, "minvalue"))  findmin(p);
  else if(!strcmp(operator, "maxvalue"))  findmax(p);
  else if(!strcmp(operator, "min")
          || !strcmp(operator, "max")
          || !strcmp(operator, "average")
          || !strcmp(operator, "median")) alloppixs(p, operator);
  else if(!strcmp(operator, "not"))       notfunc(p);
  else if(!strcmp(operator, "isblank"))   opisblank(p);
  else if(!strcmp(operator, "where"))     where(p);
#endif

    default:
      error(EXIT_FAILURE, 0, "the argument \"%d\" could not be "
            "interpretted as an operator", operator);
    }

  /* End the variable argument structure and return. */
  va_end(va);
  return o;
}
