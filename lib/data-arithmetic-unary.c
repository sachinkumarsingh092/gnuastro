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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/data.h>









/* Change input data structure type. */
gal_data_t *
data_arithmetic_change_type(gal_data_t *data, int operator,
                            unsigned char flags)
{
  int type=0;
  gal_data_t *out;

  /* Set the output type. */
  switch(operator)
    {
    case GAL_DATA_OPERATOR_TO_UCHAR:    type=GAL_DATA_TYPE_UCHAR;    break;
    case GAL_DATA_OPERATOR_TO_CHAR:     type=GAL_DATA_TYPE_UCHAR;    break;
    case GAL_DATA_OPERATOR_TO_USHORT:   type=GAL_DATA_TYPE_USHORT;   break;
    case GAL_DATA_OPERATOR_TO_SHORT:    type=GAL_DATA_TYPE_SHORT;    break;
    case GAL_DATA_OPERATOR_TO_UINT:     type=GAL_DATA_TYPE_UINT;     break;
    case GAL_DATA_OPERATOR_TO_INT:      type=GAL_DATA_TYPE_INT;      break;
    case GAL_DATA_OPERATOR_TO_ULONG:    type=GAL_DATA_TYPE_ULONG;    break;
    case GAL_DATA_OPERATOR_TO_LONG:     type=GAL_DATA_TYPE_LONG;     break;
    case GAL_DATA_OPERATOR_TO_LONGLONG: type=GAL_DATA_TYPE_LONGLONG; break;
    case GAL_DATA_OPERATOR_TO_FLOAT:    type=GAL_DATA_TYPE_FLOAT;    break;
    case GAL_DATA_OPERATOR_TO_DOUBLE:   type=GAL_DATA_TYPE_DOUBLE;   break;

    default:
      error(EXIT_FAILURE, 0, "operator value of %d not recognized in "
            "`data_arithmetic_change_type'", operator);
    }

  /* Copy to the new type. */
  out=gal_data_copy_to_new_type(data, type);

  /* Delete the input structure if the user asked for it. */
  if(flags & GAL_DATA_ARITH_FREE)
    gal_data_free(data);

  /* Return */
  return out;
}





/* Return an array of value 1 for any zero valued element and zero for any
   non-zero valued element. */
#define TYPE_CASE_FOR_NOT(TYPE, IN, IN_FINISH) {                        \
    case TYPE:                                                          \
      do *o++ = !*IN; while(++IN<IN_FINISH);                            \
      break;                                                            \
  }

gal_data_t *
data_arithmetic_not(gal_data_t *data, unsigned char flags)
{
  gal_data_t *out;

  /* 'value' will only be read from one of these based on the
     datatype. Which the caller assigned. If there is any problem, it is
     their responsability, not this function's.*/
  unsigned char     *uc = data->array,   *ucf = data->array + data->size, *o;
  char               *c = data->array,    *cf = data->array + data->size;
  unsigned short    *us = data->array,   *usf = data->array + data->size;
  short              *s = data->array,    *sf = data->array + data->size;
  unsigned int      *ui = data->array,   *uif = data->array + data->size;
  int               *in = data->array,   *inf = data->array + data->size;
  unsigned long     *ul = data->array,   *ulf = data->array + data->size;
  long               *l = data->array,    *lf = data->array + data->size;
  LONGLONG           *L = data->array,    *Lf = data->array + data->size;
  float              *f = data->array,    *ff = data->array + data->size;
  double             *d = data->array,    *df = data->array + data->size;


  /* Allocate the output array. */
  out=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, data->ndim, data->dsize,
                     data->wcs, 0, data->minmapsize);
  o=out->array;


  /* Go over the pixels and set the output values. */
  switch(data->type)
    {

    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_UCHAR,    uc,  ucf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_CHAR,     c,   cf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_LOGICAL,  c,   cf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_USHORT,   us,  usf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_SHORT,    s,   sf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_UINT,     ui,  uif)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_INT,      in,  inf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_ULONG,    ul,  ulf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_LONG,     l,   lf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_LONGLONG, L,   Lf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_FLOAT,    f,   ff)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_DOUBLE,   d,   df)

    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support bit "
            "datatype, please get in touch with us to implement it.");

    default:
      error(EXIT_FAILURE, 0, "type value (%d) not recognized "
            "in `data_arithmetic_not'", data->type);
    }

  /* Delete the input structure if the user asked for it. */
  if(flags & GAL_DATA_ARITH_FREE)
    gal_data_free(data);

  /* Return */
  return out;
}
