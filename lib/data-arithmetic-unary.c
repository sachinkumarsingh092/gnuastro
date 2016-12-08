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




/* Return an array of value 1 for any zero valued element and zero for any
   non-zero valued element. */
#define TYPE_CASE_FOR_NOT(TYPE, IN, IN_FINISH) {                        \
    case TYPE:                                                          \
      do *o++ = !*IN; while(++IN<IN_FINISH);                            \
      break;                                                            \
  }

gal_data_t *
data_arithmetic_not(gal_data_t *data)
{
  gal_data_t *out;

  /* 'value' will only be read from one of these based on the
     datatype. Which the caller assigned. If there is any problem, it is
     their responsability, not this function's.*/
  void *A=data->array;
  size_t S=data->size;
  unsigned char     *uc = A,   *ucf = A+S, *o;
  char               *c = A,    *cf = A+S;
  unsigned short    *us = A,   *usf = A+S;
  short              *s = A,    *sf = A+S;
  unsigned int      *ui = A,   *uif = A+S;
  int               *in = A,   *inf = A+S;
  unsigned long     *ul = A,   *ulf = A+S;
  long               *l = A,    *lf = A+S;
  LONGLONG           *L = A,    *Lf = A+S;
  float              *f = A,    *ff = A+S;
  double             *d = A,    *df = A+S;


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

  /* Return */
  return out;
}
