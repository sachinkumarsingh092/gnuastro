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

#include <gnuastro/data.h>





/* Generic function for all types. */
#define COPY_OTYPE_ITYPE_SET(otype, itype) {                            \
    itype *ia=in->array;                                                \
    otype *oa=out->array, *of=oa+out->size;                             \
    do *oa=*ia++; while(++oa<of);                                       \
  }





/* Output type is set, now choose the input type. */
#define COPY_OTYPE_SET(otype)                                           \
  switch(in->type)                                                      \
    {                                                                   \
    case GAL_DATA_TYPE_UCHAR:                                           \
      COPY_OTYPE_ITYPE_SET(otype, unsigned char);                       \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_CHAR:                                            \
      COPY_OTYPE_ITYPE_SET(otype, char);                                \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_USHORT:                                          \
      COPY_OTYPE_ITYPE_SET(otype, unsigned short);                      \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_SHORT:                                           \
      COPY_OTYPE_ITYPE_SET(otype, short);                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_ULONG:                                           \
      COPY_OTYPE_ITYPE_SET(otype, unsigned long);                       \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_LONG:                                            \
      COPY_OTYPE_ITYPE_SET(otype, long);                                \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      COPY_OTYPE_ITYPE_SET(otype, LONGLONG);                            \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_FLOAT:                                           \
      COPY_OTYPE_ITYPE_SET(otype, float);                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      COPY_OTYPE_ITYPE_SET(otype, double);                              \
      break;                                                            \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized for "              \
            "for newtype in COPY_OTYPE_SET", in->type);                 \
    }





gal_data_t *
gal_data_copy_to_new_type(gal_data_t *in, int newtype)
{
  gal_data_t *out;

  /* Allocate space for the output type */
  out=gal_data_alloc(NULL, newtype, in->ndim, in->dsize, in->wcs,
                     0, in->minmapsize, in->name, in->unit, in->comment);

  /* Fill in the output array: */
  switch(newtype)
    {
    case GAL_DATA_TYPE_UCHAR:
      COPY_OTYPE_SET(unsigned char);
      break;

    case GAL_DATA_TYPE_CHAR:
      COPY_OTYPE_SET(char);
      break;

    case GAL_DATA_TYPE_USHORT:
      COPY_OTYPE_SET(unsigned short);
      break;

    case GAL_DATA_TYPE_SHORT:
      COPY_OTYPE_SET(short);
      break;

    case GAL_DATA_TYPE_ULONG:
      COPY_OTYPE_SET(unsigned long);
      break;

    case GAL_DATA_TYPE_LONG:
      COPY_OTYPE_SET(long);
      break;

    case GAL_DATA_TYPE_LONGLONG:
      COPY_OTYPE_SET(LONGLONG);
      break;

    case GAL_DATA_TYPE_FLOAT:
      COPY_OTYPE_SET(float);
      break;

    case GAL_DATA_TYPE_DOUBLE:
      COPY_OTYPE_SET(double);
      break;

    default:
      error(EXIT_FAILURE, 0, "type %d not recognized for "
            "for newtype in gal_data_copy_to_new_type", newtype);
    }

  /* Return the created array */
  return out;
}
